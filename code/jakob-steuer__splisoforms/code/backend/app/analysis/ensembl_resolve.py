"""Resolve Ensembl IDs (ENST/ENSP/ENSG) that are NOT in our DB to their gene.

Powers the search "not in database — did you mean…" fallback. When a user pastes
an Ensembl ID from a release we didn't ingest, we map it to its gene via the
Ensembl REST API and, if that gene exists in our DB, point the user at it.

Design:
  * IDs are matched VERSION-STRIPPED (our own IDs are unversioned).
  * Answers are cached in `ensembl_id_cache`, negatives included, so a repeated
    miss never re-hits Ensembl.
  * We match our own `genes.ensembl_gene_id` directly by ENSG — so the common
    "gene is in our DB" case costs at most one or two REST calls, and the
    Ensembl display name is only fetched when the gene is NOT in our DB.
"""
import hashlib
import logging
import re
from datetime import datetime
from difflib import SequenceMatcher
from typing import Optional

import httpx
from sqlalchemy.orm import Session

from app.models import bio as m

log = logging.getLogger("ensembl_resolve")

_ENSEMBL_REST = "https://rest.ensembl.org"
_TIMEOUT = 5.0
# Above this length we skip the O(n^2) difflib pass (a titin-scale query would
# stall the request); the exact-md5 check still runs, so identical proteins are
# always caught — we just don't report an approximate % for giant sequences.
_MAX_DIFFLIB_LEN = 12000


def _identity(a: str, b: str) -> float:
    """Normalised similarity ∈ [0,1] (difflib, length-prefiltered) — the same
    alignment-free measure backfill_canonical_match.py uses for canonical linking."""
    if not a or not b:
        return 0.0
    min_len, max_len = min(len(a), len(b)), max(len(a), len(b))
    if max_len == 0:
        return 0.0
    if min_len / max_len < 0.2:
        return min_len / max_len
    return SequenceMatcher(a=a, b=b, autojunk=False).ratio()

# ENSG/ENST/ENSP (+ optional species infix like ENSMUS...) with a numeric body and
# an optional .N version. We only care about human here but the pattern is generic.
_ENS_RE = re.compile(r"^(ENS[A-Z]*[GTP])(\d{6,})(?:\.\d+)?$", re.IGNORECASE)


def looks_like_ensembl(raw: str) -> bool:
    return bool(_ENS_RE.match((raw or "").strip()))


def strip_version(raw: str) -> str:
    s = (raw or "").strip()
    mt = _ENS_RE.match(s)
    return (mt.group(1) + mt.group(2)).upper() if mt else s


def id_type(raw: str) -> Optional[str]:
    mt = _ENS_RE.match((raw or "").strip())
    if not mt:
        return None
    prefix = mt.group(1).upper()
    if prefix.endswith("G"):
        return "ensg"
    if prefix.endswith("T"):
        return "enst"
    if prefix.endswith("P"):
        return "ensp"
    return None


def _rest_get(client: httpx.Client, path: str) -> Optional[dict]:
    try:
        r = client.get(f"{_ENSEMBL_REST}{path}",
                        headers={"Accept": "application/json"}, timeout=_TIMEOUT)
        if r.status_code == 200:
            return r.json()
        return None
    except (httpx.HTTPError, ValueError):
        return None


def _ensg_from_ensembl(client: httpx.Client, base: str, itype: str) -> tuple[Optional[str], Optional[str]]:
    """Return (ensg, ensembl_symbol) for an ID, via the fewest REST hops.

    ENSG → itself; ENST → its Parent gene; ENSP → transcript → gene.
    Falls back to /archive/id for retired IDs (uses the possible replacement).
    """
    obj = _rest_get(client, f"/lookup/id/{base}")
    if obj is None:
        # Retired / unknown current ID — consult the archive for a replacement.
        arc = _rest_get(client, f"/archive/id/{base}")
        repl = (arc or {}).get("possible_replacement") or []
        new_id = None
        for entry in repl:
            new_id = entry.get("stable_id") or entry.get("id")
            if new_id:
                break
        if not new_id:
            return None, None
        base = strip_version(new_id)
        itype = id_type(base) or itype
        obj = _rest_get(client, f"/lookup/id/{base}")
        if obj is None:
            return None, None

    ot = (obj.get("object_type") or "").lower()
    if ot == "gene" or itype == "ensg":
        return obj.get("id") or base, obj.get("display_name")
    if ot == "transcript" or itype == "enst":
        return obj.get("Parent"), None
    if ot == "translation" or itype == "ensp":
        enst = obj.get("Parent")
        if not enst:
            return None, None
        tobj = _rest_get(client, f"/lookup/id/{enst}")
        return (tobj or {}).get("Parent"), None
    return None, None


def _match_our_gene(db: Session, ensg: Optional[str], symbol: Optional[str]) -> Optional[m.Gene]:
    gene = None
    if ensg:
        gene = db.query(m.Gene).filter(m.Gene.ensembl_gene_id == ensg).first()
    if not gene and symbol:
        gene = db.query(m.Gene).filter(m.Gene.symbol == symbol).first()
    return gene


def _fetch_protein_seq(client: httpx.Client, base: str, itype: str) -> tuple[Optional[str], bool]:
    """Fetch the protein sequence for an ENST/ENSP from Ensembl.

    Returns (seq, definitive). `definitive` is True when Ensembl gave a real answer
    (a sequence, or a clear "no protein" for a non-coding id) — so the caller can
    stop re-fetching. It is False on transient failures (rate-limit / 5xx / timeout)
    so those are retried on a later request rather than cached as an empty answer.
    """
    if itype not in ("enst", "ensp"):
        return None, True
    try:
        r = client.get(f"{_ENSEMBL_REST}/sequence/id/{base}?type=protein",
                        headers={"Accept": "application/json"}, timeout=_TIMEOUT)
    except httpx.HTTPError:
        return None, False  # transient (timeout / network)
    if r.status_code == 200:
        try:
            seq = (r.json() or {}).get("seq")
        except ValueError:
            return None, True
        if not seq:
            return None, True
        seq = re.sub(r"[^A-Za-z]", "", seq).upper()
        return (seq or None), True
    if r.status_code == 400 or r.status_code == 404:
        return None, True   # Ensembl: no protein for this id (e.g. non-coding)
    if r.status_code == 429 or r.status_code >= 500:
        return None, False  # rate-limited / server error → retry later
    return None, False


def closest_isoform(db: Session, gene_id: Optional[int], query_seq: Optional[str]) -> dict:
    """Find the DB isoform most similar to `query_seq`, computed LIVE against the
    current DB (never cached, so newly-ingested isoforms are always reflected).

    Strategy:
      1. Exact protein match by md5 across ALL isoforms — catches the same protein
         under a transcript id we didn't ingest, including one added by a future
         dataset → reported as identity 1.0 / exact.
      2. Otherwise, the highest-identity isoform WITHIN the resolved gene (a small
         candidate set), via the repo's difflib measure.
    Returns {} when nothing usable is found.
    """
    if not query_seq:
        return {}

    qmd5 = hashlib.md5(query_seq.encode()).hexdigest()
    exact = (
        db.query(m.Isoform)
        .filter(m.Isoform.protein_md5 == qmd5)
        .join(m.Gene, m.Isoform.gene_id == m.Gene.id)
        .first()
    )
    if exact is not None:
        return {
            "closest_isoform_id": exact.id,
            "closest_transcript_id": exact.transcript_id,
            "closest_transcript_name": exact.transcript_name,
            "closest_gene_id": exact.gene_id,
            "closest_gene_symbol": exact.gene.symbol if exact.gene else None,
            "closest_identity": 1.0,
            "closest_exact": True,
        }

    if not gene_id or len(query_seq) > _MAX_DIFFLIB_LEN:
        return {}

    candidates = (
        db.query(m.Isoform)
        .filter(m.Isoform.gene_id == gene_id, m.Isoform.sequence.isnot(None))
        .all()
    )
    best, best_r = None, 0.0
    for iso in candidates:
        if len(iso.sequence) > _MAX_DIFFLIB_LEN:
            continue
        r = _identity(query_seq, iso.sequence)
        if r > best_r:
            best_r, best = r, iso
    if best is None:
        return {}
    return {
        "closest_isoform_id": best.id,
        "closest_transcript_id": best.transcript_id,
        "closest_transcript_name": best.transcript_name,
        "closest_gene_id": best.gene_id,
        "closest_gene_symbol": best.gene.symbol if best.gene else None,
        "closest_identity": round(best_r, 4),
        "closest_exact": False,
    }


def resolve(db: Session, raw_id: str) -> dict:
    """Resolve an Ensembl ID to our gene (cached). Returns a plain dict:

    status ∈ {resolved_gene, gene_not_in_db, not_found, not_ensembl}
    plus gene_id / gene_symbol / resolved_ensg / isoform_count / message.
    """
    base = strip_version(raw_id)
    itype = id_type(raw_id)
    if not itype:
        return {"status": "not_ensembl", "normalized_id": base, "id_type": None,
                "message": "Not an Ensembl ID."}

    cached = db.query(m.EnsemblIdCache).filter(m.EnsemblIdCache.ensembl_id == base).first()

    if cached is None:
        ensg = None
        symbol = None
        status = "not_found"
        with httpx.Client() as client:
            ensg, symbol = _ensg_from_ensembl(client, base, itype)
            if ensg:
                status = "found"
                gene = _match_our_gene(db, ensg, symbol)
                # Only spend a call on the Ensembl display name when we can't get a
                # symbol from our own gene row.
                if gene is None and not symbol:
                    gobj = _rest_get(client, f"/lookup/id/{ensg}")
                    symbol = (gobj or {}).get("display_name")
        # Upsert cache (best-effort; a race just means a duplicate lookup).
        gene = _match_our_gene(db, ensg, symbol) if ensg else None
        cached = m.EnsemblIdCache(
            ensembl_id=base, id_type=itype, resolved_ensg=ensg,
            gene_symbol=(gene.symbol if gene else symbol),
            resolved_gene_id=(gene.id if gene else None),
            status=status, checked_at=datetime.utcnow(),
        )
        try:
            db.add(cached)
            db.commit()
        except Exception:
            db.rollback()
            cached = db.query(m.EnsemblIdCache).filter(m.EnsemblIdCache.ensembl_id == base).first()
            if cached is None:
                raise

    # Resolve to our gene LIVE from the stable external key (the ENSG), so a gene
    # newly ingested — or re-attributed in a future rebuild — is always reflected.
    # The cached resolved_gene_id is only a fallback hint.
    gene = None
    if cached.resolved_ensg:
        gene = _match_our_gene(db, cached.resolved_ensg, cached.gene_symbol)
    if gene is None and cached.resolved_gene_id:
        gene = db.query(m.Gene).filter(m.Gene.id == cached.resolved_gene_id).first()

    if cached.status == "not_found":
        return {"status": "not_found", "normalized_id": base, "id_type": itype,
                "resolved_ensg": None,
                "message": "This ID is not in the database and could not be resolved via Ensembl."}

    # Lazily fetch + cache the query protein sequence (external, stable). Backfills
    # rows created before this feature, and skips ids already checked (incl. ones
    # Ensembl has no protein for). Never blocks the resolution if the fetch fails.
    if not cached.protein_seq_checked and itype in ("enst", "ensp"):
        try:
            with httpx.Client() as client:
                qseq, definitive = _fetch_protein_seq(client, base, itype)
            if qseq is not None:
                cached.query_protein_seq = qseq
                cached.query_protein_md5 = hashlib.md5(qseq.encode()).hexdigest()
                cached.protein_seq_checked = True
                db.commit()
            elif definitive:
                # Ensembl has no protein for this id — mark checked so we don't retry.
                cached.protein_seq_checked = True
                db.commit()
            # else: transient failure — leave unchecked to retry on a later request.
        except Exception:
            db.rollback()

    if gene is not None:
        n_iso = db.query(m.Isoform).filter(m.Isoform.gene_id == gene.id).count()
        # Closest isoform is ALWAYS computed live against the current DB from the
        # cached query sequence — never persisted — so a future ingest of this (or a
        # closer) isoform is reflected on the very next request.
        closest = closest_isoform(db, gene.id, cached.query_protein_seq)
        out = {
            "status": "resolved_gene", "normalized_id": base, "id_type": itype,
            "resolved_ensg": cached.resolved_ensg, "gene_id": gene.id,
            "gene_symbol": gene.symbol, "isoform_count": n_iso,
            "message": (f"{base} isn't stored directly, but it belongs to "
                        f"{gene.symbol} — showing that gene's {n_iso} isoform(s)."),
        }
        out.update(closest)
        return out

    return {
        "status": "gene_not_in_db", "normalized_id": base, "id_type": itype,
        "resolved_ensg": cached.resolved_ensg, "gene_symbol": cached.gene_symbol,
        "message": (f"{base} maps to gene "
                    f"{cached.gene_symbol or cached.resolved_ensg or 'unknown'}, "
                    "which is not in this database."),
    }

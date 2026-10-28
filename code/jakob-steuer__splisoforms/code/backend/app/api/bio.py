from fastapi import APIRouter, Depends, HTTPException, Query, Response
from sqlalchemy import or_, func, text, case
from sqlalchemy.orm import Session, joinedload, selectinload
from typing import List, Optional
import os
import re
import time
import tempfile
import threading
import hashlib
import httpx

from app.database import get_db
from app.models import bio as bio_models
from app.schemas import bio as bio_schemas
from app.analysis import ensembl_resolve

router = APIRouter()

# On-disk cache for proxied STRING network images. We serve STRING's own
# node-edge SVG same-origin (the site CSP blocks external img/connect), and
# cache each gene's graph so we don't re-hit the STRING API on every view.
# Lives under the system temp dir so it is always writable at request time —
# the prod backend runs on a read-only rootfs, so /app/... would reject writes.
# Ephemeral (cleared on restart) is fine: entries simply re-fetch on demand.
_STRING_NETWORK_CACHE_DIR = os.path.join(tempfile.gettempdir(), "splisoforms_string_network_cache")
_STRING_NETWORK_TTL_S = 30 * 24 * 3600  # 30 days — STRING networks are stable

# In-process cache for the /neoantigens/top total_count + stats block. Both depend only
# on the filter set (not on limit/offset/sort), so without it every page flip of the
# same atlas view re-aggregates the ~2.6M-row representative set (~0.5 s with the
# partial covering indexes, seconds without). Entries are a dozen ints keyed by the
# exact WHERE clause + bound params. Staleness only matters when a neoantigen backfill
# runs against a live server — deploys restart the process anyway. 0 disables.
_ATLAS_STATS_TTL_S = int(os.getenv("ATLAS_STATS_TTL_S", "600"))
_ATLAS_STATS_MAX_ENTRIES = 256
_atlas_stats_cache: dict = {}          # key -> (expires_at [monotonic], total_count, stats)
_atlas_stats_lock = threading.Lock()


def _atlas_stats_get(key):
    now = time.monotonic()
    with _atlas_stats_lock:
        hit = _atlas_stats_cache.get(key)
        if hit is None:
            return None
        expires_at, total_count, stats = hit
        if expires_at <= now:
            _atlas_stats_cache.pop(key, None)
            return None
        return total_count, dict(stats)


def _atlas_stats_put(key, total_count, stats):
    if _ATLAS_STATS_TTL_S <= 0:
        return
    with _atlas_stats_lock:
        if len(_atlas_stats_cache) >= _ATLAS_STATS_MAX_ENTRIES:
            # insertion order ≈ age: drop the older half rather than track LRU
            for k in list(_atlas_stats_cache)[: _ATLAS_STATS_MAX_ENTRIES // 2]:
                _atlas_stats_cache.pop(k, None)
        _atlas_stats_cache[key] = (time.monotonic() + _ATLAS_STATS_TTL_S, total_count, dict(stats))

@router.get("/genes/search", response_model=bio_schemas.SearchResponse)
def search_genes(query: str = Query(..., min_length=1), source_id: int = Query(None), limit: int = Query(25, le=100), offset: int = Query(0, ge=0), db: Session = Depends(get_db)):
    """Search by Gene Symbol, ENSG, ENST, ENSP, and UniProt IDs with Smart Redirects."""
    search_term = f"%{query}%"

    # Ensembl IDs are stored UNVERSIONED (e.g. ENST00000437963). A user may paste a
    # versioned id from any release (ENST00000437963.9) — strip the .N suffix so it
    # still matches. Applies to the ID-shaped fields only (ENSG/ENST/ENSP/transcript).
    _ens = re.match(r'^(ENS[A-Z]*\d{6,})\.\d+$', query.strip(), re.IGNORECASE)
    id_query = _ens.group(1) if _ens else query.strip()
    id_term = f"%{id_query}%"

    # Base Join Structure for Genes
    base_gene_query = db.query(bio_models.Gene).outerjoin(
        bio_models.Isoform, bio_models.Gene.id == bio_models.Isoform.gene_id
    ).outerjoin(
        bio_models.TranscriptMetadata, bio_models.Isoform.id == bio_models.TranscriptMetadata.isoform_id
    )

    # Resolve Genes
    gene_query = base_gene_query
    if query.startswith("NM_"):
        gene_query = gene_query.filter(bio_models.TranscriptMetadata.refseq_id.ilike(search_term))
    elif query.upper().startswith("HGNC:"):
        hgnc_clean = query.upper().replace("HGNC:", "").strip()
        gene_query = gene_query.filter(bio_models.Gene.hgnc_symbol.ilike(f"%{hgnc_clean}%"))
    else:
        gene_query = gene_query.filter(
            or_(
                bio_models.Gene.symbol.ilike(search_term),
                bio_models.Gene.hgnc_symbol.ilike(search_term),
                bio_models.Gene.uniprot_id.ilike(search_term),
                # ENSG lookup (version-stripped) — the docstring promised it but the
                # filter never included it.
                bio_models.Gene.ensembl_gene_id.ilike(id_term),
                # Protein name (UniProt "recommended name"), so "transferrin receptor"
                # finds TFRC. Aliases/CD antigens (TfR1, CD71) are not stored yet.
                bio_models.Gene.protein_name.ilike(search_term),
            )
        )
    
    if source_id:
        gene_query = gene_query.filter(bio_models.Isoform.source_id == source_id)

    total_genes = gene_query.distinct().count()
    # Relevance order — exact symbol, then symbol prefix, then exact accession/ENSG,
    # then anything else; alphabetical within a rank. Without it a 2-letter symbol
    # such as "F2" (the paper's thrombin example) is buried under 400+ substring
    # matches (OR4F29, RNF223, …) and never reaches the first page. GROUP BY the PK
    # instead of DISTINCT so the CASE can be an ORDER BY expression.
    gene_rank = case(
        (bio_models.Gene.symbol.ilike(query), 0),
        (bio_models.Gene.hgnc_symbol.ilike(query), 0),
        (bio_models.Gene.symbol.ilike(f"{query}%"), 1),
        (bio_models.Gene.uniprot_id.ilike(query), 2),
        (bio_models.Gene.ensembl_gene_id.ilike(id_query), 2),
        else_=3,
    )
    genes_raw = (gene_query.group_by(bio_models.Gene.id)
                 .order_by(gene_rank, bio_models.Gene.symbol.asc())
                 .limit(limit).offset(offset).all())

    # Search Transcripts (Isoforms)
    trans_query = db.query(bio_models.Isoform).join(
        bio_models.Gene, bio_models.Isoform.gene_id == bio_models.Gene.id
    ).outerjoin(
        bio_models.TranscriptMetadata, bio_models.Isoform.id == bio_models.TranscriptMetadata.isoform_id
    )
    
    if source_id:
        trans_query = trans_query.filter(bio_models.Isoform.source_id == source_id)
        
    transcripts_query = trans_query.filter(
        or_(
            bio_models.Isoform.transcript_id.ilike(id_term),
            bio_models.Isoform.transcript_name.ilike(search_term),
            bio_models.Isoform.ensembl_transcript_id.ilike(id_term),
            bio_models.TranscriptMetadata.uniprot_id.ilike(search_term),
            bio_models.TranscriptMetadata.enst_id.ilike(id_term),
            bio_models.TranscriptMetadata.ensp_id.ilike(id_term),
        )
    )
    total_transcripts = transcripts_query.distinct().count()
    transcript_rank = case(
        (bio_models.Isoform.transcript_name.ilike(query), 0),
        (bio_models.Isoform.transcript_id.ilike(id_query), 0),
        (bio_models.Isoform.ensembl_transcript_id.ilike(id_query), 0),
        (bio_models.Isoform.transcript_name.ilike(f"{query}%"), 1),
        (bio_models.Isoform.transcript_name.ilike(f"{query}-%"), 1),
        else_=2,
    )
    transcripts_raw = (transcripts_query.group_by(bio_models.Isoform.id)
                       .order_by(transcript_rank, bio_models.Isoform.transcript_name.asc())
                       .limit(limit).offset(offset).all())

    # Format Results
    genes_out = []
    for g in genes_raw:
        genes_out.append({
            "match_type": "gene",
            "id": g.id,
            "gene_symbol": g.symbol,
            "uniprot_id": g.uniprot_id,
            "description": g.description
        })

    transcripts_out = []
    for t in transcripts_raw:
        source_name = getattr(t.source_dataset, 'name', 'Unknown')
        transcripts_out.append({
            "match_type": "transcript",
            "id": t.gene_id,
            "isoform_id": t.id,
            "transcript_id": t.transcript_id,
            "transcript_name": t.transcript_name,
            "gene_symbol": t.gene.symbol if t.gene else None,
            "source_name": source_name
        })

    # Drugs — resolved as TARGET–DRUG PAIRS, not bare drug rows: a drug is only meaningful
    # here in the context of the protein it was co-folded against, and the row links to
    # /drugs/{gene_id}. Matches from either side (drug name / PubChem CID, or gene symbol),
    # so "abiraterone" and "AR" both surface the AR + abiraterone pair.
    pair_query = (
        db.query(
            bio_models.Gene.id.label("gene_id"),
            bio_models.Gene.symbol.label("gene_symbol"),
            bio_models.Drug.name.label("drug_name"),
            bio_models.Drug.pubchem_cid.label("pubchem_cid"),
            bio_models.Drug.moa.label("moa"),
            func.count(bio_models.DrugComplex.id).label("n_complexes"),
        )
        .join(bio_models.DrugComplex, bio_models.DrugComplex.drug_id == bio_models.Drug.id)
        .join(bio_models.Gene, bio_models.Gene.id == bio_models.DrugComplex.gene_id)
        .filter(
            or_(
                bio_models.Drug.name.ilike(search_term),
                bio_models.Drug.pubchem_cid == query,
                # Gene symbols are short, so an EXACT match — a substring match on "AR"
                # would drag in ARRB1, ARID1A, … and bury the real hit.
                bio_models.Gene.symbol.ilike(query),
            )
        )
        .group_by(
            bio_models.Gene.id,
            bio_models.Gene.symbol,
            bio_models.Drug.name,
            bio_models.Drug.pubchem_cid,
            bio_models.Drug.moa,
        )
        # Relevance, most specific first: an exact drug name, then an exact gene symbol (a
        # bare "AR" means the receptor, not the compound "AR-12"), then a drug-name prefix,
        # then a loose substring. Ties break on evidence volume. Without this a 2-char query
        # ranks "flunARizine" alongside the real hit.
        .order_by(
            case(
                (bio_models.Drug.name.ilike(query), 0),
                (bio_models.Gene.symbol.ilike(query), 1),
                (bio_models.Drug.name.ilike(f"{query}%"), 2),
                else_=3,
            ),
            func.count(bio_models.DrugComplex.id).desc(),
            bio_models.Drug.name.asc(),
        )
    )
    total_drugs = pair_query.count()
    drugs_out = [
        {
            "match_type": "drug",
            "id": r.gene_id,
            "gene_symbol": r.gene_symbol,
            "drug_name": r.drug_name,
            "pubchem_cid": r.pubchem_cid,
            "moa": r.moa,
            "n_complexes": r.n_complexes,
        }
        for r in pair_query.limit(limit).offset(offset).all()
    ]

    return {
        "genes": genes_out,
        "transcripts": transcripts_out,
        "drugs": drugs_out,
        "total_genes": total_genes,
        "total_transcripts": total_transcripts,
        "total_drugs": total_drugs,
    }

# Stored protein sequences are uppercase, whitespace-free, and carry no stop char;
# `protein_md5 = md5(sequence.encode())` on exactly that form (verified against the DB).
# We normalise user input the same way so an exact paste hashes to the same digest.
_PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBXZUO")  # 20 aa + ambiguity/rare (B Z X U O)

def _sanitize_protein(raw: str) -> str:
    """FASTA header lines dropped, all whitespace + stop markers removed, uppercased."""
    lines = [ln for ln in (raw or "").splitlines() if not ln.lstrip().startswith(">")]
    seq = "".join(lines).upper()
    seq = re.sub(r"\s", "", seq)
    seq = seq.replace("*", "").replace(".", "").replace("-", "")
    return seq

@router.post("/search/sequence", response_model=bio_schemas.SequenceSearchResponse)
def search_by_sequence(payload: bio_schemas.SequenceSearchRequest, db: Session = Depends(get_db)):
    """Exact protein-sequence search (Phase 1): normalise → md5 → protein_md5 index.

    Returns every isoform whose stored protein sequence is byte-identical to the query.
    (Phase 2 will add MMseqs2 similarity hits under a separate `similar` field.)
    """
    seq = _sanitize_protein(payload.sequence)
    if len(seq) < 10:
        raise HTTPException(status_code=400, detail="Sequence too short — at least 10 residues required.")
    if len(seq) > 40000:
        raise HTTPException(status_code=400, detail="Sequence too long — maximum 40,000 residues.")
    bad = sorted(set(seq) - _PROTEIN_ALPHABET)
    if bad:
        raise HTTPException(status_code=400,
                            detail=f"Not a protein sequence — unexpected characters: {''.join(bad)[:20]}")

    md5 = hashlib.md5(seq.encode()).hexdigest()

    rows = (
        db.query(bio_models.Isoform)
        .filter(bio_models.Isoform.protein_md5 == md5)
        .join(bio_models.Gene, bio_models.Isoform.gene_id == bio_models.Gene.id)
        .all()
    )

    exact = [
        bio_schemas.SequenceHit(
            match_type="exact",
            isoform_id=iso.id,
            gene_id=iso.gene_id,
            transcript_id=iso.transcript_id,
            transcript_name=iso.transcript_name,
            gene_symbol=iso.gene.symbol if iso.gene else None,
            source_name=getattr(iso.source_dataset, "name", None),
            seq_length=len(iso.sequence) if iso.sequence else len(seq),
            identity=1.0,
        )
        for iso in rows
    ]

    msg = None
    if not exact:
        msg = ("No isoform has this exact sequence. Similarity search "
               "(nearest matches) is coming soon.")
    return bio_schemas.SequenceSearchResponse(
        query_length=len(seq), query_md5=md5, exact=exact, message=msg,
    )

@router.get("/search/resolve", response_model=bio_schemas.ResolveResponse)
def resolve_ensembl_id(id: str = Query(..., min_length=1), db: Session = Depends(get_db)):
    """Resolve an Ensembl ID that isn't stored directly to its gene (Phase 3).

    The frontend calls this only when a normal ID search returns nothing AND the
    query looks like an Ensembl ID. We version-strip, check the DB once more, then
    fall back to a cached Ensembl REST lookup → "belongs to gene X" suggestion.
    """
    raw = id.strip()
    if not ensembl_resolve.looks_like_ensembl(raw):
        return bio_schemas.ResolveResponse(
            query=raw, status="not_ensembl", message="Not an Ensembl ID.")

    base = ensembl_resolve.strip_version(raw)
    itype = ensembl_resolve.id_type(raw)
    like = f"%{base}%"

    # In-DB check (defensive — the base id may be present even if the versioned form
    # missed). Resolve to the owning gene so the row still links somewhere useful.
    iso = (
        db.query(bio_models.Isoform)
        .outerjoin(bio_models.TranscriptMetadata,
                   bio_models.Isoform.id == bio_models.TranscriptMetadata.isoform_id)
        .filter(
            or_(
                bio_models.Isoform.transcript_id.ilike(like),
                bio_models.Isoform.ensembl_transcript_id.ilike(like),
                bio_models.TranscriptMetadata.enst_id.ilike(like),
                bio_models.TranscriptMetadata.ensp_id.ilike(like),
            )
        )
        .first()
    )
    if iso is not None:
        return bio_schemas.ResolveResponse(
            query=raw, normalized_id=base, id_type=itype, status="in_db",
            gene_id=iso.gene_id, gene_symbol=(iso.gene.symbol if iso.gene else None),
            message="Already in the database.")
    gene_direct = db.query(bio_models.Gene).filter(
        bio_models.Gene.ensembl_gene_id.ilike(like)).first()
    if gene_direct is not None:
        n_iso = db.query(bio_models.Isoform).filter(
            bio_models.Isoform.gene_id == gene_direct.id).count()
        return bio_schemas.ResolveResponse(
            query=raw, normalized_id=base, id_type=itype, status="in_db",
            gene_id=gene_direct.id, gene_symbol=gene_direct.symbol,
            resolved_ensg=gene_direct.ensembl_gene_id, isoform_count=n_iso,
            message="Already in the database.")

    # Not in DB → cached Ensembl REST resolution.
    try:
        result = ensembl_resolve.resolve(db, raw)
    except Exception as exc:  # network / parse / DB — degrade gracefully
        return bio_schemas.ResolveResponse(
            query=raw, normalized_id=base, id_type=itype, status="error",
            message=f"Could not resolve this ID right now ({type(exc).__name__}).")

    return bio_schemas.ResolveResponse(query=raw, **result)

@router.get("/genes/{gene_id}", response_model=bio_schemas.GeneDetailResponse)
def get_gene(gene_id: str, db: Session = Depends(get_db)):
    """
    Detailed view for a single Gene, including its overview metadata and isoforms.

    Accepts either the integer DB id (legacy, e.g. `/api/genes/12345`) OR the
    HGNC symbol (e.g. `/api/genes/BRCA1`). The frontend route `/gene/[id]`
    passes whatever the user typed/navigated to straight through, so the
    backend must handle both.
    """
    q = db.query(bio_models.Gene).options(
        joinedload(bio_models.Gene.isoforms).joinedload(bio_models.Isoform.structures),
        joinedload(bio_models.Gene.isoforms).joinedload(bio_models.Isoform.metadata_record),
        joinedload(bio_models.Gene.isoforms).joinedload(bio_models.Isoform.source_dataset),
        joinedload(bio_models.Gene.isoforms).selectinload(bio_models.Isoform.pfam_domains),
        joinedload(bio_models.Gene.isoforms).joinedload(bio_models.Isoform.nmd_prediction),
        joinedload(bio_models.Gene.isoforms).selectinload(
            bio_models.Isoform.source_memberships
        ).joinedload(bio_models.IsoformSourceMembership.source),
    )
    if gene_id.isdigit():
        gene = q.filter(bio_models.Gene.id == int(gene_id)).first()
    else:
        gene = q.filter(bio_models.Gene.symbol == gene_id).first()
    if not gene:
        raise HTTPException(status_code=404, detail="Gene not found")

    for isoform in gene.isoforms:
        # Canonical flag fallback
        if isoform.is_canonical is None:
            if gene.canonical_transcript_id and isoform.id == gene.canonical_transcript_id:
                isoform.is_canonical = True
            elif gene.canonical_sequence and isoform.sequence and isoform.sequence == gene.canonical_sequence:
                isoform.is_canonical = True
            else:
                isoform.is_canonical = False
        # Hydrate into a non-ORM attribute so Pydantic can read it without SQLAlchemy
        # intercepting the assignment (source_memberships is an ORM relationship).
        setattr(isoform, "dataset_presence", [
            bio_schemas.IsoformSourceMembershipResponse(
                source_id=m.source_id,
                source_name=m.source.name if m.source else "Unknown",
                in_normal=m.in_normal,
                in_primary_tumor=m.in_primary_tumor,
                in_metastasis=m.in_metastasis,
            )
            for m in isoform.source_memberships
        ])

    return gene

@router.get("/genes/{gene_id}/isoforms", response_model=List[bio_schemas.IsoformResponse])
def get_gene_isoforms(gene_id: str, db: Session = Depends(get_db)):
    """List all isoforms for a specific gene. Accepts int DB id or HGNC symbol."""
    if gene_id.isdigit():
        gene_pk = int(gene_id)
    else:
        gene = db.query(bio_models.Gene).filter(bio_models.Gene.symbol == gene_id).first()
        if not gene:
            raise HTTPException(status_code=404, detail="Gene not found")
        gene_pk = gene.id
    isoforms = db.query(bio_models.Isoform).options(joinedload(bio_models.Isoform.gene)).filter(bio_models.Isoform.gene_id == gene_pk).all()
    if not isoforms:
        raise HTTPException(status_code=404, detail="No isoforms found for this gene")
    
    for isoform in isoforms:
        if isoform.is_canonical is None:
            if isoform.gene and isoform.gene.canonical_sequence and isoform.sequence:
                isoform.is_canonical = isoform.sequence == isoform.gene.canonical_sequence
            else:
                isoform.is_canonical = False
            
    return isoforms

@router.get("/genes/{gene_id}/string_network")
def get_gene_string_network(gene_id: str, db: Session = Depends(get_db)):
    """Proxy + cache STRING's own network SVG for a gene, served same-origin so
    it passes the site CSP (which blocks external images). Centred on the gene,
    STRING draws its high-confidence functional neighbourhood."""
    gene = db.query(bio_models.Gene).filter(bio_models.Gene.symbol == gene_id).first()
    if not gene:
        raise HTTPException(status_code=404, detail="Gene not found")

    identifier = gene.symbol
    safe_key = re.sub(r"[^A-Za-z0-9_.-]", "_", identifier)
    cache_path = os.path.join(_STRING_NETWORK_CACHE_DIR, f"{safe_key}.svg")
    headers = {"Cache-Control": "public, max-age=86400"}

    def _fresh_cache() -> Optional[bytes]:
        if os.path.exists(cache_path) and (time.time() - os.path.getmtime(cache_path)) < _STRING_NETWORK_TTL_S:
            with open(cache_path, "rb") as fh:
                return fh.read()
        return None

    cached = _fresh_cache()
    if cached is not None:
        return Response(content=cached, media_type="image/svg+xml", headers=headers)

    try:
        resp = httpx.get(
            "https://string-db.org/api/svg/network",
            params={
                "identifiers": identifier,
                "species": 9606,
                "network_flavor": "confidence",
                "required_score": 700,
                "caller_identity": "splisoforms.org",
            },
            timeout=20.0,
        )
        resp.raise_for_status()
        svg = resp.content
    except Exception:
        # Fall back to any stale cache before giving up.
        if os.path.exists(cache_path):
            with open(cache_path, "rb") as fh:
                return Response(content=fh.read(), media_type="image/svg+xml", headers=headers)
        raise HTTPException(status_code=502, detail="STRING network image unavailable")

    try:
        os.makedirs(_STRING_NETWORK_CACHE_DIR, exist_ok=True)
        with open(cache_path, "wb") as fh:
            fh.write(svg)
    except OSError:
        pass  # caching is best-effort; still serve the image

    return Response(content=svg, media_type="image/svg+xml", headers=headers)


@router.get("/genes/{gene_id}/string_interactions")
def get_gene_string_interactions(gene_id: str, db: Session = Depends(get_db)):
    """Gene-level STRING partners, deduped across the gene's isoforms. STRING
    evidence is gene-level (identical for every isoform), so we collapse to one
    row per partner keeping the best combined score."""
    if gene_id.isdigit():
        gene = db.query(bio_models.Gene).filter(bio_models.Gene.id == int(gene_id)).first()
    else:
        gene = db.query(bio_models.Gene).filter(bio_models.Gene.symbol == gene_id).first()
    if not gene:
        raise HTTPException(status_code=404, detail="Gene not found")

    rows = (
        db.query(bio_models.StringInteraction)
        .join(bio_models.Isoform, bio_models.Isoform.id == bio_models.StringInteraction.isoform_id)
        .filter(bio_models.Isoform.gene_id == gene.id)
        .all()
    )
    best: dict = {}
    for r in rows:
        name = r.partner_gene_symbol
        if not name:
            continue
        score = r.string_score or 0
        if name not in best or score > best[name]:
            best[name] = score
    partners = sorted(
        ({"partner_gene_symbol": k, "string_score": v} for k, v in best.items()),
        key=lambda p: p["string_score"],
        reverse=True,
    )
    return {"gene_symbol": gene.symbol, "partners": partners}


@router.patch("/isoform/{isoform_id}/flag", response_model=dict)
def toggle_flag(isoform_id: int, db: Session = Depends(get_db)):
    """Toggle the is_flagged field on an isoform. Returns {id, is_flagged}."""
    isoform = db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")
    isoform.is_flagged = not bool(isoform.is_flagged)
    db.commit()
    return {"id": isoform_id, "is_flagged": isoform.is_flagged}

@router.get("/isoforms/flagged", response_model=List[bio_schemas.IsoformDetailResponse])
def get_flagged_isoforms(db: Session = Depends(get_db)):
    """Return all isoforms marked as of interest by curators."""
    return db.query(bio_models.Isoform).options(
        joinedload(bio_models.Isoform.gene),
        joinedload(bio_models.Isoform.structures),
        joinedload(bio_models.Isoform.metadata_record),
    ).filter(bio_models.Isoform.is_flagged == True).all()


from fastapi import APIRouter, Depends, HTTPException, Query, BackgroundTasks

@router.get("/isoform/{isoform_id}", response_model=bio_schemas.IsoformDetailResponse)
def get_isoform(isoform_id: int, background_tasks: BackgroundTasks, db: Session = Depends(get_db)):
    """Detailed view including sequence, delta_length, URLs for cif/json structures, and dynamically cached PDB entries."""
    isoform = db.query(bio_models.Isoform).options(
        joinedload(bio_models.Isoform.metadata_record),
        joinedload(bio_models.Isoform.gene),
        joinedload(bio_models.Isoform.source_dataset),
        selectinload(bio_models.Isoform.source_memberships).joinedload(
            bio_models.IsoformSourceMembership.source
        ),
        selectinload(bio_models.Isoform.pfam_domains),
        selectinload(bio_models.Isoform.ptm_accessibility),
        joinedload(bio_models.Isoform.nmd_prediction),
        selectinload(bio_models.Isoform.structural_impacts),
        selectinload(bio_models.Isoform.neoantigen_predictions),
        selectinload(bio_models.Isoform.expression_data),
        joinedload(bio_models.Isoform.homodimer_analysis),
        selectinload(bio_models.Isoform.heterodimer_analyses),
        selectinload(bio_models.Isoform.string_interactions),
        selectinload(bio_models.Isoform.ted_domains),
        selectinload(bio_models.Isoform.junction_events),
        joinedload(bio_models.Isoform.disorder_score),
        joinedload(bio_models.Isoform.structural_metrics),
    ).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")
        
    # On-demand calculation if missing and sequence is small enough (to avoid long waits)
    if not isoform.disorder_score and isoform.sequence and len(isoform.sequence) < 1500:
        def _run_disorder():
            try:
                # We need a new session since the current one will close after response returns
                from app.database import SessionLocal
                from app.services.disorder_service import DisorderService
                bg_db = SessionLocal()
                try:
                    service = DisorderService()
                    service.analyze_isoform(bg_db, isoform_id)
                finally:
                    bg_db.close()
            except Exception as e:
                print(f"Lazy disorder analysis failed: {e}")
        background_tasks.add_task(_run_disorder)

    # Dynamically hydrate RCSB PDB targets if Uniprot ID is present but pdb_ids is empty
    meta = isoform.metadata_record
    if meta and meta.uniprot_id and not meta.uniprot_id.startswith("PB.") and not meta.pdb_ids:
        def _fetch_rcsb():
            try:
                from app.database import SessionLocal
                bg_db = SessionLocal()
                try:
                    local_meta = bg_db.query(bio_models.TranscriptMetadata).filter(bio_models.TranscriptMetadata.id == meta.id).first()
                    if not local_meta: return
                    search_query = {
                        "query": {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                                "operator": "exact_match",
                                "value": local_meta.uniprot_id
                            }
                        },
                        "return_type": "entry",
                        "request_options": {
                            "results_content_type": ["experimental"],
                            "sort": [{"sort_by": "score", "direction": "desc"}],
                            "scoring_strategy": "combined"
                        }
                    }
                    with httpx.Client(timeout=2.0) as client:
                        res = client.post("https://search.rcsb.org/rcsbsearch/v2/query", json=search_query)
                        if res.status_code == 200:
                            result_set = res.json().get("result_set", [])
                            if result_set:
                                local_meta.pdb_ids = [r["identifier"] for r in result_set[:5]]
                                bg_db.commit()
                finally:
                    bg_db.close()
            except Exception as e:
                print(f"Failed to fetch PDBs for {meta.uniprot_id}: {e}")
        background_tasks.add_task(_fetch_rcsb)

    # Resolve the canonical isoform for this gene (single query, used for multiple purposes below)
    canonical_iso = None
    if isoform.gene:
        canonical_iso = db.query(bio_models.Isoform).filter(
            bio_models.Isoform.gene_id == isoform.gene_id,
            bio_models.Isoform.is_canonical == True
        ).first()

    # Override gene.canonical_sequence with the canonical isoform's sequence (in-memory only).
    # gene.canonical_sequence is populated from UniProt at ingestion time and can drift from
    # the actual AFDB canonical isoform stored in the DB.  Using the isoform's sequence as the
    # reference ensures correct identity scores and a correct is_canonical flag.
    if canonical_iso and canonical_iso.sequence and isoform.gene:
        isoform.gene.canonical_sequence = canonical_iso.sequence

    # Hydrate gene.canonical_length from the actual canonical isoform length (preferred)
    # or fall back to transcript_metadata if unavailable.
    if isoform.gene:
        if canonical_iso and canonical_iso.sequence:
            setattr(isoform.gene, "canonical_length", len(canonical_iso.sequence))
        else:
            meta_can = db.query(bio_models.TranscriptMetadata).join(bio_models.Isoform).filter(
                bio_models.Isoform.gene_id == isoform.gene_id,
                bio_models.TranscriptMetadata.canonical_protein_length.isnot(None)
            ).first()
            if meta_can and meta_can.canonical_protein_length:
                setattr(isoform.gene, "canonical_length", int(meta_can.canonical_protein_length))

    # Determine whether this isoform is the canonical one
    if isoform.gene and isoform.gene.canonical_sequence and isoform.sequence:
        setattr(isoform, "is_canonical", isoform.sequence == isoform.gene.canonical_sequence)
    else:
        setattr(isoform, "is_canonical", isoform.is_canonical or False)

    # Hydrate into a non-ORM attribute (avoids SQLAlchemy intercepting the assignment)
    setattr(isoform, "dataset_presence", [
        bio_schemas.IsoformSourceMembershipResponse(
            source_id=m.source_id,
            source_name=m.source.name if m.source else "Unknown",
            in_normal=m.in_normal,
            in_primary_tumor=m.in_primary_tumor,
            in_metastasis=m.in_metastasis,
        )
        for m in isoform.source_memberships
    ])

    # Load canonical isoform's Pfam domains for reference (only for non-canonical isoforms)
    canonical_pfam_domains = []
    if isoform.gene and not getattr(isoform, 'is_canonical', False):
        canonical_iso_id = canonical_iso.id if canonical_iso else isoform.gene.canonical_transcript_id
        if canonical_iso_id:
            canonical_pfam_domains = db.query(bio_models.PfamDomain).filter(
                bio_models.PfamDomain.isoform_id == canonical_iso_id,
                bio_models.PfamDomain.domain_status == 'present'
            ).order_by(bio_models.PfamDomain.seq_start).all()
    setattr(isoform, 'canonical_pfam_domains', canonical_pfam_domains)

    return isoform

@router.get("/isoform/{isoform_id}/pfam", response_model=List[bio_schemas.PfamDomainResponse])
def get_pfam_domains(isoform_id: int, db: Session = Depends(get_db)):
    """Get Pfam domain annotations for a specific isoform."""
    isoform = db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")

    domains = db.query(bio_models.PfamDomain).filter(
        bio_models.PfamDomain.isoform_id == isoform_id
    ).order_by(bio_models.PfamDomain.seq_start).all()

    return domains

@router.post("/isoform/{isoform_id}/pfam/analyze", response_model=List[bio_schemas.PfamDomainResponse])
def trigger_pfam_analysis(
    isoform_id: int,
    force: bool = Query(False, description="Re-analyze even if results exist"),
    db: Session = Depends(get_db),
):
    """Trigger on-demand Pfam domain analysis for an isoform."""
    isoform = db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")

    try:
        from app.analysis.pfam_domains import analyze_isoform
        results = analyze_isoform(db, isoform, force=force)
        return results
    except FileNotFoundError as e:
        raise HTTPException(status_code=503, detail=f"Pfam database not available: {str(e)}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Pfam analysis failed: {str(e)}")

@router.get("/isoform/{isoform_id}/ptm", response_model=List[bio_schemas.PtmAccessibilityResponse])
def get_ptm_accessibility(isoform_id: int, db: Session = Depends(get_db)):
    """Get PTM site accessibility analysis for a specific isoform."""
    isoform = db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")

    sites = db.query(bio_models.PtmAccessibility).filter(
        bio_models.PtmAccessibility.isoform_id == isoform_id
    ).order_by(bio_models.PtmAccessibility.residue_position).all()

    return sites

@router.post("/isoform/{isoform_id}/ptm/analyze", response_model=List[bio_schemas.PtmAccessibilityResponse])
def trigger_ptm_analysis(
    isoform_id: int,
    force: bool = Query(False, description="Re-analyze even if results exist"),
    db: Session = Depends(get_db),
):
    """Trigger on-demand PTM site accessibility analysis for an isoform."""
    isoform = db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")

    try:
        from app.analysis.ptm_accessibility import analyze_isoform
        results = analyze_isoform(db, isoform, force=force)
        return results
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"PTM analysis failed: {str(e)}")

@router.get("/isoform/{isoform_id}/nmd", response_model=bio_schemas.NmdPredictionResponse)
def get_nmd_prediction(isoform_id: int, db: Session = Depends(get_db)):
    """Get NMD prediction for a specific isoform."""
    isoform = db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")

    prediction = db.query(bio_models.NmdPrediction).filter(
        bio_models.NmdPrediction.isoform_id == isoform_id
    ).first()

    if not prediction:
        raise HTTPException(status_code=404, detail="No NMD prediction found for this isoform")

    return prediction

@router.post("/isoform/{isoform_id}/nmd/analyze", response_model=bio_schemas.NmdPredictionResponse)
def trigger_nmd_analysis(
    isoform_id: int,
    force: bool = Query(False, description="Re-analyze even if results exist"),
    db: Session = Depends(get_db),
):
    """Trigger on-demand NMD prediction analysis for an isoform."""
    isoform = db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")

    try:
        from app.analysis.nmd_prediction import analyze_isoform
        result = analyze_isoform(db, isoform, force=force)
        if not result:
            raise HTTPException(status_code=422, detail="Isoform has no sequence, cannot predict NMD")
        return result
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"NMD analysis failed: {str(e)}")

@router.get("/isoform/{isoform_id}/structural-impact", response_model=List[bio_schemas.ExonStructuralImpactResponse])
def get_structural_impact(isoform_id: int, db: Session = Depends(get_db)):
    """Get per-residue structural impact (DSSP, phi/psi, pLDDT) for an isoform."""
    isoform = db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")

    impacts = db.query(bio_models.ExonStructuralImpact).filter(
        bio_models.ExonStructuralImpact.isoform_id == isoform_id
    ).order_by(bio_models.ExonStructuralImpact.residue_number).all()

    return impacts

@router.post("/isoform/{isoform_id}/structural-impact/analyze", response_model=List[bio_schemas.ExonStructuralImpactResponse])
def trigger_structural_impact_analysis(
    isoform_id: int,
    force: bool = Query(False, description="Re-analyze even if results exist"),
    db: Session = Depends(get_db),
):
    """Trigger on-demand per-residue structural impact analysis for an isoform."""
    isoform = db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")

    try:
        from app.analysis.structural_impact import analyze_isoform
        results = analyze_isoform(db, isoform, force=force)
        if not results:
            raise HTTPException(status_code=422, detail="No structure found or no residues extracted")
        return results
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Structural impact analysis failed: {str(e)}")

@router.get("/isoform/{isoform_id}/junctions", response_model=List[bio_schemas.IsoformJunctionResponse])
def get_isoform_junctions(isoform_id: int, db: Session = Depends(get_db)):
    """Get splice junction events for a specific isoform."""
    isoform = db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")

    junctions = db.query(bio_models.IsoformJunction).filter(
        bio_models.IsoformJunction.isoform_id == isoform_id
    ).order_by(bio_models.IsoformJunction.transcript_pos).all()

    return junctions

@router.get("/gene/{gene_id}/cell-type", response_model=List[bio_schemas.CellTypeExpressionResponse])
def get_cell_type_expression(gene_id: int, db: Session = Depends(get_db)):
    """Get cell-type expression data from HPA for a gene."""
    gene = db.query(bio_models.Gene).filter(bio_models.Gene.id == gene_id).first()
    if not gene:
        raise HTTPException(status_code=404, detail="Gene not found")

    results = db.query(bio_models.CellTypeExpression).filter(
        bio_models.CellTypeExpression.gene_id == gene_id
    ).order_by(bio_models.CellTypeExpression.expression_ncpm.desc()).all()

    return results

@router.post("/gene/{gene_id}/cell-type/analyze", response_model=List[bio_schemas.CellTypeExpressionResponse])
def trigger_cell_type_analysis(
    gene_id: int,
    force: bool = Query(False, description="Re-fetch even if results exist"),
    db: Session = Depends(get_db),
):
    """Trigger on-demand HPA cell-type expression lookup for a gene."""
    gene = db.query(bio_models.Gene).filter(bio_models.Gene.id == gene_id).first()
    if not gene:
        raise HTTPException(status_code=404, detail="Gene not found")

    try:
        from app.analysis.cell_type_specificity import analyze_gene
        results = analyze_gene(db, gene, force=force)
        return results
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Cell-type analysis failed: {str(e)}")

@router.get("/neoantigens/alleles")
def get_neoantigen_alleles(db: Session = Depends(get_db)):
    """Return all distinct HLA alleles present in the neoantigen predictions, ordered by count."""
    from app.models.bio import NeoantigenPrediction
    from sqlalchemy import func
    rows = (
        db.query(NeoantigenPrediction.hla_allele, func.count(NeoantigenPrediction.id).label("n"))
        .group_by(NeoantigenPrediction.hla_allele)
        .order_by(NeoantigenPrediction.hla_allele)
        .all()
    )
    return [{"allele": r.hla_allele, "count": r.n} for r in rows]


@router.get("/neoantigens/top")
def get_top_neoantigens(
    limit: int = Query(50, le=200),
    offset: int = 0,
    binding_category: Optional[str] = None,   # "strong_binder" | "weak_binder"
    hla_allele: Optional[str] = None,          # e.g. "HLA-A*02:01"
    gene: Optional[str] = None,                # gene symbol substring
    junction_type: Optional[str] = None,       # "novel_insertion" | "novel_substitution"
    priority_tier: Optional[str] = None,       # ranking v2: "P1" | "P2" | "P3" | "P4"
    flag: Optional[str] = None,                # decorating: immunogenic|splice_enabled|expressed|ptm_disrupting|high_fidelity
    min_plddt: Optional[float] = None,
    sort_by: str = "priority",    # priority | consensus | immunogenicity | affinity | plddt | gene | isoform | peptide | hla | junction
    sort_order: str = "asc",      # asc | desc  (asc = best-first for priority/affinity)
    db: Session = Depends(get_db),
):
    """
    Return globally ranked neoantigen predictions across all isoforms (ranking v2).

    Ranking is driven by immunological evidence: `priority_tier` (P1–P4) is the
    conjunction of multi-tool PRESENTATION consensus and BigMHC IMMUNOGENICITY.
    Splice origin / PTM disruption / expression / fold quality are DECORATING flags
    (filter with `flag`), not ranking gates. See app/analysis/neoantigen_tiers.py.

      priority_tier: P1 (strong presentation + immunogenic) · P2 (strong presentation)
                     · P3 (moderate presentation + immunogenic) · P4 (candidate)
      flag:          immunogenic · splice_enabled · expressed · ptm_disrupting · high_fidelity
                     · proteome_unique (peptide in no GENCODE reference protein)
                     · all_tools (NetMHCpan scored the row — 3-tool consensus)

    Deduplication: one row per (peptide_sequence, isoform_id), keeping the best-tier
    allele hit. Rows not presented by any tool (priority_tier IS NULL) are excluded.
    """
    # ── Build WHERE clauses ───────────────────────────────────────────────────
    # is_representative is the pre-computed dedup (one best row per peptide+isoform,
    # set by backfill_priority_tier.py) and implies priority_tier IS NOT NULL. Filtering
    # it replaces a DISTINCT ON over ~7.6M rows that ran 3× per request (page/count/stats).
    where_parts = ["n.is_representative"]
    params: dict = {}

    if binding_category:
        where_parts.append("n.binding_category = :binding_category")
        params["binding_category"] = binding_category

    if hla_allele:
        where_parts.append("n.hla_allele = :hla_allele")
        params["hla_allele"] = hla_allele

    if junction_type:
        where_parts.append("n.junction_type = :junction_type")
        params["junction_type"] = junction_type

    if min_plddt is not None:
        where_parts.append("n.avg_plddt >= :min_plddt")
        params["min_plddt"] = min_plddt

    if priority_tier in ("P1", "P2", "P3", "P4"):
        where_parts.append("n.priority_tier = :priority_tier")
        params["priority_tier"] = priority_tier

    _flag_sql = {
        "immunogenic":   "n.is_immunogenic = TRUE",
        "splice_enabled": "n.ptm_impact_status LIKE 'Splice-Enabled%'",
        "expressed":     "n.is_expressed = TRUE",
        "ptm_disrupting": "n.ptm_canonical_residue_lost = TRUE",
        "high_fidelity": "(n.avg_plddt > 70 AND n.min_plddt > 50)",
        # peptide absent from EVERY GENCODE reference protein (not just the canonical)
        "proteome_unique": "n.is_proteome_unique = TRUE",
        "all_tools":      "n.tools_agreeing IS NOT NULL AND n.netmhcpan_rank_el IS NOT NULL",
    }
    if flag in _flag_sql:
        where_parts.append(_flag_sql[flag])

    if gene:
        where_parts.append("(g.symbol ILIKE :gene OR g.symbol = :gene_exact)")
        params["gene"] = f"%{gene}%"
        params["gene_exact"] = gene.upper()

    where_sql = " AND ".join(where_parts)

    # priority_tier rank for sorting (lower = better; NULL handled by WHERE). Kept in the
    # page SELECT for display; ORDER BY uses n.priority_tier directly ('P1'<'P2'<… sorts
    # the same) so the partial index ix_neo_representative(priority_tier, consensus_rank)
    # WHERE is_representative (NeoantigenPrediction.__table_args__) serves the default
    # sort as an index scan instead of sorting the whole set.
    priority_rank_expr = """
        CASE n.priority_tier
            WHEN 'P1' THEN 1 WHEN 'P2' THEN 2 WHEN 'P3' THEN 3 WHEN 'P4' THEN 4 ELSE 5
        END
    """

    dir_ = "DESC" if sort_order.lower() == "desc" else "ASC"
    _col_map = {
        "priority":       f"n.priority_tier {dir_}, n.consensus_rank ASC NULLS LAST",
        "consensus":      f"n.consensus_rank {dir_} NULLS LAST, n.priority_tier ASC",
        "immunogenicity": f"n.bigmhc_im_score {dir_} NULLS LAST, n.priority_tier ASC",
        "affinity":       f"n.binding_affinity_nm {dir_}, n.priority_tier ASC",
        "plddt":          f"n.avg_plddt {dir_} NULLS LAST, n.priority_tier ASC",
        "gene":           f"g.symbol {dir_}, n.priority_tier ASC",
        "isoform":        f"i.transcript_id {dir_}, n.priority_tier ASC",
        "peptide":        f"n.peptide_sequence {dir_}",
        "hla":            f"n.hla_allele {dir_}",
        "junction":       f"n.junction_type {dir_}, n.priority_tier ASC",
    }
    outer_order = _col_map.get(sort_by, "n.priority_tier ASC, n.consensus_rank ASC NULLS LAST")

    # count/stats only need the genes join when a gene filter is present (it's the only
    # clause that references g.*). Stats always joins isoforms for the unique-gene count.
    needs_gene_join = gene is not None
    gene_join = ("JOIN isoforms i ON i.id = n.isoform_id "
                 "JOIN genes g ON g.id = i.gene_id") if needs_gene_join else ""

    # ── Paginated results — one lean query, index-ordered, LIMIT-bounded ───────
    paginated_sql = f"""
        SELECT n.id, n.isoform_id, n.peptide_sequence, n.peptide_start, n.peptide_end,
               n.hla_allele, n.binding_affinity_nm AS affinity_nm, n.binding_category,
               n.percentile_rank, n.junction_type, n.presentation_score, n.processing_score,
               n.ptm_impact_status, n.avg_plddt, n.min_plddt, n.is_ptm_sensitive,
               n.ptm_canonical_residue_lost, n.netmhcpan_rank_el, n.netmhcpan_rank_ba,
               n.bigmhc_el_score, n.bigmhc_im_score, n.consensus_rank, n.is_ba_rescue,
               n.presentation_tier, n.priority_tier, n.is_immunogenic, n.is_expressed,
               n.proteome_uniqueness, n.is_proteome_unique, n.tools_agreeing,
               (n.percentile_rank IS NOT NULL)::int + (n.netmhcpan_rank_el IS NOT NULL)::int
                 + (n.bigmhc_el_score IS NOT NULL)::int AS tools_available,
               {priority_rank_expr}             AS priority_rank,
               i.transcript_id, i.transcript_name,
               g.id AS gene_id, g.symbol AS gene_symbol
          FROM neoantigen_predictions n
          JOIN isoforms i ON i.id = n.isoform_id
          JOIN genes    g ON g.id = i.gene_id
         WHERE {where_sql}
         ORDER BY {outer_order}
         LIMIT :limit OFFSET :offset
    """
    params["limit"] = limit
    params["offset"] = offset
    rows = db.execute(text(paginated_sql), params).mappings().all()

    count_params = {k: v for k, v in params.items() if k not in ("limit", "offset")}

    # ── Total count + summary stats — per filter set, cached for a few minutes ──
    cache_key = (where_sql, tuple(sorted(count_params.items())))
    cached = _atlas_stats_get(cache_key)
    if cached is not None:
        total_count, stats = cached
    else:
        total_count, stats = _atlas_count_and_stats(db, where_sql, gene_join, count_params)
        _atlas_stats_put(cache_key, total_count, stats)

    return {
        "items": [dict(r) for r in rows],
        "total_count": total_count,
        "stats": stats,
    }


def _atlas_count_and_stats(db: Session, where_sql: str, gene_join: str, count_params: dict):
    """total_count + the stats block of /neoantigens/top for one filter set.

    Three statements, each shaped to stay on the partial indexes over the
    representative subset (declared on NeoantigenPrediction.__table_args__):

    * FILTER counts — one pass; every counted column is INCLUDEd in ix_neo_repr_stats,
      so this is an index-only scan whenever the filters are on included columns too
      (the default view, priority_tier, flag=proteome_unique/immunogenic/expressed/…).
      `total` doubles as the endpoint's total_count (same WHERE, same join).
    * unique_peptides — its OWN statement on purpose. Folded into the FILTER query
      (as it was before 2026-10) the planner must sort all ~2.6M rows by peptide to
      serve COUNT(DISTINCT) next to the other aggregates: an on-disk merge sort,
      ~3.5 s. Alone, it reads ix_neo_repr_stats presorted by peptide_sequence (~0.2 s).
    * unique_genes — distinct isoform_ids come from ix_neo_repr_isoform, then map to
      genes; far cheaper than joining ~2.6M rows to isoforms for a distinct count.
    """
    stats_sql = f"""
        SELECT
            COUNT(*)                                                           AS total,
            COUNT(*) FILTER (WHERE n.priority_tier = 'P1')                     AS p1,
            COUNT(*) FILTER (WHERE n.priority_tier = 'P2')                     AS p2,
            COUNT(*) FILTER (WHERE n.priority_tier = 'P3')                     AS p3,
            COUNT(*) FILTER (WHERE n.priority_tier = 'P4')                     AS p4,
            COUNT(*) FILTER (WHERE n.presentation_tier = 'STRONG')             AS strong_presentation,
            COUNT(*) FILTER (WHERE n.is_immunogenic = TRUE)                    AS immunogenic,
            COUNT(*) FILTER (WHERE n.is_expressed = TRUE)                      AS expressed,
            COUNT(*) FILTER (WHERE n.ptm_impact_status LIKE 'Splice-Enabled%') AS splice_enabled,
            COUNT(*) FILTER (WHERE n.is_ba_rescue = TRUE)                      AS ba_rescue,
            COUNT(*) FILTER (WHERE n.is_proteome_unique = TRUE)                AS proteome_unique,
            COUNT(*) FILTER (WHERE n.netmhcpan_rank_el IS NULL)                AS missing_netmhcpan
          FROM neoantigen_predictions n {gene_join}
         WHERE {where_sql}
    """
    stats_row = db.execute(text(stats_sql), count_params).mappings().first()
    stats = dict(stats_row) if stats_row else {}
    total_count = stats.get("total") or 0

    unique_peptides_sql = f"""
        SELECT COUNT(DISTINCT n.peptide_sequence)
          FROM neoantigen_predictions n {gene_join}
         WHERE {where_sql}
    """
    stats["unique_peptides"] = db.execute(text(unique_peptides_sql), count_params).scalar() or 0

    unique_genes_sql = f"""
        SELECT COUNT(DISTINCT gene_id) FROM isoforms
         WHERE id IN (SELECT n.isoform_id FROM neoantigen_predictions n {gene_join}
                       WHERE {where_sql})
    """
    stats["unique_genes"] = db.execute(text(unique_genes_sql), count_params).scalar() or 0
    return total_count, stats


@router.get("/isoform/{isoform_id}/neoantigen", response_model=List[bio_schemas.NeoantigenPredictionResponse])
def get_neoantigen_predictions(isoform_id: int, db: Session = Depends(get_db)):
    """Get neoantigen binding predictions for a specific isoform."""
    isoform = db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")

    from sqlalchemy import case as sa_case, nullslast
    predictions = (
        db.query(bio_models.NeoantigenPrediction)
        .filter(bio_models.NeoantigenPrediction.isoform_id == isoform_id)
        .order_by(
            # Prefer consensus_rank when available (lower = better); fall back to affinity
            nullslast(bio_models.NeoantigenPrediction.consensus_rank.asc()),
            bio_models.NeoantigenPrediction.binding_affinity_nm.asc(),
        )
        .all()
    )
    return predictions

@router.post("/isoform/{isoform_id}/neoantigen/analyze", response_model=List[bio_schemas.NeoantigenPredictionResponse])
def trigger_neoantigen_analysis(
    isoform_id: int,
    force: bool = Query(False, description="Re-analyze even if results exist"),
    db: Session = Depends(get_db),
):
    """Trigger on-demand neoantigen prediction for an isoform."""
    isoform = db.query(bio_models.Isoform).filter(bio_models.Isoform.id == isoform_id).first()
    if not isoform:
        raise HTTPException(status_code=404, detail="Isoform not found")

    try:
        from app.analysis.neoantigen import analyze_isoform
        results = analyze_isoform(db, isoform, force=force)
        return results
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Neoantigen analysis failed: {str(e)}")

@router.get("/sources", response_model=List[bio_schemas.SourceResponse])
def get_sources(db: Session = Depends(get_db)):
    """List all available data sources."""
    return db.query(bio_models.Source).order_by(bio_models.Source.name).all()

@router.get("/stats")
def get_db_stats(db: Session = Depends(get_db)):
    """Return database metrics for the landing page."""
    isoform_count = db.query(bio_models.Isoform).count()
    gene_count = db.query(bio_models.Gene).count()
    return {
        "isoform_count": isoform_count,
        "gene_count": gene_count
    }

@router.get("/external/uniprot/{uniprot_id}")
def get_uniprot_metadata(uniprot_id: str):
    """Proxy UniProt API to avoid CORS issues on the frontend."""
    try:
        with httpx.Client(timeout=5.0) as client:
            res = client.get(f"https://rest.uniprot.org/uniprotkb/search?query=accession:{uniprot_id}")
            res.raise_for_status()
            return res.json()
    except Exception as e:
        raise HTTPException(status_code=502, detail=f"Failed to fetch from UniProt: {str(e)}")

@router.get("/isoforms/browse", response_model=bio_schemas.BrowseResponse)
def browse_isoforms(
    limit: int = 24,
    offset: int = 0,
    sort_by: str = "impact_score",
    sort_order: str = "desc",
    min_plddt: Optional[float] = None,
    min_isoforms: Optional[int] = None,
    source_id: Optional[int] = None,
    impact_filter: Optional[str] = None,
    is_canonical: Optional[bool] = None,
    search: Optional[str] = None,
    db: Session = Depends(get_db)
):
    """Ultra-fast paginated global browse using the Discovery Index."""
    from sqlalchemy import func, desc, asc, case
    from app.models.bio import IsoformDiscoveryIndex

    # Subquery: count novel isoforms per gene (for min_isoforms filter and display badge)
    iso_count_subq = (
        db.query(
            bio_models.Isoform.gene_id,
            func.count(bio_models.Isoform.id).label("isoform_count")
        )
        .filter(bio_models.Isoform.is_canonical == False)
        .group_by(bio_models.Isoform.gene_id)
        .subquery()
    )

    # Core query: join Isoform with its Discovery Index
    q = db.query(
        bio_models.Isoform,
        bio_models.Gene.symbol.label("gene_symbol"),
        bio_models.Gene.uniprot_id.label("uniprot_id"),
        IsoformDiscoveryIndex.plddt_mean,
        bio_models.Source.name.label("source_name"),
        IsoformDiscoveryIndex.num_domains_lost,
        IsoformDiscoveryIndex.num_domains_gained,
        IsoformDiscoveryIndex.has_neoantigen,
        IsoformDiscoveryIndex.neoantigen_count,
        IsoformDiscoveryIndex.neoantigen_p1_count,
        IsoformDiscoveryIndex.impact_score,
        IsoformDiscoveryIndex.delta_length.label("idx_delta_length"),
        func.coalesce(iso_count_subq.c.isoform_count, 1).label("isoform_count"),
        IsoformDiscoveryIndex.nmd_status,
        IsoformDiscoveryIndex.has_ptm_exposure_change,
        IsoformDiscoveryIndex.has_ptm_buried_change,
        IsoformDiscoveryIndex.has_ptm_loss,
    ).join(
        IsoformDiscoveryIndex, bio_models.Isoform.id == IsoformDiscoveryIndex.isoform_id
    ).join(
        bio_models.Gene, bio_models.Isoform.gene_id == bio_models.Gene.id
    ).outerjoin(
        bio_models.Source, bio_models.Isoform.source_id == bio_models.Source.id
    ).outerjoin(
        iso_count_subq, bio_models.Gene.id == iso_count_subq.c.gene_id
    )

    # Text search filter
    if search:
        search_term = f"%{search}%"
        q = q.filter(
            bio_models.Gene.symbol.ilike(search_term) |
            bio_models.Isoform.transcript_id.ilike(search_term) |
            bio_models.Gene.uniprot_id.ilike(search_term)
        )

    # Filters (Using indexed columns)
    if min_plddt is not None:
        q = q.filter(IsoformDiscoveryIndex.plddt_mean >= min_plddt)
    if min_isoforms is not None and min_isoforms > 1:
        q = q.filter(func.coalesce(iso_count_subq.c.isoform_count, 1) >= min_isoforms)
    if source_id is not None:
        q = q.filter(bio_models.Isoform.source_id == source_id)
    if is_canonical is not None:
        q = q.filter(bio_models.Isoform.is_canonical == is_canonical)
    elif not search:
        # By default in browse mode, only show novel isoforms (non-canonical)
        q = q.filter(bio_models.Isoform.is_canonical == False)

    if impact_filter == "domain_loss":
        q = q.filter(IsoformDiscoveryIndex.num_domains_lost > 0)
    elif impact_filter == "neoantigen":
        q = q.filter(IsoformDiscoveryIndex.has_neoantigen == True)
    elif impact_filter == "ptm_impact":
        q = q.filter(
            or_(
                IsoformDiscoveryIndex.has_ptm_exposure_change == True,
                IsoformDiscoveryIndex.has_ptm_buried_change == True,
                IsoformDiscoveryIndex.has_ptm_loss == True
            )
        )
    elif impact_filter == "nmd_target":
        q = q.filter(IsoformDiscoveryIndex.nmd_status == "nmd_target")

    # Sorting — respects sort_order (asc/desc)
    order_fn = desc if sort_order != "asc" else asc
    if sort_by == "impact_score":
        q = q.order_by(order_fn(IsoformDiscoveryIndex.impact_score))
    elif sort_by == "plddt":
        q = q.order_by(order_fn(IsoformDiscoveryIndex.plddt_mean))
    elif sort_by == "isoform_count":
        q = q.order_by(order_fn(func.coalesce(iso_count_subq.c.isoform_count, 1)))
    elif sort_by == "neoantigen":
        q = q.order_by(order_fn(IsoformDiscoveryIndex.neoantigen_count))
    elif sort_by == "neoantigen_p1":
        q = q.order_by(order_fn(IsoformDiscoveryIndex.neoantigen_p1_count))
    elif sort_by == "nmd":
        nmd_rank = case(
            (IsoformDiscoveryIndex.nmd_status == "nmd_target", 2),
            (IsoformDiscoveryIndex.nmd_status == "nmd_escape", 1),
            else_=0,
        )
        q = q.order_by(order_fn(nmd_rank))
    elif sort_by == "ptm":
        ptm_rank = case(
            (IsoformDiscoveryIndex.has_ptm_loss == True, 3),
            (IsoformDiscoveryIndex.has_ptm_buried_change == True, 2),
            (IsoformDiscoveryIndex.has_ptm_exposure_change == True, 1),
            else_=0,
        )
        q = q.order_by(order_fn(ptm_rank))
    else:  # length / sequence delta — sort by signed value so desc = large insertions first
        q = q.order_by(order_fn(IsoformDiscoveryIndex.delta_length))

    # Fast Count and Fetch
    total_count = q.count()
    results = q.offset(offset).limit(limit).all()

    items = []
    for row in results:
        (iso, symbol, uniprot, plddt, source, lost, gained, neo, neo_count, neo_p1, score, idx_delta,
         iso_cnt, nmd_status, ptm_exp, ptm_bur, ptm_loss) = row
        items.append({
            "id": iso.id,
            "transcript_id": iso.transcript_id,
            "transcript_name": iso.transcript_name,
            "gene_id": iso.gene_id,
            "gene_symbol": symbol,
            "uniprot_id": uniprot,
            "plddt_mean": plddt,
            "delta_length": idx_delta if idx_delta is not None else iso.delta_length,
            "source_name": source,
            "has_neoantigen": neo,
            "neoantigen_count": neo_count or 0,
            "neoantigen_p1_count": neo_p1 or 0,
            "num_domains_lost": lost,
            "num_domains_gained": gained,
            "nmd_status": nmd_status,
            "has_ptm_exposure_change": ptm_exp,
            "has_ptm_buried_change": ptm_bur,
            "has_ptm_loss": ptm_loss,
            "isoform_count": iso_cnt,
            "impact_score": score
        })

    # Optimized Stats (Using Index)
    dist_stats = db.query(
        func.count(case((IsoformDiscoveryIndex.plddt_mean < 50, 1), else_=None)).label("p50"),
        func.count(case((IsoformDiscoveryIndex.plddt_mean.between(50, 70), 1), else_=None)).label("p70"),
        func.count(case((IsoformDiscoveryIndex.plddt_mean.between(70, 90), 1), else_=None)).label("p90"),
        func.count(case((IsoformDiscoveryIndex.plddt_mean >= 90, 1), else_=None)).label("p100"),
        func.count(case((IsoformDiscoveryIndex.num_domains_lost > 0, 1), else_=None)).label("lost"),
        func.count(case((IsoformDiscoveryIndex.num_domains_gained > 0, 1), else_=None)).label("gained"),
        func.count(case((IsoformDiscoveryIndex.has_neoantigen == True, 1), else_=None)).label("neo"),
        # PTM counts
        func.count(case((IsoformDiscoveryIndex.has_ptm_exposure_change == True, 1), else_=None)).label("ptm_exp"),
        func.count(case((IsoformDiscoveryIndex.has_ptm_buried_change == True, 1), else_=None)).label("ptm_bur"),
        func.count(case((IsoformDiscoveryIndex.has_ptm_loss == True, 1), else_=None)).label("ptm_loss"),
        # NMD counts — buckets match the values the NMD pipeline actually writes
        # (nmd_prediction.py): nmd_target / nmd_escape / predicted_translatable / unknown.
        # target & escape require the CDS 50-nt PTC analysis; until cds_sequence is
        # populated they are 0 and every isoform is translatable-or-unknown.
        func.count(case((IsoformDiscoveryIndex.nmd_status == "nmd_target", 1), else_=None)).label("nmd_target"),
        func.count(case((IsoformDiscoveryIndex.nmd_status == "nmd_escape", 1), else_=None)).label("nmd_escape"),
        func.count(case((IsoformDiscoveryIndex.nmd_status == "predicted_translatable", 1), else_=None)).label("nmd_translatable"),
        func.count(case((IsoformDiscoveryIndex.nmd_status == "unknown", 1), else_=None)).label("nmd_unknown"),
        # Additional aggregate stats
        func.avg(func.abs(IsoformDiscoveryIndex.delta_length)).label("avg_delta"),
        func.avg(
            case((IsoformDiscoveryIndex.num_domains_lost > 0, IsoformDiscoveryIndex.num_domains_lost), else_=None)
        ).label("avg_domains_lost"),
    ).first()

    stats = {
        "plddt_distribution": [
            {"bin": "0–50", "count": dist_stats.p50 or 0},
            {"bin": "50–70", "count": dist_stats.p70 or 0},
            {"bin": "70–90", "count": dist_stats.p90 or 0},
            {"bin": "90–100", "count": dist_stats.p100 or 0},
        ],
        "impact_categories": [
            {"category": "Domain Loss", "count": dist_stats.lost or 0},
            {"category": "Domain Gained", "count": dist_stats.gained or 0},
            {"category": "Neoantigens", "count": dist_stats.neo or 0},
        ],
        "ptm_distribution": [
            {"impact": "Exposed", "count": dist_stats.ptm_exp or 0},
            {"impact": "Buried", "count": dist_stats.ptm_bur or 0},
            {"impact": "Loss", "count": dist_stats.ptm_loss or 0},
        ],
        "nmd_distribution": [
            {"status": "Target", "count": dist_stats.nmd_target or 0},
            {"status": "Escape", "count": dist_stats.nmd_escape or 0},
            {"status": "Translatable", "count": dist_stats.nmd_translatable or 0},
            {"status": "Unknown", "count": dist_stats.nmd_unknown or 0},
        ],
        "source_distribution": [],
        # Preserve None when SQL AVG returns NULL (no data yet) so the
        # frontend "—" fallback renders rather than "0 aa" / "0.0".
        "avg_delta_length": round(float(dist_stats.avg_delta), 1) if dist_stats.avg_delta is not None else None,
        "avg_domains_lost": round(float(dist_stats.avg_domains_lost), 2) if dist_stats.avg_domains_lost is not None else None,
    }

    return {"items": items, "total_count": total_count, "stats": stats}


# ─── Drug–isoform layer ──────────────────────────────────────────────────────
def _resolve_gene(gene_id: str, db: Session):
    q = db.query(bio_models.Gene)
    g = (q.filter(bio_models.Gene.id == int(gene_id)).first()
         if gene_id.isdigit() else
         q.filter(bio_models.Gene.symbol == gene_id).first())
    if not g:
        raise HTTPException(status_code=404, detail="Gene not found")
    return g


@router.get("/drugs/stats", response_model=bio_schemas.DrugStatsResponse)
def drug_stats(db: Session = Depends(get_db)):
    """Proteome-wide drug-lens overview."""
    n_drugs = db.query(bio_models.Drug).count()
    n_complexes = db.query(bio_models.DrugComplex).count()
    n_genes = db.query(bio_models.DrugComplex.gene_id).distinct().count()
    fate_rows = (db.query(bio_models.DrugPocketFate.fate, func.count())
                 .group_by(bio_models.DrugPocketFate.fate).all())
    fate_counts = {f or "unknown": c for f, c in fate_rows}
    # gene-level subtractive tally (decidable rows only)
    sub_sq = (db.query(
                bio_models.DrugComplex.gene_id.label("gid"),
                func.bool_or(bio_models.DrugPocketFate.fate.in_(("lost", "disrupted"))).label("sub"))
              .join(bio_models.DrugComplex,
                    bio_models.DrugComplex.id == bio_models.DrugPocketFate.drug_complex_id)
              .filter(bio_models.DrugPocketFate.decidable.is_(True))
              .group_by(bio_models.DrugComplex.gene_id).subquery())
    n_total = db.query(sub_sq).count()
    n_sub = db.query(sub_sq).filter(sub_sq.c.sub.is_(True)).count()
    return bio_schemas.DrugStatsResponse(
        n_drugs=n_drugs, n_complexes=n_complexes, n_drugged_genes=n_genes,
        fate_counts=fate_counts, n_subtractive_genes=n_sub, n_total_scored_genes=n_total)


@router.get("/drugs/genes", response_model=List[bio_schemas.DrugGeneSummary])
def drug_genes(
    limit: int = Query(5000, le=10000),
    offset: int = 0,
    sort: str = Query("subtractive", pattern="^(subtractive|n_drugs)$"),
    db: Session = Depends(get_db),
):
    """Ranked table of drugged genes for the /explore/drugs lens.

    Two grouped queries (drug/complex counts per gene; per-(gene,isoform)
    subtractive/preserved rollup) merged in Python — no per-gene N+1.
    """
    dc, f, g = bio_models.DrugComplex, bio_models.DrugPocketFate, bio_models.Gene

    # (a) drug + complex counts per gene
    counts = dict()
    for gid, symbol, n_drugs, n_complexes in (
        db.query(dc.gene_id, g.symbol,
                 func.count(func.distinct(dc.drug_id)),
                 func.count(func.distinct(dc.id)))
          .join(g, g.id == dc.gene_id)
          .group_by(dc.gene_id, g.symbol).all()):
        counts[gid] = {"symbol": symbol, "n_drugs": n_drugs, "n_complexes": n_complexes,
                       "sub": 0, "pres": 0, "n_iso": 0}

    # (b) per-(gene, isoform) decidable rollup: an isoform is subtractive if ANY
    # decidable pocket is lost/disrupted, preserved if its calls are all preserved.
    iso = (db.query(
                dc.gene_id.label("gid"), f.isoform_id.label("iid"),
                func.bool_or(f.fate.in_(("lost", "disrupted"))).label("sub"),
                func.bool_or(f.fate == "preserved").label("pres"))
           .join(dc, dc.id == f.drug_complex_id)
           .filter(f.decidable.is_(True))
           .group_by(dc.gene_id, f.isoform_id).subquery())
    for gid, sub, pres in db.query(iso.c.gid, iso.c.sub, iso.c.pres).all():
        if gid not in counts:
            continue
        counts[gid]["n_iso"] += 1
        if sub:
            counts[gid]["sub"] += 1
        elif pres:
            counts[gid]["pres"] += 1

    results = []
    for gid, r in counts.items():
        denom = r["sub"] + r["pres"]
        results.append(bio_schemas.DrugGeneSummary(
            gene_id=gid, symbol=r["symbol"], n_drugs=r["n_drugs"],
            n_complexes=r["n_complexes"], n_isoforms_scored=r["n_iso"],
            n_subtractive=r["sub"], n_preserved=r["pres"],
            subtractive_frac=round(r["sub"] / denom, 3) if denom else 0.0))
    key = (lambda r: (r.subtractive_frac, r.n_drugs)) if sort == "subtractive" else (lambda r: r.n_drugs)
    results.sort(key=key, reverse=True)
    return results[offset:offset + limit]


def _cif_url(c) -> Optional[str]:
    """/static path to a drug complex's co-fold model.cif (served via the
    drug_complexes symlink), or None when the structure wasn't produced.

    LOWERCASED deliberately. `job_name` carries the UniProt accession in upper case
    (2-hydroxyflutamide--P10275--0e8126a2), but every directory AlphaFold 3 produced is
    lower case (…--p10275--…) — in our copy and in the source dataset alike. On macOS the
    filesystem is case-insensitive, so the mismatch resolved fine in local dev and went
    unnoticed; on the Linux server it is case-sensitive and every drug CIF 404s. Do not
    "simplify" this back to a bare job_name.
    """
    if not c.has_structure:
        return None
    slug = c.job_name.lower()
    return f"/static/data/structures/drug_complexes/{slug}/{slug}_model.cif"


@router.get("/gene/{gene_id}/drugs", response_model=bio_schemas.GeneDrugsResponse)
def gene_drugs(gene_id: str, db: Session = Depends(get_db)):
    """All drug complexes + per-isoform pocket fates for one gene."""
    g = _resolve_gene(gene_id, db)
    dc, f = bio_models.DrugComplex, bio_models.DrugPocketFate
    complexes = (db.query(dc).options(joinedload(dc.drug))
                 .filter(dc.gene_id == g.id).all())
    for c in complexes:
        setattr(c, "cif_url", _cif_url(c))
    cid_to_drug = {c.id: c.drug for c in complexes}
    fates = (db.query(f).join(dc, dc.id == f.drug_complex_id)
             .filter(dc.gene_id == g.id).all())
    for fate in fates:  # hydrate drug identity for the matrix
        d = cid_to_drug.get(fate.drug_complex_id)
        if d:
            setattr(fate, "drug_name", d.name)
            setattr(fate, "pubchem_cid", d.pubchem_cid)
    return bio_schemas.GeneDrugsResponse(
        gene_id=g.id, symbol=g.symbol,
        complexes=complexes, fates=fates)


@router.get("/isoform/{isoform_id}/drug-fate",
            response_model=List[bio_schemas.DrugPocketFateResponse])
def isoform_drug_fate(isoform_id: int, db: Session = Depends(get_db)):
    """Per-isoform drug-pocket fates (what splicing does to each drug's pocket)."""
    dc, f = bio_models.DrugComplex, bio_models.DrugPocketFate
    fates = (db.query(f).options(joinedload(f.complex).joinedload(dc.drug))
             .filter(f.isoform_id == isoform_id).all())
    for fate in fates:
        d = fate.complex.drug if fate.complex else None
        if d:
            setattr(fate, "drug_name", d.name)
            setattr(fate, "pubchem_cid", d.pubchem_cid)
    return fates

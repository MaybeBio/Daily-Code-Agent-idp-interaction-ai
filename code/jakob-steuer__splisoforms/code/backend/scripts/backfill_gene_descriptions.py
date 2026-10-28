#!/usr/bin/env python3
"""
backfill_gene_descriptions.py
=============================
Populate `genes.protein_name` and `genes.function_summary` from UniProt REST.

For every gene with a `uniprot_id` and no `protein_name` (or `--force`), this
script GETs `https://rest.uniprot.org/uniprotkb/{uniprot_id}.json` and stores:

  protein_name        proteinDescription.recommendedName.fullName.value
  function_summary    comments[type=FUNCTION].texts[0].value
                      (truncated to ~1.2 kB to keep gene rows lightweight)

Replaces the auto-generated junk Gene.description strings (e.g.
"Auto-generated for PacBio isoform PB.124576.52") for display purposes.

Serial by design: UniProt REST is rate-limited per IP; the script paces at
~5 req/s with retry on 429 / 5xx.

Usage
-----
    cd backend
    python scripts/backfill_gene_descriptions.py --dry-run
    python scripts/backfill_gene_descriptions.py
    python scripts/backfill_gene_descriptions.py --limit 50 --force
"""

import argparse
import json
import logging
import os
import sys
import time
import urllib.error
import urllib.request
from typing import Optional

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sqlalchemy as sa

from _backfill_io import Progress, print_summary, setup_logging
from app.database import SessionLocal
from app.models.bio import Gene

log = logging.getLogger("backfill_gene_descriptions")

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"
USER_AGENT = "splisoforms/1.0"
REQ_DELAY = 0.2                       # 5 req/s — well under UniProt's published limit
FUNCTION_SUMMARY_CHAR_LIMIT = 1200    # truncate over-long FUNCTION text


def fetch_uniprot(acc: str, retries: int = 3) -> Optional[dict]:
    url = f"{UNIPROT_API}/{acc}.json"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
            )
            with urllib.request.urlopen(req, timeout=15) as resp:
                return json.loads(resp.read())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                log.warning("UniProt 404 for %s (likely retired accession)", acc)
                return None
            if e.code in (429, 500, 502, 503) and attempt < retries - 1:
                wait = 2 ** attempt
                log.warning("UniProt %d for %s; retry in %ds", e.code, acc, wait)
                time.sleep(wait)
                continue
            log.error("UniProt HTTPError %s for %s", e.code, acc)
            return None
        except Exception as exc:
            log.error("UniProt fetch failed for %s: %s", acc, exc)
            return None
    return None


def parse_protein_name(data: dict) -> Optional[str]:
    """Best-effort extraction of the canonical recommended protein name."""
    pd = data.get("proteinDescription") or {}
    rec = pd.get("recommendedName") or {}
    full = (rec.get("fullName") or {}).get("value")
    if full:
        return full
    # Fallbacks: submittedName (TrEMBL) or first alternativeName
    sub = pd.get("submissionNames") or []
    if sub:
        sn = (sub[0].get("fullName") or {}).get("value")
        if sn:
            return sn
    alt = pd.get("alternativeNames") or []
    if alt:
        an = (alt[0].get("fullName") or {}).get("value")
        if an:
            return an
    return None


def parse_function_summary(data: dict) -> Optional[str]:
    """First text from the FUNCTION comment, falling back to SIMILARITY.

    Many well-characterised-enough proteins carry no CC-FUNCTION comment but do carry a
    SIMILARITY one ("Belongs to the X family."). Extracting FUNCTION alone left 2,460 genes
    with an empty description even though UniProt had usable text — e.g. C1orf174 (Q8IYL3),
    CCNJ (Q5T5M9), AADACL4 (Q5VUY2). The isoform Overview panel has always rendered either
    type, so accept both here, preferring FUNCTION.

    Truncated to a sane display length.
    """
    by_type: dict[str, str] = {}
    for comment in data.get("comments") or []:
        ctype = comment.get("commentType")
        if ctype not in ("FUNCTION", "SIMILARITY") or ctype in by_type:
            continue
        texts = comment.get("texts") or []
        if not texts:
            continue
        value = texts[0].get("value")
        if value:
            by_type[ctype] = value

    value = by_type.get("FUNCTION") or by_type.get("SIMILARITY")
    if not value:
        return None
    if len(value) > FUNCTION_SUMMARY_CHAR_LIMIT:
        value = value[:FUNCTION_SUMMARY_CHAR_LIMIT].rsplit(". ", 1)[0] + " …"
    return value


def main():
    parser = argparse.ArgumentParser(
        description="Backfill UniProt protein_name and function_summary onto genes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--dry-run", action="store_true",
                        help="Fetch a sample of N and print what would be stored; no DB writes.")
    parser.add_argument("--force",   action="store_true",
                        help="Re-fetch even if protein_name is already populated.")
    parser.add_argument("--limit",   type=int, default=None,
                        help="Process at most N genes.")
    args = parser.parse_args()

    logfile = setup_logging("backfill_gene_descriptions")

    db = SessionLocal()
    try:
        q = db.query(Gene).filter(Gene.uniprot_id.isnot(None))
        if not args.force:
            # Fetch a gene if EITHER field is missing. Filtering on protein_name alone
            # silently stranded 2,457 genes that had a name but no function_summary (a
            # UniProt entry can carry recommendedName without a CC-FUNCTION comment on the
            # first pass) — their gene pages then render an empty description.
            q = q.filter(
                sa.or_(Gene.protein_name.is_(None), Gene.function_summary.is_(None))
            )
        if args.limit:
            q = q.limit(args.limit)
        genes = q.all()

        total = len(genes)
        scope: list[str] = []
        if args.force: scope.append("force=True")
        if args.limit: scope.append(f"limit={args.limit}")
        print("[backfill_gene_descriptions] UniProt protein_name + function_summary backfill")
        print(f"  scope: {total:,} genes with UniProt"
              + (f"  |  {' | '.join(scope)}" if scope else "")
              + ("  |  DRY RUN" if args.dry_run else ""))

        if total == 0:
            print("Nothing to process.")
            return

        ok = no_data = errors = 0
        err_list: list[tuple[str, str]] = []

        with Progress(total, desc="UniProt") as pbar:
            for gene in genes:
                acc = gene.uniprot_id
                data = fetch_uniprot(acc)
                if data is None:
                    no_data += 1
                else:
                    pname = parse_protein_name(data)
                    func  = parse_function_summary(data)
                    if not pname and not func:
                        no_data += 1
                        log.info("gene_id=%d %s — no protein_name / function in UniProt response",
                                 gene.id, acc)
                    else:
                        if not args.dry_run:
                            gene.protein_name     = pname
                            gene.function_summary = func
                        log.info("gene_id=%d %s name=%r func_len=%s",
                                 gene.id, acc, (pname or "")[:60],
                                 len(func) if func else 0)
                        ok += 1
                pbar.update(1, status=f"ok={ok} empty={no_data} err={errors}")
                time.sleep(REQ_DELAY)

                # Periodic flush so we don't lose work on a crash
                if not args.dry_run and ok and ok % 200 == 0:
                    db.commit()

        if not args.dry_run:
            db.commit()

        print_summary(
            "GENE DESCRIPTIONS BACKFILL COMPLETE",
            counts=[
                ("Genes processed",           total),
                ("✓ With protein_name/func",  ok),
                ("– No UniProt data",         no_data),
                ("✗ Errors",                  errors),
            ],
            total=total,
            logfile=logfile,
            errors=err_list,
        )
        if args.dry_run:
            print("  DRY RUN — no DB writes performed.")

    except Exception:
        db.rollback()
        raise
    finally:
        db.close()


if __name__ == "__main__":
    main()

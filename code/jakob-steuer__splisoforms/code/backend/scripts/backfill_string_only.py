#!/usr/bin/env python3
"""
backfill_string_only.py
=======================
Phase 6f STRING interactions backfill (homodimer already done in 4.8).

Why a separate script?
----------------------
`run_ppi_analysis.py` runs both homodimer + STRING per isoform; with
126 k isoforms × per-API latency that came out at ~2.8/s = 12 h. Since
isoforms of the same gene share their canonical's UniProt accession,
every per-isoform STRING call was redundant. This script:

  1. Collects the SET of unique UniProt accessions referenced by any
     gene that has at least one isoform (~25 k accessions, ~5× fewer).
  2. Hits the STRING REST API in parallel (ThreadPoolExecutor=8) for
     each accession ONCE, caches partners in memory.
  3. Writes StringInteraction rows for every isoform of every gene that
     resolved to an accession with hits.

Runtime: ~30-60 min for ~25 k accessions at STRING's ~1-2 req/s effective
rate per thread.
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from app.database import SessionLocal
from app.models.bio import Gene, Isoform, StringInteraction

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("backfill_string_only")

# interaction_partners (not network): returns the query protein's DIRECT partners,
# with the query always in preferredName_A and the partner in preferredName_B. This
# is alias-proof — STRING's preferredName for our UniProt accession is often an old
# gene symbol (e.g. Q96FB5/METTL25B resolves to "RRNAD1"), which broke the previous
# `network` + gene_symbol self-matching (it silently dropped or mislabelled partners).
STRING_API_URL = "https://string-db.org/api/json/interaction_partners"
DEFAULT_WORKERS = 8
# 400 = STRING "medium confidence" (its default web view). 700 (high) hid genes whose
# partners are all sub-0.7 (e.g. METTL25B), so the site showed "no interactions" for
# proteins that STRING clearly connects. The score is stored per row for UI filtering.
REQUIRED_SCORE = 400


def fetch_partners(uniprot_id: str) -> tuple[str, list[dict] | None, str | None]:
    """One STRING request. Returns (acc, partners_list, error)."""
    try:
        r = requests.get(
            STRING_API_URL,
            params={
                "identifiers": uniprot_id,
                "species": 9606,
                "required_score": REQUIRED_SCORE,
            },
            timeout=30,
        )
        r.raise_for_status()
        return (uniprot_id, r.json(), None)
    except Exception as e:
        return (uniprot_id, None, f"{type(e).__name__}: {e}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                        help=f"Concurrent STRING requests (default {DEFAULT_WORKERS})")
    parser.add_argument("--limit", type=int,
                        help="Cap accessions processed (testing)")
    parser.add_argument("--force", action="store_true",
                        help="Re-fetch even where StringInteraction rows exist")
    args = parser.parse_args()

    db = SessionLocal()

    # Collect distinct (uniprot_id, gene_id) of genes whose canonical exists
    log.info("collecting accessions to fetch ...")
    rows = (db.query(Gene.uniprot_id, Gene.id, Gene.symbol)
              .filter(Gene.uniprot_id.isnot(None))
              .all())
    acc_to_genes: dict[str, list[tuple[int, str]]] = {}
    for acc, gid, sym in rows:
        acc_to_genes.setdefault(acc, []).append((gid, sym))
    accs = sorted(acc_to_genes.keys())
    log.info("%d distinct accessions, %d genes total", len(accs), len(rows))

    # Drop accessions already done (any StringInteraction for any of their isoforms)
    if not args.force:
        existing = set(
            row[0] for row in
            (db.query(Isoform.gene_id)
              .join(StringInteraction, StringInteraction.isoform_id == Isoform.id)
              .distinct().all())
        )
        before = len(accs)
        accs = [a for a in accs
                if any(gid not in existing for gid, _sym in acc_to_genes[a])]
        log.info("skipping %d accessions already done; %d remaining",
                 before - len(accs), len(accs))

    if args.limit:
        accs = accs[:args.limit]

    if not accs:
        log.info("nothing to fetch.")
        db.close()
        return 0

    # --- Fetch in parallel ---
    log.info("fetching STRING partners with %d workers ...", args.workers)
    results: dict[str, list[dict]] = {}
    errors = 0
    t0 = time.time()
    n = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(fetch_partners, a) for a in accs]
        for fut in as_completed(futures):
            acc, partners, err = fut.result()
            n += 1
            if err:
                errors += 1
                if errors <= 5:
                    log.error("  %s: %s", acc, err)
            else:
                results[acc] = partners or []
            if n % 500 == 0 or n == len(accs):
                el = time.time() - t0
                rate = n / el if el > 0 else 0
                eta = (len(accs) - n) / rate if rate > 0 else 0
                log.info("  [%d/%d] hits=%d err=%d (%.1f/s · eta %.0fs)",
                         n, len(accs), sum(1 for v in results.values() if v),
                         errors, rate, eta)

    log.info("fetch phase done in %.0fs. inserting StringInteraction rows ...",
             time.time() - t0)

    # --- Apply rows ---
    n_rows = 0
    n_genes_with_partners = 0
    for acc, partners in results.items():
        if not partners:
            continue
        n_genes_with_partners += 1
        gene_entries = acc_to_genes.get(acc, [])
        for gene_id, gene_symbol in gene_entries:
            isos = db.query(Isoform).filter(Isoform.gene_id == gene_id).all()
            if not isos:
                continue
            # Clear existing for these isoforms (--force or first run)
            iso_ids = [i.id for i in isos]
            if args.force:
                (db.query(StringInteraction)
                   .filter(StringInteraction.isoform_id.in_(iso_ids))
                   .delete(synchronize_session=False))
            for iso in isos:
                seen = set()
                for p in partners:
                    # interaction_partners always puts the query in _A and the
                    # partner in _B, so the partner is preferredName_B directly —
                    # no self-name matching (which failed on alias mismatches).
                    partner = p.get("preferredName_B", "") or p.get("preferredName_A", "")
                    score = p.get("score", 0)
                    if not partner or partner == gene_symbol:
                        continue
                    if partner in seen:
                        continue
                    seen.add(partner)
                    # STRING's JSON API returns score as a 0-1 float; string_score is an
                    # INTEGER column and the UI treats it as STRING's 0-1000 combined score
                    # (bar width = score/10). Store on the 0-1000 scale, not the raw float
                    # (which rounded every score to 0/1 — the whole table was ~all 1s).
                    db.add(StringInteraction(
                        isoform_id=iso.id,
                        partner_gene_symbol=partner,
                        string_score=int(round(float(score) * 1000)) if score else 0,
                    ))
                    n_rows += 1
        if n_genes_with_partners % 500 == 0:
            db.commit()
            log.info("  ... committed (%d genes, %d rows)",
                     n_genes_with_partners, n_rows)
    db.commit()

    log.info("✅ summary:")
    log.info("  accessions queried:           %d", len(accs))
    log.info("  accessions with partners:     %d", n_genes_with_partners)
    log.info("  StringInteraction rows added: %d", n_rows)
    log.info("  errors:                       %d", errors)
    db.close()
    return 0 if errors == 0 else 2


if __name__ == "__main__":
    sys.exit(main())

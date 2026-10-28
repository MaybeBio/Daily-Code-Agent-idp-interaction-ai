"""
Decorate neoantigen_predictions with canonical-PTM impact — the step
load_neoantigen_results.py defers when run with --skip-ptm.

For every predicted neoepitope this fills the *decorating* fields
(`ptm_impact_status`, `is_ptm_sensitive`, `ptm_canonical_residue_lost`) that drive
the SPL/PTM chips on the neoantigen browser.  It does NOT touch priority_tier or any
presentation/immunogenicity field — those are set at load and are independent.

Why a separate script: the canonical PTM sites come from live UniProt/EBI
(`get_ptm_sites`, ~2 calls per accession).  Doing that inline in the loader serialised
~25k network calls into the critical path.  Here it is:
  • concurrent  — a thread pool fetches accessions in parallel;
  • cached      — results persist to data/ptm_cache/canonical_ptms.json, so re-runs and
                  the prod load are instant (ship the cache);
  • resumable   — cached accessions are skipped; interrupt and restart freely.

Logic mirrors load_neoantigen_results.py (the direct/flanking-window classification),
so a --skip-ptm load followed by this script is equivalent to a full load.

    python scripts/backfill_neoantigen_ptm.py                 # full run
    python scripts/backfill_neoantigen_ptm.py --limit-accessions 20   # smoke test
"""
import os
import sys
import json
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy import text
from psycopg2.extras import execute_values

from app.database import SessionLocal, engine
from app.analysis.ptm_data import get_ptm_sites

CACHE_PATH_DEFAULT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "ptm_cache", "canonical_ptms.json",
)


def _load_cache(path):
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    return {}


def _save_cache(path, cache):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(cache, fh)
    os.replace(tmp, path)


def _classify(start, end, canonical_ptms, isoform_seq):
    """Return (ptm_impact_status, is_sensitive, residue_lost). Mirrors the loader."""
    direct = [p for p in canonical_ptms if start <= p[0] <= end]
    flanking = [p for p in canonical_ptms
                if (start - 3 <= p[0] < start) or (end < p[0] <= end + 3)]

    if direct:
        any_lost = False
        for pos, restype in direct:
            if 1 <= pos <= len(isoform_seq):
                if isoform_seq[pos - 1] != restype:
                    any_lost = True
            else:
                any_lost = True
        status = ("Direct PTM-Neoantigen (residue lost)" if any_lost
                  else "Direct PTM-Neoantigen (residue retained)")
        return status, True, any_lost

    if flanking:
        found = lost = False
        for pos, restype in flanking:
            if 1 <= pos <= len(isoform_seq) and isoform_seq[pos - 1] == restype:
                found = True
            else:
                lost = True
        if found:
            return "Processing-Blocked", True, False
        if lost:
            return "Splice-Enabled", True, False

    return "Standard", False, False


def main():
    ap = argparse.ArgumentParser(description="Decorate neoantigen_predictions with canonical-PTM impact")
    ap.add_argument("--workers", type=int, default=16, help="Concurrent UniProt fetchers")
    ap.add_argument("--cache", default=CACHE_PATH_DEFAULT, help="On-disk PTM-site cache (JSON)")
    ap.add_argument("--refresh-cache", action="store_true", help="Ignore cached accessions and refetch")
    ap.add_argument("--limit-accessions", type=int, help="Smoke test: only this many accessions")
    args = ap.parse_args()

    db = SessionLocal()

    # ── 1. gene→uniprot, and the isoforms (with a uniprot gene) that have predictions ──
    print("Collecting accessions to fetch …", flush=True)
    gene_uniprot = dict(db.execute(text(
        "SELECT id, uniprot_id FROM genes WHERE uniprot_id IS NOT NULL")).all())

    iso_rows = db.execute(text("""
        SELECT DISTINCT i.id, i.gene_id, i.sequence
          FROM neoantigen_predictions n
          JOIN isoforms i ON i.id = n.isoform_id
         WHERE i.gene_id IN (SELECT id FROM genes WHERE uniprot_id IS NOT NULL)
    """)).all()
    iso_seq = {r[0]: (r[2] or "") for r in iso_rows}
    iso_uniprot = {r[0]: gene_uniprot.get(r[1]) for r in iso_rows}
    accessions = sorted({u for u in iso_uniprot.values() if u})
    print(f"  {len(iso_rows):,} predicted isoforms map to {len(accessions):,} unique accessions")

    # ── 2. fetch canonical PTM sites (concurrent + cached) ──
    cache = {} if args.refresh_cache else _load_cache(args.cache)
    todo = [a for a in accessions if a not in cache]
    if args.limit_accessions:
        todo = todo[:args.limit_accessions]
    print(f"  cached: {len(accessions) - len(todo):,}  |  to fetch: {len(todo):,}", flush=True)

    def fetch(acc):
        try:
            sites = get_ptm_sites(acc)
            return acc, [[s.residue_position, s.residue_type] for s in sites]
        except Exception as e:
            return acc, {"__error__": str(e)}

    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(fetch, a): a for a in todo}
        for fut in as_completed(futures):
            acc, sites = fut.result()
            if isinstance(sites, dict):          # error → don't poison cache, retry next run
                continue
            cache[acc] = sites
            done += 1
            if done % 250 == 0:
                _save_cache(args.cache, cache)
                print(f"    fetched {done:,}/{len(todo):,}", flush=True)
    _save_cache(args.cache, cache)
    print(f"  fetched {done:,}; cache now holds {len(cache):,} accessions", flush=True)

    # ── 3. classify each prediction row and collect the non-Standard updates ──
    # Only isoforms whose accession actually carries PTM sites can change status.
    ptm_isoforms = [iid for iid, u in iso_uniprot.items() if cache.get(u)]
    print(f"Classifying predictions for {len(ptm_isoforms):,} PTM-bearing isoforms …", flush=True)

    updates = []            # (id, status, is_sensitive, residue_lost)
    B = 2000
    for i in range(0, len(ptm_isoforms), B):
        batch = ptm_isoforms[i:i + B]
        rows = db.execute(text("""
            SELECT id, isoform_id, peptide_start, peptide_end
              FROM neoantigen_predictions
             WHERE isoform_id = ANY(:ids)
        """), {"ids": batch}).all()
        for pid, iso_id, start, end in rows:
            sites = [tuple(s) for s in cache.get(iso_uniprot[iso_id], [])]
            status, sens, lost = _classify(start, end, sites, iso_seq.get(iso_id, ""))
            if status != "Standard":
                updates.append((pid, status, sens, lost))

    print(f"  {len(updates):,} rows change to a non-Standard PTM status", flush=True)

    # ── 4. bulk write (raw psycopg2 for speed) ──
    if updates:
        raw = engine.raw_connection()
        try:
            cur = raw.cursor()
            execute_values(cur, """
                UPDATE neoantigen_predictions AS n SET
                    ptm_impact_status          = v.status,
                    is_ptm_sensitive           = v.sens,
                    ptm_canonical_residue_lost = v.lost
                FROM (VALUES %s) AS v(id, status, sens, lost)
                WHERE n.id = v.id
            """, updates, template="(%s, %s, %s, %s)", page_size=10000)
            raw.commit()
        finally:
            raw.close()

    dist = db.execute(text("""
        SELECT ptm_impact_status, count(*) FROM neoantigen_predictions
        GROUP BY ptm_impact_status ORDER BY count(*) DESC
    """)).all()
    print("ptm_impact_status distribution:")
    for st, n in dist:
        print(f"  {st:42s} {n:,}")
    db.close()
    print("✅ PTM decoration complete.")


if __name__ == "__main__":
    main()

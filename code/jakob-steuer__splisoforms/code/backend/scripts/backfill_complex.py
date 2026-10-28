"""
Backfill complex analysis using the correct AFDB summary endpoint.

The correct endpoint is:
    https://alphafold.ebi.ac.uk/api/uniprot/summary/{uniprot_id}.json
which lists ALL structures (monomer + NVIDIA Multimer homodimers/heterodimers) with an
explicit oligomeric_state field ("HOMODIMER" / "HETERODIMER" / "MONOMER").

Usage:
    python backend/scripts/backfill_complex.py [--gene MAOB] [--reset-flags] [--force] [--limit N]

Options:
    --gene        Process a single gene symbol only (for testing)
    --reset-flags Clear afdb_complex_checked for genes before running
    --force       Re-analyze even if a complex record already exists
    --limit       Cap number of genes processed (useful for dry-runs)
"""
import os
import sys
import time
import argparse

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy import text
from app.database import SessionLocal
from app.models import bio as bio_models
from app.services.ppi_service import PPIService


def reset_stale_flags(db, gene_symbol=None):
    """
    Clear afdb_complex_checked on genes.
    """
    q = """
        UPDATE genes
           SET afdb_complex_checked = FALSE
         WHERE afdb_complex_checked = TRUE
           AND uniprot_id IS NOT NULL
    """
    if gene_symbol:
        q += f" AND symbol = '{gene_symbol}'"
    result = db.execute(text(q))
    db.commit()
    print(f"[reset] Cleared afdb_complex_checked on {result.rowcount} genes.")


def get_canonical_isoforms(db, gene_symbol=None, force=False, limit=None):
    """
    Return (isoform_id, gene_symbol, uniprot_id) for canonical isoforms of genes
    that have a uniprot_id and still need complex analysis.
    """
    gene_filter = f"AND g.symbol = '{gene_symbol}'" if gene_symbol else ""

    if force:
        q = f"""
            SELECT i.id, g.symbol, g.uniprot_id
              FROM isoforms i
              JOIN genes g ON g.id = i.gene_id
             WHERE i.is_canonical = TRUE
               AND g.uniprot_id IS NOT NULL
               {gene_filter}
             ORDER BY g.symbol
             {"LIMIT " + str(limit) if limit else ""}
        """
    else:
        q = f"""
            SELECT i.id, g.symbol, g.uniprot_id
              FROM isoforms i
              JOIN genes g ON g.id = i.gene_id
             WHERE i.is_canonical = TRUE
               AND g.uniprot_id IS NOT NULL
               {gene_filter}
               AND g.afdb_complex_checked = FALSE
             ORDER BY g.symbol
             {"LIMIT " + str(limit) if limit else ""}
        """
    rows = db.execute(text(q)).fetchall()
    return [(r[0], r[1], r[2]) for r in rows]


def main():
    parser = argparse.ArgumentParser(description="Backfill AFDB complex analysis")
    parser.add_argument("--gene", type=str, help="Single gene symbol (e.g. MAOB)")
    parser.add_argument("--reset-flags", action="store_true", default=True,
                        help="Reset stale afdb_complex_checked flags before running (default: on)")
    parser.add_argument("--no-reset", dest="reset_flags", action="store_false",
                        help="Skip flag reset")
    parser.add_argument("--force", action="store_true",
                        help="Re-analyze genes that already have a complex analysis record")
    parser.add_argument("--limit", type=int, default=None,
                        help="Process at most N genes (for testing)")
    args = parser.parse_args()

    db = SessionLocal()
    ppi = PPIService(db)

    # Step 1: reset stale flags so the service will re-check them
    if args.reset_flags:
        reset_stale_flags(db, gene_symbol=args.gene)

    # Step 2: collect canonical isoforms that need analysis
    targets = get_canonical_isoforms(db, gene_symbol=args.gene, force=args.force, limit=args.limit)
    total = len(targets)
    print(f"[*] {total} canonical isoforms queued for complex analysis.")

    found = 0
    not_found = 0
    errors = 0
    t0 = time.time()

    for idx, (iso_id, sym, uniprot_id) in enumerate(targets, 1):
        try:
            isoform = db.query(bio_models.Isoform).get(iso_id)
            if not isoform:
                continue

            isoform.gene.afdb_complex_checked = False

            complex_summaries = ppi._find_complex_entries(uniprot_id)

            if complex_summaries:
                print(f"  [{idx}/{total}] {sym} ({uniprot_id}): {len(complex_summaries)} complexes found")
                
                # If force mode, we just let it overwrite since the service uses db.delete
                ppi._run_complex_analysis(isoform, uniprot_id)
                found += 1
            else:
                isoform.gene.afdb_complex_checked = True
                db.commit()
                not_found += 1
                if idx % 50 == 0 or total <= 20:
                    print(f"  [{idx}/{total}] {sym} ({uniprot_id}): no complex in AFDB")

        except Exception as e:
            db.rollback()
            errors += 1
            print(f"  [{idx}/{total}] ERROR {sym} ({uniprot_id}): {e}")

    elapsed = time.time() - t0
    print(f"\n[+] Done in {elapsed:.1f}s")
    print(f"    Complexes found : {found}")
    print(f"    No complexes    : {not_found}")
    print(f"    Errors          : {errors}")

    db.close()


if __name__ == "__main__":
    main()

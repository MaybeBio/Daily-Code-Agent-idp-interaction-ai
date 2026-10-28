#!/usr/bin/env python3
"""
backfill_neoantigen_uniqueness.py — proteome-wide novelty of every candidate peptide.

A neoepitope in SPLISOFORMS is novel RELATIVE TO ITS GENE'S CANONICAL by
construction. That does not make it tumour-specific: the same 8–11-mer can occur
in another reference isoform of the same gene, or in a reference protein of a
different gene, in which case it is a self-peptide the immune system is tolerant
to. The preprint names a proteome-wide uniqueness filter as planned; the 2026-10
review (F-4) measured 22 % same-gene reference matches among sampled P1 peptides.

Method
------
1. Collect the distinct candidate peptides from `neoantigen_predictions`
   (~1.9 M, lengths 8–11) into a set.
2. Scan every GENCODE v47 reference protein (source_id = 1) once: for each
   position and each k in 8..11, if the k-mer is a candidate, record the gene.
   This is O(proteome × 4) set lookups (~150 M) and needs no k-mer index.
3. Classify each (peptide, isoform gene) pair:
     unique      — the peptide occurs in no reference protein
     self_gene   — it occurs only in reference isoform(s) of the same gene
     other_gene  — it occurs in a reference protein of another gene
   The canonical isoform itself never contains the peptide (novelty filter), so a
   `self_gene` hit is another annotated isoform of the gene.
4. Set-based UPDATE via a temp table (peptide, gene_id → class).

Notes
-----
* The reference used is the set of reference PROTEINS in the DB (GENCODE v47
  protein-coding, md5-deduplicated). Non-coding-derived ORFs and the long-read
  isoforms are NOT part of the background — "novel vs the normal proteome".
* Peptides carried into the DB from a later HPC batch get NULL until this script
  is re-run; it is idempotent and fast (minutes), so run it after every load.

Usage:
    python backend/scripts/backfill_neoantigen_uniqueness.py [--dry-run]
"""
import argparse
import os
import sys
import time
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from sqlalchemy import text

from _backfill_io import exit_status
from app.database import SessionLocal

KMER_LENGTHS = (8, 9, 10, 11)


def main() -> int:
    ap = argparse.ArgumentParser(description="Proteome-wide uniqueness of neoepitope peptides")
    ap.add_argument("--dry-run", action="store_true", help="classify and report, no writes")
    args = ap.parse_args()
    db = SessionLocal()
    t0 = time.time()

    # 1. candidate peptides with the gene of each isoform they come from
    print("Loading candidate peptides …", flush=True)
    rows = db.execute(text("""
        SELECT DISTINCT n.peptide_sequence, i.gene_id
          FROM neoantigen_predictions n
          JOIN isoforms i ON i.id = n.isoform_id
    """)).fetchall()
    pep_genes: dict[str, set[int]] = defaultdict(set)
    for pep, gid in rows:
        pep_genes[pep].add(gid)
    wanted = set(pep_genes)
    print(f"  {len(rows):,} (peptide, gene) pairs; {len(wanted):,} distinct peptides "
          f"({time.time() - t0:.0f}s)", flush=True)
    if not wanted:
        print("Nothing to do.")
        return 0

    # 2. scan the reference proteome once
    print("Scanning GENCODE reference proteins …", flush=True)
    ref = db.execute(text("""
        SELECT i.gene_id, i.sequence FROM isoforms i
         WHERE i.source_id = (SELECT id FROM sources WHERE name LIKE 'GENCODE%' ORDER BY id LIMIT 1)
           AND i.sequence IS NOT NULL
    """)).fetchall()
    hits: dict[str, set[int]] = defaultdict(set)
    n_lookups = 0
    for gid, seq in ref:
        seq = seq.upper()
        L = len(seq)
        for k in KMER_LENGTHS:
            for i in range(0, L - k + 1):
                kmer = seq[i:i + k]
                if kmer in wanted:
                    hits[kmer].add(gid)
            n_lookups += max(0, L - k + 1)
    print(f"  {len(ref):,} reference proteins, {n_lookups:,} k-mer lookups, "
          f"{len(hits):,} candidate peptides found in the reference proteome "
          f"({time.time() - t0:.0f}s)", flush=True)

    # 3. classify per (peptide, gene)
    cls_rows: list[tuple[str, int, str]] = []
    counts = defaultdict(int)
    for pep, genes in pep_genes.items():
        found = hits.get(pep)
        for gid in genes:
            if not found:
                c = "unique"
            elif found <= {gid}:
                c = "self_gene"
            else:
                c = "other_gene"
            cls_rows.append((pep, gid, c))
            counts[c] += 1
    print("  classification over (peptide, gene) pairs: " +
          ", ".join(f"{k}={v:,}" for k, v in sorted(counts.items())), flush=True)

    if args.dry_run:
        print("(dry-run — no writes)")
        return 0

    # 4. set-based update through a temp table
    print("Writing …", flush=True)
    db.execute(text("DROP TABLE IF EXISTS _neo_uniq"))
    db.execute(text("CREATE TEMP TABLE _neo_uniq (peptide varchar(30), gene_id int, cls varchar(20))"))
    conn = db.connection().connection          # raw psycopg2 for COPY
    cur = conn.cursor()
    import io
    buf = io.StringIO()
    for pep, gid, c in cls_rows:
        buf.write(f"{pep}\t{gid}\t{c}\n")
    buf.seek(0)
    cur.copy_from(buf, "_neo_uniq", columns=("peptide", "gene_id", "cls"))
    db.execute(text("CREATE INDEX ON _neo_uniq (peptide, gene_id)"))
    res = db.execute(text("""
        UPDATE neoantigen_predictions n
           SET proteome_uniqueness = u.cls,
               is_proteome_unique  = (u.cls = 'unique')
          FROM _neo_uniq u, isoforms i
         WHERE i.id = n.isoform_id AND u.peptide = n.peptide_sequence AND u.gene_id = i.gene_id
    """))
    db.commit()
    print(f"  updated {res.rowcount:,} prediction rows ({time.time() - t0:.0f}s)")
    dist = db.execute(text(
        "SELECT proteome_uniqueness, count(*), count(DISTINCT peptide_sequence) "
        "FROM neoantigen_predictions GROUP BY 1 ORDER BY 1")).fetchall()
    for c, n, p in dist:
        print(f"    {c or 'NULL':12s} rows={n:>10,}  peptides={p:>9,}")
    p1 = db.execute(text(
        "SELECT proteome_uniqueness, count(DISTINCT peptide_sequence) FROM neoantigen_predictions "
        "WHERE priority_tier='P1' GROUP BY 1 ORDER BY 1")).fetchall()
    print("  P1 peptides by class: " + ", ".join(f"{c}={n:,}" for c, n in p1))
    db.close()
    return exit_status(len(cls_rows), 0 if res.rowcount else len(cls_rows))


if __name__ == "__main__":
    sys.exit(main())

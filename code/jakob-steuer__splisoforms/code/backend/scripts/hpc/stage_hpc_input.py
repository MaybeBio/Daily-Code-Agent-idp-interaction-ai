#!/usr/bin/env python3
"""
Export isoform + canonical sequences from the local DB for HPC neoantigen prediction.

Writes a TSV with one row per candidate isoform:
  isoform_id, transcript_id, gene_id, sequence, canonical_sequence

Skips:
  - Isoforms that are their gene's canonical (self-comparison → no predictions).
  - Isoforms with no amino-acid sequence.
  - Isoforms whose gene has no canonical sequence (nothing to diff against).

Usage (from project root):
    python backend/scripts/hpc/stage_hpc_input.py -o hpc_input.tsv
"""
import argparse
import csv
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

from sqlalchemy import func
from app.database import SessionLocal
from app.models.bio import Gene, Isoform


def get_canonical_sequences(db) -> dict[int, tuple[str, int | None]]:
    """Build gene_id → (canonical_sequence, canonical_isoform_id) map.

    Priority:
      1. The sequence of the gene's CANONICAL ISOFORM (Gene.canonical_transcript_id)
      2. Gene.canonical_sequence (UniProt reviewed SwissProt)
      3. Longest AFDB isoform sequence
      4. Longest isoform sequence overall

    The canonical isoform comes FIRST deliberately. It is the MANE-first canonical that
    the rest of the DB (and the UI) treats as authoritative, and it is what a neoepitope
    is defined against: a peptide is only "novel" if it is ABSENT from the canonical
    protein. Diffing against anything else (e.g. a UniProt sequence that corresponds to no
    isoform, or a "longest isoform" fallback) produces peptides that ARE present in the
    canonical — i.e. not neoepitopes at all. That mismatch is exactly the bug that put
    111,295 bogus peptides in the DB; keep this aligned with the validation query:

      strpos(canonical_isoform.sequence, peptide) > 0   =>   NOT a neoepitope
    """
    genes = db.query(Gene).all()
    result: dict[int, tuple[str, int | None]] = {}

    gene_ids = [g.id for g in genes]

    isoforms_by_gene: dict[int, list] = {}
    rows = (
        db.query(
            Isoform.id,
            Isoform.gene_id,
            Isoform.transcript_id,
            Isoform.sequence,
            func.length(Isoform.sequence).label("seq_len"),
        )
        .filter(Isoform.gene_id.in_(gene_ids), Isoform.sequence.isnot(None))
        .all()
    )
    for r in rows:
        isoforms_by_gene.setdefault(r.gene_id, []).append(r)

    seq_by_iso_id = {r.id: r.sequence for r in rows}

    for gene in genes:
        gid = gene.id
        canonical_iso_id = gene.canonical_transcript_id

        # Canonical sequence (for novel-region diffing) — canonical isoform wins.
        canon_seq = seq_by_iso_id.get(canonical_iso_id) if canonical_iso_id else None

        if not canon_seq:
            canon_seq = gene.canonical_sequence

        if not canon_seq:
            isos = isoforms_by_gene.get(gid, [])
            afdb = [i for i in isos if i.transcript_id and i.transcript_id.startswith("AFDB_")]
            if afdb:
                afdb.sort(key=lambda x: x.seq_len or 0, reverse=True)
                canon_seq = afdb[0].sequence
            elif isos:
                isos.sort(key=lambda x: x.seq_len or 0, reverse=True)
                canon_seq = isos[0].sequence
            else:
                continue

        if canon_seq:
            result[gid] = (canon_seq, canonical_iso_id)

    return result


def main():
    parser = argparse.ArgumentParser(description="Stage isoform data for HPC neoantigen prediction")
    parser.add_argument("-o", "--output", required=True, help="Output TSV path")
    parser.add_argument("--limit", type=int, help="Limit number of isoforms (for testing)")
    parser.add_argument(
        "--ids",
        help="Path to a file of isoform_ids (one per line). Stage ONLY these isoforms. "
             "Used for targeted re-runs (e.g. isoforms whose canonical basis changed, or "
             "that were never submitted) instead of re-predicting the whole DB.",
    )
    args = parser.parse_args()

    db = SessionLocal()

    print("Loading canonical sequences …")
    canon_map = get_canonical_sequences(db)
    print(f"  {len(canon_map)} genes with canonical sequences")

    wanted: set[int] | None = None
    if args.ids:
        with open(args.ids) as fh:
            wanted = {int(line) for line in fh if line.strip()}
        print(f"Restricting to {len(wanted):,} isoform_ids from {args.ids}")

    print("Querying candidate isoforms …")
    query = (
        db.query(Isoform)
        .filter(Isoform.sequence.isnot(None), func.length(Isoform.sequence) > 0)
    )
    if wanted:
        query = query.filter(Isoform.id.in_(wanted))
    if args.limit:
        query = query.limit(args.limit)
    isoforms = query.all()

    if wanted:
        missing = wanted - {i.id for i in isoforms}
        if missing:
            print(f"  WARNING: {len(missing):,} requested ids have no usable sequence / do not exist")

    written = 0
    skipped_no_canonical = 0
    skipped_is_canonical = 0

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["isoform_id", "transcript_id", "gene_id", "sequence", "canonical_sequence"])

        for iso in isoforms:
            entry = canon_map.get(iso.gene_id)
            if not entry:
                skipped_no_canonical += 1
                continue

            canon_seq, canon_iso_id = entry

            if canon_iso_id and canon_iso_id == iso.id:
                skipped_is_canonical += 1
                continue

            if iso.sequence == canon_seq:
                skipped_is_canonical += 1
                continue

            writer.writerow([iso.id, iso.transcript_id, iso.gene_id, iso.sequence, canon_seq])
            written += 1

    db.close()

    print(f"\nStaged {written:,} isoforms → {args.output}")
    print(f"  Skipped (is canonical / identical seq): {skipped_is_canonical:,}")
    print(f"  Skipped (gene has no canonical):        {skipped_no_canonical:,}")
    file_mb = os.path.getsize(args.output) / (1024 * 1024)
    print(f"  File size: {file_mb:.1f} MB")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
build_md5_masterfile.py — catalogue every protein already predicted (by protein_md5)
so new AF3 batches only predict genuinely-novel sequences.

protein_md5 = md5(sequence.encode()).hexdigest()  on the raw uppercase protein seq.

"Already present" = the UNION of:
  * on-disk AF3 structures: structures/{gencode_v47,ccrcc,bc}/<md5>/<md5>_model.cif
    (the physical ground truth — a superset of what's ingested), and
  * DB isoforms that carry a Structure row.
(drug_complexes are protein+ligand co-folds, job-named not md5-named — excluded.)

Outputs (data/):
  * md5_masterfile.tsv   — protein_md5, structure_on_disk, in_db, seq_len,
                           n_transcripts, sources, example_gene, example_transcript
  * predicted_md5.txt    — sorted unique md5 list (for fast comm/grep dedup)

New-dataset workflow:
  1. for each new protein seq: md5=$(python -c "import hashlib,sys;print(hashlib.md5(sys.argv[1].encode()).hexdigest())" "$SEQ")
  2. novel = md5 NOT in predicted_md5.txt   (comm -23 <(sort new_md5s) predicted_md5.txt)
  3. feed the novel sequences to build_predict_queue.py (applies AF3 filters:
     30 ≤ len ≤ 2700, no X/B/Z/U/* residues).

Usage: python backend/scripts/build_md5_masterfile.py
"""
import os
import re
import sys
import csv
from collections import defaultdict
from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from sqlalchemy import func
from app.database import SessionLocal
from app.models.bio import Isoform, Gene, Structure

REPO = Path(__file__).resolve().parent.parent.parent
STRUCT_ROOT = REPO / "backend" / "static" / "data" / "structures"
PROTEIN_SOURCES = ["gencode_v47", "ccrcc", "bc"]
OUT_DIR = REPO / "data"
_MD5 = re.compile(r"^[0-9a-f]{32}$")


def scan_ondisk():
    """md5 -> set(source) for dirs that actually contain a model.cif."""
    present = defaultdict(set)
    for s in PROTEIN_SOURCES:
        d = STRUCT_ROOT / s
        if not d.exists():
            continue
        for name in os.listdir(d):
            if _MD5.match(name) and (d / name / f"{name}_model.cif").exists():
                present[name].add(s)
    return present


def main():
    db = SessionLocal()
    print("scanning on-disk structures …", flush=True)
    ondisk = scan_ondisk()
    print(f"  on-disk predicted md5s: {len(ondisk)}", flush=True)

    # DB metadata per md5
    print("aggregating DB metadata …", flush=True)
    has_struct = set(r[0] for r in db.query(Isoform.protein_md5)
                     .join(Structure, Structure.isoform_id == Isoform.id)
                     .filter(Isoform.protein_md5.isnot(None)).distinct().all())
    meta = {}  # md5 -> [seq_len, n_tx, gene, transcript]
    for md5, slen, sym, tx in db.query(
            Isoform.protein_md5, func.length(Isoform.sequence), Gene.symbol, Isoform.transcript_id)\
            .join(Gene, Gene.id == Isoform.gene_id)\
            .filter(Isoform.protein_md5.isnot(None)).all():
        e = meta.get(md5)
        if e is None:
            meta[md5] = [slen, 1, sym, tx]
        else:
            e[1] += 1
    db.close()

    all_md5 = set(ondisk) | has_struct | set(meta)
    OUT_DIR.mkdir(exist_ok=True)
    mf = OUT_DIR / "md5_masterfile.tsv"
    lst = OUT_DIR / "predicted_md5.txt"
    n_present = 0
    with open(mf, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein_md5", "structure_on_disk", "in_db", "seq_len",
                    "n_transcripts", "sources", "example_gene", "example_transcript"])
        for md5 in sorted(all_md5):
            m = meta.get(md5)
            on_disk = md5 in ondisk
            present = on_disk or (md5 in has_struct)
            if present:
                n_present += 1
            w.writerow([
                md5, int(on_disk), int(md5 in meta),
                m[0] if m else "", m[1] if m else "",
                ";".join(sorted(ondisk.get(md5, []))),
                m[2] if m else "", m[3] if m else "",
            ])
    # predicted set = anything with a structure (on disk or DB)
    predicted = sorted(set(ondisk) | has_struct)
    with open(lst, "w") as fh:
        fh.write("\n".join(predicted) + "\n")

    print(f"\n─── md5 masterfile ──────────────────────────────────────")
    print(f"total catalogued md5s        : {len(all_md5)}")
    print(f"  with a structure (present) : {len(predicted)}   ← skip these")
    print(f"  on disk                    : {len(ondisk)}")
    print(f"  in DB (any)                : {len(meta)}")
    print(f"  in DB but NO structure     : {len(set(meta) - has_struct)}")
    print(f"\nwrote:\n  {mf}\n  {lst}  ({len(predicted)} md5s)")


if __name__ == "__main__":
    main()

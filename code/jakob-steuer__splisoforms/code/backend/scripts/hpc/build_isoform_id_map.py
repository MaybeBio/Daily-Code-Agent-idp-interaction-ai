#!/usr/bin/env python3
"""
build_isoform_id_map.py — translate isoform ids between two builds of the DB.

A full rebuild renumbers `isoforms.id`. HPC neoantigen parquet shards are keyed
by the numeric isoform_id of the build they were staged from, so to load them
into a rebuilt DB each old id must be mapped to the new id of the SAME protein
(`protein_md5` is the stable identity).

A mapping is emitted only when the prediction is still valid in the new build:

  * the protein (md5) still exists as an Isoform in the new DB;
  * it is not the canonical isoform of its gene in the new DB (canonicals carry
    no neoepitopes by definition);
  * its gene's canonical PROTEIN is unchanged between the builds — neoepitopes
    are "novel relative to the canonical", so a changed canonical invalidates
    the peptide set and the isoform must be re-staged for the HPC.

Usage:
    python backend/scripts/hpc/build_isoform_id_map.py \
        --old postgresql://…/splisoforms --new postgresql://…/splisoforms_v2 \
        --out data/neoantigen_rebuild/id_map_old_to_v2.tsv
Prints a summary; writes <out> (old_isoform_id, new_isoform_id) and
<out>.restage.txt (new isoform ids that need a fresh HPC prediction).
"""
import argparse
import sys

from sqlalchemy import create_engine, text

Q = """
SELECT i.id, i.protein_md5, i.is_canonical, c.protein_md5 AS canonical_md5
  FROM isoforms i
  JOIN genes g ON g.id = i.gene_id
  LEFT JOIN isoforms c ON c.id = g.canonical_transcript_id
"""


def load(url: str):
    eng = create_engine(url)
    with eng.connect() as con:
        rows = con.execute(text(Q)).fetchall()
    return {r[1]: (r[0], bool(r[2]), r[3]) for r in rows}      # md5 → (id, is_canonical, canonical_md5)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", required=True, help="SQLAlchemy URL of the build the parquet was staged from")
    ap.add_argument("--new", required=True, help="SQLAlchemy URL of the rebuilt DB")
    ap.add_argument("--out", required=True, help="output TSV (old_isoform_id, new_isoform_id)")
    args = ap.parse_args()

    old = load(args.old)
    new = load(args.new)
    print(f"old build: {len(old):,} isoforms | new build: {len(new):,} isoforms")

    mapped = dropped_missing = dropped_canonical = dropped_canon_changed = 0
    restage = []
    with open(args.out, "w") as fh:
        fh.write("old_isoform_id\tnew_isoform_id\n")
        for md5, (old_id, old_canon, old_cmd5) in old.items():
            if old_canon:
                continue                      # never had predictions
            entry = new.get(md5)
            if entry is None:
                dropped_missing += 1
                continue
            new_id, new_canon, new_cmd5 = entry
            if new_canon:
                dropped_canonical += 1
                continue
            if old_cmd5 != new_cmd5:
                dropped_canon_changed += 1
                restage.append(new_id)
                continue
            fh.write(f"{old_id}\t{new_id}\n")
            mapped += 1

    # New non-canonical isoforms that did not exist in the old build also need staging
    old_md5s = set(old)
    new_only = [v[0] for m, v in new.items() if m not in old_md5s and not v[1]]
    restage.extend(new_only)
    with open(args.out + ".restage.txt", "w") as fh:
        for i in sorted(set(restage)):
            fh.write(f"{i}\n")

    print(f"mapped                          : {mapped:,}")
    print(f"dropped — protein gone          : {dropped_missing:,}")
    print(f"dropped — now canonical         : {dropped_canonical:,}")
    print(f"dropped — canonical changed     : {dropped_canon_changed:,}  (need re-prediction)")
    print(f"new non-canonical isoforms      : {len(new_only):,}  (need prediction)")
    print(f"→ {args.out}\n→ {args.out}.restage.txt ({len(set(restage)):,} new isoform ids to stage)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

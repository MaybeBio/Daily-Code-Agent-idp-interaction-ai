#!/usr/bin/env python3
"""
fix_ccrcc_structure_naming.py — normalise ccRCC structure dir naming to md5.

Background (2026-06 audit): the ccRCC AF3 export named structure directories
inconsistently — ~23k by the full protein md5, but ~10k by the AF3 *job_name*
(`GENE_md5prefix`, lower-cased), e.g. `a4galt_2ab1f08db0ce`. The DB and the
serving layer address structures by full md5, so the job_name-named dirs appear
as "missing structures" (and the real structures appear "orphaned").

The af3_inputs manifest (job_name → full md5) is the authoritative bridge. This
script renames every job_name-style dir (and its 3 inner files) to the full md5,
making all structures addressable by md5. Full-md5 dirs are left untouched.

It is idempotent and non-destructive: it only renames job_name dirs that resolve
to a manifest md5 and whose target md5 dir does not already exist. Run with
--dry-run first.

Usage:
    python scripts/fix_ccrcc_structure_naming.py --dry-run
    python scripts/fix_ccrcc_structure_naming.py
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CCRCC = REPO_ROOT / "data" / "sources" / "ccrcc"
STRUCT_DIR = CCRCC / "structures"
MANIFEST = CCRCC / "leonardo" / "manifest.tsv"

MD5_RE = re.compile(r"^[0-9a-f]{32}$")
INNER_SUFFIXES = ("_model.cif", "_confidences.json.gz",
                  "_confidences.json", "_summary_confidences.json")


def load_job2md5(manifest: Path) -> dict:
    """{job_name.lower(): full_md5} from the af3_inputs manifest."""
    out = {}
    with manifest.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_md5, i_job = header.index("md5"), header.index("job_name")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_md5, i_job):
                continue
            out[f[i_job].lower()] = f[i_md5]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not MANIFEST.exists():
        sys.exit(f"manifest not found: {MANIFEST}")

    job2md5 = load_job2md5(MANIFEST)
    print(f"manifest: {len(job2md5)} job_name→md5 entries")

    dirs = [p for p in STRUCT_DIR.iterdir() if p.is_dir()]
    job_dirs = [p for p in dirs if not MD5_RE.match(p.name)]
    print(f"structure dirs: {len(dirs)} total, {len(job_dirs)} job_name-style")

    renamed = skipped_exists = unmapped = 0
    for d in job_dirs:
        md5 = job2md5.get(d.name.lower())
        if md5 is None:
            print(f"  UNMAPPED (no manifest entry): {d.name}", file=sys.stderr)
            unmapped += 1
            continue
        target = STRUCT_DIR / md5
        if target.exists():
            print(f"  SKIP (target md5 dir exists): {d.name} → {md5}",
                  file=sys.stderr)
            skipped_exists += 1
            continue
        # rename inner files job_name_* → md5_*
        for child in d.iterdir():
            for suf in INNER_SUFFIXES:
                if child.name.endswith(suf):
                    new = d / f"{md5}{suf}"
                    if not args.dry_run:
                        child.rename(new)
                    break
        # rename the directory
        if not args.dry_run:
            d.rename(target)
        renamed += 1

    print(f"\n{'[dry-run] would rename' if args.dry_run else 'renamed'}: "
          f"{renamed}")
    print(f"skipped (target exists): {skipped_exists}")
    print(f"unmapped:                {unmapped}")
    if not args.dry_run:
        n_md5 = sum(1 for p in STRUCT_DIR.iterdir()
                    if p.is_dir() and MD5_RE.match(p.name))
        print(f"full-md5 dirs after fix: {n_md5}")


if __name__ == "__main__":
    main()

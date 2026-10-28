from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from idrfeat.io import read_fasta
from idrfeat.labels import load_residue_table


def coordinate_concordance(df: pd.DataFrame, seqs: dict[str, str]) -> dict:
    checked = matched = 0
    missing = set()
    for r in df.itertuples(index=False):
        s = seqs.get(r.accession)
        if s is None:
            missing.add(r.accession)
            continue
        if 1 <= r.aa_loc <= len(s):
            checked += 1
            matched += int(s[r.aa_loc - 1] == r.aa)
    return {
        "checked": checked,
        "matched": matched,
        "match_frac": matched / checked if checked else float("nan"),
        "missing_accessions": sorted(missing),
    }


def crux_table(df: pd.DataFrame) -> pd.DataFrame:
    covered = df[df["peptide_count"] > 0]
    groups = [("IDR", covered[covered["idr"] == 1]), ("ordered", covered[covered["idr"] == 0])]
    rows = []
    for label, sub in groups:
        rows.append(
            {
                "region": label,
                "n": int(len(sub)),
                "pct_depleted": float(sub["depleted"].mean() * 100) if len(sub) else float("nan"),
                "mean_lfc": float(sub["lfc"].mean()) if len(sub) else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="analyze_screen")
    parser.add_argument("--screen", required=True, help="per-residue screen table")
    parser.add_argument("--fasta", required=True, help="reference proteome FASTA (plain or .gz)")
    parser.add_argument("--outdir", default="data/processed/screen_qc")
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    df = load_residue_table(args.screen)
    seqs = read_fasta(args.fasta)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    genes = df["accession"].nunique()
    cov = df["peptide_count"] > 0
    conc = coordinate_concordance(df, seqs)
    crux = crux_table(df)
    crux.to_csv(out / "crux_idr_vs_ordered.csv", index=False)

    print(f"residues: {len(df)} | proteins: {genes}")
    print(f"coverage: {int(cov.sum())}/{len(df)} = {cov.mean() * 100:.0f}%")
    print(f"IDR residues: {int(df['idr'].sum())} ({df['idr'].mean() * 100:.0f}%)")
    covered = df[cov]
    print(f"depleted (covered): {int(covered['depleted'].sum())}/{len(covered)} "
          f"= {covered['depleted'].mean() * 100:.1f}%")
    print(f"coordinate match: {conc['matched']}/{conc['checked']} "
          f"= {conc['match_frac'] * 100:.1f}%")
    if conc["missing_accessions"]:
        print(f"proteins not in reference: {conc['missing_accessions']}")
    print("crux (depletion by region):")
    for r in crux.itertuples(index=False):
        print(f"  {r.region}: {r.pct_depleted:.1f}% depleted, mean LFC {r.mean_lfc:.3f} (n={r.n})")
    print(f"wrote {out}/crux_idr_vs_ordered.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

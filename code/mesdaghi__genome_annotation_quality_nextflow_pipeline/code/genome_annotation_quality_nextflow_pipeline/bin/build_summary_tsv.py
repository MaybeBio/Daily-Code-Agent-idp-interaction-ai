#!/usr/bin/env python3
"""
build_summary_tsv.py

Write a one-row TSV summarising the query proteome with the headline
numbers shown on the pipeline's plots:

    prop_ge_70           — fraction of proteins with mean pLDDT >= 70
                           (pLDDT branch — y-axis of GMM scatter)
    gmm_upper_prop       — GMM upper-component weight from a 2-component
                           GMM fitted with guided init from the reference
                           species (pLDDT branch — x-axis of GMM scatter)
    disorder_prop_gt_0p5 — proportion of proteins with mean metapredict
                           disorder score > 0.5 (metapredict branch)
    psauron_high_pct     — % of proteins with psauron in-frame_score >= 0.9
                           (psauron branch — "High" bar)
    ipr_hit_pct          — % of proteins with at least one InterProScan
                           IPR hit (InterPro branch — left-hand bar)

Every branch is optional: any metric whose inputs are not supplied is
left empty (the column is still emitted so downstream parsers see a
stable schema). total_proteins comes from the pLDDT PKL when it is given,
otherwise from counting the records in --fasta.

Usage:
    python3 build_summary_tsv.py \\
        --name DATASET \\
        [--fasta DATASET.fasta] \\
        [--plddt-pkl plddt_all_values_DATASET.pkl \\
         --plddt-reference-csv reference/plddt_model_organisms.csv] \\
        [--metapredict-csv DATASET_metapredict.csv] \\
        [--psauron-csv DATASET_psauron.csv] \\
        [--ipr-fasta DATASET.fasta --ipr-xml DATASET.xml \\
         --plot-interpro-dir bin] \\
        --out DATASET_summary.tsv
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture


# Same reference species used by plot_plddt.py for guided GMM init.
REFERENCE_SPECIES = [
    "Drosophila_melanogaster",
    "Homo_sapiens",
    "Saccharomyces_cerevisiae",
    "Mus_musculus",
    "Arabidopsis_thaliana",
]


def fit_reference_gmm_init(reference_csv: Path):
    """Replicate plot_plddt.py's guided init: average the means/weights
    from 2-component GMMs fitted to each reference species' pLDDT values."""
    df = pd.read_csv(reference_csv)
    df["Species"] = df["Species"].replace({"HS": "Homo_sapiens"})

    means_list, weights_list = [], []
    for sp in REFERENCE_SPECIES:
        sub = df[df["Species"] == sp]
        if len(sub) < 20:
            print(f"  WARNING: reference species '{sp}' has <20 rows, skipping init contribution")
            continue
        X = sub["Mean_pLDDT"].to_numpy().reshape(-1, 1)
        gmm = GaussianMixture(n_components=2, random_state=0).fit(X)
        means = gmm.means_.flatten()
        order = np.argsort(means)
        means_list.append(means[order])
        weights_list.append(gmm.weights_[order])

    if not means_list:
        return None, None
    return np.mean(means_list, axis=0), np.mean(weights_list, axis=0)


def compute_plddt_metrics(pkl_path: Path, reference_csv: Path) -> dict:
    with open(pkl_path, "rb") as f:
        species_plddt = pickle.load(f)

    # Pool all mean-pLDDT values from the PKL. A query PKL normally has a
    # single species key, but if there are multiple they're pooled here so
    # the TSV always carries one row per dataset.
    species_in_pkl = list(species_plddt.keys())
    if len(species_in_pkl) > 1:
        print(f"  NOTE: PKL contains {len(species_in_pkl)} species, pooling: {species_in_pkl}")

    values = []
    for _, protein_dict in species_plddt.items():
        values.extend(protein_dict.values())
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        raise ValueError(f"No pLDDT values in {pkl_path}")

    prop_ge_70 = float((values >= 70).mean())

    ref_means, ref_weights = fit_reference_gmm_init(reference_csv)
    X = values.reshape(-1, 1)
    if ref_means is not None:
        gmm = GaussianMixture(
            n_components=2,
            means_init=ref_means.reshape(-1, 1),
            weights_init=ref_weights,
            random_state=0,
        ).fit(X)
    else:
        print("  WARNING: no reference GMM init available, falling back to default")
        gmm = GaussianMixture(n_components=2, random_state=0).fit(X)

    means = gmm.means_.flatten()
    order = np.argsort(means)
    gmm_upper_prop = float(gmm.weights_[order][1])

    return {
        "total_proteins": int(values.size),
        "prop_ge_70": round(prop_ge_70, 4),
        "gmm_upper_prop": round(gmm_upper_prop, 4),
    }


def compute_disorder_prop_gt_0p5(metapredict_csv: Path) -> float:
    """Proportion of proteins with mean metapredict disorder score > 0.5.

    Per-row format (matches plot_metapredict.py): field 0 is the header,
    field 1 is the sequence, fields 2+ are per-residue disorder scores.
    The mean of fields 2+ is the per-protein mean disorder; we then
    compute the fraction of those means strictly greater than 0.5.
    """
    means = []
    skipped = 0
    with open(metapredict_csv) as fh:
        for line in fh:
            parts = line.strip().split(",")
            try:
                scores = [float(x) for x in parts[2:] if x.strip()]
            except ValueError:
                skipped += 1
                continue
            if scores:
                means.append(float(np.mean(scores)))
            else:
                skipped += 1
    if not means:
        return float("nan")
    arr = np.asarray(means, dtype=float)
    if skipped:
        print(f"  (skipped {skipped} malformed/empty rows)")
    return round(float((arr > 0.5).mean()), 4)


def compute_psauron_high_pct(psauron_csv: Path) -> float:
    """% of proteins with in-frame_score >= 0.9.

    Matches plot_psauron_distribution.py: read with skiprows=2, use the
    'in-frame_score' column.
    """
    df = pd.read_csv(psauron_csv, skiprows=2)
    if "in-frame_score" not in df.columns:
        raise ValueError(
            f"'in-frame_score' column not found in {psauron_csv}. "
            f"Found columns: {list(df.columns)}"
        )
    scores = df["in-frame_score"].dropna().to_numpy()
    if scores.size == 0:
        return float("nan")
    return round(float((scores >= 0.9).mean()) * 100, 3)


def count_fasta_records(fasta_path: Path) -> int:
    """Number of '>' header lines in a FASTA file."""
    with open(fasta_path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def compute_ipr_hit_pct(fasta_path: Path, xml_path: Path,
                        plot_interpro_dir: Path) -> float:
    """% of FASTA proteins with at least one IPR-bearing match.

    Reuses parse_fasta_lengths / parse_xml from plot_interpro.py so the
    number always matches what the InterPro plot's bar chart shows.
    """
    sys.path.insert(0, str(plot_interpro_dir))
    try:
        from plot_interpro import parse_fasta_lengths, parse_xml
    except ImportError as e:
        sys.exit(
            f"ERROR: cannot import plot_interpro from {plot_interpro_dir}: {e}"
        )

    fasta_lengths = parse_fasta_lengths(fasta_path)
    total = len(fasta_lengths)
    if total == 0:
        raise ValueError(f"No proteins in {fasta_path}")

    _, with_ipr, _, _ = parse_xml(xml_path)
    return round(with_ipr / total * 100, 3)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--name", required=True,
                    help="Dataset name (e.g. the query identifier)")
    ap.add_argument("--fasta", type=Path, default=None,
                    help="Query FASTA; used for total_proteins when the "
                         "pLDDT branch did not run")
    ap.add_argument("--plddt-pkl", type=Path, default=None,
                    help="Path to plddt_all_values_<dataset>.pkl (only "
                         "needed if the Protenix branch ran)")
    ap.add_argument("--plddt-reference-csv", type=Path, default=None,
                    help="Path to reference/plddt_model_organisms.csv "
                         "(needed for guided GMM init; required with "
                         "--plddt-pkl)")
    ap.add_argument("--psauron-csv", type=Path, default=None,
                    help="Path to the query psauron CSV (with the 2-line "
                         "header that plot_psauron_distribution.py skips; "
                         "only needed if the psauron branch ran)")
    ap.add_argument("--metapredict-csv", type=Path, default=None,
                    help="Path to the query metapredict CSV (per-row "
                         "per-residue disorder scores from "
                         "METAPREDICT_DISORDER; only needed if the "
                         "metapredict branch ran)")
    ap.add_argument("--ipr-fasta", type=Path, default=None,
                    help="Query FASTA (only needed if InterPro branch ran)")
    ap.add_argument("--ipr-xml", type=Path, default=None,
                    help="Query InterProScan XML (only needed if InterPro "
                         "branch ran)")
    ap.add_argument("--plot-interpro-dir", type=Path, default=None,
                    help="Directory containing plot_interpro.py (required "
                         "with --ipr-xml/--ipr-fasta)")
    ap.add_argument("--out", required=True, type=Path,
                    help="Output TSV path")
    args = ap.parse_args()

    print(f"=== build_summary_tsv :: {args.name} ===")

    # pLDDT branch (optional)
    if args.plddt_pkl:
        if args.plddt_reference_csv is None:
            sys.exit("--plddt-reference-csv is required when --plddt-pkl is given")
        print("\npLDDT metrics:")
        plddt = compute_plddt_metrics(args.plddt_pkl, args.plddt_reference_csv)
        print(f"  total_proteins = {plddt['total_proteins']:,}")
        print(f"  prop_ge_70     = {plddt['prop_ge_70']}")
        print(f"  gmm_upper_prop = {plddt['gmm_upper_prop']}")
    else:
        print("\npLDDT metrics: skipped (no PKL provided)")
        plddt = {"total_proteins": "", "prop_ge_70": "", "gmm_upper_prop": ""}
        if args.fasta:
            plddt["total_proteins"] = count_fasta_records(args.fasta)
            print(f"  total_proteins = {plddt['total_proteins']:,} (from FASTA)")

    # psauron branch (optional)
    if args.psauron_csv:
        print("\npsauron metric:")
        psauron_high_pct = compute_psauron_high_pct(args.psauron_csv)
        print(f"  psauron_high_pct = {psauron_high_pct}")
    else:
        print("\npsauron metric: skipped (no CSV provided)")
        psauron_high_pct = ""

    # metapredict branch (optional)
    if args.metapredict_csv:
        print("\nmetapredict metric:")
        disorder_prop_gt_0p5 = compute_disorder_prop_gt_0p5(args.metapredict_csv)
        print(f"  disorder_prop_gt_0p5 = {disorder_prop_gt_0p5}")
    else:
        print("\nmetapredict metric: skipped (no CSV provided)")
        disorder_prop_gt_0p5 = ""

    # InterPro branch (optional)
    has_ipr_inputs = bool(args.ipr_fasta) and bool(args.ipr_xml)
    if has_ipr_inputs:
        if args.plot_interpro_dir is None:
            sys.exit(
                "--plot-interpro-dir is required when --ipr-xml/--ipr-fasta "
                "are given"
            )
        print("\nInterPro metric:")
        ipr_hit_pct = compute_ipr_hit_pct(
            args.ipr_fasta, args.ipr_xml, args.plot_interpro_dir
        )
        print(f"  ipr_hit_pct = {ipr_hit_pct}")
    else:
        print("\nInterPro metric: skipped (no XML/FASTA provided)")
        ipr_hit_pct = ""  # leave column blank, keep schema stable

    # Emit TSV (header + one row)
    cols = ["dataset", "total_proteins", "prop_ge_70", "gmm_upper_prop",
            "disorder_prop_gt_0p5", "psauron_high_pct", "ipr_hit_pct"]
    row = [args.name, plddt["total_proteins"], plddt["prop_ge_70"],
           plddt["gmm_upper_prop"], disorder_prop_gt_0p5,
           psauron_high_pct, ipr_hit_pct]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        f.write("\t".join(cols) + "\n")
        f.write("\t".join(str(v) for v in row) + "\n")

    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()

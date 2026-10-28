from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

DATA = Path("data/depmap")
SUPP = Path("results/depmap_essentiality/supp")
SUPP.mkdir(parents=True, exist_ok=True)

CONFOUNDS = [
    "protein_length", "total_idr_len", "n_idr_segments", "idr_fraction",
    "disorder_mean", "disorder_frac", "metapredict_mean", "metapredict_frac",
    "aiupred_disorder_mean", "plddt_mean",
]
ANNOTATION_PREFIXES = ("elm_", "ptm_", "phylop", "phastcons")


def strip_entrez(name: str) -> str:
    return re.sub(r"\s*\(\d+\)$", "", str(name)).strip()


def load_control_set(path: Path) -> set[str]:
    genes = set()
    with open(path) as f:
        next(f)
        for line in f:
            s = line.strip()
            if s:
                genes.add(strip_entrez(s))
    return genes


AA_NAME = {
    "A": "alanine", "C": "cysteine", "D": "aspartate", "E": "glutamate", "F": "phenylalanine",
    "G": "glycine", "H": "histidine", "I": "isoleucine", "K": "lysine", "L": "leucine",
    "M": "methionine", "N": "asparagine", "P": "proline", "Q": "glutamine", "R": "arginine",
    "S": "serine", "T": "threonine", "V": "valine", "W": "tryptophan", "Y": "tyrosine",
}

DESC = {
    "seg_length": ("IDR geometry", "metapredict", "Mean length of the disordered segments"),
    "frac_charged": ("Composition", "localCIDER", "Fraction of charged residues (D,E,K,R)"),
    "frac_polar": ("Composition", "localCIDER", "Fraction of polar uncharged residues"),
    "frac_aromatic": ("Composition", "localCIDER", "Fraction of aromatic residues (F,W,Y)"),
    "frac_proline": ("Composition", "localCIDER", "Fraction of proline"),
    "frac_glycine": ("Composition", "localCIDER", "Fraction of glycine"),
    "frac_disorder_promoting": ("Composition", "localCIDER", "Fraction of disorder-promoting residues"),
    "fcr": ("Charge", "localCIDER", "Fraction of charged residues"),
    "ncpr": ("Charge", "localCIDER", "Net charge per residue"),
    "net_charge": ("Charge", "localCIDER", "Total net charge of the IDR"),
    "kappa": ("Charge", "localCIDER", "Charge patterning: how blocky + and - charges are (0 mixed, 1 segregated)"),
    "omega": ("Charge", "localCIDER", "Patterning of charged plus proline residues vs the rest"),
    "scd": ("Charge", "localCIDER", "Sequence charge decoration, a charge-patterning index"),
    "frac_positive": ("Charge", "localCIDER", "Fraction of positively charged residues (K,R)"),
    "frac_negative": ("Charge", "localCIDER", "Fraction of negatively charged residues (D,E)"),
    "hydropathy_mean": ("Hydropathy", "localCIDER", "Mean Kyte-Doolittle hydropathy of the IDR"),
    "lowcomplexity_frac": ("Low complexity", "NCBI SEG", "Fraction of the IDR in low-complexity sequence"),
    "longest_single_run": ("Low complexity", "sequence", "Longest single-residue run in the IDR"),
    "disorder_mean": ("Disorder", "metapredict", "Mean disorder score across the protein"),
    "disorder_max": ("Disorder", "metapredict", "Maximum disorder score"),
    "disorder_frac": ("Disorder", "metapredict", "Fraction of the protein predicted disordered"),
    "metapredict_mean": ("Disorder", "metapredict", "Mean metapredict disorder score"),
    "metapredict_max": ("Disorder", "metapredict", "Maximum metapredict disorder score"),
    "metapredict_frac": ("Disorder", "metapredict", "Fraction disordered by metapredict"),
    "elm_motif_residues": ("Motifs", "ELM", "Residues covered by short linear motifs"),
    "elm_motif_density": ("Motifs", "ELM", "Short linear motifs per residue"),
    "elm_classes_distinct": ("Motifs", "ELM", "Number of distinct motif classes present"),
    "ptm_count": ("Modification sites", "dbPTM", "Total modification sites in the IDR"),
    "ptm_density": ("Modification sites", "dbPTM", "Modification sites per residue"),
    "ptm_phospho_count": ("Modification sites", "dbPTM", "Phosphorylation sites"),
    "ptm_phospho_density": ("Modification sites", "dbPTM", "Phosphorylation sites per residue"),
    "ptm_acetyl_count": ("Modification sites", "dbPTM", "Acetylation sites"),
    "ptm_acetyl_density": ("Modification sites", "dbPTM", "Acetylation sites per residue"),
    "ptm_ubiq_count": ("Modification sites", "dbPTM", "Ubiquitination sites"),
    "ptm_ubiq_density": ("Modification sites", "dbPTM", "Ubiquitination sites per residue"),
    "ptm_methyl_count": ("Modification sites", "dbPTM", "Methylation sites"),
    "ptm_methyl_density": ("Modification sites", "dbPTM", "Methylation sites per residue"),
    "n_idr_segments": ("Size and disorder baseline", "metapredict", "Number of disordered segments"),
    "total_idr_len": ("Size and disorder baseline", "metapredict", "Total disordered length"),
    "protein_length": ("Size and disorder baseline", "sequence", "Protein length"),
    "idr_fraction": ("Size and disorder baseline", "metapredict", "Fraction of the protein that is disordered"),
}


def classify(col: str):
    if col.startswith("aa_frac_"):
        aa = col.split("_")[-1]
        return ("Composition", "localCIDER", f"Fraction of {AA_NAME.get(aa, aa)} ({aa})")
    return DESC.get(col, ("Other", "sequence", col))


def main() -> None:
    gl = pd.read_csv(DATA / "gene_level_idr_features.csv", index_col=0)
    feats = [c for c in gl.columns if gl[c].dtype != object]
    keep = [c for c in feats if gl[c].notna().mean() >= 0.5]
    confound = [c for c in CONFOUNDS if c in keep]
    character = [c for c in keep if c not in confound]
    seq_only = [c for c in keep if not c.startswith(ANNOTATION_PREFIXES)]

    ess = load_control_set(DATA / "AchillesCommonEssentialControls.csv")
    non = load_control_set(DATA / "AchillesNonessentialControls.csv")
    label = pd.Series(index=gl.index, dtype="float64")
    label[gl.index.isin(ess)] = 1.0
    label[gl.index.isin(non)] = 0.0
    train = gl[label.notna()]
    ytr = label[label.notna()].to_numpy().astype(int)

    rows = []
    for c in keep:
        group, tool, desc = classify(c)
        x = train[c].to_numpy(dtype=float)
        ok = np.isfinite(x)
        if ok.sum() >= 50 and len(np.unique(x[ok])) >= 3:
            a = roc_auc_score(ytr[ok], x[ok])
            uauc = round(float(max(a, 1 - a)), 3)
            direction = "higher in essential" if a >= 0.5 else "higher in non-essential"
        else:
            uauc, direction = float("nan"), "NA"
        rows.append({
            "feature": c, "group": group, "tool": tool, "description": desc,
            "role": "baseline" if c in confound else "character",
            "source": "annotation" if c.startswith(ANNOTATION_PREFIXES) else "sequence-derived",
            "coverage": round(float(gl[c].notna().mean()), 3),
            "univariate_auc": uauc,
            "direction": direction,
        })
    df = pd.DataFrame(rows)
    df["univariate_auc"] = pd.to_numeric(df["univariate_auc"], errors="coerce")
    df_sorted = df.sort_values("univariate_auc", ascending=False, na_position="last")
    df.to_csv(SUPP / "feature_table.csv", index=False)
    df_sorted.to_csv(SUPP / "feature_table_by_auc.csv", index=False)

    print(f"features kept {len(keep)} | character {len(character)} | baseline {len(confound)} | "
          f"sequence-only {len(seq_only)} | annotation {len(keep) - len(seq_only)}")
    print(f"wrote {SUPP}/feature_table.csv and feature_table_by_auc.csv")
    print("\nTop 12 features by univariate ROC-AUC:")
    print(f"{'feature':24s} {'group':22s} {'tool':12s} {'auc':>5s}")
    for _, r in df_sorted.head(12).iterrows():
        print(f"{r['feature']:24s} {r['group']:22s} {r['tool']:12s} {str(r['univariate_auc']):>5s}")


if __name__ == "__main__":
    main()

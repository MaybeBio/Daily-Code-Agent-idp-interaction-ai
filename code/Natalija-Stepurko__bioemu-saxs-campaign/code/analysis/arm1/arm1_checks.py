"""Checks on the Arm 1 replication list: does a large raw chi2 reflect the data's precision, and are the
twelve listed systems related to one another?

(i) Precision. For every entry, the experimental relative error is the median sigma/I over the lower half
    of the measured q range (q <= (q_min + q_max) / 2, points with I > 0 and sigma > 0). Over the clean
    BioEmu-1 entries: Spearman correlation of log raw chi2 with the relative error; median raw chi2 by
    SASBDB curve type and by source type (Kruskal-Wallis across groups of at least five; Mann-Whitney for
    in-house against synchrotron). For the twelve Arm 1 systems against the eligible pool: the median
    relative error, NRMSD and raw chi2 of the twelve, the pool medians, each system's percentile in the
    pool, the percentile of the twelve's median in the pool, and a Mann-Whitney test of the twelve against
    the rest of the pool.
(ii) Relatedness. Organism (SASBDB record, counted at species level, and the organism group of the
    coverage analysis's UniProt mapping), UniProt accession, CATH-Gene3D superfamilies, the near-duplicate
    cluster and the cross-validation sequence group (results/clusters.csv), and the largest pairwise 5-mer
    containment between any two of the twelve.

Inputs: the measured curves (BSC_ROOT/data, from `bsc fetch`), and the tracked results/resolvability.csv,
results/candidates.csv, results/campaign_arms.json, results/clusters.csv, results/entries.csv,
results/coverage/sasbdb_metadata.csv (analysis/coverage/18_sasbdb_metadata.py) and
results/coverage/per_protein_mapping.csv.
Outputs (tracked): results/arm1/relative_error.csv, results/arm1/precision.json,
results/arm1/relatedness.csv, results/arm1/relatedness.json.

    BSC_ROOT=<where bsc fetch put data/> uv run python analysis/arm1/arm1_checks.py
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from bsc import config as C
from bsc import report, robust
from bsc.score import read_experiment

REPO = Path(__file__).resolve().parents[2]
RES = REPO / "results"
OUT = RES / "arm1"
REL_ERR_DEF = "median sigma/I over q <= (q_min + q_max) / 2, I > 0 and sigma > 0"


def relative_error(q: np.ndarray, I: np.ndarray, sigma: np.ndarray) -> tuple[float, int]:
    """Median sigma/I over the lower half of the measured q range (I > 0, sigma > 0): (value, points)."""
    ok = (I > 0) & (sigma > 0)
    q, I, sigma = q[ok], I[ok], sigma[ok]
    if len(q) == 0:
        return float("nan"), 0
    k = q <= (q.min() + q.max()) / 2
    return float(np.median(sigma[k] / I[k])), int(k.sum())


def percentile_in(pool: np.ndarray, x: float) -> float:
    """Share of the pool at or below x, in percent (the 'weak' percentile of scipy)."""
    return float(stats.percentileofscore(pool, x, kind="weak"))


def _spearman(x: pd.Series, y: pd.Series) -> dict:
    r = stats.spearmanr(x, y)
    return {"rho": float(r.statistic), "p": float(r.pvalue)}


def precision(d: pd.DataFrame, meta: pd.DataFrame, cand: pd.DataFrame, arm1: list[str]) -> dict:
    m = d[(d["model"] == C.MODEL) & d["clean"]]
    m = m.merge(meta[["label", "curve", "source_type", "sec"]], on="label")
    out = {"n": int(len(m)), "rel_err_definition": REL_ERR_DEF,
           "spearman_log_chi2_vs_rel_err": _spearman(np.log(m["chi2_raw"]), m["rel_err"]),
           "spearman_nrmsd_vs_rel_err": _spearman(m["nrmsd_raw"], m["rel_err"]),
           "rel_err_median": float(m["rel_err"].median())}
    for col in ("curve", "source_type"):
        g = {k: v["chi2_raw"] for k, v in m.groupby(col)}
        out[f"chi2_raw_by_{col}"] = {k: {"n": int(len(v)), "chi2_raw_median": float(v.median()),
                                         "rel_err_median": float(m.loc[v.index, "rel_err"].median())}
                                     for k, v in sorted(g.items(), key=lambda kv: -len(kv[1]))}
        big = [v for v in g.values() if len(v) >= 5]
        out[f"chi2_raw_by_{col}_kruskal_p"] = float(stats.kruskal(*big).pvalue) if len(big) > 1 else None
    ih = m.loc[m["source_type"] == "X-ray in house", "chi2_raw"]
    sy = m.loc[m["source_type"] == "X-ray synchrotron", "chi2_raw"]
    out["in_house_vs_synchrotron_chi2_raw_mannwhitney_p"] = float(stats.mannwhitneyu(ih, sy).pvalue)
    # the twelve Arm 1 systems against the eligible pool
    pool = m[m["label"].isin(report.arm1_eligible(cand)["label"])].set_index("label")
    twelve = pool.loc[arm1]
    rest = pool.drop(index=arm1)
    rec = {"n_pool": int(len(pool)), "n_listed": int(len(twelve))}
    for col in ("rel_err", "nrmsd_raw", "chi2_raw"):
        p = pool[col].to_numpy()
        rec[col] = {"listed_median": float(twelve[col].median()), "pool_median": float(pool[col].median()),
                    "rest_median": float(rest[col].median()),
                    "listed_median_percentile_in_pool": percentile_in(p, float(twelve[col].median())),
                    "mannwhitney_listed_vs_rest_p": float(stats.mannwhitneyu(twelve[col], rest[col]).pvalue),
                    "per_system_percentile": {lab: percentile_in(p, float(twelve.at[lab, col]))
                                              for lab in arm1}}
    rec["listed_curve_types"] = twelve["curve"].value_counts().to_dict()
    rec["listed_source_types"] = twelve["source_type"].value_counts().to_dict()
    rec["pool_source_types"] = pool["source_type"].value_counts().to_dict()
    out["arm1_against_pool"] = rec
    return out


def relatedness(ent: pd.DataFrame, meta: pd.DataFrame, mp: pd.DataFrame, cl: pd.DataFrame,
                arm1: list[str]) -> tuple[pd.DataFrame, dict]:
    e, md, mpi, c = (x.set_index("label") for x in (ent, meta, mp, cl))
    t = pd.DataFrame([{"label": lab, "length": int(e.at[lab, "length"]),
                       "sasbdb_organism": md.at[lab, "organisms"],
                       "molecule": md.at[lab, "mol_names"], "uniprot": mpi.at[lab, "uniprot"],
                       "org_group": mpi.at[lab, "org_group"],
                       "gene3d_superfamilies": mpi.at[lab, "gene3d_superfamilies"],
                       "dup_cluster": int(c.at[lab, "dup_cluster"]),
                       "sequence_group": int(c.at[lab, "group"])}
                      for lab in arm1])
    t["species"] = t["sasbdb_organism"].str.split().str[:2].str.join(" ")     # strain designations dropped
    sets = {lab: robust.kmer_set(e.at[lab, "sequence"]) for lab in arm1}
    pairs = [(a, b, robust.containment(sets[a], sets[b])) for a, b in itertools.combinations(arm1, 2)]
    top = max(pairs, key=lambda x: x[2])
    sf = t["gene3d_superfamilies"].fillna("").str.split(";").map(lambda x: {y for y in x if y})
    shared_sf = {s: [lab for lab, ss in zip(t["label"], sf, strict=True) if s in ss]
                 for s in set().union(*sf) if sum(s in ss for ss in sf) > 1}
    summary = {"n": len(arm1), "n_species_recorded": int(t["species"].nunique()),
               "n_without_organism": int(t["species"].isna().sum()),
               "species": t["species"].value_counts().to_dict(),
               "organism_groups": t["org_group"].value_counts().to_dict(),
               "n_uniprot": int(t["uniprot"].nunique()),
               "n_sequence_groups": int(t["sequence_group"].nunique()),
               "n_dup_clusters": int(t["dup_cluster"].nunique()),
               "max_pairwise_5mer_containment": {"pair": [top[0], top[1]], "value": float(top[2])},
               "near_duplicate_threshold_jaccard": robust.NEAR_DUPLICATE,
               "group_threshold_containment": robust.GROUP_SIMILARITY,
               "superfamilies_shared_by_two_or_more": shared_sf}
    return t, summary


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ent = pd.read_csv(RES / "entries.csv")
    rows = []
    for lab in ent["label"]:
        x = read_experiment(lab)
        v, n = relative_error(x["q"].to_numpy(float), x["I"].to_numpy(float), x["sigma"].to_numpy(float))
        rows.append({"label": lab, "rel_err": v, "n_points_lower_half": n})
    rel = pd.DataFrame(rows)
    rel.to_csv(OUT / "relative_error.csv", index=False)
    d = pd.read_csv(RES / "resolvability.csv").merge(rel, on="label")
    d["clean"] = d["clean"].astype(bool)
    meta = pd.read_csv(RES / "coverage" / "sasbdb_metadata.csv")
    cand = pd.read_csv(RES / "candidates.csv")
    arm1 = [r["label"] for r in json.load(open(RES / "campaign_arms.json"))["arm1"]["replication"]]
    top = report.arm1_eligible(cand).sort_values("chi2_raw", ascending=False).head(report.ARM1_LISTED)
    assert arm1 == top["label"].tolist(), "the eligibility rule no longer gives the listed twelve"
    counts = json.load(open(RES / "analysis.json"))["bioemu"]["resolvability_counts"]
    assert d.loc[d["model"] == C.MODEL, "clean"].sum() == sum(counts.values()), "clean BioEmu-1 count"
    P = precision(d, meta, cand, arm1)
    json.dump(P, open(OUT / "precision.json", "w"), indent=1)
    t, S = relatedness(ent, meta, pd.read_csv(RES / "coverage" / "per_protein_mapping.csv"),
                       pd.read_csv(RES / "clusters.csv"), arm1)
    t.to_csv(OUT / "relatedness.csv", index=False)
    json.dump(S, open(OUT / "relatedness.json", "w"), indent=1)
    a = P["arm1_against_pool"]
    print(f"rho(log chi2, rel err) = {P['spearman_log_chi2_vs_rel_err']['rho']:+.2f} (n = {P['n']}); "
          f"Arm 1 rel err median {a['rel_err']['listed_median']:.4f} "
          f"(pool {a['rel_err']['pool_median']:.4f}); "
          f"{S['n_species_recorded']} species, {S['n_sequence_groups']} sequence groups")


if __name__ == "__main__":
    main()

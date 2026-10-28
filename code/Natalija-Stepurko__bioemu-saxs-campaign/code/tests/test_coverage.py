"""The coverage stage: structure groups, intervals and the stratified test, on synthetic data."""
import numpy as np
import pandas as pd
from scipy import stats

from bsc import coverage


def test_structure_groups_from_the_sifts_flags():
    mp = pd.DataFrame({"map_method": ["SASBDB", "SASBDB", "SASBDB", "unmapped"],
                       "has_pdb_50": [True, False, "False", np.nan],
                       "has_pdb_any_overlap": [True, True, False, np.nan]})
    g = coverage.structure_group(mp)
    assert g.iloc[:3].tolist() == ["pdb50", "partial", "none"] and pd.isna(g.iloc[3])


def test_wilson_interval_matches_known_values_and_stays_in_range():
    lo, hi = coverage.wilson(118, 290)
    assert abs(lo - 0.352) < 0.002 and abs(hi - 0.464) < 0.002
    assert coverage.wilson(0, 11)[0] == 0.0 and 0 < coverage.wilson(0, 11)[1] < 0.3


def test_van_elteren_with_one_stratum_is_the_rank_sum_test():
    rng = np.random.default_rng(3)
    x = np.r_[rng.normal(0, 1, 40), rng.normal(0.6, 1, 50)]
    grp = np.r_[np.ones(40, bool), np.zeros(50, bool)]
    ve = coverage.van_elteren(x, grp, np.zeros(90))
    mw = stats.mannwhitneyu(x[grp], x[~grp], alternative="two-sided", method="asymptotic",
                            use_continuity=False)
    assert abs(ve["p_two_sided"] - mw.pvalue) < 1e-9


def test_a_confounded_pooled_difference_vanishes_when_stratified():
    """Structure availability tracks the class, the class sets the error, structure has no effect."""
    rng = np.random.default_rng(7)
    rows = []
    for cls, n, p_struct, mu in (("folded", 300, 0.85, 0.3), ("disordered", 120, 0.35, 1.0)):
        for i in range(n):
            rows.append({"label": f"{cls[0]}{i}", "disorder_class": cls, "has": rng.random() < p_struct,
                         "chi2_raw": 10 ** rng.normal(mu, 0.3)})
    d = pd.DataFrame(rows)
    d["model"], d["clean"] = "bioemu", True
    d["rg_ratio"], d["nrmsd_raw"] = 1.0 + rng.normal(0, 0.05, len(d)), 0.03
    d["resolvability"] = np.where(d["chi2_raw"] <= 2, "raw fit", "strongly reweightable")
    mp = pd.DataFrame({"label": d["label"], "map_method": "SASBDB", "has_pdb_50": d["has"],
                       "has_pdb_any_overlap": d["has"]})
    out = coverage.structure_split(d.drop(columns="has"), mp)
    assert out["pooled_tests"]["chi2_raw_mannwhitney_p"] < 1e-6
    assert out["stratified_by_class"]["chi2_raw"]["p_two_sided"] > 0.01
    assert out["pooled"]["pdb50"]["n"] + out["pooled"]["under50"]["n"] == len(d)
    adj = [out["by_class"][c]["tests"]["chi2_raw_holm_p"] for c in ("folded", "disordered")]
    raw = [out["by_class"][c]["tests"]["chi2_raw_mannwhitney_p"] for c in ("folded", "disordered")]
    assert all(a >= r for a, r in zip(adj, raw, strict=True))


def test_holm_adjustment():
    adj = coverage.holm({"a": 0.01, "b": 0.04, "c": 0.03})
    assert adj == {"a": 0.03, "c": 0.06, "b": 0.06}


def test_identity_bins_close_on_the_left():
    b = coverage.identity_bin(pd.Series([0.0, 0.299, 0.3, 0.5, 0.899, 0.9, 1.0]))
    assert b.tolist() == ["<30", "<30", "30-50", "50-90", "50-90", ">=90", ">=90"]


def _synthetic_closeness(seed=11):
    """Identity tracks the class (folded proteins have close training relatives), the class sets the error,
    identity itself has no effect within a class."""
    rng = np.random.default_rng(seed)
    rows, nn = [], []
    for cls, n, mu_id, mu in (("folded", 200, 0.8, 0.3), ("partly disordered", 80, 0.5, 0.6),
                              ("disordered", 100, 0.25, 1.0)):
        for i in range(n):
            lab = f"{cls[0]}{i}"
            ident = float(np.clip(rng.normal(mu_id, 0.2), 0, 1))
            rows.append({"label": lab, "disorder_class": cls, "chi2_raw": 10 ** rng.normal(mu, 0.3),
                         "rg_ratio": 1.0 + rng.normal(0, 0.05), "nrmsd_raw": 0.03 + rng.normal(0, 0.005)})
            sf = rng.choice(["present", "absent", "none"], p=[0.6, 0.2, 0.2])
            nn.append({"label": lab, "nearest_identity": ident, "nearest_identity_md": ident / 2,
                       "id_afdb_member": ident, "sf_group": sf, "sf_group_md": sf})
    d = pd.DataFrame(rows)
    d["model"], d["clean"] = "bioemu", True
    d["resolvability"] = np.where(d["chi2_raw"] <= 2, "raw fit", "strongly reweightable")
    return d, pd.DataFrame(nn)


def test_closeness_confounded_by_class_is_null_within_class():
    d, nn = _synthetic_closeness()
    out = coverage.training_closeness(d, nn)
    assert out["spearman_pooled"]["chi2_raw"]["rho"] < -0.3 and out["spearman_pooled"]["chi2_raw"]["p"] < 1e-6
    assert sum(out["pooled"][b]["n"] for b in coverage.ID_LABELS) == len(d)
    for cls in ("folded", "partly disordered", "disordered"):
        sp = out["spearman_by_class"][cls]["chi2_raw"]
        assert sp["holm_p"] >= sp["p"] and sp["holm_p"] > 0.01
    assert out["stratified_identity_ge50"]["chi2_raw"]["p_two_sided"] > 0.01
    assert set(out["sf_group"]["pooled"]) == {"present", "absent", "none"}


def test_category_evidence_splits_every_profile_in_or_out():
    d, nn = _synthetic_closeness()
    nn["mp_dis_cat"] = np.where(d["disorder_class"] == "disordered", "mostly disordered", "neither")
    folded = d["disorder_class"] == "folded"
    mp = pd.DataFrame({"label": d["label"], "map_method": "SASBDB", "has_pdb_50": folded,
                       "has_pdb_any_overlap": d["disorder_class"] != "disordered",
                       "gene3d_superfamilies": np.where(folded, "3.40.50.300;1.10.8.10", ""),
                       "org_group": "Human"})
    rk = pd.DataFrame({"sf": ["3.40.50.300", "1.10.8.10"], "rank": [1, 250]})
    out = coverage.category_evidence(d.assign(length=120), mp, nn, rk)
    cat = out["categories"]
    for c in cat.values():
        assert c["in"]["n"] + c["rest"]["n"] == len(d)
    assert cat["alpha/beta enzyme folds"]["in"]["n"] == 200
    assert cat["superfamily beyond the top 200"]["in"]["n"] == 200
    assert cat["no experimental structure"]["in"]["n"] == 100
    assert cat["construct over 500 residues"]["in"]["n"] == 0


def test_quality_checks_by_class_counts_every_entry_once_per_flag():
    ent = pd.DataFrame({"disorder_class": ["folded"] * 4 + ["partly disordered"] * 2 + ["disordered"] * 3,
                        "flag_aggregation": [False, True, False, False, False, False, False, False, False],
                        "flag_negative_I": [False, True, False, False, False, False, False, False, False],
                        "flag_bad_errors": [False] * 9,
                        "flag_oligomer": [False, False, False, False, False, True, True, True, False],
                        "guinier_valid": [True, True, False, True, True, True, True, True, True]})
    q = coverage.quality_checks_by_class(ent)
    f, d = q["by_class"]["folded"], q["by_class"]["disordered"]
    assert (f["n"], f["pass"], d["n"], d["pass"]) == (4, 2, 3, 1)
    assert f["failed_by_check"]["aggregation upturn"] == 1
    assert f["failed_by_check"]["no valid Guinier region"] == 1
    assert d["failed_by_check"]["mass above a monomer"] == 2
    assert q["fisher_disordered_vs_folded"]["table_pass_fail"] == [[1, 2], [2, 2]]

"""The score stage on a synthetic archive laid out like the PeptoneBench predictions."""
import json

import numpy as np
import pandas as pd
import pytest


def sphere_curve(q, R):
    x = q * R
    return (3 * (np.sin(x) - x * np.cos(x)) / x**3) ** 2


@pytest.fixture
def synthetic_root(tmp_path, monkeypatch):
    from bsc import config as C
    data, results = tmp_path / "data", tmp_path / "results"
    saxs_dir = data / "PeptoneDB-SAXS" / "sasbdb-clean_data"
    pred = data / "Predictions" / "PeptoneDB-SAXS" / "bioemu"
    saxs_dir.mkdir(parents=True)
    pred.mkdir(parents=True)
    rng = np.random.default_rng(0)
    q = np.linspace(0.01, 0.3, 120)
    rows = []
    for i, (label, R_true, R_model, dis) in enumerate(
            [("SASDX01", 20.0, 20.0, 0.1), ("SASDX02", 30.0, 22.0, 0.7), ("SASDX03", 25.0, 25.0, 0.4)]):
        # the "experiment" is itself an ensemble average, as a real solution measurement is
        truth = np.mean([sphere_curve(q, r) for r in R_true + rng.normal(0, 1.5, 200)], axis=0)
        sigma = 0.02 * truth + 1e-4
        exp = truth + rng.normal(0, sigma)
        pd.DataFrame({"q": q, "I": exp, "s": sigma}).to_csv(saxs_dir / f"{label}-bift.dat", sep="\t",
                                                                index=False, header=["# q", "I(q)", "sigma"])
        radii = R_model + rng.normal(0, 1.5, 60)
        curves = pd.DataFrame(np.array([sphere_curve(q, r) for r in radii]), columns=[f"{v:.5f}" for v in q])
        curves.to_csv(pred / f"Pepsi-{label}.csv")
        rows.append({"label": label, "sequence": "A" * (50 + 30 * i), "length": 50 + 30 * i, "pH": 7.0,
                     "mean_gscore_adopt2": dis, "gscores_adopt2": "[]"})
    pd.DataFrame(rows).to_csv(data / "PeptoneDB-SAXS" / "PeptoneDB-SAXS.csv", index=False)
    monkeypatch.setattr(C, "DATA", data)
    monkeypatch.setattr(C, "RESULTS", results)
    monkeypatch.setattr(C, "SAXS_TABLE", data / "PeptoneDB-SAXS" / "PeptoneDB-SAXS.csv")
    monkeypatch.setattr(C, "SAXS_DIR", saxs_dir)
    monkeypatch.setattr(C, "PRED_DIR", pred.parent)
    monkeypatch.setattr(C, "N_JOBS", 1)
    monkeypatch.setattr(C, "THETA_GRID", 10.0 ** np.linspace(-1, 5, 9))
    return tmp_path


def test_score_stage_end_to_end(synthetic_root):
    from bsc import score
    score.main(["bioemu"])
    ent = pd.read_csv(synthetic_root / "results" / "entries.csv")
    sc = pd.read_csv(synthetic_root / "results" / "scores.csv")
    assert list(ent["disorder_class"]) == ["folded", "disordered", "partly disordered"]
    assert ent["guinier_valid"].all()
    s = sc.set_index("label")
    # the correct model fits; the mis-sized one does not, and reweighting cannot rescue it fully
    assert s.loc["SASDX01", "chi2_raw"] < 3
    assert s.loc["SASDX02", "chi2_raw"] > 20
    assert s.loc["SASDX02", "chi2_best"] > s.loc["SASDX01", "chi2_best"]
    assert (s["phi_at_best"] <= 1.0).all() and (s["n_conformers"] == 60).all()
    paths = json.load(open(synthetic_root / "results" / "paths" / "bioemu.json"))
    assert {p["label"] for p in paths} == set(s.index)
    # the ensemble Rg tracks the model radius, and the experimental Rg the true one
    assert abs(s.loc["SASDX02", "rg_model"] - 22 * np.sqrt(3 / 5)) < 1.5
    assert abs(ent.set_index("label").loc["SASDX02", "rg_exp"] - 30 * np.sqrt(3 / 5)) < 1.5


def test_analyse_stage_on_synthetic_scores(tmp_path, monkeypatch):
    """The analysis stage on a hand-made scores table with a known structure."""
    from bsc import analyse
    from bsc import config as C
    rng = np.random.default_rng(3)
    results = tmp_path / "results"
    results.mkdir()
    monkeypatch.setattr(C, "RESULTS", results)
    n = 90
    cls = np.repeat(["folded", "partly disordered", "disordered"], n // 3)
    length = rng.integers(60, 400, n)
    ent = pd.DataFrame({"label": [f"S{i:03d}" for i in range(n)], "length": length, "ph": 7.0,
                        "disorder_mean": np.select([cls == "folded", cls == "disordered"], [0.1, 0.8], 0.4),
                        "disorder_class": cls, "length_bin": "x", "n_q": 200, "q_min": 0.01, "q_max": 0.3,
                        "rg_exp": 20.0, "rg_exp_err": 0.3, "guinier_points": 30, "guinier_valid": True,
                        "upturn": 0.0, "flag_aggregation": False, "flag_negative_I": False,
                        "flag_bad_errors": False, "flag_oligomer": False, "mw_ratio": 1.0,
                        "sasbdb_title": "t", "sequence": "A"})
    ent.loc[0, "flag_aggregation"] = True
    # a flagged entry with no computed curves (no rows in scores.csv) must not count as clean
    extra = ent.iloc[[1]].assign(label="S999", flag_bad_errors=True)
    pd.concat([ent, extra]).to_csv(results / "entries.csv", index=False)
    # disordered worse, folded chi2 rising with length, model too compact for disordered
    base = np.where(cls == "disordered", 8.0, np.where(cls == "folded", 1.5, 3.0))
    chi2_raw = base * np.exp(rng.normal(0, 0.3, n)) * np.where(cls == "folded", length / 150, 1.0)
    rows = []
    for model, factor in (("bioemu", 1.0), ("alphafold2", 1.6)):
        c2 = chi2_raw * factor
        phi2 = np.where(c2 <= 2, 1.0, np.clip(1.5 / c2, 0.02, 0.9))
        phi2[(cls == "disordered") & (np.arange(n) % 7 == 0)] = np.nan
        rows.append(pd.DataFrame({"model": model, "label": ent["label"], "n_conformers": 100, "n_dropped": 0,
                                  "chi2_raw": c2, "chi2_raw_noback": c2 * 1.1,
                                  "rg_model": np.where(cls == "disordered", 16.0, 20.0),
                                  "rg_conformer_median": 20.0, "rg_conformer_iqr": 1.0,
                                  "chi2_best": np.minimum(c2, 1.1), "phi_at_best": phi2 * 0.5,
                                  "phi_at_chi2_1": phi2 * 0.8, "phi_at_chi2_2": phi2,
                                  "chi2_at_phi_0.5": np.maximum(c2 / 2, 1.0), "chi2_at_phi_0.1": 1.2,
                                  "chi2_at_phi_bench": 1.05}))
    pd.concat(rows).to_csv(results / "scores.csv", index=False)
    analyse.main()
    A = json.load(open(results / "analysis.json"))
    assert A["n_clean"] == n - 1
    assert A["n_clean"] == int(analyse.clean_mask(pd.read_csv(results / "entries.csv")).sum())
    E = A["expectations"]
    assert E["E1"]["met"] and E["E2"]["met"] and E["E4"]["met"]
    assert set(A["bioemu"]["resolvability_counts"]) <= set(analyse.KINDS)
    assert A["model_comparison"][0]["model"] == "alphafold2"
    assert A["model_comparison"][0]["median_log10_chi2_ratio_vs_bioemu"] > 0
    lo, hi = A["model_comparison"][0]["ci95_median_log10_chi2_ratio"]
    assert lo <= A["model_comparison"][0]["median_log10_chi2_ratio_vs_bioemu"] <= hi
    cand = pd.read_csv(results / "candidates.csv")
    assert cand["priority"].is_monotonic_decreasing
    assert "S000" not in set(cand["label"])          # the flagged entry is excluded
    assert (cand["resolvability"].iloc[0] in {analyse.KIND_STRONG, analyse.KIND_MODEST})
    # the additions after the first run: intervals cover the point estimates, the kind counts at the
    # pre-specified threshold equal the headline counts, and the examples are one entry per class
    R = A["rg_robustness"]
    lo, hi = R["curve_rg_ratio_by_class"]["disordered"]["ci95"]
    assert lo <= E["E4"]["median_rg_ratio_disordered"] <= hi
    H = A["headline_robustness"]
    lo, hi = H["share_raw_fit"]["ci95"]
    assert lo <= A["bioemu"]["share_raw_fit_chi2_2"] <= hi
    assert H["kind_counts_by_phi_threshold"]["0.5"] == {k: A["bioemu"]["resolvability_counts"].get(k, 0)
                                                        for k in analyse.KINDS}
    assert set(A["examples"]) == {"folded", "partly disordered", "disordered"}
    rep = pd.read_csv(results / "resolvability.csv").query("model == 'bioemu'").set_index("label")
    assert rep.loc[A["representative_path"], "resolvability"] == analyse.KIND_STRONG
    assert (results / "clusters.csv").exists()


def test_report_builds_from_analysis(tmp_path, monkeypatch):
    """Figures and the campaign proposal build from the synthetic analysis output."""
    from bsc import analyse, report
    from bsc import config as C
    test_analyse_stage_on_synthetic_scores(tmp_path, monkeypatch)
    results = tmp_path / "results"
    (results / "paths").mkdir()
    sc = pd.read_csv(results / "scores.csv")
    paths = [{"label": lab, "theta": th, "chi2": max(1.0, c2 / (1 + 10 / th)), "phi": min(1.0, th / 100)}
             for lab, c2 in zip(sc["label"], sc["chi2_raw"], strict=True) for th in (1000, 100, 10, 1, 0.1)]
    json.dump(paths, open(results / "paths" / "bioemu.json", "w"))
    docs = tmp_path / "docs"
    docs.mkdir()
    monkeypatch.setattr(C, "REPO", tmp_path)
    monkeypatch.setattr(report, "OUT", results / "figures")
    analyse.main()
    report.main()
    for f in ("fig_error_map", "fig_rg", "fig_rg_robustness", "fig_resolvability", "fig_models"):
        assert (results / "figures" / f"{f}.png").stat().st_size > 10_000
    # without DENSS maps and the validation stages the examples and prediction figures are left alone
    assert not (results / "figures" / "fig_examples.png").exists()
    assert not (results / "figures" / "fig_predict.png").exists()
    text = (docs / "CAMPAIGN.md").read_text()
    assert "## 3. Arm 1" in text and "## 4. Arm 2" in text and "| SASDB |" not in text
    arms = json.load(open(results / "campaign_arms.json"))
    assert arms["arm1"]["predictor_r2"] is None and "has not been run" in text
    assert all(r["ratio_model_over_law"] > 1.1 for r in arms["arm2"]["test"])
    assert all(abs(r["ratio_model_over_law"] - 1) <= 0.1 for r in arms["arm2"]["controls"])
    A = json.load(open(results / "analysis.json"))
    assert f"{A['bioemu']['chi2_raw_median']:.1f}" in text

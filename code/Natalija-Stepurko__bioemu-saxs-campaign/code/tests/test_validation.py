"""The two validation experiments on toy pools where the answer is known."""
import numpy as np
import pandas as pd

from bsc import analyse, features, predict, select_eval


def toy_pool(n: int, rng: np.random.Generator, signal: bool) -> tuple[pd.DataFrame, list[str], np.ndarray]:
    """A pool whose raw error is (signal=True) a clean function of two features, or (False) noise."""
    length = rng.integers(50, 500, n)
    dis = rng.random(n)
    x_noise = rng.normal(size=(n, 3))
    if signal:
        log_chi2 = -0.3 + 1.5 * dis + 0.4 * (np.log10(length) - 2) + rng.normal(0, 0.05, n)
    else:
        log_chi2 = rng.normal(0.5, 0.5, n)
    chi2 = 10.0 ** log_chi2
    # the kind follows the error: large errors are strongly reweightable or not fit
    phi2 = np.where(chi2 <= 2, 1.0, np.clip(2.0 / chi2, 0.02, 0.95))
    phi2 = np.where(log_chi2 > 1.2, np.nan, phi2)
    phi1 = np.where(np.isfinite(phi2), phi2 * 0.6, np.nan)
    cls = np.where(dis < 0.2, "folded", np.where(dis < 0.6, "partly disordered", "disordered"))
    t = pd.DataFrame({"label": [f"T{i:03d}" for i in range(n)], "length": length.astype(float),
                      "disorder_mean": dis, "class_folded": (cls == "folded").astype(float),
                      "class_partly": (cls == "partly disordered").astype(float),
                      "class_disordered": (cls == "disordered").astype(float),
                      "noise_a": x_noise[:, 0], "noise_b": x_noise[:, 1], "noise_c": x_noise[:, 2],
                      "chi2_raw": chi2, "phi_at_chi2_1": phi1, "phi_at_chi2_2": phi2, "disorder_class": cls})
    t["kind"] = [analyse.kind(r) for _, r in t.iterrows()]
    cols = features.feature_columns(t)
    groups = np.arange(n)          # every entry its own group
    return t, cols, groups


def test_predictor_cv_is_near_zero_r2_on_pure_noise():
    t, cols, groups = toy_pool(300, np.random.default_rng(1), signal=False)
    out, oof = predict.cross_validate(t, cols, groups, n_repeats=2)
    for k in ("constant", "class_mean", "ridge", "hgb"):
        assert out["regression"][k]["r2"]["mean"] < 0.1
    for k in ("rate", "class_rate", "logistic", "hgb"):
        assert out["classification"][k]["auc"]["mean"] < 0.65
    assert len(oof) == 300 and set(oof["label"]) == set(t["label"])


def test_predictor_cv_recovers_a_planted_signal():
    t, cols, groups = toy_pool(300, np.random.default_rng(2), signal=True)
    out, _ = predict.cross_validate(t, cols, groups, n_repeats=2)
    assert out["regression"]["ridge"]["r2"]["mean"] > 0.8
    assert out["regression"]["hgb"]["r2"]["mean"] > 0.8
    assert out["classification"]["hgb"]["auc"]["mean"] > 0.9
    top = {r["feature"] for r in out["permutation_importance_hgb"][:2]}
    assert "disorder_mean" in top
    assert sum(c["n"] for c in out["calibration_hgb"]) == 300


def test_policy_evaluation_on_a_toy_pool_with_a_known_answer():
    """When the error is predictable, the predicted policies approach the ceiling and beat random;
    the ceiling (observed quantities) is the best of all."""
    t, cols, groups = toy_pool(300, np.random.default_rng(3), signal=True)
    out = select_eval.evaluate(t, cols, groups, n_repeats=2, top_k=(10,), n_random=50)
    pol = out["policies"]["10"]
    rnd = pol["random"]["strong_or_notfit"]["mean"]
    assert abs(rnd - out["pool_mean"]["strong_or_notfit"]) < 0.1
    assert pol["predicted error"]["strong_or_notfit"]["mean"] > rnd + 0.3
    assert pol["predicted priority"]["strong_or_notfit"]["mean"] > rnd + 0.3
    ceiling = pol["observed priority (ceiling)"]
    assert ceiling["strong_or_notfit"]["mean"] >= pol["predicted priority"]["strong_or_notfit"]["mean"] - 1e-9
    assert ceiling["phi_reduction_chi2_2"]["mean"] > rnd
    assert out["random_single_draw_sd"]["10"]["strong_or_notfit"] > 0
    # the two components of the yield add up to it, policy by policy
    for p in pol.values():
        both = p["strongly_reweightable"]["mean"] + p["target_not_reached"]["mean"]
        assert abs(both - p["strong_or_notfit"]["mean"]) < 1e-9


def test_priority_rule_matches_the_candidate_ranking():
    log_chi2 = np.array([1.0, 1.0, 0.0, -0.5])
    w = np.array([1.0, 0.6, 1.0, 1.0])
    tract = np.array([1.0, 1.0, 1.0, 1.0])
    p = select_eval.priority(log_chi2, w, tract)
    assert p.tolist() == [1.0, 0.6, 0.0, 0.0]


def test_features_use_nothing_from_the_experimental_curve():
    cols = features.feature_columns(pd.DataFrame(columns=[
        "label", "length", "rg_model", "rg_exp", "chi2_raw", "kind", "disorder_class", "phi_at_chi2_1"]))
    assert "rg_exp" in cols            # the helper only drops known targets and identifiers ...
    ent = pd.DataFrame({"label": ["A"], "sequence": ["MKKLLPTAAAGLLLLAAQPAMA" * 4]})
    d = pd.DataFrame({"model": ["bioemu"], "clean": [True], "label": ["A"], "length": [88],
                      "disorder_mean": [0.3], "disorder_class": ["partly disordered"], "rg_model": [20.0],
                      "rg_conformer_iqr": [2.0], "rg_conformer_median": [19.0], "chi2_raw": [3.0],
                      "phi_at_chi2_1": [0.5], "phi_at_chi2_2": [0.8], "kind": ["modestly reweightable"],
                      "rg_exp": [21.0]})
    table = features.build(d, ent)
    # ... and the builder never copies an experimental quantity into the feature table
    assert not any(c.startswith("rg_exp") or "chi2" in c for c in features.feature_columns(table))
    assert abs(table["rg_law_disordered"].iloc[0] - 1.927 * 88 ** 0.598) < 1e-9

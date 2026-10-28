"""Stage `select_eval` (DESIGN §7, D1): a retrospective evaluation of the campaign's selection rule.

The 399 clean profiles are treated as a pool. On each held-out fold (5-fold, grouped by sequence
cluster), three acquisition policies rank the held-out entries from quantities known before
measuring (features.build, fitted on the training folds) and the top k are "acquired". The objective
is what those profiles would have done to the ensemble: the mean observed phi-reduction (1 − phi at
the chi2 target, 1 when the target is never reached) and the high-discrepancy target yield, the share
of entries that are strongly reweightable or whose target fit is not reached along the path, with
each of those two components separately. Policies: random; largest predicted raw error; the
campaign's priority rule with predicted quantities in place of the observed ones. The same rule with
the observed quantities is the ceiling.

This is retrospective: a profile that already exists carries no new information for the model, and
the evaluation says only whether the rule, fed with predictions, ranks existing profiles better than
chance. Output: results/select_eval.json.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor

from bsc import analyse, features, predict
from bsc import config as C

TOP_K = (10, 20)
N_RANDOM = 200
KIND_WEIGHT = {analyse.KIND_STRONG: 1.0, analyse.KIND_MODEST: 0.6, analyse.KIND_NOTFIT: 0.3,
               analyse.KIND_RAW: 0.0}


def tractability(table: pd.DataFrame) -> np.ndarray:
    cls_w = np.where(table["disorder_class"] == "disordered", 0.6, 1.0)
    return np.where(table["length"] <= 350, 1.0, 0.5) * cls_w


def priority(log_chi2: np.ndarray, kind_weight: np.ndarray, tract: np.ndarray) -> np.ndarray:
    """The campaign's rule: log10 raw chi2 (clipped at 0) × reweighting-kind weight × tractability."""
    return np.clip(log_chi2, 0, None) * kind_weight * tract


def objectives(table: pd.DataFrame) -> pd.DataFrame:
    """Observed outcomes per entry: phi-reduction at chi2 = 1 and at chi2 = 2 (1 when never reached),
    the high-discrepancy indicator (strongly reweightable or target fit not reached) and its two
    components separately."""
    out = pd.DataFrame(index=table.index)
    for t in (1, 2):
        phi = table[f"phi_at_chi2_{t}"].to_numpy(float)
        out[f"phi_reduction_chi2_{t}"] = np.where(np.isfinite(phi), 1 - phi, 1.0)
    out["strong_or_notfit"] = table["kind"].isin(predict.POSITIVE_KINDS).astype(float)
    out["strongly_reweightable"] = (table["kind"] == analyse.KIND_STRONG).astype(float)
    out["target_not_reached"] = (table["kind"] == analyse.KIND_NOTFIT).astype(float)
    return out


def fit_predictors(Xtr: pd.DataFrame, tr: pd.DataFrame, seed: int) -> dict:
    """Predictors of the quantities the rule needs, from pre-acquisition features only."""
    reg = HistGradientBoostingRegressor(random_state=seed, **predict.HGB).fit(Xtr, np.log10(tr["chi2_raw"]))
    kinds = tr["kind"].to_numpy()
    clf = HistGradientBoostingClassifier(random_state=seed, **predict.HGB).fit(Xtr, kinds)
    return {"log_chi2": reg, "kind": clf}


def evaluate(table: pd.DataFrame, cols: list[str], groups: np.ndarray, n_repeats: int = predict.N_REPEATS,
             top_k: tuple[int, ...] = TOP_K, n_random: int = N_RANDOM) -> dict:
    obj = objectives(table)
    tract = tractability(table)
    policies = ("random", "predicted error", "predicted priority", "observed priority (ceiling)")
    res = {k: {p: {o: [] for o in obj.columns} for p in policies} for k in top_k}
    draw_sd = {k: {o: [] for o in obj.columns} for k in top_k}      # spread of a single random pick of k
    rng = np.random.default_rng(C.SEED)
    for rep in range(n_repeats):
        per_fold = {k: {p: {o: [] for o in obj.columns} for p in policies} for k in top_k}
        for tr_idx, te_idx in predict.folds(groups, rep):
            tr, te = table.iloc[tr_idx], table.iloc[te_idx]
            models = fit_predictors(tr[cols], tr, C.SEED + rep)
            pred_log = models["log_chi2"].predict(te[cols])
            proba = models["kind"].predict_proba(te[cols])
            w_pred = proba @ np.array([KIND_WEIGHT[c] for c in models["kind"].classes_])
            scores = {
                "predicted error": pred_log,
                "predicted priority": priority(pred_log, w_pred, tract[te_idx]),
                "observed priority (ceiling)": priority(np.log10(te["chi2_raw"].to_numpy()),
                                                        te["kind"].map(KIND_WEIGHT).to_numpy(),
                                                        tract[te_idx]),
            }
            o_te = obj.iloc[te_idx]
            for k in top_k:
                kk = min(k, len(te))
                for p, s in scores.items():
                    top = np.argsort(-s, kind="stable")[:kk]
                    for o in obj.columns:
                        per_fold[k][p][o].append(float(o_te[o].to_numpy()[top].mean()))
                draws = pd.DataFrame([o_te.iloc[rng.choice(len(te), kk, replace=False)].mean()
                                      for _ in range(n_random)])
                for o in obj.columns:
                    per_fold[k]["random"][o].append(float(draws[o].mean()))
                    draw_sd[k][o].append(float(draws[o].std(ddof=1)))
        for k in top_k:
            for p in policies:
                for o in obj.columns:
                    res[k][p][o].append(float(np.mean(per_fold[k][p][o])))
    out = {"n": int(len(table)), "n_groups": int(len(np.unique(groups))), "n_splits": predict.N_SPLITS,
           "n_repeats": n_repeats, "n_random_draws_per_fold": n_random, "top_k": list(top_k),
           "objectives": list(obj.columns), "pool_mean": {o: float(obj[o].mean()) for o in obj.columns},
           "policies": {}}
    for k in top_k:
        out["policies"][str(k)] = {}
        out.setdefault("random_single_draw_sd", {})[str(k)] = {o: float(np.mean(draw_sd[k][o]))
                                                                for o in obj.columns}
        for p in policies:
            rec = {}
            for o in obj.columns:
                v = np.asarray(res[k][p][o])
                r = np.asarray(res[k]["random"][o])
                rec[o] = {"mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
                          "minus_random_mean": float((v - r).mean()),
                          "minus_random_sd": float((v - r).std(ddof=1)) if len(v) > 1 else 0.0,
                          "per_repeat": [float(x) for x in v]}
            out["policies"][str(k)][p] = rec
    return out


def main() -> None:
    ent, d = analyse.load()
    d["kind"] = d.apply(analyse.kind, axis=1)
    table = features.build(d, ent)
    groups = predict.sequence_groups(table)
    out = evaluate(table, features.feature_columns(table), groups)
    json.dump(out, open(C.RESULTS / "select_eval.json", "w"), indent=1)
    for k, pol in out["policies"].items():
        line = "; ".join(f"{p}: {v['phi_reduction_chi2_1']['mean']:.2f} / {v['strong_or_notfit']['mean']:.2f}"
                         for p, v in pol.items())
        print(f"  top {k} (phi-reduction at chi2=1 / high-discrepancy share) {line}",
              flush=True)


if __name__ == "__main__":
    main()

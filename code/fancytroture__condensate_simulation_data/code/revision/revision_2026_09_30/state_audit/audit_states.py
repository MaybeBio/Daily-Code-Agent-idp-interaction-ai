#!/usr/bin/env python3
"""Audit saved VAMPNet ensembles without altering submitted inputs or figures."""
from __future__ import annotations

import csv
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
REPO = ROOT.parent
os.environ["JPCL_DATA_ROOT"] = str(ROOT)
sys.path.insert(0, str(REPO / "figures"))
from supplementary import shared as sh
import make_figure3_msm_story as original

CACHE = REPO / "source_data/supplementary/derived_model_ensembles"
VARIANTS = (
    "submitted_edge_mask", "unfiltered_display_union",
    "common_submatrix", "common_conditional", "common_plus_exit",
)


def write_table(name, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (OUT / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def distance(tl, th, pl, ph):
    return np.sqrt(np.sum(0.5 * (pl + ph)[..., :, None] * (th-tl)**2, axis=(-2, -1)))


def conditional(t):
    sums = t.sum(axis=-1, keepdims=True)
    if np.any(sums <= 0):
        raise ValueError("Empty common-state transition row; conditional distance undefined")
    return t / sums


def metric(name, tl, th, pl, ph, common, ml, mh):
    if name == "submitted_edge_mask":
        return distance(tl * ml, th * mh, pl, ph)
    if name == "unfiltered_display_union":
        return distance(tl, th, pl, ph)
    a, b = tl[..., common, :][..., common], th[..., common, :][..., common]
    p, q = pl[..., common], ph[..., common]
    p, q = p / p.sum(axis=-1, keepdims=True), q / q.sum(axis=-1, keepdims=True)
    if name == "common_conditional":
        a, b = conditional(a), conditional(b)
    elif name == "common_plus_exit":
        # Retain escape probability; exits include all original noncommon states.
        a = np.concatenate((a, np.clip(1-a.sum(axis=-1, keepdims=True), 0, 1)), axis=-1)
        b = np.concatenate((b, np.clip(1-b.sum(axis=-1, keepdims=True), 0, 1)), axis=-1)
    return distance(a, b, p, q)


def load_condition(system, family, role, states):
    path = CACHE / f"{system}_{family}_{sh.TEMPS[system][role]}_n5_lag50.npz"
    raw = dict(np.load(path, allow_pickle=False))
    assert raw["t"].shape[1:] == (5, 5)
    assert np.allclose(raw["t"].sum(axis=-1), 1, atol=1e-7)
    assert np.allclose(raw["mu"].sum(axis=-1), 1, atol=1e-5)
    ts, ps = [], []
    paths = [path]
    for t, p, attempt in zip(raw["t"], raw["mu"], raw["attempts"]):
        map_path = sh.model_dir(system, family, role) / f"state_map_att{int(attempt)}_raw_to_final.npy"
        mapping = np.load(map_path).astype(int)
        assert sorted(mapping.tolist()) == list(range(5)), str(map_path)
        inv = np.argsort(mapping)
        ts.append(t[np.ix_(inv, inv)])
        ps.append(p[inv])
        paths.append(map_path)
    local_t, local_p = np.array(ts), np.array(ps)
    rows = sh.state_rows(system, family, role)
    mapping = {int(r["local_state"]): int(r.get("display_state") or r["global_state"]) for r in rows}
    assert len(set(mapping.values())) == len(mapping), "Non-bijective display mapping"
    index = {s: i for i, s in enumerate(states)}
    t = np.zeros((len(ts), len(states), len(states)))
    p = np.zeros((len(ps), len(states)))
    for i, di in mapping.items():
        p[:, index[di]] = local_p[:, i]
        for j, dj in mapping.items():
            t[:, index[di], index[dj]] = local_t[:, i, j]
    mask = np.zeros((len(states), len(states)))
    source_t = np.zeros_like(mask)
    for r in original.transition_network_rows(system):
        if (r["family"], r["role"]) != (family, role):
            continue
        i, j = original.row_display_states(r)
        mask[index[i], index[j]] = 1
        source_t[index[i], index[j]] += float(r["transition_probability"])
    return dict(t=t, p=p, local_t=local_t, local_p=local_p, raw=raw,
                rows=rows, mapping=mapping, mask=mask, source_t=source_t,
                paths=paths, retained=set(mapping.values()))


def main():
    condition_rows, matrix_rows, aligned_rows, model_rows = [], [], [], []
    results, decomposition, population_rows, contrast_rows, manifest = [], [], [], [], []
    boot_by_condition = {}
    published = {(r["system"], r["family"]): float(r["weighted_transition_distance"])
                 for r in sh.read_csv(ROOT / "project_figure/Figure3_transition_distance_source.csv")}
    rng = np.random.default_rng(20260930)
    n_boot = 4000
    for system in sh.SYSTEMS:
        # Reproduce the exact original universe of displayed network state labels.
        states = original.transition_states(system)
        for family in sh.FAMILIES:
            pair = [load_condition(system, family, role, states) for role in sh.ROLES]
            common = np.array([s in pair[0]["retained"] & pair[1]["retained"] for s in states])
            for role, data in zip(sh.ROLES, pair):
                prefix = dict(system=system, family=family, role=role, temperature=sh.TEMPS[system][role])
                raw = data["raw"]
                mean_t, mean_p = data["local_t"].mean(0), data["local_p"].mean(0)
                source_mu = {int(r["local_state"]): float(r["mu"]) for r in data["rows"]}
                population_sum = sum(source_mu.values())
                condition_rows.append(dict(**prefix, n_models=len(raw["attempts"]),
                    lag_ns=float(raw["base_lag"] * raw["dt"]),
                    retained_local_states=" ".join(map(str, data["mapping"])),
                    displayed_population_sum=population_sum,
                    ensemble_retained_population_sum=data["p"].mean(0).sum(),
                    source_vs_ensemble_population_max_abs=max(abs(source_mu[i]-mean_p[i]) for i in source_mu),
                    max_mean_population=mean_p.max(), mean_model_max_population=data["local_p"].max(1).mean(),
                    mean_n_states_gt_1pct=(data["local_p"] > 0.01).sum(1).mean(),
                    n_mean_states_gt_1pct=int((mean_p > 0.01).sum()),
                    recorded_network_entries=int(data["mask"].sum()),
                    recorded_network_diagonals=int(np.trace(data["mask"])),
                    local_row_sum_min=mean_t.sum(1).min(), local_row_sum_max=mean_t.sum(1).max()))
                for k, attempt in enumerate(raw["attempts"]):
                    for i in range(5):
                        model_rows.append(dict(**prefix, attempt=int(attempt), local_state=i,
                            display_state=data["mapping"].get(i, ""), population=data["local_p"][k, i]))
                for i in range(5):
                    for j in range(5):
                        matrix_rows.append(dict(**prefix, source_local_state=i, target_local_state=j,
                            source_display_state=data["mapping"].get(i, ""),
                            target_display_state=data["mapping"].get(j, ""),
                            population=mean_p[i], transition_probability=mean_t[i, j],
                            model_sd=data["local_t"][:, i, j].std(ddof=1)))
                for i, si in enumerate(states):
                    for j, sj in enumerate(states):
                        present = si in data["retained"] and sj in data["retained"]
                        aligned_rows.append(dict(**prefix, source_display_state=si, target_display_state=sj,
                            present=present, population=data["p"][:, i].mean() if si in data["retained"] else "",
                            transition_probability=data["t"][:, i, j].mean() if present else "",
                            included_by_submitted_mask=bool(data["mask"][i, j])))
                for path in data["paths"] + [sh.STATE_SOURCE[system][family]]:
                    manifest.append(dict(path=str(path.relative_to(ROOT)), size_bytes=path.stat().st_size,
                        sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            tl, th = [d["t"].mean(0) for d in pair]
            pl, ph = [d["p"].mean(0) for d in pair]
            ml, mh = [d["mask"] for d in pair]
            boot = []
            for d in pair:
                n = len(d["t"])
                ids = rng.integers(0, n, (n_boot, n))
                boot.append((d["t"][ids].mean(1), d["p"][ids].mean(1)))
            for name in VARIANTS:
                point = float(metric(name, tl, th, pl, ph, common, ml, mh))
                draws = metric(name, boot[0][0], boot[1][0], boot[0][1], boot[1][1], common, ml, mh)
                low, high = np.quantile(draws, [0.025, 0.975])
                boot_by_condition[(system, family, name)] = (point, draws)
                results.append(dict(system=system, family=family, definition=name, distance=point,
                    model_bootstrap_low=low, model_bootstrap_high=high,
                    common_states=" ".join(str(s) for s, ok in zip(states, common) if ok),
                    common_mass_Tl=pl[common].sum(), common_mass_Th=ph[common].sum(),
                    submitted_point=published[(system, family)]))
                if name == "submitted_edge_mask":
                    assert abs(point-published[(system, family)]) < 1e-7, (system, family, point)
            source_p = [np.array([original.state_vector(system, family, role).get(s, 0) for s in states]) for role in sh.ROLES]
            source_distance = float(distance(pair[0]["source_t"], pair[1]["source_t"], *source_p))
            results.append(dict(system=system, family=family, definition="released_network_table_only",
                distance=source_distance, submitted_point=published[(system, family)]))
            for name, a, b in (("submitted_population_table", *source_p), ("retained_ensemble_population", pl, ph)):
                a, b = a/a.sum(), b/b.sum()
                m = (a+b)/2
                js2 = sum(np.sum(x[x>0]*np.log2(x[x>0]/m[x>0])) for x in (a,b))/2
                population_rows.append(dict(system=system, family=family, definition=name,
                    total_variation=0.5*np.abs(a-b).sum(), js_distance=np.sqrt(max(js2, 0)),
                    population_normalization="conditional on retained states"))
            for name, a, b in (("submitted_edge_mask", tl*ml, th*mh), ("unfiltered_display_union", tl, th)):
                contrib = 0.5*(pl+ph)[:,None]*(b-a)**2
                cc = common[:,None] & common[None,:]
                d2 = contrib.sum()
                decomposition.append(dict(system=system, family=family, definition=name,
                    squared_distance=d2, common_common_squared=contrib[cc].sum(),
                    involving_noncommon_squared=contrib[~cc].sum(),
                    fraction_involving_noncommon=contrib[~cc].sum()/d2))
    for system in sh.SYSTEMS:
        for name in VARIANTS:
            cp, cb = boot_by_condition[(system,"condensate",name)]
            fp, fb = boot_by_condition[(system,"free",name)]
            lo, hi = np.quantile(cb-fb,[0.025,0.975])
            contrast_rows.append(dict(system=system, definition=name, condensate_minus_free=cp-fp,
                model_bootstrap_low=lo, model_bootstrap_high=hi))
    for name, rows in (("conditions.csv", condition_rows), ("full_local_transition_matrices.csv", matrix_rows),
                       ("aligned_transition_matrices_masked.csv", aligned_rows), ("model_populations.csv", model_rows),
                       ("distance_comparison.csv", results), ("distance_contrasts.csv", contrast_rows),
                       ("population_distances.csv", population_rows), ("distance_decomposition.csv", decomposition)):
        write_table(name, rows)
    write_table("input_manifest.csv", list({r["path"]:r for r in manifest}.values()))
    (OUT/"run_metadata.json").write_text(json.dumps(dict(seed=20260930, n_model_bootstrap=n_boot,
        states="zero-based local and display IDs", uncertainty="model seeds only, not sampling CI",
        source="saved VAMPNet soft-count transition matrices", python=sys.version), indent=2)+"\n")
    print("Original six published points reproduced to <1e-7; all local matrices stochastic.")
    for row in contrast_rows:
        print(row)


if __name__ == "__main__":
    main()

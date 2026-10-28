#!/usr/bin/env python3
"""Read-only saved-VAMPNet audit; all new artifacts stay beside this script."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT / "project_figure/scripts/supplementary_figures"))
import shared as s
sys.path.insert(0, str(OUT.parent / "state_audit"))
import audit_states as audit

CACHE = ROOT.parent / 'source_data/supplementary/derived_model_ensembles'
UNIT = 1000


def write_csv(name, rows):
    if not rows:
        return
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        writer.writeheader()
        writer.writerows(rows)


def transition(count):
    count = (count + count.swapaxes(-1, -2)) * 0.5
    den = count.sum(axis=-1, keepdims=True)
    result = count / np.maximum(den, 1e-12)
    result += np.eye(count.shape[-1]) * (den[..., 0] <= 1e-12)[..., :, None]
    return result


def condition(system, family, role):
    return dict(system=system, family=family, role=role, temperature=s.TEMPS[system][role])


def audit_condition(system, family, role, stride):
    tag = f"{system}_{family}_{s.TEMPS[system][role]}"
    base = condition(system, family, role)
    directory = s.model_dir(system, family, role)
    inter = directory / s.intermediate_name(system)
    source = CACHE / f"{tag}_n5_lag50.npz"
    with np.load(source) as z:
        ensemble = {k: z[k] for k in z.files}
    attempts = ensemble["attempts"]
    lag, dt = int(ensemble["base_lag"]), float(ensemble["dt"])
    assert lag == 200 and np.isclose(dt, 0.01)
    with np.load(inter / "combined_meta.npz") as z:
        lengths = z["lengths"].astype(int)
    n, nt, nb = len(attempts), len(lengths), int(np.ceil(max(lengths) / UNIT))
    counts = np.zeros((n, nt, nb, 5, 5))
    masses = np.zeros((n, nt, nb, 5))
    frame_counts = np.zeros((nt, nb), dtype=int)
    pair_counts = np.zeros((nt, nb), dtype=int)
    sample_q = []
    inventory, validation = [], []
    for ai, att in enumerate(attempts):
        mp = directory / f"state_map_att{att}_raw_to_final.npy"
        mapping = np.load(mp).astype(int)
        assert sorted(mapping.tolist()) == list(range(5)), str(mp)
        perm = np.argsort(mapping)
        qdir = inter / f"qcache_combined_n5_att{att}_lag50"
        files = sorted(qdir.glob("traj*.npy"), key=s.numeric_path_key)
        assert [p.name for p in files] == [f"traj{i}.npy" for i in range(nt)], str(qdir)
        qs = []
        max_norm_error = 0.0
        for ti, (path, length) in enumerate(zip(files, lengths)):
            raw = np.load(path, mmap_mode="r")
            assert raw.shape == (length, 5), str(path)
            assert np.isfinite(raw).all() and raw.min() >= -1e-7, str(path)
            max_norm_error = max(max_norm_error, float(abs(raw.sum(axis=1) - 1).max()))
            q = np.asarray(raw[:, perm], dtype=np.float64)
            qs.append(q[::stride].astype(np.float32))
            for bi, start in enumerate(range(0, length, UNIT)):
                end = min(start + UNIT, length)
                stop = min(end, length - lag)
                masses[ai, ti, bi] = q[start:end].sum(axis=0)
                frame_counts[ti, bi] = end - start
                pair_counts[ti, bi] = max(0, stop - start)
                if stop > start:
                    counts[ai, ti, bi] = q[start:stop].T @ q[start + lag:stop + lag]
            inventory.append(dict(**base, attempt=int(att), trajectory=ti, frames=int(length),
                                  path=str(path.relative_to(ROOT)), bytes=path.stat().st_size,
                                  mtime_ns=path.stat().st_mtime_ns,
                                  mapping=','.join(map(str, mapping)),
                                  mapping_sha256=hashlib.sha256(mp.read_bytes()).hexdigest()))
        sample_q.append(np.concatenate(qs))
        t = transition(counts[ai].sum(axis=(0, 1)))
        mu = masses[ai].sum(axis=(0, 1)) / lengths.sum()
        err_t = float(abs(t - ensemble["t"][ai][np.ix_(perm, perm)]).max())
        err_mu = float(abs(mu - ensemble["mu"][ai][perm]).max())
        assert err_t < 1e-8 and err_mu < 1e-8, (tag, att, err_t, err_mu)
        assert max_norm_error < 1e-5
        validation.append(dict(**base, attempt=int(att), max_t_cache_error=err_t,
                               max_mu_cache_error=err_mu, max_q_norm_error=max_norm_error))
    np.savez_compressed(OUT / f"{tag}_sufficient_statistics.npz", counts=counts, masses=masses,
                        frame_counts=frame_counts, pair_counts=pair_counts, lengths=lengths,
                        attempts=attempts, lag=lag, dt=dt, unit_frames=UNIT)
    q = np.stack(sample_q)
    hard = q.argmax(axis=-1)
    pair_rows, state_rows = [], []
    per_state = {i: [] for i in range(5)}
    for a in range(n):
        for b in range(a + 1, n):
            mass_a, mass_b = q[a].sum(axis=0, dtype=np.float64), q[b].sum(axis=0, dtype=np.float64)
            intersection = np.minimum(q[a], q[b]).sum(axis=0, dtype=np.float64)
            agree = hard[a] == hard[b]
            ha, hb = np.bincount(hard[a], minlength=5), np.bincount(hard[b], minlength=5)
            hi = np.bincount(hard[a][agree], minlength=5)
            soft_dice = 2 * intersection / np.maximum(mass_a + mass_b, 1e-12)
            hard_dice = 2 * hi / np.maximum(ha + hb, 1)
            pair_rows.append(dict(**base, attempt_a=int(attempts[a]), attempt_b=int(attempts[b]),
                                  sample_frames=q.shape[1], frame_stride=stride,
                                  soft_overlap=float(intersection.sum() / q.shape[1]),
                                  hard_overlap=float(agree.mean()),
                                  soft_dot=float(np.einsum('ij,ij->', q[a], q[b], dtype=np.float64) / q.shape[1])))
            for state in range(5):
                per_state[state].append((soft_dice[state], hard_dice[state]))
    full_mu = masses.sum(axis=(1, 2)) / lengths.sum()
    display = s.display_state_mapping(system, family, role)
    for state in range(5):
        vals = np.asarray(per_state[state])
        state_rows.append(dict(**base, local_state=state, display_state=display.get(state, "unmapped"),
                               mean_population=float(full_mu[:, state].mean()),
                               population_min=float(full_mu[:, state].min()),
                               population_max=float(full_mu[:, state].max()),
                               major_state=bool(full_mu[:, state].mean() >= 0.1),
                               soft_dice_mean=float(vals[:, 0].mean()),
                               soft_dice_p05=float(np.quantile(vals[:, 0], .05)),
                               hard_dice_mean=float(vals[:, 1].mean()),
                               hard_dice_p05=float(np.quantile(vals[:, 1], .05))))
    summary = dict(**base, models=n, trajectories=nt, frames=int(lengths.sum()),
                   sampled_frames=q.shape[1], frame_stride=stride,
                   soft_overlap_mean=float(np.mean([v["soft_overlap"] for v in pair_rows])),
                   hard_overlap_mean=float(np.mean([v["hard_overlap"] for v in pair_rows])),
                   soft_overlap_p05=float(np.quantile([v["soft_overlap"] for v in pair_rows], .05)),
                   hard_overlap_p05=float(np.quantile([v["hard_overlap"] for v in pair_rows], .05)))
    print(json.dumps(summary), flush=True)
    return inventory, validation, pair_rows, state_rows, summary


def units(system, family, role, mode, block):
    path = OUT / f"{system}_{family}_{s.TEMPS[system][role]}_sufficient_statistics.npz"
    with np.load(path) as z:
        c, m, lengths = z["counts"], z["masses"], z["lengths"]
    if mode == "trajectory":
        return c.sum(axis=2), m.sum(axis=2), []
    # Selection windows have different starts: synchronize by absolute slab time.
    absolute = OUT / f"{system}_{family}_{s.TEMPS[system][role]}_absolute_statistics.npz"
    with np.load(absolute) as z:
        c, m, excluded = z["counts"], z["masses"], z["excluded_trajectories"].tolist()
    factor = block // UNIT
    assert block % UNIT == 0 and c.shape[1] % factor == 0
    c = c.reshape(c.shape[0], -1, factor, 5, 5).sum(axis=2)
    m = m.reshape(m.shape[0], -1, factor, 5).sum(axis=2)
    return c, m, excluded


def draw_models(c, m, nboot, rng):
    nu = c.shape[1]
    weights = rng.multinomial(nu, np.full(nu, 1 / nu), size=nboot)
    bc = np.einsum('bu,muij->bmij', weights, c, optimize=True)
    bm = np.einsum('bu,mui->bmi', weights, m, optimize=True)
    # Draw trajectories once per replicate, shared across all fixed model seeds.
    t = transition(bc).mean(axis=1)
    p = (bm / bm.sum(axis=-1, keepdims=True)).mean(axis=1)
    pc, pm = c.sum(axis=1), m.sum(axis=1)
    return np.concatenate([transition(pc).mean(axis=0)[None], t]), np.concatenate([
        (pm / pm.sum(axis=-1, keepdims=True)).mean(axis=0)[None], p])


def mapped(t, p, mapping, states):
    out_t = np.zeros((len(t), len(states), len(states)))
    out_p = np.zeros((len(t), len(states)))
    for i, display_i in mapping.items():
        if display_i not in states:
            continue
        ii = states.index(display_i)
        out_p[:, ii] += p[:, i]
        for j, display_j in mapping.items():
            if display_j in states:
                out_t[:, ii, states.index(display_j)] = t[:, i, j]
    return out_t, out_p


def bootstrap(nboot, seed):
    rng = np.random.default_rng(seed)
    rows, diagnostics, contrast, draws = [], [], [], {}
    for block in (1000, 5000, 10000):
        for system in s.SYSTEMS:
            for family in s.FAMILIES:
                mode = "trajectory" if family == "free" else "absolute_slab_time_block"
                arrays, maps, excluded, n_units = [], [], [], []
                for role in s.ROLES:
                    c, m, omitted = units(system, family, role, mode, block)
                    arrays.append(draw_models(c, m, nboot, rng))
                    maps.append(s.display_state_mapping(system, family, role))
                    excluded.append(omitted)
                    n_units.append(c.shape[1])
                states = audit.original.transition_states(system)
                reference = [audit.load_condition(system, family, role, states) for role in s.ROLES]
                common = np.asarray([v in set(maps[0].values()) & set(maps[1].values()) for v in states])
                for metric in audit.VARIANTS:
                    tl, pl = mapped(*arrays[0], maps[0], states)
                    th, ph = mapped(*arrays[1], maps[1], states)
                    dist = audit.metric(metric, tl, th, pl, ph, common,
                                        reference[0]["mask"], reference[1]["mask"])
                    key = f"{system}_{family}_{block}_{metric}"
                    draws[key] = dist
                    rows.append(dict(system=system, family=family, metric=metric, states=','.join(map(str, states)),
                                     common_states=','.join(str(v) for v, ok in zip(states, common) if ok),
                                     method=mode, block_frames=block if family == "condensate" else 0,
                                     block_ns=block * .01 if family == "condensate" else 0,
                                     comparison_block_frames=block, units_Tl=n_units[0], units_Th=n_units[1],
                                     excluded_Tl=str(excluded[0]), excluded_Th=str(excluded[1]),
                                     estimate=float(dist[0]), bootstrap_mean=float(dist[1:].mean()),
                                     sampling_sd=float(dist[1:].std(ddof=1)),
                                     ci_low=float(np.quantile(dist[1:], .025)), ci_high=float(np.quantile(dist[1:], .975)),
                                     n_boot=nboot, fixed_models=20))
                    for role, tt, pp in (("Tl", tl, pl), ("Th", th, ph)):
                        diagnostics.append(dict(system=system, family=family, role=role, metric=metric,
                                                comparison_block_frames=block,
                                                retained_population=float(pp[0].sum()),
                                                common_population=float(pp[0, common].sum()),
                                                common_row_mass_min=float(tt[0][np.ix_(common, common)].sum(axis=1).min()),
                                                common_row_mass_max=float(tt[0][np.ix_(common, common)].sum(axis=1).max()),
                                                min_retained_row_mass=float(tt[0].sum(axis=1).min()),
                                                max_retained_row_mass=float(tt[0].sum(axis=1).max())))
            for metric in audit.VARIANTS:
                a = draws[f"{system}_condensate_{block}_{metric}"] - draws[f"{system}_free_{block}_{metric}"]
                contrast.append(dict(system=system, metric=metric, block_frames=block, estimate=float(a[0]),
                                     ci_low=float(np.quantile(a[1:], .025)), ci_high=float(np.quantile(a[1:], .975)),
                                     bootstrap_fraction_gt_zero=float(np.mean(a[1:] > 0)), n_boot=nboot))
    write_csv("sampling_distances.csv", rows)
    write_csv("sampling_contrasts.csv", contrast)
    write_csv("distance_support_diagnostics.csv", diagnostics)
    np.savez_compressed(OUT / "bootstrap_draws.npz", **draws)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stride", type=int, default=10)
    parser.add_argument("--nboot", type=int, default=3000)
    parser.add_argument("--seed", type=int, default=20260930)
    parser.add_argument("--bootstrap-only", action="store_true")
    args = parser.parse_args()
    assert args.stride >= 1 and args.nboot >= 2
    if not args.bootstrap_only:
        tables = [[] for _ in range(5)]
        missing = []
        for system in s.SYSTEMS:
            for family in s.FAMILIES:
                for role in s.ROLES:
                    try:
                        result = audit_condition(system, family, role, args.stride)
                    except FileNotFoundError as err:
                        missing.append(dict(**condition(system, family, role), missing_path=str(err.filename)))
                        continue
                    for table, values in zip(tables, result):
                        table.extend(values if isinstance(values, list) else [values])
                    for name, table in zip(("inventory.csv", "cache_validation.csv", "model_pairs.csv",
                                             "state_agreement.csv", "condition_agreement.csv"), tables):
                        write_csv(name, table)
        (OUT / "missing_inputs.json").write_text(json.dumps(missing, indent=2) + "\n")
        if missing:
            raise RuntimeError("Missing inputs recorded; complete-condition tables retained")
    bootstrap(args.nboot, args.seed)
    (OUT / "run_metadata.json").write_text(json.dumps(dict(argv=sys.argv, seed=args.seed,
        nboot=args.nboot, agreement_stride=args.stride, numpy=np.__version__, python=sys.version,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        metric_script_sha256=hashlib.sha256(Path(audit.__file__).read_bytes()).hexdigest(),
        shared_sha256=hashlib.sha256(Path(s.__file__).read_bytes()).hexdigest()), indent=2) + "\n")


if __name__ == "__main__":
    main()

"""Paired snapshot summaries for transient 0.1-20 ns temperature controls."""
import csv
import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from extract_full20ns import CASES, ROOT, OUT, pure_function, sha, write_json

RG_SOURCE = ROOT / "fus/fus_protein_single_simulation/thermo_figure1_rg_ree_relaxation.py"
INTRA_SOURCE = ROOT / "ppt_5_11_figures/compute_intrachain_contact_relaxation.py"
PARTNER_SOURCE = ROOT / "ppt_5_11_figures/compute_interchain_partner_survival.py"
PILOT_SOURCE = ROOT / "revision_2026_10_08/FUS_slab_temperature_jump/analyze_pilot_acf.py"
STYLE_SOURCE = ROOT / "project_figure/scripts/make_figure2_dynamics_story.py"
acf_rg = pure_function(RG_SOURCE, "autocorr_fft")
acf_unbiased = pure_function(INTRA_SOURCE, "autocorr_fft")
partner_curve = pure_function(PILOT_SOURCE, "partner_curve")
direct_partner = pure_function(PARTNER_SOURCE, "segment_partner_survival")
WINDOWS = {"full": (.1, 20.), "early": (.1, 10.), "late": (10., 20.)}
METRICS = ("Rg", "N_intra", "N_inter", "N_partner", "partner")
LABELS = {"Rg": r"$R_g$ ACF", "N_intra": r"$N_{\mathrm{intra}}$ ACF",
          "N_inter": r"$N_{\mathrm{inter}}$ ACF", "N_partner": "Partner-number ACF",
          "partner": r"$S_{\mathrm{partner}}$"}
MEAN_LABELS = {"density": "Core density (g/mL)", "Rg": r"$R_g$ (nm)",
               "N_intra": r"$N_{\mathrm{intra}}$", "N_inter": r"$N_{\mathrm{inter}}$",
               "N_partner": "Partner number", "load_per_partner": "Contacts per partner",
               "core_cohort_fraction": "Initial cohort in core"}
DISPLAY = {"FUS": "FUS-LCD", "A1": "A1-LCD", "TDP43": "TDP-43-LCD"}
NBOOT = 4000
RNG = np.random.default_rng(20261009)
BOOT_COUNTS = np.stack([np.bincount(row, minlength=8) / 8
                        for row in RNG.integers(0, 8, size=(NBOOT, 8))])


def table(name, rows):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def interval(values):
    return np.quantile(values, [.025, .975], method="inverted_cdf")


def crossing(curves, threshold, dt=.01):
    curves = np.atleast_2d(curves)
    hit = curves <= threshold
    k = hit.argmax(1)
    has = hit.any(1)
    before = np.maximum(k - 1, 0)
    row = np.arange(len(curves))
    result = np.full(len(curves), np.inf)
    valid = has & (k > 0)
    y0, y1 = curves[row[valid], before[valid]], curves[row[valid], k[valid]]
    result[valid] = (before[valid] + (threshold - y0) / (y1 - y0)) * dt
    result[has & (k == 0)] = 0
    return result


def acf(values, metric, mode, maxlag):
    x = np.asarray(values, float).copy()
    if mode == "detrended":
        t = np.arange(len(x)) - (len(x) - 1) / 2
        x -= t[:, None] * (t @ x / (t @ t))[None]
    fn = acf_rg if metric == "Rg" and mode != "unbiased" else acf_unbiased
    assert np.all(x.var(0) > 1e-12), "Constant chain series: needs explicit handling"
    values = np.stack([fn(x[:, c])[:maxlag + 1] for c in range(x.shape[1])])
    assert np.isfinite(values).all()
    np.testing.assert_allclose(values[:, 0], 1., atol=1e-12)
    return values.mean(0)


def test_math():
    rng = np.random.default_rng(53)
    x = rng.normal(size=(80, 4))
    d = x - x.mean(0)
    for metric in ("Rg", "N_inter", "N_partner"):
        direct = np.array([((d[:80-k] * d[k:]).mean(0) / (d*d).mean(0)).mean() for k in range(31)])
        if metric == "Rg":
            direct *= (80 - np.arange(31)) / 80
        np.testing.assert_allclose(acf(x, metric, "original", 30), direct, atol=1e-12)
    contacts = rng.integers(0, 2, size=(80, 4, 10))
    contacts[:, :, 0] = 1
    direct = np.mean([direct_partner([set(np.flatnonzero(v)) for v in contacts[:, c]], 30)[0]
                      for c in range(4)], axis=0)
    np.testing.assert_allclose(partner_curve(contacts, 30), direct, atol=1e-12)
    np.testing.assert_allclose(crossing([[1., .4, .1]], .2, 1), [5/3])
    assert np.isinf(crossing([[1., .9]], .2)[0])
    assert interval([1, 2, np.inf])[1] == np.inf


def branch_summary(case, snap, temp):
    path = OUT / "observables" / f"{case}_snapshot{snap:02d}_{temp}K.npz"
    meta = json.loads(path.with_suffix(".json").read_text())
    assert sha(path) == meta["output_sha256"]
    with np.load(path) as loaded:
        data = dict(loaded)
    curve_sets, rows, trends, per_snapshot = {}, [], [], []
    time = data["time_ns"]
    mean_series = {m: data[m].mean(1) for m in METRICS if m != "partner"}
    mean_series["density"] = data["density_core_3_5_7"][:, 1]
    per_partner = np.divide(data["N_inter"], data["N_partner"],
        out=np.zeros_like(data["N_inter"], dtype=float), where=data["N_partner"] > 0)
    assert np.all(data["N_partner"] > 0), "Empty partner sets need separate reporting"
    mean_series["load_per_partner"] = per_partner.mean(1)
    mean_series["core_cohort_fraction"] = (abs(data["chain_center_z_nm"]) < 5).mean(1)
    for window, (begin, end) in WINDOWS.items():
        keep = (time > begin + 1e-8) & (time <= end + 1e-8)
        maxlag = int(keep.sum()) // 2
        lag = np.arange(maxlag + 1) * .01
        for metric in METRICS:
            modes = ["original"]
            if window == "full" and metric != "partner":
                modes += ["detrended"]
                if metric == "Rg":
                    modes += ["unbiased"]
            for mode in modes:
                curve = (partner_curve(data["partner_contacts"][keep], maxlag) if metric == "partner"
                         else acf(data[metric][keep], metric, mode, maxlag))
                curve_sets[window, metric, mode] = curve
                threshold = .5 if metric == "partner" else .2
                tau = crossing(curve, threshold)[0]
                per_snapshot.append(dict(case=case, snapshot=snap, temperature_K=temp, window=window,
                    metric=metric, mode=mode, threshold=threshold, crossing_ns=tau,
                    censored=not np.isfinite(tau), horizon_ns=lag[-1],
                    integral_0_2ns=np.trapz(curve[:201], lag[:201]),
                    integral_0_horizon_ns=np.trapz(curve, lag), value_at_2ns=curve[200],
                    value_at_horizon=curve[-1]))
        for metric, series in mean_series.items():
            rows.append(dict(case=case, snapshot=snap, temperature_K=temp, window=window,
                metric=metric, value=series[keep].mean(), n_chains=len(data["chain_ids"])))
    for start in np.arange(0, 20, .5):
        keep = (time > max(.1, start) + 1e-8) & (time <= start + .5 + 1e-8)
        for metric, series in mean_series.items():
            trends.append(dict(case=case, snapshot=snap, temperature_K=temp,
                               time_ns=start + .25, metric=metric, value=series[keep].mean()))
    return curve_sets, rows, trends, per_snapshot


def summarize_curves(all_curves):
    curves_rows, timescale_rows, contrast_rows = [], [], []
    plot_cache = {}
    for case, (_, low, high) in CASES.items():
        keys = sorted({key[2:] for key in all_curves if key[0] == case})
        for window, metric, mode in keys:
            values = [np.stack(all_curves[case, t, window, metric, mode]) for t in (low, high)]
            lag = np.arange(values[0].shape[1]) * .01
            boot = [BOOT_COUNTS @ v for v in values]
            mean = [v.mean(0) for v in values]
            threshold = .5 if metric == "partner" else .2
            for t, m, b in zip((low, high), mean, boot):
                ci = np.quantile(b, [.025, .975], axis=0)
                plot_cache[case, t, window, metric, mode] = (lag, m, ci)
                curves_rows.extend(dict(case=case, temperature_K=t, window=window, metric=metric,
                    mode=mode, lag_ns=l, mean=v, ci_low=a, ci_high=z)
                    for l, v, a, z in zip(lag, m, ci[0], ci[1]))
                tau = crossing(m, threshold)[0]
                tau_boot = crossing(b, threshold)
                lo, hi = interval(tau_boot)
                timescale_rows.append(dict(case=case, temperature_K=t, window=window, metric=metric,
                    mode=mode, threshold=threshold, crossing_ns=tau, ci_low=lo, ci_high=hi,
                    boot_censored_fraction=np.mean(~np.isfinite(tau_boot)), horizon_ns=lag[-1]))
            measures = [("threshold_time", threshold)]
            if window == "full" and mode == "original":
                measures += [("threshold_sensitivity", th) for th in ((.4, .6) if metric == "partner" else (.1, .3))]
            measures += [("integral_0_2ns", None), ("integral_0_horizon", None),
                         ("value_at_2ns", None), ("value_at_horizon", None)]
            if lag[-1] >= 5:
                measures += [("integral_0_5ns", None), ("value_at_5ns", None)]
            for measure, th in measures:
                def compute(a):
                    a = np.atleast_2d(a)
                    if th is not None:
                        return crossing(a, th)
                    if measure.startswith("integral"):
                        end = 201 if measure == "integral_0_2ns" else 501 if measure == "integral_0_5ns" else len(lag)
                        return np.trapz(a[:, :end], lag[:end], axis=1)
                    index = 200 if measure == "value_at_2ns" else 500 if measure == "value_at_5ns" else -1
                    return a[:, index]
                point = [compute(m)[0] for m in mean]
                draws = [compute(b) for b in boot]
                finite_point = np.isfinite(point).all()
                valid = np.isfinite(draws[0]) & np.isfinite(draws[1])
                ratio_ok = valid.all() and np.all(draws[0] > 0)
                difference_ci = interval(draws[1] - draws[0]) if valid.all() else [None, None]
                ratio_ci = interval(draws[1] / draws[0]) if ratio_ok else [None, None]
                contrast_rows.append(dict(case=case, window=window, metric=metric, mode=mode,
                    measure=measure, threshold=th, low_value=point[0], high_value=point[1],
                    high_minus_low=point[1] - point[0] if finite_point else None,
                    difference_ci_low=difference_ci[0], difference_ci_high=difference_ci[1],
                    high_over_low=point[1] / point[0] if finite_point and point[0] > 0 else None,
                    ratio_ci_low=ratio_ci[0], ratio_ci_high=ratio_ci[1],
                    bootstrap_both_uncensored_fraction=valid.mean(),
                    ci_note="All 4000 paired draws used" if ratio_ok else "Ratio CI not reported: censored draws or nonpositive denominator",
                    horizon_ns=lag[-1]))
    return curves_rows, timescale_rows, contrast_rows, plot_cache


def summarize_means(rows):
    summary = []
    for case, (_, low, high) in CASES.items():
        for window in WINDOWS:
            for metric in MEAN_LABELS:
                values = [np.array([next(r["value"] for r in rows if r["case"] == case
                    and r["snapshot"] == s and r["temperature_K"] == t and r["window"] == window
                    and r["metric"] == metric) for s in range(1, 9)]) for t in (low, high)]
                b = [BOOT_COUNTS @ v for v in values]
                diff_ci = interval(b[1] - b[0])
                ratio_ci = interval(b[1] / b[0])
                summary.append(dict(case=case, window=window, metric=metric,
                    low_mean=values[0].mean(), high_mean=values[1].mean(),
                    high_minus_low=(values[1] - values[0]).mean(),
                    difference_ci_low=diff_ci[0], difference_ci_high=diff_ci[1],
                    high_over_low=values[1].mean() / values[0].mean(),
                    ratio_ci_low=ratio_ci[0], ratio_ci_high=ratio_ci[1],
                    pairs_high_gt_low=int(np.sum(values[1] > values[0]))))
    return summary


def style_module():
    spec = importlib.util.spec_from_file_location("maintext_style", STYLE_SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.setup_style()
    return module


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=200)
    plt.close(fig)


def plot_curves(cache, style, metrics, name, mode="original"):
    fig, axes = plt.subplots(len(metrics), 3, figsize=(15.5, 3.6 * len(metrics)), squeeze=False)
    for row, metric in enumerate(metrics):
        row_min = min(float(cache[case, t, "full", metric, mode][2].min())
                      for case, (_, low, high) in CASES.items() for t in (low, high))
        for col, (case, (_, low, high)) in enumerate(CASES.items()):
            ax = axes[row, col]
            for t, role in ((low, "Tl"), (high, "Th")):
                lag, mean, ci = cache[case, t, "full", metric, mode]
                ax.plot(lag, mean, color=style.TEMP_COLORS[role], lw=2, label=style.temp_label(role))
                ax.fill_between(lag, ci[0], ci[1], color=style.TEMP_COLORS[role], alpha=.16, linewidth=0)
            ax.axhline(.5 if metric == "partner" else .2, color=".5", ls=":", lw=.9)
            ax.set(xlim=(0, 10), ylabel=LABELS[metric], xlabel="Lag time (ns)")
            ax.set_ylim((0, 1.04) if metric == "partner" else (min(-.25, row_min-.02), 1.04))
            ax.set_xticks(np.arange(0, 11, 2))
            if row == 0:
                ax.set_title(DISPLAY[case], pad=12)
            if row == 0 and col == 0:
                ax.legend(ncol=2, loc="upper right", columnspacing=.8, handlelength=1.5)
            style.box_axis(ax)
            ax.text(-.17, 1.045, chr(65+row*3+col), transform=ax.transAxes, fontsize=20, fontweight="bold")
    fig.subplots_adjust(left=.08, right=.985, bottom=.09 if len(metrics)>2 else .12,
                        top=.92 if len(metrics)==2 else .94, wspace=.32, hspace=.42)
    save(fig, name)


def plot_trends(trends, style):
    metrics = ("density", "Rg", "N_intra", "N_inter", "N_partner", "core_cohort_fraction")
    fig, axes = plt.subplots(len(metrics), 3, figsize=(15.5, 19))
    for row, metric in enumerate(metrics):
        for col, (case, (_, low, high)) in enumerate(CASES.items()):
            ax = axes[row, col]
            x = np.arange(.25, 20, .5)
            for t, role in ((low, "Tl"), (high, "Th")):
                vals = np.array([[r["value"] for r in trends if r["case"]==case and r["snapshot"]==s
                    and r["temperature_K"]==t and r["metric"]==metric] for s in range(1, 9)])
                ci = np.quantile(BOOT_COUNTS @ vals, [.025, .975], axis=0)
                ax.plot(x, vals.mean(0), color=style.TEMP_COLORS[role], lw=2, label=style.temp_label(role))
                ax.fill_between(x, ci[0], ci[1], color=style.TEMP_COLORS[role], alpha=.16)
            ax.set(xlim=(0, 20), ylabel=MEAN_LABELS[metric], xlabel="Time after branching (ns)")
            if row == 0:
                ax.set_title(DISPLAY[case], pad=12)
            if row == 0 and col == 0:
                ax.legend(ncol=2, loc="best", handlelength=1.4, columnspacing=.7)
            style.box_axis(ax)
    fig.subplots_adjust(left=.09, right=.985, top=.96, bottom=.05, wspace=.36, hspace=.5)
    save(fig, "Fig_control_time_evolution")


def plot_ratios(rows, name, measure, style):
    metrics = ("Rg", "N_intra", "N_inter", "N_partner", "partner")
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.7))
    displayed = [r for r in rows if r["window"]=="full" and r["mode"]=="original"
                 and r["measure"]==measure and r["metric"] in metrics]
    upper = max(1.1, 1.07 * max(r["ratio_ci_high"] or r["high_over_low"] or 0 for r in displayed))
    labels = [r"$R_g$", r"$N_{\mathrm{intra}}$", r"$N_{\mathrm{inter}}$", "Partner\nnumber", r"$S_{\mathrm{partner}}$"]
    for ax, case in zip(axes, CASES):
        ax.axhline(1, color=".5", ls=":", lw=1)
        for i, metric in enumerate(metrics):
            r = next(r for r in rows if r["case"]==case and r["window"]=="full"
                     and r["mode"]=="original" and r["metric"]==metric and r["measure"]==measure)
            point, lo, hi = r["high_over_low"], r["ratio_ci_low"], r["ratio_ci_high"]
            if point is not None:
                ax.plot(i, point, "o", color="#4a4a4a", ms=7)
                if lo is not None and hi is not None:
                    ax.vlines(i, lo, hi, color="#4a4a4a", lw=1.5)
                    ax.hlines([lo, hi], i-.08, i+.08, color="#4a4a4a", lw=1.5)
                else:
                    ax.text(i, point, "*", ha="left", va="bottom", fontsize=18)
            else:
                ax.text(i, .07, "Not\nreached", transform=ax.get_xaxis_transform(),
                        ha="center", va="bottom", fontsize=10.5, color=".4")
        ax.set_title(DISPLAY[case], pad=12)
        ax.set_xticks(range(5), labels, rotation=25, ha="right", fontsize=13)
        ax.set_ylabel(r"$\tau(T_{\mathrm{H}})/\tau(T_{\mathrm{L}})$" if measure=="threshold_time"
                      else r"$I_5(T_{\mathrm{H}})/I_5(T_{\mathrm{L}})$")
        ax.set_ylim(0, upper)
        ax.set_xlim(-.5, 4.5)
        style.box_axis(ax)
    fig.subplots_adjust(left=.075, right=.99, bottom=.29, top=.88, wspace=.32)
    save(fig, name)


def plot_checks(cache, style):
    for metric in ("Rg", "N_inter", "partner"):
        fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.5))
        for ax, (case, (_, low, high)) in zip(axes, CASES.items()):
            for window, ls in (("early", "-"), ("late", "--")):
                for t, role in ((low, "Tl"), (high, "Th")):
                    lag, mean, _ = cache[case, t, window, metric, "original"]
                    ax.plot(lag, mean, ls=ls, color=style.TEMP_COLORS[role], lw=2,
                            label=style.temp_label(role)+(" early" if window=="early" else " late"))
            ax.set(title=DISPLAY[case], xlabel="Lag time (ns)", ylabel=LABELS[metric], xlim=(0, 5),
                   ylim=(0, 1.04) if metric=="partner" else (-.3, 1.04))
            if case == next(iter(CASES)):
                ax.legend(fontsize=12, loc="upper right", ncol=2, handlelength=1.3, columnspacing=.7)
            style.box_axis(ax)
        fig.subplots_adjust(left=.08, right=.99, bottom=.18, top=.88, wspace=.32)
        save(fig, f"Fig_early_late_{metric}")


def plot_paired_means(rows, style):
    metrics = ("density", "Rg", "N_intra", "N_inter", "N_partner", "load_per_partner")
    fig, axes = plt.subplots(3, 2, figsize=(12, 12))
    for ax, metric in zip(axes.flat, metrics):
        for i, (case, (_, low, high)) in enumerate(CASES.items()):
            for snap in range(1, 9):
                values = [next(r["value"] for r in rows if r["case"]==case and r["snapshot"]==snap
                    and r["temperature_K"]==t and r["window"]=="full" and r["metric"]==metric)
                    for t in (low, high)]
                x = np.array([i-.17, i+.17]) + (snap-4.5)*.013
                ax.plot(x, values, color=".75", lw=.8, zorder=1)
                for j, role in enumerate(("Tl", "Th")):
                    ax.plot(x[j], values[j], "o", color=style.TEMP_COLORS[role], ms=4,
                            label=style.temp_label(role) if i==0 and snap==1 else None)
        ax.set_xticks(range(3), [DISPLAY[c] for c in CASES], fontsize=13)
        ax.set_ylabel(MEAN_LABELS[metric])
        ax.set_xlim(-.5, 2.5)
        if metric == metrics[0]:
            ax.legend(ncol=2, loc="best", handlelength=1, columnspacing=.7)
        style.box_axis(ax)
    fig.subplots_adjust(left=.1, right=.985, bottom=.065, top=.97, hspace=.32, wspace=.3)
    save(fig, "Fig_paired_mean_observables")


def check_pilot_curves(all_curves):
    path = PILOT_SOURCE.parent / "acf_snapshot01_20ns/acf_curves.csv"
    with path.open() as handle:
        old = list(csv.DictReader(handle))
    errors = {}
    for temp in (280, 325):
        for metric in ("Rg", "N_intra", "N_inter", "partner"):
            selected = [r for r in old if int(r["temperature_K"])==temp and r["metric"]==metric
                        and r["window"]=="full" and r["mode"]=="original"]
            selected.sort(key=lambda r: float(r["lag_ns"]))
            reference = np.array([float(r["value"]) for r in selected])
            got = all_curves["FUS", temp, "full", metric, "original"][0]
            np.testing.assert_allclose(got, reference, atol=1e-12, rtol=1e-12)
            errors[f"{temp}_{metric}"] = float(np.max(abs(got-reference)))
    write_json(OUT / "pilot_curve_crosscheck.json", dict(max_abs_errors=errors, source_sha256=sha(path)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--from-cache", action="store_true",
                        help="Reuse validated per-snapshot curves and tables; recompute aggregation and figures")
    args = parser.parse_args()
    test_math()
    progress = json.loads((OUT / "extraction_progress.json").read_text())
    assert progress["status"] == "complete" and len(progress["completed"]) == 48
    all_curves, mean_rows, trends, per_snapshot = {}, [], [], []
    if args.from_cache:
        manifest = json.loads((OUT / "summary_cache_manifest.json").read_text())
        for name, expected in manifest.items():
            assert sha(OUT / name) == expected, f"Summary cache changed: {name}"
        with np.load(OUT / "per_snapshot_curves.npz") as cache:
            for key in cache.files:
                case, temp, window, metric, mode = key.split("__")
                all_curves[case, int(temp), window, metric, mode] = list(cache[key])
        def read_rows(name):
            with (OUT / name).open() as handle:
                rows = list(csv.DictReader(handle))
            for row in rows:
                for key in ("snapshot", "temperature_K", "n_chains"):
                    if key in row:
                        row[key] = int(row[key])
                for key in ("value", "time_ns"):
                    if key in row:
                        row[key] = float(row[key])
            return rows
        mean_rows = read_rows("per_snapshot_means.csv")
        trends = read_rows("per_snapshot_time_evolution.csv")
        per_snapshot = read_rows("per_snapshot_dynamics.csv")
        print("Reusing validated per-snapshot summaries", flush=True)
    else:
        for case, (_, low, high) in CASES.items():
            for snap in range(1, 9):
                for temp in (low, high):
                    curves, rows, trend, stats = branch_summary(case, snap, temp)
                    for key, values in curves.items():
                        all_curves.setdefault((case, temp, *key), []).append(values)
                    mean_rows.extend(rows)
                    trends.extend(trend)
                    per_snapshot.extend(stats)
                print(f"Analyzed {case} snapshot {snap}/8", flush=True)
    check_pilot_curves(all_curves)
    np.savez_compressed(OUT / "per_snapshot_curves.npz", **{
        "__".join(map(str, key)): np.stack(values) for key, values in all_curves.items()})
    curve_rows, time_rows, contrasts, cache = summarize_curves(all_curves)
    means = summarize_means(mean_rows)
    for name, rows in (("acf_ensemble.csv", curve_rows), ("timescales.csv", time_rows),
        ("paired_dynamic_contrasts.csv", contrasts), ("per_snapshot_dynamics.csv", per_snapshot),
        ("per_snapshot_means.csv", mean_rows), ("paired_mean_contrasts.csv", means),
        ("per_snapshot_time_evolution.csv", trends)):
        table(name, rows)
    np.savez_compressed(OUT / "bootstrap_pair_weights.npz", weights=BOOT_COUNTS)
    write_json(OUT / "summary_cache_manifest.json", {name: sha(OUT / name) for name in
        ("per_snapshot_curves.npz", "per_snapshot_means.csv", "per_snapshot_time_evolution.csv", "per_snapshot_dynamics.csv")})
    style = style_module()
    plot_curves(cache, style, ("Rg", "partner"), "Fig_primary_Rg_partner")
    plot_curves(cache, style, ("N_intra", "N_inter"), "Fig_contact_ACFs")
    plot_curves(cache, style, ("Rg", "N_inter", "N_partner"), "Fig_detrended_ACFs", mode="detrended")
    plot_trends(trends, style)
    plot_ratios(contrasts, "Fig_threshold_time_ratios", "threshold_time", style)
    plot_ratios(contrasts, "Fig_integral_0_5ns_ratios", "integral_0_5ns", style)
    plot_checks(cache, style)
    plot_paired_means(mean_rows, style)
    write_json(OUT / "analysis_provenance.json", dict(status="complete", cases=CASES,
        window_ns=[.1, 20.], actual_first_frame_ns=.11, dt_ns=.01, frames_per_branch=1990,
        max_lag_ns=9.95, bootstrap_draws=NBOOT, seed=20261009,
        bootstrap_unit="Eight paired source snapshots; chains averaged within each branch first; equal snapshot weights",
        uncertainty="Conditional on these eight snapshots, fixed force field and definitions; source-trajectory snapshot independence not established",
        density="Centered |z|<5 nm; all beads. No density matching or filtering imposed during analysis",
        cohort="Same initial center chains at both temperatures, retained even if they later leave the core",
        rg="Equal bead weights; bonded minimum-image unwrapping; original biased FFT normalization",
        contacts="0.66 nm periodic bead-pair cutoff; intra sequence separation >2; inter retains all 200 partner identities",
        partner="Endpoint overlap |P(t) intersection P(t+lag)|/|P(t)|; returns allowed, not continuous survival",
        tau="First interpolated crossing, 0.2 for ACFs and 0.5 for partner retention; infinity means right-censored at lag horizon",
        integrals="Signed finite-horizon areas under normalized curves; not infinite-time relaxation estimates",
        censored_bootstrap="Censored draws retained; no ratio CI if any paired bootstrap crossing is censored. Never discard censored draws",
        plots="Mean and pointwise paired-snapshot bootstrap 95% intervals, not simultaneous bands; full negative ACF tails shown",
        ratio_plot_asterisk="Point estimate available but ratio CI omitted because at least one bootstrap draw is censored",
        source_sha256={str(p): sha(p) for p in (Path(__file__), RG_SOURCE, INTRA_SOURCE, PARTNER_SOURCE, PILOT_SOURCE, STYLE_SOURCE)}))
    print("All analyses and figures complete", flush=True)


if __name__ == "__main__":
    main()

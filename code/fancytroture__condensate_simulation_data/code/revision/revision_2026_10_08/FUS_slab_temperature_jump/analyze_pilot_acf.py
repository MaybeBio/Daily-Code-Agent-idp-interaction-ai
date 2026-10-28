#!/usr/bin/env python3
"""Original-observable ACFs on the completed transient pair, with drift diagnostics."""
import ast
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree
from MDAnalysis.coordinates.DCD import DCDReader

from protocol import HERE, ROOT, BOX, whole, validate, write_json

OUT = HERE / "acf_snapshot01_20ns"
CONTACT_SOURCE = ROOT / "revision_2026_09_30/interchain_recompute/compute.py"
INTRA_SOURCE = ROOT / "ppt_5_11_figures/compute_intrachain_contact_relaxation.py"
RG_SOURCE = ROOT / "fus/fus_protein_single_simulation/thermo_figure1_rg_ree_relaxation.py"
spec = importlib.util.spec_from_file_location("contact_operator", CONTACT_SOURCE)
contact_operator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(contact_operator)
DT = .01


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_function(path, name):
    # Execute only the established pure function, never an original production script.
    node = next(n for n in ast.parse(path.read_text()).body if isinstance(n, ast.FunctionDef) and n.name == name)
    namespace = {"np": np}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), namespace)
    return namespace[name]


acf_unbiased = source_function(INTRA_SOURCE, "autocorr_fft")
acf_rg_original = source_function(RG_SOURCE, "autocorr_fft")


def table(name, rows):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def intra_count(x):
    pairs = cKDTree(np.asarray(x, float) % BOX, boxsize=BOX).query_pairs(.66, output_type="ndarray")
    if len(pairs) == 0:
        return 0
    pairs = pairs[(pairs[:, 1] - pairs[:, 0]) > 2]
    displacement = x[pairs[:, 1]] - x[pairs[:, 0]]
    displacement -= np.rint(displacement / BOX) * BOX
    return int(np.count_nonzero(np.sum(displacement**2, axis=1) < .66**2))


def intra_direct(x):
    i, j = np.triu_indices(163, k=3)
    delta = x[j] - x[i]
    delta -= np.rint(delta / BOX) * BOX
    return int(np.count_nonzero(np.sum(delta**2, axis=1) < .66**2))


def partner_curve(contacts, maxlag):
    present = np.asarray(contacts) > 0
    ntime, nchains, _ = present.shape
    curves = []
    for chain in range(nchains):
        values = present[:, chain].astype(float)
        size = values.sum(1)
        valid = size > 0
        weighted = np.divide(values, size[:, None], out=np.zeros_like(values), where=valid[:, None])
        fa = np.fft.rfft(weighted, n=2*ntime, axis=0)
        fb = np.fft.rfft(values, n=2*ntime, axis=0)
        numerator = np.fft.irfft(np.conjugate(fa)*fb, n=2*ntime, axis=0)[:maxlag+1].sum(1)
        denominator = np.array([valid[:ntime-k].sum() for k in range(maxlag+1)])
        curves.append(np.divide(numerator, denominator, out=np.full(maxlag+1, np.nan), where=denominator>0))
    curves = np.array(curves)
    assert np.isfinite(curves).all()
    assert np.all(curves >= -1e-10) and np.all(curves <= 1+1e-10)
    np.testing.assert_allclose(curves[:, 0], 1.)
    return curves.mean(0)


def centered_acf(values, metric, maxlag, detrend=False, force_unbiased=False):
    x = np.asarray(values, float).copy()
    if detrend:
        time_centered = np.arange(len(x)) - (len(x)-1)/2
        slopes = time_centered @ x / (time_centered @ time_centered)
        x -= time_centered[:, None] * slopes[None]
    function = acf_rg_original if metric == "Rg" and not force_unbiased else acf_unbiased
    curves = [function(x[:, chain]) for chain in range(x.shape[1])]
    assert all(c is not None and np.isfinite(c[:maxlag+1]).all() for c in curves)
    curve = np.mean([c[:maxlag+1] for c in curves], axis=0)
    np.testing.assert_allclose(curve[0], 1.)
    return curve


def first_crossing(lag, curve, threshold):
    hit = np.flatnonzero(curve <= threshold)
    if len(hit) == 0:
        return None
    index = int(hit[0])
    if index == 0:
        return float(lag[0])
    return float(lag[index-1] + (threshold-curve[index-1])*(lag[index]-lag[index-1])/(curve[index]-curve[index-1]))


def validate_math():
    rng = np.random.default_rng(33)
    x = rng.normal(size=(60, 3))
    for metric in ("Rg", "N_intra", "N_inter"):
        got = centered_acf(x, metric, 20)
        values = x-x.mean(0)
        expected = []
        for lag in range(21):
            cov = (values[:len(x)-lag] * values[lag:]).mean(0) / (values**2).mean(0)
            if metric == "Rg":
                cov *= (len(x)-lag)/len(x)
            expected.append(cov.mean())
        np.testing.assert_allclose(got, expected, atol=1e-12)
    contacts = rng.integers(0, 4, (60, 3, 5))
    contacts[:, :, 0] = 1
    got = partner_curve(contacts, 20)
    present = contacts > 0
    expected = [np.mean((present[:60-k]&present[k:]).sum(2)/present[:60-k].sum(2)) for k in range(21)]
    np.testing.assert_allclose(got, expected, atol=1e-12)
    xyz = rng.uniform([0, 0, 0], BOX, size=(163, 3))
    xyz[:5] = [[.1, 1, 1], [19.9, 1, 1], [.3, 1, 1], [.4, 1, 1], [.5, 1, 1]]
    assert intra_count(xyz) == intra_direct(xyz)
    shifted = xyz.copy(); shifted[::2] += BOX
    assert intra_count(shifted) == intra_direct(xyz)
    assert first_crossing(np.array([0., 1., 2.]), np.array([1., .4, .1]), .2) == 1+2/3
    assert first_crossing(np.array([0., 1.]), np.array([1., .4]), .2) is None


def extract(temperature):
    validate()
    initial = dict(np.load(HERE / "inputs/snapshot01.npz"))
    cohort = initial["initial_center_chain_ids"].astype(int)
    root = HERE / "runs" / f"snapshot01_{temperature}K"
    parts = sorted(root.glob("part[0-9][0-9][0-9][0-9]"))
    assert [p.name for p in parts] == [f"part{i:04d}" for i in range(1, 41)]
    identity = dict(source_manifests=[sha(p/"complete.json") for p in parts],
                    source_initial=sha(HERE / "inputs/snapshot01.npz"), contact_source=sha(CONTACT_SOURCE),
                    extraction_version=1, dt_ns=DT, cohort=cohort.tolist(), cutoff_nm=.66,
                    intra_sequence_exclusion=2, rg="Equal-bead-weight, matches mdtraj.compute_rg(masses=None)")
    digest = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    cache = OUT / f"observables_full10ps_{temperature}K.npz"
    if cache.exists():
        data = dict(np.load(cache, allow_pickle=False))
        assert str(data["identity_sha256"]) == digest
        return data
    times, rg, intra, partner = [], [], [], []
    checks = []
    started = time.monotonic()
    for index, part in enumerate(parts):
        record = json.loads((part / "complete.json").read_text())
        assert sha(part / "trajectory.dcd") == record["files"]["trajectory.dcd"]
        with DCDReader(str(part / "trajectory.dcd")) as reader:
            assert len(reader) == 50 and reader.n_atoms == 32600
            np.testing.assert_allclose(reader.dt, 10., rtol=1e-6)
            for local, frame in enumerate(reader):
                tick = index * 50 + local + 1
                np.testing.assert_allclose(frame.time / 1000, tick * DT, atol=1e-5)
                np.testing.assert_allclose(frame.dimensions, [200, 200, 2800, 90, 90, 90], atol=.001)
                x = frame.positions.astype(float) / 10
                chain = whole(x)[cohort]
                rg_value = np.sqrt(np.mean(np.sum((chain-chain.mean(1)[:, None])**2, axis=2), axis=1))
                intra_value = np.array([intra_count(x[c*163:(c+1)*163]) for c in cohort])
                result = contact_operator.partner_contact_counts(x, BOX, 163, cohort, cutoff=.66)
                counts = np.stack([result[int(c)] for c in cohort])
                if tick in (10, 200, 1000, 2000):
                    for ci in (0, len(cohort)-1):
                        assert intra_value[ci] == intra_direct(x[cohort[ci]*163:(cohort[ci]+1)*163])
                    checks.append(dict(time_ns=tick*DT, intra_direct_match=True))
                times.append(tick*DT); rg.append(rg_value); intra.append(intra_value); partner.append(counts)
        if (index+1) % 5 == 0:
            print(f"{temperature} K: full-resolution observables through {(index+1)*.5:g} ns; {time.monotonic()-started:.1f} s", flush=True)
    data = dict(time_ns=np.array(times), Rg=np.array(rg), N_intra=np.array(intra),
                partner_contacts=np.array(partner, dtype=np.int32), chain_ids=cohort,
                identity_sha256=np.asarray(digest))
    data["N_inter"] = data["partner_contacts"].sum(2)
    old = dict(np.load(HERE / "analysis_snapshot01_20ns_detailed" / f"observables_{temperature}K.npz"))
    selected = old["time_ns"] > 0
    index = np.rint(old["time_ns"][selected] / DT).astype(int)-1
    np.testing.assert_array_equal(data["partner_contacts"][index], old["partner_contacts"][selected])
    np.savez_compressed(cache, **data)
    write_json(OUT / f"extraction_{temperature}K.json", dict(identity=identity,
        independent_intra_checks=checks, old_inter_observables_exact_match=True,
        rg_note="Original Rg is equal-weight; previous descriptive pilot Rg was mass-weighted. This output follows the original.",
        original_dt_ns=.01, n_frames=len(times)))
    return data


def main():
    OUT.mkdir(exist_ok=True)
    validate_math()
    branches = {temperature: extract(temperature) for temperature in (325, 280)}
    curves, statistics, drifts = [], [], []
    cache = {}
    windows = [("full", .1, 20.), ("early", .1, 10.), ("late", 10., 20.)]
    for temperature, data in branches.items():
        for label, begin, end in windows:
            keep = (data["time_ns"] > begin+1e-8) & (data["time_ns"] <= end+1e-8)
            n = int(keep.sum())
            maxlag = min(1000, n//2)
            lag = np.arange(maxlag+1)*DT
            for metric in ("Rg", "N_intra", "N_inter", "partner"):
                modes = ["original"] if metric == "partner" else ["original", "linear_detrended"]
                if metric == "Rg":
                    modes.append("unbiased_normalization")
                for mode in modes:
                    if metric == "partner":
                        curve = partner_curve(data["partner_contacts"][keep], maxlag)
                    else:
                        curve = centered_acf(data[metric][keep], metric, maxlag,
                                             detrend=mode == "linear_detrended", force_unbiased=mode == "unbiased_normalization")
                    cache[temperature, label, metric, mode] = (lag, curve)
                    threshold = .5 if metric == "partner" else .2
                    crossing = first_crossing(lag, curve, threshold)
                    stats = dict(temperature_K=temperature, window=label, start_ns=begin, end_ns=end,
                        metric=metric, mode=mode, n_frames=n, n_chains=data["Rg"].shape[1], max_lag_ns=float(lag[-1]),
                        threshold=threshold, crossing_ns=crossing, crossing_censored=crossing is None,
                        integral_0_2ns=float(np.trapz(curve[lag<=2.0000001], lag[lag<=2.0000001])),
                        integral_0_5ns=float(np.trapz(curve[lag<=5.0000001], lag[lag<=5.0000001])) if lag[-1]>=5 else None)
                    statistics.append(stats)
                    curves.extend(dict(temperature_K=temperature, window=label, metric=metric, mode=mode,
                                       lag_ns=float(t), value=float(v)) for t, v in zip(lag, curve))
                if metric != "partner":
                    x = data[metric][keep].astype(float)
                    axis = np.arange(n)-(n-1)/2
                    slope = axis@x/(axis@axis)
                    fluct = x-x.mean(0)
                    residual = fluct-axis[:, None]*slope
                    drifts.append(dict(temperature_K=temperature, window=label, metric=metric,
                        mean=float(x.mean()), first_half_mean=float(x[:n//2].mean()), second_half_mean=float(x[n//2:].mean()),
                        fraction_within_chain_variance_removed_by_linear_detrend=float(1-(residual**2).sum()/(fluct**2).sum())))
    ratios = []
    for low in [r for r in statistics if r["temperature_K"] == 280]:
        high = next(r for r in statistics if r["temperature_K"] == 325
                    and all(r[k] == low[k] for k in ("window", "metric", "mode")))
        valid_ratio = low["crossing_ns"] is not None and high["crossing_ns"] is not None and low["crossing_ns"] > 0
        ratios.append(dict(window=low["window"], metric=low["metric"], mode=low["mode"],
            crossing_325_ns=high["crossing_ns"], crossing_280_ns=low["crossing_ns"],
            crossing_ratio_325_over_280=high["crossing_ns"]/low["crossing_ns"] if valid_ratio else None,
            integral2_ratio_325_over_280=high["integral_0_2ns"]/low["integral_0_2ns"] if low["integral_0_2ns"] != 0 else None))
    table("acf_curves.csv", curves); table("acf_statistics.csv", statistics)
    table("acf_ratios.csv", ratios); table("drift_diagnostics.csv", drifts)
    write_json(OUT / "summary.json", dict(single_starting_snapshot=True, chain_count=52, dt_ns=DT,
        full_window_original=[r for r in ratios if r["window"] == "full" and r["mode"] == "original"],
        full_window_detrended=[r for r in ratios if r["window"] == "full" and r["mode"] == "linear_detrended"],
        window_sensitivity_original=[r for r in ratios if r["window"] != "full" and r["mode"] == "original"],
        drift_diagnostics=drifts, pure_function_and_fft_checks_passed=True,
        source_hashes={str(p): sha(p) for p in (CONTACT_SOURCE, INTRA_SOURCE, RG_SOURCE, Path(__file__))},
        limitations=["Transient finite-window ACFs, not demonstrated equilibrium relaxation times.",
            "Original equal-weight Rg and original Rg biased normalization retained; unbiased normalization also exported.",
            "N_intra and N_inter use origin-count-corrected covariance normalized per chain, then average over chains.",
            "Fixed initial cohort, no future-residence selection; not the original 500ns selection protocol.",
            "Linear detrending is a diagnostic, not a uniquely correct removal of nonequilibrium response.",
            "No independent-snapshot confidence intervals from one paired slab; no chain-bootstrap significance claims.",
            "Density differs modestly and evolves; this is not an exact fixed-density equilibrium control."]))
    plot(cache)
    print(json.dumps([r for r in ratios if r["window"] == "full"], indent=2), flush=True)


def plot(cache):
    plt.rcParams.update({"font.size": 11, "pdf.fonttype": 42})
    colors = {325: "#b44562", 280: "#177b78"}
    titles = {"Rg": r"$R_g$", "N_intra": r"$N_{\mathrm{intra}}$", "N_inter": r"$N_{\mathrm{inter}}$", "partner": "Partner retention"}
    for view in ("full_window", "window_sensitivity"):
        fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
        for ax, metric, letter in zip(axes.flat, titles, "ABCD"):
            for temperature, color in colors.items():
                combinations = [("full", "original", "-", "full")]
                if view == "full_window" and metric != "partner":
                    combinations += [("full", "linear_detrended", "--", "detrended")]
                if view == "window_sensitivity":
                    combinations += [("early", "original", "--", "0.1-10 ns"), ("late", "original", ":", "10-20 ns")]
                for window, mode, style, label in combinations:
                    lag, curve = cache[temperature, window, metric, mode]
                    ax.plot(lag, curve, color=color, ls=style, lw=1.6, label=f"{temperature} K, {label}")
            ax.axhline(.5 if metric == "partner" else .2, color=".5", ls=":", lw=.8)
            ax.axhline(0, color=".8", lw=.6)
            ax.set(title=f"{letter}   {titles[metric]}", xlabel="Lag (ns)",
                ylabel="Partner-set overlap" if metric == "partner" else "Normalized ACF",
                xlim=(0, 10 if view == "full_window" else 5), ylim=(-.3, 1.03))
            ax.legend(frameon=False, fontsize=9)
        fig.suptitle("FUS temperature-jump pilot | Transient ACFs, one starting snapshot", fontsize=14)
        for extension in ("png", "pdf"):
            fig.savefig(OUT / f"acf_{view}.{extension}", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    main()

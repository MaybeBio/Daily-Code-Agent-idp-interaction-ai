#!/usr/bin/env python3

import argparse
import csv
import glob
import os
import re

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-paper-thermo")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from paper_style import SYSTEM_COLORS, TEMP_COLORS, add_panel_label, save_figure, setup_style
from thermo_paper_utils import ROOT, SYSTEM_DIRS, out_dir, selected_temps, temp_display_label


DT_NS = 0.01
THRESHOLD = 0.1


def parse_args():
    parser = argparse.ArgumentParser(
        description="Publication-grade thermodynamic overview: Rg, Ree, and relaxation times."
    )
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--free-chunk", type=int, default=2000)
    parser.add_argument("--cond-chunk", type=int, default=2000)
    parser.add_argument("--cond-rg-ree-stride", type=int, default=10)
    parser.add_argument("--n-bootstrap", type=int, default=400)
    parser.add_argument("--outdir", default=out_dir("figure1_rg_ree_relaxation"))
    return parser.parse_args()


def free_paths(temp):
    base = SYSTEM_DIRS["free"][temp]
    return {
        "top": os.path.join(base, "starts", "seed_000.pdb"),
        "traj_glob": os.path.join(base, "outputs_{}{}_cut_stride1".format(temp.replace("K", ""), "K"), "fus_{}_s*_r*.dcd".format(temp)),
        "source_raw": os.path.join(base, "analysis_results_{}".format(temp), "raw_data.npz"),
    }


def cond_paths(temp):
    base = SYSTEM_DIRS["condensate"][temp]
    return {
        "top": os.path.join(base, "starts", "seed_000.pdb"),
        "traj_glob": os.path.join(base, "extracted_500ns_bulk_chains", "chain*_bulk_500ns.dcd"),
        "source_raw": os.path.join(base, "analysis_results_bulk", "raw_data_rg_ree.npz"),
    }


def system_paths(system_name, temp):
    return free_paths(temp) if system_name == "free" else cond_paths(temp)


def natural_sort_key(path):
    name = os.path.basename(path)
    nums = re.findall(r"\d+", name)
    return tuple(int(x) for x in nums) if nums else (name,)


def autocorr_fft(values):
    arr = np.asarray(values, dtype=np.float64)
    arr = arr - np.mean(arr)
    n = len(arr)
    f = np.fft.fft(arr, n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[:n].real
    return acf / acf[0]


def find_precise_t(time_axis, acf, threshold):
    idx = np.where(acf <= threshold)[0]
    if len(idx) == 0:
        return np.nan
    i = int(idx[0])
    if i == 0:
        return float(time_axis[0])
    x0, x1 = float(time_axis[i - 1]), float(time_axis[i])
    y0, y1 = float(acf[i - 1]), float(acf[i])
    return x0 + (threshold - y0) * (x1 - x0) / (y1 - y0)


def bootstrap_relaxation_ci(acf_by_traj, time_axis, threshold, n_boot=400, seed=123):
    acf_by_traj = np.asarray(acf_by_traj, dtype=np.float64)
    if acf_by_traj.ndim != 2 or acf_by_traj.shape[0] == 0:
        return np.nan, np.nan
    if acf_by_traj.shape[0] == 1:
        t_one = find_precise_t(time_axis, acf_by_traj[0], threshold)
        return float(t_one), float(t_one)

    rng = np.random.default_rng(seed)
    n_traj = int(acf_by_traj.shape[0])
    boot_vals = []
    for _ in range(int(n_boot)):
        idx = rng.integers(0, n_traj, size=n_traj)
        acf_mean = np.mean(acf_by_traj[idx], axis=0)
        t_val = find_precise_t(time_axis, acf_mean, threshold)
        if np.isfinite(t_val):
            boot_vals.append(float(t_val))

    if not boot_vals:
        return np.nan, np.nan
    return float(np.percentile(boot_vals, 2.5)), float(np.percentile(boot_vals, 97.5))


def compute_rg_ree_and_relaxation(system_name, temp, free_chunk, cond_chunk, cond_rg_ree_stride, n_bootstrap):
    import mdtraj as md

    paths = system_paths(system_name, temp)
    dcd_files = sorted(glob.glob(paths["traj_glob"]), key=natural_sort_key)
    if not dcd_files:
        raise FileNotFoundError("No trajectories matched {}".format(paths["traj_glob"]))

    chunk = free_chunk if system_name == "free" else cond_chunk
    rg_ree_stride = 1 if system_name == "free" else cond_rg_ree_stride

    all_rg = []
    all_ree = []
    traj_rg = []
    traj_ree = []

    for dcd in dcd_files:
        rg_series = []
        ree_series = []
        for chunk_traj in md.iterload(dcd, top=paths["top"], chunk=chunk, stride=1):
            rg_full = md.compute_rg(chunk_traj).astype(np.float64)
            ree_full = md.compute_distances(chunk_traj, [[0, chunk_traj.n_atoms - 1]])[:, 0].astype(np.float64)
            rg_series.append(rg_full)
            ree_series.append(ree_full)
            all_rg.append(rg_full[::rg_ree_stride])
            all_ree.append(ree_full[::rg_ree_stride])

        rg_concat = np.concatenate(rg_series, axis=0)
        ree_concat = np.concatenate(ree_series, axis=0)
        traj_rg.append(rg_concat)
        traj_ree.append(ree_concat)

    min_lag = min(len(arr) for arr in traj_rg) // 2
    acf_rg_by_traj = np.stack([autocorr_fft(arr)[:min_lag] for arr in traj_rg], axis=0)
    acf_ree_by_traj = np.stack([autocorr_fft(arr)[:min_lag] for arr in traj_ree], axis=0)
    avg_rg_acf = np.mean(acf_rg_by_traj, axis=0)
    avg_ree_acf = np.mean(acf_ree_by_traj, axis=0)
    time_axis = np.arange(min_lag, dtype=np.float64) * DT_NS
    t_rg = find_precise_t(time_axis, avg_rg_acf, THRESHOLD)
    t_ree = find_precise_t(time_axis, avg_ree_acf, THRESHOLD)
    t_rg_ci_low, t_rg_ci_high = bootstrap_relaxation_ci(
        acf_rg_by_traj,
        time_axis,
        THRESHOLD,
        n_boot=n_bootstrap,
        seed=1234,
    )
    t_ree_ci_low, t_ree_ci_high = bootstrap_relaxation_ci(
        acf_ree_by_traj,
        time_axis,
        THRESHOLD,
        n_boot=n_bootstrap,
        seed=5678,
    )

    return {
        "rg": np.concatenate(all_rg, axis=0).astype(np.float32),
        "ree": np.concatenate(all_ree, axis=0).astype(np.float32),
        "acf_rg": avg_rg_acf.astype(np.float32),
        "acf_ree": avg_ree_acf.astype(np.float32),
        "time_axis": time_axis.astype(np.float32),
        "t_rg_ns": np.float32(t_rg),
        "t_ree_ns": np.float32(t_ree),
        "t_rg_ci_low_ns": np.float32(t_rg_ci_low),
        "t_rg_ci_high_ns": np.float32(t_rg_ci_high),
        "t_ree_ci_low_ns": np.float32(t_ree_ci_low),
        "t_ree_ci_high_ns": np.float32(t_ree_ci_high),
        "n_traj": np.int32(len(dcd_files)),
    }


def load_source_raw(system_name, temp):
    path = system_paths(system_name, temp)["source_raw"]
    if not os.path.exists(path):
        return None
    bundle = np.load(path)
    if system_name == "free":
        return np.asarray(bundle["rg"], dtype=np.float32), np.asarray(bundle["ree"], dtype=np.float32)
    return np.asarray(bundle["rg"], dtype=np.float32), np.asarray(bundle["ree"], dtype=np.float32)


def cache_path(outdir, system_name, temp):
    return os.path.join(outdir, "cache_{}_{}.npz".format(system_name, temp))


def load_or_build_cache(outdir, system_name, temp, force_recompute, free_chunk, cond_chunk, cond_rg_ree_stride, n_bootstrap):
    path = cache_path(outdir, system_name, temp)
    required_keys = {
        "rg",
        "ree",
        "acf_rg",
        "acf_ree",
        "time_axis",
        "t_rg_ns",
        "t_ree_ns",
        "t_rg_ci_low_ns",
        "t_rg_ci_high_ns",
        "t_ree_ci_low_ns",
        "t_ree_ci_high_ns",
        "n_traj",
    }
    if os.path.exists(path) and not force_recompute:
        data = np.load(path)
        if required_keys.issubset(set(data.files)):
            return {key: data[key] for key in data.files}

    payload = compute_rg_ree_and_relaxation(system_name, temp, free_chunk, cond_chunk, cond_rg_ree_stride, n_bootstrap)
    source = load_source_raw(system_name, temp)
    if source is not None:
        payload["rg"] = source[0]
        payload["ree"] = source[1]
    np.savez(path, **payload)
    return payload


def kde_curve(values, n_points=300):
    values = np.asarray(values, dtype=np.float64)
    density = gaussian_kde(values)
    lo, hi = float(np.min(values)), float(np.max(values))
    pad = 0.08 * (hi - lo) if hi > lo else 0.1
    xs = np.linspace(lo - pad, hi + pad, n_points)
    ys = density(xs)
    return xs, ys


def linestyle_for_temp(temp):
    styles = {
        "280K": "-",
        "290K": (0, (5, 1.6)),
        "303K": "--",
        "315K": (0, (3, 1.4, 1.2, 1.4)),
        "325K": ":",
    }
    return styles[temp]


def alpha_for_temp(temp):
    alpha_map = {
        "280K": 0.98,
        "290K": 0.90,
        "303K": 0.84,
        "315K": 0.78,
        "325K": 0.72,
    }
    return alpha_map[temp]


def write_summary_csv(path, cache, temps):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "system",
                "temperature",
                "n_samples_rg_ree",
                "rg_mean_nm",
                "rg_std_nm",
                "ree_mean_nm",
                "ree_std_nm",
                "rg_relax_ns",
                "rg_relax_ci_low_ns",
                "rg_relax_ci_high_ns",
                "ree_relax_ns",
                "ree_relax_ci_low_ns",
                "ree_relax_ci_high_ns",
                "n_traj",
            ]
        )
        for system_name in ["free", "condensate"]:
            for temp in temps:
                item = cache[(system_name, temp)]
                writer.writerow(
                    [
                        system_name,
                        temp,
                        int(len(item["rg"])),
                        "{:.8f}".format(float(np.mean(item["rg"]))),
                        "{:.8f}".format(float(np.std(item["rg"]))),
                        "{:.8f}".format(float(np.mean(item["ree"]))),
                        "{:.8f}".format(float(np.std(item["ree"]))),
                        "{:.8f}".format(float(item["t_rg_ns"])),
                        "{:.8f}".format(float(cache_value(item, "t_rg_ci_low_ns", "t_rg_ns"))),
                        "{:.8f}".format(float(cache_value(item, "t_rg_ci_high_ns", "t_rg_ns"))),
                        "{:.8f}".format(float(item["t_ree_ns"])),
                        "{:.8f}".format(float(cache_value(item, "t_ree_ci_low_ns", "t_ree_ns"))),
                        "{:.8f}".format(float(cache_value(item, "t_ree_ci_high_ns", "t_ree_ns"))),
                        int(item["n_traj"]),
                    ]
                )


def cache_value(item, key, fallback_key=None):
    if key in item:
        return item[key]
    if fallback_key is not None and fallback_key in item:
        return item[fallback_key]
    raise KeyError(key)


def asymmetric_yerr(values, lows, highs):
    vals = np.asarray(values, dtype=np.float64)
    lo = np.asarray(lows, dtype=np.float64)
    hi = np.asarray(highs, dtype=np.float64)
    err_low = np.clip(vals - lo, 0.0, None)
    err_high = np.clip(hi - vals, 0.0, None)
    return np.vstack([err_low, err_high])


def main():
    args = parse_args()
    setup_style()
    os.makedirs(args.outdir, exist_ok=True)
    temps = selected_temps()

    cache = {}
    for system_name in ["free", "condensate"]:
        for temp in temps:
            cache[(system_name, temp)] = load_or_build_cache(
                args.outdir,
                system_name,
                temp,
                args.force_recompute,
                args.free_chunk,
                args.cond_chunk,
                args.cond_rg_ree_stride,
                args.n_bootstrap,
            )

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0), dpi=300)
    panel_labels = ["A", "B", "C", "D"]
    for ax, label in zip(axes.ravel(), panel_labels):
        add_panel_label(ax, label)

    ax_rg, ax_ree, ax_trg, ax_tree = axes.ravel()

    for system_name in ["free", "condensate"]:
        for temp in temps:
            item = cache[(system_name, temp)]
            color = SYSTEM_COLORS[system_name]
            ls = linestyle_for_temp(temp)
            alpha = alpha_for_temp(temp)
            label = "{} {}".format("Free" if system_name == "free" else "Condensate", temp_display_label(temp))

            xs, ys = kde_curve(item["rg"])
            ax_rg.plot(xs, ys, color=color, linestyle=ls, lw=2.0, alpha=alpha, label=label)

            xs, ys = kde_curve(item["ree"])
            ax_ree.plot(xs, ys, color=color, linestyle=ls, lw=2.0, alpha=alpha, label=label)

    ax_rg.set_title(r"$R_g$ distributions")
    ax_rg.set_xlabel(r"$R_g$ (nm)")
    ax_rg.set_ylabel("Density")
    ax_rg.grid(alpha=0.18, linestyle="--", linewidth=0.6)
    ax_rg.legend(ncol=2, loc="upper right")

    ax_ree.set_title(r"$R_{ee}$ distributions")
    ax_ree.set_xlabel(r"$R_{ee}$ (nm)")
    ax_ree.set_ylabel("Density")
    ax_ree.grid(alpha=0.18, linestyle="--", linewidth=0.6)
    ax_ree.legend(ncol=2, loc="upper right")

    x_pos = np.arange(len(temps), dtype=np.float64)
    x_labels = [temp_display_label(t) for t in temps]
    free_rg_relax = [float(cache[("free", t)]["t_rg_ns"]) for t in temps]
    cond_rg_relax = [float(cache[("condensate", t)]["t_rg_ns"]) for t in temps]
    free_ree_relax = [float(cache[("free", t)]["t_ree_ns"]) for t in temps]
    cond_ree_relax = [float(cache[("condensate", t)]["t_ree_ns"]) for t in temps]
    free_rg_low = [float(cache_value(cache[("free", t)], "t_rg_ci_low_ns", "t_rg_ns")) for t in temps]
    free_rg_high = [float(cache_value(cache[("free", t)], "t_rg_ci_high_ns", "t_rg_ns")) for t in temps]
    cond_rg_low = [float(cache_value(cache[("condensate", t)], "t_rg_ci_low_ns", "t_rg_ns")) for t in temps]
    cond_rg_high = [float(cache_value(cache[("condensate", t)], "t_rg_ci_high_ns", "t_rg_ns")) for t in temps]
    free_ree_low = [float(cache_value(cache[("free", t)], "t_ree_ci_low_ns", "t_ree_ns")) for t in temps]
    free_ree_high = [float(cache_value(cache[("free", t)], "t_ree_ci_high_ns", "t_ree_ns")) for t in temps]
    cond_ree_low = [float(cache_value(cache[("condensate", t)], "t_ree_ci_low_ns", "t_ree_ns")) for t in temps]
    cond_ree_high = [float(cache_value(cache[("condensate", t)], "t_ree_ci_high_ns", "t_ree_ns")) for t in temps]

    offset = 0.06 if len(temps) > 1 else 0.0
    free_rg_x = x_pos - offset
    cond_rg_x = x_pos + offset
    ax_trg.errorbar(
        free_rg_x,
        free_rg_relax,
        yerr=asymmetric_yerr(free_rg_relax, free_rg_low, free_rg_high),
        marker="o",
        lw=2.1,
        capsize=3.5,
        color=SYSTEM_COLORS["free"],
        label="Free",
    )
    ax_trg.errorbar(
        cond_rg_x,
        cond_rg_relax,
        yerr=asymmetric_yerr(cond_rg_relax, cond_rg_low, cond_rg_high),
        marker="o",
        lw=2.1,
        capsize=3.5,
        color=SYSTEM_COLORS["condensate"],
        label="Condensate",
    )
    for x, y, hi in zip(free_rg_x, free_rg_relax, free_rg_high):
        ax_trg.text(x, max(y, hi) + 0.35, "{:.2f}".format(y), fontsize=8, ha="center", va="bottom")
    for x, y, hi in zip(cond_rg_x, cond_rg_relax, cond_rg_high):
        ax_trg.text(x, max(y, hi) + 0.35, "{:.2f}".format(y), fontsize=8, ha="center", va="bottom")
    ax_trg.set_title(r"$R_g$ relaxation time from autocorrelation")
    ax_trg.set_xlabel("Temperature")
    ax_trg.set_ylabel("Relaxation time (ns)")
    ax_trg.set_xticks(x_pos)
    ax_trg.set_xticklabels(x_labels)
    ax_trg.grid(alpha=0.18, linestyle="--", linewidth=0.6)
    ax_trg.legend(loc="upper right")

    ax_tree.errorbar(
        free_rg_x,
        free_ree_relax,
        yerr=asymmetric_yerr(free_ree_relax, free_ree_low, free_ree_high),
        marker="o",
        lw=2.1,
        capsize=3.5,
        color=SYSTEM_COLORS["free"],
        label="Free",
    )
    ax_tree.errorbar(
        cond_rg_x,
        cond_ree_relax,
        yerr=asymmetric_yerr(cond_ree_relax, cond_ree_low, cond_ree_high),
        marker="o",
        lw=2.1,
        capsize=3.5,
        color=SYSTEM_COLORS["condensate"],
        label="Condensate",
    )
    for x, y, hi in zip(free_rg_x, free_ree_relax, free_ree_high):
        ax_tree.text(x, max(y, hi) + 0.35, "{:.2f}".format(y), fontsize=8, ha="center", va="bottom")
    for x, y, hi in zip(cond_rg_x, cond_ree_relax, cond_ree_high):
        ax_tree.text(x, max(y, hi) + 0.35, "{:.2f}".format(y), fontsize=8, ha="center", va="bottom")
    ax_tree.set_title(r"$R_{ee}$ relaxation time from autocorrelation")
    ax_tree.set_xlabel("Temperature")
    ax_tree.set_ylabel("Relaxation time (ns)")
    ax_tree.set_xticks(x_pos)
    ax_tree.set_xticklabels(x_labels)
    ax_tree.grid(alpha=0.18, linestyle="--", linewidth=0.6)
    ax_tree.legend(loc="upper right")

    fig.suptitle(
        "Thermodynamic overview from chain size statistics and autocorrelation relaxation",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )

    outstem = os.path.join(args.outdir, "figure1_rg_ree_relaxation")
    save_figure(fig, outstem)
    write_summary_csv(os.path.join(args.outdir, "figure1_rg_ree_relaxation_summary.csv"), cache, temps)
    print("Saved:", outstem + ".png")
    print("Saved:", outstem + ".pdf")


if __name__ == "__main__":
    main()

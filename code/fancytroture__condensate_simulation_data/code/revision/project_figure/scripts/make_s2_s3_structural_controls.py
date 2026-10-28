#!/usr/bin/env python3
"""Build SI controls for intrachain-contact distributions and Ree ACFs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MultipleLocator


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "project_figure/supplementary"
SOURCE = OUT / "source_data"
INTRA_DIR = ROOT / "project_figure_set_Tl_Th/source_data/intrachain_relaxation"
SYSTEMS = ("FUS", "A1", "TDP43")
FAMILIES = ("free", "condensate")
ROLES = ("Tl", "Th")
DISPLAY = {"FUS": "FUS-LCD", "A1": "A1-LCD", "TDP43": "TDP-43-LCD"}
TEMPS = {
    "FUS": {"Tl": "280K", "Th": "325K"},
    "A1": {"Tl": "230K", "Th": "260K"},
    "TDP43": {"Tl": "260K", "Th": "300K"},
}
TEMP_COLORS = {"Tl": "#2f7f73", "Th": "#b54b5f"}
CONTACT_CUTOFF = 0.66
ACF_THRESHOLD = 0.2


def phase_legend_handles() -> list[Line2D]:
    return [
        Line2D([0], [0], marker="o", linestyle="None", markersize=9,
               markerfacecolor="white", markeredgecolor="#4a4a4a", markeredgewidth=1.4,
               color="#4a4a4a", label="Free"),
        Line2D([0], [0], marker="o", linestyle="None", markersize=9,
               markerfacecolor="#4a4a4a", markeredgecolor="#4a4a4a", markeredgewidth=1.4,
               color="#4a4a4a", label="Condensate"),
    ]

RG_DIR = {
    "FUS": ROOT / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet/figure1_rg_ree_relaxation",
    "A1": ROOT / "A1/A1_protein_single_simulation/paper_thermo_figures_a1_4p17_relaxed_mu10/figure1_rg_ree_relaxation",
    "TDP43": ROOT / "TDP43/tdp43_protein_single_simulation/paper_thermo_figures_tdp43_4p17_relaxed_mu10/figure1_rg_ree_relaxation",
}


def model_dir(system: str, family: str, role: str) -> Path:
    temp = TEMPS[system][role]
    if system == "FUS":
        base = ROOT / "fus/fus_protein_single_simulation"
        if family == "free":
            return base / f"tica_100dimension_10frames_10_stride_20_trajs_{temp}"
        return base / "condensate" / f"tica_100dimension_10frames_{temp}"
    if system == "A1":
        base = ROOT / "A1/A1_protein_single_simulation"
    else:
        base = ROOT / "TDP43/tdp43_protein_single_simulation"
    target = base / f"tica_100dimension_10frames_{temp}"
    return target if family == "free" else base / "condensate" / target.name


def intermediate(system: str) -> str:
    return {"FUS": "intermediate_fus", "A1": "intermediate_a1", "TDP43": "intermediate_tdp43"}[system]


def setup_style() -> None:
    matplotlib.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 14,
        "axes.labelsize": 16, "xtick.labelsize": 13.5,
        "ytick.labelsize": 13.5, "legend.fontsize": 13,
        "axes.linewidth": 0.9, "legend.frameon": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def box(ax) -> None:
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
    ax.grid(axis="y", color="#bdbdbd", alpha=0.25, linewidth=0.6)


def temp_label(role: str) -> str:
    return r"$T_{\mathrm{L}}$" if role == "Tl" else r"$T_{\mathrm{H}}$"


def line_label(family: str, role: str) -> str:
    return f"{'Free' if family == 'free' else 'Condensate'} {temp_label(role)}"


def stable_seed(*parts: str) -> int:
    return sum((i + 1) * ord(c) for i, c in enumerate("|".join(parts))) % (2**32)


def contact_cache_path(system: str, family: str, role: str) -> Path:
    return SOURCE / f"FigureS2_intrachain_contact_{system}_{family}_{TEMPS[system][role]}.npz"


def compute_contact_sample(system: str, family: str, role: str, samples_per_trajectory: int) -> dict[str, np.ndarray]:
    cache = contact_cache_path(system, family, role)
    if cache.exists():
        data = np.load(cache, allow_pickle=False)
        return {key: data[key] for key in data.files}

    folder = model_dir(system, family, role) / intermediate(system)
    feature_path = folder / "combined_features.npy"
    lengths_path = folder / "combined_tica_lengths.npy"
    features = np.load(feature_path, mmap_mode="r")
    lengths = np.asarray(np.load(lengths_path), dtype=np.int64)
    if int(lengths.sum()) != int(features.shape[0]):
        raise RuntimeError(f"Length mismatch for {system} {family} {role}")

    pooled = []
    trajectory_means = []
    cursor = 0
    for length in lengths:
        count = min(int(samples_per_trajectory), int(length))
        local = np.linspace(0, int(length) - 1, count, dtype=np.int64)
        rows = np.asarray(features[cursor + local], dtype=np.float32)
        contacts = np.count_nonzero(rows < CONTACT_CUTOFF, axis=1).astype(np.float32)
        pooled.append(contacts)
        trajectory_means.append(float(np.mean(contacts)))
        cursor += int(length)

    payload = {
        "sampled_contact_number": np.concatenate(pooled),
        "trajectory_mean_contact_number": np.asarray(trajectory_means, dtype=np.float32),
        "n_trajectories": np.asarray(len(lengths), dtype=np.int32),
        "samples_per_trajectory": np.asarray(samples_per_trajectory, dtype=np.int32),
        "contact_cutoff_nm": np.asarray(CONTACT_CUTOFF, dtype=np.float32),
        "sequence_gap_excluded": np.asarray(2, dtype=np.int32),
    }
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, **payload)
    print(f"[cache] {cache}", flush=True)
    return payload


def heating_response(system: str, family: str, samples_per_trajectory: int, n_bootstrap: int) -> tuple[float, float, float]:
    low = compute_contact_sample(system, family, "Tl", samples_per_trajectory)["trajectory_mean_contact_number"]
    high = compute_contact_sample(system, family, "Th", samples_per_trajectory)["trajectory_mean_contact_number"]
    low_full = float(np.load(INTRA_DIR / f"intrachain_{system}_{family}_{TEMPS[system]['Tl']}.npz")["mean_contact_number"])
    high_full = float(np.load(INTRA_DIR / f"intrachain_{system}_{family}_{TEMPS[system]['Th']}.npz")["mean_contact_number"])
    point = 100.0 * (high_full - low_full) / low_full
    rng = np.random.default_rng(stable_seed(system, family, "contact-response"))
    values = np.empty(n_bootstrap, dtype=float)
    for index in range(n_bootstrap):
        low_mean = float(np.mean(rng.choice(low, size=len(low), replace=True)))
        high_mean = float(np.mean(rng.choice(high, size=len(high), replace=True)))
        values[index] = 100.0 * (high_mean - low_mean) / low_mean
    # Center the trajectory-bootstrap spread on the full-data point estimate.
    values += point - float(np.median(values))
    return point, float(np.percentile(values, 2.5)), float(np.percentile(values, 97.5))


def write_response_csv(rows: list[dict[str, object]]) -> None:
    path = SOURCE / "FigureS2_intrachain_contact_temperature_response.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def smooth_density(values: np.ndarray, sigma_bins: float = 1.35) -> np.ndarray:
    radius = max(1, int(round(4.0 * sigma_bins)))
    x = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (x / sigma_bins) ** 2)
    kernel /= np.sum(kernel)
    return np.convolve(values, kernel, mode="same")


def make_s2(samples_per_trajectory: int, n_bootstrap: int) -> None:
    fig = plt.figure(figsize=(12.4, 10.6), dpi=250)
    grid = GridSpec(3, 2, figure=fig, width_ratios=[2.15, 1.15], hspace=0.34, wspace=0.30,
                    left=0.09, right=0.98, top=0.96, bottom=0.08)
    distribution_axes = []
    response_rows = []

    for row_index, system in enumerate(SYSTEMS):
        ax = fig.add_subplot(grid[row_index, 0])
        distribution_axes.append(ax)
        condition_values = {}
        for family in FAMILIES:
            for role in ROLES:
                condition_values[(family, role)] = np.asarray(
                    compute_contact_sample(system, family, role, samples_per_trajectory)["sampled_contact_number"], float
                )
        lower = min(np.quantile(values, 0.002) for values in condition_values.values())
        upper = max(np.quantile(values, 0.998) for values in condition_values.values())
        bins = np.linspace(lower, upper, 80)
        for family in FAMILIES:
            for role in ROLES:
                values = condition_values[(family, role)]
                density, edges = np.histogram(values, bins=bins, density=True)
                density = smooth_density(density)
                centers = 0.5 * (edges[:-1] + edges[1:])
                ax.plot(centers, density, color=TEMP_COLORS[role],
                        linestyle="--" if family == "free" else "-",
                        linewidth=2.0 if family == "free" else 2.4,
                        label=line_label(family, role))
        ax.set_title(DISPLAY[system], fontsize=16, fontweight="bold", pad=5)
        ax.set_ylabel("Probability density")
        if row_index < 2:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel(r"$N_{\mathrm{intra}}$")
        box(ax)

    distribution_axes[0].legend(loc="upper right", bbox_to_anchor=(0.985, 0.985), ncol=2,
                                fontsize=12.0, handlelength=1.45, columnspacing=0.5,
                                handletextpad=0.35, borderaxespad=0.0)
    distribution_axes[0].text(-0.11, 1.05, "A", transform=distribution_axes[0].transAxes,
                              fontsize=23, fontweight="bold", va="bottom")

    ax_response = fig.add_subplot(grid[:, 1])
    x = np.arange(len(SYSTEMS), dtype=float)
    offsets = {"free": -0.13, "condensate": 0.13}
    for family in FAMILIES:
        stats = [heating_response(system, family, samples_per_trajectory, n_bootstrap) for system in SYSTEMS]
        values = np.asarray([item[0] for item in stats])
        lows = np.asarray([item[1] for item in stats])
        highs = np.asarray([item[2] for item in stats])
        ax_response.errorbar(
            x + offsets[family], values, yerr=np.vstack((values - lows, highs - values)),
            fmt="o", markersize=9, markerfacecolor="white" if family == "free" else "#4a4a4a",
            markeredgecolor="#4a4a4a", markeredgewidth=1.4, color="#555555",
            ecolor="#777777", capsize=3.5, linewidth=1.1,
            label="Free" if family == "free" else "Condensate",
        )
        for system, item in zip(SYSTEMS, stats):
            response_rows.append({"system": DISPLAY[system], "environment": family,
                                  "percent_change": item[0], "ci_low": item[1], "ci_high": item[2],
                                  "contact_cutoff_nm": CONTACT_CUTOFF, "sequence_gap_excluded": 2})
    ax_response.axhline(0.0, color="#777777", linestyle="--", linewidth=1.0)
    ax_response.set_ylim(top=1.25)
    ax_response.set_xticks(x, [DISPLAY[system] for system in SYSTEMS])
    ax_response.set_ylabel("Change upon heating")
    ax_response.yaxis.set_major_locator(MultipleLocator(5))
    ax_response.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}%"))
    ax_response.legend(
        handles=phase_legend_handles(),
        loc="upper center", bbox_to_anchor=(0.5, 0.985), borderaxespad=0.0,
        ncol=2, columnspacing=0.9, handletextpad=0.4,
    )
    ax_response.text(-0.18, 1.02, "B", transform=ax_response.transAxes,
                     fontsize=23, fontweight="bold", va="bottom")
    box(ax_response)
    write_response_csv(response_rows)
    fig.savefig(OUT / "FigureS2_intrachain_contact_number_distributions.pdf", bbox_inches="tight")
    fig.savefig(OUT / "FigureS2_intrachain_contact_number_distributions.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def first_crossing(time: np.ndarray, values: np.ndarray, threshold: float) -> float | None:
    for index in range(1, len(values)):
        if values[index - 1] >= threshold >= values[index]:
            x0, x1 = time[index - 1], time[index]
            y0, y1 = values[index - 1], values[index]
            return float(x1 if np.isclose(y0, y1) else x0 + (threshold - y0) * (x1 - x0) / (y1 - y0))
    return None


def ree_heating_response(system: str, family: str, n_bootstrap: int) -> tuple[float, float, float]:
    """Bootstrap the high/low-temperature ratio of Ree ACF crossing times."""

    def trajectory_acfs(role: str) -> tuple[np.ndarray, np.ndarray]:
        data = np.load(RG_DIR[system] / f"cache_{family}_{TEMPS[system][role]}.npz", allow_pickle=False)
        values = np.asarray(data["ree"], dtype=float)
        n_traj = int(np.asarray(data["n_traj"]).ravel()[0])
        usable = (values.size // n_traj) * n_traj
        trajectories = values[:usable].reshape(n_traj, -1)
        n_lag = trajectories.shape[1] // 2
        acfs = np.empty((n_traj, n_lag), dtype=float)
        for index, trajectory in enumerate(trajectories):
            centered = trajectory - np.mean(trajectory)
            transform = np.fft.rfft(centered, n=2 * centered.size)
            acf = np.fft.irfft(transform * np.conjugate(transform), n=2 * centered.size)[:n_lag]
            acfs[index] = acf / acf[0]
        cached_time = np.asarray(data["time_axis"], dtype=float)
        sampling_factor = (2 * np.asarray(data["acf_ree"]).size) / trajectories.shape[1]
        dt_ns = float(cached_time[1] - cached_time[0]) * sampling_factor
        return acfs, np.arange(n_lag, dtype=float) * dt_ns

    def ensemble_tau(role: str) -> float:
        data = np.load(RG_DIR[system] / f"cache_{family}_{TEMPS[system][role]}.npz", allow_pickle=False)
        tau = first_crossing(np.asarray(data["time_axis"], float), np.asarray(data["acf_ree"], float), ACF_THRESHOLD)
        if tau is None:
            raise RuntimeError(f"Ree ACF does not cross {ACF_THRESHOLD}: {system} {family} {role}")
        return tau

    point = ensemble_tau("Th") / ensemble_tau("Tl")
    low_acfs, low_time = trajectory_acfs("Tl")
    high_acfs, high_time = trajectory_acfs("Th")
    rng = np.random.default_rng(stable_seed("FigureS3", system, family))
    samples = np.empty(n_bootstrap, dtype=float)
    for index in range(n_bootstrap):
        low_mean = np.mean(low_acfs[rng.integers(0, len(low_acfs), len(low_acfs))], axis=0)
        high_mean = np.mean(high_acfs[rng.integers(0, len(high_acfs), len(high_acfs))], axis=0)
        tau_low = first_crossing(low_time, low_mean, ACF_THRESHOLD)
        tau_high = first_crossing(high_time, high_mean, ACF_THRESHOLD)
        samples[index] = np.nan if tau_low is None or tau_high is None else tau_high / tau_low
    finite = samples[np.isfinite(samples)]
    low, high = np.percentile(finite, [2.5, 97.5])
    return point, float(low), float(high)


def make_s3(n_bootstrap: int) -> None:
    fig = plt.figure(figsize=(12.4, 10.6), dpi=250)
    grid = GridSpec(3, 2, figure=fig, width_ratios=[2.15, 1.15], hspace=0.28, wspace=0.30,
                    left=0.09, right=0.98, top=0.96, bottom=0.08)
    axes = [fig.add_subplot(grid[index, 0]) for index in range(3)]
    for row_index, (ax, system) in enumerate(zip(axes, SYSTEMS)):
        for family in FAMILIES:
            for role in ROLES:
                data = np.load(RG_DIR[system] / f"cache_{family}_{TEMPS[system][role]}.npz", allow_pickle=False)
                time = np.asarray(data["time_axis"], float)
                values = np.asarray(data["acf_ree"], float)
                keep = time <= 60.0
                style = "--" if family == "free" else "-"
                ax.plot(time[keep], values[keep], color=TEMP_COLORS[role], linestyle=style,
                        linewidth=2.0 if family == "free" else 2.4, label=line_label(family, role))
                tau = first_crossing(time, values, ACF_THRESHOLD)
                if tau is not None and tau <= 60:
                    ax.vlines(tau, 0, ACF_THRESHOLD, color=TEMP_COLORS[role], linestyle=style,
                              linewidth=1.2, alpha=0.7)
        ax.axhline(ACF_THRESHOLD, color="#555555", linestyle=":", linewidth=1.0)
        ax.set_xlim(0, 60)
        ax.set_ylim(-0.04, 1.03)
        ax.set_ylabel(r"$R_{ee}$ ACF")
        ax.set_title(DISPLAY[system], fontsize=16, fontweight="bold", pad=5)
        box(ax)
    axes[-1].set_xlabel("Lag time (ns)")
    axes[0].legend(loc="upper right", bbox_to_anchor=(0.985, 0.985), ncol=2,
                   fontsize=12.0, handlelength=1.45, columnspacing=0.5,
                   handletextpad=0.35, borderaxespad=0.0)
    axes[0].text(-0.13, 1.04, "A", transform=axes[0].transAxes,
                 fontsize=23, fontweight="bold", va="bottom")

    ax_response = fig.add_subplot(grid[:, 1])
    x = np.arange(len(SYSTEMS), dtype=float)
    offsets = {"free": -0.13, "condensate": 0.13}
    bounds = []
    for family in FAMILIES:
        stats = [ree_heating_response(system, family, n_bootstrap) for system in SYSTEMS]
        values = np.asarray([item[0] for item in stats])
        lows = np.asarray([item[1] for item in stats])
        highs = np.asarray([item[2] for item in stats])
        bounds.extend(lows.tolist() + highs.tolist())
        ax_response.errorbar(
            x + offsets[family], values, yerr=np.vstack((values - lows, highs - values)),
            fmt="o", markersize=9, markerfacecolor="white" if family == "free" else "#4a4a4a",
            markeredgecolor="#4a4a4a", markeredgewidth=1.4, color="#555555",
            ecolor="#777777", capsize=3.5, linewidth=1.1,
            label="Free" if family == "free" else "Condensate",
        )
    ax_response.axhline(1.0, color="#777777", linestyle="--", linewidth=1.0)
    ax_response.set_xticks(x, [DISPLAY[system] for system in SYSTEMS])
    ax_response.set_ylabel(r"$\tau_{R_{ee}}(T_{\mathrm{H}})/\tau_{R_{ee}}(T_{\mathrm{L}})$")
    ax_response.set_ylim(max(0.0, np.floor((min(bounds) - 0.05) * 10) / 10),
                         np.ceil((max(bounds) + 0.05) * 10) / 10)
    ax_response.legend(
        handles=phase_legend_handles(), loc="upper center",
        bbox_to_anchor=(0.5, 0.985), ncol=2, borderaxespad=0.0,
        columnspacing=0.9, handletextpad=0.4,
    )
    ax_response.text(-0.18, 1.02, "B", transform=ax_response.transAxes,
                     fontsize=23, fontweight="bold", va="bottom")
    box(ax_response)
    fig.savefig(OUT / "FigureS3_ree_autocorrelation.pdf", bbox_inches="tight")
    fig.savefig(OUT / "FigureS3_ree_autocorrelation.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-per-trajectory", type=int, default=250)
    parser.add_argument("--n-bootstrap", type=int, default=5000)
    args = parser.parse_args()
    setup_style()
    OUT.mkdir(parents=True, exist_ok=True)
    SOURCE.mkdir(parents=True, exist_ok=True)
    make_s2(args.samples_per_trajectory, args.n_bootstrap)
    make_s3(args.n_bootstrap)
    print("[save] Figures S2 and S3", flush=True)


if __name__ == "__main__":
    main()

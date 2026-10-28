#!/usr/bin/env python3
"""Build Figure 2: a compact dynamics story from old Figures 3, 5, and 6."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(os.environ.get("JPCL_DATA_ROOT", REPOSITORY_ROOT)).resolve()
DEFAULT_OUT = REPOSITORY_ROOT / "outputs/main"

SYSTEM_ORDER = ("FUS", "A1", "TDP43")
DISPLAY_NAMES = {"FUS": "FUS-LCD", "A1": "A1-LCD", "TDP43": "TDP-43-LCD"}
FAMILY_ORDER = ("free", "condensate")
ROLE_ORDER = ("Tl", "Th")
RG_THRESHOLD = 0.2
TEMP_COLORS = {"Tl": "#2f7f73", "Th": "#b54b5f"}
FAMILY_LABELS = {"free": "Free", "condensate": "Condensate"}
OBS_COLORS = {
    "rg": "#4c78a8",
    "intra": "#2f7f73",
    "partner": "#6f78b7",
    "inter": "#c96d28",
}

RG_DIRS = {
    "FUS": PROJECT_ROOT / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet/figure1_rg_ree_relaxation",
    "A1": PROJECT_ROOT / "A1/A1_protein_single_simulation/paper_thermo_figures_a1_4p17_relaxed_mu10/figure1_rg_ree_relaxation",
    "TDP43": PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/paper_thermo_figures_tdp43_4p17_relaxed_mu10/figure1_rg_ree_relaxation",
}
TEMPS = {
    "FUS": {"Tl": "280K", "Th": "325K"},
    "A1": {"Tl": "230K", "Th": "260K"},
    "TDP43": {"Tl": "260K", "Th": "300K"},
}

INTRA_DIR = PROJECT_ROOT / "project_figure_set_Tl_Th/source_data/intrachain_relaxation"
INTER_DIR = PROJECT_ROOT / "storyline_figures_for_ppt_corrected_v3"
RESULTS = PROJECT_ROOT / "ppt_5_11_figures/results"
INTER_RATIO_BOOTSTRAP = (
    PROJECT_ROOT
    / "project_figure/source_data/interchain_ratio_bootstrap/interchain_contact_relaxation_ratio_bootstrap.csv"
)


def setup_style() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 13.5,
            "axes.titlesize": 15.5,
            "axes.labelsize": 17.5,
            "xtick.labelsize": 15.5,
            "ytick.labelsize": 15.5,
            "legend.fontsize": 14.8,
            "legend.frameon": False,
            "axes.linewidth": 0.9,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fval(value: str | float) -> float:
    return float(value)


def stable_seed(*parts: str) -> int:
    text = "|".join(parts)
    return sum((index + 1) * ord(character) for index, character in enumerate(text)) % (2**32)


def display_system(system: str) -> str:
    return DISPLAY_NAMES.get(system, system)


def display_systems() -> list[str]:
    return [display_system(system) for system in SYSTEM_ORDER]


def temp_label(role: str) -> str:
    return r"$T_{\mathrm{L}}$" if role == "Tl" else r"$T_{\mathrm{H}}$"


def clean_axis(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#bdbdbd", alpha=0.25, linewidth=0.6)




def box_axis(ax) -> None:
    """Show a complete rectangular frame while preserving the existing grid."""
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
    ax.grid(axis="y", color="#bdbdbd", alpha=0.25, linewidth=0.6)

def add_panel_labels(fig) -> None:
    labels = [
        ("A", 0.055, 0.958),
        ("B", 0.765, 0.958),
        ("C", 0.055, 0.742),
        ("D", 0.055, 0.528),
        ("E", 0.055, 0.314),
        ("F", 0.055, 0.118),
    ]
    for label, x, y in labels:
        fig.text(x, y, label, fontsize=18, fontweight="bold", ha="left", va="top")


def rg_summary_rows(system: str) -> list[dict[str, str]]:
    return read_csv(RG_DIRS[system] / "figure1_rg_ree_relaxation_summary.csv")


def rg_cache(system: str, family: str, role: str) -> dict[str, np.ndarray]:
    temp = TEMPS[system][role]
    data = np.load(RG_DIRS[system] / f"cache_{family}_{temp}.npz", allow_pickle=False)
    return {key: data[key] for key in data.files}


def load_intra_cache(system: str, family: str, role: str) -> dict[str, np.ndarray]:
    temp = TEMPS[system][role]
    data = np.load(INTRA_DIR / f"intrachain_{system}_{family}_{temp}.npz", allow_pickle=False)
    return {key: data[key] for key in data.files}


def line_label(family: str, role: str) -> str:
    return f"{FAMILY_LABELS[family]} {temp_label(role)}"


def first_crossing_time(time: np.ndarray, series: np.ndarray, threshold: float) -> float | None:
    """Return the first interpolated time where the curve falls below the threshold."""
    time = np.asarray(time, dtype=float)
    series = np.asarray(series, dtype=float)
    for idx in range(1, len(time)):
        y0 = series[idx - 1]
        y1 = series[idx]
        if y0 >= threshold and y1 <= threshold:
            if np.isclose(y0, y1):
                return float(time[idx])
            frac = (threshold - y0) / (y1 - y0)
            return float(time[idx - 1] + frac * (time[idx] - time[idx - 1]))
    return None


def draw_threshold_time(ax, tau: float | None, threshold: float, *, color: str, linestyle: str, lw: float) -> None:
    """Draw a colored vertical line to indicate the threshold-crossing time."""
    if tau is None or not np.isfinite(tau):
        return
    ax.vlines(
        tau,
        0.0,
        threshold,
        colors=color,
        linestyles=linestyle,
        lw=lw,
        alpha=0.65,
        zorder=1,
    )


def rg_relaxation_summary_values(system: str, family: str, role: str) -> tuple[float, float, float]:
    temp = TEMPS[system][role]
    row = next(
        r
        for r in rg_summary_rows(system)
        if r["system"] == family and r["temperature"] == temp
    )
    return (
        fval(row["rg_relax_ns"]),
        fval(row["rg_relax_ci_low_ns"]),
        fval(row["rg_relax_ci_high_ns"]),
    )


def rg_relaxation_time(system: str, family: str, role: str) -> float:
    """Rg ACF crossing time evaluated with the Figure 2 cutoff."""
    data = rg_cache(system, family, role)
    tau = first_crossing_time(
        np.asarray(data["time_axis"], dtype=float),
        np.asarray(data["acf_rg"], dtype=float),
        RG_THRESHOLD,
    )
    if tau is None:
        raise RuntimeError(f"Rg ACF does not cross {RG_THRESHOLD} for {system} {family} {role}")
    return tau


def rg_heating_response_percent(system: str, family: str) -> float:
    tau_tl = rg_relaxation_time(system, family, "Tl")
    tau_th = rg_relaxation_time(system, family, "Th")
    return 100.0 * (tau_th - tau_tl) / tau_tl


def rg_heating_response_interval(system: str, family: str) -> tuple[float, float, float]:
    """Bootstrap the heating response directly from independent trajectories."""

    def trajectory_acfs(role: str) -> tuple[np.ndarray, np.ndarray]:
        data = rg_cache(system, family, role)
        values = np.asarray(data["rg"], dtype=float)
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
        original_frames = 2 * np.asarray(data["acf_rg"]).size
        sampling_factor = original_frames / trajectories.shape[1]
        dt_ns = float(cached_time[1] - cached_time[0]) * sampling_factor
        return acfs, np.arange(n_lag, dtype=float) * dt_ns

    low_acfs, low_time = trajectory_acfs("Tl")
    high_acfs, high_time = trajectory_acfs("Th")
    rng = np.random.default_rng(stable_seed("figure2B", system, family))
    samples = np.empty(5000, dtype=float)
    for sample_index in range(samples.size):
        low_mean = np.mean(
            low_acfs[rng.integers(0, low_acfs.shape[0], low_acfs.shape[0])],
            axis=0,
        )
        high_mean = np.mean(
            high_acfs[rng.integers(0, high_acfs.shape[0], high_acfs.shape[0])],
            axis=0,
        )
        tau_low = first_crossing_time(low_time, low_mean, RG_THRESHOLD)
        tau_high = first_crossing_time(high_time, high_mean, RG_THRESHOLD)
        samples[sample_index] = 100.0 * (tau_high / tau_low - 1.0)

    value = rg_heating_response_percent(system, family)
    low, high = np.percentile(samples[np.isfinite(samples)], [2.5, 97.5])
    return value, float(low), float(high)


def rg_phase_legend_handles() -> list[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markersize=8.5,
            markerfacecolor="white",
            markeredgecolor="#4a4a4a",
            markeredgewidth=1.5,
            color="#4a4a4a",
            label="Free",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="None",
            markersize=8.5,
            markerfacecolor="#4a4a4a",
            markeredgecolor="#4a4a4a",
            markeredgewidth=1.5,
            color="#4a4a4a",
            label="Condensate",
        ),
    ]


def acf_legend_handles() -> list[Line2D]:
    return [
        Line2D([0], [0], color=TEMP_COLORS["Tl"], ls="--", lw=2.0, label=r"Free $T_{\mathrm{L}}$"),
        Line2D([0], [0], color=TEMP_COLORS["Th"], ls="--", lw=2.0, label=r"Free $T_{\mathrm{H}}$"),
        Line2D([0], [0], color=TEMP_COLORS["Tl"], ls="-", lw=2.2, label=r"Condensate $T_{\mathrm{L}}$"),
        Line2D([0], [0], color=TEMP_COLORS["Th"], ls="-", lw=2.2, label=r"Condensate $T_{\mathrm{H}}$"),
    ]


def plot_rg_heating_response(ax, *, show_ylabel: bool = True) -> None:
    """High-to-low temperature ratio of Rg relaxation times."""
    x = np.arange(len(SYSTEM_ORDER), dtype=float)
    offsets = {"free": -0.12, "condensate": 0.12}
    styles = {
        "free": {"mfc": "white", "mec": "#4a4a4a"},
        "condensate": {"mfc": "#4a4a4a", "mec": "#4a4a4a"},
    }

    all_bounds = []
    for family in FAMILY_ORDER:
        stats = [rg_heating_response_interval(system, family) for system in SYSTEM_ORDER]
        values = 1.0 + np.asarray([item[0] for item in stats], dtype=float) / 100.0
        lows = 1.0 + np.asarray([item[1] for item in stats], dtype=float) / 100.0
        highs = 1.0 + np.asarray([item[2] for item in stats], dtype=float) / 100.0
        all_bounds.extend(lows.tolist())
        all_bounds.extend(highs.tolist())
        ax.errorbar(
            x + offsets[family],
            values,
            yerr=np.vstack((values - lows, highs - values)),
            fmt="o",
            markersize=8.5,
            markerfacecolor=styles[family]["mfc"],
            markeredgecolor=styles[family]["mec"],
            markeredgewidth=1.5,
            color="#4a4a4a",
            ecolor="#777777",
            elinewidth=1.0,
            capsize=3.0,
            capthick=1.0,
            zorder=3,
        )

    ax.axhline(1.0, color="#777777", ls="--", lw=0.9, zorder=0)
    ax.set_xticks(x, display_systems())
    ax.set_title("")
    if show_ylabel:
        ax.set_ylabel(r"$\tau_{R_g}(T_{\mathrm{H}})/\tau_{R_g}(T_{\mathrm{L}})$", fontsize=16.5)
    else:
        ax.set_ylabel(None)
        ax.tick_params(labelleft=False)

    y_min = max(0.0, np.floor((min(all_bounds) - 0.05) / 0.1) * 0.1)
    y_max = np.ceil((max(all_bounds) + 0.05) / 0.1) * 0.1
    ax.set_ylim(y_min, y_max)
    ax.tick_params(axis="both", labelsize=15.5)
    clean_axis(ax)


def plot_rg_acf(ax, system: str, show_ylabel: bool = False, show_legend: bool = False) -> None:
    threshold = RG_THRESHOLD
    for family in FAMILY_ORDER:
        for role in ROLE_ORDER:
            cache = rg_cache(system, family, role)
            time = np.asarray(cache["time_axis"], dtype=float)
            mean = np.asarray(cache["acf_rg"], dtype=float)
            keep = time <= 60
            ls = "--" if family == "free" else "-"
            lw = 1.6 if family == "free" else 2.0
            ax.plot(
                time[keep],
                mean[keep],
                color=TEMP_COLORS[role],
                ls=ls,
                lw=lw,
                label=line_label(family, role),
                zorder=3,
            )
            tau = rg_relaxation_time(system, family, role)
            draw_threshold_time(
                ax,
                tau,
                threshold,
                color=TEMP_COLORS[role],
                linestyle=ls,
                lw=1.1 if family == "free" else 1.25,
            )
    ax.axhline(threshold, color="#555555", ls=":", lw=0.9)
    ax.set_title(display_system(system), fontsize=12.5, fontweight="bold", pad=5)
    ax.set_xlabel("Lag time (ns)")
    if show_ylabel:
        ax.set_ylabel(r"$R_g$ ACF")
    ax.set_xlim(0, 60)
    ax.set_ylim(-0.04, 1.03)
    if show_legend:
        ax.legend(ncol=2, loc="upper right", fontsize=7.5, handlelength=1.8, columnspacing=0.7)
    clean_axis(ax)




def rg_cond_free_ratio(system: str, role: str) -> tuple[float, float, float]:
    temp = TEMPS[system][role]
    rows = rg_summary_rows(system)
    free = next(row for row in rows if row["system"] == "free" and row["temperature"] == temp)
    cond = next(row for row in rows if row["system"] == "condensate" and row["temperature"] == temp)
    value = fval(cond["rg_relax_ns"]) / fval(free["rg_relax_ns"])
    low = fval(cond["rg_relax_ci_low_ns"]) / fval(free["rg_relax_ci_high_ns"])
    high = fval(cond["rg_relax_ci_high_ns"]) / fval(free["rg_relax_ci_low_ns"])
    return value, low, high


def plot_phase_slowing(ax) -> None:
    intra_ratio = read_csv(INTRA_DIR / "intrachain_condensate_free_ratio.csv")
    x = np.arange(len(SYSTEM_ORDER), dtype=float)
    specs = [
        ("rg", "Rg relaxation", "o", -0.21),
        ("intra", "Intrachain contact-number relaxation", "s", 0.21),
    ]
    role_offsets = {"Tl": -0.065, "Th": 0.065}

    for obs_key, label, marker, obs_offset in specs:
        for role in ROLE_ORDER:
            values, lows, highs = [], [], []
            for system in SYSTEM_ORDER:
                if obs_key == "rg":
                    value, low, high = rg_cond_free_ratio(system, role)
                else:
                    row = next(r for r in intra_ratio if r["system"] == system and r["role"] == role)
                    value = fval(row["cond_over_free"])
                    low = fval(row["ratio_ci_low"])
                    high = fval(row["ratio_ci_high"])
                values.append(value)
                lows.append(low)
                highs.append(high)
            values = np.asarray(values)
            lows = np.asarray(lows)
            highs = np.asarray(highs)
            ax.errorbar(
                x + obs_offset + role_offsets[role],
                values,
                yerr=np.vstack((values - lows, highs - values)),
                fmt=marker,
                ms=6.2,
                color=TEMP_COLORS[role],
                markerfacecolor=OBS_COLORS[obs_key],
                markeredgecolor="#333333",
                elinewidth=0.9,
                capsize=2.5,
                lw=0,
                label=f"{label} {temp_label(role)}",
            )
    ax.axhline(1, color="#555555", ls="--", lw=1)
    ax.set_xticks(x, display_systems())
    ax.set_ylabel(r"Condensate/free relaxation ratio")
    ax.set_ylim(0, 9.0)
    ax.text(0.98, 0.95, r"$>1$: slower in condensates", transform=ax.transAxes, ha="right", va="top", fontsize=9)
    ax.legend(ncol=2, loc="upper left", fontsize=7.8, handletextpad=0.3, columnspacing=0.8)
    clean_axis(ax)


def plot_temperature_response(ax) -> None:
    intra_rows = read_csv(INTER_DIR / "intra_inter_contact_relaxation/intrachain_contact_relaxation_ratio.csv")
    partner_rows = read_csv(RESULTS / "interchain_partner_survival_ratio.csv")
    inter_rows = read_csv(INTER_RATIO_BOOTSTRAP)
    sources = [
        ("intra", "Intrachain contacts", intra_rows),
        ("inter", "Interchain contact number", inter_rows),
        ("partner", "Partner survival", partner_rows),
    ]
    x = np.arange(len(SYSTEM_ORDER), dtype=float)
    width = 0.23
    for idx, (key, label, rows) in enumerate(sources):
        values = np.asarray([fval(next(row for row in rows if row["system"] == system)["ratio_th_over_tl"]) for system in SYSTEM_ORDER])
        lows = np.asarray([fval(next(row for row in rows if row["system"] == system)["ratio_ci_low"]) for system in SYSTEM_ORDER])
        highs = np.asarray([fval(next(row for row in rows if row["system"] == system)["ratio_ci_high"]) for system in SYSTEM_ORDER])
        ax.bar(
            x + (idx - 1) * width,
            values,
            width,
            color=OBS_COLORS[key],
            edgecolor="#333333",
            linewidth=0.65,
            yerr=np.vstack((values - lows, highs - values)),
            error_kw={"lw": 0.8, "capsize": 2.5, "ecolor": "#444444"},
            label=label,
        )
    ax.axhline(1, color="#555555", ls="--", lw=1)
    ax.set_xticks(x, display_systems())
    ax.set_ylabel(r"$\tau(T_{\mathrm{H}})/\tau(T_{\mathrm{L}})$")
    ax.set_ylim(0, 3.3)
    ax.text(0.98, 0.95, r"$<1$: faster at $T_{\mathrm{H}}$; $>1$: slower at $T_{\mathrm{H}}$", transform=ax.transAxes, ha="right", va="top", fontsize=9)
    ax.legend(loc="upper left", fontsize=8.2)
    clean_axis(ax)


def plot_interchain_acf(ax, system: str, show_ylabel: bool = False, show_legend: bool = False) -> None:
    threshold = 0.2
    curves = read_csv(INTER_DIR / "interchain_contact_relaxation_curves_corrected.csv")
    for role in ROLE_ORDER:
        rows = [row for row in curves if row["system"] == system and row["role"] == role]
        time = np.asarray([fval(row["time_ns"]) for row in rows])
        mean = np.asarray([fval(row["acf_mean"]) for row in rows])
        sem = np.asarray([fval(row["acf_sem"]) for row in rows])
        ax.plot(time, mean, color=TEMP_COLORS[role], lw=1.8, label=temp_label(role), zorder=3)
        ax.fill_between(time, mean - sem, mean + sem, color=TEMP_COLORS[role], alpha=0.14, zorder=2)
        tau = first_crossing_time(time, mean, threshold)
        draw_threshold_time(ax, tau, threshold, color=TEMP_COLORS[role], linestyle="-", lw=1.15)
    ax.axhline(threshold, color="#555555", ls=":", lw=0.9)
    ax.set_title(display_system(system), fontsize=12.5, fontweight="bold", pad=5)
    ax.set_xlabel("Lag time (ns)")
    if show_ylabel:
        ax.set_ylabel("Interchain contact number ACF")
    ax.set_xlim(0, 25)
    ax.set_ylim(0.0, 1.03)
    if show_legend:
        ax.legend(ncol=2, loc="upper right", fontsize=8)
    clean_axis(ax)


def plot_intrachain_contact_acf(ax, system: str, show_ylabel: bool = False, show_legend: bool = False) -> None:
    threshold = 0.2
    for family in FAMILY_ORDER:
        for role in ROLE_ORDER:
            cache = load_intra_cache(system, family, role)
            time = np.asarray(cache["time_ns"], dtype=float)
            mean = np.asarray(cache["acf_mean"], dtype=float)
            sem = np.asarray(cache["acf_sem"], dtype=float)
            keep = time <= 6
            ls = "--" if family == "free" else "-"
            lw = 1.6 if family == "free" else 2.0
            ax.plot(
                time[keep],
                mean[keep],
                color=TEMP_COLORS[role],
                ls=ls,
                lw=lw,
                label=line_label(family, role),
                zorder=3,
            )
            ax.fill_between(
                time[keep],
                mean[keep] - sem[keep],
                mean[keep] + sem[keep],
                color=TEMP_COLORS[role],
                alpha=0.10,
                zorder=2,
            )
            tau = first_crossing_time(time[keep], mean[keep], threshold)
            draw_threshold_time(
                ax,
                tau,
                threshold,
                color=TEMP_COLORS[role],
                linestyle=ls,
                lw=1.1 if family == "free" else 1.25,
            )
    ax.axhline(threshold, color="#555555", ls=":", lw=0.9)
    ax.set_title(display_system(system), fontsize=12.5, fontweight="bold", pad=5)
    ax.set_xlabel("Lag time (ns)")
    if show_ylabel:
        ax.set_ylabel("Intrachain contact number ACF")
    ax.set_xlim(0, 6)
    ax.set_ylim(0.0, 1.03)
    if show_legend:
        ax.legend(ncol=2, loc="upper right", fontsize=7.5, handlelength=1.8, columnspacing=0.7)
    clean_axis(ax)


def plot_partner_survival(ax, system: str, show_ylabel: bool = False, show_legend: bool = False) -> None:
    threshold = 0.5
    curves = read_csv(RESULTS / "interchain_partner_survival_curves.csv")
    for role in ROLE_ORDER:
        rows = [row for row in curves if row["system"] == system and row["role"] == role]
        time = np.asarray([fval(row["time_ns"]) for row in rows])
        mean = np.asarray([fval(row["survival_mean"]) for row in rows])
        sem = np.asarray([fval(row["survival_sem"]) for row in rows])
        ax.plot(time, mean, color=TEMP_COLORS[role], lw=1.8, label=temp_label(role), zorder=3)
        ax.fill_between(time, mean - sem, mean + sem, color=TEMP_COLORS[role], alpha=0.14, zorder=2)
        tau = first_crossing_time(time, mean, threshold)
        draw_threshold_time(ax, tau, threshold, color=TEMP_COLORS[role], linestyle='-', lw=1.15)
    ax.axhline(threshold, color="#555555", ls=":", lw=0.9)
    ax.set_title(display_system(system), fontsize=12.5, fontweight="bold", pad=5)
    ax.set_xlabel("Lag time (ns)")
    if show_ylabel:
        ax.set_ylabel("Partner survival")
    ax.set_xlim(0, 50)
    ax.set_ylim(0.28, 1.03)
    if show_legend:
        ax.legend(ncol=2, loc="upper right", fontsize=8)
    clean_axis(ax)


def plot_multiscale_response_points(ax) -> None:
    intra_rows = read_csv(
        INTER_DIR
        / "intra_inter_contact_relaxation/intrachain_contact_relaxation_ratio.csv"
    )
    partner_rows = read_csv(RESULTS / "interchain_partner_survival_ratio.csv")
    inter_rows = read_csv(INTER_RATIO_BOOTSTRAP)

    # Use a monochrome summary style, but avoid reusing the filled black circle
    # that already denotes "Free" in panel B.
    sources = [
        (r"$N_{\mathrm{intra}}$", intra_rows, "D", -0.20, "#4a4a4a"),
        (r"$N_{\mathrm{inter}}$", inter_rows, "^", 0.00, "white"),
        (r"$S_{\mathrm{partner}}$", partner_rows, "s", 0.20, "white"),
    ]

    x = np.arange(len(SYSTEM_ORDER), dtype=float)

    for label, rows, marker, offset, facecolor in sources:
        selected = [next(row for row in rows if row["system"] == system) for system in SYSTEM_ORDER]
        values = np.asarray(
            [fval(row["ratio_th_over_tl"]) for row in selected],
            dtype=float,
        )
        lows = np.asarray([fval(row["ratio_ci_low"]) for row in selected], dtype=float)
        highs = np.asarray([fval(row["ratio_ci_high"]) for row in selected], dtype=float)

        ax.errorbar(
            x + offset,
            values,
            yerr=np.vstack((values - lows, highs - values)),
            fmt=marker,
            markersize=8.5,
            markerfacecolor=facecolor,
            markeredgecolor="#4a4a4a",
            markeredgewidth=1.25,
            color="#4a4a4a",
            ecolor="#777777",
            elinewidth=1.05,
            capsize=3.0,
            capthick=1.05,
            label=label,
            zorder=3,
        )

    ax.axhline(1.0, color="#777777", ls="--", lw=1.0)
    ax.set_xticks(x, display_systems())
    ax.set_ylabel(r"$\tau(T_{\mathrm{H}})/\tau(T_{\mathrm{L}})$")
    ax.set_ylim(0.0, 4.5)
    ax.set_title("")
    ax.legend(
        handles=[
            Line2D([0], [0], marker="D", linestyle="None", markersize=8.5,
                   markerfacecolor="#4a4a4a", markeredgecolor="#4a4a4a",
                   color="#4a4a4a", label=r"$N_{\mathrm{intra}}$"),
            Line2D([0], [0], marker="^", linestyle="None", markersize=8.5,
                   markerfacecolor="white", markeredgecolor="#4a4a4a",
                   color="#4a4a4a", label=r"$N_{\mathrm{inter}}$"),
            Line2D([0], [0], marker="s", linestyle="None", markersize=8.5,
                   markerfacecolor="white", markeredgecolor="#4a4a4a",
                   color="#4a4a4a", label=r"$S_{\mathrm{partner}}$"),
        ],
        ncol=3,
        loc="upper left",
        bbox_to_anchor=(0.0, 0.965),
        fontsize=15.0,
        columnspacing=1.25,
        handletextpad=0.50,
        borderaxespad=0.0,
    )
    clean_axis(ax)


def build_figure(output_dir: Path) -> None:
    fig = plt.figure(figsize=(16.8, 13.1), dpi=300)
    outer = GridSpec(4, 4, figure=fig, height_ratios=[1.0, 1.0, 1.0, 1.15], hspace=0.34, wspace=0.40, left=0.075, right=0.975, top=0.958, bottom=0.055)

    first_axes = {}
    row_centers = {}

    # Create the 3x4 matrix for A/B/D/E (columns) x systems (rows)
    column_specs = [
        ("A", plot_rg_acf, r"$R_g$ ACF"),
        ("B", plot_intrachain_contact_acf, r"$N_{\mathrm{intra}}$ ACF"),
        ("D", plot_interchain_acf, r"$N_{\mathrm{inter}}$ ACF"),
        ("E", plot_partner_survival, r"$S_{\mathrm{partner}}$"),
    ]

    axes_matrix = {}
    for col_idx, (panel_label, plotter, metric_title) in enumerate(column_specs):
        col_axes = []
        for row_idx, system in enumerate(SYSTEM_ORDER):
            ax = fig.add_subplot(outer[row_idx, col_idx])
            col_axes.append(ax)
            if row_idx == 0:
                first_axes[panel_label] = ax
            plotter(ax, system, show_ylabel=False, show_legend=False)
            box_axis(ax)
            # Place system identities consistently as subplot titles.
            ax.set_title(display_system(system), fontsize=16.0,
                         fontweight="bold", pad=5)
            label_size = 17.5
            ax.set_ylabel(metric_title, fontsize=label_size)
            if row_idx < len(SYSTEM_ORDER) - 1:
                ax.set_xlabel("")
                ax.tick_params(labelbottom=False)
        axes_matrix[panel_label] = col_axes

    # A defines the shared free/condensate and temperature encoding for A/B.
    axes_matrix["A"][0].legend(
        handles=acf_legend_handles(), loc="upper right",
        bbox_to_anchor=(0.985, 0.985), ncol=2, frameon=False,
        fontsize=12.5, handlelength=1.5, columnspacing=0.45,
        handletextpad=0.25, borderaxespad=0.0,
    )

    # Bottom row: C under A; F spanning BDE.
    # C and F now use the exact same vertical extent, so their panel labels,
    # titles, upper edges, and lower edges are aligned. Extra headroom is kept
    # inside B's y-range for the legend between the title and the data points.
    ax_f = fig.add_subplot(outer[3, 1:4])
    first_axes['F'] = ax_f
    plot_multiscale_response_points(ax_f)
    box_axis(ax_f)

    ax_b = fig.add_subplot(outer[3, 0])
    first_axes['C'] = ax_b
    plot_rg_heating_response(ax_b, show_ylabel=True)
    box_axis(ax_b)

    f_pos = ax_f.get_position()
    b_pos = ax_b.get_position()
    ax_b.set_position([b_pos.x0, f_pos.y0, b_pos.width, f_pos.height])

    ax_b.set_title("")
    # Reserve a dedicated strip above the highest confidence interval so the
    # phase legend never obscures an error bar.
    ax_b.set_ylim(ax_b.get_ylim()[0], max(1.65, ax_b.get_ylim()[1]))
    ax_b.legend(
        handles=rg_phase_legend_handles(),
        labels=["Free", "Condensate"],
        loc="upper right",
        bbox_to_anchor=(0.985, 0.965),
        ncol=2,
        frameon=False,
        fontsize=14.5,
        handletextpad=0.40,
        columnspacing=0.75,
        borderaxespad=0.0,
    )

    # Panel letters
    fig.canvas.draw()
    for label in ('A', 'B', 'C', 'D', 'E', 'F'):
        box = first_axes[label].get_position()
        fig.text(max(0.012, box.x0 - 0.03), min(box.y1 + 0.020, 0.992), label, fontsize=24, fontweight='bold', ha='left', va='bottom')

    output_dir.mkdir(parents=True, exist_ok=True)
    name = 'Figure2_dynamics_story_columnwise_layout_v22'
    for ext in ('png', 'pdf'):
        fig.savefig(output_dir / f'{name}.{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'[save] {output_dir / name}.png')
    print(f'[save] {output_dir / name}.pdf')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    setup_style()
    build_figure(args.output_dir)


if __name__ == "__main__":
    main()

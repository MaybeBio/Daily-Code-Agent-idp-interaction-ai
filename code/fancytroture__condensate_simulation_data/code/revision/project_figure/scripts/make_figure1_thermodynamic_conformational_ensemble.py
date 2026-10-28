#!/usr/bin/env python3
"""Build Figure 1 from old Figure 2A-C.

Panels (combined mode):
    A. Rg distributions for free and condensate chains.
    B. Ree distributions for free and condensate chains.
    C. Temperature-induced Rg/Ree response.
    D. Temperature-induced intrachain contact-map difference.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter, FuncFormatter, MultipleLocator


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = Path(__file__).resolve().parents[1]

SYSTEM_ORDER = ("FUS", "A1", "TDP43")
DISPLAY_NAMES = {"FUS": "FUS-LCD", "A1": "A1-LCD", "TDP43": "TDP-43-LCD"}
FAMILY_ORDER = ("free", "condensate")
ROLE_ORDER = ("Tl", "Th")

TEMP_COLORS = {"Tl": "#2f7f73", "Th": "#b54b5f"}
FAMILY_LABELS = {"free": "Free", "condensate": "Condensate"}

# Exclude local sequence-neighbor contacts from panel D before plotting and
# color scaling. A gap of 8 masks the main diagonal and the first eight
# off-diagonals, i.e. residue pairs with |i-j| <= 8.
CONTACT_SEQUENCE_GAP = {
    "FUS": 8,
    "A1": 8,
    "TDP43": 8,
}

SYSTEMS = {
    "FUS": {
        "temps": {"Tl": "280K", "Th": "325K"},
        "rg_dir": PROJECT_ROOT
        / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet/figure1_rg_ree_relaxation",
        "state": {
            "free": PROJECT_ROOT
            / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet/figure3_state_population_contactmaps/figure3_free_state_population_contactmaps_source.csv",
            "condensate": PROJECT_ROOT
            / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet/figure3_state_population_contactmaps/figure3_condensate_state_population_contactmaps_source.csv",
        },
        "contact": {
            ("free", "280K"): PROJECT_ROOT
            / "fus/fus_protein_single_simulation/tica_100dimension_10frames_10_stride_20_trajs_280K",
            ("free", "325K"): PROJECT_ROOT
            / "fus/fus_protein_single_simulation/tica_100dimension_10frames_10_stride_20_trajs_325K",
            ("condensate", "280K"): PROJECT_ROOT
            / "fus/fus_protein_single_simulation/condensate/tica_100dimension_10frames_280K",
            ("condensate", "325K"): PROJECT_ROOT
            / "fus/fus_protein_single_simulation/condensate/tica_100dimension_10frames_325K",
        },
    },
    "A1": {
        "temps": {"Tl": "230K", "Th": "260K"},
        "rg_dir": PROJECT_ROOT
        / "A1/A1_protein_single_simulation/paper_thermo_figures_a1_4p17_relaxed_mu10/figure1_rg_ree_relaxation",
        "state": {
            "free": PROJECT_ROOT
            / "ppt_5_11_figures/results/a1_230_260_state_mapping/figure3_free_state_population_contactmaps_source.csv",
            "condensate": PROJECT_ROOT
            / "ppt_5_11_figures/results/a1_230_260_state_mapping/figure3_condensate_state_population_contactmaps_source.csv",
        },
        "contact": {
            ("free", "230K"): PROJECT_ROOT / "A1/A1_protein_single_simulation/tica_100dimension_10frames_230K",
            ("free", "260K"): PROJECT_ROOT / "A1/A1_protein_single_simulation/tica_100dimension_10frames_260K",
            ("condensate", "230K"): PROJECT_ROOT
            / "A1/A1_protein_single_simulation/condensate/tica_100dimension_10frames_230K",
            ("condensate", "260K"): PROJECT_ROOT
            / "A1/A1_protein_single_simulation/condensate/tica_100dimension_10frames_260K",
        },
    },
    "TDP43": {
        "temps": {"Tl": "260K", "Th": "300K"},
        "rg_dir": PROJECT_ROOT
        / "TDP43/tdp43_protein_single_simulation/paper_thermo_figures_tdp43_4p17_relaxed_mu10/figure1_rg_ree_relaxation",
        "state": {
            "free": PROJECT_ROOT
            / "TDP43/tdp43_protein_single_simulation/paper_thermo_figures_tdp43_4p17_relaxed_mu10/figure3_state_population_contactmaps/figure3_free_state_population_contactmaps_source.csv",
            "condensate": PROJECT_ROOT
            / "TDP43/tdp43_protein_single_simulation/paper_thermo_figures_tdp43_4p17_relaxed_mu10/figure3_state_population_contactmaps/figure3_condensate_state_population_contactmaps_source.csv",
        },
        "contact": {
            ("free", "260K"): PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/tica_100dimension_10frames_260K",
            ("free", "300K"): PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/tica_100dimension_10frames_300K",
            ("condensate", "260K"): PROJECT_ROOT
            / "TDP43/tdp43_protein_single_simulation/condensate/tica_100dimension_10frames_260K",
            ("condensate", "300K"): PROJECT_ROOT
            / "TDP43/tdp43_protein_single_simulation/condensate/tica_100dimension_10frames_300K",
        },
    },
}


def contact_cmap():
    cmap = matplotlib.colormaps["RdBu_r"].copy()
    cmap.set_bad("white")
    return cmap


def nice_symmetric_limit(vmax: float) -> float:
    """Round a positive limit up to a clean scientific-notation bound."""
    vmax = float(abs(vmax))
    if vmax == 0.0:
        return 1.0
    exponent = int(np.floor(np.log10(vmax)))
    scaled = vmax / (10.0 ** exponent)
    if scaled <= 1.0:
        nice = 1.0
    elif scaled <= 2.0:
        nice = 2.0
    elif scaled <= 5.0:
        nice = 5.0
    else:
        nice = 10.0
    return nice * (10.0 ** exponent)


def colorbar_ticks_for_limit(vlim: float) -> np.ndarray:
    return np.linspace(-vlim, vlim, 5)


def scientific_scale_exponent(vlim: float) -> int:
    """Return the base-10 exponent used to scale colorbar tick labels."""
    vlim = float(abs(vlim))
    if vlim == 0.0:
        return 0
    return int(np.floor(np.log10(vlim)))


def apply_scaled_colorbar_labels(colorbar, vlim: float) -> None:
    """Show integer tick labels plus a single shared scientific multiplier."""
    exponent = int(np.floor(np.log10(abs(vlim)))) if vlim != 0 else 0
    scale = 10.0 ** exponent

    raw_ticks = colorbar_ticks_for_limit(vlim)
    scaled_ticks = raw_ticks / scale

    labels = []
    for value in scaled_ticks:
        if np.isclose(value, 0.0):
            labels.append("0")
        elif np.isclose(value, round(value)):
            labels.append(str(int(round(value))))
        else:
            labels.append(f"{value:g}")

    colorbar.set_ticks(raw_ticks)
    colorbar.set_ticklabels(labels)

    # Disable any automatic offset/scientific text so only our multiplier remains.
    colorbar.ax.yaxis.get_offset_text().set_visible(False)

    colorbar.ax.text(
        0.5,
        1.02,
        rf"$\times 10^{{{exponent}}}$",
        transform=colorbar.ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=11,
    )


def setup_style() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 15,
            "axes.titlesize": 17,
            "axes.labelsize": 16,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 13,
            "legend.frameon": False,
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


def row_temperature(row: dict[str, str]) -> str:
    return row.get("temperature", row.get("temp", ""))


def display_system(system: str) -> str:
    return DISPLAY_NAMES.get(system, system)


def temp_label(role: str, include_temp: bool = False, system: str | None = None) -> str:
    base = r"$T_{\mathrm{L}}$" if role == "Tl" else r"$T_{\mathrm{H}}$"
    if include_temp and system is not None:
        return f"{base} ({SYSTEMS[system]['temps'][role]})"
    return base


def panel_label(ax, label: str, x: float = -0.12, y: float = 1.08) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        ha="left",
        va="top",
        clip_on=False,
    )


def clean_axis(ax) -> None:
    # Keep a complete rectangular frame for panels A, B, and C.
    for spine in ax.spines.values():
        spine.set_visible(True)
    ax.grid(axis="y", color="#bdbdbd", alpha=0.25, linewidth=0.6)


def observable_label(observable: str) -> str:
    return r"$R_g$" if observable == "rg" else r"$R_{ee}$"


def smooth_density(density: np.ndarray, sigma_bins: float = 1.15) -> np.ndarray:
    """Apply light Gaussian smoothing to histogram densities."""
    radius = max(1, int(np.ceil(3 * sigma_bins)))
    x = np.arange(-radius, radius + 1, dtype=float)
    kernel = np.exp(-0.5 * (x / sigma_bins) ** 2)
    kernel /= kernel.sum()
    padded = np.pad(np.asarray(density, dtype=float), (radius, radius), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def rg_cache(system: str, family: str, role: str) -> dict[str, np.ndarray]:
    temp = SYSTEMS[system]["temps"][role]
    data = np.load(SYSTEMS[system]["rg_dir"] / f"cache_{family}_{temp}.npz", allow_pickle=False)
    return {key: data[key] for key in data.files}


def rg_summary_rows(system: str) -> list[dict[str, str]]:
    return read_csv(SYSTEMS[system]["rg_dir"] / "figure1_rg_ree_relaxation_summary.csv")


def relative_response_percent(system: str, family: str, observable: str) -> float:
    """Mean-based heating response: 100 * (high - low) / low."""
    rows = rg_summary_rows(system)
    low_temp = SYSTEMS[system]["temps"]["Tl"]
    high_temp = SYSTEMS[system]["temps"]["Th"]
    low_row = next(r for r in rows if r["system"] == family and row_temperature(r) == low_temp)
    high_row = next(r for r in rows if r["system"] == family and row_temperature(r) == high_temp)
    low_mean = fval(low_row[f"{observable}_mean_nm"])
    high_mean = fval(high_row[f"{observable}_mean_nm"])
    return 100.0 * (high_mean - low_mean) / low_mean


def trajectory_means(data: dict[str, np.ndarray], key: str) -> np.ndarray:
    values = np.asarray(data[key], dtype=float)
    n_traj = int(np.asarray(data["n_traj"]).ravel()[0])
    usable = (values.size // n_traj) * n_traj
    if usable == 0:
        raise ValueError(f"No usable trajectory blocks for {key}")
    return values[:usable].reshape(n_traj, -1).mean(axis=1)


def stable_seed(*parts: str) -> int:
    text = "|".join(parts)
    return sum((idx + 1) * ord(char) for idx, char in enumerate(text)) % (2**32)


def response_bootstrap_interval(
    system: str,
    family: str,
    observable: str,
    n_boot: int = 5000,
) -> tuple[float, float, float]:
    point = relative_response_percent(system, family, observable)
    low = trajectory_means(rg_cache(system, family, "Tl"), observable)
    high = trajectory_means(rg_cache(system, family, "Th"), observable)
    rng = np.random.default_rng(stable_seed("figure1C", system, family, observable))
    samples = np.empty(n_boot, dtype=float)
    for sample_idx in range(n_boot):
        low_sample = low[rng.integers(0, len(low), len(low))]
        high_sample = high[rng.integers(0, len(high), len(high))]
        samples[sample_idx] = 100.0 * (np.mean(high_sample) - np.mean(low_sample)) / np.mean(low_sample)
    ci_low, ci_high = np.percentile(samples, [2.5, 97.5])
    return point, float(ci_low), float(ci_high)


def plot_response_summary(ax, observable: str, *, ylim: tuple[float, float] | None = None) -> None:
    xs = np.arange(len(SYSTEM_ORDER), dtype=float)
    offsets = {"free": -0.12, "condensate": 0.12}
    family_style = {
        "free": {"marker": "o", "mfc": "white", "mec": "#4a4a4a", "color": "#4a4a4a"},
        "condensate": {"marker": "o", "mfc": "#4a4a4a", "mec": "#4a4a4a", "color": "#4a4a4a"},
    }

    # Plot free and condensate as paired symbols without connecting lines.
    for family in FAMILY_ORDER:
        stats = [response_bootstrap_interval(system, family, observable) for system in SYSTEM_ORDER]
        vals = np.asarray([item[0] for item in stats], dtype=float)
        ci_low = np.asarray([item[1] for item in stats], dtype=float)
        ci_high = np.asarray([item[2] for item in stats], dtype=float)
        yerr = np.vstack([
            np.maximum(vals - ci_low, 0.0),
            np.maximum(ci_high - vals, 0.0),
        ])
        style = family_style[family]
        ax.errorbar(
            xs + offsets[family],
            vals,
            yerr=yerr,
            fmt=style["marker"],
            markersize=10.5,
            markerfacecolor=style["mfc"],
            markeredgecolor=style["mec"],
            markeredgewidth=1.5,
            color=style["color"],
            ecolor=style["color"],
            elinewidth=1.2,
            capsize=3.0,
            capthick=1.2,
            linestyle="None",
            zorder=3,
        )

    ax.axhline(0.0, color="#9c9c9c", lw=1.0, ls="--", zorder=0)
    ax.set_xticks(xs, [display_system(system) for system in SYSTEM_ORDER])
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}%"))
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.tick_params(axis="both", labelsize=12)
    clean_axis(ax)


def response_legend_handles() -> list[Line2D]:
    return [
        Line2D([0], [0], marker="o", linestyle="None", markersize=10.5,
               markerfacecolor="white", markeredgecolor="#4a4a4a", markeredgewidth=1.5,
               color="#4a4a4a", label="Free"),
        Line2D([0], [0], marker="o", linestyle="None", markersize=10.5,
               markerfacecolor="#4a4a4a", markeredgecolor="#4a4a4a", markeredgewidth=1.5,
               color="#4a4a4a", label="Condensate"),
    ]


def state_rows(system: str, family: str) -> list[dict[str, str]]:
    return read_csv(SYSTEMS[system]["state"][family])


def ensemble_contact_map(system: str, family: str, role: str) -> np.ndarray:
    temp = SYSTEMS[system]["temps"][role]
    rows = [row for row in state_rows(system, family) if row_temperature(row) == temp]
    if not rows:
        raise RuntimeError(f"No state rows for {system} {family} {temp}")

    contacts = np.load(SYSTEMS[system]["contact"][(family, temp)] / "avg_contacts.npy")
    numerator = np.zeros(contacts.shape[1:], dtype=np.float64)
    denominator = 0.0

    for row in rows:
        local_state = int(row["local_state"])
        mu = max(fval(row["mu"]), 0.0)
        if local_state < 0 or local_state >= contacts.shape[0]:
            raise IndexError(
                f"local_state {local_state} out of range for {system} {family} {temp}; "
                f"avg_contacts has {contacts.shape[0]} states"
            )
        numerator += mu * np.asarray(contacts[local_state], dtype=np.float64)
        denominator += mu

    if denominator <= 0.0:
        raise RuntimeError(f"Zero state population weight for {system} {family} {temp}")
    return numerator / denominator


def mask_contact_map(contact: np.ndarray, gap: int) -> np.ndarray:
    """Mask local sequence-neighbor contacts with |i-j| <= gap."""
    result = np.asarray(contact, dtype=np.float64).copy()
    n_residues = result.shape[0]
    if gap < 0:
        raise ValueError(f"gap must be non-negative, got {gap}")
    for delta in range(min(gap, n_residues - 1) + 1):
        indices = np.arange(n_residues - delta)
        result[indices, indices + delta] = np.nan
        result[indices + delta, indices] = np.nan
    return result


def temperature_delta_contact_map(system: str, family: str) -> np.ndarray:
    low = ensemble_contact_map(system, family, "Tl")
    high = ensemble_contact_map(system, family, "Th")
    n_residues = min(low.shape[0], high.shape[0])
    delta = high[:n_residues, :n_residues] - low[:n_residues, :n_residues]
    gap = CONTACT_SEQUENCE_GAP.get(system, 8)
    return mask_contact_map(delta, gap=gap)


def plot_size_distribution(
    ax, system: str, observable: str, *, show_ylabel: bool = True
) -> None:
    all_values = []
    sampled: dict[tuple[str, str], np.ndarray] = {}
    for role in ROLE_ORDER:
        for family in FAMILY_ORDER:
            values = np.asarray(rg_cache(system, family, role)[observable], dtype=float)
            sample = values[:: max(1, values.size // 200000)]
            sampled[(role, family)] = sample
            all_values.append(sample)

    lo = min(float(np.quantile(values, 0.002)) for values in all_values)
    hi = max(float(np.quantile(values, 0.998)) for values in all_values)
    bins = np.linspace(lo, hi, 120)

    for role in ROLE_ORDER:
        for family in FAMILY_ORDER:
            values = sampled[(role, family)]
            density, edges = np.histogram(values, bins=bins, density=True)
            density = smooth_density(density, sigma_bins=1.2 if observable == "rg" else 1.35)
            centers = 0.5 * (edges[:-1] + edges[1:])
            ax.plot(
                centers,
                density,
                color=TEMP_COLORS[role],
                lw=2.4 if family == "condensate" else 2.0,
                ls="-" if family == "condensate" else "--",
                label=f"{FAMILY_LABELS[family]} {temp_label(role)}",
            )
    ax.set_title(display_system(system), fontsize=17, fontweight="bold")
    ax.set_xlabel(f"{observable_label(observable)} (nm)")
    if show_ylabel:
        ax.set_ylabel("Probability density")
    else:
        ax.set_ylabel("")
        ax.tick_params(labelleft=False)
    if observable == "rg":
        ax.set_xlim(1.5, 6.2)
        ax.set_ylim(0.0, 1.08)
        ax.xaxis.set_major_locator(MultipleLocator(1))
    else:
        ax.set_xlim(0.8, 18.0)
        ax.set_ylim(0.0, 0.19)
        ax.xaxis.set_major_locator(MultipleLocator(2))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    clean_axis(ax)


def build_figure(output_dir: Path, observable: str) -> None:
    """Build a two-panel figure for a single size observable plus contact maps."""
    delta_maps = {}
    delta_abs = []
    for system in SYSTEM_ORDER:
        for family in FAMILY_ORDER:
            delta = temperature_delta_contact_map(system, family)
            delta_maps[(system, family)] = delta
            delta_abs.extend(np.abs(delta[np.isfinite(delta)]).ravel())
    vmax = nice_symmetric_limit(float(np.quantile(np.asarray(delta_abs), 0.995)))

    fig = plt.figure(figsize=(16.5, 8.8), dpi=300)
    gs = GridSpec(3, 6, figure=fig, height_ratios=[1.0, 1.25, 1.25], hspace=0.52, wspace=0.42)
    panel_x = 0.082
    fig.text(panel_x, 0.925, "A", fontsize=18, fontweight="bold", ha="left", va="top")
    fig.text(panel_x, 0.585, "B", fontsize=18, fontweight="bold", ha="left", va="top")

    for system_index, system in enumerate(SYSTEM_ORDER):
        ax = fig.add_subplot(gs[0, system_index * 2 : system_index * 2 + 2])
        plot_size_distribution(ax, system, observable)
        if system_index == 0:
            ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.02), ncol=1)

    image = None
    contact_axes = []
    for family_index, family in enumerate(FAMILY_ORDER):
        for system_index, system in enumerate(SYSTEM_ORDER):
            ax = fig.add_subplot(gs[1 + family_index, system_index * 2 : system_index * 2 + 2])
            contact_axes.append(ax)
            image = ax.imshow(
                delta_maps[(system, family)],
                origin="lower",
                cmap=contact_cmap(),
                norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax),
                interpolation="nearest",
            )
            if family_index == 0:
                ax.text(
                    0.5, 1.06, display_system(system), transform=ax.transAxes,
                    ha="center", va="bottom", fontsize=13, fontweight="bold",
                )
            if system_index == 0:
                ax.text(
                    -0.30, 0.5, FAMILY_LABELS[family], transform=ax.transAxes,
                    rotation=90, ha="center", va="center", fontsize=12, fontweight="bold",
                )
            ax.set_xlabel("Residue index")
    colorbar = fig.colorbar(image, ax=contact_axes, fraction=0.018, pad=0.015, ticks=colorbar_ticks_for_limit(vmax))
    colorbar.set_label(r"$\Delta P_{\mathrm{intra}} = P(T_{\mathrm{H}}) - P(T_{\mathrm{L}})$")
    apply_scaled_colorbar_labels(colorbar, vmax)

    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = "ree" if observable == "ree" else "rg"
    name = f"Figure1_thermodynamic_conformational_ensemble_{suffix}"
    for ext in ("png", "pdf"):
        fig.savefig(output_dir / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {output_dir / name}.png")
    print(f"[save] {output_dir / name}.pdf")


def build_combined_figure(output_dir: Path) -> None:
    """Build the final three-panel figure: Rg, Ree, and contact-map changes."""
    delta_maps = {}
    delta_abs = []
    for system in SYSTEM_ORDER:
        for family in FAMILY_ORDER:
            delta = temperature_delta_contact_map(system, family)
            delta_maps[(system, family)] = delta
            delta_abs.extend(np.abs(delta[np.isfinite(delta)]).ravel())
    vmax = nice_symmetric_limit(float(np.quantile(np.asarray(delta_abs), 0.995)))

    # Four panels: A (Rg), B (Ree), C (heating-response summary), D (contact maps).
    fig = plt.figure(figsize=(17.0, 14.4), dpi=300)
    outer = GridSpec(
        4,
        1,
        figure=fig,
        height_ratios=[0.95, 0.95, 1.18, 1.95],
        hspace=0.54,
        left=0.060,
        right=0.945,
        top=0.955,
        bottom=0.055,
    )

    gs_a = outer[0].subgridspec(1, 3, wspace=0.10)
    gs_b = outer[1].subgridspec(1, 3, wspace=0.10)
    gs_csum = outer[2].subgridspec(1, 2, wspace=0.12)
    # No colorbar column here: all three contact maps receive the full width.
    gs_d = outer[3].subgridspec(1, 3, wspace=0.07)

    first_axes = {}
    summary_axes = []

    for row, (observable, grid) in enumerate((("rg", gs_a), ("ree", gs_b))):
        for system_index, system in enumerate(SYSTEM_ORDER):
            ax = fig.add_subplot(grid[0, system_index])
            if system_index == 0:
                first_axes["A" if row == 0 else "B"] = ax
            plot_size_distribution(ax, system, observable, show_ylabel=(system_index == 0))

            ax.set_ylabel("Probability density" if system_index == 0 else None, fontsize=16)

            ax.tick_params(axis="both", labelsize=15)
            ax.xaxis.label.set_size(20)
            if row == 0 and system_index == 0:
                ax.legend(
                    loc="upper right", bbox_to_anchor=(0.985, 0.985), ncol=2,
                    fontsize=12.5, frameon=False, columnspacing=0.7,
                    handlelength=1.55, handletextpad=0.35, borderaxespad=0.0,
                )

    # Panel C: mean-based temperature response summaries for Rg and Ree.
    response_values = []
    for response_observable in ("rg", "ree"):
        for system in SYSTEM_ORDER:
            for family in FAMILY_ORDER:
                response_values.extend(response_bootstrap_interval(system, family, response_observable))
    resp_min = min(response_values)
    resp_max = max(response_values)
    lower = min(-1.0, np.floor(resp_min - 0.8))
    # Reserve clear headroom for the in-panel observable labels and legend.
    upper = max(26.0, np.ceil(resp_max + 2.0))
    summary_ylim = (lower, upper)

    for summary_index, response_observable in enumerate(("rg", "ree")):
        ax = fig.add_subplot(gs_csum[0, summary_index])
        summary_axes.append(ax)
        if summary_index == 0:
            first_axes["C"] = ax
        plot_response_summary(ax, response_observable, ylim=summary_ylim)
        ax.set_ylabel("Change upon heating" if summary_index == 0 else None, fontsize=16)
        response_title = (
            r"$\mathbfit{R}_{\mathbfit{g}}$"
            if response_observable == "rg"
            else r"$\mathbfit{R}_{\mathbfit{ee}}$"
        )
        ax.set_title(
            response_title,
            fontsize=21,
            fontweight="bold",
        )
        ax.tick_params(axis="both", labelsize=15)
        if summary_index == 0:
            ax.legend(
                handles=response_legend_handles(),
                labels=["Free", "Condensate"],
                loc="upper right",
                ncol=2,
                frameon=False,
                fontsize=14,
                handletextpad=0.5,
            )
    image = None
    contact_axes = []
    for system_index, system in enumerate(SYSTEM_ORDER):
        ax = fig.add_subplot(gs_d[0, system_index])
        if system_index == 0:
            first_axes["D"] = ax
        contact_axes.append(ax)

        free_delta = delta_maps[(system, "free")]
        condensate_delta = delta_maps[(system, "condensate")]
        n_residues = min(free_delta.shape[0], condensate_delta.shape[0])
        free_delta = free_delta[:n_residues, :n_residues]
        condensate_delta = condensate_delta[:n_residues, :n_residues]
        row_index, column_index = np.indices((n_residues, n_residues))
        # With origin="lower", row >= column is the visually upper-left
        # triangle and row < column is the visually lower-right triangle.
        combined_delta = np.where(row_index >= column_index, free_delta, condensate_delta)

        image = ax.imshow(
            combined_delta,
            origin="lower",
            cmap=contact_cmap(),
            norm=TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax),
            interpolation="nearest",
            aspect="equal",
        )
        ax.set_title(display_system(system), fontsize=19, fontweight="bold", pad=9)
        ax.set_ylabel("Residue index" if system_index == 0 else None, fontsize=20)
        ax.set_xlabel("Residue index", fontsize=20)
        ax.text(
            0.035, 0.955, "Free", transform=ax.transAxes,
            ha="left", va="top", fontsize=17, fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.68, "pad": 1.5},
        )
        ax.text(
            0.965, 0.045, "Condensate", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=17, fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.68, "pad": 1.5},
        )
        ax.xaxis.set_major_locator(MultipleLocator(25))
        ax.yaxis.set_major_locator(MultipleLocator(25))
        ax.tick_params(axis="both", labelsize=14)


    # Place panel letters from the actual axes positions rather than fixed
    # page coordinates. This keeps them above each panel and away from y labels.
    fig.canvas.draw()
    for label in ("A", "B", "C", "D"):
        box = first_axes[label].get_position()
        fig.text(
            0.012,
            min(box.y1 + 0.018, 0.988),
            label,
            fontsize=24,
            fontweight="bold",
            ha="left",
            va="bottom",
        )

    # Floating colorbar matched exactly to the merged contact-map row.
    contact_boxes = [ax.get_position() for ax in contact_axes]
    c_bottom = min(box.y0 for box in contact_boxes)
    c_top = max(box.y1 for box in contact_boxes)
    c_right = max(box.x1 for box in contact_boxes)
    cax = fig.add_axes([c_right + 0.012, c_bottom, 0.014, c_top - c_bottom])
    colorbar = fig.colorbar(image, cax=cax, ticks=colorbar_ticks_for_limit(vmax))
    colorbar.set_label(
        r"$\Delta P_{\mathrm{intra}}$",
        fontsize=17,
        labelpad=12,
    )
    apply_scaled_colorbar_labels(colorbar, vmax)
    colorbar.ax.tick_params(labelsize=14)

    output_dir.mkdir(parents=True, exist_ok=True)
    name = "Figure1_thermodynamic_conformational_ensemble_v17"
    for ext in ("png", "pdf"):
        fig.savefig(output_dir / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {output_dir / name}.png")
    print(f"[save] {output_dir / name}.pdf")

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--observable", choices=("combined", "rg", "ree"), default="combined")
    args = parser.parse_args()

    setup_style()
    if args.observable == "combined":
        build_combined_figure(args.output_dir)
    else:
        build_figure(args.output_dir, args.observable)


if __name__ == "__main__":
    main()

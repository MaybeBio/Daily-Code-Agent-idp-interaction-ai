#!/usr/bin/env python3
"""Compact data-driven TOC based on the user-selected draft layout."""

from pathlib import Path
import csv
import os

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from scipy.stats import gaussian_kde
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(os.environ.get("JPCL_DATA_ROOT", REPOSITORY_ROOT)).resolve()
OUT = REPOSITORY_ROOT / "outputs/toc"
TL = "#2f7f73"
TH = "#b54b5f"
INK = "#25292e"
LIGHT = "#cbd0d4"


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


SNAPSHOTS = {
    "Tl": DATA_ROOT / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet_supplement/supp_figure8_representative_structures/representative_pdbs/condensate_280K_s2_Cond.pdb",
    "Th": DATA_ROOT / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet_supplement/supp_figure8_representative_structures/representative_pdbs/condensate_325K_s2_Cond.pdb",
}


def snapshot_xy(path):
    xyz = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    xyz = np.asarray(xyz, dtype=float)
    xyz -= xyz.mean(axis=0)
    _, _, vh = np.linalg.svd(xyz, full_matrices=False)
    return xyz @ vh[:2].T


def condensate(ax, role, xy, shared_scale):
    color = TL if role == "Tl" else TH
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1.05, 1.05)
    ax.set_aspect("equal")
    ax.axis("off")
    # A soft field identifies the condensate environment without enclosing the
    # conformation in a heavy schematic boundary.
    ax.add_patch(Circle((0, 0), 0.82, fc=color, ec="none", alpha=0.045))
    # The highlighted chain is a real representative condensate conformation.
    chain = xy / shared_scale * 0.94
    ax.plot(chain[:, 0], chain[:, 1], color=color, lw=2.8, alpha=0.95, zorder=4)
    bead_idx = np.unique(np.r_[np.arange(0, len(chain), 7), len(chain) - 1])
    ax.scatter(chain[bead_idx, 0], chain[bead_idx, 1], s=11.5, color=color,
               edgecolor="white", linewidth=0.25, zorder=5)
    ax.scatter(chain[[0, -1], 0], chain[[0, -1], 1], s=27, color=color,
               edgecolor="white", linewidth=0.55, zorder=6)


def network_edges():
    path = DATA_ROOT / "project_figure/Figure3_transition_contribution_source.csv"
    selected = [r for r in rows(path) if r["system"] == "FUS" and r["family"] == "condensate"]
    output = {"Tl": [], "Th": []}
    for row in selected:
        i, j = int(row["source_state"]), int(row["target_state"])
        if i == j:
            continue
        for role, field in (("Tl", "transition_probability_tl"), ("Th", "transition_probability_th")):
            value = float(row[field])
            if value > 0.10:
                output[role].append((i, j, value))
    return output


def rendered_structure(ax, role):
    temp = "280" if role == "Tl" else "325"
    image_path = DATA_ROOT / "toc_assets" / f"toc_fus_{temp}K_condensate.png"
    if not image_path.exists():
        image_path = REPOSITORY_ROOT / "toc_assets" / image_path.name
    image = mpimg.imread(image_path)
    rgb = image[..., :3]
    # Crop by the saturated tagged chain, not by the full gray environment.
    # The nearby gray chain fragments remain visible inside this much tighter
    # field of view, while the colored chain nearly reaches the frame edges.
    if role == "Tl":
        mask = (rgb[..., 1] > rgb[..., 0] * 1.18) & (rgb[..., 1] > rgb[..., 2] * 1.08) & (rgb[..., 1] < 0.92)
    else:
        mask = (rgb[..., 0] > rgb[..., 1] * 1.18) & (rgb[..., 0] > rgb[..., 2] * 1.08) & (rgb[..., 0] < 0.94)
    yy, xx = np.where(mask)
    # Use the same wide rectangular crop for both temperatures.  A tight crop
    # lets the tagged chain and its local condensate environment fill the box.
    pad = 35 if role == "Tl" else 90
    y0, y1 = max(0, yy.min() - pad), min(image.shape[0], yy.max() + pad)
    x0, x1 = max(0, xx.min() - pad), min(image.shape[1], xx.max() + pad)
    cy, cx = (y0 + y1) // 2, (x0 + x1) // 2
    target_aspect = 1.77
    height, width = y1 - y0, x1 - x0
    if width / height < target_aspect:
        width = int(np.ceil(height * target_aspect))
    else:
        height = int(np.ceil(width / target_aspect))
    y0, y1 = max(0, cy - height // 2), min(image.shape[0], cy + (height + 1) // 2)
    x0, x1 = max(0, cx - width // 2), min(image.shape[1], cx + (width + 1) // 2)
    image = image[y0:y1, x0:x1]
    ax.imshow(image, interpolation="lanczos", aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.0)
        spine.set_edgecolor(INK)


def contact_map(role, display_state):
    source = DATA_ROOT / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet/figure3_state_population_contactmaps/figure3_condensate_state_population_contactmaps_source.csv"
    temp = "280K" if role == "Tl" else "325K"
    match = next(r for r in rows(source)
                 if r["temperature"] == temp and int(r["display_state"]) == display_state)
    local_state = int(match["local_state"])
    avg_path = DATA_ROOT / f"fus/fus_protein_single_simulation/condensate/tica_100dimension_10frames_{temp}/avg_contacts.npy"
    data = np.asarray(np.load(avg_path, mmap_mode="r")[local_state], dtype=float)
    idx = np.arange(data.shape[0])
    data[np.abs(idx[:, None] - idx[None, :]) <= 4] = np.nan
    return np.clip(data, 1e-4, None), float(match["mu"])


def rg_distribution(ax):
    base = DATA_ROOT / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet/figure1_rg_ree_relaxation"
    all_rg = []
    for role, temp, color in (("Tl", "280K", TL), ("Th", "325K", TH)):
        values = np.asarray(np.load(base / f"cache_condensate_{temp}.npz")["rg"], dtype=float)
        values = values[::max(1, values.size // 50000)]
        all_rg.append(values)
    plot_lo, plot_hi = 1.45, 7.05
    grid = np.linspace(plot_lo, plot_hi, 240)
    for values, color in zip(all_rg, (TL, TH)):
        kde = gaussian_kde(values[::max(1, values.size // 30000)], bw_method=0.18)
        density = kde(grid)
        ax.plot(grid, density, color=color, lw=1.75, solid_capstyle="round")
    ax.set_xlim(plot_lo, plot_hi)
    ax.set_ylim(bottom=0)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.spines[:].set_visible(False)
    ax.add_patch(FancyArrowPatch(
        (0, 0), (1.04, 0), transform=ax.transAxes, shrinkA=0, shrinkB=0,
        arrowstyle="-|>", mutation_scale=7, color=INK, lw=0.75,
        clip_on=False, zorder=10))
    ax.add_patch(FancyArrowPatch(
        (0, 0), (0, 1.05), transform=ax.transAxes, shrinkA=0, shrinkB=0,
        arrowstyle="-|>", mutation_scale=7, color=INK, lw=0.75,
        clip_on=False, zorder=10))
    ax.text(0.0, 1.075, r"$P$", transform=ax.transAxes, fontsize=8.5,
            ha="center", va="bottom")
    ax.text(1.00, -0.085, r"$R_g\,(R_{ee})$", transform=ax.transAxes, fontsize=8.0,
            ha="right", va="top")


def network(ax, role, edges):
    color = TL if role == "Tl" else TH
    pos = {0: np.array([0.50, 0.82]), 1: np.array([0.12, 0.12]), 2: np.array([0.88, 0.12])}
    ax.set_xlim(-0.16, 1.18)
    ax.set_ylim(-0.18, 1.12)
    ax.axis("off")
    state_maps = {}
    populations = {}
    maps = {}
    for state in pos:
        cmap, population = contact_map(role, state)
        maps[state] = cmap
        populations[state] = population
    max_population = max(populations.values())
    for state in pos:
        cmap, population = maps[state], populations[state]
        # Encode population directly in the contact-map side length.  This is
        # deliberately stronger than area scaling so population changes remain
        # legible at TOC size.
        size = max(0.18, 0.54 * population / max_population)
        state_maps[state] = (cmap, population, size)
    for i, j, probability in edges[role]:
        pair = tuple(sorted((i, j)))
        canonical = pos[pair[1]] - pos[pair[0]]
        unit = canonical / np.linalg.norm(canonical)
        offset_a = 0.5 * state_maps[pair[0]][2] / np.max(np.abs(unit)) + 0.018
        offset_b = 0.5 * state_maps[pair[1]][2] / np.max(np.abs(unit)) + 0.018
        base_a = pos[pair[0]] + offset_a * unit
        base_b = pos[pair[1]] - offset_b * unit
        # Diagonal reverse paths are separated only in x, so their start and
        # arrival y coordinates match exactly.  The horizontal pair is
        # separated only in y.
        # Keep opposing directions visibly separated at TOC size.
        lane = np.array([0.042, 0.0]) if pair != (1, 2) else np.array([0.0, 0.042])
        if role == "Th" and pair == (0, 2):
            # Use the right edge of the upper map and the top edge of the
            # lower-right map.  These endpoints create a long, open route on
            # the outside of the network rather than a short corner-to-corner
            # connection between two large maps.
            upper_half = 0.5 * state_maps[0][2]
            lower_half = 0.5 * state_maps[2][2]
            route_a = pos[0] + np.array([upper_half + 0.065, -0.045])
            route_b = pos[2] + np.array([0.085, lower_half + 0.065])
            # Separate the two directions only along x.  Their endpoints
            # therefore share identical y coordinates at both maps.
            route_lane = np.array([0.035, 0.0])
            if i == 0:
                start, finish = route_a + route_lane, route_b + route_lane
            else:
                start, finish = route_b - route_lane, route_a - route_lane
        elif i == pair[0]:
            start, finish = base_a + lane, base_b + lane
        else:
            start, finish = base_b - lane, base_a - lane
        arrow_width = max(1.35, 0.85 * (probability / 0.10) ** 2)
        head_width = max(5.5, 1.55 * arrow_width)
        head_length = max(6.5, 1.90 * arrow_width)
        arrow_style = (
            f"Simple,tail_width={arrow_width:.3f},"
            f"head_width={head_width:.3f},head_length={head_length:.3f}"
        )
        ax.add_patch(FancyArrowPatch(
            start, finish,
            arrowstyle=arrow_style, mutation_scale=1.0,
            linewidth=0.0, facecolor=color, edgecolor=color, alpha=0.92,
            connectionstyle="arc3,rad=0",
        ))
    # Contact-map side length encodes state population.
    for state, point in pos.items():
        cmap, population, size = state_maps[state]
        x0, y0 = point - size / 2
        ax.imshow(cmap, origin="lower", cmap="viridis",
                  norm=LogNorm(vmin=1e-4, vmax=1e-1), interpolation="nearest",
                  extent=(x0, x0 + size, y0, y0 + size), zorder=5)
        ax.add_patch(Rectangle((x0, y0), size, size, fc="none", ec=color,
                               lw=1.4, zorder=6))


def measured_acf(kind, role):
    if kind == "intra":
        temp = "280K" if role == "Tl" else "325K"
        path = DATA_ROOT / f"project_figure_set_Tl_Th/source_data/intrachain_relaxation/intrachain_FUS_condensate_{temp}.npz"
        data = np.load(path, allow_pickle=False)
        return np.asarray(data["time_ns"]), np.asarray(data["acf_mean"])
    temp = "280K" if role == "Tl" else "325K"
    path = DATA_ROOT / f"project_figure/source_data/interchain_ratio_bootstrap/interchain_contact_relaxation_FUS_{temp}_parts081-150.npz"
    data = np.load(path, allow_pickle=False)
    return np.asarray(data["time_ns"]), np.asarray(data["acf_mean"])


def acf_axis(ax, kind, title):
    for role, color in (("Tl", TL), ("Th", TH)):
        time, mean = measured_acf(kind, role)
        ax.plot(time, mean, color=color, lw=2.2)
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 1.03)
    ax.set_xticks([0, 5, 10, 15])
    ax.set_yticks([0, 1], ["0", "1"])
    ax.tick_params(axis="both", labelsize=8.5)
    ax.tick_params(axis="y", length=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.1)
    ax.set_ylabel("Autocorrelation", fontsize=8.2, rotation=90, labelpad=7)
    ax.set_xlabel("Lag time", fontsize=8.2, labelpad=2)


def main():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.linewidth": 1.0,
        "mathtext.fontset": "dejavusans",
    })
    fig = plt.figure(figsize=(7.2, 4.15), dpi=300, facecolor="white")

    # Real representative FUS-LCD conformations rendered with VMD.
    fig.text(0.335, 0.912, "Global structural buffering", ha="center", va="center",
             fontsize=10.7, fontweight="bold", color=INK)
    ax_sl = fig.add_axes([0.0075, 0.565, 0.235, 0.225])
    ax_rg = fig.add_axes([0.273, 0.627, 0.130, 0.110])
    ax_sh = fig.add_axes([0.4275, 0.565, 0.235, 0.225])
    rendered_structure(ax_sl, "Tl")
    rg_distribution(ax_rg)
    rendered_structure(ax_sh, "Th")
    fig.text(0.125, 0.822, r"Low $\mathbf{\mathit{T}}$", color=TL, fontsize=10.2,
             fontweight="bold", ha="center", va="center")
    fig.text(0.545, 0.822, r"High $\mathbf{\mathit{T}}$", color=TH, fontsize=10.2,
             fontweight="bold", ha="center", va="center")
    # Actual FUS-condensate transition networks.
    edges = network_edges()
    ax_nl = fig.add_axes([-0.0225, 0.072, 0.295, 0.385])
    ax_nh = fig.add_axes([0.3975, 0.072, 0.295, 0.385])
    network(ax_nl, "Tl", edges)
    network(ax_nh, "Th", edges)
    fig.text(0.335, 0.037, "Microscopic conformational remodeling", ha="center", va="center",
             fontsize=10.2, fontweight="bold", color=INK)
    # Each temperature column links its condensate snapshot to the associated
    # conformational network with a vertical downward arrow.
    fig.patches.append(FancyArrowPatch(
        (0.125, 0.557), (0.125, 0.475), transform=fig.transFigure,
        arrowstyle="-|>", mutation_scale=10, color=INK, lw=1.1, clip_on=False))
    fig.patches.append(FancyArrowPatch(
        (0.545, 0.557), (0.545, 0.475), transform=fig.transFigure,
        arrowstyle="-|>", mutation_scale=10, color=INK, lw=1.1, clip_on=False))
    # Divider and the two measured ACF responses.
    fig.lines.append(Line2D([0.690, 0.690], [0.0, 1.0], transform=fig.transFigure,
                            color=LIGHT, lw=1.0, ls=(0, (1.5, 2.5))))
    fig.text(0.870, 0.900, "Intrachain contact", ha="center", va="center",
             fontsize=10.7, fontweight="bold", color=INK)
    fig.text(0.870, 0.435, "Interchain contact", ha="center", va="center",
             fontsize=10.7, fontweight="bold", color=INK)
    ax_intra = fig.add_axes([0.765, 0.555, 0.210, 0.290])
    ax_inter = fig.add_axes([0.765, 0.105, 0.210, 0.290])
    acf_axis(ax_intra, "intra", "Intrachain contact number")
    acf_axis(ax_inter, "inter", "Interchain contact number")

    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "TOC_graphic_data_v3.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT / "TOC_graphic_data_v3.png", dpi=900, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


if __name__ == "__main__":
    main()

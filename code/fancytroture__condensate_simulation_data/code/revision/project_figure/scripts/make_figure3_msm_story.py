#!/usr/bin/env python3
"""Build Figure 3: Markov-state model and state-dynamics story."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "ppt_5_11_figures/results"
PAPER = PROJECT_ROOT / "ppt_5_11_figures/paper_story_figures"

SYSTEM_ORDER = ("FUS", "A1", "TDP43")
DISPLAY_NAMES = {"FUS": "FUS-LCD", "A1": "A1-LCD", "TDP43": "TDP-43-LCD"}
FAMILY_ORDER = ("free", "condensate")
ROLE_ORDER = ("Tl", "Th")
TEMP_COLORS = {"Tl": "#2f7f73", "Th": "#b54b5f"}
FAMILY_COLORS = {"free": "#4a4a4a", "condensate": "#4a4a4a"}
FAMILY_LABELS = {"free": "Free", "condensate": "Condensate"}
STATE_COLORS = {
    0: "#3b82b8",
    1: "#f28e2b",
    2: "#4daf4a",
    3: "#8f63a8",
    4: "#d95f5f",
    5: "#9c6b5d",
    6: "#d77fb3",
    7: "#7f7f7f",
}
STATE_LAYOUT_TEMPLATES = {
    # Make each motif visibly larger while preserving open outer lanes for arrows.
    3: [(0.50, 0.84), (0.23, 0.24), (0.77, 0.24)],
    4: [(0.50, 0.85), (0.23, 0.54), (0.50, 0.16), (0.77, 0.54)],
    5: [(0.50, 0.86), (0.21, 0.65), (0.28, 0.16), (0.72, 0.16), (0.79, 0.65)],
}
THUMB_SIZES = {3: 0.32, 4: 0.295, 5: 0.255}
TEMPS = {
    "FUS": {"Tl": "280K", "Th": "325K"},
    "A1": {"Tl": "230K", "Th": "260K"},
    "TDP43": {"Tl": "260K", "Th": "300K"},
}
TICA_IMAGES = {
    "FUS": RESULTS / "tica_state_projection_fus_free_condensate.png",
    "A1": RESULTS / "tica_state_projection_a1_free_condensate.png",
    "TDP43": RESULTS / "tica_state_projection_tdp43_free_condensate.png",
}
CONTACTMAP_IMAGES = {
    "FUS": PROJECT_ROOT / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet/figure3_state_population_contactmaps/figure3_condensate_state_population_contactmaps.png",
    "A1": PROJECT_ROOT / "A1/A1_protein_single_simulation/paper_thermo_figures_a1_4p17_relaxed_mu10/figure3_state_population_contactmaps/figure3_condensate_state_population_contactmaps.png",
    "TDP43": PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/paper_thermo_figures_tdp43_4p17_relaxed_mu10/figure3_state_population_contactmaps/figure3_condensate_state_population_contactmaps.png",
}
AVG_CONTACTS = {
    "FUS": {
        "Tl": PROJECT_ROOT / "fus/fus_protein_single_simulation/condensate/tica_100dimension_10frames_280K/avg_contacts.npy",
        "Th": PROJECT_ROOT / "fus/fus_protein_single_simulation/condensate/tica_100dimension_10frames_325K/avg_contacts.npy",
    },
    "A1": {
        "Tl": PROJECT_ROOT / "A1/A1_protein_single_simulation/condensate/tica_100dimension_10frames_230K/avg_contacts.npy",
        "Th": PROJECT_ROOT / "A1/A1_protein_single_simulation/condensate/tica_100dimension_10frames_260K/avg_contacts.npy",
    },
    "TDP43": {
        "Tl": PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/condensate/tica_100dimension_10frames_260K/avg_contacts.npy",
        "Th": PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/condensate/tica_100dimension_10frames_300K/avg_contacts.npy",
    },
}
CONTACT_FINGERPRINT_STATES = {
    "FUS": {"Tl": 0, "Th": 2},
    "A1": {"Tl": 0, "Th": 5},
    "TDP43": {"Tl": 0, "Th": 5},
}
STATE_SOURCE = {
    "FUS": {
        "free": PROJECT_ROOT / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet/figure3_state_population_contactmaps/figure3_free_state_population_contactmaps_source.csv",
        "condensate": PROJECT_ROOT / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet/figure3_state_population_contactmaps/figure3_condensate_state_population_contactmaps_source.csv",
    },
    "A1": {
        "free": RESULTS / "a1_230_260_state_mapping/figure3_free_state_population_contactmaps_source.csv",
        "condensate": RESULTS / "a1_230_260_state_mapping/figure3_condensate_state_population_contactmaps_source.csv",
    },
    "TDP43": {
        "free": PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/paper_thermo_figures_tdp43_4p17_relaxed_mu10/figure3_state_population_contactmaps/figure3_free_state_population_contactmaps_source.csv",
        "condensate": PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/paper_thermo_figures_tdp43_4p17_relaxed_mu10/figure3_state_population_contactmaps/figure3_condensate_state_population_contactmaps_source.csv",
    },
}
MFPT_SOURCE = {
    "FUS": {
        "free": PROJECT_ROOT / "fus/fus_protein_single_simulation/figs_final_mfpt_crossT_global/mfpt_summary_combined_n5_lag50.csv",
        "condensate": PROJECT_ROOT / "fus/fus_protein_single_simulation/condensate/figs_final_mfpt_crossT_global/mfpt_summary_combined_n5_lag50.csv",
    },
    "A1": {
        "free": PROJECT_ROOT / "A1/A1_protein_single_simulation/figs_final_mfpt_crossT_global/mfpt_summary_combined_n5_lag50.csv",
        "condensate": PROJECT_ROOT / "A1/A1_protein_single_simulation/condensate/figs_final_mfpt_crossT_global/mfpt_summary_combined_n5_lag50.csv",
    },
    "TDP43": {
        "free": PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/figs_final_mfpt_crossT_global/mfpt_summary_combined_n5_lag50.csv",
        "condensate": PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/condensate/figs_final_mfpt_crossT_global/mfpt_summary_combined_n5_lag50.csv",
    },
}
LIFETIME_SOURCE = PAPER / "figureS_state_lifetime_compact_source.csv"


def setup_style() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11.0,
            "axes.labelsize": 11.5,
            "xtick.labelsize": 10.0,
            "ytick.labelsize": 10.0,
            "legend.fontsize": 9.5,
            "legend.frameon": False,
            "axes.linewidth": 0.85,
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


def display_system(system: str) -> str:
    return DISPLAY_NAMES.get(system, system)


def display_systems() -> list[str]:
    return [display_system(system) for system in SYSTEM_ORDER]


def temp_label(role: str) -> str:
    return r"$T_{\mathrm{L}}$" if role == "Tl" else r"$T_{\mathrm{H}}$"


def row_temperature(row: dict[str, str]) -> str:
    return row.get("temperature", row.get("temp", ""))


def clean_axis(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#bdbdbd", alpha=0.25, linewidth=0.6)




def box_axis(ax) -> None:
    """Use a full rectangular frame, as requested for panel B."""
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
        spine.set_color("#333333")
    ax.grid(axis="y", color="#bdbdbd", alpha=0.22, linewidth=0.6)
    ax.tick_params(direction="out", length=3.5, width=0.8)

def crop_image(image: np.ndarray, pad: int = 12) -> np.ndarray:
    rgb = image[..., :3]
    alpha = image[..., 3] if image.shape[-1] == 4 else np.ones(image.shape[:2])
    mask = (alpha > 0.05) & np.any(rgb < 0.985, axis=2)
    if not np.any(mask):
        return image
    ys, xs = np.where(mask)
    return image[
        max(0, ys.min() - pad) : min(image.shape[0], ys.max() + pad),
        max(0, xs.min() - pad) : min(image.shape[1], xs.max() + pad),
    ]


def contactmap_state_order(system: str) -> list[int]:
    states = sorted(set(state_vector(system, "condensate", "Tl")) | set(state_vector(system, "condensate", "Th")))
    return states


def contactmap_fingerprint_panel(system: str, role: str) -> np.ndarray:
    image = mpimg.imread(CONTACTMAP_IMAGES[system])
    height, width = image.shape[:2]
    states = contactmap_state_order(system)
    state = CONTACT_FINGERPRINT_STATES[system][role]
    col = states.index(state)
    ncols = len(states)

    # These source figures have a population summary above two contact-map rows.
    if ncols == 3:
        x_start, x_step, panel_width = 0.195, 0.305, 0.180
    elif ncols == 4:
        x_start, x_step, panel_width = 0.130, 0.225, 0.172
    elif ncols == 5:
        x_start, x_step, panel_width = 0.117, 0.167, 0.135
    else:
        x_start = 0.13
        panel_width = 0.78 / max(ncols, 1)
        x_step = panel_width
    x0 = int(width * (x_start + x_step * col))
    x1 = int(width * (x_start + x_step * col + panel_width))
    if role == "Tl":
        y0, y1 = int(height * 0.335), int(height * 0.625)
    else:
        y0, y1 = int(height * 0.705), int(height * 0.982)
    panel = image[y0:y1, x0:x1]
    return crop_image(panel, pad=2)


def tica_condensate_panel(system: str, role: str) -> np.ndarray:
    image = mpimg.imread(TICA_IMAGES[system])
    height, width = image.shape[:2]
    row = 1
    col = 0 if role == "Tl" else 1
    y0 = row * (height // 2)
    y1 = height
    x0 = col * (width // 2)
    x1 = (col + 1) * (width // 2) if col == 0 else width
    panel = image[y0:y1, x0:x1]
    panel = panel[265 : panel.shape[0] - 24, 18 : panel.shape[1] - 18]
    return crop_image(panel, pad=10)


def state_vector(system: str, family: str, role: str) -> dict[int, float]:
    temp = TEMPS[system][role]
    vector: dict[int, float] = {}
    for row in read_csv(STATE_SOURCE[system][family]):
        if row_temperature(row) != temp:
            continue
        state = int(row.get("display_state") or row.get("global_state", "0"))
        vector[state] = vector.get(state, 0.0) + max(fval(row["mu"]), 0.0)
    total = sum(vector.values())
    if total <= 0.0:
        return {}
    return {state: value / total for state, value in vector.items()}


def local_state_for_display(system: str, family: str, role: str, display_state: int) -> int:
    temp = TEMPS[system][role]
    for row in read_csv(STATE_SOURCE[system][family]):
        if row_temperature(row) != temp:
            continue
        state = int(row.get("display_state") or row.get("global_state", "0"))
        if state == display_state:
            return int(row["local_state"])
    raise KeyError(f"No local state for {system} {family} {role} display state {display_state}")


def display_state_for_local(system: str, family: str, role: str, local_state: int) -> int | None:
    temp = TEMPS[system][role]
    for row in read_csv(STATE_SOURCE[system][family]):
        if row_temperature(row) != temp:
            continue
        if int(row["local_state"]) == int(local_state):
            return int(row.get("display_state") or row.get("global_state", "0"))
    return None


def row_display_states(row: dict[str, str]) -> tuple[int, int]:
    system = row["system"]
    family = row["family"]
    role = row["role"]
    source = display_state_for_local(system, family, role, int(row["source_local_state"]))
    target = display_state_for_local(system, family, role, int(row["target_local_state"]))
    if source is None:
        source = int(row["source_display_state"])
    if target is None:
        target = int(row["target_display_state"])
    return source, target


def mask_contact_map(contact_map: np.ndarray, gap: int = 4) -> np.ndarray:
    masked = np.array(contact_map, dtype=float, copy=True)
    n_res = masked.shape[0]
    for offset in range(gap + 1):
        idx = np.arange(n_res - offset)
        jdx = idx + offset
        masked[idx, jdx] = np.nan
        masked[jdx, idx] = np.nan
    return np.clip(masked, 1e-4, None)


def contact_fingerprint_map(system: str, role: str) -> np.ndarray:
    display_state = CONTACT_FINGERPRINT_STATES[system][role]
    return contact_map_for_state(system, role, display_state)


def contact_map_for_state(system: str, role: str, display_state: int) -> np.ndarray:
    local_state = local_state_for_display(system, "condensate", role, display_state)
    contacts = np.load(AVG_CONTACTS[system][role], mmap_mode="r")
    return mask_contact_map(np.asarray(contacts[local_state], dtype=float))


def js_distance(p_dict: dict[int, float], q_dict: dict[int, float]) -> float:
    states = sorted(set(p_dict) | set(q_dict))
    if not states:
        return float("nan")
    p = np.array([p_dict.get(state, 0.0) for state in states], dtype=float)
    q = np.array([q_dict.get(state, 0.0) for state in states], dtype=float)
    if p.sum() <= 0.0 or q.sum() <= 0.0:
        return float("nan")
    p /= p.sum()
    q /= q.sum()
    m = 0.5 * (p + q)

    def kl(a: np.ndarray, b: np.ndarray) -> float:
        mask = a > 0.0
        return float(np.sum(a[mask] * np.log2(a[mask] / b[mask])))

    return math.sqrt(max(0.0, 0.5 * kl(p, m) + 0.5 * kl(q, m)))


def state_redistribution(system: str, family: str) -> float:
    return js_distance(state_vector(system, family, "Tl"), state_vector(system, family, "Th"))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def transition_network_rows(system: str) -> list[dict[str, str]]:
    path = RESULTS / f"transition_networks_{system.lower()}_free_condensate_global_aligned_source.csv"
    return read_csv(path)


def transition_states(system: str) -> list[int]:
    rows = transition_network_rows(system)
    states_seen = set()
    for row in rows:
        source, target = row_display_states(row)
        states_seen.add(source)
        states_seen.add(target)
    states = sorted(states_seen)
    return states


def transition_positions(system: str) -> dict[int, tuple[float, float]]:
    states = transition_states(system)
    if system == "TDP43":
        preferred = {
            0: (0.25, 3.85),
            1: (-3.35, 1.55),
            2: (-2.35, -3.65),
            5: (1.25, -3.65),
            6: (3.45, 0.70),
        }
        return {state: preferred[state] for state in states if state in preferred}
    if len(states) == 1:
        return {states[0]: (0.0, 0.0)}
    angles = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, len(states), endpoint=False)
    return {state: (2.85 * np.cos(angle), 3.45 * np.sin(angle)) for state, angle in zip(states, angles)}


def transition_matrix(system: str, family: str, role: str, states: list[int]) -> np.ndarray:
    index = {state: i for i, state in enumerate(states)}
    matrix = np.zeros((len(states), len(states)), dtype=float)
    for row in transition_network_rows(system):
        if row["family"] != family or row["role"] != role:
            continue
        source, target = row_display_states(row)
        if source not in index or target not in index:
            continue
        matrix[index[source], index[target]] += fval(row["transition_probability"])
    return matrix


def weighted_transition_distance(system: str, family: str) -> float:
    states = transition_states(system)
    tl = transition_matrix(system, family, "Tl", states)
    th = transition_matrix(system, family, "Th", states)
    pi_tl = state_vector(system, family, "Tl")
    pi_th = state_vector(system, family, "Th")
    pi_bar = np.array([0.5 * (pi_tl.get(state, 0.0) + pi_th.get(state, 0.0)) for state in states], dtype=float)
    diff = th - tl
    return float(np.sqrt(np.sum(pi_bar[:, None] * diff * diff)))


def transition_contribution_matrix(system: str, family: str = "condensate") -> tuple[list[int], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    states = transition_states(system)
    tl = transition_matrix(system, family, "Tl", states)
    th = transition_matrix(system, family, "Th", states)
    pi_tl = state_vector(system, family, "Tl")
    pi_th = state_vector(system, family, "Th")
    pi_bar = np.array([0.5 * (pi_tl.get(state, 0.0) + pi_th.get(state, 0.0)) for state in states], dtype=float)
    delta = th - tl
    contribution = pi_bar[:, None] * delta * delta
    return states, contribution, delta, tl, th


def transition_contribution_source_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for system in SYSTEM_ORDER:
        states, contribution, delta, tl, th = transition_contribution_matrix(system, "condensate")
        pi_tl = state_vector(system, "condensate", "Tl")
        pi_th = state_vector(system, "condensate", "Th")
        for i, source_state in enumerate(states):
            pi_bar = 0.5 * (pi_tl.get(source_state, 0.0) + pi_th.get(source_state, 0.0))
            for j, target_state in enumerate(states):
                rows.append(
                    {
                        "system": system,
                        "family": "condensate",
                        "source_state": str(source_state),
                        "target_state": str(target_state),
                        "transition_probability_tl": f"{tl[i, j]:.8f}",
                        "transition_probability_th": f"{th[i, j]:.8f}",
                        "delta_transition_probability": f"{delta[i, j]:.8f}",
                        "pi_bar_source": f"{pi_bar:.8f}",
                        "weighted_contribution": f"{contribution[i, j]:.10f}",
                    }
                )
    return rows


def plot_transition_contribution_matrix(ax, system: str, vmax: float):
    states, contribution, _delta, _tl, _th = transition_contribution_matrix(system, "condensate")
    im = ax.imshow(contribution, cmap="magma", vmin=0.0, vmax=vmax, interpolation="nearest")
    labels = [rf"$s_{{{state + 1}}}$" for state in states]
    ax.set_xticks(np.arange(len(states)), labels)
    ax.set_yticks(np.arange(len(states)), labels)
    ax.tick_params(length=0, labelsize=8.5, pad=1)
    ax.set_title(display_system(system), fontsize=10.8, fontweight="bold", pad=4)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color("#333333")
    return im


def median_mfpt_ratio(system: str, family: str) -> float:
    rows = read_csv(MFPT_SOURCE[system][family])
    by_temp = {row["temp"]: row for row in rows}
    tl = fval(by_temp[TEMPS[system]["Tl"]]["median"])
    th = fval(by_temp[TEMPS[system]["Th"]]["median"])
    return th / tl


def lifetime_source_rows() -> list[dict[str, str]]:
    return read_csv(LIFETIME_SOURCE)


def lifetime_value(system: str, family: str, role: str, state: int) -> float | None:
    temp = TEMPS[system][role]
    for row in lifetime_source_rows():
        if row["system"] != system or row["family"] != family or row["role"] != role:
            continue
        if row["temperature"] != temp:
            continue
        if int(row["global_state"]) == int(state):
            return fval(row["tau_lag"])
    return None


def lifetime_states(system: str, family: str) -> list[int]:
    states: set[int] = set()
    for row in lifetime_source_rows():
        if row["system"] == system and row["family"] == family:
            states.add(int(row["global_state"]))
    return sorted(states)


def plot_lifetime_family_panel(ax, family: str, show_ylabel: bool = False, show_legend: bool = False) -> None:
    xticks: list[float] = []
    xticklabels: list[str] = []
    x_centers: list[float] = []
    current_x = 0.0
    all_values: list[float] = []
    labeled_tl = False
    labeled_th = False
    for system in SYSTEM_ORDER:
        states = lifetime_states(system, family)
        group_positions: list[float] = []
        for state in states:
            x_pos = current_x
            group_positions.append(x_pos)
            tl = lifetime_value(system, family, "Tl", state)
            th = lifetime_value(system, family, "Th", state)
            if tl is not None and th is not None:
                ax.plot([x_pos, x_pos], [tl, th], color="#9a9a9a", lw=1.15, alpha=0.75, zorder=1)
            if tl is not None:
                all_values.append(tl)
                ax.scatter(
                    x_pos - 0.045,
                    tl,
                    s=58,
                    facecolor=TEMP_COLORS["Tl"],
                    edgecolor="#333333",
                    linewidth=0.7,
                    marker="o",
                    label=r"$T_{\mathrm{L}}$" if show_legend and not labeled_tl else None,
                    zorder=4,
                )
                labeled_tl = True
            if th is not None:
                all_values.append(th)
                ax.scatter(
                    x_pos + 0.045,
                    th,
                    s=58,
                    facecolor="#bd5268",
                    edgecolor="#333333",
                    linewidth=0.7,
                    marker="s",
                    label=r"$T_{\mathrm{H}}$" if show_legend and not labeled_th else None,
                    zorder=4,
                )
                labeled_th = True
            xticks.append(x_pos)
            xticklabels.append(rf"$s_{{{state + 1}}}$")
            current_x += 0.46
        if group_positions:
            x_centers.append(0.5 * (group_positions[0] + group_positions[-1]))
        current_x += 0.56

    ax.set_xticks(xticks, xticklabels)
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_title(FAMILY_LABELS[family], fontsize=11.5, pad=4)
    if show_ylabel:
        ax.set_ylabel(r"State lifetime, $\tau_i$ (lag units)")
    for x_center, system in zip(x_centers, SYSTEM_ORDER):
        ax.text(
            x_center,
            -0.25,
            display_system(system),
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=9.5,
            fontweight="bold",
            clip_on=False,
        )
    clean_axis(ax)
    if show_legend:
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles))
        ax.legend(unique.values(), unique.keys(), loc="upper right", ncol=2, handletextpad=0.35, columnspacing=0.7)


def plot_transition_network(ax, system: str, family: str, role: str, positions: dict[int, tuple[float, float]]) -> None:
    rows = [row for row in transition_network_rows(system) if row["family"] == family and row["role"] == role]
    states = sorted(
        set(row_display_states(row)[0] for row in rows)
        | set(row_display_states(row)[1] for row in rows)
    )
    populations = state_vector(system, family, role)
    strongest = [
        row
        for row in sorted(rows, key=lambda row: fval(row["transition_probability"]), reverse=True)
        if row_display_states(row)[0] != row_display_states(row)[1]
    ][:9]
    for row in strongest:
        source, target = row_display_states(row)
        probability = fval(row["transition_probability"])
        start = np.asarray(positions[source])
        end = np.asarray(positions[target])
        rad = 0.18 if source < target else -0.18
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=12,
                linewidth=0.9 + 4.1 * probability,
                color="#555555",
                alpha=0.72,
                connectionstyle=f"arc3,rad={rad}",
                shrinkA=23,
                shrinkB=23,
                zorder=2,
            )
        )
    for state in states:
        position = positions[state]
        population = populations.get(state, 0.0)
        radius = 0.27 + 0.31 * math.sqrt(max(population, 0.0))
        ax.add_patch(
            Circle(
                position,
                radius,
                facecolor=STATE_COLORS.get(state, "#999999"),
                edgecolor="#333333",
                linewidth=1.0,
                zorder=4,
            )
        )
        ax.text(position[0], position[1], rf"$s_{{{state + 1}}}$", color="white", ha="center", va="center", fontsize=11.5, fontweight="bold", zorder=5)
    ax.set_xlim(-1.42, 1.42)
    ax.set_ylim(-1.34, 1.34)
    ax.set_aspect("equal")
    ax.axis("off")


def data_to_axes_fraction(ax, x: float, y: float) -> tuple[float, float]:
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    return (x - x0) / (x1 - x0), (y - y0) / (y1 - y0)


def ordered_states_for_layout(states: list[int]) -> list[int]:
    states = sorted(states)
    if 0 in states:
        rest = [state for state in states if state != 0]
        if len(states) == 3:
            return [0] + rest
        if len(states) == 4:
            # top, left, bottom, right
            return [0, rest[0], rest[1], rest[2]]
        if len(states) == 5:
            # top, upper-left, lower-left, lower-right, upper-right
            return [0, rest[0], rest[1], rest[2], rest[3]]
    return states


def structural_axes_positions(states: list[int]) -> dict[int, tuple[float, float]]:
    ordered = ordered_states_for_layout(states)
    template = STATE_LAYOUT_TEMPLATES.get(len(ordered))
    if template is None:
        angles = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, len(ordered), endpoint=False)
        template = [(0.50 + 0.28 * np.cos(a), 0.50 + 0.28 * np.sin(a)) for a in angles]
    return {state: template[i] for i, state in enumerate(ordered)}


def plot_structural_transition_network(ax, system: str, role: str, positions: dict[int, tuple[float, float]]) -> None:
    family = "condensate"
    rows = [row for row in transition_network_rows(system) if row["family"] == family and row["role"] == role]
    states = sorted(
        set(row_display_states(row)[0] for row in rows)
        | set(row_display_states(row)[1] for row in rows)
    )
    populations = state_vector(system, family, role)

    # Larger contact maps and more open state motifs. Keep extra space below
    # the lower states so exterior arrows and state labels do not run into the
    # next protein row.
    if len(states) == 3:
        ordered = ordered_states_for_layout(states)
        coords = [(0.50, 0.84), (0.18, 0.26), (0.82, 0.26)]
        thumb_size = 0.355   # ≈1.25× the v11 three-state thumbnails
    elif len(states) == 4:
        ordered = ordered_states_for_layout(states)
        coords = [(0.50, 0.86), (0.16, 0.57), (0.50, 0.24), (0.84, 0.57)]
        thumb_size = 0.292   # ≈1.25× the v11 four-state thumbnails
    else:
        ordered = ordered_states_for_layout(states)
        angles = np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, len(ordered), endpoint=False)
        coords = [(0.50 + 0.36 * np.cos(a), 0.55 + 0.33 * np.sin(a)) for a in angles]
        thumb_size = 0.27

    # Keep both temperature columns centered in their own axes so the low- and high-temperature labels
    # are aligned across all protein rows.
    x_shift = 0.0
    state_positions = {
        state: (float(np.clip(coords[i][0] + x_shift, 0.08, 0.92)), coords[i][1])
        for i, state in enumerate(ordered)
    }

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")

    def center(state: int) -> np.ndarray:
        return np.asarray(state_positions[state], dtype=float)

    def normalize(v: np.ndarray) -> np.ndarray:
        n = float(np.hypot(v[0], v[1]))
        return v / n if n > 1e-12 else np.zeros(2, dtype=float)

    def square_edge_point(c: np.ndarray, direction: np.ndarray, gap: float = 0.018) -> np.ndarray:
        u = normalize(direction)
        half = 0.5 * thumb_size
        d = half / max(abs(u[0]), abs(u[1]), 1e-12)
        return c + (d + gap) * u

    panel_center = np.mean(np.asarray(list(state_positions.values()), dtype=float), axis=0)

    visible_rows = [
        row for row in rows
        if row_display_states(row)[0] != row_display_states(row)[1]
        and fval(row["transition_probability"]) >= 0.008
    ]
    best_by_direction: dict[tuple[int, int], dict[str, str]] = {}
    for row in visible_rows:
        key = row_display_states(row)
        if key not in best_by_direction or fval(row["transition_probability"]) > fval(best_by_direction[key]["transition_probability"]):
            best_by_direction[key] = row

    # Transition edges and edge labels are intentionally omitted in this draft.

    for state in states:
        x_frac, y_frac = state_positions[state]
        x0 = x_frac - thumb_size / 2
        x1 = x_frac + thumb_size / 2
        y0 = y_frac - thumb_size / 2
        y1 = y_frac + thumb_size / 2
        ax.imshow(
            contact_map_for_state(system, role, state),
            origin="lower", cmap="viridis",
            norm=LogNorm(vmin=1e-4, vmax=1e-1), interpolation="nearest",
            extent=(x0, x1, y0, y1),
            transform=ax.transAxes,
            clip_on=False,
            zorder=4,
        )
        ax.add_patch(
            Rectangle(
                (x0, y0),
                thumb_size,
                thumb_size,
                transform=ax.transAxes,
                facecolor="none",
                edgecolor="#2f2f2f",
                linewidth=1.05,
                clip_on=False,
                zorder=5,
            )
        )
        population = populations.get(state, 0.0)
        ax.text(
            x_frac, y0 - 0.018, rf"$s_{{{state + 1}}}$, {100.0 * population:.0f}%",
            transform=ax.transAxes, ha="center", va="top",
            fontsize=11.0, fontweight="bold", color="#222222",
            zorder=8, clip_on=False,
        )

def plot_state_occupancy_fingerprint(ax, system: str) -> None:
    states = sorted(set(state_vector(system, "condensate", "Tl")) | set(state_vector(system, "condensate", "Th")))
    y_positions = np.arange(len(ROLE_ORDER), dtype=float)
    for y, role in zip(y_positions, ROLE_ORDER):
        left = 0.0
        vector = state_vector(system, "condensate", role)
        for state in states:
            value = 100.0 * vector.get(state, 0.0)
            ax.barh(
                y,
                value,
                left=left,
                height=0.46,
                color=STATE_COLORS.get(state, "#999999"),
                edgecolor="white",
                linewidth=0.7,
            )
            if value >= 8.0:
                ax.text(left + value / 2.0, y, rf"$s_{{{state + 1}}}$", ha="center", va="center", color="white", fontsize=8.5, fontweight="bold")
            left += value
    ax.set_xlim(0, 100)
    ax.set_yticks(y_positions, [temp_label(role) for role in ROLE_ORDER])
    ax.invert_yaxis()
    ax.set_xlabel("Condensate state population (%)")
    ax.set_title(display_system(system), fontsize=13, fontweight="bold", pad=4)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#bdbdbd", alpha=0.25, linewidth=0.6)
    ax.tick_params(axis="y", length=0)


def plot_contact_fingerprint(ax, system: str, role: str) -> None:
    contact_map = contact_fingerprint_map(system, role)
    ax.imshow(contact_map, origin="lower", cmap="viridis", norm=LogNorm(vmin=1e-4, vmax=1e-1), interpolation="nearest")
    state = CONTACT_FINGERPRINT_STATES[system][role]
    vector = state_vector(system, "condensate", role)
    population = 100.0 * vector.get(state, 0.0)
    ax.set_title(rf"{temp_label(role)}  $s_{{{state + 1}}}$, {population:.1f}%", fontsize=10.5, pad=3)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)


def add_panel_labels(fig) -> None:
    fig.text(0.045, 0.955, "A", fontsize=18, fontweight="bold", ha="left", va="top")


def label_axis(ax, label: str, x: float = -0.18, y: float = 1.10) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=18, fontweight="bold", ha="left", va="top", clip_on=False)


def build_figure(output_dir: Path) -> None:
    fig = plt.figure(figsize=(15.4, 18.4), dpi=300)
    fig.subplots_adjust(left=0.055, right=0.985, top=0.965, bottom=0.055)

    # Three protein rows x two temperature columns, plus B spanning the bottom.
    outer = GridSpec(
        4,
        2,
        figure=fig,
        height_ratios=[1.20, 1.20, 1.20, 0.56],
        hspace=0.11,
        wspace=-0.04,
    )

    first_a_ax = None
    a_axes: dict[tuple[int, int], object] = {}
    for system_index, system in enumerate(SYSTEM_ORDER):
        for role_index, role in enumerate(ROLE_ORDER):
            ax = fig.add_subplot(outer[system_index, role_index])
            a_axes[(system_index, role_index)] = ax
            if first_a_ax is None:
                first_a_ax = ax
            plot_structural_transition_network(ax, system, role, transition_positions(system))
            if role_index == 0:
                ax.text(
                    0.020,
                    0.50,
                    display_system(system),
                    transform=ax.transAxes,
                    ha="right",
                    va="center",
                    fontsize=17.0,
                    fontweight="bold",
                )
            if system_index == 0:
                ax.set_title(temp_label(role), fontsize=19.0, fontweight="bold", pad=18)

    # B. Weighted transition distance, spanning the full bottom width.
    ax_td = fig.add_subplot(outer[3, :])
    x = np.arange(len(SYSTEM_ORDER), dtype=float)
    offsets = {"free": -0.12, "condensate": 0.12}
    sizes = {"free": 95, "condensate": 105}
    td_source_rows = []
    for family in FAMILY_ORDER:
        values = [weighted_transition_distance(system, family) for system in SYSTEM_ORDER]
        for system, value in zip(SYSTEM_ORDER, values):
            td_source_rows.append(
                {
                    "system": system,
                    "family": family,
                    "weighted_transition_distance": f"{value:.8f}",
                }
            )
        ax_td.scatter(
            x + offsets[family],
            values,
            s=sizes[family],
            facecolors="white" if family == "free" else FAMILY_COLORS[family],
            edgecolors=FAMILY_COLORS[family],
            linewidths=1.8,
            marker="o",
            label=FAMILY_LABELS[family],
            zorder=4,
        )
    ax_td.set_xticks(x, display_systems())
    ax_td.set_ylabel(r"Transition distance, $D_T^{(w)}$", labelpad=8)
    ax_td.legend(
        loc="upper left",
        ncol=2,
        fontsize=10.5,
        handletextpad=0.45,
        columnspacing=1.0,
        borderaxespad=0.2,
    )
    box_axis(ax_td)

    fig.canvas.draw()

    # Match panel B exactly to the outer left/right edges of the two A columns.
    # This makes B narrower than the full page and visually column-aligned.
    left_edge = min(a_axes[(row, 0)].get_position().x0 for row in range(3))
    right_edge = max(a_axes[(row, 1)].get_position().x1 for row in range(3))
    b_pos0 = ax_td.get_position()
    # Narrow B while keeping it centered under the two A columns.
    inset = 0.065 * (right_edge - left_edge)
    ax_td.set_position([
        left_edge + inset,
        b_pos0.y0,
        (right_edge - left_edge) - 2.0 * inset,
        b_pos0.height,
    ])
    fig.canvas.draw()

    panel_label_x = None
    if first_a_ax is not None:
        a_pos = first_a_ax.get_position()
        panel_label_x = max(0.006, a_pos.x0 - 0.080)
        fig.text(
            panel_label_x,
            min(0.985, a_pos.y1 + 0.020),
            "A",
            fontsize=22,
            fontweight="bold",
            ha="left",
            va="top",
        )
    b_pos = ax_td.get_position()
    fig.text(
        panel_label_x if panel_label_x is not None else max(0.006, b_pos.x0 - 0.080),
        min(0.985, b_pos.y1 + 0.020),
        "B",
        fontsize=22,
        fontweight="bold",
        ha="left",
        va="top",
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        output_dir / "Figure3_transition_distance_source.csv",
        td_source_rows,
        ["system", "family", "weighted_transition_distance"],
    )
    name = "Figure3_msm_state_dynamics_v12"
    for ext in ("png", "pdf"):
        fig.savefig(output_dir / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {output_dir / name}.png")
    print(f"[save] {output_dir / name}.pdf")

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    setup_style()
    build_figure(args.output_dir)


if __name__ == "__main__":
    main()

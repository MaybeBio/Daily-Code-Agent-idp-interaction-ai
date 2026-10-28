#!/usr/bin/env python3
"""Build Figure 3: Markov-state model and state-dynamics story."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
from supplementary.shared import mapped_matrix_stack, state_rows


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(os.environ.get("JPCL_DATA_ROOT", REPOSITORY_ROOT)).resolve()
DEFAULT_OUT = REPOSITORY_ROOT / "outputs/main"
RESULTS = PROJECT_ROOT / "ppt_5_11_figures/results"
PAPER = PROJECT_ROOT / "ppt_5_11_figures/paper_story_figures"

SYSTEM_ORDER = ("FUS", "A1", "TDP43")
DISPLAY_NAMES = {"FUS": "FUS-LCD", "A1": "A1-LCD", "TDP43": "TDP-43-LCD"}
FAMILY_ORDER = ("free", "condensate")
ROLE_ORDER = ("Tl", "Th")
TEMP_COLORS = {"Tl": "#2f7f73", "Th": "#b54b5f"}
FAMILY_COLORS = {"free": "#4c78a8", "condensate": "#e8874f"}
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
        "free": {
            "Tl": PROJECT_ROOT / "fus/fus_protein_single_simulation/tica_100dimension_10frames_10_stride_20_trajs_280K/avg_contacts.npy",
            "Th": PROJECT_ROOT / "fus/fus_protein_single_simulation/tica_100dimension_10frames_10_stride_20_trajs_325K/avg_contacts.npy",
        },
        "condensate": {
            "Tl": PROJECT_ROOT / "fus/fus_protein_single_simulation/condensate/tica_100dimension_10frames_280K/avg_contacts.npy",
            "Th": PROJECT_ROOT / "fus/fus_protein_single_simulation/condensate/tica_100dimension_10frames_325K/avg_contacts.npy",
        },
    },
    "A1": {
        "free": {
            "Tl": PROJECT_ROOT / "A1/A1_protein_single_simulation/tica_100dimension_10frames_230K/avg_contacts.npy",
            "Th": PROJECT_ROOT / "A1/A1_protein_single_simulation/tica_100dimension_10frames_260K/avg_contacts.npy",
        },
        "condensate": {
            "Tl": PROJECT_ROOT / "A1/A1_protein_single_simulation/condensate/tica_100dimension_10frames_230K/avg_contacts.npy",
            "Th": PROJECT_ROOT / "A1/A1_protein_single_simulation/condensate/tica_100dimension_10frames_260K/avg_contacts.npy",
        },
    },
    "TDP43": {
        "free": {
            "Tl": PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/tica_100dimension_10frames_260K/avg_contacts.npy",
            "Th": PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/tica_100dimension_10frames_300K/avg_contacts.npy",
        },
        "condensate": {
            "Tl": PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/condensate/tica_100dimension_10frames_260K/avg_contacts.npy",
            "Th": PROJECT_ROOT / "TDP43/tdp43_protein_single_simulation/condensate/tica_100dimension_10frames_300K/avg_contacts.npy",
        },
    },
}
CONTACT_FINGERPRINT_STATES = {
    "FUS": {"Tl": 0, "Th": 2},
    "A1": {"Tl": 0, "Th": 3},
    "TDP43": {"Tl": 0, "Th": 3},
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
            "font.size": 12.5,
            "axes.labelsize": 13.0,
            "xtick.labelsize": 11.5,
            "ytick.labelsize": 11.5,
            "legend.fontsize": 11.0,
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


def state_vector(
    system: str,
    family: str,
    role: str,
    *,
    normalize: bool = True,
) -> dict[int, float]:
    temp = TEMPS[system][role]
    vector: dict[int, float] = {}
    for row in read_csv(STATE_SOURCE[system][family]):
        if row_temperature(row) != temp:
            continue
        state = int(row.get("display_state") or row.get("global_state", "0"))
        vector[state] = vector.get(state, 0.0) + max(fval(row["mu"]), 0.0)
    if not normalize:
        return vector
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


def contact_map_for_state(system: str, role: str, display_state: int, family: str = "condensate") -> np.ndarray:
    local_state = local_state_for_display(system, family, role, display_state)
    contacts = np.load(AVG_CONTACTS[system][family][role], mmap_mode="r")
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
            3: (1.25, -3.65),
            4: (3.45, 0.70),
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


def ensemble_mean_transition_distance(
    system: str,
    family: str,
    n_bootstrap: int = 10000,
) -> tuple[float, float, float, int]:
    """Distance between ensemble-mean networks with independent model-bootstrap uncertainty."""
    states = transition_states(system)
    cache_dir = DEFAULT_OUT / "supplementary"
    tl_t, tl_p = mapped_matrix_stack(system, family, "Tl", states, cache_dir)
    th_t, th_p = mapped_matrix_stack(system, family, "Th", states, cache_dir)
    n_models = min(len(tl_t), len(th_t))
    tl_t, tl_p = tl_t[:n_models], tl_p[:n_models]
    th_t, th_p = th_t[:n_models], th_p[:n_models]
    index = {state: idx for idx, state in enumerate(states)}
    masks = {}
    for role in ROLE_ORDER:
        mask = np.zeros((len(states), len(states)), dtype=float)
        for row in transition_network_rows(system):
            if row["family"] != family or row["role"] != role:
                continue
            source, target = row_display_states(row)
            if source in index and target in index:
                mask[index[source], index[target]] = 1.0
        masks[role] = mask

    def distance(indices_tl: np.ndarray, indices_th: np.ndarray) -> float:
        mean_tl = np.mean(tl_t[indices_tl], axis=0) * masks["Tl"]
        mean_th = np.mean(th_t[indices_th], axis=0) * masks["Th"]
        mean_pi = 0.5 * (
            np.mean(tl_p[indices_tl], axis=0)
            + np.mean(th_p[indices_th], axis=0)
        )
        return float(
            np.sqrt(np.sum(mean_pi[:, None] * (mean_th - mean_tl) ** 2))
        )

    all_models = np.arange(n_models)
    estimate = distance(all_models, all_models)
    seed = {
        ("FUS", "free"): 3101,
        ("FUS", "condensate"): 3102,
        ("A1", "free"): 3201,
        ("A1", "condensate"): 3202,
        ("TDP43", "free"): 3301,
        ("TDP43", "condensate"): 3302,
    }[(system, family)]
    rng = np.random.default_rng(seed)
    bootstrap = np.empty(n_bootstrap, dtype=float)
    for idx in range(n_bootstrap):
        bootstrap[idx] = distance(
            rng.integers(0, n_models, n_models),
            rng.integers(0, n_models, n_models),
        )
    low, high = np.percentile(bootstrap, [2.5, 97.5])
    return estimate, float(low), float(high), n_models


def transition_distance_contrast(
    system: str,
    n_bootstrap: int = 20000,
) -> tuple[float, float, float, float]:
    """Bootstrap the condensate-minus-free difference in ensemble-mean distance."""
    states = transition_states(system)
    cache_dir = DEFAULT_OUT / "supplementary"
    family_data = {}
    index = {state: idx for idx, state in enumerate(states)}
    for family in FAMILY_ORDER:
        role_data = {}
        for role in ROLE_ORDER:
            matrices, populations = mapped_matrix_stack(
                system, family, role, states, cache_dir
            )
            mask = np.zeros((len(states), len(states)), dtype=float)
            for row in transition_network_rows(system):
                if row["family"] != family or row["role"] != role:
                    continue
                source, target = row_display_states(row)
                if source in index and target in index:
                    mask[index[source], index[target]] = 1.0
            role_data[role] = (matrices, populations, mask)
        family_data[family] = role_data

    def family_distance(family: str, indices_tl: np.ndarray, indices_th: np.ndarray) -> float:
        tl_t, tl_p, tl_mask = family_data[family]["Tl"]
        th_t, th_p, th_mask = family_data[family]["Th"]
        mean_tl = np.mean(tl_t[indices_tl], axis=0) * tl_mask
        mean_th = np.mean(th_t[indices_th], axis=0) * th_mask
        mean_pi = 0.5 * (
            np.mean(tl_p[indices_tl], axis=0)
            + np.mean(th_p[indices_th], axis=0)
        )
        return float(np.sqrt(np.sum(mean_pi[:, None] * (mean_th - mean_tl) ** 2)))

    all_indices = {
        family: np.arange(
            min(
                len(family_data[family]["Tl"][0]),
                len(family_data[family]["Th"][0]),
            )
        )
        for family in FAMILY_ORDER
    }
    point = (
        family_distance("condensate", all_indices["condensate"], all_indices["condensate"])
        - family_distance("free", all_indices["free"], all_indices["free"])
    )
    rng = np.random.default_rng({"FUS": 4101, "A1": 4201, "TDP43": 4301}[system])
    bootstrap = np.empty(n_bootstrap, dtype=float)
    for idx in range(n_bootstrap):
        sampled = {}
        for family in FAMILY_ORDER:
            n_models = len(all_indices[family])
            sampled[(family, "Tl")] = rng.integers(0, n_models, n_models)
            sampled[(family, "Th")] = rng.integers(0, n_models, n_models)
        bootstrap[idx] = (
            family_distance(
                "condensate",
                sampled[("condensate", "Tl")],
                sampled[("condensate", "Th")],
            )
            - family_distance(
                "free",
                sampled[("free", "Tl")],
                sampled[("free", "Th")],
            )
        )
    low, high = np.percentile(bootstrap, [2.5, 97.5])
    probability_positive = float(np.mean(bootstrap > 0.0))
    return point, float(low), float(high), probability_positive


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

    # Pull the two temperature columns toward the center without changing the
    # internal geometry of either motif.
    x_shift = 0.075 if role == "Tl" else -0.075
    if system == "TDP43" and role == "Th":
        x_shift = -0.075
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

    def square_edge_point(c: np.ndarray, direction: np.ndarray, gap: float = 0.014) -> np.ndarray:
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

    grouped: dict[tuple[int, int], list[dict[str, str]]] = {}
    for row in best_by_direction.values():
        s, t = row_display_states(row)
        grouped.setdefault(tuple(sorted((s, t))), []).append(row)

    for pair, pair_rows in grouped.items():
        a, b = pair
        ca, cb = center(a), center(b)
        midpoint = 0.5 * (ca + cb)
        chord = cb - ca
        chord_unit = normalize(chord)
        normal = np.array([-chord_unit[1], chord_unit[0]], dtype=float)

        # Choose the side that points away from the motif center. For the lower
        # horizontal edge, explicitly use the region below the contact maps.
        outward = midpoint - panel_center
        if abs(chord[0]) > 1.7 * abs(chord[1]):
            outward = np.array([0.0, -1.0 if midpoint[1] < panel_center[1] else 1.0])
        elif np.hypot(outward[0], outward[1]) < 1e-10:
            outward = normal
        side_sign = 1.0 if float(np.dot(normal, outward)) >= 0.0 else -1.0
        outward_normal = side_sign * normal

        # Bidirectional pairs are shown as separated straight arrows. Certain
        # special one-way edges (for example TDP-43 TH S1→S6) are also forced
        # to straight lines for clarity.
        pair_rows = sorted(pair_rows, key=lambda r: row_display_states(r))
        if len(states) == 3:
            parallel_offsets = [-0.055, 0.055] if len(pair_rows) > 1 else [0.0]
            for parallel_offset, row in zip(parallel_offsets, pair_rows):
                source, target = row_display_states(row)
                probability = fval(row["transition_probability"])
                cs = center(source) + parallel_offset * normal
                ct = center(target) + parallel_offset * normal
                direct = normalize(ct - cs)
                if np.hypot(direct[0], direct[1]) < 1e-12:
                    continue
                start_pt = square_edge_point(cs, direct, gap=0.034)
                end_pt = square_edge_point(ct, -direct, gap=0.034)
                ax.add_patch(
                    FancyArrowPatch(
                        start_pt,
                        end_pt,
                        arrowstyle="-|>",
                        mutation_scale=18,
                        linewidth=2.2,
                        color="#505050",
                        alpha=0.95,
                        transform=ax.transAxes,
                        clip_on=False,
                        zorder=6,
                    )
                )
                if len(pair_rows) > 1:
                    offset_sign = 1.0 if parallel_offset > 0 else -1.0
                    label_xy = 0.5 * (start_pt + end_pt) + 0.040 * offset_sign * normal
                else:
                    label_xy = 0.5 * (start_pt + end_pt) + 0.035 * outward_normal
                ax.text(
                    label_xy[0],
                    label_xy[1],
                    f"{probability:.2f}",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=10.8,
                    color="#303030",
                    zorder=7,
                    clip_on=False,
                )
            continue

        tdp43_th_straight_pairs = {
            frozenset((2, 4)),  # S5 -> S3
            frozenset((0, 3)),  # S1 -> S4
            frozenset((2, 3)),  # S4 -> S3
            frozenset((0, 2)),  # S1 -> S3
        }
        force_tdp43_th_straight = (
            system == "TDP43"
            and role == "Th"
            and frozenset(pair) in tdp43_th_straight_pairs
        )
        use_straight_pair = len(pair_rows) > 1 or force_tdp43_th_straight
        if use_straight_pair:
            if len(pair_rows) == 1:
                parallel_offsets = [0.0]
                along_t_values = [0.52]
            else:
                parallel_offsets = [-0.085, 0.085]
                along_t_values = [0.40, 0.60]

            for lane_index, (parallel_offset, along_t, row) in enumerate(zip(parallel_offsets, along_t_values, pair_rows)):
                source, target = row_display_states(row)
                probability = fval(row["transition_probability"])
                cs, ct = center(source), center(target)
                direct = normalize(ct - cs)
                if np.hypot(direct[0], direct[1]) < 1e-12:
                    continue

                line_normal = np.array([-direct[1], direct[0]], dtype=float)
                offset_vec = parallel_offset * line_normal
                start_pt = square_edge_point(cs, direct, gap=0.014) + offset_vec
                end_pt = square_edge_point(ct, -direct, gap=0.018) + offset_vec
                if system == "TDP43" and role == "Th":
                    half = 0.5 * thumb_size
                    if (source, target) == (0, 2):
                        # S1 -> S3: leave the midpoint of S1's left edge and
                        # enter the midpoint of S3's top edge.
                        start_pt = np.array([cs[0] - half - 0.018, cs[1]])
                        end_pt = np.array([ct[0], ct[1] + half + 0.018])
                    elif (source, target) == (3, 2):
                        # S4 -> S3: leave the midpoint of S4's left edge and
                        # enter the lower-right part of S3. The angled entry
                        # keeps the arrowhead clear of the centered state label.
                        start_pt = np.array([cs[0] - half - 0.018, cs[1]])
                        end_pt = np.array([ct[0] + 0.090, ct[1] - half - 0.018])
                    elif (source, target) == (0, 3):
                        # Start below the S1 population label so the vertical
                        # S1 -> S4 arrow does not run through the text.
                        start_pt = np.array([cs[0], cs[1] - half - 0.075])
                        end_pt = np.array([ct[0], ct[1] + half + 0.030])
                    direct = normalize(end_pt - start_pt)
                    line_normal = np.array([-direct[1], direct[0]], dtype=float)
                arrow = FancyArrowPatch(
                    start_pt,
                    end_pt,
                    arrowstyle="-|>",
                    mutation_scale=23,
                    linewidth=1.65 + 5.0 * probability,
                    color="#505050",
                    alpha=0.95,
                    transform=ax.transAxes,
                    clip_on=False,
                    zorder=6,
                )
                ax.add_patch(arrow)

                if len(pair_rows) == 1:
                    label_xy = 0.5 * (start_pt + end_pt) + 0.008 * line_normal
                else:
                    offset_sign = 1.0 if parallel_offset >= 0 else -1.0
                    label_xy = start_pt + along_t * (end_pt - start_pt)
                    label_xy = label_xy + 0.020 * offset_sign * line_normal
                if system == "TDP43" and role == "Th" and len(pair_rows) == 1:
                    # Keep each value close to its own arrow, but away from the
                    # crowded network center.
                    label_t = {
                        (4, 2): 0.18,  # S5 -> S3
                        (3, 2): 0.52,  # S4 -> S3
                        (0, 2): 0.42,  # S1 -> S3
                        (0, 3): 0.05,  # S1 -> S4
                    }.get((source, target), 0.48)
                    label_side = {
                        (4, 2): -1.0,
                        (3, 2): 1.0,
                        (0, 2): -1.0,
                        (0, 3): -1.0,
                    }.get((source, target), 1.0)
                    label_offset = 0.022
                    if (source, target) == (0, 3):
                        label_offset = 0.045
                    elif (source, target) == (4, 2):
                        label_offset = 0.028
                    elif (source, target) == (0, 2):
                        label_offset = 0.038
                    elif (source, target) == (3, 2):
                        label_offset = 0.050
                    label_xy = (
                        start_pt
                        + label_t * (end_pt - start_pt)
                        + label_offset * label_side * line_normal
                    )
                ax.text(
                    label_xy[0], label_xy[1], f"{probability:.2f}",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=10.8, color="#303030",
                    zorder=7, clip_on=False,
                )
            continue

        # Single-direction pairs keep a clean exterior lane.
        lane_offsets = [0.21]

        for lane_offset, row in zip(lane_offsets, pair_rows):
            source, target = row_display_states(row)
            probability = fval(row["transition_probability"])
            cs, ct = center(source), center(target)

            # Aim the entry/exit directions slightly outward so the lines leave
            # the outer corners rather than cutting through the gap between maps.
            source_dir = normalize((ct - cs) + 0.65 * outward_normal)
            target_dir = normalize((cs - ct) + 0.65 * outward_normal)
            start_pt = square_edge_point(cs, source_dir, gap=0.012)
            end_pt = square_edge_point(ct, target_dir, gap=0.016)
            control = 0.5 * (start_pt + end_pt) + lane_offset * outward_normal

            path = matplotlib.path.Path(
                [start_pt, control, end_pt],
                [matplotlib.path.Path.MOVETO, matplotlib.path.Path.CURVE3, matplotlib.path.Path.CURVE3],
            )
            arrow = FancyArrowPatch(
                path=path,
                arrowstyle="-|>",
                mutation_scale=23,
                linewidth=1.65 + 5.0 * probability,
                color="#505050",
                alpha=0.95,
                transform=ax.transAxes,
                clip_on=False,
                zorder=6,
            )
            ax.add_patch(arrow)

            # Label near the apex of its own Bézier lane.
            label_xy = 0.25 * start_pt + 0.50 * control + 0.25 * end_pt
            label_xy = label_xy + 0.007 * outward_normal
            ax.text(
                label_xy[0], label_xy[1], f"{probability:.2f}",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=10.8, color="#303030",
                zorder=7, clip_on=False,
            )

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
            fontsize=13.0, fontweight="bold", color="#222222",
            zorder=8, clip_on=False,
        )

def plot_combined_structural_transition_network(
    ax,
    system: str,
    family: str,
    show_triangle_temperature_labels: bool = False,
) -> None:
    """Overlay low- and high-temperature dynamics on one structural network."""
    # Display the original ensemble populations.  Do not renormalize after
    # omitting states that were not retained by cross-condition alignment;
    # this keeps the labels identical to the full state-resolved maps in SI.
    populations = {
        role: state_vector(system, family, role, normalize=False)
        for role in ROLE_ORDER
    }
    states = sorted(set(populations["Tl"]) | set(populations["Th"]))
    ordered = ordered_states_for_layout(states)
    if len(states) == 3:
        coords, thumb_size = [(0.50, 0.84), (0.16, 0.20), (0.84, 0.20)], 0.225
    elif len(states) == 4:
        coords, thumb_size = [(0.50, 0.87), (0.13, 0.56), (0.50, 0.13), (0.87, 0.56)], 0.198
    else:
        coords = [(0.50, 0.88), (0.12, 0.64), (0.25, 0.14), (0.75, 0.14), (0.88, 0.64)]
        thumb_size = 0.176
    if system == "A1" and family == "condensate":
        # A square layout avoids the former central crossing: s2<->s1 is
        # vertical at left, s2<->s4 is horizontal below, s3->s1 is horizontal
        # above, and the remaining s4->s1 path occupies the diagonal.
        coords = [(0.25, 0.80), (0.25, 0.20), (0.75, 0.80), (0.75, 0.20)]
        thumb_size = 0.187
    positions = {state: np.asarray(coords[i], dtype=float) for i, state in enumerate(ordered)}

    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal"); ax.axis("off")

    def normalize(vector: np.ndarray) -> np.ndarray:
        length = float(np.hypot(vector[0], vector[1]))
        return vector / length if length > 1e-12 else np.zeros(2)

    def edge_point(center: np.ndarray, direction: np.ndarray) -> np.ndarray:
        unit = normalize(direction)
        distance = (0.5 * thumb_size) / max(abs(unit[0]), abs(unit[1]), 1e-12)
        # Connections that leave a node through its lower edge need extra
        # clearance for the state label and populations. Side and upper-edge
        # connections can be longer without entering that annotation zone.
        annotation_gap = 0.080 if unit[1] < -0.25 else 0.043
        return center + (distance + annotation_gap) * unit

    visible: dict[tuple[str, int, int], float] = {}
    for role in ROLE_ORDER:
        for row in transition_network_rows(system):
            if row["family"] != family or row["role"] != role:
                continue
            source, target = row_display_states(row)
            probability = fval(row["transition_probability"])
            if source == target or probability <= 0.10:
                continue
            key = (role, source, target)
            visible[key] = max(probability, visible.get(key, -np.inf))

    pair_groups: dict[tuple[int, int], list[tuple[str, int, int, float]]] = {}
    for (role, source, target), probability in visible.items():
        if source in positions and target in positions:
            pair_groups.setdefault(tuple(sorted((source, target))), []).append(
                (role, source, target, probability)
            )

    routed_edges: list[tuple[str, int, int, float, float, float]] = []
    lane_spacing = 0.065
    for pair, edges in sorted(pair_groups.items()):
        low_state, _high_state = pair
        edges.sort(key=lambda item: (
            ROLE_ORDER.index(item[0]),
            0 if item[1] == low_state else 1,
        ))
        center = 0.5 * (len(edges) - 1)
        for lane_index, edge in enumerate(edges):
            lane_rank = lane_index - center
            routed_edges.append((*edge, lane_rank * lane_spacing, lane_rank * 0.180))

    for role, source, target, probability, lateral, label_shift in routed_edges:
        color = TEMP_COLORS[role]
        direct = normalize(positions[target] - positions[source])
        # Only visible arrows receive lanes. Their fixed spacing is centered
        # on the state-pair axis, so omitted sub-threshold transitions leave
        # neither a blank lane nor a one-sided bundle.
        low_state, high_state = sorted((source, target))
        canonical = normalize(positions[high_state] - positions[low_state])
        canonical_normal = np.array([-canonical[1], canonical[0]])
        start = edge_point(positions[source], direct) + lateral * canonical_normal
        end = edge_point(positions[target], -direct) + lateral * canonical_normal

        # Use a common geometric construction for every lane in a state pair,
        # so reverse and temperature-resolved arrows have matching lengths.
        if system == "TDP43" and family == "condensate" and source == 1 and target == 0:
            # Extend the sixth panel's s2->s1 arrow toward both map edges while
            # retaining a small clear gap from each contact map.
            start = start - 0.040 * direct
            end = end + 0.040 * direct
        arrow_width = max(1.60, 0.55 * (probability / 0.10) ** 2)
        head_width = max(5.5, 1.55 * arrow_width)
        head_length = max(6.5, 1.90 * arrow_width)
        arrow_style = (
            f"Simple,tail_width={arrow_width:.3f},"
            f"head_width={head_width:.3f},head_length={head_length:.3f}"
        )
        ax.add_patch(FancyArrowPatch(
            start, end, arrowstyle=arrow_style, mutation_scale=1.0,
            linewidth=0.0, facecolor=color, edgecolor=color, alpha=0.95,
            capstyle="butt", joinstyle="miter",
            connectionstyle="arc3,rad=0", transform=ax.transAxes,
            clip_on=False, zorder=3,
        ))
        # Stagger values substantially along their corresponding arrows.
        # They remain directly on the line, but stay away from both the
        # common midpoint and the arrowhead/node regions.
        label_fraction = 0.5 + label_shift
        if system == "TDP43" and family == "condensate" and {
            source, target
        } in ({1, 2}, {0, 2}):
            # Using the same directed fraction places reverse-arrow labels on
            # opposite physical sides of the pair, separating 0.14/0.18 and
            # 0.19/0.18, respectively.
            label_fraction = 0.34
        on_line = start + label_fraction * (end - start)
        label = on_line
        ax.text(label[0], label[1], f"{probability:.2f}", transform=ax.transAxes,
                ha="center", va="center", fontsize=10.4, fontweight="bold", color=color,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.94, "pad": 0.45},
                zorder=4, clip_on=False)

    cmap = matplotlib.colormaps["viridis"].copy(); cmap.set_bad("white")
    for state in states:
        center = positions[state]
        x0, x1 = center[0] - thumb_size / 2, center[0] + thumb_size / 2
        y0, y1 = center[1] - thumb_size / 2, center[1] + thumb_size / 2
        maps = {role: contact_map_for_state(system, role, state, family)
                for role in ROLE_ORDER if state in populations[role]}
        template = next(iter(maps.values()))
        low = maps.get("Tl", np.full_like(template, np.nan, dtype=float))
        high = maps.get("Th", np.full_like(template, np.nan, dtype=float))
        n = min(low.shape[0], high.shape[0])
        low, high = np.asarray(low[:n, :n], float), np.asarray(high[:n, :n], float)
        row_index, column_index = np.indices((n, n))
        combined = np.where(row_index >= column_index, low, high)
        np.fill_diagonal(combined, np.nan)
        ax.imshow(combined, origin="lower", cmap=cmap, norm=LogNorm(vmin=1e-4, vmax=1e-1),
                  interpolation="nearest", extent=(x0, x1, y0, y1),
                  transform=ax.transAxes, clip_on=False, zorder=5)
        ax.add_patch(Rectangle((x0, y0), thumb_size, thumb_size, transform=ax.transAxes,
                               facecolor="none", edgecolor="#2f2f2f", linewidth=1.1,
                               clip_on=False, zorder=6))
        if show_triangle_temperature_labels and state == 0:
            ax.text(x0 + 0.012, y1 - 0.012, r"$T_{\mathrm{L}}$", transform=ax.transAxes,
                    ha="left", va="top", fontsize=11.5, fontweight="bold",
                    color=TEMP_COLORS["Tl"], zorder=7,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 0.35})
            ax.text(x1 - 0.012, y0 + 0.012, r"$T_{\mathrm{H}}$", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=11.5, fontweight="bold",
                    color=TEMP_COLORS["Th"], zorder=7,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.88, "pad": 0.35})
        ax.text(center[0], y0 - 0.018, rf"$s_{{{state + 1}}}$", transform=ax.transAxes,
                ha="center", va="top", fontsize=15.2, fontweight="bold", color="#202020", zorder=7,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.94, "pad": 0.4})
        low_text = rf"{100 * populations['Tl'].get(state, 0):.0f}%" if state in populations["Tl"] else "--"
        high_text = rf"{100 * populations['Th'].get(state, 0):.0f}%" if state in populations["Th"] else "--"
        population_above = (
            state == 0
            or (system == "A1" and family == "condensate" and state == 2)
            or (system == "TDP43" and family == "condensate" and state in (1, 4))
        )
        population_y = y1 + 0.012 if population_above else y0 - 0.061
        population_va = "bottom" if population_above else "top"
        ax.text(center[0] - 0.020, population_y, low_text, transform=ax.transAxes,
                ha="right", va=population_va, fontsize=10.5, fontweight="bold", color=TEMP_COLORS["Tl"], zorder=7,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.94, "pad": 0.3})
        ax.text(center[0] + 0.020, population_y, high_text, transform=ax.transAxes,
                ha="left", va=population_va, fontsize=10.5, fontweight="bold", color=TEMP_COLORS["Th"], zorder=7,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.94, "pad": 0.3})

    ax.plot([], [], color=TEMP_COLORS["Tl"], lw=3.0, label=r"$T_{\mathrm{L}}$")
    ax.plot([], [], color=TEMP_COLORS["Th"], lw=3.0, label=r"$T_{\mathrm{H}}$")


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
    fig = plt.figure(figsize=(15.6, 12.4), dpi=300)
    fig.subplots_adjust(left=0.070, right=0.985, top=0.940, bottom=0.060)

    # Free and condensate networks are paired by system, with B below.
    outer = GridSpec(
        3,
        3,
        figure=fig,
        height_ratios=[1.0, 1.0, 0.43],
        hspace=0.10,
        wspace=0.12,
        right=0.98,
    )

    first_a_ax = None
    a_axes: dict[tuple[int, int], object] = {}
    for family_index, family in enumerate(FAMILY_ORDER):
        for system_index, system in enumerate(SYSTEM_ORDER):
            ax = fig.add_subplot(outer[family_index, system_index])
            a_axes[(family_index, system_index)] = ax
            if first_a_ax is None:
                first_a_ax = ax
            plot_combined_structural_transition_network(
                ax,
                system,
                family,
                show_triangle_temperature_labels=(family == "free" and system == "FUS"),
            )
            if family_index == 0:
                ax.set_title(display_system(system), fontsize=19.0,
                             fontweight="bold", pad=13)
            if system_index == 0:
                ax.text(-0.055, 0.50, FAMILY_LABELS[family], transform=ax.transAxes,
                        rotation=90, ha="center", va="center", fontsize=17.0,
                        fontweight="bold", clip_on=False)
    fig.canvas.draw()
    first_box = a_axes[(0, 0)].get_position()
    second_box = a_axes[(0, 1)].get_position()
    # Align the legend center with the upper state's y coordinate and the
    # lower-right state's x coordinate in the first FUS network.
    legend_x = first_box.x0 + 0.84 * (first_box.x1 - first_box.x0)
    legend_y = first_box.y0 + 0.84 * (first_box.y1 - first_box.y0)
    fig.legend(
        handles=[
            Line2D([0], [0], color=TEMP_COLORS["Tl"], lw=3.0, label=r"$T_{\mathrm{L}}$"),
            Line2D([0], [0], color=TEMP_COLORS["Th"], lw=3.0, label=r"$T_{\mathrm{H}}$"),
        ],
        loc="center", bbox_to_anchor=(legend_x, legend_y), ncol=1,
        bbox_transform=fig.transFigure,
        fontsize=15.0, handlelength=1.7, handletextpad=0.45,
        borderaxespad=0.0,
    )
    # B. Weighted transition distance, spanning the full bottom width.
    ax_td = fig.add_subplot(outer[2, :])
    x = np.arange(len(SYSTEM_ORDER), dtype=float)
    offsets = {"free": -0.12, "condensate": 0.12}
    sizes = {"free": 95, "condensate": 105}
    td_source_rows = []
    td_contrast_rows = []
    td_plot_bounds = []
    for family in FAMILY_ORDER:
        values = []
        lows = []
        highs = []
        for system in SYSTEM_ORDER:
            value, low, high, n_models = ensemble_mean_transition_distance(
                system,
                family,
            )
            values.append(value)
            lows.append(low)
            highs.append(high)
            td_plot_bounds.extend([low, high])
            td_source_rows.append(
                {
                    "system": system,
                    "family": family,
                    "weighted_transition_distance": f"{value:.8f}",
                    "percentile_2.5": f"{low:.8f}",
                    "percentile_97.5": f"{high:.8f}",
                    "n_models": str(n_models),
                }
            )
        ax_td.errorbar(
            x + offsets[family],
            values,
            yerr=[
                np.asarray(values) - np.asarray(lows),
                np.asarray(highs) - np.asarray(values),
            ],
            fmt="o",
            markersize=np.sqrt(sizes[family]),
            markerfacecolor="white" if family == "free" else "#4a4a4a",
            markeredgecolor="#4a4a4a",
            markeredgewidth=1.8,
            color="#4a4a4a",
            linewidth=1.5,
            capsize=4,
            label=FAMILY_LABELS[family],
            zorder=4,
        )
    for system in SYSTEM_ORDER:
        difference, low, high, probability_positive = transition_distance_contrast(system)
        td_contrast_rows.append(
            {
                "system": system,
                "condensate_minus_free": f"{difference:.8f}",
                "percentile_2.5": f"{low:.8f}",
                "percentile_97.5": f"{high:.8f}",
                "bootstrap_probability_positive": f"{probability_positive:.6f}",
                "n_bootstrap": "20000",
                "resampling": "independent across temperature and environment",
            }
        )
    ax_td.set_xticks(x, display_systems())
    ax_td.set_ylabel(
        r"$D_T^{(w)}$",
        fontsize=17.0,
        labelpad=10,
    )
    ax_td.tick_params(axis="both", labelsize=16.0)
    td_min = min(td_plot_bounds)
    td_max = max(td_plot_bounds)
    td_span = max(td_max - td_min, 0.01)
    ax_td.set_ylim(td_min - 0.10 * td_span, td_max + 0.30 * td_span)
    ax_td.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="None", markersize=np.sqrt(sizes["free"]),
                   markerfacecolor="white", markeredgecolor="#4a4a4a", markeredgewidth=1.8,
                   color="#4a4a4a", label=FAMILY_LABELS["free"]),
            Line2D([0], [0], marker="o", linestyle="None", markersize=np.sqrt(sizes["condensate"]),
                   markerfacecolor="#4a4a4a", markeredgecolor="#4a4a4a", markeredgewidth=1.8,
                   color="#4a4a4a", label=FAMILY_LABELS["condensate"]),
        ],
        loc="upper left",
        ncol=2,
        fontsize=15.5,
        handletextpad=0.45,
        columnspacing=1.0,
        borderaxespad=0.2,
    )
    box_axis(ax_td)

    fig.canvas.draw()

    # Match panel B to the merged network width.
    left_edge = min(ax.get_position().x0 for ax in a_axes.values())
    right_edge = max(ax.get_position().x1 for ax in a_axes.values())
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

    if first_a_ax is not None:
        a_pos = first_a_ax.get_position()
        fig.text(
            0.055,
            min(0.985, a_pos.y1 + 0.020),
            "A",
            fontsize=24,
            fontweight="bold",
            ha="left",
            va="top",
        )
    b_pos = ax_td.get_position()
    fig.text(
        0.055,
        min(0.985, b_pos.y1 + 0.045),
        "B",
        fontsize=24,
        fontweight="bold",
        ha="left",
        va="top",
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        output_dir / "Figure3_transition_distance_source.csv",
        td_source_rows,
        [
            "system",
            "family",
            "weighted_transition_distance",
            "percentile_2.5",
            "percentile_97.5",
            "n_models",
        ],
    )
    write_csv(
        output_dir / "Figure3_transition_distance_contrast_source.csv",
        td_contrast_rows,
        [
            "system",
            "condensate_minus_free",
            "percentile_2.5",
            "percentile_97.5",
            "bootstrap_probability_positive",
            "n_bootstrap",
            "resampling",
        ],
    )
    name = "Figure3_msm_state_dynamics_v19"
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

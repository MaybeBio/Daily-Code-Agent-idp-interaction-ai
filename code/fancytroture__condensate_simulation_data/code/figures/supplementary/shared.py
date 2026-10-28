#!/usr/bin/env python3
"""Build the manuscript supporting-figure set for Figures 1-3.

The script only uses analysis products already generated from the selected
temperature pairs. It does not invent missing cutoff scans. Every panel is
written as PDF and 400 dpi PNG, with compact source tables where appropriate.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ.get("JPCL_DATA_ROOT", REPOSITORY_ROOT)).resolve()
DEFAULT_OUT = REPOSITORY_ROOT / "outputs/supplementary"
RESULTS = ROOT / "ppt_5_11_figures/results"
STORY = ROOT / "storyline_figures_for_ppt_corrected_v3"
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
PHASE_COLORS = {"free": "#555555", "condensate": "#c96d28"}
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

RG_DIR = {
    "FUS": ROOT / "fus/fus_protein_single_simulation/paper_thermo_figures_2temps_280_325_vampnet/figure1_rg_ree_relaxation",
    "A1": ROOT / "A1/A1_protein_single_simulation/paper_thermo_figures_a1_4p17_relaxed_mu10/figure1_rg_ree_relaxation",
    "TDP43": ROOT / "TDP43/tdp43_protein_single_simulation/paper_thermo_figures_tdp43_4p17_relaxed_mu10/figure1_rg_ree_relaxation",
}

DENSITY_DIR = {
    "FUS": ROOT / "fus/fus_condensate_simulation",
    "A1": ROOT / "A1/A1_condensate_simulation",
    "TDP43": ROOT / "TDP43/tdp43_condensate_simulation",
}
DENSITY_PREFIX = {"FUS": "condensate_", "A1": "condensate_", "TDP43": "condensate_"}

CONTACT_MAP_IMAGES = {
    "FUS": {
        "free": STORY / "state_population_contactmaps_full/fus_free_state_aligned_contactmaps.png",
        "condensate": STORY / "state_population_contactmaps_full/fus_condensate_aligned_contactmaps.png",
    },
    "A1": {
        "free": STORY / "state_population_contactmaps_full/a1_free_state_aligned_contactmaps.png",
        "condensate": STORY / "state_population_contactmaps_full/a1_condensate_aligned_contactmaps.png",
    },
    "TDP43": {
        "free": STORY / "state_population_contactmaps_full/tdp43_free_state_aligned_contactmaps.png",
        "condensate": STORY / "state_population_contactmaps_full/tdp43_condensate_aligned_contactmaps.png",
    },
}

TICA_IMAGES = {
    "FUS": RESULTS / "tica_state_projection_fus_free_condensate.png",
    "A1": RESULTS / "tica_state_projection_a1_free_condensate.png",
    "TDP43": RESULTS / "tica_state_projection_tdp43_free_condensate.png",
}

STATE_SOURCE = {
    "FUS": {
        "free": RG_DIR["FUS"].parent / "figure3_state_population_contactmaps/figure3_free_state_population_contactmaps_source.csv",
        "condensate": RG_DIR["FUS"].parent / "figure3_state_population_contactmaps/figure3_condensate_state_population_contactmaps_source.csv",
    },
    "A1": {
        "free": RESULTS / "a1_230_260_state_mapping/figure3_free_state_population_contactmaps_source.csv",
        "condensate": RESULTS / "a1_230_260_state_mapping/figure3_condensate_state_population_contactmaps_source.csv",
    },
    "TDP43": {
        "free": RG_DIR["TDP43"].parent / "figure3_state_population_contactmaps/figure3_free_state_population_contactmaps_source.csv",
        "condensate": RG_DIR["TDP43"].parent / "figure3_state_population_contactmaps/figure3_condensate_state_population_contactmaps_source.csv",
    },
}



def model_dir(system: str, family: str, role: str) -> Path:
    temp = TEMPS[system][role]
    if system == "FUS":
        base = ROOT / "fus/fus_protein_single_simulation"
        name = f"tica_100dimension_10frames_10_stride_20_trajs_{temp}" if family == "free" else f"tica_100dimension_10frames_{temp}"
        return base / name if family == "free" else base / "condensate" / name
    if system == "A1":
        base = ROOT / "A1/A1_protein_single_simulation"
        name = f"tica_100dimension_10frames_{temp}"
        return base / name if family == "free" else base / "condensate" / name
    base = ROOT / "TDP43/tdp43_protein_single_simulation"
    name = f"tica_100dimension_10frames_{temp}"
    return base / name if family == "free" else base / "condensate" / name


def intermediate_name(system: str) -> str:
    return {"FUS": "intermediate_fus", "A1": "intermediate_a1", "TDP43": "intermediate_tdp43"}[system]


def setup_style() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10.5,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.fontsize": 9,
            "legend.frameon": False,
            "axes.linewidth": 0.85,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    for row in rows[1:]:
        for field in row:
            if field not in fieldnames:
                fieldnames.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def save(fig, out: Path, name: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    # 400 dpi keeps the large contact-map/TICA composites publication-ready
    # without exhausting memory during the raster export. PDF is also saved.
    fig.savefig(out / f"{name}.png", dpi=400, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {name}")


def panel(ax, label: str) -> None:
    ax.text(-0.12, 1.05, label, transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom", clip_on=False)


def clean(ax, grid: bool = True) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    if grid:
        ax.grid(axis="y", color="#bdbdbd", alpha=0.25, lw=0.6)


def temp_label(role: str) -> str:
    return r"$T_{\mathrm{L}}$" if role == "Tl" else r"$T_{\mathrm{H}}$"


def rg_cache(system: str, family: str, role: str) -> dict[str, np.ndarray]:
    path = RG_DIR[system] / f"cache_{family}_{TEMPS[system][role]}.npz"
    data = np.load(path, allow_pickle=False)
    return {key: data[key] for key in data.files}


def intra_cache(system: str, family: str, role: str) -> dict[str, np.ndarray]:
    path = INTRA_DIR / f"intrachain_{system}_{family}_{TEMPS[system][role]}.npz"
    data = np.load(path, allow_pickle=False)
    return {key: data[key] for key in data.files}


def crossing(time: np.ndarray, values: np.ndarray, threshold: float) -> float:
    time = np.asarray(time, float)
    values = np.asarray(values, float)
    for idx in range(1, len(values)):
        if values[idx - 1] >= threshold >= values[idx]:
            y0, y1 = values[idx - 1], values[idx]
            if np.isclose(y0, y1):
                return float(time[idx])
            return float(time[idx - 1] + (threshold - y0) * (time[idx] - time[idx - 1]) / (y1 - y0))
    return float("nan")


def density_paths(system: str, role: str) -> tuple[Path, Path]:
    folder = DENSITY_DIR[system] / f"{DENSITY_PREFIX[system]}{TEMPS[system][role]}_200chains"
    return folder / "overall_rhoz_gml.dat", folder / "block_rhoz_mgml.dat"


def slab_snapshot(system: str, role: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    import mdtraj as md

    folder = density_paths(system, role)[0].parent
    prefix = {"FUS": "fus", "A1": "a1", "TDP43": "tdp43"}[system]
    trajectories = sorted(folder.glob(f"{prefix}_slab_vampnet_part*.dcd"), key=numeric_path_key)
    topology_candidates = [folder / f"{prefix}_slab_vampnet_start.pdb", folder / "start_L300.pdb"]
    topology = next(path for path in topology_candidates if path.exists())
    with md.formats.DCDTrajectoryFile(str(trajectories[-1])) as handle:
        last_index = len(handle) - 1
    frame = md.load_frame(str(trajectories[-1]), last_index, top=str(topology))
    xyz = np.asarray(frame.xyz[0], dtype=float)
    chain_len = xyz.shape[0] // 200
    chain_index = np.repeat(np.arange(200), chain_len)[: xyz.shape[0]]
    center = 0.5 * np.asarray(frame.unitcell_lengths[0], dtype=float)
    xyz = xyz - center
    return xyz[:, 0], xyz[:, 2], chain_index


def trajectory_means(data: dict[str, np.ndarray], key: str) -> np.ndarray:
    values = np.asarray(data[key], float)
    n_traj = int(np.asarray(data["n_traj"]).ravel()[0])
    usable = (values.size // n_traj) * n_traj
    return values[:usable].reshape(n_traj, -1).mean(axis=1)


def bootstrap_response(tl: np.ndarray, th: np.ndarray, rng: np.random.Generator, n_boot: int = 3000) -> np.ndarray:
    values = np.empty(n_boot)
    for idx in range(n_boot):
        a = tl[rng.integers(0, len(tl), len(tl))]
        b = th[rng.integers(0, len(th), len(th))]
        values[idx] = 100.0 * (np.mean(b) - np.mean(a)) / np.mean(a)
    return values


def cropped_image(path: Path) -> np.ndarray:
    image = mpimg.imread(path)
    rgb = image[..., :3]
    mask = np.any(rgb < 0.985, axis=2)
    if not np.any(mask):
        return image
    y, x = np.where(mask)
    pad = 8
    return image[max(0, y.min()-pad):min(image.shape[0], y.max()+pad), max(0, x.min()-pad):min(image.shape[1], x.max()+pad)]


def curve_from_csv(path: Path, system: str, role: str, ykey: str) -> tuple[np.ndarray, np.ndarray]:
    rows = [r for r in read_csv(path) if r["system"] == system and r["role"] == role]
    return np.asarray([float(r["time_ns"]) for r in rows]), np.asarray([float(r[ykey]) for r in rows])


def ratio_vs_threshold(metric: str, system: str, thresholds: np.ndarray, family: str = "condensate") -> np.ndarray:
    curves = []
    for role in ROLES:
        if metric == "rg":
            d = rg_cache(system, family, role); t, y = d["time_axis"], d["acf_rg"]
        elif metric == "intra":
            d = intra_cache(system, family, role); t, y = d["time_ns"], d["acf_mean"]
        elif metric == "inter":
            t, y = curve_from_csv(STORY / "interchain_contact_relaxation_curves_corrected.csv", system, role, "acf_mean")
        else:
            t, y = curve_from_csv(RESULTS / "interchain_partner_survival_curves.csv", system, role, "survival_mean")
        curves.append((t, y))
    tau_l = np.asarray([crossing(*curves[0], th) for th in thresholds])
    tau_h = np.asarray([crossing(*curves[1], th) for th in thresholds])
    return tau_h / tau_l


def summary_rows(metric: str) -> list[dict[str, str]]:
    if metric == "intra":
        return read_csv(STORY / "intra_inter_contact_relaxation/intrachain_contact_relaxation_summary.csv")
    if metric == "inter":
        return read_csv(STORY / "interchain_contact_relaxation_summary_corrected.csv")
    return read_csv(RESULTS / "interchain_partner_survival_summary.csv")


def parse_window_cache(path: Path) -> tuple[float, float]:
    data = np.load(path, allow_pickle=False)
    return float(np.asarray(data["relax_time_ns"])), float(np.asarray(data["last_part"]) - np.asarray(data["first_part"]) + 1)


def loss_gain_cache(system: str, role: str) -> dict[str, np.ndarray]:
    temp = TEMPS[system][role]
    path = next(iter(sorted(RESULTS.glob(f"partner_loss_gain_{system}_{temp}_parts*_pstride020.npz"))))
    data = np.load(path, allow_pickle=False)
    return {k: data[k] for k in data.files}


def numeric_path_key(path: Path) -> int:
    digits = "".join(ch for ch in path.stem if ch.isdigit())
    return int(digits or 0)


def estimate_t_from_qfiles(qfiles: list[Path], lag: int) -> np.ndarray:
    """Reproduce analyze_single.py's reversible soft-count transition matrix."""
    count = np.zeros((5, 5), dtype=np.float64)
    for path in qfiles:
        q = np.load(path, mmap_mode="r")
        if len(q) <= lag:
            continue
        a = np.asarray(q[:-lag], dtype=np.float64)
        b = np.asarray(q[lag:], dtype=np.float64)
        count += a.T @ b
    count = 0.5 * (count + count.T)
    row_sum = count.sum(axis=1, keepdims=True)
    transition = count / np.maximum(row_sum, 1e-12)
    bad = row_sum[:, 0] <= 1e-12
    transition[bad] = 0.0
    transition[bad, bad] = 1.0
    return transition


def model_ensemble(system: str, family: str, role: str, out: Path) -> dict[str, np.ndarray]:
    """Load or derive ITS, CK, T, and populations for all saved attempts."""
    cache_dir = out / "source_data/derived_model_ensembles"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{system}_{family}_{TEMPS[system][role]}_n5_lag50.npz"
    if cache_path.exists():
        data = np.load(cache_path, allow_pickle=False)
        return {key: data[key] for key in data.files}

    directory = model_dir(system, family, role)
    inter = directory / intermediate_name(system)
    h5 = inter / "data_combined_att0.hdf5"
    with h5py.File(h5, "r") as handle:
        group = handle["combined/0/5/50"]
        lags = np.asarray(group["lags_for_its"], dtype=int)
        dt = float(group.attrs.get("dt", 0.01))
        base_lag = int(group.attrs.get("ck_base_lag", 200))
        ck_steps = int(group.attrs.get("ck_steps", 4))

    attempt_dirs = sorted(
        inter.glob("qcache_combined_n5_att*_lag50"),
        key=lambda p: int(p.name.split("att")[-1].split("_")[0]),
    )
    transitions, populations, its_stack, cke_stack, ckp_stack, attempts = [], [], [], [], [], []
    for qdir in attempt_dirs:
        attempt = int(qdir.name.split("att")[-1].split("_")[0])
        map_path = directory / f"state_map_att{attempt}_raw_to_final.npy"
        qfiles = sorted(qdir.glob("traj*.npy"), key=numeric_path_key)
        if not map_path.exists() or not qfiles:
            continue

        total = np.zeros(5, dtype=np.float64)
        n_frames = 0
        for path in qfiles:
            q = np.load(path, mmap_mode="r")
            total += np.asarray(q, dtype=np.float64).sum(axis=0)
            n_frames += len(q)
        mu = total / max(n_frames, 1)

        t_by_lag = {int(lag): estimate_t_from_qfiles(qfiles, int(lag)) for lag in sorted(set(lags) | {base_lag * k for k in range(1, ck_steps + 1)})}
        its = np.empty((4, len(lags)), dtype=np.float64)
        for lag_idx, lag in enumerate(lags):
            eigenvalues = np.abs(np.linalg.eigvals(t_by_lag[int(lag)]))
            stationary = int(np.argmin(np.abs(eigenvalues - 1.0)))
            eigenvalues = np.delete(eigenvalues, stationary)
            eigenvalues = np.sort(np.clip(eigenvalues, 1e-12, 1 - 1e-12))[::-1]
            its[:, lag_idx] = -float(lag) * dt / np.log(eigenvalues[:4])

        t_base = t_by_lag[base_lag]
        ckp = np.asarray([np.linalg.matrix_power(t_base, k) for k in range(1, ck_steps + 1)])
        cke = np.asarray([t_by_lag[base_lag * k] for k in range(1, ck_steps + 1)])
        transitions.append(t_base)
        populations.append(mu)
        its_stack.append(its)
        cke_stack.append(cke)
        ckp_stack.append(ckp)
        attempts.append(attempt)

    payload = {
        "attempts": np.asarray(attempts, dtype=np.int32),
        "lags": lags.astype(np.int32),
        "dt": np.asarray(dt, dtype=np.float64),
        "base_lag": np.asarray(base_lag, dtype=np.int32),
        "t": np.asarray(transitions),
        "mu": np.asarray(populations),
        "its": np.asarray(its_stack),
        "cke": np.asarray(cke_stack),
        "ckp": np.asarray(ckp_stack),
    }
    np.savez_compressed(cache_path, **payload)
    print(f"[derive] {system} {family} {TEMPS[system][role]}: {len(attempts)} model attempts")
    return payload


def aligned_transition_stack(system: str, family: str, role: str, out: Path) -> tuple[np.ndarray, np.ndarray]:
    directory = model_dir(system, family, role)
    ensemble = model_ensemble(system, family, role, out)
    matrices, populations = [], []
    for raw_t, raw_mu, attempt in zip(ensemble["t"], ensemble["mu"], ensemble["attempts"]):
        attempt = int(attempt)
        map_path = directory / f"state_map_att{attempt}_raw_to_final.npy"
        if not map_path.exists():
            continue
        mapping = np.load(map_path).astype(int)
        final_t = np.zeros_like(raw_t); final_mu = np.zeros_like(raw_mu)
        for raw_i, final_i in enumerate(mapping):
            final_mu[final_i] = raw_mu[raw_i]
            for raw_j, final_j in enumerate(mapping):
                final_t[final_i, final_j] = raw_t[raw_i, raw_j]
        matrices.append(final_t); populations.append(final_mu)
    return np.asarray(matrices), np.asarray(populations)


def aligned_ck_stack(
    system: str,
    family: str,
    role: str,
    out: Path,
) -> tuple[np.ndarray, np.ndarray]:
    """Align estimated and predicted CK matrices to final within-condition states."""
    directory = model_dir(system, family, role)
    ensemble = model_ensemble(system, family, role, out)
    estimated, predicted = [], []
    for raw_est, raw_pred, attempt in zip(
        ensemble["cke"],
        ensemble["ckp"],
        ensemble["attempts"],
    ):
        map_path = directory / f"state_map_att{int(attempt)}_raw_to_final.npy"
        if not map_path.exists():
            continue
        mapping = np.load(map_path).astype(int)
        final_est = np.zeros_like(raw_est)
        final_pred = np.zeros_like(raw_pred)
        for raw_i, final_i in enumerate(mapping):
            for raw_j, final_j in enumerate(mapping):
                final_est[:, final_i, final_j] = raw_est[:, raw_i, raw_j]
                final_pred[:, final_i, final_j] = raw_pred[:, raw_i, raw_j]
        estimated.append(final_est)
        predicted.append(final_pred)
    return np.asarray(estimated), np.asarray(predicted)


def display_state_mapping(system: str, family: str, role: str) -> dict[int, int]:
    """Map final within-condition state ranks to cross-condition display states."""
    return {
        int(row["local_state"]): int(row.get("display_state") or row.get("global_state", 0))
        for row in state_rows(system, family, role)
    }


def state_rows(system: str, family: str, role: str) -> list[dict[str, str]]:
    temp = TEMPS[system][role]
    return [r for r in read_csv(STATE_SOURCE[system][family]) if r.get("temperature", r.get("temp", "")) == temp]


def mapped_matrix_stack(system: str, family: str, role: str, states: list[int], out: Path) -> tuple[np.ndarray, np.ndarray]:
    stack, pops = aligned_transition_stack(system, family, role, out)
    mapping = {int(r["local_state"]): int(r.get("display_state") or r.get("global_state", 0)) for r in state_rows(system, family, role)}
    index = {state: idx for idx, state in enumerate(states)}
    out_t = np.zeros((len(stack), len(states), len(states)))
    out_p = np.zeros((len(stack), len(states)))
    for attempt in range(len(stack)):
        for local_i, global_i in mapping.items():
            if global_i not in index: continue
            out_p[attempt, index[global_i]] += pops[attempt, local_i]
            for local_j, global_j in mapping.items():
                if global_j in index:
                    out_t[attempt, index[global_i], index[global_j]] = stack[attempt, local_i, local_j]
    return out_t, out_p


def mapped_distance_distribution(system: str, family: str, out: Path) -> tuple[np.ndarray, list[int]]:
    states = sorted({int(r.get("display_state") or r.get("global_state", 0)) for role in ROLES for r in state_rows(system, family, role)})
    tl_t, tl_p = mapped_matrix_stack(system, family, "Tl", states, out)
    th_t, th_p = mapped_matrix_stack(system, family, "Th", states, out)
    n = min(len(tl_t), len(th_t))
    pi_bar = 0.5 * (tl_p[:n] + th_p[:n])
    dist = np.sqrt(np.sum(pi_bar[:, :, None] * (th_t[:n] - tl_t[:n]) ** 2, axis=(1, 2)))
    return dist, states


def write_readme(out: Path) -> None:
    text = """# Supporting figures

Batch entry point: `figures/make_supplementary_figures.py`.
Implementation modules are in `figures/supplementary/`.

Examples:

```bash
python figures/make_supplementary_figures.py --only 1,3,4
python figures/make_s2_s3_structural_controls.py
```

1. S1: condensate density profiles at low and high temperature.
2. S2: intrachain contact-number distributions and heating response.
3. S3: full state-resolved intrachain contact maps.
4. S4: low-dimensional state-space projections and assignments.
5. S5: VAMPNet implied-timescale validation across 20 models.
6. S6: state-resolved Chapman-Kolmogorov tests.
7. S7: end-to-end-distance autocorrelations and heating response.
8. S8: partner turnover, partner-number change, and compensation.

Important scope notes:

- TDP-43-LCD uses 260 K and 300 K throughout.
- S6 uses a per-system pooled global TICA basis. TIC1 and TIC2 are the first
  two slow coordinates after PCA compression of intrachain contact features;
  they are projection coordinates, not physical spatial axes.
- S5 and S8 use the same cross-condition state identities as Figure 3.
- Figure 3B reports distances between ensemble-mean transition matrices.
  Model ensembles are resampled independently across temperatures and
  environments. Direct 95% bootstrap intervals for the condensate-minus-free
  contrast exclude zero for FUS-LCD, A1-LCD, and TDP-43-LCD.
- In S2, heating raises fractional partner turnover but modestly lowers the
  compensation index. The defensible statement is that substantial loss-gain
  compensation persists at both temperatures, not that heating strengthens it.
- S3 uses the same 0.66 nm interchain bead-contact observable and matched
  trajectory windows as Figure 2D. Panel F reports independent segment-bootstrap
  confidence intervals for the high-minus-low temperature differences rather
  than comparing the absolute scales of load-state persistence and partner
  survival. The load-state point estimates are positive throughout, although
  the FUS-LCD interval includes zero at 20 ns; partner-survival differences are
  negative throughout. These results explain coexistence but do not establish
  redistribution as a strict causal origin of the slower connected ACF.
"""
    (out / "README.md").write_text(text, encoding="utf-8")

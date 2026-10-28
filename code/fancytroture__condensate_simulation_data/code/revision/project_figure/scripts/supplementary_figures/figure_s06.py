"""Generate supplementary Figure S6: global TICA state-space projections."""

from __future__ import annotations

import re
from scipy.ndimage import gaussian_filter

from .shared import *  # noqa: F403


SOURCE_CSV = {
    "FUS": RESULTS / "tica_state_projection_fus_free_condensate_source.csv",
    "A1": RESULTS / "tica_state_projection_a1_free_condensate_source.csv",
    "TDP43": RESULTS / "tica_state_projection_tdp43_free_condensate_source.csv",
}

FREE_STATE_COLORS = ["#1f77b4", "#4c78a8", "#72b7b2", "#54a24b", "#9ecae9", "#b279a2"]
COND_STATE_COLORS = ["#dd8452", "#ff9d76", "#f2cf5b", "#8ab17d", "#c44e52", "#b07aa1"]
STATE_SUP = {"free": "Free", "condensate": "Cond"}
SAMPLE_PER_CONDITION = 2200
RNG_SEED = 20240511


def state_color(family: str, display_state: int) -> str:
    palette = FREE_STATE_COLORS if family == "free" else COND_STATE_COLORS
    return palette[int(display_state) % len(palette)]


def state_label(family: str, display_state: int) -> str:
    return r"$s_{{{}}}$".format(int(display_state) + 1)


def read_projection_rows(system: str) -> list[dict[str, str]]:
    rows = read_csv(SOURCE_CSV[system])
    compact = {
        "FUS": {0: 0, 1: 1, 2: 2},
        "A1": {0: 0, 1: 1, 2: 2, 5: 3},
        "TDP43": {0: 0, 1: 1, 2: 2, 5: 3, 6: 4},
    }[system]
    for row in rows:
        raw = int(row["display_state"])
        row["display_state"] = str(compact[raw])
    return rows


def rows_for_condition(rows: list[dict[str, str]], family: str, role: str) -> list[dict[str, str]]:
    selected = [row for row in rows if row["family"] == family and row["role"] == role]
    if not selected:
        raise RuntimeError(f"No projection rows for {family} {role}")
    return selected


def load_global_model(model_dir: Path) -> dict[str, np.ndarray]:
    return {
        "feature_mean": np.load(model_dir / "feature_mean.npy"),
        "feature_scale": np.load(model_dir / "feature_scale.npy"),
        "pca_mean": np.load(model_dir / "pca_mean.npy"),
        "pca_components": np.load(model_dir / "pca_components.npy"),
        "tica_eigenvectors": np.load(model_dir / "tica_eigenvectors.npy"),
    }


def transform_features_to_tic(
    feature_path: Path,
    flat_indices: np.ndarray,
    model: dict[str, np.ndarray],
    chunk_size: int = 1500,
) -> np.ndarray:
    features = np.load(feature_path, mmap_mode="r")
    coords = np.empty((flat_indices.size, 2), dtype=np.float32)
    mean = np.asarray(model["feature_mean"], dtype=np.float32)
    scale = np.asarray(model["feature_scale"], dtype=np.float32)
    pca_mean = np.asarray(model["pca_mean"], dtype=np.float32)
    pca_components = np.asarray(model["pca_components"], dtype=np.float32)
    tica_vectors = np.asarray(model["tica_eigenvectors"][:, :2], dtype=np.float32)

    for start in range(0, flat_indices.size, chunk_size):
        stop = min(start + chunk_size, flat_indices.size)
        block = np.asarray(features[flat_indices[start:stop]], dtype=np.float32)
        block -= mean
        block /= scale
        block -= pca_mean
        z_pca = block @ pca_components.T
        coords[start:stop] = z_pca @ tica_vectors
    return coords


def sampled_flat_indices(feature_path: Path, system: str, family: str, role: str) -> np.ndarray:
    n_total = int(np.load(feature_path, mmap_mode="r").shape[0])
    sample_size = min(SAMPLE_PER_CONDITION, n_total)
    seed_offset = 1000 * SYSTEMS.index(system) + 200 * FAMILIES.index(family) + 31 * ROLES.index(role)
    rng = np.random.default_rng(RNG_SEED + seed_offset)
    return np.sort(rng.choice(n_total, size=sample_size, replace=False)).astype(np.int64)


def raw_to_final_map(tica_dir: Path, att: int, n_states: int) -> np.ndarray:
    path = tica_dir / f"state_map_att{att}_raw_to_final.npy"
    if not path.exists():
        return np.arange(n_states, dtype=np.int32)
    mapper = np.asarray(np.load(path), dtype=np.int32).ravel()
    if mapper.size < n_states:
        return np.arange(n_states, dtype=np.int32)
    return mapper[:n_states]


def condition_display_labels(
    flat_indices: np.ndarray,
    qcache_dir: Path,
    local_to_display: dict[int, int],
) -> np.ndarray:
    inter_dir = qcache_dir.parent
    lengths_path = inter_dir / "combined_tica_lengths.npy"
    if lengths_path.exists():
        lengths = np.asarray(np.load(lengths_path), dtype=np.int64)
    else:
        q_files = sorted(qcache_dir.glob("traj*.npy"))
        lengths = np.asarray([np.load(path, mmap_mode="r").shape[0] for path in q_files], dtype=np.int64)
    starts = np.concatenate([[0], np.cumsum(lengths)])

    match = re.search(r"_att(\d+)_lag", qcache_dir.name)
    att = int(match.group(1)) if match else 0
    first_q = next(iter(sorted(qcache_dir.glob("traj*.npy"))))
    n_states = int(np.load(first_q, mmap_mode="r").shape[1])
    mapper = raw_to_final_map(qcache_dir.parent.parent, att, n_states)

    labels = np.full(flat_indices.size, -1, dtype=np.int32)
    traj_ids = np.searchsorted(starts[1:], flat_indices, side="right")
    for traj_id in np.unique(traj_ids):
        mask = traj_ids == traj_id
        local_idx = flat_indices[mask] - starts[int(traj_id)]
        q = np.load(qcache_dir / f"traj{int(traj_id)}.npy", mmap_mode="r")
        raw = np.argmax(q[local_idx], axis=1).astype(np.int32)
        final = mapper[raw]
        labels[mask] = np.asarray([local_to_display.get(int(state), -1) for state in final], dtype=np.int32)
    return labels


def load_condition_payload(system: str, family: str, role: str, rows: list[dict[str, str]]) -> dict[str, object]:
    condition_rows = rows_for_condition(rows, family, role)
    first = condition_rows[0]
    feature_path = Path(first["features_file"])
    qcache_dir = Path(first["qcache_dir"])
    local_to_display = {int(row["local_state"]): int(row["display_state"]) for row in condition_rows}
    local_to_mu = {int(row["local_state"]): float(row["mu"]) for row in condition_rows}
    cache_dir = DEFAULT_OUT / "source_data/s6_state_space_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{system}_{family}_{role}_sample{SAMPLE_PER_CONDITION}_seed{RNG_SEED}_compact_v2.npz"

    if cache_path.exists():
        with np.load(cache_path, allow_pickle=False) as cached:
            coords = np.asarray(cached["coords_all"], dtype=np.float32)
            labels = np.asarray(cached["labels"], dtype=np.int32)
        keep = labels >= 0
        print(f"[cache] S6 {system} {family} {role}", flush=True)
    else:
        print(f"[project] S6 {system} {family} {role}", flush=True)
        model = load_global_model(Path(first["global_tica_model_dir"]))
        flat_indices = sampled_flat_indices(feature_path, system, family, role)
        coords = transform_features_to_tic(feature_path, flat_indices, model)
        labels = condition_display_labels(flat_indices, qcache_dir, local_to_display)
        np.savez_compressed(cache_path, coords_all=coords, labels=labels)
        keep = labels >= 0

    return {
        "system": system,
        "family": family,
        "role": role,
        "temp": first["temp"],
        "coords_all": coords,
        "coords_state": coords[keep],
        "labels": labels[keep],
        "local_to_display": local_to_display,
        "local_to_mu": local_to_mu,
        "source": first,
    }


def density_grid(
    coords: np.ndarray,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    bins: int = 135,
    sigma: float = 2.2,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    hist, xedges, yedges = np.histogram2d(coords[:, 0], coords[:, 1], bins=bins, range=[xlim, ylim])
    return gaussian_filter(hist.astype(np.float64), sigma=sigma).T, xedges, yedges


def draw_background(ax, coords: np.ndarray, xlim: tuple[float, float], ylim: tuple[float, float]) -> None:
    hist, xedges, yedges = density_grid(coords, xlim, ylim, bins=135, sigma=2.6)
    vmax = np.percentile(hist[hist > 0], 98.5) if np.any(hist > 0) else 1.0
    ax.imshow(
        hist,
        origin="lower",
        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
        cmap="Greys",
        vmin=0.0,
        vmax=vmax,
        alpha=0.30,
        interpolation="bilinear",
        aspect="auto",
        zorder=0,
    )


def draw_state_contours(
    ax,
    payload: dict[str, object],
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> list[dict[str, object]]:
    coords = np.asarray(payload["coords_state"], dtype=np.float64)
    labels = np.asarray(payload["labels"], dtype=np.int32)
    family = str(payload["family"])
    rows = []
    for display_state in sorted(np.unique(labels).tolist()):
        points = coords[labels == display_state]
        if points.shape[0] < 35:
            continue
        hist, xedges, yedges = density_grid(points, xlim, ylim, bins=135, sigma=2.15)
        peak = float(np.nanmax(hist))
        if peak <= 0:
            continue
        xcenters = 0.5 * (xedges[:-1] + xedges[1:])
        ycenters = 0.5 * (yedges[:-1] + yedges[1:])
        color = state_color(family, int(display_state))
        ax.contourf(
            xcenters,
            ycenters,
            hist,
            levels=[peak * 0.18, peak * 0.42, peak],
            colors=[color, color],
            alpha=0.18,
            antialiased=True,
            zorder=1,
        )
        ax.contour(
            xcenters,
            ycenters,
            hist,
            levels=[peak * 0.18, peak * 0.42],
            colors=[color, color],
            linewidths=[0.9, 1.45],
            alpha=0.95,
            zorder=2,
        )
        cx = float(np.median(points[:, 0]))
        cy = float(np.median(points[:, 1]))
        ax.scatter([cx], [cy], s=24, color=color, edgecolor="white", linewidth=0.55, zorder=3)
        ax.text(
            cx,
            cy,
            state_label(family, int(display_state)),
            fontsize=16.0,
            fontweight="bold",
            ha="center",
            va="center",
            color="#202020",
            zorder=4,
        )
        rows.append(
            {
                "system": payload["system"],
                "family": family,
                "role": payload["role"],
                "temperature": payload["temp"],
                "display_state": int(display_state),
                "n_points": int(points.shape[0]),
                "tic1_center": f"{cx:.6f}",
                "tic2_center": f"{cy:.6f}",
            }
        )
    return rows


def axis_limits(payloads: list[dict[str, object]]) -> dict[str, tuple[tuple[float, float], tuple[float, float]]]:
    limits = {}
    for system in SYSTEMS:
        coords = np.concatenate(
            [np.asarray(payload["coords_all"], dtype=np.float64) for payload in payloads if payload["system"] == system],
            axis=0,
        )
        xmin, xmax = np.nanpercentile(coords[:, 0], [0.5, 99.5])
        ymin, ymax = np.nanpercentile(coords[:, 1], [0.5, 99.5])
        xpad = max((float(xmax) - float(xmin)) * 0.12, 0.2)
        ypad = max((float(ymax) - float(ymin)) * 0.12, 0.2)
        limits[system] = ((float(xmin) - xpad, float(xmax) + xpad), (float(ymin) - ypad, float(ymax) + ypad))
    return limits


def figure_s6_tica(out: Path) -> None:
    payloads = []
    for system in SYSTEMS:
        rows = read_projection_rows(system)
        for family, role in (("free", "Tl"), ("free", "Th"), ("condensate", "Tl"), ("condensate", "Th")):
            payloads.append(load_condition_payload(system, family, role, rows))

    limits = axis_limits(payloads)
    payload_index = {(p["system"], p["family"], p["role"]): p for p in payloads}
    conditions = (("free", "Tl"), ("free", "Th"), ("condensate", "Tl"), ("condensate", "Th"))

    fig, axes = plt.subplots(3, 4, figsize=(15.4, 10.9), sharex=False, sharey=False)
    source_rows = []
    for row_idx, system in enumerate(SYSTEMS):
        xlim, ylim = limits[system]
        for col_idx, (family, role) in enumerate(conditions):
            ax = axes[row_idx, col_idx]
            payload = payload_index[(system, family, role)]
            draw_background(ax, np.asarray(payload["coords_all"], dtype=np.float64), xlim, ylim)
            source_rows.extend(draw_state_contours(ax, payload, xlim, ylim))
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
            ax.text(
                0.96,
                0.95,
                temp_label(role),
                transform=ax.transAxes,
                fontsize=18.0,
                fontweight="bold",
                ha="right",
                va="top",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.2},
                zorder=6,
            )
            ax.grid(color="#d9d9d9", alpha=0.28, linestyle="--", linewidth=0.55)
            ax.set_axisbelow(True)
            for spine in ax.spines.values():
                spine.set_color("#3a3a3a")
                spine.set_linewidth(0.75)
            ax.tick_params(labelsize=14.0, length=4.2, width=0.9, pad=2.2)
            if col_idx == 0:
                ax.set_ylabel("TIC2", fontsize=16.0)
            else:
                ax.set_yticklabels([])
            if row_idx == 2:
                ax.set_xlabel("TIC1", fontsize=16.0)
            else:
                ax.set_xticklabels([])

        axes[row_idx, 0].text(
            -0.24,
            0.5,
            DISPLAY[system],
            transform=axes[row_idx, 0].transAxes,
            rotation=90,
            ha="center",
            va="center",
            fontsize=21.0,
            fontweight="bold",
        )

    fig.subplots_adjust(left=0.075, right=0.990, bottom=0.055, top=0.900, hspace=0.22, wspace=0.105)
    fig.canvas.draw()

    group_titles = [("Free", 0, 1), ("Condensate", 2, 3)]
    for title, first_col, last_col in group_titles:
        box_left = axes[0, first_col].get_position()
        box_right = axes[0, last_col].get_position()
        fig.text(
            0.5 * (box_left.x0 + box_right.x1),
            0.947,
            title,
            ha="center",
            va="bottom",
            fontsize=22.0,
            fontweight="bold",
        )

    write_csv(out / "source_data/FigureS6_state_space_projection_summary.csv", source_rows)
    save(fig, out, "FigureS6_state_space_projections_and_assignments")

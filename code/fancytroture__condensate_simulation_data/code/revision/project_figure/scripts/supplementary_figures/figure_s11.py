"""Generate supplementary Figure S11: contact-definition sensitivity."""

from .shared import *  # noqa: F403


CUTOFFS = (0.60, 0.66, 0.72)
CUTOFF_CACHE_DIRS = {
    0.60: DEFAULT_OUT / "source_data/contact_cutoff_sensitivity/cutoff_060",
    0.66: DEFAULT_OUT / "source_data/s4_contact_load_cache",
    0.72: DEFAULT_OUT / "source_data/contact_cutoff_sensitivity/cutoff_072",
}
CUTOFF_ANALYSIS_DIRS = {
    0.60: DEFAULT_OUT / "source_data/contact_cutoff_sensitivity/cutoff_060_analysis",
    0.66: DEFAULT_OUT / "source_data/s4_contact_load_analysis",
    0.72: DEFAULT_OUT / "source_data/contact_cutoff_sensitivity/cutoff_072_analysis",
}


def cutoff_cache(cutoff: float, system: str, role: str) -> dict[str, np.ndarray]:
    temp = TEMPS[system][role]
    matches = sorted(
        CUTOFF_CACHE_DIRS[cutoff].glob(
            f"partner_contact_load_{system}_{temp}_parts*_pstride005.npz"
        )
    )
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one cutoff cache for {cutoff:.2f} nm {system} {role}; found {len(matches)}"
        )
    data = np.load(matches[0], allow_pickle=False)
    return {key: data[key] for key in data.files}


def connected_acf(values: np.ndarray) -> np.ndarray | None:
    array = np.asarray(values, dtype=float)
    if len(array) < 4:
        return None
    array = array - np.mean(array)
    variance = float(np.dot(array, array) / len(array))
    if not np.isfinite(variance) or variance <= 1e-12:
        return None
    transform = np.fft.fft(array, n=2 * len(array))
    correlation = np.fft.ifft(transform * np.conjugate(transform))[:len(array)].real
    return correlation / (variance * np.arange(len(array), 0, -1, dtype=float))


def segment_curves(
    cache: dict[str, np.ndarray],
    metric: str,
    max_lag_ns: float,
) -> tuple[np.ndarray, np.ndarray]:
    loads = np.asarray(cache["contact_load"], dtype=np.int32)
    valid = np.asarray(cache["valid"], dtype=bool)
    dt = float(np.asarray(cache["dt_ns"]).reshape(-1)[0])
    max_lag = int(round(max_lag_ns / dt))
    curves = []
    for segment_load, segment_valid in zip(loads, valid):
        values = segment_load[segment_valid]
        if len(values) <= max_lag:
            continue
        if metric == "inter":
            curve = connected_acf(np.sum(values, axis=1))
            if curve is None:
                continue
            curves.append(curve[:max_lag + 1])
            continue
        present = values > 0
        survival = np.full(max_lag + 1, np.nan)
        survival[0] = 1.0
        for lag in range(1, max_lag + 1):
            initial = present[:-lag]
            future = present[lag:]
            initial_count = np.sum(initial, axis=1)
            usable = initial_count > 0
            if np.any(usable):
                retained = np.sum(initial & future, axis=1)
                survival[lag] = float(np.mean(retained[usable] / initial_count[usable]))
        curves.append(survival)
    return np.arange(max_lag + 1) * dt, np.asarray(curves)


def bootstrap_ratio(
    tl_time: np.ndarray,
    tl_curves: np.ndarray,
    th_time: np.ndarray,
    th_curves: np.ndarray,
    threshold: float,
    rng: np.random.Generator,
    n_boot: int = 2000,
) -> tuple[float, float, float]:
    tl_cross = crossing(tl_time, np.nanmean(tl_curves, axis=0), threshold)
    th_cross = crossing(th_time, np.nanmean(th_curves, axis=0), threshold)
    point = th_cross / tl_cross
    samples = []
    for _ in range(n_boot):
        tl = np.nanmean(tl_curves[rng.integers(0, len(tl_curves), len(tl_curves))], axis=0)
        th = np.nanmean(th_curves[rng.integers(0, len(th_curves), len(th_curves))], axis=0)
        low_cross = crossing(tl_time, tl, threshold)
        high_cross = crossing(th_time, th, threshold)
        if np.isfinite(low_cross) and np.isfinite(high_cross) and low_cross > 0:
            samples.append(high_cross / low_cross)
    low, high = np.percentile(samples, [2.5, 97.5])
    return float(point), float(low), float(high)


def stable_load_ratio(cutoff: float, system: str) -> tuple[float, float, float]:
    path = CUTOFF_ANALYSIS_DIRS[cutoff] / "FigureS4_contact_load_mechanism_summary.csv"
    rows = [
        row for row in read_csv(path)
        if row["system"] == system
        and row["metric"] == "conditional_replacement_eps0.05"
        and np.isclose(float(row["lag_ns"]), 2.0)
    ]
    by_role = {row["role"]: row for row in rows}
    tl = by_role["Tl"]
    th = by_role["Th"]
    point = float(th["value"]) / float(tl["value"])
    low = float(th["ci_low"]) / float(tl["ci_high"])
    high = float(th["ci_high"]) / float(tl["ci_low"])
    return point, low, high


def figure_s11_contact_sensitivity(out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.5))
    system_colors = {"FUS": "#3f6c8e", "A1": "#b56a32", "TDP43": "#69619b"}
    markers = {"FUS": "o", "A1": "s", "TDP43": "^"}
    rng = np.random.default_rng(20260723)
    source: list[dict[str, object]] = []

    for system in SYSTEMS:
        inter_values = []
        partner_values = []
        replacement_values = []
        for cutoff in CUTOFFS:
            tl_cache = cutoff_cache(cutoff, system, "Tl")
            th_cache = cutoff_cache(cutoff, system, "Th")
            inter_tl = segment_curves(tl_cache, "inter", 30.0)
            inter_th = segment_curves(th_cache, "inter", 30.0)
            partner_tl = segment_curves(tl_cache, "partner", 50.0)
            partner_th = segment_curves(th_cache, "partner", 50.0)
            inter = bootstrap_ratio(*inter_tl, *inter_th, 0.20, rng)
            partner = bootstrap_ratio(*partner_tl, *partner_th, 0.50, rng)
            replacement = stable_load_ratio(cutoff, system)
            inter_values.append(inter)
            partner_values.append(partner)
            replacement_values.append(replacement)
            source.append({
                "system": DISPLAY[system],
                "cutoff_nm": cutoff,
                "interchain_t0.2_ratio": inter[0],
                "interchain_ci_low": inter[1],
                "interchain_ci_high": inter[2],
                "partner_t50_ratio": partner[0],
                "partner_ci_low": partner[1],
                "partner_ci_high": partner[2],
                "stable_load_replacement_ratio_2ns": replacement[0],
                "replacement_ci_low": replacement[1],
                "replacement_ci_high": replacement[2],
            })

        for ax, values in zip(axes, (inter_values, partner_values, replacement_values)):
            point = np.asarray([value[0] for value in values])
            low = np.asarray([value[1] for value in values])
            high = np.asarray([value[2] for value in values])
            ax.fill_between(CUTOFFS, low, high, color=system_colors[system], alpha=0.10)
            ax.plot(CUTOFFS, point, marker=markers[system], ms=5.5, lw=1.5,
                    color=system_colors[system], label=DISPLAY[system])

    for ax in axes:
        ax.axhline(1, color="#777777", ls="--", lw=0.9)
        ax.axvline(0.66, color="#999999", ls=":", lw=0.9)
        ax.set_xticks(CUTOFFS, [f"{cutoff:.2f}" for cutoff in CUTOFFS])
        ax.set_xlabel("Contact cutoff (nm)")
        clean(ax)
    axes[0].set_ylabel(r"Timescale ratio  $\tau(T_{\mathrm{H}})/\tau(T_{\mathrm{L}})$")
    axes[0].set_title("Interchain contact number", fontweight="bold")
    axes[1].set_title("Partner identity", fontweight="bold")
    axes[2].set_title("Stable-load partner replacement", fontweight="bold")
    axes[2].set_ylabel(r"Replacement ratio  $T_{\mathrm{H}}/T_{\mathrm{L}}$")
    axes[2].legend(loc="upper right", fontsize=9)
    panel(axes[0], "A")
    panel(axes[1], "B")
    panel(axes[2], "C")
    fig.tight_layout(w_pad=1.7)
    write_csv(out / "source_data/FigureS11_contact_definition_sensitivity.csv", source)
    save(fig, out, "FigureS11_contact_definition_sensitivity")

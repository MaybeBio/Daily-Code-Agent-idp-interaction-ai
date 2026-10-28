"""Generate supplementary Figure S10: timescale-threshold robustness."""

from .shared import *  # noqa: F403


def inter_cache(system: str, role: str) -> dict[str, np.ndarray]:
    path = next(
        iter(sorted(
            (ROOT / "project_figure/source_data/interchain_ratio_bootstrap").glob(
                f"interchain_contact_relaxation_{system}_{TEMPS[system][role]}_parts*.npz"
            )
        ))
    )
    data = np.load(path, allow_pickle=False)
    return {key: data[key] for key in data.files}


def partner_cache(system: str, role: str) -> dict[str, np.ndarray]:
    path = next(
        iter(sorted(
            RESULTS.glob(
                f"interchain_partner_survival_{system}_{TEMPS[system][role]}_parts*_pstride020_lag50p0ns.npz"
            )
        ))
    )
    data = np.load(path, allow_pickle=False)
    return {key: data[key] for key in data.files}


def ratio_curve(
    tl_time: np.ndarray,
    tl_curves: np.ndarray,
    th_time: np.ndarray,
    th_curves: np.ndarray,
    thresholds: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = 1500,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    point = np.asarray([
        crossing(th_time, np.mean(th_curves, axis=0), threshold)
        / crossing(tl_time, np.mean(tl_curves, axis=0), threshold)
        for threshold in thresholds
    ])
    bootstrap = np.full((n_boot, len(thresholds)), np.nan)
    for sample in range(n_boot):
        tl = np.mean(tl_curves[rng.integers(0, len(tl_curves), len(tl_curves))], axis=0)
        th = np.mean(th_curves[rng.integers(0, len(th_curves), len(th_curves))], axis=0)
        for idx, threshold in enumerate(thresholds):
            tl_cross = crossing(tl_time, tl, threshold)
            th_cross = crossing(th_time, th, threshold)
            if np.isfinite(tl_cross) and np.isfinite(th_cross) and tl_cross > 0:
                bootstrap[sample, idx] = th_cross / tl_cross
    low = np.nanpercentile(bootstrap, 2.5, axis=0)
    high = np.nanpercentile(bootstrap, 97.5, axis=0)
    return point, low, high


def figure_s10_threshold_robustness(out: Path) -> None:
    metrics = (
        ("intra", "Intrachain contact number", np.linspace(0.15, 0.30, 9), 0.20),
        ("partner", "Partner identity", np.linspace(0.425, 0.60, 8), 0.50),
        ("inter", "Interchain contact number", np.linspace(0.15, 0.30, 9), 0.20),
    )
    system_colors = {"FUS": "#3f6c8e", "A1": "#b56a32", "TDP43": "#69619b"}
    markers = {"FUS": "o", "A1": "s", "TDP43": "^"}
    fig, axes = plt.subplots(1, 3, figsize=(13.3, 4.4))
    source: list[dict[str, object]] = []
    rng = np.random.default_rng(20260723)

    for ax, (metric, title, thresholds, main_threshold) in zip(axes, metrics):
        for system in SYSTEMS:
            if metric == "intra":
                tl = intra_cache(system, "condensate", "Tl")
                th = intra_cache(system, "condensate", "Th")
                ratio = np.asarray([
                    crossing(th["time_ns"], th["acf_mean"], threshold)
                    / crossing(tl["time_ns"], tl["acf_mean"], threshold)
                    for threshold in thresholds
                ])
                low = np.full_like(ratio, np.nan)
                high = np.full_like(ratio, np.nan)
            else:
                tl = partner_cache(system, "Tl") if metric == "partner" else inter_cache(system, "Tl")
                th = partner_cache(system, "Th") if metric == "partner" else inter_cache(system, "Th")
                curve_key = "segment_survival" if metric == "partner" else "segment_acfs"
                time_key = "time_ns"
                ratio, low, high = ratio_curve(
                    tl[time_key], tl[curve_key], th[time_key], th[curve_key],
                    thresholds, rng,
                )
                ax.fill_between(thresholds, low, high, color=system_colors[system], alpha=0.10)
            ax.plot(
                thresholds,
                ratio,
                marker=markers[system],
                ms=5.2,
                lw=1.5,
                color=system_colors[system],
                label=DISPLAY[system],
            )
            for threshold, value, ci_low, ci_high in zip(thresholds, ratio, low, high):
                source.append({
                    "metric": metric,
                    "system": DISPLAY[system],
                    "threshold": float(threshold),
                    "timescale_ratio_th_over_tl": float(value),
                    "bootstrap_ci_low": float(ci_low),
                    "bootstrap_ci_high": float(ci_high),
                })
        ax.axhline(1, color="#777777", ls="--", lw=0.9)
        ax.axvline(main_threshold, color="#9a9a9a", ls=":", lw=0.9)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.set_xlabel("Crossing threshold")
        clean(ax)

    axes[0].set_ylabel(r"Timescale ratio  $\tau(T_{\mathrm{H}})/\tau(T_{\mathrm{L}})$")
    axes[0].set_ylim(0, 1.15)
    axes[1].set_ylim(0, 1.15)
    inter_high = [
        row["bootstrap_ci_high"]
        for row in source
        if row["metric"] == "inter" and np.isfinite(row["bootstrap_ci_high"])
    ]
    axes[2].set_ylim(0, max(3.6, 1.05 * max(inter_high, default=3.4)))
    axes[2].legend(loc="upper right", fontsize=9)
    panel(axes[0], "A")
    panel(axes[1], "B")
    panel(axes[2], "C")
    fig.tight_layout(w_pad=1.6)
    write_csv(out / "source_data/FigureS10_threshold_robustness.csv", source)
    save(fig, out, "FigureS10_relaxation_threshold_robustness")

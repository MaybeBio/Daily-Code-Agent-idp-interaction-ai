"""Generate supplementary Figure S7: implied-timescale validation."""

from .shared import *  # noqa: F403


def _nice_log_floor(value: float) -> float:
    if not np.isfinite(value) or value <= 0:
        return 1e-2
    exponent = np.floor(np.log10(value))
    base = 10.0 ** exponent
    for step in (10.0, 7.0, 5.0, 3.0, 2.0, 1.0):
        candidate = step * base
        if candidate <= value:
            return candidate
    return base


def _nice_log_ceil(value: float) -> float:
    if not np.isfinite(value) or value <= 0:
        return 1.0
    exponent = np.floor(np.log10(value))
    base = 10.0 ** exponent
    for step in (1.0, 2.0, 3.0, 5.0, 7.0, 10.0):
        candidate = step * base
        if candidate >= value:
            return candidate
    return 10.0 * base


def figure_s7_its(out: Path) -> None:
    fig, axes = plt.subplots(3, 4, figsize=(15.2, 10.2), sharex=True)
    panel_data = {}
    y_values = []
    for system in SYSTEMS:
        for family, role in (("free", "Tl"), ("free", "Th"), ("condensate", "Tl"), ("condensate", "Th")):
            ensemble = model_ensemble(system, family, role, out)
            stack = ensemble["its"]
            x = ensemble["lags"] * float(ensemble["dt"])
            selected_lag = int(ensemble["base_lag"]) * float(ensemble["dt"])
            mean = np.nanmean(stack, axis=0)
            lo = np.nanpercentile(stack, 5, axis=0)
            hi = np.nanpercentile(stack, 95, axis=0)
            panel_data[(system, family, role)] = (x, selected_lag, mean, lo, hi)
            finite = np.r_[x[np.isfinite(x) & (x > 0)], mean[np.isfinite(mean) & (mean > 0)], lo[np.isfinite(lo) & (lo > 0)], hi[np.isfinite(hi) & (hi > 0)]]
            if finite.size:
                y_values.append(finite)
    y_all = np.concatenate(y_values)
    y_limits = (_nice_log_floor(float(np.nanmin(y_all)) * 0.9), _nice_log_ceil(float(np.nanmax(y_all)) * 1.08))

    rows = []
    for i, system in enumerate(SYSTEMS):
        for j, (family, role) in enumerate((("free", "Tl"), ("free", "Th"), ("condensate", "Tl"), ("condensate", "Th"))):
            ax = axes[i, j]
            x, selected_lag, mean, lo, hi = panel_data[(system, family, role)]
            colors = plt.cm.viridis(np.linspace(0.12, 0.88, mean.shape[0]))
            ax.plot(x, x, "--", color="#777777", lw=1.1)
            ax.axvline(selected_lag, color="#222222", ls="--", lw=1.4, alpha=0.85, zorder=1)
            for mode, color in enumerate(colors):
                valid = np.isfinite(mean[mode]) & (mean[mode] > 0)
                ax.fill_between(x[valid], lo[mode, valid], hi[mode, valid], color=color, alpha=0.16)
                ax.plot(x[valid], mean[mode, valid], "o-", color=color, ms=3.4, lw=1.35)
                for lag, value, low, high in zip(x[valid], mean[mode, valid], lo[mode, valid], hi[mode, valid]):
                    rows.append({"system": DISPLAY[system], "family": family, "role": role, "mode": mode + 2, "lag_ns": lag, "its_mean_ns": value, "p05": low, "p95": high, "selected_lag_ns": selected_lag})
            ax.set_yscale("log")
            ax.set_ylim(*y_limits)
            ax.text(
                0.95,
                0.07,
                temp_label(role),
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=18.0,
                fontweight="bold",
                zorder=10,
            )
            if i == 2: ax.set_xlabel("Lag time (ns)", fontsize=17.0)
            if j == 0: ax.set_ylabel("Implied timescale (ns)", fontsize=17.0)
            ax.tick_params(axis="both", which="major", labelsize=15.0, width=1.1, length=5.0)
            ax.tick_params(axis="both", which="minor", width=0.8, length=2.8)
            if j > 0:
                ax.tick_params(axis="y", which="both", labelleft=False)
            ax.grid(True, which="both", alpha=0.18, lw=0.65)
        axes[i, 0].text(
            -0.43,
            0.5,
            DISPLAY[system],
            transform=axes[i, 0].transAxes,
            rotation=90,
            ha="center",
            va="center",
            fontsize=17.0,
            fontweight="bold",
        )
    fig.subplots_adjust(left=0.105, right=0.990, bottom=0.070, top=0.900, hspace=0.27, wspace=0.235)
    fig.canvas.draw()

    for title, first_col, last_col in (("Free", 0, 1), ("Condensate", 2, 3)):
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

    write_csv(out / "source_data/FigureS7_its_model_ensemble.csv", rows)
    save(fig, out, "FigureS7_vampnet_implied_timescale_validation")

"""Generate supplementary Figure S9: transition-distance model ensemble."""

from .shared import *  # noqa: F403


def figure_s9_transition_robustness(out: Path) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(6.6, 4.8))
    rows: list[dict[str, object]] = []
    x = np.arange(len(SYSTEMS), dtype=float)

    for idx, system in enumerate(SYSTEMS):
        free, free_states = mapped_distance_distribution(system, "free", out)
        condensate, condensate_states = mapped_distance_distribution(system, "condensate", out)
        n = min(len(free), len(condensate))
        jitter = np.linspace(-0.045, 0.045, n)
        ax.plot(
            np.full(n, x[idx] - 0.12) + jitter,
            free[:n],
            "o",
            ms=3.2,
            mfc="white",
            mec="#4a4a4a",
            mew=0.75,
            color="#4a4a4a",
            alpha=0.38,
        )
        ax.plot(
            np.full(n, x[idx] + 0.12) + jitter,
            condensate[:n],
            "o",
            ms=3.2,
            mfc="#4a4a4a",
            mec="#4a4a4a",
            mew=0.75,
            color="#4a4a4a",
            alpha=0.38,
        )
        for attempt in range(n):
            ax.plot(
                [x[idx] - 0.12 + jitter[attempt], x[idx] + 0.12 + jitter[attempt]],
                [free[attempt], condensate[attempt]],
                color="#b8b8b8",
                lw=0.55,
                alpha=0.45,
                zorder=0,
            )
        for values, offset, family in ((free[:n], -0.12, "free"), (condensate[:n], 0.12, "condensate")):
            mean = float(np.mean(values))
            low, high = np.percentile(values, [5, 95])
            ax.errorbar(
                x[idx] + offset,
                mean,
                yerr=[[mean - low], [high - mean]],
                fmt="o",
                ms=7,
                capsize=3,
                color="#1f1f1f",
                mec="#1f1f1f",
                mew=1.2,
                mfc="white" if family == "free" else "#1f1f1f",
                zorder=4,
            )
        difference = condensate[:n] - free[:n]
        for attempt in range(n):
            rows.append({
                "system": DISPLAY[system],
                "attempt_pair": attempt,
                "free_states": " ".join(str(state + 1) for state in free_states),
                "condensate_states": " ".join(str(state + 1) for state in condensate_states),
                "free_distance": float(free[attempt]),
                "condensate_distance": float(condensate[attempt]),
                "condensate_minus_free": float(difference[attempt]),
            })

    ax.set_xticks(x, [DISPLAY[system] for system in SYSTEMS])
    ax.set_ylabel("Weighted transition distance")
    clean(ax)
    ax.spines[["top", "right"]].set_visible(True)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="None", markersize=7.5,
                   markerfacecolor="white", markeredgecolor="#1f1f1f", markeredgewidth=1.2,
                   color="#1f1f1f", label="Free"),
            Line2D([0], [0], marker="o", linestyle="None", markersize=7.5,
                   markerfacecolor="#1f1f1f", markeredgecolor="#1f1f1f", markeredgewidth=1.2,
                   color="#1f1f1f", label="Condensate"),
        ],
        ncol=2,
        loc="upper left",
        handletextpad=0.45,
        columnspacing=1.0,
    )
    fig.tight_layout()
    write_csv(out / "source_data/FigureS9_transition_distance_model_ensemble.csv", rows)
    save(fig, out, "FigureS9_transition_distance_model_ensemble")

"""Generate the partner-turnover summary used as supplementary Figure S4."""

from .shared import *  # noqa: F403


def figure_s3_compensation(out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.8))
    rows = []
    summary = read_csv(RESULTS / "partner_loss_gain_compensation_summary.csv")
    metrics = [
        ("mean_turnover_fraction", r"$\mathcal{T}$", "Fractional turnover", (0.52, 0.76)),
        ("mean_abs_delta_fraction", r"$\mathcal{D}$", "Fractional change", (0.13, 0.212)),
        ("compensation_index", r"$\mathcal{C}$", "Compensation index", (0.695, 0.805)),
    ]
    x = np.arange(len(SYSTEMS), dtype=float)

    for label, ax, (key, symbol, ylabel, ylim) in zip(("A", "B", "C"), axes, metrics):
        for role, offset in (("Tl", -0.12), ("Th", 0.12)):
            values = [
                float(next(row for row in summary if row["system"] == system and row["role"] == role)[key])
                for system in SYSTEMS
            ]
            ax.plot(
                x + offset,
                values,
                "o",
                ms=9.0,
                color=TEMP_COLORS[role],
                label=temp_label(role),
            )
            for system, value in zip(SYSTEMS, values):
                rows.append({"system": DISPLAY[system], "role": role, "metric": key, "value": value})

        ax.set_xticks(x, [DISPLAY[system] for system in SYSTEMS])
        ax.set_ylabel(f"{ylabel}, {symbol}", fontsize=14.0, labelpad=8)
        ax.set_ylim(*ylim)
        if label == "C":
            ax.set_yticks([0.70, 0.75, 0.80], ["0.70", "0.75", "0.80"])
        ax.tick_params(axis="y", labelsize=11.5, width=1.0, length=5)
        ax.tick_params(axis="x", labelsize=10.5, width=1.0, length=5)
        clean(ax)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1.0)
        panel(ax, label)

    axes[0].legend(
        ncol=2,
        loc="upper right",
        borderaxespad=0.5,
        fontsize=11.5,
        handletextpad=0.5,
        columnspacing=1.2,
    )
    fig.tight_layout(w_pad=2.2)
    write_csv(out / "source_data/FigureS2_loss_gain_compensation.csv", rows)
    save(fig, out, "FigureS2_partner_loss_gain_compensation")

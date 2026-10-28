"""Generate supplementary Figure S4 in the same visual format as Figure S3."""

from .shared import *  # noqa: F403


ANALYSIS_DIR = (
    ROOT
    / "project_figure/supplementary/source_data/s4_contact_load_analysis"
)
SUMMARY_PATH = ANALYSIS_DIR / "FigureS4_contact_load_mechanism_summary.csv"
DISTRIBUTION_PATH = ANALYSIS_DIR / "FigureS4_contact_load_distributions.npz"
SYSTEM_MARKERS = {"FUS": "o", "A1": "s", "TDP43": "^"}
COMPONENTS = (
    ("turnover_fraction_exit_loss", "Exit loss", "#c47a3b"),
    ("turnover_fraction_enter_gain", "Entry gain", "#e1ad74"),
    ("turnover_fraction_retained_loss", "Retained loss", "#557f9b"),
    ("turnover_fraction_retained_gain", "Retained gain", "#9bb8c9"),
)
DELTA_COLORS = {
    "temperature_difference_load_state_persistence_q": "#3f6c8e",
    "temperature_difference_partner_survival": "#a45b5b",
}


def _metric_rows(
    rows: list[dict[str, str]],
    system: str,
    role: str,
    metric: str,
) -> list[dict[str, str]]:
    selected = [
        row
        for row in rows
        if row["system"] == system
        and row["role"] == role
        and row["metric"] == metric
    ]
    return sorted(selected, key=lambda row: float(row["lag_ns"]))


def _condition_curve(
    ax,
    rows: list[dict[str, str]],
    system: str,
    role: str,
    metric: str,
) -> None:
    selected = _metric_rows(rows, system, role, metric)
    x = np.asarray([float(row["lag_ns"]) for row in selected])
    y = np.asarray([float(row["value"]) for row in selected])
    low = np.asarray([float(row["ci_low"]) for row in selected])
    high = np.asarray([float(row["ci_high"]) for row in selected])
    color = TEMP_COLORS[role]
    ax.fill_between(x, low, high, color=color, alpha=0.09, lw=0)
    ax.plot(
        x,
        y,
        color=color,
        marker=SYSTEM_MARKERS[system],
        ms=6.2,
        lw=2.0,
    )


def _difference_curve(
    ax,
    rows: list[dict[str, str]],
    system: str,
    metric: str,
) -> None:
    selected = _metric_rows(rows, system, "Th-Tl", metric)
    x = np.r_[0.0, [float(row["lag_ns"]) for row in selected]]
    y = np.r_[0.0, [float(row["value"]) for row in selected]]
    low = np.r_[0.0, [float(row["ci_low"]) for row in selected]]
    high = np.r_[0.0, [float(row["ci_high"]) for row in selected]]
    color = DELTA_COLORS[metric]
    ax.fill_between(x, low, high, color=color, alpha=0.09, lw=0)
    ax.plot(
        x,
        y,
        color=color,
        marker=SYSTEM_MARKERS[system],
        ms=4.6,
        lw=1.5,
    )


def _decomposition_panel(ax, rows: list[dict[str, str]]) -> None:
    centers = np.arange(len(SYSTEMS), dtype=float)
    offsets = {"Tl": -0.18, "Th": 0.18}
    width = 0.32
    for system_index, system in enumerate(SYSTEMS):
        for role in ROLES:
            x = centers[system_index] + offsets[role]
            bottom = 0.0
            for metric, label, color in COMPONENTS:
                selected = _metric_rows(rows, system, role, metric)
                value = float(selected[0]["value"])
                ax.bar(
                    x,
                    value,
                    width=width,
                    bottom=bottom,
                    color=color,
                    edgecolor="white",
                    lw=0.45,
                    label=label if system_index == 0 and role == "Tl" else None,
                )
                ax.text(
                    x,
                    bottom + 0.5 * value,
                    f"{100 * value:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=10.0,
                    color="white" if color in ("#c47a3b", "#557f9b") else "#303030",
                )
                bottom += value
            ax.text(
                x,
                -0.052,
                temp_label(role),
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="top",
                fontsize=11.0,
            )
    ax.set_xlim(-0.56, len(SYSTEMS) - 0.44)
    ax.set_ylim(0, 1.19)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_xticks(centers, [DISPLAY[system] for system in SYSTEMS])
    ax.tick_params(axis="x", pad=20, labelsize=11.5)
    ax.tick_params(axis="y", labelsize=11.5)
    ax.set_ylabel("Fraction of contact-load turnover", fontsize=13)
    ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.99, 0.99),
        ncol=2,
        fontsize=10.0,
        columnspacing=0.8,
    )


def _restore_box(ax) -> None:
    clean(ax)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)


def figure_s4_contact_load(out: Path) -> None:
    rows = read_csv(SUMMARY_PATH)
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.2, 4.7))

    _decomposition_panel(ax_a, rows)
    _restore_box(ax_a)
    panel(ax_a, "A")

    for system in SYSTEMS:
        for role in ROLES:
            _condition_curve(
                ax_b,
                rows,
                system,
                role,
                "conditional_replacement_eps0.05",
            )
    ax_b.set_xlim(0, 20.5)
    ax_b.set_ylim(0.35, 0.71)
    ax_b.set_xlabel("Lag time (ns)", fontsize=13)
    ax_b.set_ylabel("Partner replacement fraction", fontsize=13)
    ax_b.tick_params(axis="both", labelsize=11.5)
    ax_b.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color=TEMP_COLORS["Tl"],
                lw=1.8,
                label=r"$T_{\mathrm{L}}$",
            ),
            Line2D(
                [0],
                [0],
                color=TEMP_COLORS["Th"],
                lw=1.8,
                label=r"$T_{\mathrm{H}}$",
            ),
            *[
            Line2D(
                [0],
                [0],
                marker=SYSTEM_MARKERS[system],
                color="#555555",
                lw=0,
                ms=5.3,
                label=DISPLAY[system],
            )
            for system in SYSTEMS
            ],
        ],
        loc="lower right",
        ncol=1,
        fontsize=10.2,
        handlelength=1.7,
        labelspacing=0.45,
    )
    _restore_box(ax_b)
    panel(ax_b, "B")
    fig.tight_layout(w_pad=3.2)

    source_rows = [
        {
            "system": DISPLAY[row["system"]],
            "role": row["role"],
            "lag_ns": row["lag_ns"],
            "metric": row["metric"],
            "value": row["value"],
            "ci_low": row["ci_low"],
            "ci_high": row["ci_high"],
            "n_segments": row["n_segments"],
        }
        for row in rows
        if row["metric"] in {metric for metric, _, _ in COMPONENTS}
        or row["metric"] == "conditional_replacement_eps0.05"
    ]
    write_csv(out / "source_data/FigureS3_contact_load_mechanism.csv", source_rows)
    save(fig, out, "FigureS3_contact_load_redistribution_and_state_persistence")

"""Generate supplementary Figure S8: state-resolved CK tests."""

from .shared import *  # noqa: F403


def _state_label(display_state: int) -> str:
    return rf"$s_{{{display_state + 1}}}$"


def _small_ck_axis(ax) -> None:
    ax.set_xlim(1.6, 8.4)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xticks([2, 4, 6, 8])
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.tick_params(labelsize=15.5, length=4.0, width=1.0, pad=2.0)
    ax.grid(True, color="#d9d9d9", alpha=0.28, lw=0.55)
    ax.spines[["top", "right"]].set_visible(False)
    for spine in ax.spines.values():
        spine.set_linewidth(0.85)


def figure_s8_ck(out: Path) -> None:
    source: list[dict[str, object]] = []
    panel_letters = ("A", "B", "C")
    conditions = (("free", "Tl"), ("free", "Th"), ("condensate", "Tl"), ("condensate", "Th"))

    for page_idx, system in enumerate(SYSTEMS):
        fig = plt.figure(figsize=(15.8, 11.8))
        outer = fig.add_gridspec(
            2, 2, left=0.095, right=0.985, bottom=0.075, top=0.855,
            hspace=0.300, wspace=0.145,
        )

        for j, (family, role) in enumerate(conditions):
            slot = outer[j // 2, j % 2]
            block_ax = fig.add_subplot(slot)
            block_ax.set_xticks([])
            block_ax.set_yticks([])
            for spine in block_ax.spines.values():
                spine.set_visible(False)
            phase_label = "Free" if family == "free" else "Condensate"
            block_ax.set_title(
                f"{phase_label}, {temp_label(role)}", fontsize=22.5,
                fontweight="bold", color="#111111", pad=12,
            )

            estimated, predicted = aligned_ck_stack(system, family, role, out)
            ensemble = model_ensemble(system, family, role, out)
            x = np.arange(1, estimated.shape[1] + 1) * int(ensemble["base_lag"]) * float(ensemble["dt"])
            mapping = display_state_mapping(system, family, role)
            ordered = sorted(mapping.items(), key=lambda item: item[1])
            subgrid = slot.subgridspec(len(ordered), len(ordered), hspace=0.125, wspace=0.145)

            for row_idx, (local_i, display_i) in enumerate(ordered):
                for col_idx, (local_j, display_j) in enumerate(ordered):
                    ax = fig.add_subplot(subgrid[row_idx, col_idx])
                    est = estimated[:, :, local_i, local_j]
                    pred = predicted[:, :, local_i, local_j]
                    est_mean = np.mean(est, axis=0)
                    pred_mean = np.mean(pred, axis=0)
                    pred_lo, pred_hi = np.percentile(pred, [5, 95], axis=0)
                    curve_color = TEMP_COLORS[role]
                    ax.fill_between(x, pred_lo, pred_hi, color=curve_color, alpha=0.14, linewidth=0)
                    ax.plot(x, pred_mean, color=curve_color, lw=2.1, ls="-")
                    ax.plot(x, est_mean, color=curve_color, lw=1.9, ls="--")
                    if col_idx == 0:
                        ax.set_ylabel(f"from {_state_label(display_i)}", fontsize=17.0, labelpad=2.5)
                    else:
                        ax.set_yticklabels([])
                    if row_idx != len(ordered) - 1:
                        ax.set_xticklabels([])
                    else:
                        ax.set_xlabel(f"to {_state_label(display_j)}", fontsize=16.5, labelpad=2.0)
                    _small_ck_axis(ax)
                    is_diagonal = local_i == local_j
                    for step, time, estimate, prediction, low, high in zip(range(1, len(x) + 1), x, est_mean, pred_mean, pred_lo, pred_hi):
                        source.append({
                            "system": DISPLAY[system], "family": family, "role": role,
                            "from_state": display_i + 1, "to_state": display_j + 1,
                            "ck_step": step, "time_ns": time, "estimated": estimate,
                            "predicted": prediction, "predicted_p05": low,
                            "predicted_p95": high, "self_transition": bool(is_diagonal),
                        })

        style_handles = [
            Line2D([0], [0], color="#444444", lw=2.1, label="Predicted"),
            Line2D([0], [0], color="#444444", lw=1.9, ls="--", label="Estimated"),
        ]
        fig.legend(handles=style_handles, loc="upper center", ncol=2, bbox_to_anchor=(0.54, 0.915), fontsize=17.5)
        fig.text(0.018, 0.965, panel_letters[page_idx], ha="left", va="top", fontsize=26.0, fontweight="bold")
        fig.text(0.025, 0.465, DISPLAY[system], ha="center", va="center", rotation=90, fontsize=20.0, fontweight="bold")
        save(fig, out, f"FigureS8{panel_letters[page_idx]}_state_resolved_chapman_kolmogorov_tests")

    write_csv(out / "source_data/FigureS8_state_resolved_ck.csv", source)

"""Generate supplementary Figure S5: free-chain structural transition networks."""

from .shared import *  # noqa: F403
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyArrowPatch, Rectangle


def mask_contact_map(contact_map: np.ndarray, gap: int = 4) -> np.ndarray:
    masked = np.asarray(contact_map, dtype=float).copy()
    size = masked.shape[0]
    for offset in range(gap + 1):
        index = np.arange(size - offset)
        shifted = index + offset
        masked[index, shifted] = np.nan
        masked[shifted, index] = np.nan
    return np.clip(masked, 1e-4, None)


def free_contact_map(system: str, role: str, display_state: int) -> np.ndarray:
    mapping = display_state_mapping(system, "free", role)
    local_state = next(local for local, display in mapping.items() if display == display_state)
    contacts = np.load(model_dir(system, "free", role) / "avg_contacts.npy", mmap_mode="r")
    return mask_contact_map(contacts[local_state])


def square_edge(center: np.ndarray, direction: np.ndarray, size: float, gap: float) -> np.ndarray:
    unit = direction / max(float(np.linalg.norm(direction)), 1e-12)
    distance = (0.5 * size) / max(abs(unit[0]), abs(unit[1]), 1e-12)
    return center + (distance + gap) * unit


def draw_structural_network(
    ax,
    system: str,
    role: str,
    out: Path,
    source: list[dict[str, object]],
) -> None:
    states = sorted(
        {
            int(row.get("display_state") or row.get("global_state", 0))
            for current_role in ROLES
            for row in state_rows(system, "free", current_role)
        }
    )
    matrices, populations = mapped_matrix_stack(system, "free", role, states, out)
    transition = np.mean(matrices, axis=0)
    population = np.mean(populations, axis=0)
    index = {state: idx for idx, state in enumerate(states)}

    coords = [(0.50, 0.84), (0.18, 0.26), (0.82, 0.26)]
    x_shift = 0.075 if role == "Tl" else -0.075
    positions = {
        state: np.asarray((coords[idx][0] + x_shift, coords[idx][1]), dtype=float)
        for idx, state in enumerate(states)
    }
    thumb_size = 0.355

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")

    for pair_index, (state_a, state_b) in enumerate(
        ((states[0], states[1]), (states[0], states[2]), (states[1], states[2]))
    ):
        center_a = positions[state_a]
        center_b = positions[state_b]
        direction = center_b - center_a
        unit = direction / max(float(np.linalg.norm(direction)), 1e-12)
        normal = np.asarray((-unit[1], unit[0]))
        for source_state, target_state, lane in (
            (state_a, state_b, -0.055),
            (state_b, state_a, 0.055),
        ):
            probability = float(transition[index[source_state], index[target_state]])
            if probability < 0.008:
                continue
            source_center = positions[source_state] + lane * normal
            target_center = positions[target_state] + lane * normal
            direct = target_center - source_center
            start = square_edge(source_center, direct, thumb_size, 0.034)
            end = square_edge(target_center, -direct, thumb_size, 0.034)
            label = 0.5 * (start + end) + (0.040 if lane > 0 else -0.040) * normal
            ax.add_patch(
                FancyArrowPatch(
                    start,
                    end,
                    arrowstyle="-|>",
                    mutation_scale=18,
                    linewidth=2.2,
                    color="#505050",
                    alpha=0.95,
                    transform=ax.transAxes,
                    clip_on=False,
                    zorder=2,
                )
            )
            ax.text(
                label[0],
                label[1],
                f"{probability:.2f}",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=9.1,
                color="#303030",
                zorder=3,
                clip_on=False,
            )
            source.append(
                {
                    "system": DISPLAY[system],
                    "role": role,
                    "source_state": source_state + 1,
                    "target_state": target_state + 1,
                    "transition_probability": probability,
                    "population": "",
                }
            )

    for state in states:
        x, y = positions[state]
        x0 = x - thumb_size / 2
        y0 = y - thumb_size / 2
        ax.imshow(
            free_contact_map(system, role, state),
            origin="lower",
            cmap="viridis",
            norm=LogNorm(vmin=1e-4, vmax=1e-1),
            interpolation="nearest",
            extent=(x0, x0 + thumb_size, y0, y0 + thumb_size),
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
        pop = float(population[index[state]])
        ax.text(
            x,
            y0 - 0.018,
            rf"$s_{{{state + 1}}}$, {100 * pop:.0f}%",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=11,
            fontweight="bold",
            color="#222222",
            zorder=6,
            clip_on=False,
        )
        source.append(
            {
                "system": DISPLAY[system],
                "role": role,
                "source_state": state + 1,
                "target_state": "",
                "transition_probability": "",
                "population": pop,
            }
        )


def figure_s5_free_networks(out: Path) -> None:
    fig = plt.figure(figsize=(15.4, 13.8), dpi=300)
    grid = GridSpec(3, 2, figure=fig, hspace=0.11, wspace=-0.14)
    axes = {}
    source: list[dict[str, object]] = []
    for row, system in enumerate(SYSTEMS):
        for col, role in enumerate(ROLES):
            ax = fig.add_subplot(grid[row, col])
            axes[row, col] = ax
            draw_structural_network(ax, system, role, out, source)
            if col == 0:
                ax.text(
                    0.020,
                    0.50,
                    DISPLAY[system],
                    transform=ax.transAxes,
                    ha="right",
                    va="center",
                    fontsize=17,
                    fontweight="bold",
                )
    fig.subplots_adjust(left=0.055, right=0.985, top=0.975, bottom=0.035)
    fig.canvas.draw()
    # Align each heading with the visual center of its inward-shifted motif.
    for col, role in enumerate(ROLES):
        axis_box = axes[0, col].get_position()
        motif_x = 0.575 if role == "Tl" else 0.425
        fig.text(
            axis_box.x0 + motif_x * axis_box.width,
            min(0.994, axis_box.y1 + 0.016),
            temp_label(role),
            ha="center",
            va="bottom",
            fontsize=19,
            fontweight="bold",
        )
    write_csv(out / "source_data/FigureS5_free_transition_networks.csv", source)
    save(fig, out, "FigureS5_free_state_transition_networks")

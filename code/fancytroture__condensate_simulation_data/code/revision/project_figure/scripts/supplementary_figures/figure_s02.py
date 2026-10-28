"""Generate supplementary Figure S2 directly from state contact-map arrays."""

from .shared import *  # noqa: F403


MASK_GAP = 4


def _mapped_contact_maps(
    system: str,
    family: str,
    role: str,
) -> dict[int, dict[str, object]]:
    """Return aligned contact maps and populations keyed by display state."""
    contacts = np.load(
        model_dir(system, family, role) / "avg_contacts.npy",
        mmap_mode="r",
    )

    mapped: dict[int, dict[str, object]] = {}
    for row in state_rows(system, family, role):
        local_state = int(row["local_state"])
        display_state = int(row.get("display_state") or row.get("global_state", 0))

        mapped[display_state] = {
            "contact_map": np.asarray(contacts[local_state], dtype=float),
            "population": float(row["mu"]),
        }

    return mapped


def _masked_for_log(contact_map: np.ndarray, vmin: float) -> np.ndarray:
    """Mask near-diagonal contacts and clip finite values for logarithmic display."""
    masked = np.array(contact_map, dtype=float, copy=True)

    residue = np.arange(masked.shape[0])
    masked[np.abs(residue[:, None] - residue[None, :]) <= MASK_GAP] = np.nan

    finite = np.isfinite(masked)
    masked[finite] = np.clip(masked[finite], vmin, None)

    return masked


def figure_s2_contact_maps(out: Path) -> None:
    """Build Figure S2 with family titles centered over the top visible blocks."""

    payload = {
        (system, family, role): _mapped_contact_maps(system, family, role)
        for system in SYSTEMS
        for family in FAMILIES
        for role in ROLES
    }

    values: list[np.ndarray] = []
    for states in payload.values():
        for state in states.values():
            masked = _masked_for_log(state["contact_map"], 1e-4)
            values.append(masked[np.isfinite(masked)])

    all_values = np.concatenate(values)
    vmin = max(1e-4, float(np.percentile(all_values, 2.0)))
    vmax = float(np.percentile(all_values, 99.5))

    norm = LogNorm(vmin=vmin, vmax=vmax)

    cmap = matplotlib.colormaps["viridis"].copy()
    cmap.set_bad("white")

    fig = plt.figure(figsize=(16.8, 13.2))

    outer = fig.add_gridspec(
        3,
        2,
        left=0.075,
        right=0.99,
        bottom=0.035,
        top=0.925,
        hspace=0.16,
        wspace=0.07,
        width_ratios=[3.07, 5.14],
    )

    row_centers = (0.80, 0.495, 0.19)

    # Store the top visible axes for each family. The title sits above the
    # first row, so its center should follow the visible block directly below it
    # rather than hidden/empty state slots in lower rows.
    family_title_axes = {family: [] for family in FAMILIES}

    for row_index, system in enumerate(SYSTEMS):
        fig.text(
            0.022,
            row_centers[row_index],
            DISPLAY[system],
            rotation=90,
            ha="center",
            va="center",
            fontsize=20,
            fontweight="bold",
        )

        for family_index, family in enumerate(FAMILIES):
            state_slots = 3 if family == "free" else 5

            states = sorted(
                set(payload[(system, family, "Tl")])
                | set(payload[(system, family, "Th")])
            )

            inner = outer[row_index, family_index].subgridspec(
                2,
                state_slots,
                hspace=0.10,
                wspace=0.035,
            )

            for role_index, role in enumerate(ROLES):
                maps = payload[(system, family, role)]

                for state_index in range(state_slots):
                    ax = fig.add_subplot(inner[role_index, state_index])
                    empty_state = False

                    if row_index == 0:
                        family_title_axes[family].append(ax)

                    if state_index >= len(states):
                        ax.set_axis_off()
                        continue

                    state = states[state_index]

                    if state in maps:
                        state_data = maps[state]
                        contact_map = state_data["contact_map"]

                        ax.imshow(
                            _masked_for_log(contact_map, vmin),
                            origin="lower",
                            cmap=cmap,
                            norm=norm,
                            interpolation="nearest",
                            rasterized=True,
                        )

                        ax.text(
                            0.965,
                            0.035,
                            f"{100.0 * state_data['population']:.1f}%",
                            transform=ax.transAxes,
                            ha="right",
                            va="bottom",
                            fontsize=13.5,
                            fontweight="bold",
                            color="white",
                        )

                    else:
                        empty_state = True

                    if empty_state:
                        ax.set_axis_off()
                        continue

                    reference_maps = maps
                    if not reference_maps:
                        other_role = "Th" if role == "Tl" else "Tl"
                        reference_maps = payload[(system, family, other_role)]

                    n_res = next(iter(reference_maps.values()))["contact_map"].shape[0]
                    ticks = [0, n_res // 2, n_res - 1]

                    ticklabels = [0, 74, 147] if system == "TDP43" else ticks

                    ax.set_xticks(ticks)
                    ax.set_yticks(ticks)
                    ax.set_xticklabels(ticklabels)
                    ax.set_yticklabels(ticklabels)
                    ax.tick_params(
                        axis="both",
                        labelsize=11,
                        length=2.8,
                        pad=1,
                    )

                    if role_index == 0 or state not in payload[(system, family, "Tl")]:
                        ax.set_title(
                            rf"$s_{{{state + 1}}}$",
                            fontsize=19,
                            fontweight="bold",
                            pad=7,
                        )
                        ax.tick_params(labelbottom=False)

                    if role_index == 1 and state_index != 0:
                        ax.tick_params(labelbottom=False)

                    if state_index == 0 and family == "free":
                        temperature_label = (
                            r"$\mathbfit{T}_{\mathrm{L}}$"
                            if role == "Tl"
                            else r"$\mathbfit{T}_{\mathrm{H}}$"
                        )

                        ax.set_ylabel(
                            temperature_label,
                            fontsize=21,
                            rotation=0,
                            labelpad=19,
                            va="center",
                            color="#111111",
                            fontweight="bold",
                        )
                    else:
                        ax.tick_params(labelleft=False)

                    for spine in ax.spines.values():
                        spine.set_linewidth(0.65)

    # Resolve the final subplot positions before placing family titles.
    fig.canvas.draw()

    family_titles = {
        "free": "Free",
        "condensate": "Condensate",
    }

    for family in FAMILIES:
        boxes = [ax.get_position() for ax in family_title_axes[family]]
        family_left = min(box.x0 for box in boxes)
        family_right = max(box.x1 for box in boxes)
        family_center = 0.5 * (family_left + family_right)

        fig.text(
            family_center,
            0.982,
            family_titles[family],
            ha="center",
            va="top",
            fontsize=22,
            fontweight="bold",
        )

    save(
        fig,
        out,
        "FigureS4_full_state_resolved_intrachain_contact_maps",
    )

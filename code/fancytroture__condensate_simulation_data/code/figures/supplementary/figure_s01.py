"""Generate supplementary Figure S1."""

from .shared import *  # noqa: F403


def figure_s1_density(out: Path) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(8.2, 9.4), sharex=True)
    source = []
    for row, system in enumerate(SYSTEMS):
        ax_profile = axes[row]
        for role in ROLES:
            overall_path, block_path = density_paths(system, role)
            overall = np.loadtxt(overall_path)
            z, rho = overall[:, 0], overall[:, 1]
            block_data = np.loadtxt(block_path)
            block_z = block_data[:, 0]
            # Block-density files are stored in mg/mL, whereas the overall
            # profile is stored in g/mL.
            block_profiles = block_data[:, 1:] / 1000.0
            if not np.allclose(z, block_z):
                raise ValueError(f"Density grids do not match for {system} {role}")
            block_mean = np.mean(block_profiles, axis=1)
            block_sd = np.std(block_profiles, axis=1, ddof=1)
            ax_profile.fill_between(
                z,
                block_mean - block_sd,
                block_mean + block_sd,
                color=TEMP_COLORS[role],
                alpha=0.22,
                linewidth=0,
            )
            ax_profile.plot(z, rho, color=TEMP_COLORS[role], lw=2.8, label=temp_label(role))
            for z_value, rho_value, mean_value, sd_value in zip(z, rho, block_mean, block_sd):
                source.append(
                    {
                        "system": DISPLAY[system],
                        "role": role,
                        "temperature": TEMPS[system][role],
                        "z_nm": z_value,
                        "density_g_ml": rho_value,
                        "block_mean_density_g_ml": mean_value,
                        "block_sd_density_g_ml": sd_value,
                        "n_time_blocks": block_profiles.shape[1],
                    }
                )
        ax_profile.text(
            0.025,
            0.94,
            DISPLAY[system],
            transform=ax_profile.transAxes,
            fontsize=16,
            fontweight="bold",
            ha="left",
            va="top",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.5},
        )
        ax_profile.set_ylabel(r"$\rho(z)$ (g/mL)", fontsize=14)
        if row == len(SYSTEMS) - 1:
            ax_profile.set_xlabel("z (nm)", fontsize=15)
        ax_profile.set_xlim(-40, 40)
        ax_profile.tick_params(axis="both", labelsize=12.5, width=1.0, length=5)
        clean(ax_profile)
        ax_profile.spines[["top", "right"]].set_visible(True)
        for spine in ax_profile.spines.values():
            spine.set_linewidth(1.05)
    axes[0].legend(
        ncol=2,
        loc="upper right",
        fontsize=14,
        handlelength=2.4,
        columnspacing=1.6,
        borderaxespad=0.7,
    )
    fig.tight_layout(h_pad=1.0)
    write_csv(out / "source_data/FigureS1_density_profiles.csv", source)
    save(fig, out, "FigureS1_phase_separation_and_density_stability")

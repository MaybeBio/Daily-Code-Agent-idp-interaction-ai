#!/usr/bin/env python3
"""Supplementary Figure 3: |dG_ii - dG_jj| for mixing vs demixing systems,
split by interaction type.

Mirrors the overlaid-histogram style of Figure3.py's Panel D, but broken out
into one panel per chemical interaction category (Charged/Polar/Hydrophobic
pairings) instead of pooling all pairs together. As of this version, the
mix/demix split uses the same compound criteria as Figure3.py's Panel D:
    Mix   = D_norm <= 0.13 AND dGij < -3
    Demix = D_norm >  0.13 AND (dG_ii < -3 OR dG_jj < -3)
Pairs matching neither condition are excluded from both histograms (but kept
in the exported CSV, labeled 'Excluded'), same as Panel D.

Data come from the rebuilt master CSV plus gin_samples.csv, merged the same
way Figure3.py does it, so the GIN categories and homotypic dG values are
taken from one source rather than being duplicated per-script.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from matplotlib.ticker import MaxNLocator
from pathlib import Path

plt.rcParams.update({
    'font.size': 24, 'axes.linewidth': 2, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.size': 8, 'xtick.major.width': 2, 'ytick.major.size': 8, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300
})

DATA_DIR = Path("supplementary3_data")
HET_CSV = "../DDM_Heterotypic/heterotypic_main.csv"
GIN_CSV = "../DATASETS/gin_samples.csv"

MIN_DENSE_FRAC = 0.05  # two phases only if >= 5% of voxels sit at or above c_den
D_CUTOFF = 0.13        # boundary between mix (<=) and demix (>), matches Figure3.py Panel D
DG_CUTOFF = -3         # dG threshold for both legs of the mix/demix criteria, matches Panel D
COMP_MAX = 95          # drop near-pure compositions

# GIN cluster -> general category (Supplementary Table S1)
GIN_CATEGORY = {
    **{k: 'Charged' for k in [3, 7, 8, 9, 17, 18, 19, 23, 24, 25, 26, 29]},
    **{k: 'Polar' for k in [0, 1, 2, 4, 5, 6, 11, 12, 13, 14, 15, 16, 21, 22, 27]},
    **{k: 'Hydrophobic' for k in [10, 20, 28]},
}

INTERACTIONS = ['Charged-Charged', 'Charged-Hydrophobic', 'Charged-Polar',
                'Hydrophobic-Hydrophobic', 'Hydrophobic-Polar', 'Polar-Polar']

# Same base palette as Supplementary_Figure4.py (sns husl, one hue per panel),
# so a given interaction type reads as the same color across both figures.
# Within each panel here, a lighter tint marks Mix and a darker shade marks Demix.
PANEL_BASE_COLOR = dict(zip(INTERACTIONS, sns.color_palette("husl", 6)))
MIX_TINT = 0.45     # blend toward white, for Mix
DEMIX_SHADE = 0.30  # blend toward black, for Demix


def tint(color, frac):
    """Blend a color toward white by frac (0=unchanged, 1=white)."""
    r, g, b = mcolors.to_rgb(color)
    return (r + (1 - r) * frac, g + (1 - g) * frac, b + (1 - b) * frac)


def shade(color, frac):
    """Blend a color toward black by frac (0=unchanged, 1=black)."""
    r, g, b = mcolors.to_rgb(color)
    return (r * (1 - frac), g * (1 - frac), b * (1 - frac))


PANEL_COLORS = {interaction: (tint(base, MIX_TINT), shade(base, DEMIX_SHADE))
                for interaction, base in PANEL_BASE_COLOR.items()}


def set_mirrored_yaxis(ax, pos_max, neg_max, target_ticks=5):
    """Give a bidirectional bar chart (positive values above 0, negative below)
    a consistent tick step on both sides, then label every tick with its
    magnitude only (abs value) so it reads as a single shared 'Count' axis.

    Relying on the default locator alone can leave one side with no non-zero
    tick at all when its range is much smaller than the other's -- e.g. demix
    bars topping out at ~15 while mix bars reach -60 with the default step of
    20, which skips right over the demix side and leaves it unlabeled. This
    picks one 'nice' step from the larger of the two ranges, then lays ticks
    out independently on each side using that step, so every panel gets at
    least one labeled tick near each side's own max.
    """
    overall_max = max(pos_max, neg_max, 1)
    nice_ticks = MaxNLocator(nbins=target_ticks).tick_values(0, overall_max)
    step = nice_ticks[1] - nice_ticks[0] if len(nice_ticks) > 1 else max(overall_max, 1)

    pos_ticks = np.arange(step, pos_max + step, step) if pos_max > 0 else np.array([])
    neg_ticks = -np.arange(step, neg_max + step, step) if neg_max > 0 else np.array([])
    yticks = np.concatenate([neg_ticks[::-1], [0], pos_ticks])

    ax.set_yticks(yticks)
    ax.set_yticklabels([f'{abs(int(t))}' for t in yticks])
    ax.set_ylim(neg_ticks.min() - step * 0.3 if len(neg_ticks) else -step * 0.3,
                pos_ticks.max() + step * 0.3 if len(pos_ticks) else step * 0.3)


def save_csv(df, filename):
    DATA_DIR.mkdir(exist_ok=True)
    df.to_csv(DATA_DIR / filename, index=False)


def apply_two_phase_gate(df):
    """A voxel map resolves two phases only when at least MIN_DENSE_FRAC of
    its voxels sit at or above c_den. Otherwise c_dil is the mean over all
    voxels and no dense value (and so no dGij) is defined."""
    df = df.copy()
    df['two_phase'] = df['frac_vox_ge_cden'] >= MIN_DENSE_FRAC
    one = ~df['two_phase']
    for dil in ('cdil', 'cdil1', 'cdil2'):
        df.loc[one, dil] = df.loc[one, f'{dil}_onephase']
    df.loc[one, ['cden', 'cden1', 'cden2', 'dGij', 'dGij_err']] = np.nan
    return df


def get_interaction_type(gin1, gin2):
    if pd.isna(gin1) or pd.isna(gin2):
        return None
    cat1, cat2 = GIN_CATEGORY.get(int(gin1)), GIN_CATEGORY.get(int(gin2))
    if not cat1 or not cat2:
        return None
    return '-'.join(sorted([cat1, cat2]))


def load_data():
    """Heterotypic results with per-species GIN group and homotypic dG attached."""
    df = apply_two_phase_gate(pd.read_csv(HET_CSV))
    df = df[df['two_phase']]          # gated-out (single-phase) maps are not analysed
    df_prot = pd.read_csv(GIN_CSV)

    for name_col, suf in [('seq_name1', '_p1'), ('seq_name2', '_p2')]:
        df = pd.merge(df, df_prot, left_on=name_col, right_on='seq_name',
                      how='left', suffixes=('', suf))
        df = df.drop(columns=['seq_name'])
        df = df.rename(columns={c: c + suf for c in df_prot.columns
                                if c != 'seq_name' and c in df.columns
                                and not c.endswith(('_p1', '_p2'))})

    # Keep both mixing and demixing systems here; the mix/demix split happens
    # per interaction-type panel below, same as Figure3.py Panel D.
    df = df[(df['comp1'] < COMP_MAX) & (df['comp2'] < COMP_MAX)].copy()
    df['dG_diff'] = (df['dG_p1'] - df['dG_p2']).abs()
    df = df.dropna(subset=['dGij', 'dG_diff', 'D_norm', 'dG_p1', 'dG_p2'])
    df['interaction'] = df.apply(
        lambda r: get_interaction_type(r['GIN_group_p1'], r['GIN_group_p2']), axis=1)
    return df[df['interaction'].notna()].copy()


def main():
    df = load_data()
    print(f"{len(df)} rows across {df['interaction'].nunique()} interaction types "
          f"(mix: D_norm <= {D_CUTOFF} & dGij < {DG_CUTOFF}; "
          f"demix: D_norm > {D_CUTOFF} & (dG_ii < {DG_CUTOFF} or dG_jj < {DG_CUTOFF}))")

    fig, axes = plt.subplots(2, 3, figsize=(24, 16))
    bins = np.arange(0, 10.5, 0.5)
    centers = (bins[:-1] + bins[1:]) / 2

    demix_group_all = pd.Series('Excluded', index=df.index)

    for idx, (interaction, ax) in enumerate(zip(INTERACTIONS, axes.flatten())):
        ax.text(-0.1, 1.12, chr(65 + idx), transform=ax.transAxes,
                fontsize=28, fontweight='bold', va='top', ha='right')
        df_int = df[df['interaction'] == interaction]

        mix_mask = (df_int['D_norm'] <= D_CUTOFF) & (df_int['dGij'] < DG_CUTOFF)
        demix_mask = (df_int['D_norm'] > D_CUTOFF) & ((df_int['dG_p1'] < DG_CUTOFF) | (df_int['dG_p2'] < DG_CUTOFF))
        mix_df = df_int[mix_mask]
        demix_df = df_int[demix_mask]

        demix_group_all.loc[mix_df.index] = 'Mix'
        demix_group_all.loc[demix_df.index] = 'Demix'

        if len(mix_df) + len(demix_df) < 3:
            ax.text(0.5, 0.5, f'{interaction}\n(n={len(mix_df) + len(demix_df)}, insufficient data)',
                    ha='center', va='center', transform=ax.transAxes)
        else:
            mix_color, demix_color = PANEL_COLORS[interaction]

            demix_counts, _ = np.histogram(demix_df['dG_diff'], bins=bins)
            mix_counts, _ = np.histogram(mix_df['dG_diff'], bins=bins)

            ax.bar(centers, demix_counts, width=0.45, color=demix_color, edgecolor='black',
                   label=f'Demix, n={len(demix_df)}')
            ax.bar(centers, -mix_counts, width=0.45, color=mix_color, edgecolor='black',
                   label=f'Mix, n={len(mix_df)}')
            ax.axhline(0, color='black', linewidth=1)
            set_mirrored_yaxis(ax, demix_counts.max(), mix_counts.max())
            ax.legend(loc='lower right', framealpha=0.95, edgecolor='black')

        ax.set_xlabel('|$\\Delta$G$_{ii}$ - $\\Delta$G$_{jj}$| ($k_BT$)')
        ax.set_ylabel('Count')
        ax.set_title(interaction, fontweight='bold', fontsize=24)
        ax.set_xlim(-0.1, 10)
        ax.set_box_aspect(1)

    plt.tight_layout()
    plt.savefig('Supplementary_Figure3.pdf', bbox_inches='tight')
    save_csv(df.assign(demix_group=demix_group_all)[
                 ['seq_name1', 'seq_name2', 'comp1', 'comp2', 'interaction',
                  'dG_p1', 'dG_p2', 'dG_diff', 'dGij', 'D_norm', 'demix_group']],
             'supplementary3_data.csv')
    print(f"Saved Supplementary_Figure3.pdf, data in {DATA_DIR}/")


if __name__ == '__main__':
    main()
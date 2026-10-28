#!/usr/bin/env python3
"""Supplementary Figure 4: Heterotypic dGij vs the average homotypic dG,
split by interaction type.

Data come from the rebuilt master CSV plus gin_samples.csv, merged the same
way Figure3.py does it, so the GIN categories and homotypic dG values are
taken from one source rather than being duplicated per-script.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import linregress

plt.rcParams.update({
    'font.size': 24, 'axes.linewidth': 2, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.size': 8, 'xtick.major.width': 2, 'ytick.major.size': 8, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300
})

DATA_DIR = Path("supplementary4_data")
HET_CSV = "../DDM_Heterotypic/heterotypic_main.csv"
GIN_CSV = "../DATASETS/gin_samples.csv"

MIN_DENSE_FRAC = 0.05 # two phases only if >= 5% of voxels sit at or above c_den
D_CUTOFF = 0.13       # well-mixed only
COMP_MAX = 95         # drop near-pure compositions

# GIN cluster -> general category (Supplementary Table S1)
GIN_CATEGORY = {
    **{k: 'Charged' for k in [3, 7, 8, 9, 17, 18, 19, 23, 24, 25, 26, 29]},
    **{k: 'Polar' for k in [0, 1, 2, 4, 5, 6, 11, 12, 13, 14, 15, 16, 21, 22, 27]},
    **{k: 'Hydrophobic' for k in [10, 20, 28]},
}

INTERACTIONS = ['Charged-Charged', 'Charged-Hydrophobic', 'Charged-Polar',
                'Hydrophobic-Hydrophobic', 'Hydrophobic-Polar', 'Polar-Polar']


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

    df = df[(df['D_norm'] < D_CUTOFF)
            & (df['comp1'] < COMP_MAX) & (df['comp2'] < COMP_MAX)].copy()
    df['avg_dG'] = (df['dG_p1'] + df['dG_p2']) / 2
    df = df.dropna(subset=['dGij', 'avg_dG'])
    df['interaction'] = df.apply(
        lambda r: get_interaction_type(r['GIN_group_p1'], r['GIN_group_p2']), axis=1)
    return df[df['interaction'].notna()].copy()


def main():
    df = load_data()
    print(f"{len(df)} well-mixed rows (D_norm < {D_CUTOFF}) across "
          f"{df['interaction'].nunique()} interaction types")

    fig, axes = plt.subplots(2, 3, figsize=(24, 16))
    colors = sns.color_palette("husl", 6)
    fits = []

    for idx, (interaction, ax) in enumerate(zip(INTERACTIONS, axes.flatten())):
        ax.text(-0.1, 1.12, chr(65 + idx), transform=ax.transAxes,
                fontsize=28, fontweight='bold', va='top', ha='right')
        df_int = df[df['interaction'] == interaction]

        if len(df_int) < 3:
            ax.text(0.5, 0.5, f'{interaction}\n(n={len(df_int)}, insufficient data)',
                    ha='center', va='center', transform=ax.transAxes)
        else:
            sns.histplot(data=df_int, x='avg_dG', y='dGij', binwidth=(0.3, 0.3),
                         cbar=True, ax=ax, cbar_kws={'shrink': 0.7, 'pad': 0.02, 'format': '%d'},
                         color=colors[idx])
            cbar_ax = ax.figure.axes[-1]
            cbar_ax.yaxis.set_label_position('left')
            cbar_ax.set_ylabel('Count', rotation=90, labelpad=8)

            slope, intercept, r, *_ = linregress(df_int['avg_dG'], df_int['dGij'])
            x_fit = np.array([-10, 0])
            ax.plot(x_fit, slope * x_fit + intercept, 'k--', lw=3, zorder=10)
            ax.plot([-10, 0], [-10, 0], color='grey', ls='--', lw=2, alpha=0.5)
            ax.text(0.05, 0.95,
                    f'Slope = {slope:.2f}\n$R^2$ = {r**2:.2f}\nn = {len(df_int)}',
                    transform=ax.transAxes, fontsize=22, va='top',
                    bbox=dict(boxstyle='round', fc='white', alpha=0.8, ec='black'))
            fits.append({'interaction': interaction, 'slope': slope,
                         'intercept': intercept, 'r2': r ** 2, 'n': len(df_int)})

        ax.set_xlabel('($\\Delta$G$_{ii}$ + $\\Delta$G$_{jj}$)/2 ($k_BT$)')
        ax.set_ylabel('Heterotypic $\\Delta$G$_{ij}$ ($k_BT$)')
        ax.set_title(interaction, fontsize=24, fontweight='bold')
        ax.set_xlim(0.5, -10)
        ax.set_ylim(0.5, -10)
        ax.set_box_aspect(1)

    plt.tight_layout()
    plt.savefig('Supplementary_Figure4.pdf', bbox_inches='tight')

    save_csv(df[['seq_name1', 'seq_name2', 'comp1', 'comp2', 'interaction',
                 'dG_p1', 'dG_p2', 'avg_dG', 'dGij', 'D_norm']],
             'supplementary4_data.csv')
    if fits:
        save_csv(pd.DataFrame(fits), 'supplementary4_fits.csv')
        print(pd.DataFrame(fits).to_string(index=False))
    print(f"Saved Supplementary_Figure4.pdf, data in {DATA_DIR}/")


if __name__ == '__main__':
    main()
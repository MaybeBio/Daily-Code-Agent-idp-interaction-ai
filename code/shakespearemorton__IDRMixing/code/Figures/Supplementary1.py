#!/usr/bin/env python3
"""Supplementary Figure: Violin comparison"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import seaborn as sns

plt.rcParams.update({
    'font.size': 20, 'axes.linewidth': 2, 'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.size': 8, 'xtick.major.width': 2, 'ytick.major.size': 8, 'ytick.major.width': 2,
    'legend.frameon': False, 'figure.dpi': 300
})

# Retaining only the colors used for the violin plots
COLORS = {'ds1': '#ffffff', 'ds2': '#969696'}

def create_violin_figure():
    df_training = pd.read_csv('../DATASETS/df_training.csv')
    gin_samples = pd.read_csv('../DATASETS/gin_samples.csv')
    
    # Create figure scaled for a single row of 6 plots
    fig = plt.figure(figsize=(30, 7))
    gs_violin = gridspec.GridSpec(1, 6, figure=fig, wspace=0.5)
    
    params = ['N', 'mean_lambda', 'dG', 'faro', 'fcr', 'ncpr']
    ylabels = {'N': 'N', 'mean_lambda': r'$\langle\lambda\rangle$', 'dG': r'$\Delta G$',
               'faro': r'$f_{\mathrm{aro}}$', 'fcr': 'FCR', 'ncpr': 'NCPR'}
    titles = {'N': 'Sequence Length', 'mean_lambda': 'Mean Residue\nHydrophobicity', 'dG': 'Transfer Free Energy',
              'faro': 'Fraction Aromatic\nResidues', 'fcr': 'Fraction Charged\nResidues', 'ncpr': 'Net Charge Per\nResidue'}
    palette = {'Our Dataset': COLORS['ds1'], 'von Bülow Dataset': COLORS['ds2']}
    
    # Labels A through F for the new individual plots
    panel_labels = ['A', 'B', 'C', 'D', 'E', 'F']
    
    for idx, param in enumerate(params):
        ax_v = fig.add_subplot(gs_violin[0, idx])
        
        # Add the A-F panel label
        ax_v.text(-0.1, 1.1, panel_labels[idx], transform=ax_v.transAxes, fontsize=24, fontweight='bold', va='top', ha='right')
        
        data = pd.concat([pd.DataFrame({'value': gin_samples[param].values, 'source': 'Our Dataset'}),
                          pd.DataFrame({'value': df_training[param].values, 'source': 'von Bülow Dataset'})], ignore_index=True)
        data = data.dropna(subset=['value'])
        
        sns.violinplot(x='source', y='value', data=data, ax=ax_v, palette=palette, inner='quartile',
                       hue='source', split=True, legend=False, density_norm='width')
        sns.stripplot(x='source', y='value', data=data, ax=ax_v, color='black', size=8, alpha=0.5, jitter=True)
        
        ax_v.set_ylabel(ylabels[param], fontsize=18)
        ax_v.set_xlabel(titles[param], fontsize=18)
        ax_v.set_xticks([])
        
        ax_v.spines['top'].set_visible(False)
        ax_v.spines['right'].set_visible(False)
    
    # Add legend above the figure
    fig.legend(handles=[mpatches.Patch(fc=c, ec='black', lw=1.5, label=l) for l, c in palette.items()],
               loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=2, fontsize=18)
    
    plt.tight_layout()
    plt.savefig('Supplementary_Figure_Violins.pdf', bbox_inches='tight')
    print("Supplementary Figure saved")

if __name__ == '__main__':
    create_violin_figure()
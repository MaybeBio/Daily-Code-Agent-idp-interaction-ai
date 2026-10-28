import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import MDAnalysis as mda
from MDAnalysis.analysis import contacts, distances
import io
from PIL import Image
from collections import defaultdict
import itertools
import pickle
import seaborn as sns
import glob
import matplotlib.ticker as ticker
import os

# High-quality figure settings — Science Advances style
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
rcParams['font.size'] = 8
rcParams['axes.labelsize'] = 8
rcParams['axes.titlesize'] = 8
rcParams['xtick.labelsize'] = 7
rcParams['ytick.labelsize'] = 7
rcParams['figure.dpi'] = 600
rcParams['savefig.dpi'] = 600
rcParams['figure.facecolor'] = 'white'
rcParams['axes.facecolor'] = 'white'

fig_dir = "/Users/adelielouet/Documents/science/Reports/MSM_paper/mv_copies_march_8/figures_hi_res"
os.makedirs(fig_dir, exist_ok=True)


xtc = '/Users/adelielouet/Documents/science/AB_G5_original_simu_analysis/trajectories/Gabis_paper/traj_all-skip-0-noW_G5.xtc'
pdb = '/Users/adelielouet/Documents/science/AB_G5_original_simu_analysis/trajectories/Gabis_paper/template_G5.pdb'

abeta="A\u03B2-42"
alpha_syn_c="\u03B1-syn-C-term"
alpha_full_length="\u03B1-syn-full"

protein_name=abeta

def contacts_within_cutoff(u, group_a, group_b, radius=3.0):
    timeseries = []
    for ts in u.trajectory:
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([n_contacts])
    return np.array(timeseries)



u = mda.Universe(pdb, xtc)

resid_list = ['resid ' + str(i) for i in range(1, len(u.residues))]

ca_df = pd.DataFrame(index=range(len(u.trajectory)))

for y in resid_list:
    #ligand = u.select_atoms('resname liga')
    ligand=u.select_atoms('resid 43')
    pocket = u.select_atoms(y)
    ca = contacts_within_cutoff(u, ligand, pocket, radius=3.0)
    ca_df[y] = ca.flatten()

ca_df.rename(columns={x: y for x, y in zip(ca_df.columns, range(0, len(ca_df.columns)))}, inplace=True)

new_column_names = [f"Aromatic Residue {i}" for i in range(0, 5)]
if len(new_column_names) == len(ca_df.columns):
    ca_df.columns = new_column_names

ca_df = ca_df.applymap(lambda x: 1 if x > 0 else x)
ca_df = ca_df.replace(0, np.nan)

ca_df.rename(columns={x:y for x,y in zip(ca_df.columns,range(0,len(ca_df.columns)))})


new_column_names = [f"Residue {i}" for i in range(1, 43)]
if len(new_column_names) == len(ca_df.columns):
    ca_df.columns = new_column_names

ca_df = ca_df.applymap(lambda x: 1 if x > 0 else x)
ca_df = ca_df.replace(0, np.nan)

palette = sns.dark_palette("#69d", n_colors=len(ca_df.columns))


fig, ax = plt.subplots(figsize=(6, 3.5))

save_feat = []
for i, (col, color) in enumerate(zip(ca_df.columns, palette)):
    y = [i if not pd.isnull(ca_df[col][j]) else None for j in range(len(ca_df))]
    x = list(range(len(ca_df)))
    save_feat.append([x, y])
    ax.scatter(x, y, s=1.3, color=color, alpha=0.8, linewidths=0)

ax.set_yticks(range(0, len(ca_df.columns), 5))
ax.set_yticklabels(range(0, len(ca_df.columns), 5), fontsize=7)
ax.tick_params(axis='x', labelsize=7)
ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
ax.xaxis.get_offset_text().set_fontsize(7)
import matplotlib.ticker as mticker
ax.xaxis.set_major_formatter(mticker.FuncFormatter(
    lambda x, _: r'$\times10^{' + f'{int(np.log10(max(abs(x),1))):d}' + r'}$'
    if x != 0 else '0'
))
ax.xaxis.set_major_formatter(mticker.FuncFormatter(
    lambda val, _: f'{val/1e5:.1f}' if val != 0 else '0'
))
ax.set_xlabel(r'Frames ($\times10^5$)', fontsize=8)
ax.set_xlabel('Frames', fontsize=8)
ax.set_ylabel('Side Chains', fontsize=8)
ax.text(-0.08, 1.02, 'a', transform=ax.transAxes, fontsize=10, fontweight='bold', va='top')
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
plt.tight_layout()
plt.xlim(100000, 200000)  # Zoom into the desired time range

plt.savefig(f"{fig_dir}/handoff_scatter_zoom.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/handoff_scatter_zoom.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()




########## SNS Heatmap #######################
#1. Case A - handoff between contact[i] --> contact[j], irrespective of whether or not contact[i] loses its contact
ca_df = ca_df.replace(np.nan,0)
ca_df_mat=ca_df.to_numpy()

ligand = u.select_atoms('not protein')
protein = u.select_atoms('protein')

sequence=ca_df_mat
residues = ((list(atom.resname for atom in protein.residues)))
num_residues = len(residues)

transition_matrix = np.zeros((num_residues, num_residues))

for step in range(len(sequence) - 1):
    current_contact = sequence[step]
    next_contact = sequence[step + 1]

    for i, residue_i in enumerate(residues):
        for j, residue_j in enumerate(residues):
            if current_contact[i] == 1 and next_contact[j] == 1:
                transition_matrix[i][j] += 1


#2. Case B -  handoff between contact[i] --> contact[j] with contact[i] losing its contact
transition_matrix_lose = np.zeros((num_residues, num_residues))

for step in range(len(sequence) - 1):
    current_contact = sequence[step]
    next_contact = sequence[step + 1]

    for i, residue_i in enumerate(residues):
        for j, residue_j in enumerate(residues):
            if current_contact[i] == 1 and next_contact[j] == 1 and next_contact[i] == 0:
                transition_matrix_lose[i][j] += 1


# Case C - handoff between contact[i] --> contact[j] with contact[i] keeping its contact

transition_matrix_keep = np.zeros((num_residues, num_residues))

for step in range(len(sequence) - 1):
    current_contact = sequence[step]
    next_contact = sequence[step + 1]

    for i, residue_i in enumerate(residues):
        for j, residue_j in enumerate(residues):
            if current_contact[i] == 1 and next_contact[j] == 1 and next_contact[i] == 1:
                transition_matrix_keep[i][j] += 1


matrices = [
    transition_matrix,
    transition_matrix_lose,
    transition_matrix_keep
]

titles = [
    r"Residue $\mathit{i}$ loses / keeps contact",
    r"Residue $\mathit{i}$ loses contact",
    r"Residue $\mathit{i}$ keeps contact"
]

fig, axes = plt.subplots(
    1, 3,
    figsize=(7, 2.5),
    constrained_layout=True,
    sharex=True,
    sharey=True
)

for idx, (ax, mat, title) in enumerate(zip(axes, matrices, titles)):
    hm = sns.heatmap(
        mat,
        ax=ax,
        cmap="rocket_r",
        square=True,
        cbar=False
    )

    ax.set_title(title, fontsize=8, pad=6)
    ax.set_xlabel(r"To Residues $\mathit{j}$", fontsize=8)

    if idx == 0:
        ax.set_ylabel(r"From Residues $\mathit{i}$", fontsize=8, rotation=90)
    else:
        ax.set_ylabel("")

    ax.tick_params(axis='y', labelleft=True, labelsize=6, rotation=0)
    ax.tick_params(axis='x', labelsize=6)
    ax.set_aspect("equal")

cbar = fig.colorbar(
    hm.collections[0],
    ax=axes,
    location="right",
    shrink=0.9,
    pad=0.02
)
cbar.ax.tick_params(labelsize=6)

plt.savefig(f"{fig_dir}/handoff_heatmap.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/handoff_heatmap.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()
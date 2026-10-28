from collections import Counter, defaultdict, OrderedDict
import itertools
import math
import random
from itertools import combinations, islice
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
import scipy
from scipy.linalg import expm, eig, norm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
import community as c
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
import mdtraj as md
import pickle
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

protein_name='abeta_gabis'

fig_dir = "/Users/adelielouet/Documents/science/Reports/MSM_paper/mv_copies_march_8/figures_hi_res"
os.makedirs(fig_dir, exist_ok=True)

file = open(f"/Users/adelielouet/Documents/science/dd_proj/msm_full_model_final_draft_march_9/pickled_files/{protein_name}/kd_centrality_ordered.pckl","rb")
kd_centrality_ordered = pickle.load(file)
file.close()

file = open(f"/Users/adelielouet/Documents/science/dd_proj/msm_full_model_final_draft_march_9/pickled_files/{protein_name}/inv_map.pckl","rb")
inv_map = pickle.load(file)
file.close()

file = open(f"/Users/adelielouet/Documents/science/dd_proj/msm_full_model_final_draft_march_9/pickled_files/{protein_name}/dictionary_transitions.pckl", 'rb')
dictionary_transitions_sorted = pickle.load(file)
file.close()

### Calculating Linear Motifs
evaluate_number_of_clusters = 10
fig, axs = plt.subplots(evaluate_number_of_clusters, 1, figsize=(3.5, 8), sharex=True)

all_cluster_elements = []
axs = axs.ravel()

for i, x in enumerate(list(kd_centrality_ordered.keys())[:evaluate_number_of_clusters]):
    cluster = [elem for xs in inv_map[x] for elem in xs]
    frequency = Counter(cluster)
    all_cluster_elements.extend(cluster)

    x_vals = sorted(frequency.keys())
    y_vals = [frequency[num] for num in x_vals]

    axs[i].bar(x_vals, y_vals, color='#4C72B0', width=0.8)
    axs[i].set_ylabel('Frequency', fontsize=7)
    axs[i].set_title(f'Cluster {x}', fontsize=7, pad=2)
    axs[i].tick_params(labelsize=6)
    sns.despine(ax=axs[i])

plt.xlabel('Residue', fontsize=8)
plt.tight_layout()
plt.savefig(f"{fig_dir}/cluster_residue_frequency.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/cluster_residue_frequency.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()


### Accumulated frequency across top clusters
frequency_all_cluster = Counter(all_cluster_elements)
x_vals_all = sorted(frequency_all_cluster.keys())
y_vals_all = [frequency_all_cluster[num] for num in x_vals_all]

fig, ax = plt.subplots(figsize=(3.5, 3.5))
ax.bar(x_vals_all, y_vals_all, color='#4C72B0', width=0.8)
ax.set_xlabel('Residue', fontsize=8)
ax.set_ylabel('Residue Frequency', fontsize=8)
ax.tick_params(labelsize=7)
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(f"{fig_dir}/accumulated_residue_frequency.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/accumulated_residue_frequency.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()


### Gliding versus Hopping Mechanism
def jaccard(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    similarity = float(intersection) / union
    distance = 1 - similarity
    return similarity, distance


jaccard_scores_more_less_plot_1 = []

for keys, values in dictionary_transitions_sorted.items():
    similarity = (jaccard(keys[0], keys[1]))[0]
    jaccard_scores_more_less_plot_1.append([values, similarity])

jaccard_scores_more_less_sorted_plot_1 = sorted(jaccard_scores_more_less_plot_1, key=lambda x: x[1], reverse=True)
total = sum([x[0] for x in jaccard_scores_more_less_plot_1])

normalized_data = [[x[0] / total, x[1]] for x in jaccard_scores_more_less_plot_1]

aggregated = defaultdict(float)
for value, key in normalized_data:
    aggregated[key] += value

aggregated_dict = dict(sorted(aggregated.items()))
x = list(aggregated_dict.keys())
y = list(aggregated_dict.values())

colors = sns.color_palette("Blues", len(x))
fig, ax = plt.subplots(figsize=(3.5, 3.5))
ax.bar(range(len(x)), y, color=colors, alpha=0.8)
ax.set_xticks(range(len(x)))
ax.set_xticklabels([f"{val:.2f}" for val in x], rotation=45, ha='right', fontsize=6)
ax.set_xlabel('Jaccard Similarity', fontsize=8)
ax.set_ylabel('Normalized Transition Frequency', fontsize=8)
ax.tick_params(labelsize=7)
sns.despine(ax=ax)
ax.get_legend() and ax.get_legend().remove()
plt.tight_layout()
plt.savefig(f"{fig_dir}/jaccard_similarity.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/jaccard_similarity.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()





#Plot part b
jaccard_scores_more_less = {}
count_more_1 = 0
count_less_70 = 0
count_less_60=0

for keys, values in dictionary_transitions_sorted.items():
    if count_more_1 >= 10 and count_less_70 >= 10 and count_less_60>= 10:
        break  

    j_similarity, j_distancejaccard = jaccard(list(keys[0]), list(keys[1]))
    k = tuple(keys)

    if count_more_1 < 10 and j_similarity == 1.0:
        jaccard_scores_more_less[k] = j_similarity,values
        count_more_1 += 1

    if count_less_70 < 10 and j_similarity <= 0.70:
        jaccard_scores_more_less[k] = j_similarity,values
        count_less_70 += 1


    elif count_less_60 < 10 and j_similarity <= 0.60:
        jaccard_scores_more_less[k] = j_similarity,values
        count_less_60 += 1

x_jac = list(jaccard_scores_more_less.keys())
y_jac= list(jaccard_scores_more_less.values())
#y_normalized_jac = y_jac
unique_states = set()
for k in x_jac:
    unique_states.add(k[0])
    unique_states.add(k[1])
state_to_num = {state: f"{i}" for i, state in enumerate(sorted(unique_states), 1)}

y_sum=sum(dictionary_transitions_sorted.values())
y_normalized_jac = [(sim, trans / y_sum) for sim, trans in y_jac]

transition_labels = [f"T{i+1}" for i in range(len(x_jac))]

transition_map = {
    f"T{i+1}": (state_to_num[k[0]], state_to_num[k[1]]) for i, k in enumerate(x_jac)
}

bar_colors = []
for sim, trans in y_normalized_jac:
    if sim == 1.0:
        bar_colors.append("#deebf7")
    elif 0.65 <= sim <= 0.70:
        bar_colors.append("#9ecae1")
    else:
        bar_colors.append("#3182bd")

bar_values = [trans for sim, trans in y_normalized_jac]

fig, ax = plt.subplots(figsize=(3.5, 3.5))
ax.bar(transition_labels, bar_values, color=bar_colors, width=0.7)

blue_patch_1 = mpatches.Patch(color="#deebf7", label="Similarity = 1.0")
blue_patch_2 = mpatches.Patch(color="#9ecae1", label="0.65 <= Similarity <= 0.70")
blue_patch_3 = mpatches.Patch(color="#3182bd", label="Other Similarity Values")
ax.legend(handles=[blue_patch_1, blue_patch_2, blue_patch_3], fontsize=5, frameon=False,
          loc='upper right')
ax.set_ylabel('Normalized Transition Counts', fontsize=8)
ax.set_xlabel('Transitions', fontsize=8)
ax.tick_params(labelsize=5)
ax.set_xticks(range(len(transition_labels)))
ax.set_xticklabels(transition_labels, rotation=90, fontsize=5)
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(f"{fig_dir}/jaccard_transitions.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/jaccard_transitions.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()

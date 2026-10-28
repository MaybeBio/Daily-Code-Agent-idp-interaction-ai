
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.linalg import expm
import pandas as pd
from collections import defaultdict, OrderedDict, Counter
import itertools
import pickle
import mdtraj as md
from itertools import islice
from math import nan, isnan
from itertools import combinations
from sklearn.preprocessing import normalize
from numpy import linalg as eig
from matplotlib.backends.backend_pdf import PdfPages
from numpy import random
import signal
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

distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta_gabis/d_24_t_com_avg.pkl', 'rb'))
distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta_gabis/d_24_t_closest.pkl', 'rb'))


w_com=0.9
w_closest=0.1
# w_com=0.7
# w_closest=0.3

distances_combined = (w_com * np.array(distances_com)) + (w_closest * np.array(distances_closest))

distance_threshold_combined = (0.75*w_com)+(0.45*w_closest)

distances=distances_combined

number_contact=[]
for ts, values in enumerate(distances):
    number_contact.append(sum(1 for x in values if x <= distance_threshold_combined))

print(f"Average number of contacts between c alpha and ligand com is {np.mean(number_contact)}")
uplet_types=list(range(2,6))

n_cols = 2
n_rows = int(len(uplet_types)/n_cols)
# axes = axes.flatten()  # Flatten for easier indexing

transiton_cutoff=2000

sample_size=12 #12/2

# blues = ['#b3cde0', '#6497b1', '#1f4579', '#0b2c55']  # Lighter to darker blues
# reds = ['#f7b7b7', '#f08f8f', '#d94e4e', '#b31f1f']  # Lighter to darker reds
blues = ['#d9d9d9', '#969696', '#525252', '#000000']
reds   = ['#fcaeae', '#fb6a6a', '#de2d26', '#7f0000']

# Setting up figure size and subplots
fig, axes = plt.subplots(n_rows, n_cols, figsize=(7, 3.5))
axes = axes.flatten()

for idx,uplet_type in enumerate(uplet_types):
    print(f'w_com/w_closest ratio is {w_com} /{w_closest} and using {transiton_cutoff} transition matrix size')
    # Initialize a dictionary to store counts of simultaneous contacts for each residue pair
    residue_pairs = list(combinations(range(42), uplet_type))  # all unique pairs of residues
    contact_counts_top_uplet_type_indices = {pair: 0 for pair in residue_pairs}

    num_timesteps = distances.shape[0]
    contact_counts_uplet_type_timesteps=[]

    for t in range(num_timesteps):
        close_residues = np.where(distances[t, :] < distance_threshold_combined)[0]
        close_residues_values = [x for x in distances[t] if x < distance_threshold_combined]

        paired = list(zip(close_residues_values, close_residues))
        sorted_pairs = sorted(paired, key=lambda x: x[0])
        top_uplet_type_indices = [pair[1] for pair in sorted_pairs[:uplet_type]]
    # print(close_residues)

        if len(top_uplet_type_indices) == uplet_type:
            contact_counts_uplet_type_timesteps.append(top_uplet_type_indices)


    #this sorts the values of sublists in order so that we don't have any repeats ([6,3,9) =[3,6,9]]) 
    unique_uplets_pre_process=[sorted(sublist) for sublist in contact_counts_uplet_type_timesteps]

    frequency = Counter(tuple(sorted(sublist)) for sublist in unique_uplets_pre_process)
    print(f'For {uplet_type}, there are {len(frequency)} unique pairs')

    #removes all uplets that are less than 100 in frequency
    #filtered_keys = [key for key, count in frequency.items() if count > 50] 

    #keeps top 1000 uplets
    filtered_keys = [item for item, count in frequency.most_common(transiton_cutoff)]

    #helps to transorm the sublists into usable format
    data_preprocessed = [tuple(x) for x in unique_uplets_pre_process]

    #removes all sublists that arent present in filtered_keys (eg- any sublist/uplet type that hasnt occured X times)
    data = [sublist for sublist in data_preprocessed if sublist in filtered_keys]
    
    print(f'For {uplet_type}, there were {len(data_preprocessed)} transitions, but now filtered there are {len(data)}')

    transition_matrix = np.zeros((len(filtered_keys), len(filtered_keys)), dtype=int)

    value_to_index =  {tuple(row): index for index, row in enumerate(filtered_keys)}


    for i in range(len(data) - 1):
        current_value = data[i]
        next_value = data[i + 1]
        
        sorted_current_value = tuple(sorted(current_value))
        sorted_next_value = tuple(sorted(next_value))
        
        transition_matrix[value_to_index[sorted_current_value], value_to_index[sorted_next_value]] += 1

    x_normed = normalize(transition_matrix, axis=1, norm='l1')

    flat_indices = np.unravel_index(
        np.argpartition(-transition_matrix.ravel(), len(transition_matrix))[:len(transition_matrix)],
        transition_matrix.shape
    )

    top_indices_values = [
        ((i, j), transition_matrix[i, j]) 
        for i, j in zip(flat_indices[0], flat_indices[1])]
    
    top_indices_values_sorted=sorted(top_indices_values, key=lambda x: x[1],reverse=True)  

    #now to look at radnom samples:
    int_sample_from=[random.randint(0, len(value_to_index) - 1) for x in range(2)]


    all_empirical_values = []
    all_predicted_values = []

    for state_number in (list(range(0, 2)) + int_sample_from):
        index_1 = top_indices_values_sorted[state_number][0][0]
        index_2 = top_indices_values_sorted[state_number][0][1]
        transition_quantity = top_indices_values_sorted[state_number][1]

        start_state = list(list(value_to_index.keys())[index_1])
        end_state = list(list(value_to_index.keys())[index_2])

        #print(f"Run {idx + 1}: Start State: {start_state}, End State: {end_state}")

        time_steps = np.arange(0, 20)
        empirical_values = []
        predicted_values = []

        for tau in time_steps:
            # Calculate theoretical n-step transition matrix
            theoretical_n_step_matrix = np.linalg.matrix_power(x_normed, tau)

            sorted_start_state_index = value_to_index[tuple(sorted(start_state))]
            sorted_end_state_index = value_to_index[tuple(sorted(end_state))]

            theoretical_prob = theoretical_n_step_matrix[sorted_start_state_index, sorted_end_state_index]

            # Calculate empirical n-step transition matrix
            num_states = transition_matrix.shape[0]
            empirical_n_step_matrix = np.zeros((num_states, num_states), dtype=float)

            for i in range(len(data) - tau):
                current_value = tuple(sorted(data[i]))
                next_value = tuple(sorted(data[i + tau]))

                current_index = value_to_index[current_value]
                next_index = value_to_index[next_value]

                empirical_n_step_matrix[current_index, next_index] += 1
                
            empirical_n_step_matrix = normalize(empirical_n_step_matrix, axis=1, norm='l1')
            empirical_prob = empirical_n_step_matrix[sorted_start_state_index, sorted_end_state_index]

            empirical_values.append(empirical_prob)
            predicted_values.append(theoretical_prob)

        all_empirical_values.append(empirical_values)
        all_predicted_values.append(predicted_values)

    mean_empirical = np.mean(all_empirical_values, axis=0)
    std_empirical = np.std(all_empirical_values, axis=0)

    mean_predicted = np.mean(all_predicted_values, axis=0)
    std_predicted = np.std(all_predicted_values, axis=0)

    ax = axes[idx]
    
    ax.plot(time_steps, mean_empirical, label='Mean Empirical', color=blues[3], lw=2)

    ax.fill_between(time_steps, mean_empirical - std_empirical, mean_empirical + std_empirical,
                    color=blues[1], alpha=0.4)
    ax.plot(time_steps, mean_predicted, label='Mean Predicted', color=reds[3], lw=2)

    ax.fill_between(time_steps, mean_predicted - std_predicted, mean_predicted + std_predicted,
                    color=reds[1], alpha=0.4)

    if idx >= (n_rows - 1) * n_cols:
        ax.set_xticks(time_steps)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(4))
        ax.tick_params(axis='x', labelsize=7)
    else:
        ax.tick_params(axis='x', labelbottom=False)

    ax.tick_params(axis='y', labelsize=7)
    if idx % n_cols != 0:
        ax.tick_params(axis='y', labelleft=False)

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

# One global legend
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(
    handles, labels,
    loc='upper center',
    bbox_to_anchor=(0.5, 1.05),
    ncol=2,
    fontsize=7,
    frameon=False
)

fig.supylabel("Transition Probability", fontsize=8, x=0.01)
fig.supxlabel("Time Step", fontsize=8, y=0.01)
plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig(f"{fig_dir}/ck_test.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/ck_test.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()

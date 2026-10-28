import numpy as np
from matplotlib.ticker import ScalarFormatter
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from tqdm import tqdm
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import io
from collections import defaultdict
import itertools
import pickle
import seaborn as sns
import mdtraj as md
from itertools import islice
from deeptime.util.data import TrajectoryDataset, TrajectoriesDataset
from deeptime.util.torch import MLP
from deeptime.decomposition.deep import VAMPNet
from deeptime.decomposition import VAMP
from deeptime.util.validation import implied_timescales, ck_test
from deeptime.plots import plot_implied_timescales, plot_ck_test
from torch.utils.data import DataLoader
import deeptime as dt
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics.cluster import adjusted_mutual_info_score
from deeptime.data import sqrt_model
from MDAnalysis.coordinates.base import Timestep
from MDAnalysis.coordinates.memory import MemoryReader
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

if torch.cuda.is_available():
    device = torch.device("cuda")
    torch.backends.cudnn.benchmark = True
else:
    device = torch.device("cpu")
torch.set_num_threads(12)


xtc=('/Users/adelielouet/Documents/science/AB_G5_original_simu_analysis/trajectories/Gabis_paper/traj_all-skip-0-noW_G5.xtc')
pdb=('/Users/adelielouet/Documents/science/AB_G5_original_simu_analysis/trajectories/Gabis_paper/template_G5.pdb')

t = md.load(xtc, top=pdb)
topology = t.topology

ligand_atom_indices = [636, 628, 650]
protein_ca_indices = t.top.select("protein and name CA")
table, bonds = topology.to_dataframe()

file = open('/Users/adelielouet/Documents/science/AB_G5_original_simu_analysis/trajectories/Gabis_paper/distances_residue_com_liga.pickle','rb')
distances = pickle.load(file)
file.close()

distance_t_40=(np.stack(distances, axis=1))
distance_t_40=distance_t_40.reshape(256128, 42)


## Clustering number ONE: dividing trajecotry into metastates

#input vairables
lag=10
nstates=4
lr=5e-3
epochs=30  #n_epochs=nb_epoch: Specifies the number of training epochs (iterations over the entire dataset).
batch_size=10000
#data1 = [Xi.detach().cpu().numpy()]  # Input data

normalized_features = distance_t_40#distance_features_t_120#distance_features_t_120 / distance_features_t_120.sum(axis=1, keepdims=True)

tensor_data = torch.tensor(normalized_features) #distance_features_t_120)

# Step 3: Detach from the graph and move to CPU if necessary
Xi = tensor_data.detach().cpu()

# Step 4: Convert to NumPy array and wrap it in a list
data1 = [Xi.numpy()]
dataset = TrajectoriesDataset.from_numpy(lag, data1)   # define some lag time
#dataset = dt.util.data.TrajectoryDataset(lagtime=tau, trajectory=data.astype(np.float32)) 


n_val = int(len(dataset)*.1) # Portion of data set as test set
train_data, val_data = torch.utils.data.random_split(dataset, [len(dataset) - n_val, n_val])    
d=nstates-1 # define nstates before; d=number of hidden layers in the network
dimsize = torch.ones(d+2,dtype=int) #size of each layer of the network
n_in = data1[0].shape[1] #n_in = number of inputs for the network
n_out = nstates # n_out= number of nodes in output layer = number of states
frac = 1/np.power((n_in/n_out),1/d) # low network model should be chosen
dimsize[0]= data1[0].shape[1] # n_in
dimsize[1] = int(np.ceil(frac*n_in)) # 1/d
dimsize[d+1] = nstates 

for h in range(2,d+1):
    dimsize[h] = int(np.ceil(frac*dimsize[h-1]))  # set nodes for each layer
print(n_in,n_out,dimsize)
lobe1 = MLP(units=dimsize, nonlinearity=nn.ELU,initial_batchnorm=True) # Define a lobe network
lobe = nn.Sequential(lobe1,nn.Softmax(dim=1)) # Output layer of the lobe ; fuzzy clustering
lobe = lobe.to(device=device) 
#print(lobe)
# lobe.to(device=device)
#print(lobe)

vampnet = VAMPNet(lobe=lobe, learning_rate=lr, device=device)  # lr = learning rate has to be defined; to start with 5e-3 is good

loader_train = DataLoader(train_data, batch_size=batch_size, shuffle=True)   #Batch size should be changed; better results with batchsize between 4000-10000
loader_val = DataLoader(val_data, batch_size=len(val_data), shuffle=False)

model = vampnet.fit(loader_train, n_epochs=epochs,
                   validation_loader=loader_val, progress=tqdm).fetch_model()


state_probabilities = model.transform(data1[0])   #converts input data into fuzzy clusters
assignments = state_probabilities.argmax(1) #gets discrete clusters from fuzzy probabilities
lagtimes = np.arange(5, 201 ,5, dtype=np.int32) # this range can be changed; try different possibilities
vamp_models = [VAMP(lagtime=lag, observable_transform=model).fit_fetch(data1) for lag in tqdm(lagtimes)] 
its_data = implied_timescales(vamp_models)  #gets timescales out; might be meaningless for clustering IDP binding; but worth checking how it varies over multiple trials

state_numbers=dict(pd.value_counts(assignments))
print(f'The state values are: {state_numbers}')

fig, ax = plt.subplots(figsize=(3.5, 3.5))
ax.plot(*vampnet.train_scores.T, label='Training')
ax.plot(*vampnet.validation_scores.T, label='Validation')
ax.set_xlabel('Step', fontsize=8)
ax.set_ylabel('Score', fontsize=8)
ax.tick_params(labelsize=7)
ax.legend(fontsize=7, frameon=False)
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(f"{fig_dir}/vampnet_scores_iter1.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/vampnet_scores_iter1.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()


### This is to add the ts to each assignemnt
frames_assignments_labeled = {key: [] for key in set(assignments)}

for timestep, (cluster, feature) in enumerate(zip(assignments, normalized_features)):
    frames_assignments_labeled[cluster].append((feature, timestep))


## Clustering number two: dividing trajecotry into microstates:

frames_assignments = {key: [] for key in set(assignments)}

for x, y in zip(assignments, normalized_features):
    frames_assignments[x].append(y)


lag=5
nstates=4
lr=5e-3
epochs=30  #n_epochs=nb_epoch: Specifies the number of training epochs (iterations over the entire dataset).
batch_size=500

count=0
assignments_16_states=[]
values_16_states=[]
ts_seperated=[]
for x in (range(0,nstates)):
    tensor_data = torch.tensor(frames_assignments[x]) 

    Xi = tensor_data.detach().cpu()
    data1 = [Xi.numpy()]
    dataset = TrajectoriesDataset.from_numpy(lag, data1) 

    n_val = int(len(dataset)*.1) # Portion of data set as test set
    train_data, val_data = torch.utils.data.random_split(dataset, [len(dataset) - n_val, n_val])    
    d=nstates-1 # define nstates before; d=number of hidden layers in the network
    dimsize = torch.ones(d+2,dtype=int) #size of each layer of the network
    n_in = data1[0].shape[1] #n_in = number of inputs for the network
    n_out = nstates # n_out= number of nodes in output layer = number of states
    frac = 1/np.power((n_in/n_out),1/d) # low network model should be chosen
    dimsize[0]= data1[0].shape[1] # n_in
    dimsize[1] = int(np.ceil(frac*n_in)) # 1/d
    dimsize[d+1] = nstates 

    for h in range(2,d+1):
        dimsize[h] = int(np.ceil(frac*dimsize[h-1]))  # set nodes for each layer
    lobe1 = MLP(units=dimsize, nonlinearity=nn.ELU,initial_batchnorm=True) # Define a lobe network
    lobe = nn.Sequential(lobe1,nn.Softmax(dim=1)) # Output layer of the lobe ; fuzzy clustering
    lobe = lobe.to(device=device) 
    #print(lobe)
    # lobe.to(device=device)
    #print(lobe)

    vampnet = VAMPNet(lobe=lobe, learning_rate=lr, device=device)  # lr = learning rate has to be defined; to start with 5e-3 is good

    loader_train = DataLoader(train_data, batch_size=batch_size, shuffle=True)   #Batch size should be changed; better results with batchsize between 4000-10000
    loader_val = DataLoader(val_data, batch_size=len(val_data), shuffle=False)

    model = vampnet.fit(loader_train, n_epochs=epochs,
                    validation_loader=loader_val, progress=tqdm).fetch_model()

    state_probabilities = model.transform(data1[0])   #converts input data into fuzzy clusters
    assignments = state_probabilities.argmax(1) #gets discrete clusters from fuzzy probabilities
    lagtimes = np.arange(5, 201 ,5, dtype=np.int32) # this range can be changed; try different possibilities
    vamp_models = [VAMP(lagtime=lag, observable_transform=model).fit_fetch(data1) for lag in tqdm(lagtimes)] 
    its_data = implied_timescales(vamp_models)  #gets timescales out; might be meaningless for clustering IDP binding; but worth checking how it varies over multiple trials

    assignments_plus_4 = [x + count for x in assignments]
    print(set(assignments_plus_4),set(assignments))
    print(f'count is {count}')
    count+=4

    assignments_16_states.append(assignments_plus_4)
    values_16_states.append(frames_assignments[x])

assignments_16_states = list(itertools.chain(*assignments_16_states))
values_16_states_flat = np.array(list(itertools.chain.from_iterable(values_16_states)))

#values_16_states_flat_reversed = values_16_states_flat.reshape(len(distance_matrix_t_40_3), len(distance_matrix_t_40_3[0]), len(distance_matrix_t_40_3[0][0])) # for 3 ligand points
values_16_states_flat_reversed = values_16_states_flat.reshape(len(distance_t_40), len(distance_t_40[0])) #for com

## This is to assign the corresponding timesteps to each cluster
ts_shortcut=[]
for keys,values in frames_assignments_labeled.items():
    for o in values:
        ts_shortcut.append(o[1])

frames_assignments_labeled_16 = {key: [] for key in set(assignments_16_states)}

for (cluster, ts) in zip(assignments_16_states, ts_shortcut):
    frames_assignments_labeled_16[cluster].append(ts)

state_numbers=dict(pd.value_counts(assignments_16_states))
print(f'The state values are: {state_numbers}')


### This is to add the ts to each assignemnt
frames_assignments_labeled = {key: [] for key in set(assignments)}

for timestep, (cluster, feature) in enumerate(zip(assignments, normalized_features)):
    frames_assignments_labeled[cluster].append((feature, timestep))


## Clustering number two: dividing trajecotry into microstates:

frames_assignments = {key: [] for key in set(assignments)}

for x, y in zip(assignments, normalized_features):
    frames_assignments[x].append(y)

def clean(array):
    cleaned_arrays = []

    for arr in array:
        cleaned_array = np.array(arr).tolist()
        cleaned_arrays.append(cleaned_array)

    return(cleaned_arrays)


######## Vampnet Analysis
# Training/Validation score
fig, ax = plt.subplots(figsize=(3.5, 3.5))
ax.plot(*vampnet.train_scores.T, label='Training')
ax.plot(*vampnet.validation_scores.T, label='Validation')
ax.set_xlabel('Step', fontsize=8)
ax.set_ylabel('Score', fontsize=8)
ax.tick_params(labelsize=7)
ax.legend(fontsize=7, frameon=False)
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(f"{fig_dir}/vampnet_scores_iter2.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/vampnet_scores_iter2.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()


#put the new assingmnets back in: 

# values_16_states_flat=values_16_states_flat+state_17
# assignments_16_states=assignments_16_states+assignments_3

case_assignments = {key: [] for key in set(assignments_16_states)}

for x, y in zip(assignments_16_states, distance_t_40):
    case_assignments[x].append(y)

fig, axs = plt.subplots(4, 4, figsize=(7, 7))
axs = axs.flatten()
num_states = len(set(assignments_16_states))

for idx, (key, values) in enumerate(case_assignments.items()):
    distance_cutoff_mat = []
    for x in values:
        y = np.zeros((len(x), 1))
        print(x.shape)
        for i in np.ndindex(x.shape):
            if x[i] <= 0.6:
                y[i] = 1
        distance_cutoff_mat.append(y)

    distance_cutoff_mat_flat = [y.flatten() for y in distance_cutoff_mat]
    distance_cutoff_mat_flat = clean(distance_cutoff_mat_flat)
    distance_cutoff_mat = np.array(distance_cutoff_mat_flat)

    co_occurrence_matrix_vn = np.zeros((distance_cutoff_mat.shape[1], distance_cutoff_mat.shape[1]))
    for frame in distance_cutoff_mat:
        for i in range(len(frame)):
            for j in range(i, len(frame)):
                if frame[i] == 1 and frame[j] == 1:
                    co_occurrence_matrix_vn[i, j] += 1
                    if i != j:
                        co_occurrence_matrix_vn[j, i] += 1

    sns.heatmap(co_occurrence_matrix_vn, cmap="rocket_r", square=True, ax=axs[idx], cbar=False)
    axs[idx].set_title(f'Cluster {idx+1}', fontsize=7, pad=3)
    axs[idx].tick_params(axis='x', labelsize=5)
    axs[idx].tick_params(axis='y', labelsize=5)
    xticks = axs[idx].get_xticks()
    yticks = axs[idx].get_yticks()
    axs[idx].set_xticks(xticks[::10])
    axs[idx].set_yticks(yticks[::10])
    if idx in [12, 13, 14, 15]:
        axs[idx].set_xlabel('Residue Index', fontsize=7)
    if idx in [0, 4, 8, 12]:
        axs[idx].set_ylabel('Residue Index', fontsize=7)

plt.tight_layout()
plt.savefig(f"{fig_dir}/vampnet_heatmap.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/vampnet_heatmap.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()

########## Radius of Gyration Analysis:
rg_data = {}

num_clusters = len(frames_assignments_labeled_16)
palette = sns.color_palette("Blues", n_colors=num_clusters)

fig, ax = plt.subplots(figsize=(3.5, 3.5))
for idx, (cluster_id, timesteps) in enumerate(frames_assignments_labeled_16.items()):
    selected_traj = t.slice(timesteps)
    rg_values = md.compute_rg(selected_traj)
    rg_data[cluster_id] = rg_values
    sns.kdeplot(rg_values, label=f'Cluster {cluster_id+1}', color=palette[idx], ax=ax)

ax.set_xlabel(r'Radius of Gyration $R_g$ (nm)', fontsize=8)
ax.set_ylabel('Density', fontsize=8)
ax.tick_params(labelsize=7)
ax.legend(fontsize=6, frameon=False, loc='upper right')
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(f"{fig_dir}/vampnet_kde.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/vampnet_kde.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()


# Shows the residues involved in each group
ligand_atom_indices = [628]  # You can add more ligand atom indices here
protein_ca_indices = t.top.select("protein and name CA")
num_clusters = len(frames_assignments_labeled_16)
num_cols = 4
num_rows = (num_clusters + num_cols - 1) // num_cols  # ensures no extra empty row

fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(6 * num_cols, 4.5 * num_rows))
axes = axes.flatten()

cluster_dict = {}

for i, (cluster_idx, timesteps) in enumerate(frames_assignments_labeled_16.items()):
    all_ts = []

    for timestep in timesteps:
        selected_traj = t.slice(timestep)
        pairs = list(itertools.product(ligand_atom_indices, protein_ca_indices))
        distances_ligand = md.compute_distances(selected_traj, pairs, periodic=False)
        all_ts.append(distances_ligand[0])

    cluster_dict[cluster_idx] = all_ts

# cmap = LinearSegmentedColormap.from_list("abeta_blue", ["#d1e5f0", "#2166ac"]) 
fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(5.5, 3.5))
axes = axes.flatten()

for cluster_idx, all_ts in cluster_dict.items():
    mean = np.mean(all_ts, axis=0)
    stdev = np.std(all_ts, axis=0)

    i = cluster_idx
    ax = axes[i]

    ax.plot(range(len(mean)), mean, label=f"Cluster {cluster_idx}", color="#4d8ac7", linewidth=1)
    ax.fill_between(range(len(mean)), mean - stdev, mean + stdev, alpha=0.3, color="#4d8ac7")

    ax.set_ylim(0, 10)
    ax.set_xticks(range(0, len(protein_ca_indices), 5))
    ax.tick_params(axis='y', labelsize=5)
    ax.tick_params(axis='x', labelsize=5)

    # Show x-axis ticks only on bottom row
    if i >= (num_rows - 1) * num_cols:
        ax.set_xticklabels([protein_ca_indices[j] for j in range(0, len(protein_ca_indices), 5)], rotation=45, fontsize=5)
    else:
        ax.tick_params(axis='x', bottom=True, labelbottom=False)

    # Show y-axis ticks only on leftmost column
    if i % num_cols == 0:
        ax.tick_params(axis='y', labelleft=True)
    else:
        ax.tick_params(axis='y', labelleft=False)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

# Turn off any unused axes
for j in range(i + 1, len(axes)):
    axes[j].axis("off")

fig.supxlabel(r'C$\alpha$ Index', fontsize=8, y=0.01)
fig.supylabel('Average Contact Distance', fontsize=8, x=0.01)
plt.tight_layout()
plt.subplots_adjust(hspace=0.15, wspace=0.1)
plt.savefig(f"{fig_dir}/vampnet_lineplot.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/vampnet_lineplot.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()

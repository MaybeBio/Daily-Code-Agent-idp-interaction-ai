# Vamonet refined clustering, trying to do two subsequent clustering via vmaonet
# attempt to seperate data prior to inputting in vampnet

import numpy as np
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
import itertools

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

import os
fig_dir = "/Users/adelielouet/Documents/science/Reports/MSM_paper/mv_copies_march_8/figures_hi_res"
os.makedirs(fig_dir, exist_ok=True)

if torch.cuda.is_available():
    device = torch.device("cuda")
    torch.backends.cudnn.benchmark = True
else:
    device = torch.device("cpu")
torch.set_num_threads(12)

from deeptime.data import sqrt_model


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

asssingnments_dict={}
for run_number in range(0,5):
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

    for x in sorted(frames_assignments.keys()):
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
    values_16_states_flat_reversed = values_16_states_flat.reshape(len(distance_t_40), len(distance_t_40[0])) #for com

    ts_shortcut=[]
    for keys,values in frames_assignments_labeled.items():
        for o in values:
            ts_shortcut.append(o[1])

    frames_assignments_labeled_16 = {key: [] for key in set(assignments_16_states)}

    
    frames_assignments_labeled_16 = {key: [] for key in set(assignments_16_states)}

    for (cluster, ts) in zip(assignments_16_states, ts_shortcut):
        frames_assignments_labeled_16[cluster].append(ts)

    clusters_dict=frames_assignments_labeled_16.copy()
    print(f'round {run_number} done')

    max_frame = max(frame for frames in clusters_dict.values() for frame in frames)

    assignments_16_states = [-1] * (max_frame + 1)
    for cluster, frames in clusters_dict.items():
        for frame in frames:
            assignments_16_states[frame] = int(cluster) 

    asssingnments_dict[run_number]=assignments_16_states


ami_scores = []
for i, j in itertools.combinations(range(len(asssingnments_dict)), 2):
    score = adjusted_mutual_info_score(asssingnments_dict[i], asssingnments_dict[j])
    ami_scores.append(score)

fig, ax = plt.subplots(figsize=(3.5, 3.5))
sns.kdeplot(ami_scores, ax=ax)
ax.set_xlabel('AMI Score', fontsize=8)
ax.set_ylabel('Density', fontsize=8)
ax.tick_params(labelsize=7)
ax.legend(fontsize=6, frameon=False, loc='upper right')
sns.despine(ax=ax)
plt.tight_layout()
plt.savefig(f"{fig_dir}/vampnet_ami.pdf", format='pdf', bbox_inches='tight', pad_inches=0.05)
plt.savefig(f"{fig_dir}/vampnet_ami.png", format='png', dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close()

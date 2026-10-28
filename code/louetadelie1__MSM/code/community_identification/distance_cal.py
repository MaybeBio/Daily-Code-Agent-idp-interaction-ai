import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import mdtraj as md
import itertools
import glob
import os
import sys

### The following is fitted for alpha_synuclein from DE Shaw. Made to handle dcd files
current_dir=sys.argv[1]
protein_name=sys.argv[2]
zinc_name=sys.argv[3]
simu_length=sys.argv[4]

xtc=(f'{current_dir}/{simu_length}')
pdb=(f'{current_dir}/template.pdb')
#xtc=('/home/al2108/rds/hpc-work/alpha_syn/de_analysis/combined_frames.xtc')

#protein_name=pdb.split('/')[-2]
#zinc_name=pdb.split('/')[-1]

traj = md.load(xtc, top=pdb)
topology = traj.topology
protein = traj.top.select('protein')
ligand_indices = topology.select('not protein')

# ### CALCULATING THE AVERAGE BETWEEN COM OF LIGAND AND COM OF EACH RESIUDUE && TAKING SMALLEST DISTANCE BETWEEN LIGAND AND PROTEIN AA
com_calcualtion_all=[]
clostest_calcualtion_all=[]

com_calcualtion=[]
clostest_calcualtion=[]
for resid in range(0,traj.topology.n_residues-1): #-1 beacus -1 is ligand
    res_x=topology.select(f'resid {resid}')
    pairs = list(itertools.product(res_x, ligand_indices))
    distances=md.compute_distances(traj, pairs, periodic=True, opt=True)
    com_calcualtion.append((np.average(distances,axis=1)).tolist())
    clostest_calcualtion.append(np.min(distances,axis=1).tolist())

d_24_t=np.reshape(com_calcualtion,((traj.topology.n_residues-1), len(traj))).T
d_24_t_closest=np.reshape(clostest_calcualtion,((traj.topology.n_residues-1), len(traj))).T

dir=f'/ptmp/adlouet/mass_produce_md_simulations/output_post_processed/{protein_name}/{zinc_name}/distances/{simu_length}/'
os.makedirs(dir, exist_ok=True)

with open(f'{dir}d_24_t_closest.pkl', 'wb') as file:
          pickle.dump(d_24_t_closest, file)

with open(f'{dir}d_24_t_com_avg.pkl', 'wb') as file:
          pickle.dump(d_24_t, file)

import itertools
import math
import smplotlib
import random
from collections import Counter, defaultdict, OrderedDict
from itertools import combinations, islice
from math import nan, isnan
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy.linalg import expm, eig, norm
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
from sklearn import preprocessing
from sklearn.preprocessing import MinMaxScaler, normalize
import networkx as nx
from community import community_louvain
import community as c  # You may only need one of these
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
import mdtraj as md
import pickle
import operator
import numpy as np
from numpy.linalg import eig
from tqdm import tqdm
import time

def sigmoid_magnification(val, A, k, x0):
    return A / (1 + math.exp(-k * (val - x0)))


def cleanup(d):
    return {tuple(int(k) for k in key): float(value[0]) for key, value in d.items()}

# === Protein-specific settings == #

# proteins = ['full_alpha_syn_fau']#,'alpha_syn_lig_50','alpha_syn_lig_12','gcpr_ZINC000052011510','gcpr_ZINC000052011511','gcpr_ZINC000058184512']
# proteins = ['gcpr_ZINC000052011510','gcpr_ZINC000052011511','gcpr_ZINC000058184512']
# proteins = [
# 'alpha_syn_lig_12', 'alpha_syn_lig_20', 
#    'alpha_syn_lig_26', 'alpha_syn_lig_30', 'alpha_syn_lig_4','alpha_syn_lig_40','alpha_syn_lig_50','abeta_gabis'
#     'medin_cm10','medin_cm8','medin_urea',
# 'abeta_d8','abeta_d4','abeta_g5_new_protocol'
# ,'abeta_ph5_298k','abeta_ph7_278k','abeta_ph5_298k'
# ]
proteins = [
    'medin_cm10'
]

for protein_name in proteins:
    try:
        print(f"Processing {protein_name}...")

        # === Protein-specific settings == #
        if protein_name == 'abeta_gabis':
            pdb = '/Users/adelielouet/Documents/science/AB_G5_original_simu_analysis/trajectories/Gabis_paper/template_G5.pdb'
            label="abeta_gabis"
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta_gabis/d_24_t_com_avg.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta_gabis/d_24_t_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
        
        elif protein_name == 'alpha_syn_lig_40':
            pdb = '/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_40/ligand_40_alpha_syn_c_term.pdb'
            # xtc = '/Users/adelielouet/Documents/science/medin/urea/urea_post_run/concatenate_urea.xtc'
            # resname_ligand = "1UNL"
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/alpha_syn_lig_40/distances/t_resnum_com.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/alpha_syn_lig_40/distances/t_resnum_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True

        elif protein_name == 'alpha_syn_lig_50':
            pdb = '/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_50/ligand_50_alpha_syn_c_term.pdb'
            # xtc = '/Users/adelielouet/Documents/science/medin/urea/urea_post_run/concatenate_urea.xtc'
            # resname_ligand = "1UNL"
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/alpha_syn_lig_50/distances/t_resnum_com.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/alpha_syn_lig_50/distances/t_resnum_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True

        elif protein_name == 'alpha_syn_lig_12':
            pdb = '/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_12/DESRES-Trajectory_jacs2022-5447858-no-water-glue/jacs2022-5447858-no-water-glue/lig_12_cterm.pdb'
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_12/DESRES-Trajectory_jacs2022-5447858-no-water-glue/jacs2022-5447858-no-water-glue/t_resnum_com.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_12/DESRES-Trajectory_jacs2022-5447858-no-water-glue/jacs2022-5447858-no-water-glue/t_resnum_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True

        elif protein_name == 'full_alpha_syn_fau':
            pdb = '/Users/adelielouet/Documents/science/alpha_syn/pbmetad/DE_SHAW_traj/system.pdb'
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/full_alpha_syn_fau/d_24_t.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/full_alpha_syn_fau/d_24_t_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=5
            label='Full Length \u03B1-syn Fasudil' 
            pos_scaling=True



####################### Remove asap:
        elif protein_name == 'abeta_gabis':
            pdb = '/Users/adelielouet/Documents/science/AB_G5_original_simu_analysis/trajectories/Gabis_paper/template_G5.pdb'
            label="abeta_gabis"
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta_gabis/d_24_t_com_avg.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta_gabis/d_24_t_closest.pkl', 'rb'))
            # w_closest =  0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=5 #this is only for visauls, change it back for pop equi
            # w_com = 0.45
            # w_closest = 0.55
            w_com = 0.9 #this is only for visauls, change it back for pop equi
            w_closest = 0.1 #this is only for visauls, change it back for pop equi
            
        elif protein_name == 'medin_cm8':
            pdb = '/Users/adelielouet/Documents/science/medin/cm8/cm8_post_run/protein_ligand.gro'
            xtc = '/Users/adelielouet/Documents/science/medin/cm8/cm8/concatenate_cm8.xtc'
            w_file='/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_cm8/weights/COLVAR_REWEIGHT'
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_cm8/distances/d_24_t_com_avg.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_cm8/distances/d_24_t_closest.pkl', 'rb'))
            resname_ligand = "1UNL"
            special_char = 'Medin_cm8'
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            w_com = 0.3
            w_closest = 0.7
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'medin_cm10':
            pdb = '/Users/adelielouet/Documents/science/medin/cm10/cm10_post_run/protein_ligand.gro'
            xtc = '/Users/adelielouet/Documents/science/medin/cm10/cm10/concatenate_cm10.xtc'
            w_file='/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_cm10/weights/COLVAR_REWEIGHT'
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_cm10/distances/d_24_t_com_avg.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_cm10/distances/d_24_t_closest.pkl', 'rb'))
            resname_ligand = "1UNL"
            special_char = 'Medin_cm10'
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            w_com = 0.3
            w_closest = 0.7
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'medin_urea':
            pdb = '/Users/adelielouet/Documents/science/medin/urea/urea_post_run/protein_urea.gro'
            xtc = '/Users/adelielouet/Documents/science/medin/urea/urea_post_run/concatenate_urea.xtc'
            resname_ligand = "1UNL"
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_urea/distances/d_24_t_com_avg.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_urea/distances/d_24_t_closest.pkl', 'rb'))
            w_file = None
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'alpha_syn_lig_40':
            pdb = '/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_40/ligand_40_alpha_syn_c_term.pdb'
            # xtc = '/Users/adelielouet/Documents/science/medin/urea/urea_post_run/concatenate_urea.xtc'
            # resname_ligand = "1UNL"
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/alpha_syn_lig_40/distances/t_resnum_com.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/alpha_syn_lig_40/distances/t_resnum_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True

        elif protein_name == 'alpha_syn_lig_50':
            pdb = '/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_50/ligand_50_alpha_syn_c_term.pdb'
            # xtc = '/Users/adelielouet/Documents/science/medin/urea/urea_post_run/concatenate_urea.xtc'
            # resname_ligand = "1UNL"
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/alpha_syn_lig_50/distances/t_resnum_com.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/alpha_syn_lig_50/distances/t_resnum_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'alpha_syn_lig_12':
            pdb = '/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_12/DESRES-Trajectory_jacs2022-5447858-no-water-glue/jacs2022-5447858-no-water-glue/lig_12_cterm.pdb'
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_12/DESRES-Trajectory_jacs2022-5447858-no-water-glue/jacs2022-5447858-no-water-glue/t_resnum_com.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_12/DESRES-Trajectory_jacs2022-5447858-no-water-glue/jacs2022-5447858-no-water-glue/t_resnum_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'alpha_syn_lig_20':
            pdb = '/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_20/DESRES-Trajectory_jacs2022-5447842-no-water-glue/jacs2022-5447842-no-water-glue/lig_20_cterm.pdb'
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_20/DESRES-Trajectory_jacs2022-5447842-no-water-glue/jacs2022-5447842-no-water-glue/t_resnum_com.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_20/DESRES-Trajectory_jacs2022-5447842-no-water-glue/jacs2022-5447842-no-water-glue/t_resnum_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'alpha_syn_lig_26':
            pdb = '/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_26/DESRES-Trajectory_jacs2022-5447843-no-water-glue/jacs2022-5447843-no-water-glue/lig_26_cterm.pdb'
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_26/DESRES-Trajectory_jacs2022-5447843-no-water-glue/jacs2022-5447843-no-water-glue/t_resnum_com.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_26/DESRES-Trajectory_jacs2022-5447843-no-water-glue/jacs2022-5447843-no-water-glue/t_resnum_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'alpha_syn_lig_30':
            pdb = '/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_30/DESRES-Trajectory_jacs2022-5447857-no-water-glue/jacs2022-5447857-no-water-glue/lig_30_cterm.pdb'
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_30/DESRES-Trajectory_jacs2022-5447857-no-water-glue/jacs2022-5447857-no-water-glue/t_resnum_com.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_30/DESRES-Trajectory_jacs2022-5447857-no-water-glue/jacs2022-5447857-no-water-glue/t_resnum_closest.pkl', 'rb'))
            w_com =  0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'alpha_syn_lig_4':
            pdb = '/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_4/DESRES-Trajectory_jacs2022-12293914-no-water-glue/jacs2022-12293914-no-water-glue/lig_4_cterm.pdb'
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_4/DESRES-Trajectory_jacs2022-12293914-no-water-glue/jacs2022-12293914-no-water-glue/t_resnum_com.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/lig_4/DESRES-Trajectory_jacs2022-12293914-no-water-glue/jacs2022-12293914-no-water-glue/t_resnum_closest.pkl', 'rb'))
            w_com =  0.3
            w_closest =  0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'abeta_d4':
            pdb = '/Users/adelielouet/Documents/science/shared_files/maria/trajectories/D4/complex_noW_3.gro'
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta_D4/distances/d_24_t_com_avg.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta_D4/distances/d_24_t_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'abeta_d8':
            pdb = '/Users/adelielouet/Documents/science/shared_files/maria/trajectories/D8/complex_noW_2.gro'
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta_D8/distances/d_24_t_com_avg.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta_D8/distances/d_24_t_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest =0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'abeta_g5_new_protocol':
            pdb = '/Users/adelielouet/Documents/science/shared_files/maria/trajectories/G5_new_proto/complex_noW_4.gro'
            distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta_G5_new_protocol/distances/d_24_t_com_avg.pkl', 'rb'))
            distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta_G5_new_protocol/distances/d_24_t_closest.pkl', 'rb'))
            w_com = 0.3
            w_closest = 0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'ph_5_278k':
            pdb = '/Users/adelielouet/Documents/science/dd_proj/abeta_ph_simulations/ph5_278k/Complex_AB.pdb'
            label=protein_name.split('abeta_')[1]
            distances_com = pickle.load(open(f'//Users/adelielouet/Documents/science/dd_proj/abeta_ph_simulations/pickled_files/{label}/distances_2/d_24_t_com_avg.pkl', 'rb'))
            distances_closest = pickle.load(open(f'//Users/adelielouet/Documents/science/dd_proj/abeta_ph_simulations/pickled_files/{label}/distances_2/d_24_t_closest.pkl', 'rb'))
            w_com =  0.3
            w_closest =  0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'ph_7_278k':
            pdb = '/Users/adelielouet/Documents/science/dd_proj/abeta_ph_simulations/ph5_278k/Complex_AB.pdb'
            label=protein_name.split('abeta_')[1]
            distances_com = pickle.load(open(f'//Users/adelielouet/Documents/science/dd_proj/abeta_ph_simulations/pickled_files/{label}/distances_2/d_24_t_com_avg.pkl', 'rb'))
            distances_closest = pickle.load(open(f'//Users/adelielouet/Documents/science/dd_proj/abeta_ph_simulations/pickled_files/{label}/distances_2/d_24_t_closest.pkl', 'rb'))
            w_com =  0.3
            w_closest =  0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'ph_5_298k':
            pdb = '/Users/adelielouet/Documents/science/dd_proj/abeta_ph_simulations/ph5_278k/Complex_AB.pdb'
            label=protein_name.split('abeta_')[1]
            distances_com = pickle.load(open(f'//Users/adelielouet/Documents/science/dd_proj/abeta_ph_simulations/pickled_files/{label}/distances_2/d_24_t_com_avg.pkl', 'rb'))
            distances_closest = pickle.load(open(f'//Users/adelielouet/Documents/science/dd_proj/abeta_ph_simulations/pickled_files/{label}/distances_2/d_24_t_closest.pkl', 'rb'))
            w_com =  0.3
            w_closest =  0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

        elif protein_name == 'abeta_ph5_278_v2':
            pdb = '/Users/adelielouet/Documents/science/dd_proj/abeta_ph_simulations/ph5_278k/Complex_AB.pdb'
            label=protein_name.split('abeta_')[1]
            distances_com = pickle.load(open(f'//Users/adelielouet/Documents/science/dd_proj/abeta_ph_simulations/pickled_files/{label}/distances_2/d_24_t_com_avg.pkl', 'rb'))
            distances_closest = pickle.load(open(f'//Users/adelielouet/Documents/science/dd_proj/abeta_ph_simulations/pickled_files/{label}/distances_2/d_24_t_closest.pkl', 'rb'))
            w_com =  0.3
            w_closest =  0.7
            w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0
            combined_threshold=True
            uplet_type=3

 ######################## Remove asap


        else:
            print(f"Skipping {protein_name}: configuration not defined.")
            continue 

        w_eq /= w_eq.sum()  # Normalize so weights sum to 1
        distances_combined = (w_com * np.array(distances_com)) + (w_closest * np.array(distances_closest))
        #distance_threshold_combined = (0.75 * w_com) + (0.45 * w_closest)
        distance_threshold_combined=0.55
        number_contact = np.array([np.sum(values <= distance_threshold_combined) for values in distances_combined])
        print(f"distance is {distance_threshold_combined}")
        ##########

        bound_mask = number_contact != 0
        unbound_mask = ~bound_mask

        avg_contacts_when_bound = np.sum(w_eq[bound_mask] * number_contact[bound_mask]) / np.sum(w_eq[bound_mask])

        fraction_bound_weighted = np.sum(w_eq[bound_mask])

        total_bound_weighted = np.sum(w_eq[bound_mask]) * len(number_contact)

        print(fraction_bound_weighted)

        # === Generate Transition Matrix == #
        ### This is where you decide the -uplet type
        uplet_type=3

        u = mda.Universe(pdb)#, xtc)

        protein_residues = u.select_atoms("protein").residues
        num_residues = len(protein_residues)

        residue_pairs = list(combinations(range(num_residues), uplet_type))  # all unique pairs of residues
        contact_counts_top_uplet_type_indices = {pair: 0 for pair in residue_pairs}

        distances=distances_combined
        num_timesteps = distances.shape[0]

        contact_counts_uplet_type_timesteps=[]
        contact_counts_uplet_type_timesteps_inlcuding_0=[]

        w_eq_filtered_zeros=[]
        for t,w in zip(range(num_timesteps),w_eq):
            if w!=0.0:
                close_residues = np.where(distances[t, :] < distance_threshold_combined)[0]
                close_residues_values = [x for x in distances[t] if x < distance_threshold_combined]
                
                paired = list(zip(close_residues_values, close_residues))
                sorted_pairs = sorted(paired, key=lambda x: x[0])
                top_uplet_type_indices = [pair[1] for pair in sorted_pairs[:uplet_type]]
            # print(close_residues)
                if len(top_uplet_type_indices) == uplet_type:
                    contact_counts_uplet_type_timesteps.append(top_uplet_type_indices)
                    contact_counts_uplet_type_timesteps_inlcuding_0.append(top_uplet_type_indices)
                    w_eq_filtered_zeros.append(w)
                else:
                    contact_counts_uplet_type_timesteps_inlcuding_0.append(([-1]*uplet_type))
            else:
                print('no state exists when weighted')
                contact_counts_uplet_type_timesteps_inlcuding_0.append(([-1]*uplet_type))

        #this sorts the values of sublists in order so that we don't have any repeats ([6,3,9) =[3,6,9]]) 
        unique_uplets_pre_process=[sorted(sublist) for sublist in contact_counts_uplet_type_timesteps]
        unique_uplets_pre_process_0=[sorted(sublist) for sublist in contact_counts_uplet_type_timesteps_inlcuding_0]

        frequency = Counter(tuple(sorted(sublist)) for sublist in unique_uplets_pre_process_0)
        print(f'For {uplet_type}, there are {len(frequency)} unique pairs')

        ####################### Here you need to remove some of the pairs or else you won't be abel to converge the transition matrix
        # keys = [str(k) for k in frequency.keys()]
        # values = sorted(list(frequency.values()), reverse=True)

        # del values[0]
        # del keys[0]

        # n = len(keys)
        # x = np.arange(1, n + 1)  # x positions from 1..n

        # plt.figure(figsize=(12, 6))
        # plt.bar(x, values)

        # step = max(1, n // 5)
        # xtick_positions = np.arange(1, n + 1, step)

        # plt.xticks(xtick_positions)   # numeric ticks only

        # plt.xlabel(f"{uplet_type}-tuple index")
        # plt.ylabel("Frequency")
        # plt.title(f"Frequencies of {uplet_type}-tuple combinations")
        # plt.tight_layout()
        # plt.show()

        # the above plot is to better understand the distribution to see which ones we need to keep to only keep the top performers
        sorted_x = sorted(frequency.items(), key=operator.itemgetter(1),reverse=True)

        filtered_keys=[item for item, count in frequency.items()]

        data_preprocessed = [tuple(x) for x in unique_uplets_pre_process_0]

        data = [sublist for sublist in data_preprocessed if sublist in filtered_keys]#[filtered_keys]
        value_to_index =  {tuple(row): index for index, row in enumerate(filtered_keys)}#(filtered_keys)

        print(len(filtered_keys),len(value_to_index))

        ### Need for an unbound state to avoid absorbing states- start with unbound and edn with undbound, or else risk having the leaving state absorbed
        unbound_state =[-1] * uplet_type
        data.insert(0, unbound_state)
        data.append(unbound_state)

        print(len(filtered_keys),len(value_to_index))
        unbound_count=contact_counts_uplet_type_timesteps_inlcuding_0.count(unbound_state)
        unbound_fraction=(1-(unbound_count/len(contact_counts_uplet_type_timesteps_inlcuding_0)))*100
        print(f'The ligand is bound {unbound_fraction}% of the time')

        ############## bruh ##############
        sorted_x = dict(sorted(frequency.items(), key=operator.itemgetter(1), reverse=True))
        cut_size = int(len(sorted_x) / 1)
        cut_frequency = dict(list(sorted_x.items())[:cut_size])

        transition_matrix = np.zeros(
            (cut_size, cut_size),
            dtype=np.float64
        )

        filtered_data = [tuple(sorted(x)) for x in data if tuple(sorted(x)) in cut_frequency]

        cut_keys = list(cut_frequency.keys())
        value_to_index = {k: i for i, k in enumerate(cut_keys)}

        for i in range(len(filtered_data) - 1):
            current_value = tuple(sorted(filtered_data[i]))
            next_value = tuple(sorted(filtered_data[i + 1]))
            weight = 1  # placeholder weight

            row = value_to_index[current_value]
            col = value_to_index[next_value]

            transition_matrix[row, col] += weight * 100000

        row_sums = transition_matrix.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # prevent NaN rows

        x_normed = transition_matrix / row_sums

        ### Check the nature of your system: 
        # you cant have sum of any row be 0 or else it means you go into a state but never come out of it, or vice versa, which isnt possible

        # this is row with 1 at the diagonal position and 0 everywhere else
        absorbing_states = np.where(np.all(x_normed == np.eye(x_normed.shape[0]), axis=1))[0]

        if len(absorbing_states) > 0:
            print(f"Found {len(absorbing_states)} absorbing states at indices: {absorbing_states}")
        else:
            print("No absorbing states found.")

        # print("calcaulting eigenvalues")
        eigenvalues, eigenvectors = eig(x_normed.T)
        # equilibrium_eigenvector_index=[i for i, e in enumerate(eigenvalues) if np.isclose(e, 1.00)]
        equilibrium_eigenvector_index = [np.argmin(np.abs(eigenvalues - 1.0))]
        equilibrium_eigenvector=(eigenvectors[:,equilibrium_eigenvector_index])
        P_eq = np.abs(np.real(equilibrium_eigenvector / np.sum(equilibrium_eigenvector)))


        cut_keys = list(cut_frequency.keys())
        p_eq_keys = {cut_keys[i]: P_eq[i] for i in range(len(cut_keys))}

        #### 
        P_eq_flat = P_eq.flatten()
        equilibrium_matrix = np.tile(P_eq_flat, (x_normed.shape[0], 1)) # matrix where each row is a copy of the equilibrium distribution vector P_eq_flat
        kd=list(map(lambda x: (P_eq.sum()-x)/x, P_eq))

        kd_keys = {cut_keys[i]: kd[i] for i in range(len(cut_keys))}

        # populating the dictionary
        matrix=transition_matrix

        ls_letters=cut_frequency

        coord_dict = {}
        rows = ls_letters
        cols = ls_letters

        dictionary_transitions={}
        for idx_i,i in enumerate(rows):
            for idx_j,j in enumerate(cols):
                dictionary_transitions[i,j]=matrix[idx_i][idx_j]

        dictionary_transitions_sorted={k: v for k, v in sorted(dictionary_transitions.items(), key=lambda item: item[1],reverse=True)}
        # --> run until here for the hopping v gliding mechanism
        merged_data = defaultdict(int)

        for (key, value) in dictionary_transitions_sorted.items():
            key_tuple = tuple(key)
            merged_data[key_tuple] += value
        merged_output = [(list(key), value) for key, value in merged_data.items()]

        filtered_merged_output=(list(filter(lambda x: x[1] != 0, merged_output)))

        #### Calcualting whether or not bound
        traj = md.load(pdb)
        box_lengths = traj.unitcell_lengths  
        v = np.prod(box_lengths, axis=1)  # volume per frame in nm^3
        volume_L = v.mean() * 1e-24       # average volume in liters

        unbound_key_pop_equilibrium =tuple([-1] * uplet_type)

        unbound_pop_equilibrium = p_eq_keys[unbound_key_pop_equilibrium].item()  # get scalar from np.array
        bound_pop_equilibrium = sum(v.item() for k, v in p_eq_keys.items() if k != unbound_key_pop_equilibrium)

        # print("Unbound from ts:", unbound_pop_equilibrium)
        # print("Bound from ts:", bound_pop_equilibrium)

        Kd_pop_equilibrium = (unbound_pop_equilibrium / bound_pop_equilibrium) * (1 / (scipy.constants.N_A * volume_L))
        print(f"Equilibrium solved from ts is {Kd_pop_equilibrium}")
        # Using the frames as means of calculating kd
        # Load and normalize weights

        bound_mask = np.array(number_contact) >= uplet_type
        unbound_mask = ~bound_mask

        w_bound = w_eq[bound_mask].sum()
        w_unbound = w_eq[unbound_mask].sum()

        Kd_pop_weighted = (w_unbound / w_bound) * (1 / (scipy.constants.N_A * volume_L))
        # print(f"Weighted Kd from populations: {Kd_pop_weighted:.3e} M")

        frame_time_ns = 0.02
        frame_time_s = frame_time_ns * 1e-9

        T_bound_s = (w_eq[bound_mask].sum()) * frame_time_s
        T_unbound_s = (w_eq[unbound_mask].sum()) * frame_time_s

        n_binding = sum((number_contact[i-1] == uplet_type and number_contact[i] >= uplet_type) for i in range(1, len(number_contact)))
        n_unbinding = sum((number_contact[i-1] != uplet_type and number_contact[i] < uplet_type) for i in range(1, len(number_contact)))

        ligand_concentration = 1e-3  # mol/L

        k_off = n_unbinding / T_bound_s
        k_on = n_binding / (T_unbound_s * ligand_concentration)

        Kd_kinetic_weighted = k_off / k_on
        # print(f"Weighted Kd from kinetics: {Kd_kinetic_weighted:.3e} M")

        print(f"Finished {protein_name} successfully!\n")

    except Exception as e:
        print(f"Error processing {protein_name}: {e}")
        continue


from community import community_louvain
from collections import Counter
from collections import defaultdict
import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
import mdtraj as md
from sklearn import preprocessing
from sklearn.preprocessing import MinMaxScaler
import seaborn as sns
import pickle
import networkx as nx
import random
import matplotlib.pyplot as plt
import random
import pickle
import smplotlib
import numpy as np
import os

#### Re-transforming these clusters into trajectory snippets
def convert_to_tcl(residue_dict, cluster_id):
    if cluster_id not in residue_dict:
        return ""
    residue_lists = residue_dict[cluster_id]
    tcl_res_list = "set res_list {\n"
    for res_tuple in residue_lists:
        tcl_res_list += f"    {{{' '.join(map(str, res_tuple))}}}\n"
    tcl_res_list += "}\n"

    return tcl_res_list


alpha_full_length = "alpha_syn_full_peptide"
medin_cm14 = 'medin_cm14'
medin_cm8 = 'medin_cm8'
abeta = "abeta"
medin_urea = 'medin_urea'
uplet_type=5

protein_name = abeta

if protein_name == 'abeta':
    pdb = '/Users/adelielouet/Documents/science/AB_G5_original_simu_analysis/trajectories/Gabis_paper/template_G5.pdb'
    xtc = '/Users/adelielouet/Documents/science/AB_G5_original_simu_analysis/trajectories/Gabis_paper/traj_all-skip-0-noW_G5.xtc'
    # w_file = "/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta/weights/weights_corr_G5"
    # distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta/distances/d_24_t_com_avg.pkl', 'rb'))
    # distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/abeta/distances/d_24_t_com_avg.pkl', 'rb'))
    # w_eq = np.loadtxt(w_file)
    resname_ligand = "liga"
    special_char = 'A\u03B2'
    w_com = 0.8
    w_closest = 0.2

elif protein_name == 'medin_cm14':
    pdb = '/Users/adelielouet/Documents/science/medin/cm14/simu/prepared_system/CM14/prepared_tpr_files_v1/template_complex.gro'
    xtc = '/Users/adelielouet/Documents/science/medin/cm14/simu/prepared_system/CM14/prepared_tpr_files_v1/cat_trjcat.xtc'
    resname_ligand = "1UNL"
    special_char = 'Medin_cm14'
    w_com = 0.45
    w_closest = 0.55
    w_eq=np.array([1.0] *len(distances_com))   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0

elif protein_name == 'medin_cm8':
    pdb = '/Users/adelielouet/Documents/science/medin/cm8/cm8_post_run/protein_ligand.gro'
    xtc = '/Users/adelielouet/Documents/science/medin/cm8/cm8_post_run/concatenate_cm8_skip_10.xtc'
    w_file='/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_cm8/weights/COLVAR_REWEIGHT'
    distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_cm8/distances/d_24_t_com_avg.pkl', 'rb'))
    distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_cm8/distances/d_24_t_closest.pkl', 'rb'))
    resname_ligand = "1UNL"
    special_char = 'Medin_cm8'
    w_com = 0.45
    w_closest = 0.55
    w_eq=np.array([1.0] *len(distances_com))

elif protein_name == 'alpha_syn_full_peptide':
    pdb = '/Users/adelielouet/Documents/science/alpha_syn/DE_SHAW/system.pdb'
    # xtc = ???  # xtc path appears to be missing for alpha
    resname_ligand = "*"
    special_char = "\u03B1-syn-full"
    w_com = 0.1
    w_closest = 0.9
    w_eq=[1.0] *len(distances_com)   # this is for systems that don't need to be reweighted, assign each frame a weight of 1.0

elif protein_name == 'medin_urea':
    pdb = '/Users/adelielouet/Documents/science/medin/urea/urea_post_run/protein_urea.gro'
    xtc = '/Users/adelielouet/Documents/science/medin/urea/urea_post_run/concatenate_urea.xtc'
    resname_ligand = "1UNL"
    distances_com = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_urea/distances/d_24_t_com_avg.pkl', 'rb'))
    distances_closest = pickle.load(open(f'/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_urea/distances/d_24_t_closest.pkl', 'rb'))
    w_file='/Users/adelielouet/Documents/science/dd_proj/msm_full_model_paper/pickled_files/medin_cm8/weights/COLVAR_REWEIGHT'
    special_char = 'Medin_urea'
    w_com = 0.45
    w_closest = 0.55
    w_eq=np.array([1.0] *len(distances_com))

else:
    pass

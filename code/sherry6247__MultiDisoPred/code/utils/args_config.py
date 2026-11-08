'''
Date: 2025-07-18 17:16:48
LastEditors: Sicen Liu
LastEditTime: 2026-03-02 22:05:11
FilePath: /liusicen/methods/MoSE/utils/args_config.py
Description:  
Copyright: © 2025 liusicen@smbu.edu.cn. All rights reserved.
'''

import os
import numpy as np
import torch
from tqdm import tqdm, trange


DATASET_PATH = "/home/suhanghua/MoSE/Datasets/datasets/"
FEATURE_PATH = "/home/suhanghua/MoSE/Data/"
CAID3_PATH = "/home/suhanghua/MoSE/Datasets/CAID3-analysis/"
CAID2_PATH = "/home/suhanghua/MoSE/Datasets/CAID2_datasets/"

import argparse
def seed(s):
    if s.isdigit():
        s = int(s)
        if 0 <= s <= 9999:
            return s
        else:
            raise argparse.ArgumentTypeError("Seed must be between 0 and 9999")
    elif s == "random":
        return np.random.randint(0, 9999)
    else:
        raise ValueError("Invalid seed value")

def parse_args():
    parser = argparse.ArgumentParser(description='Hybirds')
    
    # Task
    parser.add_argument("--model_name", type=str, default="MoSE")
    parser.add_argument("--pssm", type=int, default=1, help="Utilize pssm feature")
    parser.add_argument("--drBERT", type=int, default=1, help="Utilize DrBERT feature")
    parser.add_argument("--esm2", type=int, default=1, help="Utilize esm2 feature")
    parser.add_argument("--active_experts", type=str, default="uni1,uni2,uni3,syn,red",
                        help="Comma-separated list of active experts for ablation. Default: all 5 experts.")
    parser.add_argument("--no_interaction_loss", type=int, default=0, help="Disable interaction loss (1=yes)")
    parser.add_argument("--no_synergy_loss", type=int, default=0, help="Disable synergy loss (1=yes)")
    parser.add_argument("--no_redundancy_loss", type=int, default=0, help="Disable redundancy loss (1=yes)")
    
    # Dataset parameters
    parser.add_argument("--train_dataset_path", type=str, default="/home/suhanghua/MoSE/Datasets/datasets/training_dataset_25/DM4845_nr25_from_CAID3_with_labels.fasta")
    parser.add_argument("--test_dataset", type=str, default="CAID3_disorder_pdb", help="[DISORDER723_test,CASP,SL329]")
    parser.add_argument("--train_dataset", type=str, default= "DM4845_nr25_from_CAID3")
    parser.add_argument("--evo_path", type=str, default="/home/suhanghua/MoSE/Datasets/Evo_features")
    parser.add_argument("--preT_feature_path", type=str, default="/home/suhanghua/MoSE/Datasets/PreT_features")
    
    # Training parameters
    parser.add_argument("--batch_size", type=int, default=4)
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--device", type=str, default='0')
    parser.add_argument("--lr", type=float, default=0.001)
    parser.add_argument("--weight_decay", type=float, default=0.01)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--debug", type=bool, default=False)
    parser.add_argument("--seed", type=seed, default='random', help="random seed")
    
    # Architecture parameters
    parser.add_argument("--hidden_size", type=int, default=128)
    parser.add_argument("--heads", type=int, default=2)
    parser.add_argument("--num_layers", type=int, default=2)
    parser.add_argument("--interact_init", type=float, default=0.1, help="Whether to initialize the interaction module")
    parser.add_argument("--proj_dim", type=int, default=128)
    parser.add_argument("--triplet_margin", type=float, default=0.2)
    parser.add_argument("--rw_hidden", type=int, default=256)
    parser.add_argument("--rw_depth", type=int, default=2)
    parser.add_argument("--rw_tau", type=float, default=1.0)
    parser.add_argument("--rand_scale", type=float, default=1.0)
    
    # Evaluation parameters
    parser.add_argument("--eval_every", type=int, default=100)
    
    # Dropout parameters
    parser.add_argument("--dropout", type=float, default=0.3)

    # Save parameters
    parser.add_argument("--save_result", type=str, default="/home/suhanghua/MoSE/saved_models")
    parser.add_argument("--saved_embedding", type=int, default=1000)
    parser.add_argument("--saved_model_path", type=str, default="/home/suhanghua/MoSE")

    train_data_dict = {
        "Train_723_bc25" : "/home/suhanghua/MoSE/Datasets/datasets/training_dataset_25/Train_723_bc25.fasta",
        "Train_494_bc25" : "/home/suhanghua/MoSE/Datasets/datasets/training_dataset_25/Train_494_bc25.fasta",
        "Train_329_bc25" : "/home/suhanghua/MoSE/Datasets/datasets/training_dataset_25/Train_329_bc25.fasta",
        "Train_casp_bc25" : "/home/suhanghua/MoSE/Datasets/datasets/training_dataset_25/Train_casp_bc25.fasta",
        "Train_CAID3_pdb_bc25" : "/home/suhanghua/MoSE/Datasets/datasets/training_dataset_25/Train_update_caid3_pdb_bc25.fasta",
        "Train_CAID2_pdb_bc25" : "/home/suhanghua/MoSE/Datasets/datasets/training_dataset_25/Train_CAID2_pdb_bc25.fasta",
        "Train_CAID2_nox_bc25" : "/home/suhanghua/MoSE/Datasets/datasets/training_dataset_25/Train_CAID2_nox_bc25.fasta",
        "DM4229_504_training": "/home/suhanghua/MoSE/Data/datasets/DM4229_504_training.fasta",
        "DM4845_training": "/home/suhanghua/MoSE/Data/datasets/DM4845_training.fasta",
        "DM4229_training": "/home/suhanghua/MoSE/Data/datasets/DM4229_training.fasta",
        "lipid_training_converted": "/home/suhanghua/MoSE/lipid_txt/converted/lipid_training_converted.fasta",
        "CAID2_binding": "/home/suhanghua/MoSE/Data/datasets/CAID2_binding.fasta",
        "CAID3_binding": "/home/suhanghua/MoSE/Data/datasets/CAID3_binding.fasta",
        "lipid_training": "/home/suhanghua/MoSE/Data/datasets/lipid_training.fasta",
        "DM4845_nr25_from_CAID3": "/home/suhanghua/MoSE/Datasets/datasets/training_dataset_25/DM4845_nr25_from_CAID3_with_labels.fasta",
        "DM4845_training_nr25_DISORDER723_test": "/home/suhanghua/MoSE/Data/datasets/DM4845_training_nr25_DISORDER723_test.fasta",
        "DM4845_training_nr25_CAID2_disorder_pdb": "/home/suhanghua/MoSE/Data/datasets/DM4845_training_nr25_CAID2_disorder_pdb.fasta",
        "DM4845_training_nr25_CAID3_disorder_pdb": "/home/suhanghua/MoSE/Data/datasets/DM4845_training_nr25_CAID3_disorder_pdb.fasta",
        "DM4845_training_nr25_MXD494": "/home/suhanghua/MoSE/Data/datasets/DM4845_training_nr25_MXD494_fixed.fasta",
        "lipid_TR552_processed_merged": "/home/suhanghua/MoSE/Data/datasets/lipid_TR552_processed_merged.fasta",
        "TR552_idr_binding": "/home/suhanghua/MoSE/Data/datasets/TR552_idr_binding.fasta",
        "TR552": "/home/suhanghua/MoSE/Data/datasets/TR552.fasta",
        # 全量 TR552 部署训练用：fasta 同上，但特征目录名为 TR552_processed
        "TR552_processed": "/home/suhanghua/MoSE/Data/datasets/TR552.fasta",
    } 
    
    test_data_dict = {
        "DISORDER723_test": DATASET_PATH + "DISORDER723_test.fasta",
        "DisProt832_test": DATASET_PATH +  "DisProt832_test.fasta",
        "S1_test": DATASET_PATH + "S1_test.fasta",
        "CASP": DATASET_PATH + "CASP.fasta",
        "Disprot504": DATASET_PATH + "Disprot504.fasta",
        "SL329": DATASET_PATH + "SL329.fasta",
        "MXD494": DATASET_PATH + "MXD494.fasta",
        "CAID3_disorder_nox" : CAID3_PATH + "Datasets_updated/disorder_nox.fasta",
        "CAID3_disorder_pdb" : CAID3_PATH + "Datasets_updated/disorder_pdb.fasta",
        "CAID2_disorder_nox" : CAID2_PATH + "disorder_nox.fasta",
        "CAID2_disorder_pdb" : CAID2_PATH + "disorder_pdb.fasta",
        "lipid_test_converted": "/home/suhanghua/MoSE/lipid_txt/converted/lipid_test_converted.fasta",
        "lipid_validation_converted": "/home/suhanghua/MoSE/lipid_txt/converted/lipid_validation_converted.fasta",
        "CAID2_binding": "/home/suhanghua/MoSE/Data/datasets/CAID2_binding.fasta",
        "CAID3_binding": "/home/suhanghua/MoSE/Data/datasets/CAID3_binding.fasta",
        "lipid_test": "/home/suhanghua/MoSE/Data/datasets/lipid_test.fasta",
        "CAID3_disorder_pdb": DATASET_PATH + "disorder_pdb_CAID3.fasta",
        "binding_idr": "/home/suhanghua/MoSE/Data/datasets/binding_idr.fasta",
    }

    args = parser.parse_args()

    args.train_dataset_path = train_data_dict[args.train_dataset]
    
    test_dataset_path = test_data_dict[args.test_dataset]
    args.test_dataset_path = test_dataset_path
    
    return args
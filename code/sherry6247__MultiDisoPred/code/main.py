'''
Date: 2025-07-18 17:05:41
LastEditors: Sicen Liu
LastEditTime: 2026-11-01 14:43:41
FilePath: /liusicen/mygithub/MultiDisoPred/main.py
Description:  This is the main entry point for training the Prompt Pool IDP model.
'''

import os
import sys
import sys, getopt
import torch
import numpy as np
import re
import os
import random
import gc
import torch.nn as nn
from functools import partial

from src.model import MultiDisoPred, SimpleConcatMLP, SimpleAttentionFusion
from torch.optim import AdamW
from utils.args_config import parse_args
from utils.data_processing import load_dataset, proteinDataset, kfold_split, collate_fn_batch, CAID3_collate_fn_batch, collect_span_embeddings_from_dataset
from torch.utils.data import DataLoader, RandomSampler, SequentialSampler
from src.Trainer import Trainer
from torch.utils.tensorboard import SummaryWriter
# import wandb

# os.environ["WANDB_MODE"] = "offline"

def setup_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

def get_freest_gpu():
    import subprocess
    try:
        result = subprocess.run(
            ['nvidia-smi', '--query-gpu=memory.free', '--format=csv,nounits,noheader'],
            capture_output=True, text=True, check=True
        )
        memory_free = [int(x) for x in result.stdout.strip().split('\n')]
        return memory_free.index(max(memory_free))
    except Exception:
        return 0

def main():
    
    args = parse_args()

    if torch.cuda.is_available():
        torch.cuda.manual_seed(args.seed)
        use_cuda = True
        if torch.cuda.device_count() > 1:
            # DataParallel 要求主设备必须是 cuda:0
            device = torch.device('cuda:0')
            print(f"Running on CUDA:0 (DataParallel with {torch.cuda.device_count()} GPUs)")
        else:
            # 当 CUDA_VISIBLE_DEVICES 限制了可见 GPU 时，直接用 cuda:0
            if torch.cuda.device_count() == 1:
                device = torch.device('cuda:0')
            else:
                best_gpu = get_freest_gpu()
                device = torch.device(f'cuda:{best_gpu}')
            print(f"Running on CUDA:{device}")
    else:
        use_cuda = False
        device = torch.device('cpu')
        print("Running on CPU")

    setup_seed(args.seed)

    args.device = device
    
    criterion = {
        "MultiDisoPred" : nn.BCEWithLogitsLoss(reduction='none'),
    }
    args.criterion = criterion[args.model_name]

    feature_dim_dict = {
        "pssm":20,
        "esm2":1280,     # 完整ESM2特征维度
        "drBERT":768,    # 完整DR-BERT特征维度
    }
    feature_dims = []
    feature_indices = []
    if args.pssm:
        feature_dims.append(feature_dim_dict['pssm'])
        feature_indices.append(0)
    if args.esm2:
        feature_dims.append(feature_dim_dict['esm2'])
        feature_indices.append(1)
    if args.drBERT:
        feature_dims.append(feature_dim_dict['drBERT'])
        feature_indices.append(2)

    ####################################################################
    ## Load the dataset
    #####################################################################
    print("Loading training dataset...")
    train_dataset = load_dataset(args, args.train_dataset_path, dataset_type=args.train_dataset)
    print("Loading test dataset...")
    test_dataset = load_dataset(args, args.test_dataset_path, dataset_type=args.test_dataset)
    
    # 及时释放内存
    gc.collect()
    
    collate_fn_with_params = partial(collate_fn_batch, args=args)
    
    # K-Fold: return dataset lists
    train_Datas, val_Datas = kfold_split(train_dataset, seed=args.seed)
    
    # 及时释放内存
    del train_dataset
    gc.collect()
    
    for k, (train_Dataset, val_Dataset) in  enumerate(zip(train_Datas, val_Datas)):
        #############################Train#######################################
        train_Dataset = proteinDataset(train_Dataset)
        train_sampler = RandomSampler(train_Dataset)
        train_dataloader = DataLoader(train_Dataset,
                                      sampler=train_sampler,
                                      batch_size=args.batch_size,
                                      collate_fn=collate_fn_with_params)
        val_Dataset = proteinDataset(val_Dataset)
        val_dataloader = DataLoader(val_Dataset,
                                    sampler=SequentialSampler(val_Dataset),
                                    batch_size=1,
                                    collate_fn=collate_fn_with_params)
        
        #############################Test#######################################
        testDataset = proteinDataset(test_dataset)
        if args.test_dataset in ['CASP','SL329','CAID3_disorder_pdb', 'CAID3_disorder_nox', 'CAID2_disorder_nox', 'CAID2_disorder_pdb', 'binding_idr']:
            test_collate_fn_with_params = partial(CAID3_collate_fn_batch, args=args)
        else:
            test_collate_fn_with_params = partial(collate_fn_batch, args=args)
        test_dataloader = DataLoader(testDataset,
                                    sampler=SequentialSampler(testDataset),
                                    batch_size=1,
                                    collate_fn=test_collate_fn_with_params)
        

        if args.model_name == 'MultiDisoPred':
            active_expert_names = [e.strip() for e in args.active_experts.split(',')]
            model = MultiDisoPred(feature_dims=feature_dims, feature_indices=feature_indices, d_model=args.hidden_size, dropout=args.dropout, proj_dim=args.proj_dim,
            # lambda_int=args.lambda_int,
            triplet_margin=args.triplet_margin,
            rw_hidden=args.rw_hidden,
            rw_depth=args.rw_depth,
            rw_tau=args.rw_tau,
            rand_scale=args.rand_scale,
            active_expert_names=active_expert_names,
            no_synergy_loss=bool(args.no_synergy_loss),
            no_redundancy_loss=bool(args.no_redundancy_loss))
        elif args.model_name == 'SimpleConcatMLP':
            model = SimpleConcatMLP(feature_dims=feature_dims, hidden_size=args.hidden_size, dropout=args.dropout, feature_indices=feature_indices)
        elif args.model_name == 'SimpleAttentionFusion':
            model = SimpleAttentionFusion(feature_dims=feature_dims, d_model=args.hidden_size, dropout=args.dropout, proj_dim=args.proj_dim, feature_indices=feature_indices)
        

        if torch.cuda.device_count() > 1:
            print("Let's use", torch.cuda.device_count(), "GPUs!")
            # dim = 0 [30, xxx] -> [10, ...], [10, ...], [10, ...] on 3 GPUs
            model = nn.DataParallel(model)
            model.to('cuda')
        optimizer = AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay) 

        args.optimizer = optimizer
        args.k = k
        
        trainer = Trainer(args, model, train_dataloader, val_dataloader, test_dataloader)
        trainer.train()
    


if __name__ == "__main__":
    main()
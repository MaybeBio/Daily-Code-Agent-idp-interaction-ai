'''
Author: Sicen Liu
Date: 2025-05-16 15:42:32
LastEditTime: 2026-11-01 14:42:46
FilePath: /liusicen/mygithub/MultiDisoPred/src/Trainer.py
Description: 

Copyright (c) 2025 by ${liusicen_cs@outlook.com}, All Rights Reserved. 
'''

import os
import torch
import shutil
import numpy as np
from tqdm import tqdm, trange
from torch.utils.data import DataLoader, TensorDataset
# import wandb
from collections import defaultdict

from src.evaluate import eval, predict

import json
import pickle
from collections import OrderedDict

from torch.utils.tensorboard import SummaryWriter


"""adjust_learning_rate"""
def lr_poly(base_lr, iter, max_iter, power):
    # ratio_length = 1 - (float(current_length) / 30) 
    # iter = iter + ratio_length
    # iter = iter + current_length
    if iter > max_iter:
        iter = iter % max_iter
    return base_lr * ((1 - float(iter) / max_iter) ** (power))#+ (float(current_length) / 30) ** (power))
    # return base_lr * (((1 - float(iter) / max_iter) ** (power))+ 0.1*((1 - (float(current_length) / 30) ** (power))))

def adjust_learning_rate(optimizer, i_iter, args):
    lr = lr_poly(args.lr, i_iter, 100000, 0.9)
    optimizer.param_groups[0]['lr'] = np.min(np.around(lr,8))
    if len(optimizer.param_groups) > 1:
        optimizer.param_groups[1]['lr'] = lr * 10
    return lr

class Trainer:

    def __init__(
            self, args, 
            model, 
            train_dataloader:DataLoader, 
            val_dataloader:DataLoader, 
            test_dataloader
    ):
        self.args = args
        self.model = model
        self.train_dataloader = train_dataloader
        self.val_dataloader = val_dataloader
        self.test_dataloader = test_dataloader
        self.optimizer = args.optimizer
    
    def train(self):
        print("Start Training...")
        if not self.args.debug:
            writer = SummaryWriter(log_dir="runs/{}/{}_{}/k_{}-seed_{}-pssm_{}-esm2_{}-drbert_{}-valmet_{}".format(self.args.model_name,self.args.train_dataset, self.args.test_dataset, self.args.k, self.args.seed, self.args.pssm, self.args.esm2, self.args.drBERT, getattr(self.args, 'val_metric', 'AUC')))
        patients, t_patients = 0, 0
        model = self.model
        global_step = 0
        nb_tr_setp = 0
        tr_loss = 0
        val_loss = 0
        dev_auc_best, dev_mcc_best, dev_bacc_best, val_test_auc_best = 0, 0, 0, 0
        val_metric = getattr(self.args, 'val_metric', 'AUC')
        higher_better = (val_metric != 'BCE_loss')
        dev_best = float('-inf') if higher_better else float('inf')
        test_auc_best = 0
        best_results = defaultdict(list)
        cu_iter = 0
        # for epoch in trange(int(self.args.epochs), desc="Epoch:", ncols=100):
        def compute_pos_weight_from_loader(loader, device):
            pos = 0
            neg = 0
            for batch in loader:
                *_, input_mask, labels = batch
                m = input_mask.to(device).bool()
                y = labels.to(device)
                yy = y[m]
                pos += int((yy == 1).sum().item())
                neg += int((yy == 0).sum().item())
            pos = max(pos, 1)
            neg = max(neg, 1)
            return neg / pos
        pos_weight = compute_pos_weight_from_loader(self.train_dataloader, self.args.device)

        for epoch in range(int(self.args.epochs)):
            print("\n--------K: {} epoch: {} seed: {} model_name:{} train_data:{}/{} test_data:{} --------".format(self.args.k, epoch, self.args.seed, self.args.model_name, self.args.train_dataset, len(self.train_dataloader.dataset), self.args.test_dataset))
            tr_loss = 0
            nb_tr_examples, nb_tr_steps = 0, 0
            prog_iter = tqdm(self.train_dataloader, desc='Training', ncols=100, colour='blue')
            model = model.to(self.args.device)
            model.train()
            for _, batch in enumerate(prog_iter):
                self.optimizer.zero_grad()

                id, seq, ori_labels, input_feature1, input_feature2, input_feature3, input_mask, labels = batch
                # Move tensors to GPU
                input_feature1 = input_feature1.to(self.args.device)
                input_feature2 = input_feature2.to(self.args.device)
                input_feature3 = input_feature3.to(self.args.device)
                input_mask = input_mask.to(self.args.device)
                labels = labels.to(self.args.device)

                # cu_iter +=1
                # adjust_learning_rate(self.optimizer, cu_iter, self.args)
                if self.args.model_name == 'MultiDisoPred':
                    logits, aux = model(input_feature1, input_feature2, input_feature3, input_mask)
                    weights = aux['weights']
                    expert_logits = aux['expert_logits']
                    interaction_loss = aux['interaction_loss']
                    interaction_loss_terms = aux['interaction_loss_terms']
                    self.args.criterion = torch.nn.BCEWithLogitsLoss(reduction='none', pos_weight=torch.tensor(pos_weight).to(self.args.device))
                elif self.args.model_name in ['SimpleConcatMLP', 'SimpleAttentionFusion']:
                    logits, _ = model(input_feature1, input_feature2, input_feature3, input_mask)
                    self.args.criterion = torch.nn.BCEWithLogitsLoss(reduction='none', pos_weight=torch.tensor(pos_weight).to(self.args.device))

                loss = self.args.criterion(logits, labels)
                denom = input_mask.sum().clamp_min(1.0)
                loss = (loss * input_mask).sum() / denom

                if self.args.model_name == 'MultiDisoPred' and not getattr(self.args, 'no_interaction_loss', False):
                    loss = loss + self.args.interact_init *  interaction_loss
                
                #清零梯度
                self.optimizer.zero_grad()
                
                loss.backward()
                self.optimizer.step()

                # train_results = cal_batch(self.args, labels, logits, residue_masks, sequences_masks)

                tr_loss += loss.item()
                nb_tr_examples += 1
                nb_tr_steps += 1

                # Display loss
                # if self.is_local_0 == 0:
                prog_iter.set_postfix({"train_loss" :'%.4f' % (tr_loss / nb_tr_steps), "lr" : "%.6f" % (self.optimizer.param_groups[0]['lr'])})

                torch.cuda.empty_cache()
            if not self.args.debug:
                writer.add_scalar("Loss/train", tr_loss / nb_tr_steps, epoch)
        
            # Validation 
            print("=================VALIDATION =================")
            val_results = eval(self.args, model, self.val_dataloader)
            if not self.args.debug:
                
                for vk, vv in val_results.items():
                    writer.add_scalar("val/{}".format(vk), vv, epoch)
            
            feature_params = {"esm2":self.args.esm2, "drBERT":self.args.drBERT, "pssm":self.args.pssm}
            feature_name = "-".join([k for k, v in feature_params.items() if v == 1])
            if getattr(self.args, 'no_interaction_loss', False):
                feature_name += "-no-interaction"
            if getattr(self.args, 'no_synergy_loss', False):
                feature_name += "-no-synergy"
            if getattr(self.args, 'no_redundancy_loss', False):
                feature_name += "-no-redundancy"
            active_experts = getattr(self.args, 'active_experts', '')
            if active_experts and active_experts.replace(' ', '') != 'uni1,uni2,uni3,syn,red':
                feature_name += '-' + active_experts.replace(',', '-').replace(' ', '')
            feature_name = [feature_name]

            val_met = val_results[val_metric]

            is_better = (val_met > dev_best) if higher_better else (val_met < dev_best)
            if is_better:
                patients = 0
                dev_best = val_met
                save_path = os.path.join(self.args.save_result, "{}".format(self.args.model_name), "{}".format(self.args.test_dataset), "{}".format(self.args.train_dataset), "SEED_{}".format(self.args.seed), "{}".format("-".join(feature_name)), "lr_{}_in_{}_pj_{}tm_{}_rwh_{}_rwd_{}_rwt_{}_rs_{}".format(self.args.lr, self.args.interact_init, self.args.proj_dim, self.args.triplet_margin, self.args.rw_hidden, self.args.rw_depth, self.args.rw_tau, self.args.rand_scale), "val_metric_{}".format(val_metric), "K_{}".format(self.args.k))
                if not os.path.exists(save_path):
                    os.makedirs(save_path)
                torch.save({
                    "epoch" : epoch,
                    "k" : self.args.k,
                    "seed" : self.args.seed,
                    "model_state_dict" : model.state_dict()
                }, os.path.join(save_path,"model.pth"))
                try:
                    val_selected_keys = ['AUC', 'APS-PRAUC', 'F_max', 'MCC', 'Sn', 'Sp', 'F1', 'BCE_loss']
                    val_selected_results = {k: val_results[k] for k in val_selected_keys if k in val_results}
                    val_selected_results['val_metric'] = val_metric
                    val_selected_results['Epoch'] = epoch
                    val_selected_results['K'] = self.args.k
                    val_selected_results['Seed'] = self.args.seed
                    val_selected_results['Model'] = self.args.model_name
                    val_selected_results['Test_Dataset'] = self.args.test_dataset
                    val_selected_results['hidden_size'] = self.args.hidden_size
                    val_selected_results['num_layers'] = self.args.num_layers
                    val_selected_results['learning_rate'] = self.args.lr
                    val_selected_results["train_dataset_len"] = (len(self.train_dataloader.dataset)+len(self.val_dataloader.dataset))
                    val_selected_results['active_experts'] = getattr(self.args, 'active_experts', 'all')
                    
                    with open(os.path.join(save_path, 'best_val_results.json'), 'w') as f:
                        json.dump(val_selected_results, f, indent=4)
                    
                except Exception as e:
                    print(f"Error saving model: {e}")

               
            else:
                patients += 1
                if patients == 10:
                    print("Early stopping at epoch {}.".format(epoch))
                    break
                


        writer.close()
            
            
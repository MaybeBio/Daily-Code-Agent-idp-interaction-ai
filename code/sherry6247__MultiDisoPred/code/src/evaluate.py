'''
Author: Sicen Liu
Date: 2024-12-13 07:37:56
LastEditTime: 2026-11-01 14:37:54
FilePath: /liusicen/mygithub/MultiDisoPred/src/evaluate.py
Description: 

Copyright (c) 2024 by ${liusicen_cs@outlook.com}, All Rights Reserved. 
'''
import pickle
import numpy as np
from sklearn import metrics
import torch
from src.metrics import get_results, predict_get_matrics,predict_get_polt_matrics
import wandb
from collections import defaultdict
from tqdm import tqdm
import os
import json

# evaluate
def eval(args, model, eval_dataloader, type='val', save_path=None):
    model.eval()
    eval_loss = 0
    eval_loss_bce = 0
    # batch_results = defaultdict(list)
    data_ids, data_seqs, data_ori_labels, data_labels, data_logits, data_residue_masks, data_sequences_masks = [], [], [], [], [], [], []
    data_gt_masks = []
    data_gt_seqs = []
    data_gt_labels = []
    val_prog_iter = tqdm(eval_dataloader, desc='Val', ncols=100, colour='green')
    nb_val_steps = 0
    with torch.no_grad():
        for idx, batch in enumerate(val_prog_iter):
            id, seq, ori_labels, input_feature1, input_feature2, input_feature3, input_mask, labels = batch
            # Move tensors to GPU
            input_feature1 = input_feature1.to(args.device)
            input_feature2 = input_feature2.to(args.device)
            input_feature3 = input_feature3.to(args.device)
            input_mask = input_mask.to(args.device)
            labels = labels.to(args.device)
            if args.model_name == 'ResidueI2MoE':
                logits, aux = model(input_feature1, input_feature2, input_feature3, input_mask)
                weights = aux['weights']
                expert_logits = aux['expert_logits']
                interaction_loss = aux['interaction_loss']
                interaction_loss_terms = aux['interaction_loss_terms']
            elif args.model_name in ['SimpleConcatMLP', 'SimpleAttentionFusion']:
                logits, _ = model(input_feature1, input_feature2, input_feature3, input_mask)
            e_loss = args.criterion(logits, labels)
            denom = input_mask.sum().clamp_min(1.0)
            e_loss_bce = (e_loss * input_mask).sum() / denom
            e_loss = (e_loss * input_mask).mean()
            if args.model_name == 'ResidueI2MoE':
                e_loss = e_loss + args.interact_init * interaction_loss
            eval_loss += e_loss.item()
            eval_loss_bce += e_loss_bce.item()
            logits = torch.sigmoid(logits)
            
            nb_val_steps += 1
            val_prog_iter.set_postfix(val_loss='%.4f' % (eval_loss / nb_val_steps))
            # val_prog_iter.set_postfix(val_loss='%.4f' % (e_loss.item()))
            torch.cuda.empty_cache()
            
            # # 计算batch数据的auc, F1, mcc, acc
            # batch_result = cal_batch(args, labels, logits, residue_masks, sequences_masks)
            # for k, v in batch_result.items():
            #     batch_results[k].append(v)
            data_ids.append(id)
            data_seqs.append(seq)
            data_ori_labels.extend(ori_labels)
            data_labels.extend(labels.detach().cpu().numpy())
            data_logits.extend(logits.detach().cpu().numpy())
            data_gt_masks.extend(input_mask.detach().cpu().numpy())
    
    # 在验证集上使用自适应阈值，测试集使用固定阈值
    use_adaptive = (type == 'val')
    reuslts = get_results(args, data_labels, data_logits, data_gt_masks, use_adaptive_threshold=use_adaptive)
    reuslts['BCE_loss'] = float(eval_loss_bce / nb_val_steps)
    if type == 'test' or type == 'val_test':
        # 保存预测的结果
        try:
            result_path = os.path.join(save_path, '{}_predict_result.txt'.format(type))
            with open(result_path, 'w') as file:
                new_results = {}
                for  k, v in reuslts.items():
                    new_results[k] = float(v)
                json.dump(new_results, file)  # 写入字典
                file.write('\n')
                file.write('====='*5) 
                for id, seq, label, logit in zip(data_ids, data_seqs, data_labels, data_logits):
                    file.write(f'\n>{id}:\n')
                    for sq, ll , lt in zip(seq, label, logit):
                        file.write("{}\t{}\t{}\n".format(sq, ll, lt))
        except Exception as e:
            print(f"Error saving results: {e}")

    return reuslts

def preict_get_results(args, ids, seqs, labels, logits, input_mask):
    true_label, pre_label, pre_prob = [], [], []
    for i, logit in enumerate(logits):
        true_label.append(labels[i])
        pre_label.append((logit >= args.threshold).astype(int))
    results, sorted_total_samples = predict_get_matrics(ids, seqs, true_label, pre_label, logits, input_mask)
    return results, sorted_total_samples

# prediction part
def predict(args, model, eval_dataloader):
    model.eval()
    eval_loss = 0
    # batch_results = defaultdict(list)
    data_ids, data_seqs, data_ori_labels, data_labels, data_logits, data_residue_masks, data_sequences_masks = [], [], [], [], [], [], []
    data_gt_masks = []
    data_gt_seqs = []
    data_gt_labels = []
    tri_feature_list = []
    preT_feature_list = []
    fused_feature_list = []
    val_prog_iter = tqdm(eval_dataloader, desc='Predict', ncols=100, colour='green')
    nb_val_steps = 0
    with torch.no_grad():
        for idx, batch in enumerate(val_prog_iter):
            id, seq, ori_labels, input_features, input_mask, labels = batch
            # Move tensors to GPU
            input_features = input_features.to(args.device)
            input_mask = input_mask.to(args.device)
            labels = labels.to(args.device)
            if args.model_name == 'PromptIDP_Att':
                logits, prompt_orth_loss = model(input_features)
            elif args.model_name == 'MultiDisoPred' or args.model_name == 'ConditionalIDRModel':
                logits_dict = model(input_features)
                logits = logits_dict['logits']
            else:
                logits = model(input_features)
            e_loss = args.criterion(logits, labels)
            e_loss = (e_loss * input_mask).mean()
            if args.model_name == 'PromptIDP_Att':
                e_loss = e_loss + args.lambda_orth * prompt_orth_loss
            eval_loss += e_loss.item()
            logits = torch.sigmoid(logits)
            
            nb_val_steps += 1
            val_prog_iter.set_postfix(val_loss='%.4f' % (eval_loss / nb_val_steps))
            # val_prog_iter.set_postfix(val_loss='%.4f' % (e_loss.item()))
            torch.cuda.empty_cache()
            
            data_ids.append(id)
            data_seqs.append(seq)
            data_ori_labels.extend(ori_labels)
            data_labels.extend(labels.detach().cpu().numpy())
            data_logits.extend(logits.detach().cpu().numpy())
            data_gt_masks.extend(input_mask.detach().cpu().numpy())
    
    results, sorted_total_samples = preict_get_results(args, data_ids, data_seqs, data_labels, data_logits, data_gt_masks)
                
        
    return results, sorted_total_samples


def preict_get_polt_results(args, ids, seqs, labels, logits, input_mask, route_ws, route_idx, gates):
    true_label, pre_label, pre_prob = [], [], []
    for i, logit in enumerate(logits):
        true_label.append(labels[i])
        pre_label.append((logit >= args.threshold).astype(int))
    results, sorted_total_samples = predict_get_polt_matrics(ids, seqs, true_label, pre_label, logits, input_mask, route_ws, route_idx, gates)
    return results, sorted_total_samples
def predict_polt(args, model, eval_dataloader):
    model.eval()
    eval_loss = 0
    # batch_results = defaultdict(list)
    data_ids, data_seqs, data_ori_labels, data_labels, data_logits, data_residue_masks, data_sequences_masks = [], [], [], [], [], [], []
    data_gt_masks = []
    data_routes_ws = []
    data_routes_idx = []
    data_gates = []
    val_prog_iter = tqdm(eval_dataloader, desc='Predict', ncols=100, colour='green')
    nb_val_steps = 0
    with torch.no_grad():
        for idx, batch in enumerate(val_prog_iter):
            id, seq, ori_labels, input_features, input_mask, labels = batch
            # Move tensors to GPU
            input_features = input_features.to(args.device)
            input_mask = input_mask.to(args.device)
            labels = labels.to(args.device)
            if args.model_name == 'PromptIDP_Att':
                logits, prompt_orth_loss = model(input_features)
            elif args.model_name == 'MultiDisoPred':
                logits_dict = model(input_features)
                logits = logits_dict['logits']
            else:
                logits = model(input_features)
            e_loss = args.criterion(logits, labels)
            e_loss = (e_loss * input_mask).mean()
            if args.model_name == 'PromptIDP_Att':
                e_loss = e_loss + args.lambda_orth * prompt_orth_loss
            eval_loss += e_loss.item()
            logits = torch.sigmoid(logits)
            
            nb_val_steps += 1
            val_prog_iter.set_postfix(val_loss='%.4f' % (eval_loss / nb_val_steps))
            # val_prog_iter.set_postfix(val_loss='%.4f' % (e_loss.item()))
            torch.cuda.empty_cache()
            
            data_ids.append(id)
            data_seqs.append(seq)
            data_ori_labels.extend(ori_labels)
            data_labels.extend(labels.detach().cpu().numpy())
            data_logits.extend(logits.detach().cpu().numpy())
            data_gt_masks.extend(input_mask.detach().cpu().numpy())
            data_routes_ws.extend(logits_dict['route_weights'].detach().cpu().numpy())
            data_routes_idx.extend(logits_dict['route_idx'].detach().cpu().numpy())
            data_gates.extend(logits_dict['gate'].detach().cpu().numpy())
    
    results, sorted_total_samples = preict_get_polt_results(args, data_ids, data_seqs, data_labels, data_logits, data_gt_masks, data_routes_ws, data_routes_idx, data_gates)
                
        
    return results, sorted_total_samples

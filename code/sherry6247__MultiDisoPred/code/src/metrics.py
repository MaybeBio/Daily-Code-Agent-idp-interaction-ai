'''
Author: Sicen Liu
Date: 2024-09-10 06:52:43
LastEditTime: 2026-01-04 15:36:24
FilePath: /liusicen/methods/Prompt_pool_IDP/src/metrics.py
Description: 

Copyright (c) 2024 by ${liusicen_cs@outlook.com}, All Rights Reserved. 
'''
import numpy as np
from sklearn import metrics
from sklearn.metrics import matthews_corrcoef, balanced_accuracy_score, f1_score, roc_auc_score, auc, precision_recall_curve, roc_curve, average_precision_score, average_precision_score
from math import *
from collections import OrderedDict


def data2numpy(data):
    return data.cpu().detach().numpy()


def get_matrics(true_label, pre_label, pred_logit, mask):
        results = {}
        total_true_label, total_pre_logit = [], []
        TP = 0
        FN = 0
        TN = 0
        FP = 0
        other = 0

        for i in range(len(pre_label)):
            for j in range(len(pre_label[i])):
                if mask[i][j] == 1:
                    if pre_label[i][j] == 0. and true_label[i][j] == 1.:
                        FN = FN + 1
                    elif pre_label[i][j] == 1. and true_label[i][j] == 1.:
                        TP = TP + 1
                    elif pre_label[i][j] == 0. and true_label[i][j] == 0.:
                        TN = TN + 1
                    elif pre_label[i][j] == 1. and true_label[i][j] == 0.:
                        FP = FP + 1
                    elif true_label[i][j] == '2':
                        other = other + 1
                    total_true_label.append(true_label[i][j])
                    total_pre_logit.append(pred_logit[i][j])
                else: continue

        # print('TP =', TP)
        # print('FP =', FP)
        # print('TN =', TN)
        # print('FN =', FN)
        # print('other =', other)
        # print('if the test set is whole ordered, the FPR is adopted!')
        # print('the residue FPR is :', round(FP/(1.0*(FP + TN)), 5))

        if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) != 0:
            Sn = TP / float(TP + FN)  # Sensitivity
            Sp = TN / float(TN + FP)  # Specificity
            BACC = (Sn + Sp) / 2.0      # Balanced Accuracy
            MCC = ((TP * TN) - (FP * FN)) / float(sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
            F1 = (2. * TP ) / (2 * TP + FP + FN)
        else:
            Sn = 0
            Sp = 0
            BACC = 0
            MCC = 0
            F1 = 0
        Sn = round(Sn, 4)
        Sp = round(Sp, 4)
        BACC = round(BACC, 4)
        MCC = round(MCC, 4)
        F1 = round(F1, 4)
        results["Sn"] = Sn  # Sensitivity
        results["Sp"] = Sp  # Specificity
        results["BACC"] = BACC  # Balanced Accuracy
        results["MCC"] = MCC  # Matthews Correlation Coefficient
        results["F1"] = F1

        # 计算AUC，处理边界情况
        try:
            # 检查是否只有一种类别
            unique_labels = np.unique(np.array(total_true_label))
            if len(unique_labels) < 2:
                print(f"Warning: Only one class present in labels ({unique_labels[0] if len(unique_labels) > 0 else 'none'}). AUC set to 0.")
                AUC = 0.0
                APS_PRAUC = 0.0
                F_max = 0.0
            else:
                AUC = roc_auc_score(np.array(total_true_label), np.array(total_pre_logit))
                
                # 计算APS-PRAUC
                precisions, recalls, _ = precision_recall_curve(np.array(total_true_label), np.array(total_pre_logit))
                APS_PRAUC = auc(recalls, precisions)
                
                # 计算F_max
                f1_scores = 2 * (precisions * recalls) / (precisions + recalls + 1e-8)  # 避免除零
                F_max = np.nanmax(f1_scores)  # 忽略NaN值
        except Exception as e:
            print(f"Warning: Metrics calculation failed: {e}")
            AUC = 0.0
            APS_PRAUC = 0.0
            F_max = 0.0
        
        results["AUC"] = AUC
        results["APS-PRAUC"] = APS_PRAUC
        results["F_max"] = F_max
        
        # 打印关键指标
        print("AUC:{}".format(AUC))
        print("APS-PRAUC:{}".format(APS_PRAUC))
        print("F1:{}".format(F1))
        print("F_max:{}".format(F_max))
        print("Sensitivity (Sn):{}".format(Sn))
        print("Specificity (Sp):{}".format(Sp))
        print("Balanced Accuracy (BACC):{}".format(BACC))
        print("Matthews Correlation Coefficient (MCC):{}".format(MCC))
        print("TP:{} FP:{} TN:{} FN:{}".format(TP, FP, TN, FN))
        
        return results

def find_optimal_threshold(true_label, logits, input_mask, threshold_range=None):
    """
    在验证集上寻找最优阈值，使F1分数最大化
    
    Args:
        true_label: 真实标签列表
        logits: 预测概率列表
        input_mask: 掩码列表
        threshold_range: 阈值搜索范围，默认0.01到0.99，步长0.01（更精细的搜索）
    
    Returns:
        optimal_threshold: 最优阈值
        best_f1: 最佳F1分数
        best_results: 最优阈值对应的评估结果
    """
    if threshold_range is None:
        # 对于极端不平衡数据，使用更精细的阈值搜索
        threshold_range = np.arange(0.01, 0.99, 0.01)
    
    best_f1 = 0
    optimal_threshold = 0.5
    best_results = None
    
    for threshold in threshold_range:
        pre_label = []
        for i, logit in enumerate(logits):
            pre_label.append((logit >= threshold).astype(int))
        
        results = get_matrics(true_label, pre_label, logits, input_mask)
        f1 = results.get("F1", 0)
        
        if f1 > best_f1:
            best_f1 = f1
            optimal_threshold = threshold
            best_results = results
    
    print(f"\n自适应阈值搜索结果:")
    print(f"最优阈值: {optimal_threshold:.2f}")
    print(f"最佳F1: {best_f1:.4f}")
    
    return optimal_threshold, best_f1, best_results


def get_results(args, labels, logits, input_mask, use_adaptive_threshold=False):
    """
    计算评估结果
    
    Args:
        args: 参数对象
        labels: 真实标签
        logits: 预测概率
        input_mask: 掩码
        use_adaptive_threshold: 是否使用自适应阈值（仅在验证集上使用）
    
    Returns:
        results: 评估结果字典
    """
    true_label, pre_label = [], []
    
    # 如果使用自适应阈值，先搜索最优阈值
    if use_adaptive_threshold:
        optimal_threshold, _, _ = find_optimal_threshold(labels, logits, input_mask)
        threshold = optimal_threshold
        print(f"使用自适应阈值: {threshold:.2f} (原阈值: {args.threshold})")
    else:
        threshold = args.threshold
    
    for i, logit in enumerate(logits):
        true_label.append(labels[i])
        pre_label.append((logit >= threshold).astype(int))
    results = get_matrics(true_label, pre_label, logits, input_mask)
    return results

def predict_get_matrics(ids, seqs, true_label, pre_label, pred_logit, mask):
    results = {}
    total_samples = {}
    total_true_label, total_pre_logit = [], []
    TP = 0
    FN = 0
    TN = 0
    FP = 0
    other = 0

    for i in range(len(pre_label)):
        for j in range(len(pre_label[i])):
            if mask[i][j] == 1:
                if pre_label[i][j] == 0. and true_label[i][j] == 1.:
                    FN = FN + 1
                elif pre_label[i][j] == 1. and true_label[i][j] == 1.:
                    TP = TP + 1
                elif pre_label[i][j] == 0. and true_label[i][j] == 0.:
                    TN = TN + 1
                elif pre_label[i][j] == 1. and true_label[i][j] == 0.:
                    FP = FP + 1
                elif true_label[i][j] == '2':
                    other = other + 1
                total_true_label.append(true_label[i][j])
                total_pre_logit.append(pred_logit[i][j])
            else: continue

    for i in range(len(pre_label)):
            id = ids[i]
            mask_i = mask[i]
            label_i = true_label[i]
            pred_i = pred_logit[i] #[round(l,4) for l in pred_logit[i]] # 保留小数点后四位
            seq = seqs[i]
            total_samples['{}'.format(id)] = ((seq, label_i, pred_i))
        
    sorted_total_samples = {key: total_samples[key] for key in ids if key in total_samples}
    # print('TP =', TP)
    # print('FP =', FP)
    # print('TN =', TN)
    # print('FN =', FN)
    # print('other =', other)
    # print('if the test set is whole ordered, the FPR is adopted!')
    # print('the residue FPR is :', round(FP/(1.0*(FP + TN)), 5))

    if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) != 0:
        Sn = TP / float(TP + FN)
        Sp = TN / float(TN + FP )
        Bacc = (TP / float(TP + FN) + TN / float(TN + FP)) / 2.0
        MCC = ((TP * TN) - (FP * FN)) / float(sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
        ACC = (TP + TN )/((TP + TN + FP + FN)*1.0)
        F1 = (2. * TP ) / (2 * TP + FP + FN)
        Prec = TP / float(TP + FP)
    else:
        Sn = 0
        Sp = 0
        Bacc = 0
        MCC = 0
        ACC = 0
        F1 = 0
        Prec = 0
    Sn = round(Sn, 4)
    Sp =  round(Sp, 4)
    Bacc =  round(Bacc, 4)
    ACC =  round(ACC, 4)
    MCC =  round(MCC, 4)
    F1 =  round(F1, 4)
    Prec =  round(Prec, 4)
    results["Sn"] = Sn
    results["Sp"] = Sp
    results["Bacc"] = Bacc
    results["ACC"] = ACC
    results["MCC"] = MCC
    results["F1"] = F1
    results["Prec"] = Prec
    
    AP = average_precision_score(np.array(total_true_label), np.array(total_pre_logit))
    results["AP"] = AP
    
    AUC = roc_auc_score(np.array(total_true_label), np.array(total_pre_logit))
    results["AUC"] = AUC
    
    # Fmax
    # 计算不同阈值下的Precision、Recall
    precisions, recalls, thresholds = precision_recall_curve(np.array(total_true_label), np.array(total_pre_logit))

    # 计算每个阈值对应的F1分数，选择最大值
    f1_scores = 2 * (precisions * recalls) / (precisions + recalls + 1e-8)  # 避免除零
    f1_max = np.nanmax(f1_scores)  # 忽略NaN值
    optimal_threshold = thresholds[np.nanargmax(f1_scores)]
    results["F_max"] = f1_max
    results["F_max_threshold"] = optimal_threshold
    
    # 计算 PR-AUC
    pr_auc = auc(recalls, precisions)
    results["APS-PRAUC"] = pr_auc
    
    for k, v in results.items():
        print("{}:{}".format(k, v))
    
    return results, sorted_total_samples
    
def predict_get_matrics1(ids, seqs, true_label, pre_label, pred_logit, mask):
        results = {}
        total_samples = {}
        total_true_label, total_pre_logit = [], []
        TP = 0
        FN = 0
        TN = 0
        FP = 0
        other = 0

        for i in range(len(pre_label)):
            for j in range(len(pre_label[i])):
                if mask[i][j] == 1:
                    if pre_label[i][j] == 0. and true_label[i][j] == 1.:
                        FN = FN + 1
                    elif pre_label[i][j] == 1. and true_label[i][j] == 1.:
                        TP = TP + 1
                    elif pre_label[i][j] == 0. and true_label[i][j] == 0.:
                        TN = TN + 1
                    elif pre_label[i][j] == 1. and true_label[i][j] == 0.:
                        FP = FP + 1
                    elif true_label[i][j] == '2':
                        other = other + 1
                    total_true_label.append(true_label[i][j])
                    total_pre_logit.append(pred_logit[i][j])
                else: continue
        for i in range(len(pre_label)):
            id = ids[i]
            mask_i = mask[i]
            label_i = true_label[i]
            pred_i = pred_logit[i] #[round(l,4) for l in pred_logit[i]] # 保留小数点后四位
            seq = seqs[i]
            total_samples['{}'.format(id)] = ((seq, label_i, pred_i))
        
        sorted_total_samples = {key: total_samples[key] for key in ids if key in total_samples}
        

        # print('TP =', TP)
        # print('FP =', FP)
        # print('TN =', TN)
        # print('FN =', FN)
        # print('other =', other)
        # print('if the test set is whole ordered, the FPR is adopted!')
        # print('the residue FPR is :', round(FP/(1.0*(FP + TN)), 5))

        if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) != 0:
            Sn = TP / float(TP + FN)
            Sp = TN / float(TN + FP )
            Bacc = (TP / float(TP + FN) + TN / float(TN + FP)) / 2.0
            MCC = ((TP * TN) - (FP * FN)) / float(sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
            ACC = (TP + TN )/((TP + TN + FP + FN)*1.0)
            F1 = (2. * TP ) / (2 * TP + FP + FN)
            Prec = TP / float(TP + FP)
        else:
            Sn = 0
            Sp = 0
            Bacc = 0
            MCC = 0
            ACC = 0
            F1 = 0
            Prec = 0
        Sn = round(Sn, 4)
        Sp = round(Sp, 4)
        Bacc = round(Bacc, 4)
        ACC = round(ACC, 4)
        MCC = round(MCC, 4)
        F1 = round(F1, 4)
        Prec = round(Prec, 4)
        results["Sn"] = Sn
        results["Sp"] = Sp
        results["Bacc"] = Bacc
        results["ACC"] = ACC
        results["MCC"] = MCC
        results["F1"] = F1
        results["Prec"] = Prec
        
        AP = average_precision_score(np.array(total_true_label), np.array(total_pre_logit))
        results["AP"] = AP
        
        AUC = roc_auc_score(np.array(total_true_label), np.array(total_pre_logit))
        results["AUC"] = AUC
        for k, v in results.items():
            print("{}:{}".format(k, v))
        
        return results, sorted_total_samples


def predict_get_polt_matrics(ids, seqs, true_label, pre_label, pred_logit, mask, route_ws, route_idx, gates):
    results = {}
    total_samples = {}
    total_true_label, total_pre_logit = [], []
    TP = 0
    FN = 0
    TN = 0
    FP = 0
    other = 0

    for i in range(len(pre_label)):
        for j in range(len(pre_label[i])):
            if mask[i][j] == 1:
                if pre_label[i][j] == 0. and true_label[i][j] == 1.:
                    FN = FN + 1
                elif pre_label[i][j] == 1. and true_label[i][j] == 1.:
                    TP = TP + 1
                elif pre_label[i][j] == 0. and true_label[i][j] == 0.:
                    TN = TN + 1
                elif pre_label[i][j] == 1. and true_label[i][j] == 0.:
                    FP = FP + 1
                elif true_label[i][j] == '2':
                    other = other + 1
                total_true_label.append(true_label[i][j])
                total_pre_logit.append(pred_logit[i][j])
            else: continue

    for i in range(len(pre_label)):
            id = ids[i]
            mask_i = mask[i]
            label_i = true_label[i]
            pred_i = pred_logit[i] #[round(l,4) for l in pred_logit[i]] # 保留小数点后四位
            seq = seqs[i]
            route_ws_i = route_ws[i]
            route_idx_i = route_idx[i]
            gates_i = gates[i]
            total_samples['{}'.format(id)] = ((seq, label_i, pred_i, route_ws_i, route_idx_i, gates_i))
        
    sorted_total_samples = {key: total_samples[key] for key in ids if key in total_samples}
    # print('TP =', TP)
    # print('FP =', FP)
    # print('TN =', TN)
    # print('FN =', FN)
    # print('other =', other)
    # print('if the test set is whole ordered, the FPR is adopted!')
    # print('the residue FPR is :', round(FP/(1.0*(FP + TN)), 5))

    if (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) != 0:
        Sn = TP / float(TP + FN)
        Sp = TN / float(TN + FP )
        Bacc = (TP / float(TP + FN) + TN / float(TN + FP)) / 2.0
        MCC = ((TP * TN) - (FP * FN)) / float(sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
        ACC = (TP + TN )/((TP + TN + FP + FN)*1.0)
        F1 = (2. * TP ) / (2 * TP + FP + FN)
        Prec = TP / float(TP + FP)
    else:
        Sn = 0
        Sp = 0
        Bacc = 0
        MCC = 0
        ACC = 0
        F1 = 0
        Prec = 0
    Sn = round(Sn, 4)
    Sp =  round(Sp, 4)
    Bacc =  round(Bacc, 4)
    ACC =  round(ACC, 4)
    MCC =  round(MCC, 4)
    F1 =  round(F1, 4)
    Prec =  round(Prec, 4)
    results["Sn"] = Sn
    results["Sp"] = Sp
    results["Bacc"] = Bacc
    results["ACC"] = ACC
    results["MCC"] = MCC
    results["F1"] = F1
    results["Prec"] = Prec
    
    AP = average_precision_score(np.array(total_true_label), np.array(total_pre_logit))
    results["AP"] = AP
    
    AUC = roc_auc_score(np.array(total_true_label), np.array(total_pre_logit))
    results["AUC"] = AUC
    
    # Fmax
    # 计算不同阈值下的Precision、Recall
    precisions, recalls, thresholds = precision_recall_curve(np.array(total_true_label), np.array(total_pre_logit))

    # 计算每个阈值对应的F1分数，选择最大值
    f1_scores = 2 * (precisions * recalls) / (precisions + recalls + 1e-8)  # 避免除零
    f1_max = np.nanmax(f1_scores)  # 忽略NaN值
    optimal_threshold = thresholds[np.nanargmax(f1_scores)]
    results["F_max"] = f1_max
    results["F_max_threshold"] = optimal_threshold
    
    # 计算 PR-AUC
    pr_auc = auc(recalls, precisions)
    results["APS-PRAUC"] = pr_auc
    
    for k, v in results.items():
        print("{}:{}".format(k, v))
    
    return results, sorted_total_samples
  
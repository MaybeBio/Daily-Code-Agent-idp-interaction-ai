"""
This script calculates the results for table 1 of DARUMA-binding.
"""
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, precision_recall_curve, f1_score


def read_caid_result_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

        result_dict = {}
        for line in lines:
            if line.startswith("#") or not line.strip():
                continue
            
            if line.startswith(">"):
                ac = line[1:].strip()
                result_dict[ac] = {"aa": [], "probability": [], "binary_prediction": []}
            else:
                _, _aa, _probability, _binary_prediction = line.strip().split("\t")
                result_dict[ac]["aa"].append(_aa)
                result_dict[ac]["probability"].append(float(_probability))
                result_dict[ac]["binary_prediction"].append(int(_binary_prediction))

    return result_dict


def read_caid_dataset_file(file_path):
    with open(file_path, 'r') as f:
        text = f.read()

        dataset_dict = {}
        for data in text.strip(">").split("\n>"):
            lines = data.strip().split("\n")
            ac = lines[0].strip()
            aa = lines[1].strip()
            label = list(map(int, lines[2].strip()))
            dataset_dict[ac] = {"aa": aa, "label": label}
    
    return dataset_dict


def calculate_fmax(all_labels, all_probs):
    precision, recall, thresholds = precision_recall_curve(all_labels, all_probs)
    f_scores = 2 * (precision * recall) / (precision + recall + 1e-8)
    fmax = np.max(f_scores)
    best_threshold = thresholds[np.argmax(f_scores)] if len(thresholds) > 0 else 0.5
    return fmax, best_threshold


def calculate_metrics(probabilities, labels):
    auc = roc_auc_score(labels, probabilities)
    aps = average_precision_score(labels, probabilities)
    fmax = calculate_fmax(labels, probabilities)[0]
    return auc, aps, fmax


result_dict = read_caid_result_file("../results/caid3/caid3_DARUMA-binding.caid")
dataset_dict = read_caid_dataset_file("../../data/test/caid3_binding.fasta")

full_probabilities = []
full_labels = []
for ac in result_dict:
    full_probabilities.extend(result_dict[ac]["probability"])
    full_labels.extend(dataset_dict[ac]["label"])

auc, aps, fmax = calculate_metrics(np.array(full_probabilities), np.array(full_labels))
coverage = len(result_dict) / len(dataset_dict)

print(f"AUC: {auc:.4f}")
print(f"APS: {aps:.4f}")
print(f"Fmax: {fmax:.4f}")
print(f"Coverage: {coverage:.4f}")

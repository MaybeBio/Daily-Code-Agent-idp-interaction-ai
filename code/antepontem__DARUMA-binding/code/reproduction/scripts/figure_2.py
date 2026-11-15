"""
This script calculates the results for figure 2 of DARUMA-binding.
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, precision_recall_curve


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


def make_roc_curve(datasets):
    for name, y_true, y_score in datasets:
        fpr, tpr, _ = roc_curve(y_true, y_score)
        plt.plot(fpr, tpr, label=name)

    plt.plot([0, 1], [0, 1], "--")
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    plt.legend()
    plt.show()


def make_pr_curve(datasets):
    for name, y_true, y_score in datasets:
        precision, recall, _ = precision_recall_curve(y_true, y_score)
        plt.plot(recall, precision, label=name)

    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.legend()
    plt.show()


dataset_dict = read_caid_dataset_file("../../data/test/caid3_binding.fasta")

from pathlib import Path
files = Path("../results/caid3").glob("*.caid")

roc_datasets = []
for file in files:
    print(file.name)
    result_dict = read_caid_result_file(file)

    full_probabilities = []
    full_labels = []
    for ac in result_dict:
        full_probabilities.extend(result_dict[ac]["probability"])
        full_labels.extend(dataset_dict[ac]["label"])
    
    roc_datasets.append((file.stem, full_labels, full_probabilities))

make_roc_curve(roc_datasets)
make_pr_curve(roc_datasets)

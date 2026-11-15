"""
This script calculates the results for figure 3.
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score


def read_occlusion_result_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

        result_dict = {
            "seq": [],
            "y": [],
            "original_probs": [],
            "plddt_occluded_probs": [],
            "rsa_occluded_probs": [],
            "helix_occluded_probs": [],
            "sheet_occluded_probs": [],
            "coil_occluded_probs": [],
        }

        for line in lines[1:]:
            residue, y, original_probs, plddt_occluded_probs, rsa_occluded_probs, helix_occluded_probs, sheet_occluded_probs, coil_occluded_probs = line.strip().split(",")
            result_dict["seq"].append(residue)
            result_dict["y"].append(int(y))
            result_dict["original_probs"].append(float(original_probs))
            result_dict["plddt_occluded_probs"].append(float(plddt_occluded_probs))
            result_dict["rsa_occluded_probs"].append(float(rsa_occluded_probs))
            result_dict["helix_occluded_probs"].append(float(helix_occluded_probs))
            result_dict["sheet_occluded_probs"].append(float(sheet_occluded_probs))
            result_dict["coil_occluded_probs"].append(float(coil_occluded_probs))

    return result_dict


def make_roc_curve(datasets):
    for name, y_true, y_score in datasets:
        fpr, tpr, _ = roc_curve(y_true, y_score)
        plt.plot(fpr, tpr, label=name)

    plt.plot([0, 1], [0, 1], "--")
    plt.xlabel("FPR")
    plt.ylabel("TPR")
    plt.legend()
    plt.show()



from pathlib import Path
files = Path("../results/occlusion_analysis").glob("*.csv")

full_probabilities_dict = {
    "y": [],
    "original_probs": [],
    "plddt_occluded_probs": [],
    "rsa_occluded_probs": [],
    "helix_occluded_probs": [],
    "sheet_occluded_probs": [],
    "coil_occluded_probs": [],
}

for file in files:
    result_dict = read_occlusion_result_file(file)

    for key in full_probabilities_dict.keys():
        full_probabilities_dict[key].extend(result_dict[key])


roc_datasets = []
for key in full_probabilities_dict.keys():
    if key == "y":
        continue
    roc_datasets.append((key, full_probabilities_dict["y"], full_probabilities_dict[key]))
make_roc_curve(roc_datasets)


# AUC
original_auc = roc_auc_score(full_probabilities_dict["y"], full_probabilities_dict["original_probs"])
print(f"original_probs: AUC: {original_auc:.4f}")

for key in full_probabilities_dict.keys():
    if key == "y" or key == "original_probs":
        continue

    auc = roc_auc_score(full_probabilities_dict["y"], full_probabilities_dict[key])
    print(f"{key}: AUC: {auc:.4f}, delta AUC: {original_auc - auc:.4f}")
    
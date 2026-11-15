"""
This script calculates the results for figure 4.
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


def make_violin_plot(data, labels, title=None):
    plt.violinplot(data, showmedians=True)

    means = [np.mean(d) for d in data]
    positions = range(1, len(data) + 1)

    plt.scatter(positions, means, marker="*", s=100, label="Mean")

    plt.xticks(positions, labels)
    plt.ylabel("Value")

    if title is not None:
        plt.title(title)

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


for key in full_probabilities_dict.keys():
    full_probabilities_dict[key] = np.array(full_probabilities_dict[key])


full_delta_probabilities_dict = {
    "plddt": full_probabilities_dict["original_probs"] - full_probabilities_dict["plddt_occluded_probs"],
    "rsa": full_probabilities_dict["original_probs"] - full_probabilities_dict["rsa_occluded_probs"],
    "helix": full_probabilities_dict["original_probs"] - full_probabilities_dict["helix_occluded_probs"],
    "sheet": full_probabilities_dict["original_probs"] - full_probabilities_dict["sheet_occluded_probs"],
    "coil": full_probabilities_dict["original_probs"] - full_probabilities_dict["coil_occluded_probs"],
}

make_violin_plot(
    [
        full_delta_probabilities_dict["plddt"],
        full_delta_probabilities_dict["rsa"],
        full_delta_probabilities_dict["helix"],
        full_delta_probabilities_dict["sheet"],
        full_delta_probabilities_dict["coil"]
    ],
    ["PLDDT", "RSA", "Helix", "Sheet", "Coil"],
    title="Figure 4(a)"
)

TP_mask = (full_probabilities_dict["y"] == 1) & (full_probabilities_dict["original_probs"] >= 0.26)
FN_mask = (full_probabilities_dict["y"] == 1) & (full_probabilities_dict["original_probs"] < 0.26)
FP_mask = (full_probabilities_dict["y"] == 0) & (full_probabilities_dict["original_probs"] >= 0.26)
TN_mask = (full_probabilities_dict["y"] == 0) & (full_probabilities_dict["original_probs"] < 0.26)

make_violin_plot(
    [
        full_delta_probabilities_dict["plddt"][TP_mask],
        full_delta_probabilities_dict["rsa"][TP_mask],
        full_delta_probabilities_dict["helix"][TP_mask],
        full_delta_probabilities_dict["sheet"][TP_mask],
        full_delta_probabilities_dict["coil"][TP_mask]
    ],
    ["PLDDT", "RSA", "Helix", "Sheet", "Coil"],
    title="Figure 4(b) - True Positive"
)

make_violin_plot(
    [
        full_delta_probabilities_dict["plddt"][FN_mask],
        full_delta_probabilities_dict["rsa"][FN_mask],
        full_delta_probabilities_dict["helix"][FN_mask],
        full_delta_probabilities_dict["sheet"][FN_mask],
        full_delta_probabilities_dict["coil"][FN_mask]
    ],
    ["PLDDT", "RSA", "Helix", "Sheet", "Coil"],
    title="Figure 4(c) - False Negative"
)

make_violin_plot(
    [
        full_delta_probabilities_dict["plddt"][FP_mask],
        full_delta_probabilities_dict["rsa"][FP_mask],
        full_delta_probabilities_dict["helix"][FP_mask],
        full_delta_probabilities_dict["sheet"][FP_mask],
        full_delta_probabilities_dict["coil"][FP_mask]
    ],
    ["PLDDT", "RSA", "Helix", "Sheet", "Coil"],
    title="Figure 4(d) - False Positive"
)

make_violin_plot(
    [
        full_delta_probabilities_dict["plddt"][TN_mask],
        full_delta_probabilities_dict["rsa"][TN_mask],
        full_delta_probabilities_dict["helix"][TN_mask],
        full_delta_probabilities_dict["sheet"][TN_mask],
        full_delta_probabilities_dict["coil"][TN_mask]
    ],
    ["PLDDT", "RSA", "Helix", "Sheet", "Coil"],
    title="Figure 4(e) - True Negative"
)


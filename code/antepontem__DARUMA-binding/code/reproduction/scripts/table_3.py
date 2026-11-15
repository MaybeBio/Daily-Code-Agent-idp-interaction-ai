"""
This script calculates the results for table 3 of DARUMA-binding.
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score


def read_caid_result_csv_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

        result_dict = {
            "seq": [],
            "y": [],
            "af_plddt": [],
            "af_rsa": [],
            "af_ss8": [],
            "db_preds": [],
            "db_plddt": [],
            "db_rsa": [],
            "db_ss3": [],
            "db_ss8": [],
        }

        for line in lines[1:]:
            line = line.strip().split(",")
            result_dict["seq"].append(line[0])
            result_dict["y"].append(int(line[1]))
            result_dict["af_plddt"].append(float(line[2]) if line[2] != "" else None)
            result_dict["af_rsa"].append(float(line[3]) if line[3] != "" else None)
            result_dict["af_ss8"].append(line[4] if line[4] != "" else None)
            result_dict["db_preds"].append(int(line[5]))
            result_dict["db_plddt"].append(float(line[6]))
            result_dict["db_rsa"].append(float(line[7]))
            result_dict["db_ss3"].append(line[8])
            result_dict["db_ss8"].append(line[9])

    return result_dict


def calc_TPR(y, x):
    y = np.asarray(y)
    x = np.asarray(x)
    return f"{np.sum((y == 1) & (x == 1)) / np.sum(y == 1):.3f}"

def calc_FPR(y, x):
    y = np.asarray(y)
    x = np.asarray(x)
    return f"{np.sum((y == 0) & (x == 1)) / np.sum(y == 0):.3f}"

def calc_precision(y, x):
    y = np.asarray(y)
    x = np.asarray(x)
    return f"{np.sum((y == 1) & (x == 1)) / np.sum(x == 1):.3f}"


from pathlib import Path
files = Path("../results/caid3/csv/").glob("*.csv")

full_result_dict = {
    "seq": [],
    "y": [],
    "af_plddt": [],
    "af_rsa": [],
    "af_ss8": [],
    "db_preds": [],
    "db_plddt": [],
    "db_rsa": [],
    "db_ss3": [],
    "db_ss8": [],
}

for file in files:
    result_dict = read_caid_result_csv_file(file)

    for key in full_result_dict.keys():
        full_result_dict[key].extend(result_dict[key])


mask = np.array(full_result_dict["af_plddt"]) != None
for key in full_result_dict.keys():
    full_result_dict[key] = np.array(full_result_dict[key])[mask]


pLDDT_threshold = 70
RSA_threshold = 0.25

HH_mask = (full_result_dict["af_plddt"] > pLDDT_threshold) & (full_result_dict["af_rsa"] > RSA_threshold)
HL_mask = (full_result_dict["af_plddt"] > pLDDT_threshold) & (full_result_dict["af_rsa"] <= RSA_threshold)
LH_mask = (full_result_dict["af_plddt"] <= pLDDT_threshold) & (full_result_dict["af_rsa"] > RSA_threshold)
LL_mask = (full_result_dict["af_plddt"] <= pLDDT_threshold) & (full_result_dict["af_rsa"] <= RSA_threshold)

for mask, name in zip([HH_mask, HL_mask, LH_mask, LL_mask], ["HH", "HL", "LH", "LL"]):
    n = np.sum(mask)
    n_binding = np.sum(full_result_dict["y"][mask] == 1)
    tpr = calc_TPR(full_result_dict["y"][mask], full_result_dict["db_preds"][mask])
    fpr = calc_FPR(full_result_dict["y"][mask], full_result_dict["db_preds"][mask])
    precision = calc_precision(full_result_dict["y"][mask], full_result_dict["db_preds"][mask])
    print(f"{name}: n={n}, n_binding={n_binding}, TPR={tpr}, FPR={fpr}, Precision={precision}")


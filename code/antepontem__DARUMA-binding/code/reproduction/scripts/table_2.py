"""
This script calculates the results for table 2.
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


def calc_PCC(y, x):
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    return f"{np.corrcoef(x, y)[0, 1]:.3f}"

def calc_RMSE(y, x):
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    return f"{np.sqrt(np.mean((y - x) ** 2)):.3f}"

def calc_accuracy(y, x):
    y = np.asarray(y)
    x = np.asarray(x)
    return f"{np.mean(y == x):.3f}"

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

q8_to_q3 = {
    "H": "H",
    "G": "H",
    "I": "H",
    "E": "S",
    "B": "S",
    "T": "C",
    "S": "C",
    "C": "C",
}

# Table 2(a)
print("Table 2(a)")
plddt_pcc = calc_PCC(full_result_dict["af_plddt"], full_result_dict["db_plddt"])
plddt_rmse = calc_RMSE(full_result_dict["af_plddt"]/100, full_result_dict["db_plddt"]/100)
print(f"pLDDT: PCC={plddt_pcc}, RMSE={plddt_rmse}") 

rsa_pcc = calc_PCC(full_result_dict["af_rsa"], full_result_dict["db_rsa"])
rsa_rmse = calc_RMSE(full_result_dict["af_rsa"], full_result_dict["db_rsa"])
print(f"RSA: PCC={rsa_pcc}, RMSE={rsa_rmse}")

q8_acc = calc_accuracy(full_result_dict["af_ss8"], full_result_dict["db_ss8"])
q3_acc = calc_accuracy([q8_to_q3[x] for x in full_result_dict["af_ss8"]], full_result_dict["db_ss3"])
print(f"SS: 8-class accuracy={q8_acc}, 3-class accuracy={q3_acc}")
print()

# Table 2(b)
TP_mask = (np.array(full_result_dict["y"]) == 1) & (np.array(full_result_dict["db_preds"]) == 1)
FN_mask = (np.array(full_result_dict["y"]) == 1) & (np.array(full_result_dict["db_preds"]) == 0)
FP_mask = (np.array(full_result_dict["y"]) == 0) & (np.array(full_result_dict["db_preds"]) == 1)
TN_mask = (np.array(full_result_dict["y"]) == 0) & (np.array(full_result_dict["db_preds"]) == 0)

print("Table 2(b)")
for mask, name in zip([TP_mask, FN_mask, FP_mask, TN_mask], ["TP", "FN", "FP", "TN"]):
    print(name, np.sum(mask))
    plddt_pcc = calc_PCC(full_result_dict["af_plddt"][mask], full_result_dict["db_plddt"][mask])
    plddt_rmse = calc_RMSE(full_result_dict["af_plddt"][mask]/100, full_result_dict["db_plddt"][mask]/100)
    print(f"  pLDDT: PCC={plddt_pcc}, RMSE={plddt_rmse}")
    rsa_pcc = calc_PCC(full_result_dict["af_rsa"][mask], full_result_dict["db_rsa"][mask])
    rsa_rmse = calc_RMSE(full_result_dict["af_rsa"][mask], full_result_dict["db_rsa"][mask])
    print(f"  RSA: PCC={rsa_pcc}, RMSE={rsa_rmse}")
    q8_acc = calc_accuracy(full_result_dict["af_ss8"][mask], full_result_dict["db_ss8"][mask])
    q3_acc = calc_accuracy([q8_to_q3[x] for x in full_result_dict["af_ss8"][mask]], full_result_dict["db_ss3"][mask])
    print(f"  SS: 8-class accuracy={q8_acc}, 3-class accuracy={q3_acc}")


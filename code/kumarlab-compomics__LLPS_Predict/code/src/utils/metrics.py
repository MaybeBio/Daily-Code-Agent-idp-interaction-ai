# evaluation metrics
import numpy as np
from sklearn.metrics import (
    classification_report, 
    confusion_matrix, 
    roc_auc_score,
    average_precision_score, 
    matthews_corrcoef,
    f1_score
)


def classification_metrics(target, pred, threshold: float = 0.5):
    # AUROC
    auroc = roc_auc_score(target, pred)
    print(f"AUROC: {auroc:.2f}")
    # PRAUC
    prauc = average_precision_score(target, pred)
    print(f"PRAUC: {prauc:.2f}")
    # convert predicted probabilities to class labels
    pred_cls = np.where(pred < threshold, 0, 1)
    print('Optimal threshold:', threshold)
    # MMC
    mcc = matthews_corrcoef(target, pred_cls)
    print(f"MCC: {mcc:.2f}")
    # F1 score
    f1 = f1_score(target, pred_cls)
    print(f"F1 score: {f1:.2f}")
    # classification report
    report = classification_report(target, pred_cls)
    print('-'*10)
    print('Classification report:')
    print(report)
    # confusion matrix
    confusion = confusion_matrix(target, pred_cls)
    print('Confusion matrix:')
    print(confusion)
    return auroc, prauc, mcc, f1, report, confusion

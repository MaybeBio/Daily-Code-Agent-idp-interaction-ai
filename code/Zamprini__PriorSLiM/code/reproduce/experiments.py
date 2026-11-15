"""Recompute manuscript statistics from frozen sample-level measurements."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def mutation(input_dir, output_dir):
    d = pd.read_csv(input_dir / 'mutation.tsv', sep='\t')
    rows = []
    for branch in ('final_logit','intrinsic_logit','explicit_logit','motifae_logit','accessibility_logit_change'):
        delta = d['MUT_' + branch] - d['WT_' + branch]
        np.testing.assert_allclose(delta,d['delta_' + branch],atol=1e-6)
        d['correct_' + branch] = delta * d.expected_sign > 0
    for paper in ('A','B','A+B'):
        for group in d.group.unique():
            g = d[d.group.eq(group) & (d.paper.eq(paper) if paper != 'A+B' else True)]
            rows.append(dict(paper=paper,group=group,n=len(g),proteins=g.accession.nunique(),correct=int(g.correct_final_logit.sum()),agreement=g.correct_final_logit.mean()))
    pd.DataFrame(rows).to_csv(output_dir / 'mutation_summary.tsv', sep='\t',index=False)


def cv10(input_dir, output_dir):
    from sklearn.metrics import average_precision_score, roc_auc_score, brier_score_loss, precision_score, recall_score, f1_score
    d = pd.read_csv(input_dir / 'cv10.tsv.gz', sep='\t')
    rows=[]
    for (fold, operating_point, k), g in d.groupby(['fold','operating_point','K']):
        y,p=g.label.to_numpy(),g.probability.to_numpy()
        threshold=g.threshold.iloc[0]
        pred=p >= threshold
        rows.append(dict(fold=fold,operating_point=operating_point,K=k,threshold=threshold,test_ap=average_precision_score(y,p),test_auc=roc_auc_score(y,p),test_brier=brier_score_loss(y,p),test_precision=precision_score(y,pred),test_recall=recall_score(y,pred),test_f1=f1_score(y,pred)))
    pd.DataFrame(rows).to_csv(output_dir / 'cv10_test_per_fold.tsv',sep='\t',index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('experiment', choices=['cv10','functional','structural','mutation','localization','localization-test'])
    parser.add_argument('--input-dir',type=Path,default=Path(__file__).parent/'data/experiments')
    parser.add_argument('--output-dir',type=Path,default=Path('experiment_results'))
    a = parser.parse_args()
    if a.experiment == 'localization-test':
        from localization import self_test
        print(f'{self_test()} localization checks passed')
        return
    a.output_dir.mkdir(parents=True,exist_ok=False)
    if a.experiment == 'mutation':
        mutation(a.input_dir,a.output_dir)
    elif a.experiment == 'cv10':
        cv10(a.input_dir,a.output_dir)
    else:
        if a.experiment == 'structural':
            import structural as module
        elif a.experiment == 'localization':
            import localization as module
        else:
            import functional as module
        module.run(a.input_dir/a.experiment if a.experiment == "functional" else a.input_dir,a.output_dir)


if __name__ == '__main__':
    main()

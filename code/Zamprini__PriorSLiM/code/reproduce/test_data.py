"""Check the result tables used by the numerical figure panels."""
import hashlib,json,unittest
from pathlib import Path
import numpy as np
import pandas as pd
DATA=Path(__file__).resolve().parent/'data'
def table(figure,name):
    return pd.read_csv(DATA/figure/(name+'.tsv'),sep='\t')
def bh(values):
    p=np.asarray(values,float)
    order=np.argsort(p)
    q=np.empty(len(p))
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q
class FigureDataTests(unittest.TestCase):
    def test_hashes(self):
        for name,entry in json.loads((DATA/'manifest.json').read_text()).items():
            self.assertEqual(hashlib.sha256((DATA.parent/name).read_bytes()).hexdigest(),entry['sha256'],name)
    def test_multiple_testing_families(self):
        a=table('fig06','main_results');b=table('fig06','primary360_tests')
        self.assertEqual(len(a),43);self.assertEqual(len(b),360)
        self.assertEqual(int(a.p_value.isna().sum()),6)
        valid=a.p_value.notna()
        np.testing.assert_allclose(bh(a.p_value.fillna(1))[valid],a.loc[valid,'q_joint'],rtol=1e-10,atol=1e-14)
        np.testing.assert_allclose(bh(b.p_value),b.q_bh_primary360,rtol=1e-10,atol=1e-14)
        for contrast,median,count in [('file1_union_vs_file1_outside',1.2349752868,10),('file2_multi_vs_one',1.5884809433,7)]:
            d=b[(b.stratum=='Domain')&(b.min_run==3)&(b.contrast==contrast)]
            self.assertEqual(sorted(d.fold),list(range(1,11)))
            self.assertAlmostEqual(d.effect_ratio.median(),median,places=9)
            self.assertEqual(int((d.q_bh_primary360<.05).sum()),count)
    def test_mutation_cases(self):
        d=table('experiments','mutation')
        self.assertEqual(len(d),31)
        self.assertEqual(sorted(d.groupby('group').size()),[14,17])
    def test_validation_seeds(self):
        for name in ['joint_deletion_endpoint','joint_deletion_curves','retention_near600_display','d_original_absolute_values','e_training_strategy_seed_values','f_budget_envelope']:
            d=table('fig03',name)
            self.assertEqual(sorted(d.seed.unique()),list(range(2026,2031)))
    def test_external_topk(self):
        a=table('fig04','pooled_summary');b=table('fig04','pooled_per_fold')
        for method in ['PriorSLiM','Explicit_only','MotifAE_only']:
            rows=a[(a.method==method)&(a['mode']=='single_window')].sort_values('k')
            self.assertEqual(rows.k.tolist(),list(range(1,11)))
            for r in rows.itertuples():
                d=b[(b.method==method)&(b['mode']=='single_window')&(b.k==r.k)]
                self.assertEqual(len(d),10);self.assertEqual(r.loci,4286)
                self.assertAlmostEqual(d.recall.mean(),r.recall_mean,places=12)
                self.assertAlmostEqual(d.recall.std(ddof=1),r.recall_sd,places=12)
    def test_benchmark_counts(self):
        d=table('fig02','aligned_test_predictions')
        self.assertEqual(len(d),9765);self.assertEqual(int(d.label.sum()),1395)
        self.assertEqual(len(table('fig02','g_all_candidate_dependence')),5880)
        self.assertEqual(len(table('fig02','g_broad_candidate_evidence')),1459)
if __name__=='__main__':unittest.main()

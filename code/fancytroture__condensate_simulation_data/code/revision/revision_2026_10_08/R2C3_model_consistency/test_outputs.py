#!/usr/bin/env python3
"""Independent checks on the response tables, rather than their plot code."""
import csv
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
C2 = OUT.parent / 'R2C2_state_remodeling/response_package'
KEYS = [(s,f,r) for s in ('FUS','A1','TDP43') for f in ('free','condensate') for r in ('Tl','Th')]


def read(name, directory=OUT):
    with (directory/name).open() as handle:
        return list(csv.DictReader(handle))


def select(rows,key):
    return [r for r in rows if (r['system'],r['family'],r['role'])==key]


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for data in iter(lambda:handle.read(8*1024*1024),b''):
            digest.update(data)
    return digest.hexdigest()


class OutputTests(unittest.TestCase):
    def test_complete_matrices_and_populations(self):
        matrices = read('Table_C3_3_full_transition_matrices.csv')
        populations = read('Table_C3_1_all_populations.csv')
        self.assertEqual(len(matrices),300)
        self.assertEqual(len(populations),60)
        reference = read('full_transition_matrices.csv',C2)
        for key in KEYS:
            rows = select(matrices,key)
            t = np.array([float(r['transition_probability']) for r in rows]).reshape(5,5)
            mu = np.array([float(r['mean_population']) for r in select(populations,key)])
            self.assertTrue(np.isfinite(t).all() and (t>=0).all())
            np.testing.assert_allclose(t.sum(1),1,atol=1e-12)
            np.testing.assert_allclose(mu.sum(),1,atol=1e-12)
            expected = np.array([float(r['T']) for r in select(reference,key)]).reshape(5,5)
            np.testing.assert_allclose(t,expected,atol=1e-14,rtol=0)

    def test_distances_from_exported_full_matrices(self):
        matrices = read('Table_C3_3_full_transition_matrices.csv')
        populations = read('Table_C3_1_all_populations.csv')
        assignments = read('common_state_correspondence.csv')
        distances = read('Table_C3_4_six_comparisons.csv')
        self.assertEqual(len(distances),6)
        for row in distances:
            system,family = row['system'],row['family']
            pairs = [r for r in assignments if r['system']==system and r['family']==family]
            local = [[int(r[f'local_{role}'])-1 for r in pairs] for role in ('TL','TH')]
            arrays,pops = [],[]
            for role,indices in zip(('Tl','Th'),local):
                t = np.array([float(r['transition_probability']) for r in select(matrices,(system,family,role))]).reshape(5,5)
                mu = np.array([float(r['mean_population']) for r in select(populations,(system,family,role))])
                arrays.append(t[np.ix_(indices,indices)])
                mu = mu[indices]
                pops.append(mu/mu.sum())
            # Explicit double loop, independent of the production metric implementation.
            squared = sum((pops[0][i]+pops[1][i])/2*(arrays[1][i,j]-arrays[0][i,j])**2
                          for i in range(len(pairs)) for j in range(len(pairs)))
            self.assertAlmostEqual(squared**.5,float(row['full_cache_estimate']),places=13)

    def test_complete_model_pair_diagnostics(self):
        pairs = read('all_2280_model_pairs.csv')
        agreement = read('Table_C3_2_state_agreement.csv')
        self.assertEqual(len(pairs),2280)
        self.assertEqual(len(agreement),60)
        for key in KEYS:
            rows = select(pairs,key)
            self.assertEqual(len(rows),190)
            self.assertEqual(len({(r['attempt_a'],r['attempt_b']) for r in rows}),190)
            for r in rows:
                self.assertGreaterEqual(float(r['optimized_hard_agreement'])+1e-12,float(r['hard_agreement']))
            for r in select(agreement,key):
                self.assertTrue(0<=float(r['soft_dice_mean'])<=1+1e-6)
                self.assertTrue(0<=int(r['hard_dice_valid_pairs'])<=190)

    def test_reviewers_specific_numbers(self):
        summaries = read('condition_summary.csv')
        hi = select(summaries,('TDP43','condensate','Th'))[0]
        self.assertAlmostEqual(float(hi['max_mean_population']),.3699592593,places=9)
        self.assertAlmostEqual(float(hi['mean_model_max_population']),.5059238391,places=9)
        self.assertEqual(json.loads(hi['dominant_model_counts']),[11,8,1,0,0])
        free = select(summaries,('TDP43','free','Th'))[0]
        self.assertAlmostEqual(float(free['omitted_population']),.0774448652,places=9)
        self.assertEqual(sum(int(r['n_mean_states_gt_1pct'])==3 for r in summaries),8)
        self.assertEqual(sum(int(r['n_mean_states_gt_1pct'])==4 for r in summaries),4)
        values = {r['definition']:float(r['distance']) for r in read('Table_C3_6_number_reconciliation.csv')
                  if r['system']=='A1' and r['family']=='free'}
        self.assertAlmostEqual(values['submitted_edge_mask'],.128849907657,places=11)
        self.assertAlmostEqual(values['unfiltered_display_union'],.053804006268,places=11)

    def test_current_C2_intervals_not_old_versions(self):
        rows = read('Table_C3_4_temperature_distances.csv')
        ref = read('Table_R2C2_2_metrics.csv',C2)
        for r in rows:
            if r['uncertainty'] not in ('sampling','model_seeds'):
                continue
            e = next(x for x in ref if (x['system'],x['family'],x['uncertainty'],x['metric'])==
                     (r['system'],r['family'],r['uncertainty'],'D_common_pooled_weights'))
            for k in ('estimate','ci_low','ci_high'):
                self.assertAlmostEqual(float(r[k]),float(e[k]),places=13)

    def test_original_files_and_compact_inputs_unchanged(self):
        before = json.loads((OUT/'protected_original_files.json').read_text())
        for name,digest in before.items():
            self.assertEqual(sha(ROOT/name),digest,name)
        for name,digest in json.loads((OUT/'input_sha256.json').read_text()).items():
            self.assertEqual(sha(Path(name)),digest,name)

    def test_block_sensitivity_and_subset_disclosure(self):
        contrast = read('Table_C3_5_condensate_minus_free.csv')
        for system in ('FUS','A1','TDP43'):
            for r in contrast:
                if r['system']==system and r['uncertainty'].startswith('sampling'):
                    if system=='FUS':
                        self.assertLess(float(r['ci_low']),0)
                        self.assertGreater(float(r['ci_high']),0)
                    else:
                        self.assertGreater(float(r['ci_low']),0)
        for r in read('Table_C3_4_six_comparisons.csv'):
            delta = abs(float(r['sampling_eligible_estimate'])-float(r['full_cache_estimate']))
            if r['system']=='A1' and r['family']=='condensate':
                self.assertAlmostEqual(delta,.00026991793692242,places=13)
            else:
                self.assertLess(delta,1e-13)


if __name__=='__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(OutputTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    (OUT/'test_results.json').write_text(json.dumps(dict(tests=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),passed=result.wasSuccessful()),indent=2)+'\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)

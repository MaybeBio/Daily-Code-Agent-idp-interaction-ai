#!/usr/bin/env python3
"""Small independent checks of the estimators and exported results."""
import unittest
import numpy as np
import analyze_followup as a


class FollowupTests(unittest.TestCase):
    def test_transition_estimator_and_missing_row(self):
        c = np.array([[2., 4.], [0., 6.]])
        np.testing.assert_allclose(a.transition(c), [[.5,.5],[.25,.75]])
        with self.assertRaises(AssertionError):
            a.transition(np.zeros((2,2)))

    def test_population_change_is_not_transition_change(self):
        t = np.array([[[.8,.2],[.3,.7]]])
        arrays = [(t, np.array([[.8,.2]])), (t, np.array([[.2,.8]]))]
        metrics, *_ = a.metrics(arrays, [{0:0,1:1}]*2)
        self.assertEqual(float(metrics['D_common_pooled_weights'][0]), 0.)
        self.assertGreater(float(metrics['population_JS_distance_base2_conditional_C'][0]), 0.)

    def test_label_permutation_and_substochastic_rows(self):
        t = np.array([[[.7,.2,.1],[.1,.7,.2],[.3,.3,.4]]])
        p = np.array([[.4,.4,.2]])
        perm = [1,2,0]
        arrays = [(t,p),(t[:,perm][:,:,perm],p[:,perm])]
        metrics, *_ = a.metrics(arrays,[{0:0,1:1},{0:2,1:0}])
        self.assertEqual(float(metrics['D_common_pooled_weights'][0]), 0.)
        np.testing.assert_allclose(t[0,:2,:2].sum(1), [.9,.8])

    def test_contact_index_roundtrip(self):
        x,y=np.triu_indices(146,1)
        old=np.zeros((5,146,146))
        old[:,x,y]=np.arange(len(x))
        new=a.indexing_diagnostic(old)
        np.testing.assert_array_equal(new[:,x,y+2],old[:,x,y])
        self.assertEqual(new[0,0,3],old[0,0,1])

    def test_exported_entries_against_full_matrices(self):
        for row in a.old.read(a.OUT/'common_transition_entry_changes.csv'):
            if row['uncertainty']!='model_seeds':
                continue
            pair=[a.old.load_full(a.FULL,row['system'],row['family'],role) for role in a.s.ROLES]
            i,j=int(row['state'])-1,int(row['target_state'])-1
            values=[d['t'][d['mapping'][i],d['mapping'][j]] for d in pair]
            self.assertAlmostEqual(float(row['estimate']),values[1]-values[0],places=12)

    def test_no_unmapped_states_exported_as_zero(self):
        for row in a.old.read(a.OUT/'population_changes.csv'):
            common=set.intersection(*[set(m) for m in a.maps(row['system'],row['family'])])
            self.assertIn(int(row['state'])-1,common)
            self.assertGreater(float(row['population_TL']),0.)
            self.assertGreater(float(row['population_TH']),0.)


if __name__=='__main__':
    unittest.main()

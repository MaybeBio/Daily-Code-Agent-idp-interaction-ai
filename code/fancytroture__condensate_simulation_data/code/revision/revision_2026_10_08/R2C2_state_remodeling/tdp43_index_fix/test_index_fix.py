#!/usr/bin/env python3
import json
import unittest
import numpy as np
from scipy.sparse import csr_matrix
import rebuild_contacts as r
import compare_fixed as c


class IndexFixTests(unittest.TestCase):
    def test_pair_reconstruction(self):
        v=np.arange(10585,dtype=np.float32)[None,None,:]
        old=r.make_maps(v,146)
        new=r.make_maps(v,148)
        i,j=np.triu_indices(146,1)
        np.testing.assert_array_equal(old[...,i,j],new[...,i,j+2])
        self.assertTrue(np.all(new[...,np.arange(147),np.arange(1,148)]==0))
        self.assertEqual(new.shape,(1,1,148,148))

    def test_sparse_weighted_sum(self):
        rng=np.random.default_rng(1)
        contacts=rng.random((35,91))<.04
        q=rng.random((35,10))
        q/=q.sum(1,keepdims=True)
        sparse=(csr_matrix(contacts).T@q).T
        np.testing.assert_allclose(sparse,q.T@contacts,atol=1e-14)

    def test_permutation_of_sufficient_statistics(self):
        rng=np.random.default_rng(2)
        counts=rng.random((2,7,5,5));mass=rng.random((2,7,5))
        p=np.array([[1,0,2,4,3],[4,1,2,0,3]])
        cc,mm=c.permute_statistics(counts,mass,p)
        for model in range(2):
            for block in range(7):
                np.testing.assert_array_equal(cc[model,block],counts[model,block][np.ix_(p[model],p[model])])
                np.testing.assert_array_equal(mm[model,block],mass[model,block][p[model]])

    def test_trajectory_validation(self):
        records=json.loads((r.OUT/'trajectory_pair_validation.json').read_text())
        self.assertEqual(len(records),4)
        for x in records:
            self.assertEqual(x['atoms'],148)
            self.assertLessEqual(x['first_frame_distance_max_error'],1e-5)

    def test_all_legacy_alignments_reproduced(self):
        records=json.loads((r.OUT/'contact_rebuild_validation.json').read_text())
        self.assertEqual(len(records),4)
        for x in records:
            self.assertEqual(x['legacy']['changed_models'],0)
            if x['rebuild_performed']:
                self.assertLess(x['old_average_error'][0],2e-4)
            else:
                self.assertEqual(x['saved_shape'],[5,148,148])

    def test_probabilities_normalized_after_fix(self):
        records=c.previous.old.read(r.OUT/'fixed_full_local_matrices.csv')
        for family in c.s.FAMILIES:
            for role in c.s.ROLES:
                rows=[x for x in records if x['family']==family and x['role']==role]
                self.assertEqual(len(rows),25)
                for i in range(5):
                    row=[x for x in rows if int(x['local_source'])==i]
                    self.assertAlmostEqual(sum(float(x['transition_probability']) for x in row),1.)


if __name__=='__main__':unittest.main()

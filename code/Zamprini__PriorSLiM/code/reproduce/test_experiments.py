"""Numerical boundary checks for the experiment statistics."""
import unittest
import numpy as np
import functional
import localization
import robust_glm


class ExperimentTests(unittest.TestCase):
    def test_localization(self):
        self.assertEqual(localization.self_test(), 205)

    def test_bh(self):
        np.testing.assert_allclose(functional.bh([0.01, 0.04, 0.03], 4),[0.04, 0.05333333333333334, 0.05333333333333334])

    def test_dyadic_covariance(self):
        rng=np.random.default_rng(2027)
        score=rng.normal(size=(25,4))
        a=rng.integers(0,8,25).astype(str)
        b=rng.integers(0,8,25).astype(str)
        naive=np.zeros((4,4))
        for i in range(25):
            for j in range(25):
                if {a[i],b[i]} & {a[j],b[j]}:
                    naive += np.outer(score[i],score[j])
        np.testing.assert_allclose(robust_glm.dyadic_cov(score,np.eye(4),a,b),naive*25/21,atol=1e-12)


if __name__ == '__main__':
    unittest.main()

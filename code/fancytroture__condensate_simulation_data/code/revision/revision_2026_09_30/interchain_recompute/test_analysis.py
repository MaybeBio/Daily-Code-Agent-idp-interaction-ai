#!/usr/bin/env python3
import unittest
import numpy as np
from analyze import runs, correlation_values, partner_values, crossing, interval, ratio_bounds
from compute import pilot


class AnalysisTest(unittest.TestCase):
    def test_missing_frame_splits(self):
        rr=runs(np.array([0,50,150,200]),np.array([0.,.5,1.5,2.]))
        self.assertEqual([(s.start,s.stop) for s in rr],[(0,2),(2,4)])

    def test_header_jump_splits(self):
        rr=runs(np.array([0,50,100,150]),np.array([0.,.5,31.,31.5]))
        self.assertEqual([(s.start,s.stop) for s in rr],[(0,2),(2,4)])

    def test_covariance_fft_reference(self):
        legacy,_=pilot.load_legacy()
        x=np.random.default_rng(22).normal(size=300)
        direct=np.array([correlation_values(x,k).mean() for k in range(80)])
        np.testing.assert_allclose(direct,legacy.autocorr_fft(x)[:80],atol=1e-14)
        self.assertAlmostEqual(direct[0],1)

    def test_constant_no_fake_acf(self):
        self.assertEqual(len(correlation_values(np.ones(20),0)),0)

    def test_partner_return_not_continuous_survival(self):
        c=np.array([[1,0],[0,1],[1,0]])
        np.testing.assert_array_equal(partner_values(c,1),[0,0])
        np.testing.assert_array_equal(partner_values(c,2),[1])

    def test_empty_origin_is_missing(self):
        c=np.array([[0,0],[1,0]])
        self.assertTrue(np.isnan(partner_values(c,1)[0]))

    def test_crossing_interpolation(self):
        self.assertAlmostEqual(crossing(np.array([0.,.5,1.]),np.array([1.,.4,0.]),.2),.75)

    def test_censored_not_dropped(self):
        self.assertTrue(np.isinf(crossing(np.array([0.,1.]),np.array([1.,.3]),.2)))
        lo,hi=interval(np.r_[np.ones(90),np.full(10,np.inf)])
        self.assertEqual(lo,1)
        self.assertTrue(np.isinf(hi))

    def test_censored_ratio_bounds(self):
        lo,hi=ratio_bounds(np.array([np.inf,10.,np.inf]),np.array([10.,np.inf,np.inf]),30.)
        np.testing.assert_array_equal(lo,[3.,0.,0.])
        np.testing.assert_array_equal(hi,[np.inf,1/3,np.inf])


if __name__=='__main__':
    unittest.main()

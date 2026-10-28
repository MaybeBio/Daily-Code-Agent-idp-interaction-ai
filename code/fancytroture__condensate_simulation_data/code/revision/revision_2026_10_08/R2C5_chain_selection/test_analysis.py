import unittest
import numpy as np
from extract import contacts,check_counts
from analyze import eligible,unwrap,curves,all_inside,crossing


class Checks(unittest.TestCase):
    def test_contacts_periodic_and_reference(self):
        rng=np.random.default_rng(17)
        xyz=rng.uniform(0,4,(600,3));box=np.full(3,4.)
        xyz[:2]=[[.01,.2,.2],[3.99,.2,.2]]
        c=contacts(xyz,box,3)
        check_counts(xyz,box,3,c)
        np.testing.assert_array_equal(c,c.T)
        self.assertFalse(c.diagonal().any())

    def test_residence_nested(self):
        flags=np.ones((2500,200),bool);flags[700,0]=False
        lo,sh=eligible(flags)
        self.assertFalse(lo[0]);self.assertFalse(sh[1,0]);self.assertTrue(sh[0,0])
        self.assertEqual(sh[:,0].sum(),4)

    def test_unwrap_drift(self):
        box=np.full(3,10.)
        t=np.arange(100)*.1
        coords=np.zeros((100,2,3));coords[:,0,0]=9.+t;coords[:,1,0]=2.-t
        out=unwrap(coords%box,box)
        np.testing.assert_allclose(out[:,0,0],t,atol=1e-12)
        np.testing.assert_allclose(out[:,1,0],-t,atol=1e-12)

    def test_inside_endpoints(self):
        flags=np.ones((8,1),bool);flags[4]=False
        np.testing.assert_array_equal(all_inside(flags,3)[:,0],[1,1,0,0,0,1])

    def test_acf_partner_bruteforce(self):
        rng=np.random.default_rng(6)
        c=rng.integers(0,3,(200,3,200))
        ac,ps=curves(c)
        np.testing.assert_allclose(ac[:,0],1.)
        np.testing.assert_allclose(ps[:,0],1.)
        x=c.sum(axis=2).astype(float);x-=x.mean(axis=0)
        for k in [1,10,60]:
            ref=(x[:-k]*x[k:]).mean(axis=0)/(x*x).mean(axis=0)
            np.testing.assert_allclose(ac[:,k],ref)
            p=c>0
            ref=((p[:-k]&p[k:]).sum(axis=2)/p[:-k].sum(axis=2)).mean(axis=0)
            np.testing.assert_allclose(ps[:,k],ref)

    def test_crossing_censoring(self):
        self.assertTrue(np.isinf(crossing(np.ones(61),.2)))
        self.assertAlmostEqual(crossing(1-np.arange(61)/60,.2),24.)


if __name__=='__main__':unittest.main()

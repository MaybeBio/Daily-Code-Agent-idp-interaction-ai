#!/usr/bin/env python

import unittest
import os
import shutil
import sys
import subprocess
import ihm.reader
import pickle
import IMP


TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))


class Tests(unittest.TestCase):
    def test_simple(self):
        """Test model building"""
        os.chdir(os.path.join(TOPDIR, 'production_scripts'))
        p = subprocess.check_call(["python", "sample.py", "--test"])
        # todo: assert outputs, run analysis

    def test_pickle(self):
        """Test that pickled ReplicaExchange object works"""
        # Set up modeling but don't run sampling
        os.chdir(os.path.join(TOPDIR, 'production_scripts'))
        with open('sample.py') as fh:
            contents = fh.read().replace('rex.execute_macro()', '')
        sys.argv = [sys.argv[0], '--test']
        g = {}
        exec(contents, g)
        rex = g['rex']
        del g
        rex.vars['number_of_frames'] = 2

        dump = pickle.dumps((rex.model, rex))

        # Run the original ReplicaExchange and get the final score
        IMP.random_number_generator.seed(99)
        rex.execute_macro()
        rs = IMP.pmi.tools.get_restraint_set(rex.model)
        original_score = rs.evaluate(False)
        del rex, rs

        # With the same random seed, we should get the exact same trajectory
        # with the pickled object
        newm, newrex = pickle.loads(dump)
        IMP.random_number_generator.seed(99)
        newrex.execute_macro()
        rs = IMP.pmi.tools.get_restraint_set(newrex.model)
        new_score = rs.evaluate(False)
        self.assertAlmostEqual(original_score, new_score, delta=1e-4)

    def test_mmcif(self):
        """Test generation of mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'production_scripts'))
        if os.path.exists("pol_ii_g.cif"):
            os.unlink("pol_ii_g.cif")
        # Potentially override methods that need network access
        env = os.environ.copy()
        env['PYTHONPATH'] = os.path.join(TOPDIR, 'test', 'mock') \
                            + ':' + env.get('PYTHONPATH', '')
        p = subprocess.check_call(
                ["python", "sample.py", "--mmcif", "--dry-run"], env=env)
        # Check output file
        self._check_mmcif_file('pol_ii_g.cif')

    def _check_mmcif_file(self, fname):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 3)
        self.assertEqual(s.citations[0].doi, '10.1038/s41594-018-0118-5')
        self.assertEqual(len(s.software), 2)
        self.assertEqual(len(s.orphan_starting_models), 12)
        # Should be 1 state
        self.assertEqual(len(s.state_groups), 1)
        state1, = s.state_groups[0]
        # Should be 1 model
        self.assertEqual(sum(len(x) for x in state1), 1)
        # Check # of spheres and atoms in each model
        m = state1[0][0]
        self.assertEqual(len(m._spheres), 3640)
        self.assertEqual(len(m._atoms), 0)
        # Should be 1 ensemble
        self.assertEqual([e.num_models for e in s.ensembles], [1640])
        # Just one restraint - crosslinks
        xl, = s.restraints
        self.assertEqual(len(xl.experimental_cross_links), 40)
        self.assertEqual(len(xl.cross_links), 40)


if __name__ == '__main__':
    unittest.main()

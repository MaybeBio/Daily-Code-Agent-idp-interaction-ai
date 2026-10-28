#!/usr/bin/env python3
"""Focused numerical and provenance checks for saved VAMPNet resampling."""
import unittest
import numpy as np
import analyze as a
from build_absolute import accumulate


class NumericalTests(unittest.TestCase):
    def test_reversible_soft_counts(self):
        rng = np.random.default_rng(7)
        q = rng.dirichlet(np.ones(5), 4100)
        c = q[:-200].T @ q[200:]
        t = a.transition(c)
        symmetric = .5 * (c + c.T)
        pi = symmetric.sum(axis=1) / symmetric.sum()
        np.testing.assert_allclose(t.sum(axis=1), 1)
        np.testing.assert_allclose(pi[:, None] * t, pi[None, :] * t.T)

    def test_absolute_blocks_preserve_all_pairs_no_joins(self):
        rng = np.random.default_rng(3)
        q1, q2 = rng.dirichlet(np.ones(5), 4317), rng.dirichlet(np.ones(5), 1833)
        counts, masses = np.zeros((10, 5, 5)), np.zeros((10, 5))
        accumulate(q1, 5317, 5000, counts, masses)
        accumulate(q2, 6379, 5000, counts, masses)
        np.testing.assert_allclose(counts.sum(0), q1[:-200].T @ q1[200:] + q2[:-200].T @ q2[200:])
        np.testing.assert_allclose(masses.sum(0), q1.sum(0) + q2.sum(0))
        np.testing.assert_allclose(counts.sum(), len(q1) + len(q2) - 400)

    def test_identical_units_no_artificial_seed_variance(self):
        rng = np.random.default_rng(3)
        c = rng.uniform(size=(3, 1, 5, 5))
        m = rng.uniform(size=(3, 1, 5))
        t, p = a.draw_models(np.repeat(c, 8, axis=1), np.repeat(m, 8, axis=1), 50, rng)
        np.testing.assert_allclose(t, np.repeat(t[:1], len(t), axis=0))
        np.testing.assert_allclose(p, np.repeat(p[:1], len(p), axis=0))

    def test_common_exit_remains_stochastic(self):
        t = np.array([[.4, .3, .3], [.1, .4, .5], [.2, .3, .5]])
        p = np.array([.2, .5, .3])
        common = np.array([True, True, False])
        result = a.audit.metric('common_plus_exit', t, t, p, p, common, None, None)
        self.assertEqual(float(result), 0.)


class SavedDataTests(unittest.TestCase):
    def test_absolute_and_trajectory_counts_agree(self):
        for system in a.s.SYSTEMS:
            for role in a.s.ROLES:
                tag = f'{system}_condensate_{a.s.TEMPS[system][role]}'
                with np.load(a.OUT / f'{tag}_sufficient_statistics.npz') as z:
                    c, m = z['counts'], z['masses']
                with np.load(a.OUT / f'{tag}_absolute_statistics.npz') as z:
                    ac, am, excluded = z['counts'], z['masses'], z['excluded_trajectories']
                keep = np.ones(c.shape[1], dtype=bool)
                keep[excluded] = False
                np.testing.assert_allclose(ac.sum(1), c[:, keep].sum(axis=(1, 2)), rtol=1e-10)
                np.testing.assert_allclose(am.sum(1), m[:, keep].sum(axis=(1, 2)), rtol=1e-10)

    def test_sampling_point_matches_audit_except_excluded_chain(self):
        for system in a.s.SYSTEMS:
            states = a.audit.original.transition_states(system)
            for family in a.s.FAMILIES:
                if system == 'A1' and family == 'condensate':
                    continue
                data = [a.audit.load_condition(system, family, role, states) for role in a.s.ROLES]
                fitted = []
                for role in a.s.ROLES:
                    mode = 'trajectory' if family == 'free' else 'absolute_slab_time_block'
                    c, m, _ = a.units(system, family, role, mode, 5000)
                    t, p = a.draw_models(c, m, 2, np.random.default_rng(1))
                    fitted.append(a.mapped(t, p, a.s.display_state_mapping(system, family, role), states))
                common = np.array([v in data[0]['retained'] & data[1]['retained'] for v in states])
                for metric in a.audit.VARIANTS:
                    expected = a.audit.metric(metric, data[0]['t'].mean(0), data[1]['t'].mean(0),
                        data[0]['p'].mean(0), data[1]['p'].mean(0), common, data[0]['mask'], data[1]['mask'])
                    actual = a.audit.metric(metric, fitted[0][0][0], fitted[1][0][0], fitted[0][1][0],
                        fitted[1][1][0], common, data[0]['mask'], data[1]['mask'])
                    self.assertAlmostEqual(float(expected), float(actual), places=9)


if __name__ == '__main__':
    unittest.main(verbosity=2)

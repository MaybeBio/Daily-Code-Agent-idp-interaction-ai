#!/usr/bin/env python3
import unittest
import numpy as np
from compute import partner_contact_counts, pilot


class ContactsTest(unittest.TestCase):
    def test_periodic_pair_and_same_chain_exclusion(self):
        xyz = np.arange(200*2*3, dtype=float).reshape(400, 3) + 20
        xyz[0] = [.1, .1, .1]
        xyz[1] = [.15, .1, .1]
        xyz[2] = [9.9, .1, .1]
        # The other beads are isolated in a large box; use a dedicated small system
        # for the periodic check, still retaining the 200-chain output format.
        xyz = np.array([[.1, 0, 0], [.15, 0, 0], [9.9, 0, 0], [5, 5, 5]])
        counts = partner_contact_counts(xyz, np.array([10., 10., 10.]), 2, [0, 1], .3)
        self.assertEqual(counts[0][0], 0)
        self.assertEqual(counts[0][1], 2)
        self.assertEqual(counts[1][0], 2)

    def test_random_bruteforce(self):
        rng = np.random.default_rng(94)
        xyz = rng.uniform(-5, 15, size=(20, 3))
        box = np.array([10., 10., 10.])
        counts = partner_contact_counts(xyz, box, 2, list(range(10)), 2.)
        for chain in range(10):
            ref = np.zeros(200, dtype=int)
            for a in range(chain*2, chain*2+2):
                for b in range(20):
                    d = xyz[a]-xyz[b]
                    d -= box*np.round(d/box)
                    if b//2 != chain and np.dot(d, d) <= 4:
                        ref[b//2] += 1
            np.testing.assert_array_equal(counts[chain], ref)

    def test_legacy_equivalence(self):
        legacy, _ = pilot.load_legacy()
        rng = np.random.default_rng(31)
        xyz = rng.random((200*3, 3))*5
        box = np.full(3, 5.)
        ids = [0, 2, 199]
        result = partner_contact_counts(xyz, box, 3, ids)
        ref = legacy.count_contacts_for_chains(xyz, box, 3, ids, .66, np.repeat(np.arange(200), 3))
        self.assertEqual({c: int(v.sum()) for c,v in result.items()}, ref)

    def test_invalid_box(self):
        with self.assertRaises(ValueError):
            partner_contact_counts(np.zeros((2, 3)), np.zeros(3), 1, [0])


if __name__ == '__main__':
    unittest.main()

"""
test_sojourn_tracking.py

Tests for the A/B sojourn-time tracking added to
viz_util/overdamped_langevin.py (SojournTracker / track_sojourns, and the
OverdampedHarmonicLangevin._solve plumbing that drives it).

Geometry: two regions of a 1-D coordinate separated by a buffer,

    A = {r <= a}        B = {r >= b}        buffer = (a, b),  a <= b

with Schmitt-trigger (hysteresis) assignment: a walker stays "assigned" to the
region it last reached; excursions into the buffer change nothing. A sojourn
that opens at step s and is closed by a crossing detected at step i has duration
``i - s`` (the very first sojourn of each walker is taken to enter at step 0).

The handcrafted trajectories below pin down every one of those rules; the
randomized block additionally proves the vectorised implementation is bit-for-bit
identical to a transcription of the original scalar double-loop.

Run with:
    python -m pytest test/test_sojourn_tracking.py -v
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pytest
from viz_util.overdamped_langevin import (
    SojournTracker,
    track_sojourns,
    OverdampedHarmonicLangevin,
)


# ---------------------------------------------------------------------------
# Scalar reference: a faithful transcription of the ORIGINAL (pre-vectorised)
# inner double-loop, used to cross-check the vectorised version.
# ---------------------------------------------------------------------------

def scalar_reference(positions, a, b):
    positions = np.asarray(positions, dtype=float)
    if positions.ndim == 1:
        positions = positions[:, None]
    n, N = positions.shape
    A_intervals = [[] for _ in range(N)]
    B_intervals = [[] for _ in range(N)]
    in_A = positions[0] <= a
    in_B = positions[0] > a
    A_entry = np.zeros(N)
    B_entry = np.zeros(N)
    for i in range(1, n):
        r = positions[i]
        for k in range(N):
            if r[k] >= b and in_A[k]:
                B_entry[k] = i
                A_intervals[k].append(B_entry[k] - A_entry[k])
                in_A[k] = False
                in_B[k] = True
            if r[k] <= a and in_B[k]:
                A_entry[k] = i
                B_intervals[k].append(A_entry[k] - B_entry[k])
                in_B[k] = False
                in_A[k] = True
    return A_intervals, B_intervals


def as_lists(intervals):
    """Flatten a list of arrays/lists to a list of plain int lists for == comparison."""
    return [[int(v) for v in walker] for walker in intervals]


# ---------------------------------------------------------------------------
# Handcrafted single-walker trajectories
# ---------------------------------------------------------------------------

class TestSingleWalker:
    A, B = -1.0, 1.0

    def test_full_A_B_A_B_sequence(self):
        # entries:  A at step 0, B at 3, A at 5, B at 7
        traj = [-2.0, -0.5, 0.5, 1.5, 0.3, -1.5, 0.5, 1.2]
        A, B = track_sojourns(traj, self.A, self.B)
        assert as_lists(A) == [[3 - 0, 7 - 5]]   # == [[2, 1]]
        assert as_lists(B) == [[5 - 3]]              # == [[1]]

    def test_stays_in_A_no_crossing(self):
        # never reaches b -> assignment never leaves A -> nothing recorded
        traj = [-2.0, -0.9, -0.5, -3.0, -1.5]
        A, B = track_sojourns(traj, self.A, self.B)
        assert as_lists(A) == [[]]
        assert as_lists(B) == [[]]

    def test_buffer_excursion_does_not_reset_entry(self):
        # dips below a (re-entering A physically) but, already assigned to A,
        # records nothing until it finally reaches b -> single A-sojourn from 0.
        traj = [-2.0, 0.0, -1.5, 0.0, 1.5]
        A, B = track_sojourns(traj, self.A, self.B)
        assert as_lists(A) == [[4 - 0]]              # == [[3]]
        assert as_lists(B) == [[]]

    def test_touching_boundary_counts_as_crossing(self):
        # r == b and r == a are inclusive (>= b, <= a)
        traj = [-1.0, 1.0, -1.0]
        A, B = track_sojourns(traj, self.A, self.B)
        assert as_lists(A) == [[1 - 0]]              # == [[0]]
        assert as_lists(B) == [[2 - 1]]              # == [[0]]

    def test_initial_buffer_position_is_assigned_to_B(self):
        # r0 in the buffer (r0 > a) is classified as B, per init convention.
        traj = [0.0, -2.0]                                # B-sojourn closes at step 1
        A, B = track_sojourns(traj, self.A, self.B)
        assert as_lists(A) == [[]]
        assert as_lists(B) == [[1 - 0]]              # == [[0]]


# ---------------------------------------------------------------------------
# Vectorisation across walkers: each column evolves independently
# ---------------------------------------------------------------------------

class TestMultiWalker:
    A, B = -1.0, 1.0

    def test_independent_columns(self):
        # walker 0 oscillates A<->B; walker 1 never reaches b.
        positions = np.array([
            [-2.0, -2.0],
            [ 1.5,  0.0],
            [-1.5,  0.5],
            [ 1.2, -0.5],
        ])
        A, B = track_sojourns(positions, self.A, self.B)
        assert as_lists(A) == [[1, 1], []]
        assert as_lists(B) == [[1], []]

    def test_matches_per_column_single_runs(self):
        rng = np.random.default_rng(0)
        positions = rng.normal(scale=1.5, size=(200, 6))
        A_all, B_all = track_sojourns(positions, self.A, self.B)
        for k in range(positions.shape[1]):
            A_k, B_k = track_sojourns(positions[:, k], self.A, self.B)
            assert as_lists([A_all[k]]) == as_lists(A_k)
            assert as_lists([B_all[k]]) == as_lists(B_k)


# ---------------------------------------------------------------------------
# Output data structure
# ---------------------------------------------------------------------------

class TestDataStructure:
    def test_returns_list_of_int_arrays(self):
        A, B = track_sojourns([-2.0, 1.5, -1.5], -1.0, 1.0)
        assert isinstance(A, list) and isinstance(B, list)
        assert len(A) == 1 and len(B) == 1
        assert isinstance(A[0], np.ndarray)
        assert np.issubdtype(A[0].dtype, np.integer)

    def test_durations_nonnegative_when_a_strictly_less_than_b(self):
        rng = np.random.default_rng(7)
        positions = rng.normal(scale=2.0, size=(500, 4))
        A, B = track_sojourns(positions, -1.0, 1.0)
        for walker in A + B:
            assert np.all(walker >= 0)


# ---------------------------------------------------------------------------
# Randomized equivalence to the original scalar implementation
# ---------------------------------------------------------------------------

class TestEquivalenceToScalarReference:

    @pytest.mark.parametrize("seed", range(8))
    @pytest.mark.parametrize("ab", [(-1.0, 1.0), (-0.2, 0.3), (0.0, 0.0)])
    def test_random_walks(self, seed, ab):
        a, b = ab
        rng = np.random.default_rng(seed)
        # random walks wide enough to cross both boundaries many times
        steps = rng.normal(scale=0.5, size=(300, 5))
        positions = np.cumsum(steps, axis=0)
        A_ref, B_ref = scalar_reference(positions, a, b)
        A_vec, B_vec = track_sojourns(positions, a, b)
        assert as_lists(A_vec) == as_lists(A_ref)
        assert as_lists(B_vec) == as_lists(B_ref)


# ---------------------------------------------------------------------------
# End-to-end through OverdampedHarmonicLangevin._solve
# ---------------------------------------------------------------------------

def _run_with_prescribed_path(path, a, b):
    """
    Drive a real OverdampedHarmonicLangevin so its high-res coordinate follows
    `path` exactly: zero force, and noise increments injected so that
    diff*sqrt(dt) == 1 turns self.noise into the per-step displacement.
    """
    path = np.asarray(path, dtype=float)
    n = len(path)
    incr = np.diff(path)
    # gamma=kB=1, T=0.5, dt=1  ->  diff = sqrt(2*kB*T/gamma) = 1, sqrt(dt) = 1
    assert np.isclose(np.sqrt(2 * 1.0 * 0.5 / 1.0), 1.0)
    sim = OverdampedHarmonicLangevin(
        tau=float(n), dim=1, N=1, dt=1.0,
        gamma=1.0, kB=1.0, T=0.5,
        F=lambda r: 0.0 * r,
        r0=float(path[0]),
        ab=(a, b),
        seed=0,
    )
    it = iter(incr)
    def fake_noise():
        d = next(it, 0.0)
        sim.noise = np.array([[d]])          # shape (dim, N) = (1, 1)
    sim._set_noise = fake_noise
    sim.generate()
    return sim


class TestEndToEnd:
    def test_solve_matches_track_sojourns(self):
        a, b = -1.0, 1.0
        path = [-2.0, -0.5, 0.5, 1.5, 0.3, -1.5, 0.5, 1.2]
        sim = _run_with_prescribed_path(path, a, b)
        A_ref, B_ref = track_sojourns(path, a, b)
        assert as_lists(sim.A_intervals) == as_lists(A_ref)
        assert as_lists(sim.B_intervals) == as_lists(B_ref)
        # and the exact expected values
        assert as_lists(sim.A_intervals) == [[3, 2]]
        assert as_lists(sim.B_intervals) == [[2]]

    def test_ab_none_leaves_intervals_empty(self):
        sim = OverdampedHarmonicLangevin(
            tau=20.0, dim=1, N=3, dt=0.1,
            T=1.0, F=lambda r: -r, seed=1,
        )
        sim.generate()
        assert sim.A_intervals == [[], [], []]
        assert sim.B_intervals == [[], [], []]

    def test_stochastic_run_structurally_sane(self):
        N = 5
        sim = OverdampedHarmonicLangevin(
            tau=50.0, dim=1, N=N, dt=0.05,
            T=1.0, gamma=1.0, kB=1.0,
            F=lambda r: -r, r0=0.0,
            ab=(-0.5, 0.5), seed=3,
        )
        sim.generate()
        assert len(sim.A_intervals) == N and len(sim.B_intervals) == N
        for walker in list(sim.A_intervals) + list(sim.B_intervals):
            assert isinstance(walker, np.ndarray)
            assert np.all(walker >= 0)

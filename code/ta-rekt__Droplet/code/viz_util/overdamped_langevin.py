import abc
import numpy as np
from typing import Optional, Callable

class Generator(metaclass=abc.ABCMeta):
    """Abstract trajectory generator."""

    def __init__(
        self,
        T: float,
        dim: int = 1,
        N: int = 1,
        dt: float = 1.0,
        seed: Optional[int] = None,
    ):
        self.T = T                # Total time
        self.dim = dim
        self.N = N
        self.dt = dt
        self.n = int(T / dt)
        self.rng = (
            np.random.default_rng(seed) if seed is not None else np.random.default_rng()
        )

    @abc.abstractmethod
    def generate(self):
        ...

class OverdampedHarmonicLangevin(Generator):
    """
    Overdamped Langevin dynamics with an arbitrary spatial force F(r):

        gamma dr/dt = F(r) + sqrt(2 kB T gamma) eta(t),
        <eta(t) eta(t')> = delta(t - t')

    Euler-Maruyama in physical units, optionally with reflecting walls.

    Parameters
    ----------
    tau, dim, N, dt, gamma, kB, T, F, r0, seed : as before.
    r_min, r_max : float, optional
        If given, coordinates are reflected at these walls after each step.
        Use r_min to keep trajectories away from singularities like 1/r.
    ab : (a, b) tuple, optional
        If given, the durations of successive sojourns of coordinate 0 in
        A = {r <= a} and B = {r >= b} (separated by the buffer (a, b)) are
        tracked per walker. After ``generate()`` they are available, as ragged
        lists of N 1-D int arrays, on ``self.A_intervals`` and ``self.B_intervals``.
        See :class:`SojournTracker`.
    store : {"full", "light"}, optional
        "full" (default) keeps the subsampled trajectory of every walker
        (``self.r`` has shape ``(ns, dim, N)``). "light" keeps only walker 0
        (shape ``(ns, dim, 1)``), which avoids the O(ns*N) memory blowup when
        many walkers are used purely for sojourn statistics; the A/B intervals
        are still tracked over all N walkers.
    """

    def __init__(
        self,
        tau: float,
        dim: int = 1,
        N: int = 1,
        nskip: int = 1,
        dt: float = 1.0,
        gamma: float = 1.0,
        kB: float = 1.0,
        T: float = 1.0,
        F: Optional[Callable] = None,
        r0: Optional[np.ndarray] = None,
        r_min: Optional[float] = None,
        r_max: Optional[float] = None,
        ab: Optional[tuple] = None,
        store: str = "full",
        seed: Optional[int] = None,
    ):
        super().__init__(tau, dim, N, dt, seed)
        self.traj_id = "Overdamped Langevin with arbitrary force"

        self.nskip = nskip
        self.ns = (self.n - 1) // nskip + 1          # robust to n % nskip != 0

        self.T_bath = T

        if not callable(F):
            raise ValueError(
                "F is expected to be a callable F(r) -> force, with r and the "
                "return value of shape (dim, N)."
            )
        self.F = F

        self.kB = kB
        self.gamma = gamma
        self.r0 = r0

        self.r_min = r_min
        self.r_max = r_max

        if ab is not None:
            assert ab[0] <= ab[1], "a must be < b"
        
        self.ab = ab

        if store not in ("full", "light"):
            raise ValueError("store must be 'full' or 'light'")
        self.store = store

        self.A_intervals = [[] for i in range(self.N)]
        self.B_intervals = [[] for i in range(self.N)]

        self.shape = (self.ns, self.dim, 1 if store == "light" else self.N)
        self.r = np.empty(self.shape)

        self._set_init_cond()
        # self._set_noise()

    def _set_init_cond(self):
        if self.r0 is None:
            self.r[0] = 0.0

        elif np.shape(self.r0) == (self.dim, self.N) or np.ndim(self.r0) == 0:
            self.r[0] = self.r0
        else:
            raise ValueError(
                "r0 is expected to be a float or an "
                f"array of shape {(self.dim, self.N)}"
            )

    def _set_noise(self):
        self.noise = self.rng.normal(size=(self.dim, self.N))   # was (1, dim, N)

    def _solve(self):
        sqrt_dt = np.sqrt(self.dt)
        diff = np.sqrt(2.0 * self.kB * self.T_bath / self.gamma)
        inv_gamma = 1.0 / self.gamma
        r_min, r_max = self.r_min, self.r_max
        # high-res integration state is always full width (dim, N), independent of
        # the storage mode (which only controls how many walkers self.r keeps).
        r_curr = np.empty((self.dim, self.N))
        r_curr[...] = 0.0 if self.r0 is None else self.r0

        if self.ab is not None:
            a, b = self.ab
            tracker = SojournTracker(a, b, self.N)
            tracker.init_state(r_curr[0])         # classify the initial coordinate

        for i in range(1, self.n):                # i = high-res step index
            self._set_noise()
            r_curr = (
                r_curr
                + self.F(r_curr) * inv_gamma * self.dt
                + diff * sqrt_dt * self.noise
            )
            if r_min is not None:
                r_curr = np.where(r_curr < r_min, 2 * r_min - r_curr, r_curr)
            if r_max is not None:
                r_curr = np.where(r_curr > r_max, 2 * r_max - r_curr, r_curr)

            if self.ab is not None:               # vectorised over all N walkers
                tracker.update(r_curr[0], i)

            if i % self.nskip == 0:
                j = i // self.nskip
                if j < self.ns:
                    self.r[j] = r_curr[:, :1] if self.store == "light" else r_curr

        if self.ab is not None:
            self.A_intervals, self.B_intervals = tracker.as_arrays()

    def generate(self):
        self._solve()


class SojournTracker:
    """
    Vectorised tracker of sojourn durations in two regions of a 1-D coordinate,

        A = {r <= a}    and    B = {r >= b},

    separated by a buffer (a, b), for ``N`` independent walkers advanced in
    lock-step (Schmitt-trigger / hysteresis semantics).  A walker is "assigned"
    to whichever region it last reached: a sojourn in A opens when the walker
    first reaches ``r <= a`` and closes the step it first reaches ``r >= b``
    (and symmetrically for B).  Excursions into the buffer ``a < r < b`` never
    change the assignment.

    Durations (in integration steps) are stored as a *ragged* array — one list
    per walker — appended to only on the rare steps where a walker actually
    crosses, so the per-step cost is the two O(N) boolean masks below rather
    than the O(N) Python loop of the original scalar implementation.

    Convention (kept identical to that scalar implementation): an A-sojourn that
    opened at step ``s`` and is closed by a crossing detected at step ``i`` is
    recorded with duration ``i - s``; the entry step of the first sojourn of
    each walker is taken to be 0.
    """

    def __init__(self, a, b, N):
        assert a <= b, "need a <= b for the buffer [a, b]"
        self.a, self.b, self.N = a, b, N
        self.A_intervals = [[] for _ in range(N)]
        self.B_intervals = [[] for _ in range(N)]
        self.in_A = np.zeros(N, dtype=bool)        # last reached A (awaiting B)
        self.in_B = np.zeros(N, dtype=bool)        # last reached B (awaiting A)
        self.A_entry = np.zeros(N, dtype=np.int64)
        self.B_entry = np.zeros(N, dtype=np.int64)

    def init_state(self, r0):
        """Classify the initial coordinate r0 (shape (N,)) into A vs B."""
        r0 = np.asarray(r0, dtype=float).reshape(self.N)
        self.in_A = r0 <= self.a
        self.in_B = ~self.in_A                      # i.e. r0 > a, as in the original

    def update(self, r0, i):
        """Advance the trackers by one step; r0 = coordinate (shape (N,)) at step i."""
        a, b = self.a, self.b
        r0 = np.asarray(r0).reshape(self.N)

        # A -> B: walkers assigned to A that have now reached r >= b close an A-sojourn.
        cross_AB = (r0 >= b) & self.in_A
        if cross_AB.any():
            for k in np.flatnonzero(cross_AB):
                self.A_intervals[k].append(int(i - self.A_entry[k]))
            self.B_entry[cross_AB] = i
            self.in_A[cross_AB] = False
            self.in_B[cross_AB] = True

        # B -> A: walkers assigned to B that have now reached r <= a close a B-sojourn.
        # Evaluated after the A->B flip above, matching the original two-`if` order.
        cross_BA = (r0 <= a) & self.in_B
        if cross_BA.any():
            for k in np.flatnonzero(cross_BA):
                self.B_intervals[k].append(int(i - self.B_entry[k]))
            self.A_entry[cross_BA] = i
            self.in_B[cross_BA] = False
            self.in_A[cross_BA] = True

    def as_arrays(self):
        """Return (A_intervals, B_intervals) as lists of N 1-D int64 arrays."""
        A = [np.asarray(x, dtype=np.int64) for x in self.A_intervals]
        B = [np.asarray(x, dtype=np.int64) for x in self.B_intervals]
        return A, B


def track_sojourns(positions, a, b):
    """
    Offline / testable driver for :class:`SojournTracker`.

    Parameters
    ----------
    positions : array_like, shape (n, N) or (n,)
        The tracked coordinate for n steps and N walkers (a 1-D input is treated
        as a single walker).  ``positions[0]`` sets the initial A/B assignment;
        steps ``1 .. n-1`` are then fed through the tracker, exactly as ``_solve``
        feeds the live trajectory.
    a, b : float
        Buffer boundaries, ``a <= b``.

    Returns
    -------
    (A_intervals, B_intervals) : tuple of lists of N 1-D int64 arrays
        Per-walker sojourn durations, identical to what ``_solve`` stores on
        ``sim.A_intervals`` / ``sim.B_intervals``.
    """
    positions = np.asarray(positions, dtype=float)
    if positions.ndim == 1:
        positions = positions[:, None]
    n, N = positions.shape

    tr = SojournTracker(a, b, N)
    tr.init_state(positions[0])
    for i in range(1, n):
        tr.update(positions[i], i)
    return tr.as_arrays()

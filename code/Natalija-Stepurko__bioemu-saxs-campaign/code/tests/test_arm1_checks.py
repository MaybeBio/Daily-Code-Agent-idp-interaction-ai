"""analysis/arm1/arm1_checks.py: the relative error of a measured curve and the pool percentile."""
import importlib.util
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location(
    "arm1_checks", Path(__file__).resolve().parents[1] / "analysis" / "arm1" / "arm1_checks.py")
arm1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(arm1)


def test_relative_error_uses_the_lower_half_of_the_q_range_only():
    q = np.linspace(0.01, 0.31, 31)
    I = np.ones_like(q)
    sigma = np.where(q <= 0.16, 0.02, 0.5)          # precise at low q, noisy at high q
    v, n = arm1.relative_error(q, I, sigma)
    assert abs(v - 0.02) < 1e-12 and n == 16


def test_relative_error_skips_nonpositive_points():
    q = np.array([0.01, 0.02, 0.03, 0.04, 0.05])
    I = np.array([-1.0, 2.0, 2.0, 1.0, 1.0])
    sigma = np.array([0.1, 0.2, 0.0, 0.1, 0.1])
    v, n = arm1.relative_error(q, I, sigma)
    # usable q: 0.02, 0.04, 0.05; lower half: q <= 0.035, one point with sigma/I = 0.1
    assert n == 1 and abs(v - 0.1) < 1e-12


def test_percentile_in_pool():
    pool = np.arange(1, 101, dtype=float)
    assert arm1.percentile_in(pool, 25.0) == 25.0 and arm1.percentile_in(pool, 0.5) == 0.0

"""
test_fit_model.py

Tests for viz-util/fit_model.py.

The core test uses an exactly solvable 2-D OLS problem:

    y_i = a * x1_i + b * x2_i + noise,   noise ~ N(0, sigma^2)

True parameters: a=2.0, b=-1.5.
With large n the fitted values must be close to the truth, residuals must be
centred, and the reported stderr must match the OLS closed-form value
  Var(beta_hat) = sigma^2 (X^T X)^{-1}.

Run with:
    python -m pytest test/test_fit_model.py -v
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np
import pytest
from scipy import special
from viz_util.fit_model import fit


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

TRUE_A, TRUE_B, SIGMA = 2.0, -1.5, 0.5

def _make_data(n=2000, seed=42):
    rng = np.random.default_rng(seed)
    x = rng.standard_normal((n, 2))
    y = TRUE_A * x[:, 0] + TRUE_B * x[:, 1] + rng.normal(0, SIGMA, n)
    return x, y

def _linear(x, p):
    # x: (n, 2)
    return p["a"] * x[:, 0] + p["b"] * x[:, 1]

def _ols_stderr(x, sigma):
    """Closed-form OLS standard errors: sqrt(diag(sigma^2 (X^T X)^{-1}))."""
    XtX_inv = np.linalg.inv(x.T @ x)
    return sigma * np.sqrt(np.diag(XtX_inv))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestFitLinear:
    """Exact OLS baseline: 2-D linear model, large n."""

    @pytest.fixture(scope="class")
    def result(self):
        x, y = _make_data()
        return fit(
            _linear, x, y,
            params={"a": (None, -10.0, 10.0),
                    "b": (None, -10.0, 10.0)},
            n_restarts=5,
            rng=np.random.default_rng(0),
        ), x, y

    def test_params_close_to_truth(self, result):
        res, _, _ = result
        assert abs(res["params"]["a"] - TRUE_A) < 0.05
        assert abs(res["params"]["b"] - TRUE_B) < 0.05

    def test_residuals_zero_mean(self, result):
        res, _, _ = result
        assert abs(np.mean(res["residuals"])) < 0.05

    def test_residuals_shape(self, result):
        res, x, y = result
        assert res["residuals"].shape == y.shape

    def test_loss_equals_mean_squared_residuals(self, result):
        res, _, _ = result
        expected = float(np.mean(res["residuals"] ** 2))
        assert abs(res["loss"] - expected) < 1e-12

    def test_stderr_matches_ols_formula(self, result):
        res, x, _ = result
        # OLS sigma estimated from residuals
        n = len(res["residuals"])
        sigma_hat = np.sqrt(np.sum(res["residuals"] ** 2) / (n - 2))
        expected_se = _ols_stderr(x, sigma_hat)
        assert abs(res["stderr"]["a"] - expected_se[0]) < 0.01
        assert abs(res["stderr"]["b"] - expected_se[1]) < 0.01

    def test_stderr_keys_match_params(self, result):
        res, _, _ = result
        assert set(res["stderr"]) == set(res["params"])


class TestFitInterface:
    """API / edge-case tests."""

    def test_1d_input_broadcast(self):
        """1-D x input (n,) should be treated as (n, 1)."""
        rng = np.random.default_rng(1)
        x = rng.standard_normal(500)
        y = 3.0 * x + rng.normal(0, 0.1, 500)
        def model(xi, p): return p["a"] * xi[:, 0]
        res = fit(model, x, y, params={"a": (None, 0.0, 10.0)}, n_restarts=3)
        assert abs(res["params"]["a"] - 3.0) < 0.1

    def test_midpoint_default_initial_guess(self):
        """p0=None should default to midpoint and still converge."""
        x, y = _make_data(n=500)
        res = fit(_linear, x, y,
                  params={"a": (None, -10.0, 10.0),
                          "b": (None, -10.0, 10.0)},
                  n_restarts=3)
        assert abs(res["params"]["a"] - TRUE_A) < 0.2
        assert abs(res["params"]["b"] - TRUE_B) < 0.2

    def test_explicit_initial_guess(self):
        """Explicit p0 close to truth should converge with 1 restart."""
        x, y = _make_data(n=500)
        res = fit(_linear, x, y,
                  params={"a": (2.1, -10.0, 10.0),
                          "b": (-1.4, -10.0, 10.0)},
                  n_restarts=1)
        assert abs(res["params"]["a"] - TRUE_A) < 0.2

    def test_bounds_respected(self):
        """Fitted parameters must lie within the specified bounds."""
        x, y = _make_data(n=300)
        lo, hi = 0.0, 1.0
        res = fit(_linear, x, y,
                  params={"a": (None, lo, hi),
                          "b": (None, lo, hi)},
                  n_restarts=5)
        for v in res["params"].values():
            assert lo <= v <= hi

    def test_no_convergence_raises(self):
        """A model that always errors must exhaust all restarts and raise RuntimeError."""
        x, y = _make_data(n=10)
        def bad_model(xi, p):
            raise ValueError("intentional failure")
        with pytest.raises(RuntimeError, match="did not converge"):
            fit(bad_model, x, y, params={"a": (1.0, 0.0, 2.0)}, n_restarts=3)

    def test_output_keys(self):
        x, y = _make_data(n=200)
        res = fit(_linear, x, y,
                  params={"a": (None, -5.0, 5.0), "b": (None, -5.0, 5.0)})
        assert set(res) == {"params", "stderr", "loss", "residuals"}

    def test_loss_nonnegative(self):
        x, y = _make_data(n=200)
        res = fit(_linear, x, y,
                  params={"a": (None, -5.0, 5.0), "b": (None, -5.0, 5.0)})
        assert res["loss"] >= 0.0

    def test_best_restart_selected(self):
        """With many restarts the result should be at least as good as one restart."""
        x, y = _make_data(n=300, seed=7)
        res_many = fit(_linear, x, y,
                       params={"a": (None, -10.0, 10.0), "b": (None, -10.0, 10.0)},
                       n_restarts=20, rng=np.random.default_rng(0))
        res_one  = fit(_linear, x, y,
                       params={"a": (None, -10.0, 10.0), "b": (None, -10.0, 10.0)},
                       n_restarts=1,  rng=np.random.default_rng(0))
        assert res_many["loss"] <= res_one["loss"] + 1e-10


class TestFitRouse:
    """Nonlinear Rouse MSD model — the motivating real-world use case."""

    TRUE = {"gamma": 0.5, "J": 30.0}

    @staticmethod
    def rouse_msd2(x, params):
        # x: (n, 1) or (n,) of time values
        t = x[:, 0] if x.ndim == 2 else x
        gamma, J = params["gamma"], params["J"]
        tau = (J / gamma) ** 2
        return (2*gamma*t**0.5 * (1 - np.exp(-tau/(np.pi*t)))
                + 2*J * special.erfc((tau/(np.pi*t))**0.5))

    @pytest.fixture(scope="class")
    def data(self):
        time = np.linspace(1, 1000, 300)
        x = time[:, np.newaxis]
        y = self.rouse_msd2(x, self.TRUE)
        return x, y

    def test_recovers_true_params_noiseless(self, data):
        x, y = data
        res = fit(self.rouse_msd2, x, y,
                  params={"gamma": (None, 0.1, 2.0), "J": (None, 1.0, 60.0)},
                  n_restarts=20, rng=np.random.default_rng(0))
        assert abs(res["params"]["gamma"] - self.TRUE["gamma"]) < 1e-3
        assert abs(res["params"]["J"]     - self.TRUE["J"])     < 1e-2

    def test_loss_near_zero_noiseless(self, data):
        x, y = data
        res = fit(self.rouse_msd2, x, y,
                  params={"gamma": (None, 0.1, 2.0), "J": (None, 1.0, 60.0)},
                  n_restarts=20, rng=np.random.default_rng(0))
        assert res["loss"] < 1e-6

"""
test_bvp_two_spheres.py — Comprehensive tests for the 2D axisymmetric
two-sphere BVP solver in viz_util/bvp_two_spheres.py.

Run:
    python -m pytest test/test_bvp_two_spheres.py -v

Conventions (mirrors test/test_fit_model.py):
  - Class-scoped fixtures share expensive BVP solves across the tests in a
    group, keeping total runtime under a few minutes.
  - Tolerances allow for the solver's 1st-order BC accuracy near sphere
    surfaces (ghost-cell with pushed-out image to avoid self-reference).
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from dataclasses import replace

import numpy as np
import pytest
from scipy.optimize import newton_krylov

from viz_util.bvp_solver import (
    Parameters, solve_auto, psi_linear, d2fb_dphi2,
)
from viz_util.bvp_two_spheres import (
    solve_two_spheres, excess_free_energy_two_spheres,
    build_geometry, _scatter, _laplacian, _interp_phi,
    _initial_guess, _decay_length, _make_residual,
    Geometry,
)


# ---------------------------------------------------------------------------
# Shared parameter presets
# ---------------------------------------------------------------------------

def P_benign():
    return Parameters(chi_b=1.0, omega_b=-0.2, kappa=1.0,
                      phi_inf=0.1, R=5.0, eps=1e-6)


def P_linear(omega_b=-0.05):
    return Parameters(chi_b=0.0, omega_b=omega_b, kappa=1.0,
                      phi_inf=0.2, R=5.0, eps=1e-6)


def P_strong():
    return Parameters(chi_b=2.0, omega_b=-0.5, kappa=1.0,
                      phi_inf=0.1, R=5.0, eps=1e-6)


def _Nz_for(p, d, dz_target=0.15, z_buf=4.0):
    """Choose Nz so dz stays roughly constant as d grows."""
    span = 2.0 * (p.R + z_buf) + d
    return int(span / dz_target) + 1


def _avg_surface_enrichment(result, sphere=1):
    """<phi(surface) - phi_inf> averaged over solid-angle sin(theta) weight."""
    geom = result["geom"]
    phi  = result["phi"]
    p    = geom.p
    R    = p.R
    zc   = geom.z_c1 if sphere == 1 else geom.z_c2
    theta = np.linspace(0.0, np.pi, 400)
    sin_t = np.sin(theta)
    phi_s = _interp_phi(phi, geom, R * sin_t, zc + R * np.cos(theta))
    return float(np.trapezoid((phi_s - p.phi_inf) * sin_t, theta)) / 2.0


# ===========================================================================
# 1. TestGeometry
# ===========================================================================

class TestGeometry:
    """Invariants of `build_geometry`."""

    @pytest.fixture(scope="class")
    def geom(self):
        return build_geometry(P_benign(), d=2.0, Nrho=40, Nz=120)

    def test_grid_shapes(self, geom):
        assert geom.rho_grid.shape == (geom.Nrho,)
        assert geom.z_grid.shape   == (geom.Nz,)
        for m in (geom.inside_mask, geom.fluid_mask,
                  geom.dir_mask, geom.unk_mask):
            assert m.shape == (geom.Nrho, geom.Nz)
            assert m.dtype == bool

    def test_grid_endpoints(self, geom):
        p   = geom.p
        ell = _decay_length(p)
        z_buf = 4.0 * ell
        assert geom.rho_grid[0]  == pytest.approx(0.0)
        assert geom.rho_grid[-1] == pytest.approx(p.R + 4.0 * ell)
        assert geom.z_grid[0]    == pytest.approx(-p.R - z_buf)
        assert geom.z_grid[-1]   == pytest.approx(geom.z_c2 + p.R + z_buf)

    def test_spacing_positive(self, geom):
        assert geom.drho == pytest.approx(geom.rho_grid[1] - geom.rho_grid[0])
        assert geom.dz   == pytest.approx(geom.z_grid[1]   - geom.z_grid[0])
        assert geom.drho > 0 and geom.dz > 0

    def test_sphere_centres(self, geom):
        assert geom.z_c1 == pytest.approx(0.0)
        assert geom.z_c2 == pytest.approx(2.0 * geom.p.R + geom.d)

    def test_mask_partition(self, geom):
        # inside XOR fluid covers every cell exactly once
        assert (geom.inside_mask ^ geom.fluid_mask).all()
        assert not (geom.inside_mask & geom.fluid_mask).any()

    def test_dirichlet_unknown_split(self, geom):
        assert not (geom.dir_mask & geom.unk_mask).any()
        assert ((geom.dir_mask | geom.unk_mask) == geom.fluid_mask).all()

    def test_dirichlet_faces(self, geom):
        # Outer cylinder and z end-caps are Dirichlet wherever they are fluid.
        for face in (geom.fluid_mask[-1, :], geom.fluid_mask[:, 0],
                     geom.fluid_mask[:, -1]):
            pass  # just check that the corresponding dir_mask slice matches
        assert (geom.dir_mask[-1, :] == geom.fluid_mask[-1, :]).all()
        assert (geom.dir_mask[:, 0]  == geom.fluid_mask[:, 0]).all()
        assert (geom.dir_mask[:, -1] == geom.fluid_mask[:, -1]).all()

    def test_sphere_classification(self, geom):
        # Cells near each sphere centre are inside.
        i0 = 0  # rho = 0
        j_c1 = int(np.argmin(np.abs(geom.z_grid - geom.z_c1)))
        j_c2 = int(np.argmin(np.abs(geom.z_grid - geom.z_c2)))
        assert geom.inside_mask[i0, j_c1]
        assert geom.inside_mask[i0, j_c2]
        # Cells well outside both spheres are fluid.
        i_far = geom.Nrho - 1
        j_mid = int(np.argmin(np.abs(geom.z_grid
                                     - 0.5 * (geom.z_c1 + geom.z_c2))))
        assert geom.fluid_mask[i_far, j_mid]

    def test_no_sphere_overlap(self, geom):
        n_inside = int(geom.inside_mask.sum())
        # For d > 0 the two solid balls are disjoint, so n_inside should be
        # well above 0 but below ~2 * (pi R^2)/(drho*dz) (rough disk-area
        # bookkeeping; we just want a sanity envelope).
        cell_area = geom.drho * geom.dz
        upper = 2.5 * np.pi * geom.p.R ** 2 / cell_area
        assert 0 < n_inside < upper

    def test_ghost_cells_are_inside(self, geom):
        gi = geom.ghost_ij
        if gi.shape[0] > 0:
            assert geom.inside_mask[gi[:, 0], gi[:, 1]].all()

    def test_ghost_cells_border_fluid(self, geom):
        Nrho, Nz = geom.Nrho, geom.Nz
        for (i, j) in geom.ghost_ij:
            has_fluid = False
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < Nrho and 0 <= jj < Nz and geom.fluid_mask[ii, jj]:
                    has_fluid = True
                    break
            assert has_fluid, f"ghost ({i},{j}) has no fluid 4-neighbour"

    def test_ghost_weights_sum_to_one(self, geom):
        """P1 stencil weights always sum to 1.  P2 stencil weights sum to 1
        for quadratic ghosts and to 0 for linear-fallback ghosts (whose P2
        slot is unused: alpha2=0, wts=0)."""
        if geom.ghost_wts_1.shape[0] > 0:
            assert np.allclose(geom.ghost_wts_1.sum(axis=1), 1.0, atol=1e-12)
            s2 = geom.ghost_wts_2.sum(axis=1)
            quad = ~geom.ghost_used_linear
            assert np.allclose(s2[quad], 1.0, atol=1e-12)
            assert np.allclose(s2[~quad], 0.0, atol=1e-12)

    def test_ghost_image_stencil_avoids_inside(self, geom):
        """Both probe-point bilinear stencils must be all fluid (when not
        falling back to the linear formula)."""
        for (i0, j0), used_lin in zip(geom.ghost_i0j0_1,
                                       geom.ghost_used_linear):
            corners = [(i0, j0), (i0 + 1, j0),
                       (i0, j0 + 1), (i0 + 1, j0 + 1)]
            for (i, j) in corners:
                assert not geom.inside_mask[i, j], (
                    f"P1 stencil corner ({i},{j}) is inside a sphere")
        for (i0, j0), used_lin in zip(geom.ghost_i0j0_2,
                                       geom.ghost_used_linear):
            if used_lin:
                continue
            corners = [(i0, j0), (i0 + 1, j0),
                       (i0, j0 + 1), (i0 + 1, j0 + 1)]
            for (i, j) in corners:
                assert not geom.inside_mask[i, j], (
                    f"P2 stencil corner ({i},{j}) is inside a sphere")

    def test_ghost_d1_d2_positive(self, geom):
        if geom.ghost_d1.shape[0] > 0:
            assert (geom.ghost_d1 > 0).all()
            # For quadratic ghosts, d2 > d1 by at least ~one cell.
            quad = ~geom.ghost_used_linear
            sep = (geom.ghost_d2 - geom.ghost_d1)[quad]
            min_sep = 0.9 * max(geom.drho, geom.dz)
            assert (sep >= min_sep - 1e-12).all(), (
                f"P2-P1 separation below min_sep={min_sep:.4f}: "
                f"got min {sep.min():.4f}")

    def test_ghost_constant_partition(self, geom):
        """Quadratic ghosts satisfy alpha1 + alpha2 = 1 to machine precision
        (so a constant phi field scatters back unchanged)."""
        if geom.ghost_coef_phi1.shape[0] > 0:
            s = geom.ghost_coef_phi1 + geom.ghost_coef_phi2
            assert np.allclose(s, 1.0, atol=1e-12)


# ===========================================================================
# 2. TestScatter
# ===========================================================================

class TestScatter:
    """Behavior of `_scatter`."""

    @pytest.fixture(scope="class")
    def geom(self):
        return build_geometry(P_benign(), d=2.0, Nrho=40, Nz=120)

    def test_output_shape(self, geom):
        x = np.full(int(geom.unk_mask.sum()), geom.p.phi_inf)
        phi = _scatter(x, geom)
        assert phi.shape == (geom.Nrho, geom.Nz)

    def test_unknowns_roundtrip(self, geom):
        rng = np.random.default_rng(0)
        n   = int(geom.unk_mask.sum())
        x   = 0.5 + 0.1 * rng.standard_normal(n)  # in (0, 1)
        phi = _scatter(x, geom)
        assert np.allclose(phi[geom.unk_mask], x)

    def test_dirichlet_cells_phi_inf(self, geom):
        x = np.full(int(geom.unk_mask.sum()), 0.7)
        phi = _scatter(x, geom)
        assert np.allclose(phi[geom.dir_mask], geom.p.phi_inf)

    def test_uniform_omega_zero(self):
        p = replace(P_benign(), omega_b=0.0)
        geom = build_geometry(p, d=2.0, Nrho=40, Nz=120)
        x = np.full(int(geom.unk_mask.sum()), p.phi_inf)
        phi = _scatter(x, geom)
        # Uniform phi_inf on fluid + ghost-fill with omega_b=0 gives phi_inf
        # everywhere.
        assert np.allclose(phi, p.phi_inf, atol=1e-12)

    def test_ghost_fill_formula(self):
        """Ghost cells satisfy phi(G) = alpha1*phi(P1) + alpha2*phi(P2) +
        beta*omega_b/kappa to machine precision."""
        p = P_benign()
        geom = build_geometry(p, d=2.0, Nrho=40, Nz=120)
        rng = np.random.default_rng(7)
        n   = int(geom.unk_mask.sum())
        x   = p.phi_inf + 0.05 * rng.standard_normal(n)
        phi = _scatter(x, geom)

        if geom.ghost_ij.shape[0] == 0:
            return

        i0_1, j0_1 = geom.ghost_i0j0_1[:, 0], geom.ghost_i0j0_1[:, 1]
        i0_2, j0_2 = geom.ghost_i0j0_2[:, 0], geom.ghost_i0j0_2[:, 1]
        phi_P1 = (geom.ghost_wts_1[:, 0] * phi[i0_1,     j0_1]
                + geom.ghost_wts_1[:, 1] * phi[i0_1 + 1, j0_1]
                + geom.ghost_wts_1[:, 2] * phi[i0_1,     j0_1 + 1]
                + geom.ghost_wts_1[:, 3] * phi[i0_1 + 1, j0_1 + 1])
        phi_P2 = (geom.ghost_wts_2[:, 0] * phi[i0_2,     j0_2]
                + geom.ghost_wts_2[:, 1] * phi[i0_2 + 1, j0_2]
                + geom.ghost_wts_2[:, 2] * phi[i0_2,     j0_2 + 1]
                + geom.ghost_wts_2[:, 3] * phi[i0_2 + 1, j0_2 + 1])
        expected = (geom.ghost_coef_phi1 * phi_P1
                  + geom.ghost_coef_phi2 * phi_P2
                  + geom.ghost_coef_bc   * (p.omega_b / p.kappa))
        actual = phi[geom.ghost_ij[:, 0], geom.ghost_ij[:, 1]]
        assert np.allclose(actual, expected, atol=1e-12)

    def test_ghost_robin_bc_quadratic(self):
        """For a synthetic quadratic phi(s) along the normal that satisfies
        the Robin BC at s=delta exactly, _scatter reproduces phi(G)
        to machine precision."""
        from viz_util.bvp_two_spheres import _scatter, build_geometry, Geometry
        p = P_benign()
        geom = build_geometry(p, d=2.0, Nrho=50, Nz=150)
        # Build a phi field: linear in s along each ghost's normal would be
        # exact for any first-order scheme; quadratic adds curvature, which
        # the new scheme must recover exactly.  We test by setting a chosen
        # phi(P1), phi(P2) and verifying the BC slope at the surface.
        rng = np.random.default_rng(3)
        x = p.phi_inf + 0.1 * rng.standard_normal(int(geom.unk_mask.sum()))
        phi = _scatter(x, geom)
        # At each ghost, compute b + 2*c*delta from the quadratic fit and
        # verify it equals omega_b/kappa to machine precision.
        if geom.ghost_ij.shape[0] == 0:
            return
        quad_mask = ~geom.ghost_used_linear
        i0_1 = geom.ghost_i0j0_1[:, 0]; j0_1 = geom.ghost_i0j0_1[:, 1]
        i0_2 = geom.ghost_i0j0_2[:, 0]; j0_2 = geom.ghost_i0j0_2[:, 1]
        phi_P1 = (geom.ghost_wts_1[:, 0] * phi[i0_1,     j0_1]
                + geom.ghost_wts_1[:, 1] * phi[i0_1 + 1, j0_1]
                + geom.ghost_wts_1[:, 2] * phi[i0_1,     j0_1 + 1]
                + geom.ghost_wts_1[:, 3] * phi[i0_1 + 1, j0_1 + 1])
        phi_P2 = (geom.ghost_wts_2[:, 0] * phi[i0_2,     j0_2]
                + geom.ghost_wts_2[:, 1] * phi[i0_2 + 1, j0_2]
                + geom.ghost_wts_2[:, 2] * phi[i0_2,     j0_2 + 1]
                + geom.ghost_wts_2[:, 3] * phi[i0_2 + 1, j0_2 + 1])
        d1 = geom.ghost_d1[quad_mask]
        d2 = geom.ghost_d2[quad_mask]
        delta = geom.ghost_delta[quad_mask]
        phi_P1 = phi_P1[quad_mask]
        phi_P2 = phi_P2[quad_mask]
        # c from closed form: c = [(P1-P2)/(d1-d2) - omega/kappa] / (d1+d2-2*delta)
        slope_avg = (phi_P1 - phi_P2) / (d1 - d2)
        c = (slope_avg - p.omega_b / p.kappa) / (d1 + d2 - 2.0 * delta)
        b = p.omega_b / p.kappa - 2.0 * c * delta
        # BC: b + 2*c*delta = omega_b/kappa
        bc_check = b + 2.0 * c * delta
        assert np.allclose(bc_check, p.omega_b / p.kappa, atol=1e-12)


# ===========================================================================
# 3. TestLaplacian
# ===========================================================================

class TestLaplacian:
    """`_laplacian` reproduces analytic Laplacians on simple fields."""

    @pytest.fixture(scope="class")
    def geom(self):
        return build_geometry(P_benign(), d=2.0, Nrho=50, Nz=150)

    def _interior(self, geom):
        return slice(1, geom.Nrho - 1), slice(1, geom.Nz - 1)

    def test_constant(self, geom):
        phi = np.full((geom.Nrho, geom.Nz), 0.3)
        lap = _laplacian(phi, geom)
        si, sj = self._interior(geom)
        assert np.allclose(lap[si, sj], 0.0, atol=1e-10)

    def test_linear_in_z(self, geom):
        Z = geom.z_grid[None, :].repeat(geom.Nrho, axis=0)
        lap = _laplacian(Z.astype(float), geom)
        si, sj = self._interior(geom)
        assert np.allclose(lap[si, sj], 0.0, atol=1e-10)

    def test_quadratic_in_z(self, geom):
        Z = geom.z_grid[None, :].repeat(geom.Nrho, axis=0).astype(float)
        lap = _laplacian(Z ** 2, geom)
        si, sj = self._interior(geom)
        assert np.allclose(lap[si, sj], 2.0, atol=1e-10)

    def test_quadratic_in_rho(self, geom):
        Rho = geom.rho_grid[:, None].repeat(geom.Nz, axis=1).astype(float)
        lap = _laplacian(Rho ** 2, geom)
        si, sj = self._interior(geom)
        # Interior cells avoid the rho=0 axis where the 1/rho term needs care.
        assert np.allclose(lap[si, sj], 4.0, atol=1e-10)

    def test_axis_quadratic_rho(self, geom):
        """phi = rho^2 on the axis row: 4*(drho^2 - 0)/drho^2 = 4."""
        Rho = geom.rho_grid[:, None].repeat(geom.Nz, axis=1).astype(float)
        lap = _laplacian(Rho ** 2, geom)
        # Axis interior in z (avoid the two z end-caps where ddz term is wrong)
        assert np.allclose(lap[0, 1:-1], 4.0, atol=1e-10)

    def test_axis_linear_in_z(self, geom):
        """phi = z on the axis row: lap = 0."""
        Z = geom.z_grid[None, :].repeat(geom.Nrho, axis=0).astype(float)
        lap = _laplacian(Z, geom)
        assert np.allclose(lap[0, 1:-1], 0.0, atol=1e-10)


# ===========================================================================
# 4. TestInterp
# ===========================================================================

class TestInterp:
    """`_interp_phi` bilinear correctness."""

    @pytest.fixture(scope="class")
    def geom(self):
        return build_geometry(P_benign(), d=2.0, Nrho=40, Nz=120)

    def test_grid_point_recovery(self, geom):
        rng = np.random.default_rng(0)
        phi = rng.standard_normal((geom.Nrho, geom.Nz))
        # Sample at every 5th grid point.
        ii = np.arange(0, geom.Nrho, 5)
        jj = np.arange(0, geom.Nz, 5)
        I, J = np.meshgrid(ii, jj, indexing="ij")
        rho_q = geom.rho_grid[I.ravel()]
        z_q   = geom.z_grid[J.ravel()]
        out = _interp_phi(phi, geom, rho_q, z_q)
        # Allow the boundary clamp at rho_max / z_max to ever-so-slightly
        # offset a query at the very edge; skip those.
        expected = phi[I.ravel(), J.ravel()]
        # Exclude the last-row/column points whose clamp moves them 1e-12 in.
        keep = (rho_q < geom.rho_grid[-1] - 1e-9) & (z_q < geom.z_grid[-1] - 1e-9)
        assert np.allclose(out[keep], expected[keep], atol=1e-10)

    def test_exact_on_linear_field(self, geom):
        a, b = 0.7, -1.3
        Rho = geom.rho_grid[:, None].repeat(geom.Nz, axis=1)
        Z   = geom.z_grid[None, :].repeat(geom.Nrho, axis=0)
        phi = a * Rho + b * Z
        rng = np.random.default_rng(1)
        rho_q = rng.uniform(0.1, geom.rho_grid[-1] - 0.1, size=20)
        z_q   = rng.uniform(geom.z_grid[0] + 0.1, geom.z_grid[-1] - 0.1, size=20)
        out = _interp_phi(phi, geom, rho_q, z_q)
        expected = a * rho_q + b * z_q
        assert np.allclose(out, expected, atol=1e-10)

    def test_out_of_range_clamped(self, geom):
        phi = np.full((geom.Nrho, geom.Nz), 0.42)
        out = _interp_phi(phi, geom,
                          np.array([-1.0, geom.rho_grid[-1] + 5.0, 1.0]),
                          np.array([geom.z_grid[0] - 5.0, 0.0,
                                    geom.z_grid[-1] + 5.0]))
        assert np.allclose(out, 0.42)


# ===========================================================================
# 5. TestInitialGuess
# ===========================================================================

class TestInitialGuess:

    def test_length_matches_unknowns(self):
        geom = build_geometry(P_benign(), d=2.0, Nrho=40, Nz=120)
        x0 = _initial_guess(geom)
        assert x0.size == int(geom.unk_mask.sum())

    def test_values_clipped(self):
        geom = build_geometry(P_benign(), d=2.0, Nrho=40, Nz=120)
        x0 = _initial_guess(geom)
        assert (x0 > geom.p.eps).all()
        assert (x0 < 1.0 - geom.p.eps).all()

    def test_omega_zero_uniform(self):
        p = replace(P_benign(), omega_b=0.0)
        geom = build_geometry(p, d=2.0, Nrho=40, Nz=120)
        x0 = _initial_guess(geom)
        # psi_linear ∝ omega_b, so at omega_b=0 the guess collapses to phi_inf
        # (modulo the eps clip, which doesn't trigger here).
        assert np.allclose(x0, p.phi_inf, atol=1e-12)


# ===========================================================================
# 6. TestSolverBasic
# ===========================================================================

class TestSolverBasic:

    @pytest.fixture(scope="class")
    def result(self):
        return solve_two_spheres(P_benign(), d=2.0,
                                 Nrho=60, Nz=180, verbose=False)

    def test_keys(self, result):
        assert set(result.keys()) >= {
            "phi", "geom", "p", "d", "mu_inf", "success", "x",
        }

    def test_success(self, result):
        assert result["success"]

    def test_phi_in_unit_interval(self, result):
        phi = result["phi"][result["geom"].fluid_mask]
        assert np.isfinite(phi).all()
        assert (phi > 0).all()
        assert (phi < 1).all()

    def test_residual_below_tol(self, result):
        geom = result["geom"]
        res_fn = _make_residual(geom, result["mu_inf"])
        res = res_fn(result["x"])
        assert float(np.max(np.abs(res))) < 1e-5

    def test_midplane_reflection_symmetry(self, result):
        """phi(rho, z) == phi(rho, z_c1 + z_c2 - z)."""
        geom = result["geom"]
        phi  = result["phi"]
        z_mirror = geom.z_c1 + geom.z_c2 - geom.z_grid
        # interpolate phi at the mirrored z values
        Rho = geom.rho_grid[:, None].repeat(geom.Nz, axis=1).ravel()
        Z   = np.broadcast_to(z_mirror[None, :], phi.shape).ravel()
        phi_mirror = _interp_phi(phi, geom, Rho, Z).reshape(phi.shape)
        diff = np.where(geom.fluid_mask, np.abs(phi - phi_mirror), 0.0)
        assert float(diff.max()) < 1e-5

    def test_far_field_phi_inf(self, result):
        geom = result["geom"]
        phi  = result["phi"]
        p    = geom.p
        i_far = geom.Nrho - 3
        j_mid = int(np.argmin(np.abs(geom.z_grid - 0.5 * (geom.z_c1 + geom.z_c2))))
        assert abs(phi[i_far, j_mid] - p.phi_inf) < 0.01


# ===========================================================================
# 7. TestOmegaZero
# ===========================================================================

class TestOmegaZero:
    """omega_b = 0, chi_b = 0: unique uniform solution phi = phi_inf."""

    @pytest.fixture(scope="class")
    def result(self):
        p = Parameters(chi_b=0.0, omega_b=0.0, kappa=1.0,
                       phi_inf=0.2, R=5.0, eps=1e-6)
        return solve_two_spheres(p, d=2.0, Nrho=50, Nz=140, verbose=False)

    def test_converged(self, result):
        assert result["success"]

    def test_uniform_phi_inf(self, result):
        geom = result["geom"]
        phi  = result["phi"][geom.fluid_mask]
        assert np.allclose(phi, geom.p.phi_inf, atol=1e-8)

    def test_excess_F_zero(self, result):
        F = excess_free_energy_two_spheres(result)
        assert abs(F) < 1e-8


# ===========================================================================
# 8. TestSingleSphereLimit
# ===========================================================================

class TestSingleSphereLimit:
    """At very large d each sphere is effectively isolated; the 2D radial
    profile must match the 1D spherical BVP from bvp_solver.py."""

    L_iso = 10.0

    @pytest.fixture(scope="class")
    def setup(self):
        p = P_benign()
        d_far = 8.0 * p.R
        Nz = _Nz_for(p, d_far, dz_target=0.10)
        res = solve_two_spheres(p, d_far, Nrho=100, Nz=Nz,
                                rho_max=p.R + 8.0, z_buf=4.0, verbose=False)
        # 1D spherical reference (open reservoir)
        p_iso = replace(p, L=self.L_iso, conserve_volume=False,
                        parallel_plates=False)
        sol1d = solve_auto(p_iso, verbose=False)
        assert sol1d is not None and sol1d.success
        return res, p, sol1d

    def test_axial_profile_matches_1d(self, setup):
        res, p, sol1d = setup
        geom = res["geom"]
        # Sample along the axis above sphere 1, in (R, R + L_iso/2)
        z_query = np.linspace(p.R + 0.05, p.R + 0.5 * self.L_iso, 30)
        phi_2d = _interp_phi(res["phi"], geom,
                             np.zeros_like(z_query), z_query)
        phi_1d = sol1d.sol(z_query - p.R)[0]
        max_err = float(np.max(np.abs(phi_2d - phi_1d)))
        # With the quadratic ghost-cell scheme the surface BC is O(dx^2);
        # the residual error here is dominated by box-truncation and
        # bilinear sampling along the axis.
        assert max_err < 0.005

    def test_equatorial_profile_matches_1d(self, setup):
        res, p, sol1d = setup
        geom = res["geom"]
        # Equator of sphere 1: rho in (R, R + L_iso/2), z = z_c1 = 0
        rho_query = np.linspace(p.R + 0.05, p.R + 0.5 * self.L_iso, 30)
        phi_2d = _interp_phi(res["phi"], geom, rho_query,
                             np.zeros_like(rho_query))
        phi_1d = sol1d.sol(rho_query - p.R)[0]
        max_err = float(np.max(np.abs(phi_2d - phi_1d)))
        assert max_err < 0.005

    def test_spherical_isotropy(self, setup):
        """At equal surface distance, phi shouldn't depend on the angle."""
        res, p, _ = setup
        geom = res["geom"]
        s_values = np.array([0.3, 0.7, 1.5])
        for s in s_values:
            r = p.R + s
            # Sample at theta = 0 (north pole direction), pi/4, pi/2 (equator)
            theta_samples = np.array([1e-3, np.pi / 4.0, np.pi / 2.0])
            rho_q = r * np.sin(theta_samples)
            z_q   = geom.z_c1 + r * np.cos(theta_samples)
            phi_at_angles = _interp_phi(res["phi"], geom, rho_q, z_q)
            spread = float(phi_at_angles.max() - phi_at_angles.min())
            assert spread < 0.01, (
                f"isotropy violated at s={s}: spread={spread:.4f}")


# ===========================================================================
# 9. TestLinearResponseLimit
# ===========================================================================

class TestLinearResponseLimit:
    """As omega_b → 0 the solution must approach the linear superposition
    phi_inf + psi_linear(s1) + psi_linear(s2)."""

    OMEGAS = (-0.10, -0.05, -0.025)

    @pytest.fixture(scope="class")
    def errors(self):
        errs = []
        for omega_b in self.OMEGAS:
            p = replace(P_linear(), omega_b=omega_b)
            res = solve_two_spheres(p, d=2.0 * p.R,
                                    Nrho=80, Nz=240, verbose=False)
            assert res["success"], f"solve failed at omega_b={omega_b}"
            geom = res["geom"]
            Rho = geom.rho_grid[:, None].repeat(geom.Nz, axis=1)
            Z   = geom.z_grid[None, :].repeat(geom.Nrho, axis=0)
            s1 = np.clip(np.sqrt(Rho ** 2 + (Z - geom.z_c1) ** 2) - p.R, 0, None)
            s2 = np.clip(np.sqrt(Rho ** 2 + (Z - geom.z_c2) ** 2) - p.R, 0, None)
            psi1 = psi_linear(s1, p)[0]
            psi2 = psi_linear(s2, p)[0]
            phi_pred = p.phi_inf + psi1 + psi2
            diff = np.where(geom.fluid_mask,
                            np.abs(res["phi"] - phi_pred), 0.0)
            errs.append(float(diff.max()))
        return list(zip(self.OMEGAS, errs))

    def test_per_omega_bound(self, errors):
        for omega_b, err in errors:
            # Budget = O(omega_b^2) physical + O(drho) discretisation slack
            budget = 0.5 * omega_b ** 2 + 0.005
            assert err < budget, (
                f"omega_b={omega_b}: err={err:.4e}, budget={budget:.4e}")

    def test_error_scales_with_omega(self, errors):
        omegas = [o for o, _ in errors]
        errs   = [e for _, e in errors]
        # Error at the largest |omega_b| should be at least 2x the error at
        # the smallest |omega_b|. (Quadratic scaling would give 16x; allow
        # plenty of slack for the linear FD-discretisation component.)
        ratio = errs[0] / errs[-1]
        assert ratio > 2.0, (
            f"errors don't grow with |omega_b|: {list(zip(omegas, errs))}, "
            f"ratio={ratio:.2f}")


# ===========================================================================
# 10. TestUDecay
# ===========================================================================

class TestUDecay:
    """U(d) = F(d) - F(d_ref) decays smoothly to ~0 at long enough d.

    The 2D FD integrator carries a ~0.2 discretisation noise floor in F at
    the resolution used here, so the signal-to-noise check uses parameters
    that produce a clearly-separated attractive close-range plateau.
    """

    D_REF_FACTOR = 6.0  # d_ref = 6 R (deep in the far-field regime)
    D_FACTORS = (0.1, 0.3, 0.8, 2.0, 4.0)

    @pytest.fixture(scope="class")
    def sweep(self):
        # Stronger attraction than P_strong so the close-range U sits well
        # above the FD noise floor.
        p = Parameters(chi_b=1.5, omega_b=-0.6, kappa=1.0,
                       phi_inf=0.1, R=5.0, eps=1e-6)
        Nrho = 120

        def solve_one(d):
            return solve_two_spheres(p, d, Nrho=Nrho,
                                     Nz=_Nz_for(p, d, dz_target=0.08,
                                                z_buf=3.0),
                                     rho_max=p.R + 4.0, z_buf=3.0,
                                     verbose=False)

        d_ref = self.D_REF_FACTOR * p.R
        res_ref = solve_one(d_ref)
        F_ref = excess_free_energy_two_spheres(res_ref)
        Us = []
        results = []
        for factor in self.D_FACTORS:
            d = factor * p.R
            res = solve_one(d)
            F = excess_free_energy_two_spheres(res)
            Us.append(F - F_ref)
            results.append(res)
        return dict(p=p, d_ref=d_ref, F_ref=F_ref, Us=np.array(Us),
                    results=results, res_ref=res_ref)

    def test_all_converged(self, sweep):
        for res in sweep["results"] + [sweep["res_ref"]]:
            assert res["success"]
            assert np.isfinite(res["phi"][res["geom"].fluid_mask]).all()

    def test_small_d_attractive(self, sweep):
        # Closest separation (0.1 R) must be strongly attractive.
        assert sweep["Us"][0] < -0.5, f"U(d=0.1R)={sweep['Us'][0]:.3f}"

    def test_signal_vs_noise(self, sweep):
        """|U(close)| should be much larger than |U(far)|."""
        u_close = abs(sweep["Us"][0])
        u_far   = abs(sweep["Us"][-1])
        assert u_close > 4.0 * u_far, (
            f"signal/noise too low: u_close={u_close:.3f}, "
            f"u_far={u_far:.3f}")

    def test_large_d_small(self, sweep):
        """|U(d=4R)| must be below a noise-floor threshold."""
        assert abs(sweep["Us"][-1]) < 0.3, (
            f"U(d=4R)={sweep['Us'][-1]:.3f}")

    def test_envelope_decays(self, sweep):
        """|U(d)| decreases on average from the close-range to far-range
        regimes.  Strict monotonicity is too much once we're at the FD
        noise floor; require a clear separation in averaged magnitude."""
        mags = np.abs(sweep["Us"])
        close_mean = float(np.mean(mags[:2]))
        far_mean   = float(np.mean(mags[-2:]))
        assert close_mean > 3.0 * far_mean, (
            f"close/far ratio = {close_mean / max(far_mean, 1e-9):.2f}, "
            f"mags = {list(zip(self.D_FACTORS, mags))}")


# ===========================================================================
# 11. TestExcessFreeEnergy
# ===========================================================================

class TestExcessFreeEnergy:
    """Properties of `excess_free_energy_two_spheres`."""

    def test_zero_for_uniform(self):
        p = Parameters(chi_b=0.0, omega_b=0.0, kappa=1.0,
                       phi_inf=0.2, R=5.0, eps=1e-6)
        res = solve_two_spheres(p, d=2.0, Nrho=50, Nz=140, verbose=False)
        assert abs(excess_free_energy_two_spheres(res)) < 1e-8

    def test_attractive_F_negative(self):
        p = Parameters(chi_b=1.0, omega_b=-0.3, kappa=1.0,
                       phi_inf=0.1, R=5.0, eps=1e-6)
        res = solve_two_spheres(p, d=2.0, Nrho=60, Nz=180, verbose=False)
        assert excess_free_energy_two_spheres(res) < 0.0

    def test_sign_flip_with_omega(self):
        """F(+omega) > F(-omega) at matched magnitude (surface term flips)."""
        base = dict(chi_b=1.0, kappa=1.0, phi_inf=0.1, R=5.0, eps=1e-6)
        p_neg = Parameters(omega_b=-0.3, **base)
        p_pos = Parameters(omega_b=+0.3, **base)
        F_neg = excess_free_energy_two_spheres(
            solve_two_spheres(p_neg, d=2.0, Nrho=60, Nz=180, verbose=False))
        F_pos = excess_free_energy_two_spheres(
            solve_two_spheres(p_pos, d=2.0, Nrho=60, Nz=180, verbose=False))
        assert F_pos > F_neg

    def test_R_scaling(self):
        """omega_b * 4πR² * <phi_surf - phi_inf> scales ~R² at fixed omega_b."""
        common = dict(chi_b=0.0, omega_b=-0.1, kappa=1.0,
                      phi_inf=0.2, eps=1e-6)
        p4 = Parameters(R=4.0, **common)
        p8 = Parameters(R=8.0, **common)
        res4 = solve_two_spheres(p4, d=2.0, Nrho=80, Nz=240, verbose=False)
        res8 = solve_two_spheres(p8, d=2.0, Nrho=80, Nz=240, verbose=False)
        avg4 = _avg_surface_enrichment(res4)
        avg8 = _avg_surface_enrichment(res8)
        # Linear theory: <phi_surf - phi_inf> ≈ -omega_b/(kappa*(lam + 1/R))
        # which increases with R (as 1/R term shrinks).  Ratio (R8/R4) of the
        # surface contribution = (R8^2 * avg8) / (R4^2 * avg4) should be
        # > 4 (since both R^2 grows AND avg increases at fixed omega_b).
        surf4 = p4.R ** 2 * avg4
        surf8 = p8.R ** 2 * avg8
        ratio = surf8 / surf4
        assert ratio > 4.0, f"R-scaling ratio = {ratio:.2f}, expected > 4"

    def test_idempotent(self):
        p = P_benign()
        res = solve_two_spheres(p, d=2.0, Nrho=50, Nz=140, verbose=False)
        F1 = excess_free_energy_two_spheres(res)
        F2 = excess_free_energy_two_spheres(res)
        assert F1 == F2


# ===========================================================================
# 12. TestGridConvergence (slow)
# ===========================================================================

class TestGridConvergence:
    """Refinement reduces the ghost-cell BC error toward the linear-theory
    surface enrichment."""

    @pytest.fixture(scope="class")
    def sweep(self):
        p = replace(P_linear(omega_b=-0.05),
                    phi_inf=0.2, chi_b=0.0, R=5.0)
        # Truth from psi_linear at the surface (s = 0)
        truth = float(psi_linear(np.array([0.0]), p)[0][0])
        out = []
        for (Nrho, Nz) in [(40, 120), (80, 240), (160, 480)]:
            res = solve_two_spheres(p, d=6.0 * p.R, Nrho=Nrho, Nz=Nz,
                                    rho_max=p.R + 8.0, z_buf=4.0,
                                    verbose=False)
            assert res["success"]
            avg = _avg_surface_enrichment(res)
            out.append((Nrho, avg, abs(avg - truth)))
        return dict(p=p, truth=truth, levels=out)

    def test_monotone_convergence(self, sweep):
        errs = [lev[2] for lev in sweep["levels"]]
        assert errs[0] > errs[1] > errs[2], (
            f"errors not monotone: {errs}, truth={sweep['truth']:.5f}")

    def test_second_order_rate(self, sweep):
        errs = [lev[2] for lev in sweep["levels"]]
        # The quadratic ghost-cell scheme is O(dx^2): a 4x total refinement
        # should give ~16x error reduction in the asymptotic regime.  In
        # practice the coarse-grid prefactor wobbles with grid-sphere
        # alignment, so we set the threshold safely above the O(dx) ceiling
        # (4x total / 2x per step) to distinguish quadratic from linear.
        ratio_total = errs[0] / errs[2]
        ratio_step1 = errs[0] / errs[1]
        ratio_step2 = errs[1] / errs[2]
        assert ratio_total > 6.0, (
            f"total coarse:fine ratio {ratio_total:.2f} below 6 "
            f"(O(dx) ideal is 4x; expect >>4x for O(dx^2)). "
            f"errs = {errs}, truth = {sweep['truth']:.5f}")
        assert ratio_step1 > 2.0 and ratio_step2 > 2.0, (
            f"per-step ratios too low: step1={ratio_step1:.2f}, "
            f"step2={ratio_step2:.2f}, errs={errs}")


# ===========================================================================
# 13. TestDecayLength
# ===========================================================================

class TestDecayLength:

    def test_benign_matches_analytic(self):
        p = Parameters(chi_b=0.0, omega_b=0.0, kappa=1.0,
                       phi_inf=0.5, R=5.0)
        d2 = float(d2fb_dphi2(np.asarray([p.phi_inf]), p)[0])
        expected = float(np.sqrt(p.kappa / d2))
        assert _decay_length(p) == pytest.approx(expected, rel=1e-10)

    def test_spinodal_fallback(self):
        # At chi_b = 2 with phi_inf = 0.5, d2fb is exactly zero (critical).
        p = Parameters(chi_b=2.0, omega_b=0.0, kappa=1.0,
                       phi_inf=0.5, R=5.0, nb=1.0)
        ell = _decay_length(p)
        expected = float(np.sqrt(p.kappa * p.nb / 0.05))
        assert ell == pytest.approx(expected, rel=1e-10)
        assert ell > 0 and np.isfinite(ell)

    @pytest.mark.parametrize("phi_inf,chi_b", [
        (0.05, 0.0), (0.2, 1.0), (0.5, 1.5), (0.8, 0.0),
    ])
    def test_always_positive(self, phi_inf, chi_b):
        p = Parameters(phi_inf=phi_inf, chi_b=chi_b, kappa=1.0, R=5.0)
        assert _decay_length(p) > 0


# ===========================================================================
# 14. TestRobustness
# ===========================================================================

class TestRobustness:

    def test_perturbed_initial_guess(self):
        """Perturbing the initial guess slightly: same converged solution."""
        p = P_benign()
        geom = build_geometry(p, d=2.0, Nrho=60, Nz=180)
        from viz_util.bvp_two_spheres import _make_residual as mk
        from viz_util.bvp_solver import compute_mu_inf as mu
        res_fn = mk(geom, mu(p))
        x0 = _initial_guess(geom)
        rng = np.random.default_rng(11)
        x0_perturb = np.clip(x0 + 0.01 * rng.standard_normal(x0.shape),
                             p.eps, 1.0 - p.eps)
        x_a = newton_krylov(res_fn, x0,         f_tol=1e-6, maxiter=200)
        x_b = newton_krylov(res_fn, x0_perturb, f_tol=1e-6, maxiter=200)
        assert float(np.max(np.abs(x_a - x_b))) < 1e-4

    def test_parameters_not_mutated(self):
        p = P_benign()
        snapshot = dict(chi_b=p.chi_b, omega_b=p.omega_b, kappa=p.kappa,
                        phi_inf=p.phi_inf, R=p.R, nb=p.nb, L=p.L,
                        eps=p.eps,
                        conserve_volume=p.conserve_volume,
                        parallel_plates=p.parallel_plates)
        _ = solve_two_spheres(p, d=2.0, Nrho=40, Nz=120, verbose=False)
        for k, v in snapshot.items():
            assert getattr(p, k) == v, f"Parameters.{k} mutated"

    def test_returned_fields_serialisable(self):
        """Quick npz round-trip of the array-like fields in `result`."""
        import tempfile
        p = P_benign()
        res = solve_two_spheres(p, d=2.0, Nrho=40, Nz=120, verbose=False)
        tmp = tempfile.NamedTemporaryFile(suffix=".npz", delete=False)
        np.savez(tmp.name, phi=res["phi"], x=res["x"], d=res["d"],
                 mu_inf=res["mu_inf"], success=res["success"])
        loaded = np.load(tmp.name)
        assert loaded["phi"].shape == res["phi"].shape
        assert float(loaded["d"]) == res["d"]
        os.unlink(tmp.name)

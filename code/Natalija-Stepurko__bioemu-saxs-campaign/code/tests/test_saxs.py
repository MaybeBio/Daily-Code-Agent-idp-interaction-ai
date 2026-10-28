import numpy as np

from bsc import saxs


def sphere_curve(q, R):
    x = q * R
    return (3 * (np.sin(x) - x * np.cos(x)) / x**3) ** 2


def test_scale_and_background_recovers_factors():
    q = np.linspace(0.01, 0.3, 200)
    model = sphere_curve(q, 20.0)
    exp = 3.5 * model + 0.02
    sigma = np.full_like(q, 0.01)
    c, b = saxs.scale_and_background(model, exp, sigma)
    assert abs(c - 3.5) < 1e-6 and abs(b - 0.02) < 1e-6
    assert saxs.chi2(model, exp, sigma) < 1e-10


def test_chi2_near_one_for_noise_at_sigma():
    rng = np.random.default_rng(0)
    q = np.linspace(0.01, 0.3, 400)
    model = sphere_curve(q, 20.0)
    sigma = 0.02 * model + 1e-4
    exp = model + rng.normal(0, sigma)
    c2 = saxs.chi2(model, exp, sigma)
    assert 0.8 < c2 < 1.25


def test_guinier_recovers_rg_of_sphere():
    R = 25.0
    q = np.linspace(0.005, 0.25, 500)
    I = sphere_curve(q, R)
    g = saxs.guinier(q, I, 0.001 * I + 1e-6)
    assert g["valid"]
    assert abs(g["rg"] - R * np.sqrt(3 / 5)) / (R * np.sqrt(3 / 5)) < 0.03
    assert abs(g["upturn"]) < 0.02


def test_guinier_flags_aggregation_upturn():
    R = 25.0
    q = np.linspace(0.005, 0.25, 500)
    I = sphere_curve(q, R) * (1 + 0.6 * np.exp(-((q / 0.012) ** 2)))
    g = saxs.guinier(q, I, 0.001 * I + 1e-6)
    assert g["upturn"] > 0.1


def test_reweighting_moves_towards_the_data():
    rng = np.random.default_rng(1)
    q = np.linspace(0.01, 0.3, 150)
    radii = np.linspace(15, 35, 40)
    curves = np.array([sphere_curve(q, r) for r in radii])
    truth = 0.8 * sphere_curve(q, 32.0) + 0.2 * sphere_curve(q, 16.0)
    sigma = 0.01 * truth + 1e-5
    exp = truth + rng.normal(0, sigma)
    raw = saxs.chi2(saxs.ensemble_curve(curves), exp, sigma)
    path = saxs.reweighting_curve(curves, exp, sigma, 10.0 ** np.linspace(-2, 6, 17))
    best = min(p["chi2"] for p in path)
    assert raw > 10 and best < 1.5
    phis = [p["phi"] for p in path]
    assert phis[0] > 0.9 and phis[-1] < phis[0]
    w, _ = saxs.reweight(curves, exp, sigma, 1e-3)
    assert abs(w.sum() - 1) < 1e-9
    assert radii[np.argmax(w)] > 28


def test_phi_at_chi2_and_chi2_at_phi():
    path = [{"theta": 10, "chi2": 9.0, "phi": 1.0}, {"theta": 1, "chi2": 2.0, "phi": 0.5},
            {"theta": 0.1, "chi2": 0.9, "phi": 0.1}]
    assert saxs.phi_at_chi2(path, 1.0) == 0.1
    assert saxs.phi_at_chi2(path, 2.5) == 0.5
    assert np.isnan(saxs.phi_at_chi2(path, 0.5))
    assert saxs.chi2_at_phi(path, 0.4) == 2.0


def test_nrmsd_log_on_a_synthetic_case():
    """Zero for the data themselves; a 10% error in sphere radius registers; a constant rescaling
    of sigma leaves the score unchanged while chi2 moves by the square of the factor."""
    q = np.linspace(0.01, 0.3, 150)
    truth = sphere_curve(q, 25.0)
    sigma = 0.02 * truth + 1e-5
    assert saxs.nrmsd_log(truth, truth, sigma) < 1e-12
    wrong = sphere_curve(q, 27.5)
    v = saxs.nrmsd_log(wrong, truth, sigma)
    assert 0.005 < v < 0.5
    assert abs(saxs.nrmsd_log(wrong, truth, 3 * sigma) - v) < 0.1 * v
    assert saxs.chi2(wrong, truth, 3 * sigma) < saxs.chi2(wrong, truth, sigma) / 5
    # reweighting_curve can hand back the weights of every point on the path
    curves = np.array([sphere_curve(q, r) for r in (24.0, 25.0, 26.0)])
    path, ws = saxs.reweighting_curve(curves, truth, sigma, np.array([100.0, 1.0]), return_weights=True)
    assert len(path) == len(ws) == 2 and all(abs(w.sum() - 1) < 1e-9 for w in ws)


def _background_case(seed: int):
    """Data on a realistic intensity scale with non-uniform errors and a non-zero buffer offset."""
    rng = np.random.default_rng(seed)
    q = np.linspace(0.01, 0.3, 120)
    radii = np.linspace(15, 35, 30)
    curves = 4e4 * np.array([sphere_curve(q, r) for r in radii])
    truth = 0.7 * curves[25] + 0.3 * curves[3] + 150.0          # needs a background offset
    sigma = 0.02 * truth + 5.0
    exp = truth + rng.normal(0, sigma)
    return curves, exp, sigma


def test_reweighting_gradient_matches_finite_differences():
    curves, exp, sigma = _background_case(5)
    f = saxs.reweight_objective(curves, exp, sigma, theta=0.5)
    x = np.random.default_rng(6).normal(0, 0.5, len(curves))
    _, g = f(x)
    h = 1e-5
    num = np.array([(f(x + h * e)[0] - f(x - h * e)[0]) / (2 * h) for e in np.eye(len(x))])
    assert np.max(np.abs(num - g)) / np.max(np.abs(num)) < 1e-4


def test_reweighting_with_a_background_offset_improves_on_uniform_weights():
    curves, exp, sigma = _background_case(7)
    uniform = saxs.chi2(saxs.ensemble_curve(curves), exp, sigma)
    w, c2 = saxs.reweight(curves, exp, sigma, theta=0.1)
    assert c2 < 0.5 * uniform and saxs.kish_fraction(w) < 0.9

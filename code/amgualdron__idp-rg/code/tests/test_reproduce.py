"""Reproduction checks that need no MD frames: descriptors, the released network and closed form against the stored
tables and results, and the distribution's numerics.   python -m pytest tests/   (or: python tests/test_reproduce.py)"""
import hashlib
import json
import pathlib
import sys

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import idprg  # noqa: E402
from idprg import closed_form as cf, data, descriptors, distribution as D, network  # noqa: E402

RESULTS = ROOT / "results"


def test_network_file_is_the_released_one():
    spec = json.loads((ROOT / "model/network.json").read_text())
    assert hashlib.sha256((ROOT / "model/network.eqx").read_bytes()).hexdigest() == spec["sha256"]["network.eqx"]


def test_descriptors_reproduce_the_protein_table():
    P = data.proteins().sample(300, random_state=0)
    X = descriptors.table(P.sequence)
    for c in descriptors.COLUMNS + ("log_seq_length", "N"):
        assert np.abs(X[c].to_numpy() - P[c].to_numpy()).max() < 1e-12, c


def test_nonstandard_residues_are_rejected():
    try:
        descriptors.descriptors("MKXV")
    except ValueError:
        return
    raise AssertionError("accepted X")


def _check(kind, params, name):
    r = pd.read_csv(RESULTS / f"{kind}_{name}.csv")
    lnA, lnSm, u = params
    assert np.abs(lnA - r.lnA.to_numpy()).max() < 1e-10
    assert np.abs(lnSm - r.lnSm.to_numpy()).max() < 1e-10
    assert np.abs(D.nu_of(u) - r.nu.to_numpy()).max() < 1e-10


def test_validation_parameters():
    P = data.proteins()
    V = P.iloc[data.split(P, "standard")[1]]
    net, spec = network.load()
    _check("network", network.parameters(net, spec, V), "val")
    _check("closed_form", cf.parameters(cf.load(), V), "val")


def test_benchmark_from_sequences():
    B = data.benchmark()
    for kind, closed in (("network", False), ("closed_form", True)):
        p = idprg.predict(B.sequence.tolist(), closed_form=closed)
        r = pd.read_csv(RESULTS / f"{kind}_benchmark.csv")
        assert np.abs(p.mean_rg.to_numpy() - r["mean"].to_numpy()).max() < 1e-9, kind
        assert np.abs(p.nu.to_numpy() - r.nu.to_numpy()).max() < 1e-10, kind


def test_distribution_numerics():
    lnA, lnSm, nu = 1.1, 1.3, 0.6
    rg = np.exp(lnSm + np.linspace(-1.5, 1.2, 200001))
    p = D.pdf(rg, lnA, lnSm, nu)
    assert abs(np.trapezoid(p, rg) - 1) < 1e-6
    m, s, k, e = D.moments(lnA, lnSm, nu)
    assert abs(np.trapezoid(p * rg, rg) - m) < 1e-6
    assert abs(np.sqrt(np.trapezoid(p * (rg - m) ** 2, rg)) - s) < 1e-6
    assert abs(D.cdf(D.quantiles(0.3, lnA, lnSm, nu), lnA, lnSm, nu) - 0.3) < 1e-6
    assert abs(rg[np.argmax(p)] - np.exp(lnSm)) < 1e-3          # S_m is the mode


def test_paper_numbers_are_current():
    """paper/numbers.tex and the generated tables of MODEL.md and README.md match results/ (paper/make_numbers.py)."""
    sys.path.insert(0, str(ROOT / "paper"))
    import make_numbers
    new, _ = make_numbers.outputs()
    stale = [p.name for p, t in new.items() if p.read_text() != t]
    assert not stale, f"out of date, run python paper/make_numbers.py: {stale}"


def test_closed_form_tex_coefficients():
    """Every coefficient typed in paper/closed_form.tex is that of model/closed_form.json, to 6 significant figures."""
    coef, text = cf.load(), (ROOT / "paper/closed_form.tex").read_text()
    values = [coef["ce_threshold"]["c0"]] + list(coef["ln_A"].values())
    for r in ("ln_Sm", "u"):
        values += [coef[r]["intercept"]] + list(coef[r]["coef"].values())
    for v in values:
        assert f"{v:.6g}" in text, v


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("passed", name)

"""
capped_model.py — the parameters of the capped hydration forcing, the main calibrated model from 2026-10-24
(saturating_forcing.py): charged-group hydrophilic term (D, E, K, R) and TEA term, both capped at c_sat = 2 full contacts
per residue; TEA strength at its prior; hydrophilic strength fitted on the four Wuttke constructs only (curves to 320 K,
the routine that fitted S1, heldout_exploratory.fit_free), with its standard error from the Jacobian of that fit
(SE = sqrt((J^T J)^-1 chi2 / dof)). Mpipi-GG has no Wuttke runs: its terms take the mean relative size SD(term)/SD(chi)
of the two calibrated models (total_response_evidence.load), now with the capped terms.
    BENDER_REANALYSIS_DIR=... python capped_model.py      (refit and write the parameter table)
Output: heldout_temperature/capped_forcing/capped_model_params.csv
"""
import sys

sys.dont_write_bytecode = True

import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.environ.get('BENDER_REANALYSIS_DIR', HERE))
sys.path.insert(0, HERE)

BM = (__import__('os').environ['CS_DATA'] + '/BENDER-MT')
PARAMS = f'{BM}/heldout_temperature/capped_forcing/capped_model_params.csv'
CSAT = 2.0
HYDSET = 'DEKR'
GAMMA_TEA = {'C': 2.0, 'R': 0.4}            # priors, as the linear calibration
KEY = {'C': 'calvados', 'R': 'recharged'}


def fit(model):
    from scipy.optimize import least_squares
    import saturating_forcing as SF
    H = SF.HX.H
    H.OUT = SF.HA.OUT
    units = SF.wuttke_units(model, 'S1', CSAT)

    def resid(p):
        gam = SF.gam_of(model, p[0])
        r = [H.unit_resid(model, u, kind, d, gam, None) for u, (kind, d) in units.items()]
        return np.concatenate(r + [np.array([p[0] / 10.0])])

    sol = least_squares(resid, x0=[0.3 if model == 'recharged' else 1.0], method='trf', bounds=(0, np.inf))
    r = resid(sol.x)
    dof = max(1, len(r) - 1)
    cov = np.linalg.inv(sol.jac.T @ sol.jac) * (r @ r) / dof
    return float(sol.x[0]), float(np.sqrt(cov[0, 0])), float(np.sqrt(np.mean(r[:-1] ** 2)))


def params():
    """{'C': (gamma_h, se), 'R': (gamma_h, se)} and the TEA priors, from the stored table."""
    t = pd.read_csv(PARAMS).set_index('model')
    return {k: (float(t.loc[k, 'gamma_h']), float(t.loc[k, 'gamma_h_se'])) for k in ('C', 'R')}, dict(GAMMA_TEA)


def main():
    rows = []
    for k, model in KEY.items():
        g, se, rms = fit(model)
        ref = pd.read_csv(f'{BM}/heldout_temperature/results/exploratory/saturating/saturating_{model}.csv')
        ref = ref[(ref.variant == 'S1') & (ref.c_sat == CSAT)].gamma_h.iloc[0]
        if abs(g - ref) > 1e-3 * ref:
            raise ValueError(f'{model}: refit {g} differs from the scan value {ref}')
        rows.append({'model': k, 'name': model, 'c_sat': CSAT, 'hydrophilic_set': HYDSET, 'gamma_tea': GAMMA_TEA[k],
                     'gamma_h': g, 'gamma_h_se': se, 'wuttke_rms': rms})
        print(f'{model}: gamma_h {g:.4f} +- {se:.4f} (Wuttke rms {rms:.3f}; scan value {ref:.4f})')
    os.makedirs(os.path.dirname(PARAMS), exist_ok=True)
    pd.DataFrame(rows).to_csv(PARAMS, index=False)


if __name__ == '__main__':
    main()

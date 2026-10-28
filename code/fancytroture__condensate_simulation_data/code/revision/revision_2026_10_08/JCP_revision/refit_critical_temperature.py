"""Reproduce the reported three-temperature estimates from the bundled CSV."""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, minimize_scalar


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path(__file__).resolve().parent.parent /
                        'critical_temperature_slab_simulations/density_3to7us_20261012/si_figure/critical_temperature_input_densities.csv')
    args = parser.parse_args()
    with args.input.open() as stream:
        rows = list(csv.DictReader(stream))
    print('Unweighted fit: rho_dense - rho_dilute = A * (Tc - T)^0.325')
    print('Reference: Dignon et al., PLoS Comput Biol 14, e1005941 (2018).')
    print('https://doi.org/10.1371/journal.pcbi.1005941')
    for system in ('FUS', 'A1', 'TDP43'):
        selected = [r for r in rows if r['system'] == system]
        assert len(selected) == 3
        t = np.array([float(r['temperature_K']) for r in selected])
        gap = np.array([float(r['dense_g_ml']) - float(r['dilute_g_ml']) for r in selected])

        def objective(tc):
            basis = (tc - t) ** .325
            a = np.dot(basis, gap) / np.dot(basis, basis)
            return np.sum((a * basis - gap) ** 2)

        result = minimize_scalar(objective, bounds=(max(t) + .001, 650),
                                 method='bounded', options={'xatol': 1e-8})
        assert result.success
        tc = result.x
        basis = (tc - t) ** .325
        a = np.dot(basis, gap) / np.dot(basis, basis)
        check = least_squares(lambda p: p[0] * (p[1] - t) ** .325 - gap,
                              [a, tc], bounds=([0, max(t) + .001], [10, 650]),
                              xtol=1e-12, ftol=1e-12, gtol=1e-12)
        np.testing.assert_allclose(check.x, [a, tc], rtol=1e-6, atol=1e-5)
        print(f'{system}: Tc={tc:.6f} K; A={a:.9f}; RMSE={np.sqrt(result.fun/3):.9f} g/mL')


if __name__ == '__main__':
    main()

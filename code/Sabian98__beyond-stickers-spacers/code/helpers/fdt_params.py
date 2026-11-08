"""fdt_params.py: definitions copied verbatim from the unpublished TODO_T2_fdt_computation.py (BENDER-MT) by tools/extract_helpers.py.
Do not edit by hand."""
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse

from typing import Tuple

# Units: energies in kJ/mol, distances in nm, temperature in K
R_GAS = 8.314462618e-3    # kJ mol^-1 K^-1; beta = 1/(R_GAS T) is the sampled Boltzmann factor

R_CALVADOS = 8.3145e-3    # kT = 8.3145e-3 * T inside CALVADOS genParamsDH (Bjerrum length only)

EPS_AH = 0.2 * 4.184      # kJ/mol (eps_lj = 0.2 kcal/mol)

RC_AH = 2.0               # nm, Ashbaugh-Hatch cutoff (cutoff_lj)

RC_YU = 4.0               # nm, Debye-Hückel cutoff (cutoff_yu)

IONIC_STRENGTH = 0.15     # mol/L

PH = 7.5

HIS_PKA = 6.0

def eps_w(T: float) -> float:
    """
    Water relative permittivity [Derivation §10.1, Eq. 10.1], as in CALVADOS genParamsDH.
    """
    return (5321.0 / T
            + 233.76
            - 0.9297 * T
            + 0.001417 * T**2
            - 0.8292e-6 * T**3)

def deps_w_dT(T: float) -> float:
    """
    d(eps_w)/dT in K^-1 [Derivation §10.1, Eq. 10.2-10.3].
    """
    return (-5321.0 / T**2
            - 0.9297
            + 0.002834 * T
            - 2.4876e-6 * T**2)

def bjerrum_length(T: float) -> float:
    """
    Bjerrum length l_B(T) in nm [Derivation §10.2, Eq. 10.7]:

        l_B = C / (eps_w(T) * kT),   C = e^2 N_A / (4 pi eps_0) = 138.9 kJ nm/mol

    with kT = R_CALVADOS * T exactly as in CALVADOS genParamsDH.
    """
    C = (1.6021766**2 * 6.02214076 * 1000) / (4 * np.pi * 8.854188)
    return C / (eps_w(T) * R_CALVADOS * T)

def screening_kappa(T: float, ionic: float = IONIC_STRENGTH) -> float:
    """
    Inverse Debye screening length kappa(T) in nm^-1 [Derivation §10.4, Eq. 10.21]:
    kappa = sqrt(8 pi l_B I N_A / 10).
    """
    D = 8 * np.pi * 6.02214076 / 10
    return np.sqrt(D * bjerrum_length(T) * ionic)

def dkappa_dT(T: float, ionic: float = IONIC_STRENGTH) -> float:
    """
    d(kappa)/dT in nm^-1 K^-1 at fixed ionic strength [Derivation §10.4, Eq. 10.25]:
    kappa'/kappa = (1/2) l_B'/l_B = -(1/2)(1/T + eps_w'/eps_w)
    """
    return -(screening_kappa(T, ionic) / 2.0) * (1.0 / T + deps_w_dT(T) / eps_w(T))

def get_dielectric_params(T: float, ionic: float = IONIC_STRENGTH) -> dict:
    """
    All temperature-dependent electrostatic parameters at T. Call once per temperature.

    eps_yu = l_B kT = C / eps_w is the Yukawa energy prefactor (kJ nm/mol); its derivative is
    eps_yu' = -eps_yu * eps_w'/eps_w [Derivation Eq. 10.18].
    """
    ew, dew = eps_w(T), deps_w_dT(T)
    lB = bjerrum_length(T)
    eps_yu = lB * R_CALVADOS * T
    return {
        'eps_w': ew,
        'deps_w_dT': dew,
        'l_B': lB,
        'kappa': screening_kappa(T, ionic),
        'dkappa_dT': dkappa_dT(T, ionic),
        'eps_yu': eps_yu,
        'deps_yu_dT': -eps_yu * dew / ew,
        'beta': 1.0 / (R_GAS * T),
        'ionic': ionic,
        'T': float(T),
    }

def his_charge(ph: float = PH) -> float:
    """CALVADOS flexhis histidine charge, 1/(1 + 10^(pH - 6)); 0.031 at pH 7.5."""
    return 1.0 / (1.0 + 10.0 ** (ph - HIS_PKA))

def calvados_residue_params(
    sequence: str,
    residues_df,
    ph: float = PH,
    charge_termini: str = 'both',
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Per-residue CALVADOS parameters for one chain.

    Args:
        sequence:       one-letter sequence (use the topology's sequence, which is what was simulated)
        residues_df:    residues_CALVADOS2.csv indexed by one-letter code
                        (columns lambdas, sigmas, q, bondlength)
        ph:             simulation pH (prepare.py: 7.5)
        charge_termini: 'both' | 'N' | 'C' | 'none' (prepare.py: 'both')
    Returns:
        charges (N,), lambdas (N,), sigmas (N,) in nm, bond_r0 (N-1,) in nm
    """
    unknown = sorted(set(sequence) - set(residues_df.index))
    if unknown:
        raise ValueError(f'residues not in the CALVADOS table: {unknown}')
    r = residues_df.loc[list(sequence)]
    charges = r['q'].to_numpy(dtype=float).copy()
    charges[np.array(list(sequence)) == 'H'] = his_charge(ph)
    if charge_termini in ('N', 'both'):
        charges[0] += 1.0
    if charge_termini in ('C', 'both'):
        charges[-1] -= 1.0
    bl = r['bondlength'].to_numpy(dtype=float)
    return (charges, r['lambdas'].to_numpy(dtype=float), r['sigmas'].to_numpy(dtype=float),
            0.5 * (bl[:-1] + bl[1:]))

def pair_setup(charges: np.ndarray, lambdas: np.ndarray, sigmas: np.ndarray):
    """
    Nonbonded pair list (i < j, j - i >= 2; CALVADOS excludes only bonded pairs) and mixed parameters.
    Returns i, j, sigma_ij, lambda_ij, AH shift, q_i q_j.
    """
    i, j = np.triu_indices(len(charges), k=2)
    s = 0.5 * (sigmas[i] + sigmas[j])
    l = 0.5 * (lambdas[i] + lambdas[j])
    shift = (s / RC_AH) ** 12 - (s / RC_AH) ** 6
    return i, j, s, l, shift, charges[i] * charges[j]

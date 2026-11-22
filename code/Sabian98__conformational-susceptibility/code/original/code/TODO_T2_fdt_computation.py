"""
TODO-T2: Per-frame G_T computation for FDT validation
=======================================================
BENDER-MT Paper 1 — Conformational Susceptibility

WHAT WE ARE COMPUTING AND WHY
------------------------------
For every retained frame x_k of a protein simulated at temperature T we need one scalar,
the generalized thermal force (the derivative of the reduced energy at fixed coordinates):

    u(x, T)   = beta(T) * H(x, T),        beta = 1 / (R T)
    G_T(x_k)  = du/dT |_x                                   [Derivation Eq. 9.5, 12.3]

The susceptibility of any conformational observable A (Rg, g_eff, transitivity, ...) is then

    chi_T(A) = d<A>/dT = -Cov(A, G_T)                       [Derivation Eq. 9.16]

CALVADOS-2 ENERGY
-----------------
    H(x, T) = V_bond(x) + V_AH(x) + V_DH(x, T)

Only the Debye-Hückel term depends on T (through eps_w(T) and kappa(T)). Hence

    G_T = -(V_bond + V_AH) / (R T^2) + du_DH/dT,     u_DH = beta V_DH = sum q_i q_j l_B(T) B(r_ij, T)

du_DH/dT is the derivative of the REDUCED DH energy. It already contains the Boltzmann piece
-V_DH/(R T^2) (Derivation Eq. 10.41 and 10.45 are the same quantity), so V_DH must not appear
in the first term as well.

TWO EXACT DECOMPOSITIONS (each sums to chi_total)
-------------------------------------------------
(a) reweighting vs force-field change [Derivation §9, §12.4]
        chi_Boltzmann  = +Cov(A, H) / (R T^2)                 H = V_bond + V_AH + V_DH
        chi_dielectric = -beta * Cov(A, dV_DH/dT)             partial derivative at fixed x
    with dV_DH/dT = R T * du_DH/dT + V_DH / T
        f_dielectric   = |chi_dielectric| / (|chi_Boltzmann| + |chi_dielectric|)
(b) by energy term [Derivation §11]
        chi_bond = Cov(A, V_bond) / (R T^2)
        chi_AH   = Cov(A, V_AH)   / (R T^2)
        chi_DH   = -Cov(A, du_DH/dT)

FORCE FIELD (as built by the calvados package for prepare.py's Components; checked against
KULL-Centre/CALVADOS calvados/interactions.py, sequence.py, components.py, data/default_*.yaml)
----------------------------------------------------------------------------------------------
  Ashbaugh-Hatch:  eps = 0.2 kcal/mol = 0.8368 kJ/mol, rc = 2.0 nm, shifted at rc,
                   sigma_ij = (sigma_i + sigma_j)/2, lambda_ij = (lambda_i + lambda_j)/2 (arithmetic)
  Debye-Hückel:    eps_yu = l_B kT, kappa = sqrt(8 pi l_B I N_A/10), rc = 4.0 nm, shifted at rc
  Bonds:           harmonic, k = 8033 kJ/mol/nm^2, r0 = mean of the residues' bondlength (0.38 nm)
  Exclusions:      bonded (i, i+1) pairs only
  Charges:         residues table q; His = 1/(1 + 10^(pH-6)) (flexhis); charge_termini='both'
                   adds +1 to residue 1 and -1 to residue N; pH 7.5, I = 0.15 M (prepare.py)

VALIDATION
----------
  local slope    d<A>/dT from the five ensemble means (np.gradient on the non-uniform grid),
                 compared with chi_T(A) at each temperature
  integration    integral of chi_T dT over 278-320 K (trapezoid) = <A>_320 - <A>_278
                 [Derivation §13], judged by a z-score from the block standard errors

run_sanity_checks() tests the pair energies against the literal CALVADOS expressions, G_T and
dV_DH/dT against central finite differences in T at fixed coordinates, and both decompositions
and the pairwise map against chi_total.
"""

import numpy as np
from typing import Tuple

# ===========================================================================
# SECTION 1: CALVADOS-2 CONSTANTS AND DIELECTRIC PARAMETERS
# ===========================================================================

# Units: energies in kJ/mol, distances in nm, temperature in K
R_GAS = 8.314462618e-3    # kJ mol^-1 K^-1; beta = 1/(R_GAS T) is the sampled Boltzmann factor
R_CALVADOS = 8.3145e-3    # kT = 8.3145e-3 * T inside CALVADOS genParamsDH (Bjerrum length only)

EPS_AH = 0.2 * 4.184      # kJ/mol (eps_lj = 0.2 kcal/mol)
RC_AH = 2.0               # nm, Ashbaugh-Hatch cutoff (cutoff_lj)
RC_YU = 4.0               # nm, Debye-Hückel cutoff (cutoff_yu)
KB_BOND = 8033.0          # kJ mol^-1 nm^-2
IONIC_STRENGTH = 0.15     # mol/L
PH = 7.5
HIS_PKA = 6.0

MAX_ELEMS = 4e6           # frames x pairs per chunk (bounds memory for long chains)
N_BLOCKS_CHI = 10         # contiguous blocks for the jackknife standard errors


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


def dbjerrum_dT(T: float) -> float:
    """
    d(l_B)/dT in nm/K [Derivation §10.2, Eq. 10.12]:  l_B'/l_B = -1/T - eps_w'/eps_w
    """
    return bjerrum_length(T) * (-1.0 / T - deps_w_dT(T) / eps_w(T))


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


# ===========================================================================
# SECTION 2: PAIR ENERGY FUNCTIONS
# ===========================================================================

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


def ah_pair_energy(r: np.ndarray, s: np.ndarray, l: np.ndarray, shift: np.ndarray) -> np.ndarray:
    """
    Ashbaugh-Hatch pair energy (kJ/mol), temperature independent. CALVADOS expression:

        eps * select(step(r - 2^(1/6) s), 4 l ((s/r)^12 - (s/r)^6 - shift),
                                          4 ((s/r)^12 - (s/r)^6 - l shift) + (1 - l)),   r < rc
    """
    x6 = (s / r) ** 6
    lj = x6 * x6 - x6
    v = np.where(r < 2.0 ** (1.0 / 6.0) * s,
                 4.0 * EPS_AH * (lj - l * shift) + EPS_AH * (1.0 - l),
                 4.0 * EPS_AH * l * (lj - shift))
    return np.where(r < RC_AH, v, 0.0)


def dh_pair_energy(r: np.ndarray, qq: np.ndarray, params: dict) -> Tuple[np.ndarray, np.ndarray]:
    """
    Debye-Hückel pair energy V_ij = q_i q_j eps_yu B(r, T) (kJ/mol, zero beyond rc), with the shifted
    Yukawa function B = exp(-kappa r)/r - exp(-kappa rc)/rc [Derivation Eq. 10.30], and its partial
    temperature derivative at fixed r [Derivation Eq. 10.39]:

        dV_ij/dT = q_i q_j [eps_yu' B + eps_yu kappa' (-exp(-kappa r) + exp(-kappa rc))]   (kJ/mol/K)
    """
    kappa = params['kappa']
    e_r = np.exp(-kappa * r)
    e_c = np.exp(-kappa * RC_YU)
    B = e_r / r - e_c / RC_YU
    inside = r < RC_YU
    v = np.where(inside, qq * params['eps_yu'] * B, 0.0)
    dv = np.where(inside, qq * (params['deps_yu_dT'] * B
                                + params['eps_yu'] * params['dkappa_dT'] * (e_c - e_r)), 0.0)
    return v, dv


def bond_energy(xyz_frames: np.ndarray, bond_r0) -> np.ndarray:
    """Harmonic bond energy per frame, sum 0.5 k (b - r0)^2 (kJ/mol)."""
    b = np.linalg.norm(np.diff(xyz_frames, axis=1), axis=2)
    return (0.5 * KB_BOND * (b - bond_r0) ** 2).sum(axis=1)


# ===========================================================================
# SECTION 3: PER-FRAME ENERGIES AND G_T
# ===========================================================================

def per_frame_energies(
    xyz_frames: np.ndarray,
    charges: np.ndarray,
    lambdas: np.ndarray,
    sigmas: np.ndarray,
    bond_r0,
    params: dict,
    max_elems: float = MAX_ELEMS,
) -> dict:
    """
    Per-frame energies and the generalized thermal force for one protein at params['T'].

    Args:
        xyz_frames: (n_frames, N, 3) in nm
        charges, lambdas, sigmas, bond_r0: from calvados_residue_params()
        params:     from get_dielectric_params(T)
    Returns dict of (n_frames,) arrays:
        V_AH, V_DH, V_bond, H      kJ/mol
        dVdT_DH                    kJ/mol/K, partial dV_DH/dT at fixed coordinates
        dudT_DH                    K^-1, d(beta V_DH)/dT
        G_T                        K^-1, d(beta H)/dT = -(V_AH + V_bond)/(R T^2) + dudT_DH
    """
    xyz_frames = np.asarray(xyz_frames, dtype=float)
    n_frames = xyz_frames.shape[0]
    T = params['T']
    RT = R_GAS * T
    i, j, s, l, shift, qq = pair_setup(charges, lambdas, sigmas)
    ch = qq != 0
    V_AH = np.zeros(n_frames)
    V_DH = np.zeros(n_frames)
    dVdT = np.zeros(n_frames)
    step = max(1, int(max_elems // max(1, len(i))))
    for a in range(0, n_frames, step):
        x = xyz_frames[a:a + step]
        r = np.linalg.norm(x[:, i] - x[:, j], axis=2)
        V_AH[a:a + step] = ah_pair_energy(r, s, l, shift).sum(axis=1)
        if ch.any():
            v, dv = dh_pair_energy(r[:, ch], qq[ch], params)
            V_DH[a:a + step] = v.sum(axis=1)
            dVdT[a:a + step] = dv.sum(axis=1)
    V_bond = bond_energy(xyz_frames, bond_r0)
    dudT = -V_DH / (RT * T) + dVdT / RT
    return {
        'V_AH': V_AH,
        'V_DH': V_DH,
        'V_bond': V_bond,
        'H': V_AH + V_DH + V_bond,
        'dVdT_DH': dVdT,
        'dudT_DH': dudT,
        'G_T': -(V_AH + V_bond) / (RT * T) + dudT,
    }


# ===========================================================================
# SECTION 4: SUSCEPTIBILITY ESTIMATES
# ===========================================================================

def _cov(x: np.ndarray, y: np.ndarray) -> float:
    """Unbiased sample covariance (n - 1 denominator) [Derivation §12.2]."""
    return float(((x - x.mean()) * (y - y.mean())).sum() / (len(x) - 1))


def susceptibility(A_frames: np.ndarray, energies: dict, T: float, n_blocks: int = N_BLOCKS_CHI) -> dict:
    """
    chi_T(A) and its two decompositions from paired per-frame arrays.

    A_frames and every array in `energies` must come from the same frames, in time order.
    Frames where A is undefined (NaN, e.g. degree assortativity of a regular graph) are dropped
    from both. Standard errors come from a delete-one-block jackknife over n_blocks contiguous
    blocks [Derivation §12.5].

    Returns dict with chi_total, chi_total_se, chi_Boltzmann, chi_dielectric, f_dielectric,
    chi_AH, chi_bond, chi_DH (units of A per K), A_mean, A_mean_se, A_sd, n_frames.
    """
    A = np.asarray(A_frames, dtype=float)
    ok = np.isfinite(A)
    n = int(ok.sum())
    keys = ['chi_total', 'chi_total_se', 'chi_Boltzmann', 'chi_dielectric', 'f_dielectric',
            'chi_AH', 'chi_bond', 'chi_DH', 'A_mean', 'A_mean_se', 'A_sd']
    if n < 3:
        return {**{k: np.nan for k in keys}, 'n_frames': n}
    a = A[ok]
    e = {k: np.asarray(v, dtype=float)[ok] for k, v in energies.items()}
    RT = R_GAS * T

    def chis(idx):
        aa = a[idx]
        c = lambda k: _cov(aa, e[k][idx])  # noqa: E731
        return {
            'chi_total': -c('G_T'),
            'chi_Boltzmann': c('H') / (RT * T),
            'chi_dielectric': -c('dVdT_DH') / RT,
            'chi_AH': c('V_AH') / (RT * T),
            'chi_bond': c('V_bond') / (RT * T),
            'chi_DH': -c('dudT_DH'),
        }

    out = chis(np.arange(n))
    denom = abs(out['chi_Boltzmann']) + abs(out['chi_dielectric'])
    out['f_dielectric'] = abs(out['chi_dielectric']) / denom if denom > 0 else np.nan
    out['A_mean'] = float(a.mean())
    out['A_sd'] = float(a.std(ddof=1))
    if n >= 2 * n_blocks:
        blocks = np.array_split(np.arange(n), n_blocks)
        jk = np.array([chis(np.setdiff1d(np.arange(n), b))['chi_total'] for b in blocks])
        out['chi_total_se'] = float(np.sqrt((n_blocks - 1) / n_blocks * ((jk - jk.mean()) ** 2).sum()))
        bm = np.array([a[b].mean() for b in blocks])
        out['A_mean_se'] = float(bm.std(ddof=1) / np.sqrt(n_blocks))
    else:
        out['chi_total_se'] = np.nan
        out['A_mean_se'] = np.nan
    out['n_frames'] = n
    return out


# ===========================================================================
# SECTION 5: VALIDATION — TEMPERATURE SLOPES AND INTEGRATION
# ===========================================================================

TEMPS = np.array([278., 288., 300., 310., 320.])


def local_slopes(A_means: np.ndarray, temps: np.ndarray = TEMPS) -> np.ndarray:
    """
    d<A>/dT at each temperature from the ensemble means: second-order central differences on the
    non-uniform grid for interior points, one-sided at 278 and 320 K. Compare with chi_T(A; T).
    """
    return np.gradient(np.asarray(A_means, dtype=float), np.asarray(temps, dtype=float))


def linear_slope(A_means: np.ndarray, temps: np.ndarray = TEMPS) -> float:
    """Least-squares slope of <A> against T over all temperatures (the dataset's slope feature)."""
    T_c = temps - temps.mean()
    A_c = A_means - A_means.mean()
    return float((A_c * T_c).sum() / (T_c ** 2).sum())


def global_integration_check(
    chi_T_values: np.ndarray,
    chi_T_se: np.ndarray,
    A_means: np.ndarray,
    A_means_se: np.ndarray,
    temps: np.ndarray = TEMPS,
) -> dict:
    """
    integral of chi_T dT (trapezoid) vs <A>_last - <A>_first [Derivation §13, Eq. 13.3].

    Temperatures are independent simulations, so the standard errors add in quadrature.
    Returns integral, integral_se, observed_diff, observed_diff_se, z, relative_error.
    """
    t = np.asarray(temps, dtype=float)
    w = np.zeros(len(t))
    dt = np.diff(t)
    w[:-1] += dt / 2
    w[1:] += dt / 2
    integral = float((w * chi_T_values).sum())
    integral_se = float(np.sqrt((w ** 2 * np.asarray(chi_T_se) ** 2).sum()))
    observed = float(A_means[-1] - A_means[0])
    observed_se = float(np.sqrt(A_means_se[-1] ** 2 + A_means_se[0] ** 2))
    tot = np.hypot(integral_se, observed_se)
    return {
        'integral': integral,
        'integral_se': integral_se,
        'observed_diff': observed,
        'observed_diff_se': observed_se,
        'z': (integral - observed) / tot if tot > 0 else np.nan,
        'relative_error': abs(integral - observed) / abs(observed) if observed != 0 else np.nan,
    }


# ===========================================================================
# SECTION 6: PAIRWISE SUSCEPTIBILITY MAPS (TODO-T4)
# ===========================================================================

def pairwise_susceptibility_map(
    A_frames: np.ndarray,
    xyz_frames: np.ndarray,
    charges: np.ndarray,
    lambdas: np.ndarray,
    sigmas: np.ndarray,
    bond_r0,
    params: dict,
    max_elems: float = MAX_ELEMS,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Per-pair susceptibility chi_ij(A) = -Cov(A, d(beta V_ij)/dT) [Derivation Eq. 11.10]:

        AH pairs:  d(beta V_ij)/dT = -V_AH_ij / (R T^2)
        DH pairs:  d(beta V_ij)/dT = du_DH_ij/dT   (complete derivative, Eq. 10.45)

    plus the bonded term chi_0 = Cov(A, V_bond)/(R T^2). Exact sum rule [Derivation Eq. 11.11]:
    sum_{i<j} (chi_AH_ij + chi_DH_ij) + chi_0 = chi_total.

    Returns chi_map_AH (N, N), chi_map_DH (N, N), both symmetric, and chi_0.
    """
    A = np.asarray(A_frames, dtype=float)
    ok = np.isfinite(A)
    a = A[ok] - A[ok].mean()
    xyz = np.asarray(xyz_frames, dtype=float)[ok]
    n, N = len(a), xyz.shape[1]
    T = params['T']
    RT = R_GAS * T
    i, j, s, l, shift, qq = pair_setup(charges, lambdas, sigmas)
    chi_AH = np.zeros(len(i))
    chi_DH = np.zeros(len(i))
    step = max(1, int(max_elems // max(1, n)))
    for p in range(0, len(i), step):
        sl = slice(p, p + step)
        r = np.linalg.norm(xyz[:, i[sl]] - xyz[:, j[sl]], axis=2)
        chi_AH[sl] = (a @ ah_pair_energy(r, s[sl], l[sl], shift[sl])) / (n - 1) / (RT * T)
        v, dv = dh_pair_energy(r, qq[sl], params)
        chi_DH[sl] = -(a @ (-v / (RT * T) + dv / RT)) / (n - 1)
    chi_0 = float(a @ bond_energy(xyz, bond_r0)) / (n - 1) / (RT * T)
    maps = []
    for vals in (chi_AH, chi_DH):
        M = np.zeros((N, N))
        M[i, j] = vals
        M[j, i] = vals
        maps.append(M)
    return maps[0], maps[1], chi_0


# ===========================================================================
# SECTION 7: SANITY CHECKS
# ===========================================================================

def _calvados_literal_pair(r, s1, s2, l1, l2, q1, q2, T, ionic=IONIC_STRENGTH):
    """AH + DH pair energy transcribed literally from calvados/interactions.py (independent of Section 2)."""
    s = 0.5 * (s1 + s2)
    l = 0.5 * (l1 + l2)
    shift = (s / RC_AH) ** 12 - (s / RC_AH) ** 6
    if r - 2 ** (1 / 6) * s >= 0:
        ah = EPS_AH * 4 * l * ((s / r) ** 12 - (s / r) ** 6 - shift)
    else:
        ah = EPS_AH * (4 * ((s / r) ** 12 - (s / r) ** 6 - l * shift) + (1 - l))
    kT = 8.3145 * T * 1e-3
    fepsw = lambda t: 5321 / t + 233.76 - 0.9297 * t + 0.1417 * 1e-2 * t * t - 0.8292 * 1e-6 * t ** 3  # noqa: E731
    lB = 1.6021766 ** 2 / (4 * np.pi * 8.854188 * fepsw(T)) * 6.02214076 * 1000 / kT
    eps_yu = lB * kT
    k_yu = np.sqrt(8 * np.pi * lB * ionic * 6.02214076 / 10)
    dh = q1 * q2 * eps_yu * (np.exp(-k_yu * r) / r - np.exp(-k_yu * RC_YU) / RC_YU)
    return (ah if r < RC_AH else 0.0), (dh if r < RC_YU else 0.0)


def _synthetic_chain(n_res=40, n_frames=24, seed=1):
    """Random-walk chains with CALVADOS-like parameters, for derivative and sum-rule tests."""
    rng = np.random.default_rng(seed)
    steps = rng.normal(size=(n_frames, n_res - 1, 3))
    steps *= (0.38 * rng.uniform(0.95, 1.05, (n_frames, n_res - 1, 1))
              / np.linalg.norm(steps, axis=2, keepdims=True))
    xyz = np.concatenate([np.zeros((n_frames, 1, 3)), np.cumsum(steps, axis=1)], axis=1)
    q = rng.choice([-1.0, 0.0, 0.0, 1.0], size=n_res)
    q[0] += 1.0
    q[-1] -= 1.0
    lam = rng.uniform(0.0, 1.0, n_res)
    sig = rng.uniform(0.45, 0.68, n_res)
    return xyz, q, lam, sig, np.full(n_res - 1, 0.38)


def run_sanity_checks(verbose: bool = True):
    """
    Numerical checks; raise AssertionError on any failure. Run before the full pipeline.
    """
    say = print if verbose else (lambda *a, **k: None)
    say("Running sanity checks...")

    # 1-3: dielectric parameters at 300 K
    ew_300 = eps_w(300.0)
    assert 77 < ew_300 < 80, f"eps_w(300K) = {ew_300:.2f}, expected ~77.7"
    lB_300 = bjerrum_length(300.0)
    assert 0.6 < lB_300 < 0.8, f"l_B(300K) = {lB_300:.3f} nm, expected ~0.71 nm"
    kap_inv = 1.0 / screening_kappa(300.0, 0.15)
    assert 0.7 < kap_inv < 1.0, f"kappa^-1(300K) = {kap_inv:.3f} nm, expected ~0.79 nm"
    say(f"  eps_w(300K) = {ew_300:.3f}, l_B = {lB_300:.4f} nm, kappa^-1 = {kap_inv:.4f} nm ✓")

    # 4: signs. eps_w falls with T; l_B = C/(eps_w kT) rises because eps_w T falls.
    assert deps_w_dT(300.0) < 0, "deps_w/dT should be negative"
    assert dbjerrum_dT(300.0) > 0, "dl_B/dT should be positive in CALVADOS"
    say(f"  deps_w/dT(300K) = {deps_w_dT(300.0):.4f} K^-1, dl_B/dT = {dbjerrum_dT(300.0):.5f} nm/K ✓")

    # 5: analytic parameter derivatives vs central differences
    h = 1e-3
    for f, df in ((eps_w, deps_w_dT), (bjerrum_length, dbjerrum_dT), (screening_kappa, dkappa_dT)):
        fd = (f(300.0 + h) - f(300.0 - h)) / (2 * h)
        assert abs(df(300.0) - fd) < 1e-6 * max(1.0, abs(fd)), f"{df.__name__}: {df(300.0)} vs FD {fd}"
    say("  eps_w', l_B', kappa' match finite differences ✓")

    # 6: pair energies vs the literal CALVADOS expressions (both AH branches, inside/outside cutoffs)
    p300 = get_dielectric_params(300.0)
    s1, s2, l1, l2, q1, q2 = 0.60, 0.65, 0.97, 0.18, 1.0, -1.0
    s = 0.5 * (s1 + s2)
    for r in (0.55, 0.62, 0.70, 0.9, 1.9, 2.5, 3.9, 4.2):
        ah_ref, dh_ref = _calvados_literal_pair(r, s1, s2, l1, l2, q1, q2, 300.0)
        shift = (s / RC_AH) ** 12 - (s / RC_AH) ** 6
        ah = float(ah_pair_energy(np.array(r), np.array(s), np.array(0.5 * (l1 + l2)), np.array(shift)))
        dh = float(dh_pair_energy(np.array(r), np.array(q1 * q2), p300)[0])
        assert abs(ah - ah_ref) < 1e-9 and abs(dh - dh_ref) < 1e-9, f"pair energy mismatch at r={r}"
    say("  AH and DH pair energies equal the CALVADOS expressions ✓")

    # 7: G_T and dV_DH/dT vs central differences of beta H and V_DH at fixed coordinates
    xyz, q, lam, sig, b0 = _synthetic_chain()
    T0, h = 300.0, 1e-3
    E = per_frame_energies(xyz, q, lam, sig, b0, get_dielectric_params(T0))
    Ep = per_frame_energies(xyz, q, lam, sig, b0, get_dielectric_params(T0 + h))
    Em = per_frame_energies(xyz, q, lam, sig, b0, get_dielectric_params(T0 - h))
    G_fd = (Ep['H'] / (R_GAS * (T0 + h)) - Em['H'] / (R_GAS * (T0 - h))) / (2 * h)
    dV_fd = (Ep['V_DH'] - Em['V_DH']) / (2 * h)
    scale = np.abs(E['G_T']).max()
    assert np.abs(E['G_T'] - G_fd).max() < 1e-6 * scale, "G_T differs from d(beta H)/dT"
    assert np.abs(E['dVdT_DH'] - dV_fd).max() < 1e-6 * max(1e-12, np.abs(dV_fd).max()), \
        "dVdT_DH differs from dV_DH/dT"
    assert np.allclose(Ep['V_AH'], Em['V_AH']) and np.allclose(Ep['V_bond'], Em['V_bond'])
    say(f"  G_T = d(beta H)/dT and dV_DH/dT match finite differences "
        f"(max rel. diff {np.abs(E['G_T'] - G_fd).max() / scale:.1e}) ✓")

    # 8: decompositions and the pairwise sum rule
    A = np.sqrt(((xyz - xyz.mean(axis=1, keepdims=True)) ** 2).sum(axis=2).mean(axis=1))
    chi = susceptibility(A, E, T0, n_blocks=4)
    tol = 1e-9 * max(1e-30, abs(chi['chi_total']))
    assert abs(chi['chi_Boltzmann'] + chi['chi_dielectric'] - chi['chi_total']) < tol
    assert abs(chi['chi_AH'] + chi['chi_bond'] + chi['chi_DH'] - chi['chi_total']) < tol
    M_AH, M_DH, chi_0 = pairwise_susceptibility_map(A, xyz, q, lam, sig, b0, get_dielectric_params(T0))
    pair_sum = np.triu(M_AH, 2).sum() + np.triu(M_DH, 2).sum() + chi_0
    assert abs(pair_sum - chi['chi_total']) < 1e-8 * abs(chi['chi_total']), "pairwise sum rule fails"
    say("  chi_Boltzmann + chi_dielectric = chi_AH + chi_bond + chi_DH = pair-map sum = chi_total ✓")

    # 9: two-state system (Derivation §14): d<A>/dT = Cov(A, H)/(R T^2) for a fixed landscape
    beta_t = 1.0 / (R_GAS * T0)
    dH = 2.0
    p1 = 1.0 / (1.0 + np.exp(beta_t * dH))
    p0 = 1.0 - p1
    direct = (1.0 / (1.0 + np.exp(dH / (R_GAS * (T0 + h))))
              - 1.0 / (1.0 + np.exp(dH / (R_GAS * (T0 - h))))) / (2 * h)
    assert abs(direct - p0 * p1 * dH / (R_GAS * T0 ** 2)) < 1e-9, "two-state check failed"
    say(f"  Two-state validation: d<A>/dT = {direct:.6e} K^-1 ✓")

    say("\nAll sanity checks passed ✓\n")
    say("Dielectric parameters at the five simulation temperatures:")
    say(f"{'T (K)':>8} {'eps_w':>8} {'l_B (nm)':>10} {'kappa (nm^-1)':>14} {'kappa^-1 (nm)':>14}")
    say("-" * 58)
    for T in TEMPS:
        p = get_dielectric_params(float(T))
        say(f"{T:>8.0f} {p['eps_w']:>8.3f} {p['l_B']:>10.4f} {p['kappa']:>14.4f} {1 / p['kappa']:>14.4f}")


if __name__ == '__main__':
    run_sanity_checks()
    print()
    print("The runner is TODO_T2_trajectory_loader.py (modes: check, run, all, validate).")

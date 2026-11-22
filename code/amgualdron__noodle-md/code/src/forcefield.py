"""The physics of the model: residue parameters, electrostatics, and every force.

One bead per amino acid (the HPS model of Dignon et al. 2018). Units throughout:
    length nm, energy kcal/mol, mass g/mol, charge in elementary charges e.
HOOMD's unit of time follows from these: sqrt(g/mol * nm^2 / (kcal/mol)) = 0.489 ps,
so a time step dt = 0.01 is about 4.9 fs.

How HOOMD forces work, for readers who have not used it: every force (between residue
pairs, along bonds, or from a wall) is an object holding its parameters per particle
type. There is one particle type per amino acid, named by its letter ("A", "R", ...),
so `force.params["K"]` holds lysine's parameters and `force.params[("K", "E")]` those of
the lysine-glutamate pair. simulate.py hands the list of force objects returned here to
the integrator, which adds them all up at every step.

Nothing in here reads files or knows about runs; simulate.py calls these functions
with numbers taken from the run's physics.yaml.

HOOMD is imported only inside the functions that build forces: loading it takes a while
(up to minutes on a cluster's shared filesystem), and the residue table and formulas are
also used by commands that never simulate (noodle analyze, export, clean, launch's plan).
"""

import numpy as np


# ── Residue table ─────────────────────────────────────────────────────────────
# mass (g/mol), charge (e), sigma = bead diameter (nm), and lambda = hydropathy, from
# 0 (hydrophilic, no attraction) to 1 (most hydrophobic). Dignon et al. 2018.
AA_PARAMS = {
    "A": dict(mass=71.08,  charge=0.0,  sigma=0.504, hydropathy=0.730),
    "R": dict(mass=156.20, charge=1.0,  sigma=0.656, hydropathy=0.000),
    "N": dict(mass=114.10, charge=0.0,  sigma=0.568, hydropathy=0.432),
    "D": dict(mass=115.10, charge=-1.0, sigma=0.558, hydropathy=0.378),
    "C": dict(mass=103.10, charge=0.0,  sigma=0.548, hydropathy=0.595),
    "Q": dict(mass=128.10, charge=0.0,  sigma=0.602, hydropathy=0.514),
    "E": dict(mass=129.10, charge=-1.0, sigma=0.592, hydropathy=0.459),
    "G": dict(mass=57.05,  charge=0.0,  sigma=0.450, hydropathy=0.649),
    "H": dict(mass=137.10, charge=0.5,  sigma=0.608, hydropathy=0.514),
    "I": dict(mass=113.20, charge=0.0,  sigma=0.618, hydropathy=0.973),
    "L": dict(mass=113.20, charge=0.0,  sigma=0.618, hydropathy=0.973),
    "K": dict(mass=128.20, charge=1.0,  sigma=0.636, hydropathy=0.514),
    "M": dict(mass=131.20, charge=0.0,  sigma=0.618, hydropathy=0.838),
    "F": dict(mass=147.20, charge=0.0,  sigma=0.636, hydropathy=1.000),
    "P": dict(mass=97.12,  charge=0.0,  sigma=0.556, hydropathy=1.000),
    "S": dict(mass=87.08,  charge=0.0,  sigma=0.518, hydropathy=0.595),
    "T": dict(mass=101.10, charge=0.0,  sigma=0.562, hydropathy=0.676),
    "W": dict(mass=186.20, charge=0.0,  sigma=0.678, hydropathy=0.946),
    "Y": dict(mass=163.20, charge=0.0,  sigma=0.646, hydropathy=0.865),
    "V": dict(mass=99.07,  charge=0.0,  sigma=0.586, hydropathy=0.892),
}
AMINO_ACIDS = list(AA_PARAMS)  # particle type order in every simulation


# ── Electrostatics ────────────────────────────────────────────────────────────
BOLTZMANN = 0.001987204259             # kcal/mol/K
ELEMENTARY_CHARGE = 1.602176634e-19    # C
VACUUM_PERMITTIVITY = 8.8541878128e-12  # F/m
AVOGADRO = 6.02214076e23               # 1/mol
JOULES_PER_KCAL = 4184.0
NM_PER_M = 1e9
NM3_PER_LITRE = 1e24

# e^2 / (4 pi eps0): the Coulomb energy of two unit charges 1 nm apart in vacuum,
# converted from J*m to kcal/mol*nm (about 332 kcal/mol)
COULOMB_KCAL_NM = (ELEMENTARY_CHARGE**2 / (4 * np.pi * VACUUM_PERMITTIVITY)
                   * AVOGADRO / JOULES_PER_KCAL * NM_PER_M)


def thermal_energy(T_kelvin):
    """kT in kcal/mol."""
    return BOLTZMANN * T_kelvin


def water_dielectric(T_kelvin):
    """Relative permittivity of water, empirical fit (Akerlof & Oshry 1950); ~78.5 at 298 K."""
    T = T_kelvin
    return 5321 / T + 233.76 - 0.9297 * T + 1.417e-3 * T**2 - 8.292e-7 * T**3


def debye_huckel(T_kelvin, ionic_concentration_M):
    """Return (bjerrum_length_nm, kappa_per_nm): the two numbers of screened electrostatics.

    Bjerrum length l_B = e^2 / (4 pi eps0 eps_water kT): the distance at which two unit
        charges in water interact with an energy of kT (~0.71 nm at room temperature).
    kappa = 1 / Debye length: how quickly salt screens charges. For a 1:1 salt at
        concentration c, kappa^2 = 8 pi l_B n, with n = c in ions per nm^3.
    Every charged interaction in the model (residue pairs and the charged wall) uses
    these same two numbers.
    """
    bjerrum = COULOMB_KCAL_NM / (water_dielectric(T_kelvin) * thermal_energy(T_kelvin))
    ions_per_nm3 = ionic_concentration_M * AVOGADRO / NM3_PER_LITRE
    kappa = np.sqrt(8 * np.pi * bjerrum * ions_per_nm3)
    return bjerrum, kappa


# ── Residue-residue forces ────────────────────────────────────────────────────
def residue_pair_forces(physics, nlist):
    """HPS non-bonded forces between all residue pairs, plus the chain bonds.

    Short range (Ashbaugh-Hatch): a Lennard-Jones well whose attraction is scaled by
    the average hydropathy lambda of the two residues:
        U = LJ(r) + (1 - lambda) epsilon   inside the LJ minimum (r < 2^(1/6) sigma)
        U = lambda * LJ(r)                 outside it, up to lj_r_cut
    HOOMD has no such potential, so it is built from two LJ forces that add up to it:
        lj_core: LJ with epsilon*(1 - lambda), cut at the minimum and shifted up to zero
                 there, so it is purely repulsive
        lj_tail: ordinary LJ with epsilon*lambda, out to lj_r_cut
    Long range: screened Coulomb (Debye-Huckel, "Yukawa" in HOOMD),
        U = q_i q_j l_B kT exp(-kappa r) / r,   shifted to zero at yukawa_r_cut.
    Bonds: harmonic springs, U = k/2 (r - r0)^2. The `nlist` (neighbour list, HOOMD's
    list of which beads are close enough to interact) is told to skip bonded pairs.

    Speed: the pair forces take ~85% of the run time, and HOOMD evaluates each force for
    every pair within its r_cut, even when the result is zero. So the screened Coulomb
    term is switched off (r_cut = 0) for every pair that is not charged on both sides:
    its energy there is exactly zero, and skipping it makes runs ~1.5x faster. With
    lj_r_cut below yukawa_r_cut, the neighbour list of those pairs is also shorter.
    """
    import hoomd.md  # (see the top of this file)
    epsilon = physics["epsilon_lj"]
    kT = thermal_energy(physics["temperature_K"])
    bjerrum, kappa = debye_huckel(physics["temperature_K"], physics["ionic_concentration_M"])

    lj_core = hoomd.md.pair.LJ(nlist=nlist, mode="shift")
    lj_tail = hoomd.md.pair.LJ(nlist=nlist)
    yukawa = hoomd.md.pair.Yukawa(nlist=nlist, mode="shift")

    for i, a in enumerate(AMINO_ACIDS):
        for b in AMINO_ACIDS[i:]:  # every pair once: (A, A), (A, R), ... (V, V)
            sigma = (AA_PARAMS[a]["sigma"] + AA_PARAMS[b]["sigma"]) / 2
            lam = (AA_PARAMS[a]["hydropathy"] + AA_PARAMS[b]["hydropathy"]) / 2
            lj_core.params[(a, b)] = dict(epsilon=epsilon * (1 - lam), sigma=sigma)
            lj_core.r_cut[(a, b)] = 2 ** (1 / 6) * sigma
            lj_tail.params[(a, b)] = dict(epsilon=epsilon * lam, sigma=sigma)
            lj_tail.r_cut[(a, b)] = physics["lj_r_cut"]
            charge_product = AA_PARAMS[a]["charge"] * AA_PARAMS[b]["charge"]
            yukawa.params[(a, b)] = dict(epsilon=bjerrum * kT * charge_product, kappa=kappa)
            yukawa.r_cut[(a, b)] = physics["yukawa_r_cut"] if charge_product != 0 else 0.0

    bonds = hoomd.md.bond.Harmonic()
    bonds.params["AA-bond"] = dict(k=physics["harmonic_k"], r0=physics["harmonic_r0"])
    return [lj_core, lj_tail, yukawa, bonds]


# ── Wall forces (grafted geometry) ────────────────────────────────────────────
# The wall is a flat surface at height z_wall, treated as a smooth continuum.
#   lj93:  integrating LJ sites over a solid half-space gives the 9-3 potential,
#          the standard form for a featureless surface.
#   lj126: the plain 12-6 LJ form, as between two residues.
# HOOMD's Mie potential covers both (exponents n, m) and normalises the well depth to
# epsilon. Every HOOMD wall potential is shifted so that U(r_cut) = 0.
WALL_EXPONENTS = {"lj93": (9, 3), "lj126": (12, 6)}


def wall_contact_distance(sigma_residue, physics):
    """Distance from the wall of the potential minimum for a residue of this size.

    The purely repulsive wall is cut here, so a bead at this height feels no force
    at all; it is also where the grafted anchor bead is placed.
    """
    n, m = WALL_EXPONENTS[physics["wall_potential"]]
    sigma = (sigma_residue + physics["wall_sigma"]) / 2
    return (n / m) ** (1 / (n - m)) * sigma


def wall_forces(physics, z_wall, z_ceiling):
    """Residue-wall forces: the surface below the chain and a ceiling above it.

    Same Ashbaugh-Hatch split as the residue pairs, with wall_attraction playing the
    role of lambda:
        wall_attraction = 0  ->  purely repulsive wall: the default
        wall_attraction > 0  ->  attractive tail of depth epsilon*lambda_i, where
            lambda_i = wall_attraction * (residue hydropathy)   for wall_attraction_scale: hydropathy
            lambda_i = wall_attraction                          for wall_attraction_scale: uniform
    HOOMD shifts the attractive tail to zero at wall_r_cut, which makes the well
    ~1% shallower than epsilon*lambda_i (forces are unaffected).
    A charged surface (wall_charge_density != 0) adds screened electrostatics.
    The ceiling only exists to stop beads wrapping through the periodic box in z; it is
    always purely repulsive and placed out of reach of the chain.
    """
    import hoomd.md  # (see the top of this file)
    n, m = WALL_EXPONENTS[physics["wall_potential"]]
    epsilon = physics["wall_epsilon"]
    # A Plane wall acts on the beads on the side its normal points to
    surface = hoomd.wall.Plane(origin=(0, 0, z_wall), normal=(0, 0, 1))
    ceiling = hoomd.wall.Plane(origin=(0, 0, z_ceiling), normal=(0, 0, -1))

    core_surface = hoomd.md.external.wall.Mie(walls=[surface])
    core_ceiling = hoomd.md.external.wall.Mie(walls=[ceiling])
    tail_surface = hoomd.md.external.wall.Mie(walls=[surface])

    for aa in AMINO_ACIDS:
        sigma = (AA_PARAMS[aa]["sigma"] + physics["wall_sigma"]) / 2
        r_min = wall_contact_distance(AA_PARAMS[aa]["sigma"], physics)
        if physics["wall_attraction_scale"] == "hydropathy":
            lam = physics["wall_attraction"] * AA_PARAMS[aa]["hydropathy"]
        elif physics["wall_attraction_scale"] == "uniform":
            lam = physics["wall_attraction"]
        else:
            raise ValueError("wall_attraction_scale must be 'hydropathy' or 'uniform'")

        # r_extrap: a bead closer than this to a wall (or past it) is pushed back by a
        # linear extrapolation of the force, instead of feeling a near-infinite force.
        core_surface.params[aa] = dict(epsilon=epsilon * (1 - lam), sigma=sigma, n=n, m=m,
                                       r_cut=r_min, r_extrap=0.1 * sigma)
        core_ceiling.params[aa] = dict(epsilon=epsilon, sigma=sigma, n=n, m=m,
                                       r_cut=r_min, r_extrap=0.1 * sigma)
        # r_cut = 0 switches a type off entirely (all of them when wall_attraction = 0)
        tail_surface.params[aa] = dict(epsilon=epsilon * lam, sigma=sigma, n=n, m=m,
                                       r_cut=physics["wall_r_cut"] if lam > 0 else 0.0,
                                       r_extrap=0.0)

    forces = [core_surface, core_ceiling]
    if physics["wall_attraction"] > 0:
        forces.append(tail_surface)
    if physics["wall_charge_density"] != 0:
        forces.append(charged_surface(physics, surface, reach=z_ceiling - z_wall))
    return forces


def charged_surface(physics, surface, reach):
    """Screened electrostatics of a uniformly charged surface (Debye-Huckel).

    Linearised Poisson-Boltzmann for a plane carrying sigma_s charges per nm^2, with
    salt water on one side only (a solid substrate), gives a residue of charge q at
    height h above the surface the energy
        U(h) = A exp(-kappa h),    A = 4 pi l_B kT q sigma_s / kappa,
    with the same l_B and kappa as the residue pairs, so the surface is screened
    exactly like the chain. Image charges from the dielectric jump at the surface are
    neglected, as usual at this level of coarse-graining.

    HOOMD has no wall potential of exactly this form, but its Morse wall contains it:
        U_Morse(h) = D0 [exp(-2 alpha (h - r0)) - 2 exp(-alpha (h - r0))]
    With alpha = kappa and r0 = -OFFSET / kappa, a point far *below* the surface:
        second term = -2 D0 exp(-OFFSET) exp(-kappa h),  which is A exp(-kappa h)
                      when D0 = -(A / 2) exp(OFFSET);
        first term  = -(A / 2) exp(-OFFSET) exp(-2 kappa h),  smaller than the second
                      by a factor exp(-OFFSET) / 2 ~ 5e-14: zero for all purposes.
    The cutoff is the ceiling, out of reach of every bead. (HOOMD's shift to zero
    there changes the energy by a negligible constant and the force not at all.)
    check_physics.py compares the result with the formula above. Running this inside
    HOOMD, rather than as a Python force, makes charged-surface runs ~3x faster.
    """
    import hoomd.md  # (see the top of this file)
    OFFSET = 30
    kT = thermal_energy(physics["temperature_K"])
    bjerrum, kappa = debye_huckel(physics["temperature_K"], physics["ionic_concentration_M"])

    force = hoomd.md.external.wall.Morse(walls=[surface])
    for aa in AMINO_ACIDS:
        q = AA_PARAMS[aa]["charge"]
        A = 4 * np.pi * bjerrum * kT * q * physics["wall_charge_density"] / kappa
        force.params[aa] = dict(D0=-A / 2 * np.exp(OFFSET), alpha=kappa, r0=-OFFSET / kappa,
                                r_cut=reach if q != 0 else 0.0, r_extrap=0.0)
    return force

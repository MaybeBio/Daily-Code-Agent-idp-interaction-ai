"""What a finished task measured: averages with error bars, from its trajectory.

Every job calls analyse_task() right after each of its tasks finishes, on the compute
node, and writes the averages into the task's row of its job file. ./noodle analyze then
only has to gather and pool these rows.

Every column, with its exact definition, units and error bar, is described in
OBSERVABLES.md. In short (lengths in nm, times in model time units tu = 0.489 ps):

  every chain           rg_nm, rg_sq_nm2, rg_sd_nm, rh_nm, ree_nm, ree_sq_nm2, asphericity,
                        shape_anisotropy, prolateness, transverse_rms_nm, nu, r0_nm, bond_cos,
                        persistence_nm, contacts_per_residue
  grafted and brushes   com_height_nm, end_height_nm, max_height_nm, height_95_nm, rg_perp_nm,
                        rg_par_nm, end_lateral_rms_nm, tilt_cos, bond_order, adsorbed_fraction
  brushes               brush_height_nm, concentration_mg_ml, interchain_contacts_per_residue,
                        lateral_heterogeneity
  diagnostics           tau_rg_tu, tau_ree_tu, independent_samples, drift

Each observable X also has X_err, one standard error from block averaging (BLOCKS blocks).

Curves (density against height, distances against separation along the chain, ...) are
not computed here: with thousands of tasks they would be far larger than everything
else. They can be computed later from the trajectories (OBSERVABLES.md, last section).
"""

import math

import gsd.hoomd
import numpy as np
from scipy.spatial import cKDTree

from forcefield import AA_PARAMS, wall_contact_distance
from runs import DIAGNOSTICS, OBSERVABLES

BLOCKS = 10               # blocks for the error bars (and for nu's error bar)
CHUNK_FRAMES = 50         # frames read from the trajectory at a time (keeps memory small)
CONTACT_FACTOR = 1.5      # a contact: closer than this x the mean diameter of the two residues
ADSORBED_MARGIN = 0.5     # nm: a residue touches the surface within its contact distance + this
NU_MIN_SEPARATION = 5     # nu is fitted for |i - j| >= this
PAIRS_AT_ONCE = 1_000_000  # residue pairs held in memory at once (~8 MB per array)
MG_ML_PER_DALTON_NM3 = 1.66054  # 1 dalton per nm^3 is 1.66054 mg/mL


# ── Small statistics ──────────────────────────────────────────────────────────
def block_error(series):
    """Standard error of the mean of a correlated series: cut it into BLOCKS consecutive
    blocks and use the scatter of the block averages."""
    n = len(series) // BLOCKS * BLOCKS
    if n == 0:
        return math.nan
    means = np.asarray(series[:n]).reshape(BLOCKS, -1).mean(axis=1)
    return float(means.std(ddof=1) / math.sqrt(BLOCKS))


def autocorrelation(x):
    """Normalised autocorrelation function of x (frames first, then chains and, for a
    vector, its components), averaged over the chains; vectors use the dot product."""
    x = x - x.mean(axis=0)
    n = len(x)
    spectrum = np.fft.rfft(x, 2 * n, axis=0)
    covariance = np.fft.irfft(spectrum * np.conj(spectrum), axis=0)[:n]
    covariance /= (n - np.arange(n)).reshape((n,) + (1,) * (x.ndim - 1))  # the same weight at every lag
    covariance = covariance.reshape(n, -1).mean(axis=1)
    return covariance / covariance[0] if covariance[0] > 0 else np.full(n, np.nan)


def correlation_frames(x):
    """Integrated autocorrelation time in frames: 1/2 + the sum of the autocorrelation
    function up to where it first drops below zero (where noise starts to dominate)."""
    acf = autocorrelation(x)
    if np.isnan(acf[0]):
        return math.nan
    first_negative = int(np.argmax(acf < 0)) or len(acf)
    return 0.5 + float(acf[1:first_negative].sum())


def drifts(series):
    """True if the two halves of a series differ by more than 3 standard errors."""
    half = len(series) // 2
    first, second = series[:half], series[half:2 * half]
    error = math.hypot(block_error(first), block_error(second))
    return bool(error > 0 and abs(second.mean() - first.mean()) > 3 * error)


def fit_nu(separations, r2):
    """(nu, r0) of sqrt(R2) = r0 s^nu, least squares in log-log over s >= NU_MIN_SEPARATION."""
    use = (separations >= NU_MIN_SEPARATION) & (r2 > 0)
    if use.sum() < 3:
        return math.nan, math.nan
    nu, log_r0 = np.polyfit(np.log(separations[use]), 0.5 * np.log(r2[use]), 1)
    return float(nu), float(math.exp(log_r0))


# ── The trajectory ────────────────────────────────────────────────────────────
def read_chunks(path, chains, n):
    """Yield (first frame, positions) for CHUNK_FRAMES frames at a time. positions:
    [frame, chain, residue, xyz], each chain made whole again: the file holds positions
    wrapped into the periodic box, and every bond is shorter than half the box, so its true
    vector is its minimum-image vector and residue i sits at residue 0 plus the first i bonds.
    Also yields the raw (wrapped) positions, which the brush's lateral grid needs."""
    with gsd.hoomd.open(path, "r") as f:
        box = np.asarray(f[0].configuration.box[:3], dtype=np.float64)
        for start in range(0, len(f), CHUNK_FRAMES):
            stop = min(start + CHUNK_FRAMES, len(f))
            raw = np.stack([f[i].particles.position for i in range(start, stop)]).astype(np.float64)
            raw = raw.reshape(len(raw), chains, n, 3)
            bonds = np.diff(raw, axis=2)
            bonds -= box * np.round(bonds / box)
            whole = np.concatenate([raw[:, :, :1], raw[:, :, :1] + np.cumsum(bonds, axis=2)], axis=2)
            yield start, whole, raw


def frame_count(path):
    with gsd.hoomd.open(path, "r") as f:
        return len(f), np.asarray(f[0].configuration.box[:3], dtype=np.float64)


# ── One task ──────────────────────────────────────────────────────────────────
def analyse_task(trajectory, sequence, physics, geometry, frame_time):
    """The averages {column: value} of one finished task.

    geometry: simulate.system_geometry() of the task (chains, z_wall, ...).
    frame_time: model time units between two saved frames."""
    n = len(sequence)
    chains = geometry.get("chains", 1)
    grafted = physics["geometry"] == "grafted"
    brush = chains > 1
    frames, box = frame_count(trajectory)
    if frames < 2 * BLOCKS:
        raise ValueError(f"only {frames} frames: at least {2 * BLOCKS} are needed")

    masses = np.array([AA_PARAMS[aa]["mass"] for aa in sequence])
    sigmas = np.array([AA_PARAMS[aa]["sigma"] for aa in sequence])
    weights = masses / masses.sum()
    index = np.arange(n)
    separation = np.abs(index[:, None] - index[None, :])
    contact_cutoff_sq = (CONTACT_FACTOR * (sigmas[:, None] + sigmas[None, :]) / 2) ** 2
    intrachain_pair = separation >= 3
    block_of_frame = np.arange(frames) * BLOCKS // frames
    # chains whose residue pairs are handled together: about PAIRS_AT_ONCE pairs (memory)
    chain_groups = np.array_split(np.arange(chains), math.ceil(chains * n * n / PAIRS_AT_ONCE))

    # per frame and chain: series[name] is a list of [frames_in_chunk, chains] arrays
    series = {}
    # accumulated over frames (and chains), per block where an error bar needs it
    r2_by_block = np.zeros((BLOCKS, n))          # sums of R_ij^2 by separation
    pairs_by_block = np.zeros((BLOCKS, n))       # how many pairs went into them
    ree_vectors = []                             # [frames, chains, 3], for its relaxation time

    if grafted:
        free_end = n - 1 if physics["graft_terminus"] == "N" else 0
        anchor = n - 1 - free_end
        contact_height = np.array([wall_contact_distance(s, physics) for s in sigmas])
        not_anchor = index != anchor

    def add(name, values):
        series.setdefault(name, []).append(values)

    for start, x, raw in read_chunks(trajectory, chains, n):
        blocks = block_of_frame[start:start + len(x)]
        # ── shape: the gyration tensor of each chain in each frame ──
        com = np.einsum("i,fcij->fcj", weights, x)
        rel = x - com[:, :, None, :]
        gyration = np.einsum("i,fcij,fcik->fcjk", weights, rel, rel)
        moments = np.linalg.eigvalsh(gyration)   # ascending
        rg_sq = moments.sum(axis=-1)
        deviation = moments - rg_sq[..., None] / 3
        add("rg_sq_nm2", rg_sq)
        add("asphericity", (moments[..., 2] - 0.5 * (moments[..., 0] + moments[..., 1])) / rg_sq)
        add("shape_anisotropy", 1.5 * np.sum(deviation ** 2, axis=-1) / rg_sq ** 2)
        add("prolateness", 27 * np.prod(deviation, axis=-1) / rg_sq ** 3)

        # ── ends and bonds ──
        ree = x[:, :, -1] - x[:, :, 0]
        ree_length = np.linalg.norm(ree, axis=-1)
        ree_vectors.append(ree)
        add("ree_nm", ree_length)
        add("ree_sq_nm2", ree_length ** 2)
        axis = ree / np.maximum(ree_length, 1e-12)[..., None]
        from_start = x - x[:, :, :1]
        along = np.einsum("fcij,fcj->fci", from_start, axis)
        add("transverse_rms_nm", np.sqrt(np.mean(np.sum(from_start ** 2, axis=-1) - along ** 2, axis=-1).clip(0)))
        bonds = np.diff(x, axis=2)
        unit_bonds = bonds / np.linalg.norm(bonds, axis=-1, keepdims=True)
        add("bond_cos", np.mean(np.sum(unit_bonds[:, :, 1:] * unit_bonds[:, :, :-1], axis=-1), axis=-1))
        add("persistence_nm", 0.5 * (np.einsum("fcj,fcj->fc", ree, unit_bonds[:, :, 0])
                                     + np.einsum("fcj,fcj->fc", ree, unit_bonds[:, :, -1])))

        # ── pairs within each chain, one frame (and, for a brush, a group of chains) at a time ──
        inverse_mean, contacts = np.zeros(rg_sq.shape), np.zeros(rg_sq.shape)
        for k in range(len(x)):
            for c in chain_groups:
                r = rel[k, c]                                           # [chain, residue, xyz], centred
                squared = np.sum(r ** 2, axis=-1)
                d2 = (squared[:, :, None] + squared[:, None, :] - 2 * np.einsum("cix,cjx->cij", r, r)).clip(0)
                r2_by_block[blocks[k]] += np.bincount(separation.ravel(), weights=d2.sum(axis=0).ravel(), minlength=n)
                pairs_by_block[blocks[k]] += np.bincount(separation.ravel(), minlength=n) * len(r)
                distance = np.sqrt(d2)
                np.einsum("cii->ci", distance)[:] = np.inf           # leave out i = j
                inverse_mean[k, c] = np.sum(1 / distance, axis=(1, 2)) / n**2
                touching = (d2 < contact_cutoff_sq) & intrachain_pair
                contacts[k, c] = touching.sum(axis=(1, 2)) / n
        add("_inverse_distance", inverse_mean)
        add("contacts_per_residue", contacts)

        if not grafted:
            continue
        # ── heights above the surface (z is not periodic for grafted chains) ──
        height = x[..., 2] - geometry["z_wall"]                           # [frame, chain, residue]
        add("com_height_nm", com[..., 2] - geometry["z_wall"])
        add("end_height_nm", height[:, :, free_end])
        add("max_height_nm", height.max(axis=-1))
        add("height_95_nm", np.percentile(height, 95, axis=-1))
        add("_rg_perp_sq", gyration[..., 2, 2])
        add("_rg_par_sq", gyration[..., 0, 0] + gyration[..., 1, 1])
        lateral = x[:, :, free_end, :2] - x[:, :, anchor, :2]
        add("_end_lateral_sq", np.sum(lateral ** 2, axis=-1))
        up = x[:, :, free_end] - x[:, :, anchor]
        add("tilt_cos", up[..., 2] / np.maximum(np.linalg.norm(up, axis=-1), 1e-12))
        add("bond_order", np.mean(1.5 * unit_bonds[..., 2] ** 2 - 0.5, axis=-1))
        touching = height < contact_height + ADSORBED_MARGIN
        add("adsorbed_fraction", touching[:, :, not_anchor].mean(axis=-1))

        if not brush:
            continue
        # ── the brush as a whole, one frame at a time ──
        area = box[0] * box[1]
        side = math.isqrt(chains)
        brush_height, concentration = np.zeros((len(x), 1)), np.zeros((len(x), 1))
        interchain, heterogeneity = np.zeros((len(x), 1)), np.zeros((len(x), 1))
        chain_of = np.repeat(np.arange(chains), n)
        sigma_all = np.tile(sigmas, chains)
        for k in range(len(x)):
            h = height[k].ravel()
            brush_height[k] = 2 * h.mean()
            top = np.percentile(h, 95)
            concentration[k] = np.tile(masses, chains)[h < top].sum() / (area * top) * MG_ML_PER_DALTON_NM3
            # contacts between chains: a tree of all residues, periodic in x and y
            points = raw[k].reshape(-1, 3) - np.array([-box[0] / 2, -box[1] / 2, geometry["z_wall"]])
            lateral = points[:, :2]
            lateral %= box[:2]
            lateral[lateral >= box[:2]] = 0  # -1e-17 % L gives L, which the tree refuses
            points[:, 2] = points[:, 2].clip(0)
            tree = cKDTree(points, boxsize=[box[0], box[1], 10 * box[2]])
            pairs = tree.query_pairs(CONTACT_FACTOR * sigmas.max(), output_type="ndarray")
            if len(pairs):
                i, j = pairs[:, 0], pairs[:, 1]
                d = points[i] - points[j]
                d[:, :2] -= box[:2] * np.round(d[:, :2] / box[:2])
                keep = (chain_of[i] != chain_of[j]) & (np.sum(d ** 2, axis=1)
                                                         < (CONTACT_FACTOR * (sigma_all[i] + sigma_all[j]) / 2) ** 2)
                interchain[k] = 2 * keep.sum() / (chains * n)
            counts = np.histogram2d(points[:, 0], points[:, 1], bins=side, range=[[0, box[0]], [0, box[1]]])[0]
            heterogeneity[k] = counts.std() / counts.mean()
        add("brush_height_nm", brush_height)
        add("concentration_mg_ml", concentration)
        add("interchain_contacts_per_residue", interchain)
        add("lateral_heterogeneity", heterogeneity)

    series = {name: np.concatenate(parts) for name, parts in series.items()}  # [frames, chains] each

    # ── averages and error bars ──
    averages = {}

    def average(name, values):
        averages[name] = float(values.mean())
        averages[f"{name}_err"] = block_error(values.mean(axis=1))

    rg = np.sqrt(series["rg_sq_nm2"])
    average("rg_nm", rg)
    average("rg_sq_nm2", series["rg_sq_nm2"])
    averages["rg_sd_nm"], averages["rg_sd_nm_err"] = float(rg.std()), math.nan
    inverse = series["_inverse_distance"]
    averages["rh_nm"] = 1 / float(inverse.mean())
    averages["rh_nm_err"] = block_error(inverse.mean(axis=1)) / float(inverse.mean()) ** 2
    for name in ("ree_nm", "ree_sq_nm2", "asphericity", "shape_anisotropy", "prolateness", "transverse_rms_nm"):
        average(name, series[name])
    r2 = r2_by_block.sum(axis=0) / np.maximum(pairs_by_block.sum(axis=0), 1)
    averages["nu"], averages["r0_nm"] = fit_nu(index, r2)
    fits = [fit_nu(index, r2_by_block[b] / np.maximum(pairs_by_block[b], 1)) for b in range(BLOCKS)]
    averages["nu_err"] = float(np.std([f[0] for f in fits], ddof=1) / math.sqrt(BLOCKS))
    averages["r0_nm_err"] = float(np.std([f[1] for f in fits], ddof=1) / math.sqrt(BLOCKS))
    for name in ("bond_cos", "persistence_nm", "contacts_per_residue"):
        average(name, series[name])
    if grafted:
        for name in ("com_height_nm", "end_height_nm", "max_height_nm", "height_95_nm", "tilt_cos",
                     "bond_order", "adsorbed_fraction"):
            average(name, series[name])
        for name, squared in (("rg_perp_nm", "_rg_perp_sq"), ("rg_par_nm", "_rg_par_sq"),
                              ("end_lateral_rms_nm", "_end_lateral_sq")):
            mean = float(series[squared].mean())
            averages[name] = math.sqrt(mean)
            averages[f"{name}_err"] = block_error(series[squared].mean(axis=1)) / (2 * math.sqrt(mean)) if mean > 0 else math.nan
    if brush:
        for name in ("brush_height_nm", "concentration_mg_ml", "interchain_contacts_per_residue", "lateral_heterogeneity"):
            average(name, series[name])

    # ── diagnostics ──
    tau_rg = correlation_frames(rg)
    tau_ree = correlation_frames(np.concatenate(ree_vectors))
    averages.update(tau_rg_tu=tau_rg * frame_time, tau_ree_tu=tau_ree * frame_time,
                    independent_samples=frames / (2 * tau_rg) if tau_rg > 0 else math.nan,
                    drift=int(drifts(rg.mean(axis=1)) or (grafted and drifts(series["com_height_nm"].mean(axis=1)))))

    # every column, in order (None where it does not apply, e.g. heights of a free chain)
    columns = [c for name in OBSERVABLES for c in (name, f"{name}_err")] + DIAGNOSTICS
    return {c: averages.get(c) for c in columns}

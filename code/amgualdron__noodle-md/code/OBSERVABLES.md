# noodle-md observables

The exact meaning of every column in the two results tables, `manifest.csv` and
`analysis/replicates.csv`, and how each number is computed from a simulation.
The last section shows how to compute anything else (curves such as a density profile)
from the saved trajectories.

**Contents**

1. [The two tables](#1-the-two-tables)
2. [Notation](#2-notation)
3. [Which protein, which conditions](#3-which-protein-which-conditions)
4. [Sequence descriptors](#4-sequence-descriptors)
5. [Size and shape of a chain](#5-size-and-shape-of-a-chain)
6. [Chains attached to a surface](#6-chains-attached-to-a-surface)
7. [Brushes](#7-brushes)
8. [Error bars](#8-error-bars)
9. [Quality checks](#9-quality-checks)
10. [From replicates to the manifest](#10-from-replicates-to-the-manifest)
11. [Computing other quantities from the trajectories](#11-computing-other-quantities-from-the-trajectories)

---

## 1. The two tables

| table | one row per | what it is for |
|---|---|---|
| `manifest.csv` | protein × physics set | **the main result**: each observable averaged over the replicates |
| `analysis/replicates.csv` | task (one protein, one physics set, one replicate) | each simulation on its own, to see how much replicates differ |

A **task** is one simulation: one protein (or one brush of it) under one physics set.
**Replicates** are tasks that are identical except for their random seed. Every
observable `X` comes with `X_err`, its error bar (one standard error, section 8).

An **empty cell** means "does not apply" (heights of a free chain, brush columns of a
single chain) or "could not be computed" (for example ν for a chain too short to fit,
section 5). Numbers are written with 6 significant digits.

Columns of `manifest.csv`, in order: identity and conditions (section 3), sequence
descriptors (section 4), `replicates`, the observables of sections 5–7 each followed by
its `_err`, `reduced_density`, and the quality checks (section 9).

Columns of `analysis/replicates.csv`: `task_id`, `name`, `physics`, `replicate`, `seed`,
`N`, `chains`, `frames`, `run_time_h` (section 3), the observables with their `_err`,
and the quality checks.

---

## 2. Notation

**Units.** Lengths in nm, energies in kcal/mol, masses in g/mol (Da), charges in units of
the elementary charge e. Times in the model's time unit **tu** = 0.489 ps. Model times
compare tasks with each other; they are not real times (the friction is a numerical
choice, not water's).

**A chain** has N residues, numbered i = 1 … N from the N-terminus. Residue i has
position **r**ᵢ, mass mᵢ, charge qᵢ and diameter σᵢ (the model's values for its amino
acid; histidine has charge +0.5). The mass fractions are wᵢ = mᵢ / Σⱼ mⱼ. The centre of
mass is **R**cm = Σᵢ wᵢ **r**ᵢ. The distance between residues i and j is rᵢⱼ = |**r**ᵢ − **r**ⱼ|.

**Bonds.** **b**ᵢ = **r**ᵢ₊₁ − **r**ᵢ (i = 1 … N − 1), with unit vectors **b̂**ᵢ = **b**ᵢ / |**b**ᵢ|.
The end-to-end vector is **R**ee = **r**N − **r**₁.

**Whole chains.** The simulation box is periodic, and the saved positions are wrapped
into it, so a chain can appear cut in pieces. Before anything is measured each chain is
made whole: every bond is much shorter than half the box, so its true vector is the
shortest one between its two residues through the periodic boundaries, and residue i is
placed at residue 1 plus the first i − 1 bonds.

**Averages.** ⟨ … ⟩ is the average over every saved frame of the production run (the
equilibration is never saved) and, in a brush, over every chain. A quantity is first
computed per frame (and per chain), then averaged. F is the number of frames.

**Heights** (chains attached to a surface). The surface is the plane z = z_wall; the
height of residue i is hᵢ = zᵢ − z_wall. The **anchor** is the attached residue (residue 1
when grafted by the N-terminus, residue N by the C-terminus); it never moves. The **free
end** is the other terminal residue.

**Contacts.** Two residues are in contact when rᵢⱼ < 1.5 × (σᵢ + σⱼ)/2, i.e. closer than
1.5 times their mean diameter (0.7 to 1 nm, depending on the residues).

---

## 3. Which protein, which conditions

| column | meaning |
|---|---|
| `name` | the protein's name, from the sequence file |
| `physics` | the physics set: the named group of settings it was simulated with (`physics.yaml` of the run) |
| `sequence` | the amino-acid sequence, one-letter code |
| `source_csv`, `source_row` | the sequence file and its row (counted from 1, header not counted) |
| `geometry` | `bulk` (a free chain) or `grafted` (attached to a surface; a brush too) |
| `temperature_K` | temperature, K |
| `ionic_concentration_M` | salt concentration (1:1 salt), mol/L; sets the screening of charges |
| `graft_terminus` | grafted only: which end is attached, `N` or `C` |
| `wall_charge_density` | grafted only: surface charge, e/nm² (negative attracts K, R and H) |
| `wall_attraction` | grafted only: the surface's attraction, 0 (purely repulsive) to 1 |
| `grafting_density` | grafted only: chains per nm² of surface (0 for a single attached chain) |
| `chains` | chains in each simulation: 1, or the brush's number (100 = 10 × 10) |
| `replicates` | manifest only: how many replicates went into the row's averages |
| `task_id`, `replicate`, `seed` | replicates.csv only: the task's number in the run, its replicate number (1, 2, …), its random seed |
| `N` | number of residues of the protein |
| `frames` | replicates.csv only: frames the task saved (F) |
| `run_time_h` | replicates.csv only: the task's computing time, hours |

---

## 4. Sequence descriptors

Numbers computed from the sequence alone, before any simulation, with the model's own
charges q (K, R +1; D, E −1; H +0.5; all others 0) and hydropathies λ (0 = most
hydrophilic, 1 = most hydrophobic).

| column | meaning | definition |
|---|---|---|
| `N` | length | number of residues |
| `net_charge` | net charge | Σᵢ qᵢ |
| `ncpr` | net charge per residue | Σᵢ qᵢ / N |
| `fcr` | fraction of charged residues | (number of residues with qᵢ ≠ 0) / N; histidine counts |
| `f_plus` | positive charge per residue | Σ of the positive qᵢ / N |
| `f_minus` | negative charge per residue | Σ of \|qᵢ\| over the negative qᵢ / N |
| `mean_hydro` | mean hydropathy | Σᵢ λᵢ / N |
| `scd` | sequence charge decoration (Sawle & Ghosh 2015): how strongly the charges' pattern favours compaction (negative) or expansion (positive) | Σ over i < j of qᵢ qⱼ √(j − i), divided by N |
| `shd` | sequence hydropathy decoration (Zheng et al. 2020): how hydrophobic residues close along the chain are | Σ over i < j of (λᵢ + λⱼ) / (j − i), divided by N |
| `kappa` | charge patterning κ (Das & Pappu 2013): 0 = positive and negative charges well mixed, 1 = in separate blocks | below |

**κ.** Take a window of g consecutive residues. In it, let f₊ and f₋ be the positive
charge and the |negative charge| per residue, and σ = (f₊ − f₋)² / (f₊ + f₋) (σ = 0 for a
window without charges). σ_seq is the same quantity for the whole sequence. δ_g is the
mean, over all N − g + 1 windows, of (σ_window − σ_seq)². Then

    κ = (δ₅ + δ₆) / max(δ₅ + δ₆)

where the maximum is taken over the most segregated arrangements of the same residues:
all positive residues in one block and all negative in another, in either order, the
histidines at either end of the positive block, and the neutral residues split between
the two ends of the sequence or placed between the blocks. κ is empty for a sequence
without charges. For about 2% of sequences, mostly with few charges (`fcr` < 0.2), an
arrangement other than blocks is even more uneven and κ comes out slightly above 1;
there it says little.

---

## 5. Size and shape of a chain

Every task. In a brush, these are measured on each of its chains and averaged over them.

**The gyration tensor** of a chain in one frame:

    G_ab = Σᵢ wᵢ (rᵢ,a − Rcm,a)(rᵢ,b − Rcm,b)      a, b ∈ {x, y, z}

Its eigenvalues λ₁ ≤ λ₂ ≤ λ₃ are the squared extents of the chain along its three
principal axes, and Rg² = λ₁ + λ₂ + λ₃.

| column | meaning | definition | values |
|---|---|---|---|
| `rg_nm` | radius of gyration: the chain's overall size | ⟨Rg⟩, Rg = √(λ₁ + λ₂ + λ₃) (mass-weighted) | nm |
| `rg_sq_nm2` | mean square radius of gyration | ⟨Rg²⟩ | nm²; ≥ `rg_nm`² |
| `rg_sd_nm` | how much the size fluctuates | standard deviation of Rg over all frames (and chains) | nm; no error bar of its own (section 8) |
| `rh_nm` | hydrodynamic radius (Kirkwood approximation): the size that sets how fast the chain diffuses | 1 / ⟨ (1/N²) Σ over i ≠ j of 1/rᵢⱼ ⟩ (each pair counted twice, both orders) | nm; an approximation, compare trends rather than absolute values with experiments |
| `ree_nm` | end-to-end distance | ⟨\|**R**ee\|⟩ | nm |
| `ree_sq_nm2` | mean square end-to-end distance | ⟨\|**R**ee\|²⟩ | nm² |
| `asphericity` | how far from a sphere | ⟨(λ₃ − (λ₁ + λ₂)/2) / Rg²⟩ | 0 sphere, 0.25 flat disc, 1 rod |
| `shape_anisotropy` | relative shape anisotropy κ² (not the sequence κ) | ⟨(3/2) Σₖ (λₖ − Rg²/3)² / Rg⁴⟩ | 0 sphere (or any shape as symmetric as a cube), 0.25 flat disc, 1 rod |
| `prolateness` | flattened or elongated (Dima & Thirumalai 2004) | ⟨27 Πₖ (λₖ − Rg²/3) / Rg⁶⟩ | −0.25 flat disc (oblate) … 0 … 2 rod (prolate) |
| `transverse_rms_nm` | how far the chain strays from the straight line through its ends | ⟨√( (1/N) Σᵢ dᵢ² )⟩, dᵢ = distance of **r**ᵢ from the line through **r**₁ along **R**ee | nm; 0 for a straight chain |
| `nu` | scaling exponent ν: how the distance between two residues grows with their separation along the chain | fit below | 1/3 collapsed globule, 0.5 ideal chain, 0.588 swollen chain in good solvent |
| `r0_nm` | the prefactor of the same fit | fit below | nm |
| `bond_cos` | local stiffness: the mean cosine between consecutive bonds | ⟨(1/(N − 2)) Σᵢ **b̂**ᵢ · **b̂**ᵢ₊₁⟩ | −1 … 1; 0 freely jointed |
| `persistence_nm` | persistence length (Flory's definition): how far the chain keeps its initial direction | ⟨½ (**R**ee · **b̂**₁ + **R**ee · **b̂**N₋₁)⟩: the end-to-end vector projected on the first bond's and on the last bond's direction | nm |
| `contacts_per_residue` | compactness: contact partners per residue within the chain | ⟨(1/N) × number of ordered pairs (i, j) in contact with \|i − j\| ≥ 3⟩ | 0 for an extended chain; larger when compact |

**The ν fit.** For every separation s = |i − j| along the chain, R²(s) is the mean of
rᵢⱼ² over all pairs with that separation, all frames and (brush) all chains. Then a
straight line is fitted, by unweighted least squares, to

    log √R²(s) = log r₀ + ν log s        for every s from 5 to N − 1

giving `nu` and `r0_nm`. Short separations are left out because there the chain is not
yet self-similar. Fewer than 3 points (N < 8) leaves both empty. Their error bars come
from repeating the fit on each of 10 consecutive blocks of frames (section 8).

---

## 6. Chains attached to a surface

Grafted tasks, single chains and brushes alike (in a brush: per chain, averaged over
all chains). Empty for free chains. Heights hᵢ are measured from the surface (section 2);
the anchor rests on the surface, at about 0.6 nm.

| column | meaning | definition | values |
|---|---|---|---|
| `com_height_nm` | how high the chain stands | ⟨Zcm − z_wall⟩, Zcm the z of the centre of mass | nm |
| `end_height_nm` | height of the free end | ⟨h_free end⟩ | nm |
| `max_height_nm` | height of the chain's highest residue | ⟨maxᵢ hᵢ⟩ | nm |
| `height_95_nm` | height below which 95% of the chain's residues are | ⟨95th percentile of h₁ … hN⟩ (linear interpolation between residues) | nm |
| `rg_perp_nm` | extent perpendicular to the surface | √⟨G_zz⟩ | nm |
| `rg_par_nm` | extent parallel to the surface | √⟨G_xx + G_yy⟩ | nm; `rg_perp_nm`² + `rg_par_nm`² = `rg_sq_nm2` |
| `end_lateral_rms_nm` | how far sideways the free end wanders from the anchor | √⟨(x_end − x_anchor)² + (y_end − y_anchor)²⟩ | nm |
| `tilt_cos` | how upright the chain is | ⟨cos θ⟩, θ the angle between anchor → free end and the surface normal (+z) | 1 straight up, 0 lying flat |
| `bond_order` | how the bonds are oriented | ⟨(1/(N − 1)) Σᵢ P₂(b̂ᵢ,z)⟩, P₂(c) = (3c² − 1)/2 | 1 all vertical, 0 random, −0.5 all flat |
| `adsorbed_fraction` | how much of the chain touches the surface | ⟨(1/(N − 1)) × number of residues other than the anchor with hᵢ < dᵢ + 0.5 nm⟩ | 0 … 1 |

dᵢ in `adsorbed_fraction` is residue i's **contact distance**: the height where the
surface's repulsion ends and where a residue resting on the surface sits. For a surface
of diameter σ_w (`wall_sigma`), with σ = (σᵢ + σ_w)/2, it is 3^(1/6) σ ≈ 1.20 σ for the
default 9-3 surface (`wall_potential: lj93`) and 2^(1/6) σ ≈ 1.12 σ for a 12-6 surface
(`lj126`): 0.53 to 0.71 nm, depending on the residue.

---

## 7. Brushes

Brush tasks only (more than one chain). These describe the whole brush in each frame
(all C chains of N residues together, in a box Lx × Ly wide), averaged over frames.

| column | meaning | definition | values |
|---|---|---|---|
| `brush_height_nm` | thickness of the brush | ⟨2 × mean of hᵢ over all C·N residues⟩ (a brush of uniform density up to a height H gives exactly H) | nm |
| `concentration_mg_ml` | protein concentration inside the brush | ⟨M / (Lx Ly h₉₅)⟩ × 1.66054, h₉₅ = 95th percentile of all C·N heights, M = total mass (Da) of the residues with hᵢ < h₉₅ | mg/mL (1 Da/nm³ = 1.66054 mg/mL) |
| `interchain_contacts_per_residue` | how much neighbouring chains touch | ⟨2 × (number of pairs of residues on different chains in contact) / (C·N)⟩; distances through the periodic sides in x and y | 0 when chains stay apart |
| `lateral_heterogeneity` | whether the residues spread evenly over the surface or bunch together | the surface is divided into a √C × √C grid of square cells, one grafting spacing wide (one cell per chain); ⟨standard deviation / mean of the residue counts of the cells⟩, all heights included | 0 perfectly even; larger when chains bundle |
| `reduced_density` | manifest only: how crowded the brush is compared with the size of a free chain | grafting density × π × `rg_sq_nm2` of the same protein as a single chain, at the same temperature and salt, in the same run (a bulk set if there is one, else a single attached chain); empty if the run has no such single chain | < 1 chains apart ("mushrooms"); well above 1 a stretched brush |

The size and shape columns (section 5) and the surface columns (section 6) of a brush
are averages over its chains.

---

## 8. Error bars

**In `analysis/replicates.csv`**: each observable's `_err` is one standard error of its
average over the task's frames, by **block averaging**. Successive frames are correlated,
so they can't be treated as independent; blocks much longer than the correlation time can.

1. The per-frame values are averaged over the chains (a brush) to give one series of F
   values. Its first 10 × ⌊F/10⌋ values are cut into 10 consecutive blocks of equal length.
2. `X_err` = (standard deviation of the 10 block means, with n − 1) / √10.

Columns defined through a square root or an inverse of an average take the error of
that average, propagated to first order: `rh_nm` (error of ⟨1/…⟩ × Rh²), `rg_perp_nm`,
`rg_par_nm`, `end_lateral_rms_nm` (error of the mean square / (2 × the value)). For `nu`
and `r0_nm` the fit is repeated on the R²(s) of each block, and the error is the standard
deviation of the 10 fits / √10. `rg_sd_nm` has no error bar in this table.

**In `manifest.csv`**: see section 10.

---

## 9. Quality checks

They say whether a simulation ran long enough for its averages to be trusted.

| column | meaning | definition |
|---|---|---|
| `tau_rg_tu` | correlation time of the chain's size: how long the chain takes to forget its size | integrated autocorrelation time of Rg, below, in tu |
| `tau_ree_tu` | correlation time of the chain's orientation and stretch | the same for the end-to-end vector **R**ee, in tu |
| `independent_samples` | how many independent measurements of the size the task made, per chain | F / (2 τ_Rg), τ_Rg in frames. **Below 20, the averages are not trustworthy**: run longer, or add replicates |
| `drift` | whether the chain was still changing during production | 1 if the mean Rg (and, for a chain attached to a surface, the mean centre-of-mass height) of the first half of the frames differs from that of the second half by more than 3 × √(e₁² + e₂²), with e₁, e₂ the halves' block errors (section 8); else 0. A 1 means it needed a longer equilibration |

**Autocorrelation time.** For each chain, x(t) is the series of Rg over the frames
t = 0 … F − 1, minus its mean. Its autocorrelation at a lag of k frames is
C(k) = mean over t of x(t) x(t + k) (the F − k available products), averaged over the
chains and divided by C(0). For a vector (**R**ee) the product is the dot product. Then

    τ = ( ½ + Σ C(k) for k = 1 … k* − 1 ) × Δt

where k* is the first lag at which C drops below 0 (beyond it the sum is mostly noise),
and Δt is the time between two frames in tu: (production steps / frames) × the time step.

---

## 10. From replicates to the manifest

Each row of `manifest.csv` combines the n replicates of one protein under one physics
set (`replicates` = n):

- **Each observable** is the mean of the replicates' values.
- **Its error bar** `X_err` is the standard error of that mean: (standard deviation of the
  n values, with n − 1) / √n, when n ≥ 2. Replicates are fully independent simulations,
  so this is the most reliable error bar. With n = 1 it is the replicate's own block error
  (section 8). Replicates where a value is empty are left out of its mean.
- `tau_rg_tu`, `tau_ree_tu`: the mean over replicates.
- `independent_samples`: the **sum** over replicates (the independent measurements behind
  the average).
- `drift`: how many of the replicates drifted.
- `reduced_density`: computed from the row's own grafting density and the manifest row of
  the matching single chain (section 7).

---

## 11. Computing other quantities from the trajectories

Curves such as a density profile, residue-by-residue averages, or any new observable
are not computed on the cluster: for thousands of tasks they would be far larger than
the tables. Instead, compute them afterwards from the trajectories, on any computer.
You need the run's `tasks.csv`, its `physics.yaml` and the `.gsd` files from
`trajectories/` (they stay there until `noodle clean`; copy them to storage first).

**What a trajectory holds.** `trajectories/task_NNNN_<name>.gsd` is one task's production
frames, readable with the `gsd` Python package (installed with noodle-md; MDAnalysis
also reads it). Each frame has:

- `particles.position`: positions, nm, wrapped into the periodic box. Chain k (counting
  from 0) is particles k·N … k·N + N − 1, residue 1 first;
- `particles.types` and `particles.typeid`: each particle's amino acid (one-letter code);
- `bonds`: the bonded pairs;
- `configuration.box`: the box, Lx, Ly, Lz in nm, centred on 0;
- `configuration.step`: the time step. Two frames are (production steps / frames)
  steps apart; times in tu are steps × `dt`.

For a chain attached to a surface the **surface** is the plane z_wall = −Lz/2 + 0.5 nm,
and the box is not periodic in z for the chains. Masses, charges and diameters of each
amino acid are in `src/forcefield.py` (`AA_PARAMS`).

**Example: the density profile of a brush** (residues per nm³ against height)

```python
import gsd.hoomd
import numpy as np

with gsd.hoomd.open("runs/my_run/trajectories/task_0003_P12345.gsd") as f:
    Lx, Ly, Lz = f[0].configuration.box[:3]
    z_wall = -Lz / 2 + 0.5                                   # the surface
    heights = np.concatenate([frame.particles.position[:, 2] - z_wall for frame in f])
    frames = len(f)

bin_width = 0.25                                             # nm
counts, edges = np.histogram(heights, np.arange(0, heights.max() + bin_width, bin_width))
density = counts / (frames * Lx * Ly * bin_width)            # residues per nm^3
```

**Example: whole chains, then the distances between residues against their separation**

```python
import gsd.hoomd
import numpy as np

def whole_chains(frame, n):
    """[chains, n, 3] positions of chains of n residues, each made whole."""
    box = np.asarray(frame.configuration.box[:3], dtype=float)
    x = frame.particles.position.astype(float).reshape(-1, n, 3)
    bonds = np.diff(x, axis=1)
    bonds -= box * np.round(bonds / box)                     # each bond's shortest vector
    return np.concatenate([x[:, :1], x[:, :1] + np.cumsum(bonds, axis=1)], axis=1)

n = 120                                                      # the task's N (tasks.csv, column length)
with gsd.hoomd.open("runs/my_run/trajectories/task_0001_P12345.gsd") as f:
    x = np.stack([whole_chains(frame, n) for frame in f])    # [frame, chain, residue, xyz]
r2 = np.sum((x[:, :, :, None] - x[:, :, None, :]) ** 2, axis=-1).mean(axis=(0, 1))
distance = [np.sqrt(np.diagonal(r2, s).mean()) for s in range(1, n)]   # √R²(s), s = 1 … N − 1
```

(For long chains or many frames, loop over the frames instead of stacking them all, to
keep memory small.) The functions in `src/observables.py` do the same steps for every
column above and can serve as a starting point.

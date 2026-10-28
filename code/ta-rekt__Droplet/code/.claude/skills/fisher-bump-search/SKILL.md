---
name: fisher-bump-search
description: >-
  Parameter-search playbook for the surface-phase-transition Fisher-rate study.
  Use when locating the single-sphere prewetting line (solve_isolated_sphere),
  tuning the two-sphere BVP (spec/fisher_two_spheres.jl) to straddle it, or
  building cluster/parallel search scripts to find the bump in the Fisher
  information rate. Governs any work under the fisher_two_spheres /
  fisher_rate_from_fpt_solver pipeline.
---

# Fisher-bump parameter search

## 1. Scientific objective

Show that a **surface phase transition** between two beads of a Rouse polymer —
as they approach and their adsorption layers merge — puts a **discontinuity in
the pair potential `U(d)`**. That discontinuity makes the bead binding/unbinding
process **hypersensitive to the background concentration `phi_inf`**, which shows
up as a **bump (peak) in the Fisher information rate** of the binding/unbinding
dwell times with respect to `phi_inf`.

We compute the Fisher rate two ways and contrast them:

- **Linear-response (closed form):** `make_U_tot(p, lam)` in
  [viz_util/fisher_rate_from_fpt_solver.jl](https://github.com/ta-rekt/Droplet/blob/HEAD/viz_util/fisher_rate_from_fpt_solver.jl)
  — an exponential adsorption well, valid only in the linear (small-perturbation)
  regime. Cannot see the transition (it is a nonlinear/first-order effect).
- **Full nonlinear BVP:** [spec/fisher_two_spheres.jl](https://github.com/ta-rekt/Droplet/blob/HEAD/spec/fisher_two_spheres.jl)
  solves the axisymmetric Cahn–Hilliard / Flory–Huggins two-sphere BVP, captures
  the capillary-bridging / surface-transition branch, and yields `U_full(d)`.
  Its dashboard [viz_util/dashboard_fisher_two_spheres.py](https://github.com/ta-rekt/Droplet/blob/HEAD/viz_util/dashboard_fisher_two_spheres.py)
  emits a JSON (`U_full_grid`, `Phi_surf_grid`, `phi_inf`, `d_grid`, …) that is
  the **alternative input** to the fpt-solver script (in place of `make_U_tot`).

The whole point of the search: **tune the nonlinear BVP so that, at small radius
`R = 5` and fixed `chi_b = 4`, sweeping `phi_inf` drives the two beads across the
surface transition** — producing the discontinuity that the linear model misses.

**"Optimal" is loose — this is positioning, not maximization.** We do NOT do a
systematic search for the largest Fisher-rate bump. The target is simply to
choose parameters (`omega_b`, and the `phi_inf`- and `d`-sweep ranges) such that:
- the two spheres **cross** the surface transition at a gap `d` **inside** the
  `d`-sweep range, **ideally near the middle** of that range; and
- the crossing happens at a `phi_inf` **near the middle** of the `phi_inf`-sweep
  range.
Centering the crossing in both windows is what makes the discontinuity (and hence
the Fisher-rate feature) cleanly visible; that is the whole objective.

## 2. Physics of the knobs (Flory–Huggins symmetric bulk, `nb = 1`)

| symbol | meaning | role in the search |
|---|---|---|
| `phi_inf` | bulk volume fraction (the *signal* being sensed) | swept axis; higher → closer to saturation → stronger adsorption |
| `omega_b` | surface adsorption field (`< 0` attractive) | tuning knob; more negative → stronger surface attraction |
| `chi_b` | Flory interaction (fixed `= 4`, `> 2` so bulk coexistence/bridging exists) | fixed |
| `R` | sphere radius (fixed `= 5`, small) | fixed; smaller `R` ⇒ more curvature ⇒ transition needs **more negative** `omega_b` than at `R=10` |
| `kappa` | gradient (surface-tension) coefficient, `= 1` | fixed |
| `d` | surface-to-surface gap between the two spheres | integrated over via `P_eq(d)` |
| `bc_reg` | Robin-BC regularization | pins a unique finite adsorption branch near the transition (see §5) |

Dilute binodal for `chi_b=4, nb=1` is `phi_b ≈ 0.0179` (root of
`log(phi/(1-phi)) + chi*(1-2phi) = 0`). The transition lives **below** the
binodal, at small `phi_inf`.

## 3. Anchor data (known result — must reproduce first)

For a **single, isolated sphere** at `chi = 4`, `R = 10`, `nb = 1`, the surface
phase transition crosses these `(omega_b, phi_inf)` points:

| omega_b | phi_inf |
|---|---|
| -1.45 | 0.008 |
| -1.41 | 0.009 |
| -1.37 | 0.010 |

Trend: as `phi_inf` ↑, the transition `omega_b` ↑ (less negative) — nearer bulk
saturation needs less surface field. **We need the analogous line for `R = 3` and
`R = 5`,** which sit at **more negative `omega_b`** (more curvature) for the same
`phi_inf`.

**Confirmed (2026-07-06):** the §5 `phi_s`-inflection method on
`solve_isolated_sphere` (chi=4, R=10, bc_reg=0.02, N=1000-1500) reproduces the
known line's **trend and magnitude to within ~0.01-0.05 in `omega_b`** — good
enough. It is *not* an exact match (fixed-`bc_reg` inflection line is ~1.6× steeper
in slope than the very-approximate anchor points, and no single `bc_reg` hits all
three). The anchor values are approximate and the Julia solver is a minimal port
of the Python `bvp_solver.py` pipeline, so treat the 1D line as a **seed only**
(§6c) — don't tune `bc_reg`/resolution to match anchors exactly.

## 4. Tools & pipeline

**Single sphere (cheap, 1D — use to map the transition line):**
`solve_isolated_sphere(p; L_iso, N, f_tol, maxiter, bc_reg, omega_continuation)`
in [src/bvp_solver.jl](https://github.com/ta-rekt/Droplet/blob/HEAD/src/bvp_solver.jl). Returns an `IsolatedSphere1D`
with `.z`, `.phi`, `.dphi`, `.success`; **surface volume fraction = `phi[1]`**
(z = 0 is the sphere surface). It traces ONE (thin/continuation) branch — it
takes no custom `x0`, only `omega_continuation` from the uniform state. `p =
Parameters(nb, chi_b, omega_b, kappa, phi_inf, R, eps, …)`.

**Two spheres (expensive, 2D — the actual target):**
`solve_two_spheres(p, d; rho_max, z_buf, Nrho, Nz, f_tol, maxiter, x0_init,
omega_continuation, bc_reg, verbose)` and `excess_free_energy_two_spheres`,
`phi_on_sphere`, `build_geometry` in [src/bvp_two_spheres.jl](https://github.com/ta-rekt/Droplet/blob/HEAD/src/bvp_two_spheres.jl).

**Full pipeline:**
```
spec/fisher_two_spheres.jl   (BVP sweep, parallel over phi_inf via runJL/sweep.jl)
   -> data/<date>/fisher_two_spheres_.../output_cpu*.json
dashboard_fisher_two_spheres.py   (assembles U_full_grid, Phi_surf_grid, Fisher rate)
   -> media/figures/<date>/<exp>.json  +  .png
viz_util/fisher_rate_from_fpt_solver.jl   (reads that JSON as U(d) input to the FPT solver)
```

## 5. Detecting the transition  ✅ DECIDED: phi_s inflection, bc_reg > 0

The isolated solver has no free-energy function, so the transition is read off
the **surface volume fraction `phi_s = phi[1]`** as a control parameter is swept
on the **unique regularized branch (`bc_reg > 0`)**:

> **Transition = the control value where `phi_s` rises most steeply, i.e. the
> point of maximal `|d phi_s / d(control)|` (inflection of the sigmoidal
> `phi_s`-vs-control curve).**

Practically: sweep `omega_b` (fine grid) at fixed `phi_inf`, warm-started /
`omega_continuation`, record `phi_s[k] = s.phi[1]`, then locate the peak of the
finite-difference derivative (parabolic-refine the discrete argmax for sub-grid
accuracy). Do the same the other way (sweep `phi_inf` at fixed `omega_b`) as a
cross-check — both give the same transition line.

`bc_reg`: the Robin BC is `-kappa*phi'(0) + omega_b + bc_reg*log(phi_s/(1-phi_s)) = 0`.
`bc_reg > 0` pins a finite adsorption branch (the log term diverges as `phi_s→1`,
opposing runaway), turning the true first-order jump into a **steep but
continuous** sigmoid whose inflection marks the transition. Smaller `bc_reg` ⇒
sharper inflection (closer to the true discontinuity) but harder to converge;
`bc_reg ≈ 0.02` is the working value. Report the inflection at a fixed `bc_reg`
so the line is defined consistently across the sweep.

## 6. Time-efficient search strategy

1. **Map the line cheaply on the 1D solver first.** `solve_isolated_sphere` is
   ~ms–s; `solve_two_spheres` is ~0.6–2 s/solve and a full `d`-march is tens of
   seconds. Bracket the transition in `(omega_b, phi_inf)` on the single sphere
   before spending any two-sphere budget.
2. **Continuation / warm starts everywhere.** Sweep `omega_b` (or `phi_inf`)
   monotonically and warm-start each solve from the previous solution. For two
   spheres, march `d` inward (dilute branch) and outward from contact (bridged
   branch), reusing `_warm_x0` / `_dense_seed` (already in the spec).
3. **Coarse → fine bracketing.** Locate the transition to a coarse grid, then
   bisect only the bracketing interval. Do not run a fine global grid.
4. **Center the crossing, don't maximize.** Objective is positioning (§1): pick
   `omega_b` and the `phi_inf`/`d` windows so the two-sphere transition crossing
   lands **near the middle** of both sweeps. Procedure:
   a. On the 1D single sphere, find the `R=5` transition line `omega_b*(phi_inf)`
      (§5). Pick a `phi_inf_mid` for the centre of the intended `phi_inf` sweep,
      read off `omega_b* = omega_b*(phi_inf_mid)`, and **fix `omega_b` there**.
   b. With that `omega_b`, the two-sphere pair crosses the transition as `d`
      shrinks (the two adsorption layers merge). Choose the `d`-sweep range so the
      crossing sits mid-range; choose the `phi_inf`-sweep range so `phi_inf_mid`
      is its centre and both pre- and post-transition `phi_inf` are covered.
   c. The single-sphere line is only a *seed*: the two-sphere crossing shifts
      (bridging lowers the barrier), so confirm on `solve_two_spheres` and nudge
      `omega_b` / the ranges until the crossing is centred. This is a small manual
      loop, not a grid optimization.
5. **Respect the small-d wall (physics, not resolution).** In the undersaturated
   regime the inward march stalls near `d/R ≈ 0.14`; increasing `Nrho/Nz` does
   NOT reach smaller `d`. Only the bridged (dense-seeded, near-binodal) branch
   reaches contact. Don't waste resolution chasing sub-wall `d`.
6. **Resolution budget (from prior study, P_strong2, R=10):** `Nrho ≈ 240–351`,
   `dz_target ≈ 0.04–0.05`, `rho_max ≈ R+4`, `z_buf = 3`, `bc_reg ≈ 0.02`,
   warm-started. Runtime ∝ n_unknowns^~1.1; the 10-min per-target budget has
   ~10× headroom. Re-tune for `R=5` (fewer unknowns → cheaper).

## 7. Cluster execution  ✅ DECIDED: SLURM array, one task per omega_b

Structure: a **SLURM job array where each task fixes one `omega_b`** and runs the
full `phi_inf × d` two-sphere BVP sweep for that `omega_b`, writing its own JSON
shard. This matches the "fix `omega_b`, sweep `phi_inf`" objective (§6) and
load-balances naturally. (Contrast: today only local pipe-parallelism exists —
`runJL` pipes N `julia … sim/sweep.jl <spec> <cpukey> <ncpus> <outdir>`, each CPU
taking a shard of `alphas`.)

**Budget = coarse scout → refine:** Stage 1 runs a coarse `(omega_b × phi_inf)`
grid with a coarse `d`-grid to bracket where the crossing centres; Stage 2 runs a
fine grid only in the bracketing window. Minimize total solves.

⚠️ **Still needed — concrete SLURM env** (ask before writing the sbatch): partition/
queue name, cores per node & mem, walltime cap, how Julia is loaded (`module load`
vs an absolute path / juliaup), account/QOS, and output directory convention on
the cluster filesystem. Write the sbatch parameterized with these at the top.

## 8. Resolved decisions (this study)

- **Detection (§5):** `phi_s = phi[1]` inflection (max `|d phi_s/d control|`) on
  the unique `bc_reg ≈ 0.02` branch. No hysteresis / free-energy needed.
- **Objective (§1, §6):** *position*, don't maximize. Fix `omega_b` on the R=5
  transition line at the chosen `phi_inf_mid`; set `phi_inf`/`d` sweep ranges so
  the two-sphere crossing is centred in both. Confirm & nudge on `solve_two_spheres`.
- **Free variable:** `omega_b` (the tuning knob); `phi_inf` is the swept signal.
- **Cluster (§7):** SLURM array, one task per `omega_b`; coarse-scout → refine.
- **Remaining input needed:** concrete SLURM env (see §7).

## 9. Pitfalls (from prior work / memory)

- Small-d wall (`d/R ≈ 0.14`) is physics — don't fight it with resolution.
- Bridged branch only exists near saturation; dense seeds run away to `phi=1`
  when undersaturated.
- FPT solver: dense eigendecomposition loses the slow mode if the absorbing
  point sits ≳ 20 kBT above the well (MFPT overflows) — keep the barrier sane.
- Grid mismatch: the fpt script tabulates `U` on `[R, x_ridge]`; a BVP-JSON `U(d)`
  is on `d_grid` (surface gap). Convert `d → x = d + 2R` (or `d + R` per the
  script's convention) consistently when feeding `U_full` in.

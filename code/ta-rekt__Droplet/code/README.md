# Droplet

Lattice Monte Carlo simulations of polymers and biomolecular condensates. The repo also
holds a continuum pipeline for condensate-mediated enhancer–promoter sensing:
single-sphere prewetting → two-sphere BVP → 1D first-passage solver → Rouse-chain
simulation → Fisher information. The Julia core is in `src/`. Experiments are spec files
under `spec/`, and figures are made in the three `viz*.ipynb` notebooks. `CLAUDE.md` and
`context/*.md` document the lattice-MC conventions.

## Setup

`Manifest.toml` was resolved with **Julia 1.11.1**, so use exactly that version. Install
it with [juliaup](https://github.com/JuliaLang/juliaup):

```bash
curl -fsSL https://install.julialang.org | sh    # installs juliaup; open a new shell afterwards
juliaup add 1.11.1
juliaup default 1.11.1                           # runJL calls plain `julia`, so the default must be 1.11.1
julia --version                                  # -> julia version 1.11.1
```

Pin `1.11.1` rather than relying on the `release` channel, because `juliaup update` moves
`release` to the newest version. To run a single command on 1.11.1 without changing the
default, use `julia +1.11.1 ...`.

Then install the Julia and Python dependencies:

```bash
julia --project=@. -e 'using Pkg; Pkg.instantiate()'
pip install -r requirements.txt jupyter        # numpy 2.2.6, scipy, matplotlib, numdifftools, pillow
```

- **Run everything from the repo root.** `--project=@.` finds `Project.toml` from the
  working directory, and every `data/...` path in specs, `param.jl` files and notebooks
  is relative to it.
- The figures use `text.usetex = True`, so a LaTeX install is needed.
- NumPy ≥ 2.4 removed `np.trapz`. In `viz-rouse.ipynb` only cell 12 and the legacy cells
  call it. Swap in `np.trapezoid` or use the pinned 2.2.6.
- `data/` and `spec/jobs/` (cluster copies of specs) are gitignored. On a fresh clone,
  regenerate the data as described below. `media/figures/` is tracked and holds the
  reference PDFs.

## Notebooks

| Notebook | Contents |
|---|---|
| `viz.ipynb` | Lattice-MC experiments from Nov 2025 to Aug 2026: 3D polymer and density maps, filament and enhancer prewetting, polymer and heteropolymer phases, looping, transcription and enhancer search, binding curves, Gillespie mRNA accumulation and LLR, polymer MSD, and Flory–Huggins phase diagrams. Reads `data/<date>/<exp>/output_cpu*.json` from `spec/*.jl` sweeps. |
| `viz-mc.ipynb` | MC proximal-promoter transcription analysis (`data/20260417`, `20260418`, `20260823`): ON/OFF times with hysteresis, conditional local concentration, Fisher rates with bootstrap errors, decision time. |
| `viz-rouse.ipynb` | The continuum/Rouse pipeline: two-sphere pair potentials, FPT PDFs, Rouse vs 1D Fisher information, Fisher rate vs concentration and vs genomic separation, prewetting constructions, and single-sphere Fisher rates. See [Reproducing viz-rouse.ipynb](#reproducing-viz-rouseipynb). |

## Running scripts

### Parallel sweeps: `runJL` + `sim/sweep.jl`

```bash
./runJL sim/sweep.jl <specfile> <ncpus> <outdir>
```

- `runJL` starts `ncpus` copies of `julia --project=@. sim/sweep.jl <specfile> <k> <ncpus> <outdir>`
  joined by pipes, so they all run at once. Only the last process's stdout reaches the
  terminal, but stderr from every process does. **Specs must not print to stdout.**
  Once a downstream process exits, an upstream one that prints dies with `EPIPE` and
  writes no output. Use `@info`/`@warn`, which go to stderr.
- `<specfile>` must be relative to the repo root, because `sweep.jl` does `include("../" * specfile)`
  from `sim/`. An absolute path breaks.
- The spec must define `slist, alphas, betas, nparams, blurb, checkpoint_path, target!`.
  CPU `k` runs `keya = (k-1)n+1 … kn` with `n = div(nparams, ncpus)`, and each of those
  runs every `beta`. **The remainder is silently dropped**, so pick an `ncpus` that
  divides `nparams`.
- Output is `<outdir>/output_cpu<k>.json`, shaped `{keya: {keyb: record}}`.
- Several `runJL` lines can be piped together to run concurrently; see `run.sh`.

### `spec/local/<date>/`: snapshots of what ran that day

- **Lattice-MC specs** (≤ 2026-03) are dated copies of `spec/*.jl` files, run as above.
  Their contract is in `context/spec_template_guide.md`.
- **Job folders** (≥ 2026-07-20), e.g. `spec/local/20260912/test3/`, contain a
  `param.jl` and verbatim copies of the stage scripts. `param.jl` defines `PARAMS`,
  the single source of truth. Every script does `include(joinpath(@__DIR__, "param.jl"))`,
  then reads its inputs from and writes its outputs to `PARAMS.sim.output_dir`.
  To change a run, edit `param.jl`, never the scripts. `spec/local/make_grid.jl`
  generates one folder per parameter-grid cell.

| Stage script (in a job folder) | How to run | Reads → writes (in `PARAMS.sim.output_dir`) |
|---|---|---|
| `demo_bvp.jl` (single sphere) | `julia --project=@. <folder>/demo_bvp.jl` | → `demo_bvp_profiles.json` |
| `fisher_two_spheres.jl` (two-sphere BVP, sharded over gap `d`) | `./runJL sim/sweep.jl <folder>/fisher_two_spheres.jl <n_d+1> <output_dir>` | → `output_cpu*.json`, `bridged_seed_*.json` |
| `bvp_potential_interp.jl` | `julia --project=@. <folder>/bvp_potential_interp.jl` | `output_cpu*.json` → `bvp_potential_interp.json` |
| `fisher_rate_from_fpt_solver.jl` | `julia --project=@. <folder>/fisher_rate_from_fpt_solver.jl` | `bvp_potential_interp.json` → `fisher_rate_from_fpt_solver.json` |
| `rouse_mfpt_attractive.jl` | `./runJL sim/sweep.jl <folder>/rouse_mfpt_attractive.jl <ncpus> <other_dir>` | both JSONs above → `<other_dir>/output_cpu*.json` |
| `sweep_rouse_s.jl` | `julia --project=@. <folder>/sweep_rouse_s.jl` | reruns interp + FPT for each `s` → `fisher_vs_rouse_s.json` |

For the sweep specs, `sim/sweep.jl` stamps `spec_path`, `output_dir` and `ncpus` into
`PARAMS.sim` *after* the include, so every record carries its own provenance under
`record["params"]["sim"]`.

### Example: `rouse_mfpt_attractive.jl`

```bash
./runJL sim/sweep.jl spec/local/20260912/test3/rouse_mfpt_attractive.jl 30 data/20260912/rouse
```

- **Inputs** are read at include time from the `output_dir` in `param.jl`
  (`data/20260912/test3`). Run `bvp_potential_interp.jl` and
  `fisher_rate_from_fpt_solver.jl` first.
- **The outdir must differ** from that `output_dir`. Otherwise it overwrites the BVP
  `output_cpu*.json`.
- **`alphas`** are the solver's φ∞ values (`n_phi − 2 = 30`), thinned to `n_phi_run`.
  The script warns if `ncpus` does not divide `nparams`.
- **Run-length knobs** live in the spec, not in `param.jl`: `nsteps`, `n_phi_run`,
  `save_mode` and `seed0`. Size `nsteps` against `tau_rouse`.
- **Seeds** are `seed0 + keya`, so the RNG stream does not depend on `ncpus`. Output is
  still not bit-reproducible across machines or BLAS kernels, so compare statistics,
  not bits.
- **Start small** (`nsteps = 1_000_000`, few CPUs), time it, then scale up.

## Reproducing viz-rouse.ipynb

### 1. Data

Skip any step whose output already exists. The steps marked ✓ were re-run on
2026-10-07 and matched the stored files bit for bit.

```bash
# (a) Two-sphere BVP sweep: the root of the chain. 18 gaps → 18 CPUs, 4–18 min per gap.
./runJL sim/sweep.jl spec/local/20260902/test3/fisher_two_spheres.jl 18 data/20260902/test3

# (b) Fan the BVP output out. The three downstream dirs hold byte-identical copies.
for d in data/20260912/test3 data/20260912/test3_clone data/20260927/s110; do
  mkdir -p $d && cp data/20260902/test3/output_cpu*.json $d/; done

# (c) Main chain (s = 300, a = 20, R = 130): interp + FPT take ~1.5 min.   ✓
julia --project=@. spec/local/20260912/test3/bvp_potential_interp.jl
julia --project=@. spec/local/20260912/test3/fisher_rate_from_fpt_solver.jl
julia --project=@. spec/local/20260912/test3/demo_bvp.jl
julia --project=@. spec/local/20260902/test3/demo_bvp.jl      # only the binodal is read

# (d) Rouse chain: the expensive step. First edit spec/local/20260912/test3/rouse_mfpt_attractive.jl
#     to match the reference run:
#       nsteps = 700_000_000;  n_phi_run = 30;  alphas = REC_PHI[PHI_SEL]
#     (the local copy currently has alphas = REC_PHI[[21,22,23]]).
#     Reference run: 30 CPUs, ~27 h per φ on a loaded cluster, before the 2026-09-27
#     SymTridiagonal speed-up. The notebook needs all 30 φ values.
./runJL sim/sweep.jl spec/local/20260912/test3/rouse_mfpt_attractive.jl 30 data/20260912/rouse

# (e) Fisher vs genomic separation (s = 5…1000, 64 values; budget ≲1.5 h).
#     Rewrites its param.jl in place and restores it on exit. Don't kill it with -9.
#     It leaves the last s's interp/FPT JSONs in data/20260912/test3_clone.
julia --project=@. spec/local/20260912/test3_clone/sweep_rouse_s.jl

# (f) Interp at s = 110 (used for the potential figures).
julia --project=@. spec/local/20260927/s110/bvp_potential_interp.jl

# (g) Single-sphere runs (~10–20 s each).   ✓
for r in large_R medium_R small_R; do julia --project=@. spec/local/20260825/$r/demo_bvp.jl; done
for r in small_R medium_R large_R; do julia --project=@. spec/local/20261004/$r/demo_bvp.jl; done
for r in 1 2 3 4 5;              do julia --project=@. spec/local/20260825/run$r/demo_bvp.jl; done
```

**Folder-name pitfall in step (g).** The notebook needs `data/20260825/run1…run5` to hold
h0 (`omega_b`) = −0.05, −0.04, −0.0376, −0.0365 (h\*), −0.0358, in that order. The spec
folders were renamed after the original run, so `spec/local/20260825/run{5,1,2,4,3}` hold
those values. After the loop, rotate the output folders:

```bash
(cd data/20260825 && mv run5 tmp && mv run3 run5 && mv run2 run3 && mv run1 run2 && mv tmp run1)
python -c "import json;[print(r,json.load(open(f'data/20260825/{r}/demo_bvp_profiles.json'))['params']['surface']['omega_b']) for r in ['run1','run2','run3','run4','run5']]"
```

**BVP seeding in step (a).** `bvp.seed_path` warm-starts the bridged branch from
`bridged_seed_*.json` files. The stored records seeded from `data/20260902/test3` (see
`record["seed_from"]`). From an empty directory, the sweep falls back to the analytic seed
with a warning and writes new seeds. Rerunning into the same directory then
warm-starts every gap. This path has not been re-run, so it is not verified to reproduce
the stored records bit for bit.

### 2. Notebook

Cells share globals, and several names are reassigned between sections (`v_tf`,
`binodal`, `phi_inf`, `b`, `R`, `F`), so do not run the notebook top to bottom.
Restart the kernel and run the groups below in order. Cell numbers are 0-based
indices as of 2026-10-07; section headings are in quotes. Each figure is written by the
`savefig` cell right after its plot cell. Those cells overwrite the tracked PDFs, so check
the result with `git diff --stat media/`.

**Group A: two-sphere and Rouse figures** (needs data from steps a–f)

| Run | Section | Writes `media/figures/…` |
|---|---|---|
| 1, 128 | imports; "Flory Huggins energy density" helpers (`to_nM`, `omega`, …) | — |
| 83, 87 | "Two sphere BVP (*d*-major)": loads `data/20260902/test3` records, the binodal, and the `s110` interp | — |
| 89, 90 | "Interpolated potentials" | `20260917/F_attr.pdf` |
| 91, 92 | "Interpolated potentials" | `20260912/combined_potential_vs_s.pdf` |
| 94, 95 | "Minimum *s* per max *c*" | `20260912/min_sep.pdf` |
| 97, 99, 100 | "Exact FPT PDFs from solver": sets `coef = 50`, used by every later cell | `20260916/exit_return_pdfs_langevin.pdf` |
| 105 | "Comparison with Rouse simulation": loads `data/20260912/rouse` into `data`/`alphas`; warns on an `x_div` mismatch | — |
| 102, 103 | "FPT (Rouse)": run *after* 105 | `20260916/exit_return_pdfs_rouse.pdf` |
| 110, 111, 112, 113 | "Fisher rate from Rouse simulation": asserts the Rouse φ grid matches the solver's | `20260916/fisher_per_cycle.pdf` |
| 116, 117 | "Polymer Fisher rates": sets `v_tf = 120` (nb = 4000) | `20260917/allFs_tau_rouse.pdf` |
| 122, 124 | "Fisher rate versus genomic separation": reuses 116's globals | `20260831/fisher_vs_s.pdf` |

**Group B: prewetting and single-sphere figures** (needs step g)

In a fresh kernel, run cells 1 and 128 first. Then run cell 135 ("Estimating kappa"). It
sets `v_tf = 67` (Klf4; these runs use nb = 2230), along with `v_water`, `b`, `kappa`,
`h0_lo` and `phi_sat`. If group A ran in the same kernel, cell 135 must run after it.

| Run | Section | Writes `media/figures/…` |
|---|---|---|
| 135, 138, 139 | "Prewetting graphical construction" (138 only; 137 is broken) | `20260821/pw_graphical.pdf` |
| 143, 144 | "Influence of curvature" | `20260825/pw_curvature3.pdf` |
| 146, 147 | "Single sphere from numerical solver (varying h0)" | `20260912/surface_fisher_3vals.pdf` |
| 149, 150 | "Single sphere from solver (varying R)" | `20260912/surface_fisher_3vals_vsR.pdf` |

**Group C: no data needed**

| Run | Section | Writes `media/figures/…` |
|---|---|---|
| 1, 4–13 | "Recovering the Arrhenius exit time": overdamped Langevin toy (`viz_util/overdamped_langevin.py`). Cell 12 uses `np.trapz`. | `20260606/empirical_vs_boltzmann_for_arbitrary_potential_high_T.pdf` |
| 156–159 | "Parallel plates": redefines `F`, so run it last | `20260802/parallel_plates_graphical_hi_c_low_h0_1.pdf` |

**Legacy (not reproducible from the current specs):** cells 19–73 (old
Rouse/Langevin runs in `data/20260615`, `20260627`, `20260629` and JSONs under
`media/figures/20260528` and `20260628`; several of these cells already error); 74–78
(φ-major BVP, `data/20260728`); 118–120; 137; 151–154; 160–161.

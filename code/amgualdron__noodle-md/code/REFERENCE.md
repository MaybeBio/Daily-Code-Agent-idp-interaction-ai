# noodle-md reference

How noodle-md works and why: the model, the pipeline from a CSV of sequences to tables
of results, every file it writes, and how to change it. The [README](https://github.com/amgualdron/noodle-md/blob/HEAD/README.md) is for
using the tool; this document is for understanding it, checking it, and developing it
further. The exact definition of every column of the results tables is in
[OBSERVABLES.md](https://github.com/amgualdron/noodle-md/blob/HEAD/OBSERVABLES.md).

**Contents**

1. [What it does](#1-what-it-does)
2. [Words used](#2-words-used)
3. [The model](#3-the-model)
4. [The pipeline, step by step](#4-the-pipeline-step-by-step)
5. [Files](#5-files)
6. [Sequence descriptors](#6-sequence-descriptors)
7. [Observables and statistics](#7-observables-and-statistics)
8. [How it is checked](#8-how-it-is-checked)
9. [The code, its rules, and how to extend it](#9-the-code-its-rules-and-how-to-extend-it)
10. [Limitations and open items](#10-limitations-and-open-items)

---

## 1. What it does

noodle-md runs coarse-grained molecular dynamics of intrinsically disordered proteins
(IDPs) for many sequences at once, and turns the trajectories into tables of averaged
properties. Each protein can be simulated in four situations:

| situation | what is simulated |
|---|---|
| **bulk** | one free chain in a large box of implicit salt water |
| **grafted** | one chain attached by one end (N- or C-terminus) to a flat surface |
| **charged wall** | a grafted chain over a uniformly charged surface |
| **brush** | 100 copies of the chain grafted on a square lattice (10 x 10), periodic sideways like an infinite brush |

A **run** is one study: a list of sequences, one or more physics sets, a number of
replicates. Runs go from thousands of short tasks on a SLURM cluster to a handful on
a laptop, with the same commands. The tool is meant for computing: each user works in
their own copy of the repository on the cluster, launches runs, and takes the result
tables home (`noodle export`) for plotting and further analysis.

The pipeline in one picture:

```
sequences CSV + run.yaml + physics.yaml
        │  noodle launch: check, freeze, descriptors -> manifest.csv, group tasks into jobs
        ▼
SLURM job array (one core per job) ── each job, task after task:
        │     build -> equilibrate -> produce trajectory -> analyse -> results into its job file
        ▼
jobs/job_NNNN.csv
        │  noodle analyze: pool replicates
        ▼
manifest.csv · analysis/replicates.csv · analysis/report.json
        │  noodle export
        ▼
exports/<name>/   -> your own computer
```

---

## 2. Words used

| word | meaning |
|---|---|
| **run** | a folder `runs/<name>/`: one study, with its settings, tasks and results. Each run has a number |
| **physics set** | a named group of physics settings in `physics.yaml` (temperature, salt, surface, run length...) |
| **task** | one simulation: one protein, under one physics set, one replicate |
| **replicate** | an independent copy of a task with its own random seed |
| **job** | what SLURM runs: one core working through one or more tasks, one after another |
| **job array** | many jobs submitted as one SLURM submission; each element is one job |
| **frame** | one saved configuration (positions of every residue) |
| **equilibration / production** | the first part of a task, discarded (the chain forgets its start); the second, saved as frames |
| **checkpoint** | a task's full state saved to disk, so it can continue later |
| **trajectory** | a task's frames, one `.gsd` file |
| **observable** | a number measured from a trajectory (Rg, height, ...) |
| **descriptor** | a number computed from the sequence alone (charge patterning, hydropathy...) |
| **tu** | the model's unit of time (0.489 ps; see 3.1) |

---

## 3. The model

The HPS model (hydropathy scale; Dignon et al., *PLoS Comput. Biol.* 2018), run with
HOOMD-blue on CPUs. All physics lives in `src/forcefield.py` (forces) and
`src/simulate.py` (system building, dynamics). Every parameter value is in `physics.yaml`.

### 3.1 Units

Length nm, energy kcal/mol, mass g/mol (Da), charge e. The unit of time that follows
from these is √(g/mol · nm² / (kcal/mol)) = 0.489 ps. The time step `dt` 0.02 is ~10 fs.

### 3.2 Residues

One bead per amino acid, with mass, charge q, diameter σ and hydropathy λ (0 =
hydrophilic, 1 = most hydrophobic). Charges: K and R +1, D and E −1, **H +0.5**, all
others 0. The full table is `AA_PARAMS` in `src/forcefield.py`. Only the 20 standard
letters are accepted; sequences with any other letter are skipped at launch.

### 3.3 Forces between residues

- **Short range, Ashbaugh-Hatch.** For a pair with mean diameter σ = (σi + σj)/2 and
  mean hydropathy λ = (λi + λj)/2, with U_LJ the 12-6 Lennard-Jones potential of depth ε:

  U(r) = U_LJ(r) + (1 − λ) ε  for r < 2^(1/6) σ,  U(r) = λ U_LJ(r)  beyond, up to `lj_r_cut`.

  So every pair has the same repulsive core, and the attractive well scales with λ.
  HOOMD has no such potential; it is built from two LJ forces: a core with depth
  ε(1 − λ), cut and shifted at its minimum, plus a tail with depth ελ.
- **Electrostatics, Debye-Hückel.** U(r) = qi qj l_B kT e^(−κr) / r, shifted to zero at
  `yukawa_r_cut`. The Bjerrum length l_B = e²/(4π ε0 εr kT) uses the dielectric
  constant of water at the temperature (empirical fit, ~78.5 at 298 K, l_B ≈ 0.71 nm).
  The inverse Debye length κ = √(8π l_B n) for a 1:1 salt of n ions per nm³
  (Debye length 0.79 nm at 150 mM). Pairs not charged on both sides skip this force
  entirely (their energy is exactly zero; skipping it is ~1.5x faster).
- **Bonds**: harmonic, U = k/2 (r − r0)², k = 1920 kcal/mol/nm², r0 = 0.38 nm. Bonded
  pairs have no other interaction.

### 3.4 Dynamics

Langevin dynamics at the set temperature, friction γ = 1 for every residue (HOOMD's
default). The thermostat's random kicks come from the task's seed, so every replicate
is an independent trajectory. The friction is a numerical choice, not water's:
**times in this model compare tasks with each other, they are not real times.**

### 3.5 Geometries

**Bulk.** One chain in a cubic periodic box `box_size` wide (500 nm: the chain never
meets its periodic image). It starts as a straight line; equilibration makes it a coil.

**Grafted.** A flat surface at height z_wall near the bottom of the box. The box is
periodic in x and y. In z it reaches from the surface to a **ceiling** placed at the
chain's contour length (N − 1) × 0.38 nm plus `wall_headroom` (2 nm) above it, then one
interaction range of empty space. So no bead can wrap through the box in z, and none can
interact with the surface through the periodic boundary. The grafted residue (first or
last, `graft_terminus`) is the **anchor**: it sits at its contact distance from the
surface and is left out of the integrator, so it never moves. The other beads start as a
random walk growing upward from it, each bead at least 0.5 nm from every earlier one.

The **surface** interaction, per residue, uses the same split as the residue pairs, with
`wall_attraction` in the role of λ:

- a repulsive core: a 9-3 potential (`wall_potential: lj93`, the potential of a smooth
  solid half-space of LJ sites) or a 12-6 one (`lj126`), with σ = (σ_residue +
  `wall_sigma`)/2 and depth `wall_epsilon`, cut at its minimum (the **contact
  distance**), so a residue there feels no force;
- an optional attractive tail of depth `wall_epsilon` × λ_i, out to `wall_r_cut`, with
  λ_i = `wall_attraction` × λ_residue (`wall_attraction_scale: hydropathy`) or
  `wall_attraction` for every residue (`uniform`). 0 (the default) = purely repulsive.
- an optional **surface charge** `wall_charge_density` σs (e/nm²), linearised
  Poisson-Boltzmann for a plane with salt water on one side:
  U(h) = 4π l_B kT q σs / κ · e^(−κh) for a residue of charge q at height h, with the
  same l_B and κ as the residue pairs. It is accurate while the surface potential
  |e ψ / kT| = 4π l_B |σs| / κ stays below ~1 (at 150 mM: 0.35 for 0.05 e/nm², 0.71
  for 0.1, 1.41 for 0.2). `noodle launch` prints it for every charged set and warns
  above 1. Implemented as a HOOMD Morse wall whose second term is exactly this
  exponential (the derivation is in `charged_surface` in `src/forcefield.py`).

The **ceiling** is always purely repulsive and out of the chain's reach.

**Brush** (`grafting_density` > 0). `chains` chains (100 = 10 x 10) are grafted at the
centres of a square lattice of cells. With grafting density ρ (chains/nm²) the box is
√(chains/ρ) wide: 100 nm at 0.01, 44.7 at 0.05, 18.3 at 0.3, 10 at 1. If the box would
be narrower than 2 × (largest interaction range + 0.5 nm) = 8 nm, a residue could feel
the same neighbour twice through the periodic boundary, so the lattice grows to the next
square number of chains until it is wide enough (never needed with 100 chains, since ρ
is at most 1). Each chain starts grown upward inside its own lattice column, kept
0.25 nm from the column's sides, so no two beads start closer than 0.5 nm. Every other
setting works as for a single grafted chain.

### 3.6 Run length

A chain's relaxation time grows like N^2.2 (Rouse dynamics with excluded volume), so:

steps = multiplier × N^2.2 / dt,

with `rouse_multiplier` (2) for equilibration and `production_multiplier` for production.
A number in `equilibration_steps` / `production_steps` overrides the formula (the `test_*`
sets use tiny fixed lengths). `frames` frames are saved, evenly spaced over production.

| sets | production_multiplier | frames | basis |
|---|---|---|---|
| `grafted*` | 8 | 2000 | measured (below) |
| `baseline` (bulk) | 6 | 3000 | provisional: a free chain should relax faster than a grafted one |
| `brush*` | 3 | 500 | provisional: every frame averages over 100 chains (section 10) |

**Grafted chains, measured** (a `grafted` run at multiplier 6, 1500 frames, 35 proteins,
N = 22-239, one replicate each). The Rg autocorrelation time is τ_Rg = c × N^2.2 tu with
c = 0.006-0.062 (median 0.021, 90th percentile 0.042). The spread is about 1.7x at a
given N and follows how expanded the chain is: expanded chains relax slowly, compact
hydrophobic ones fast (c correlates +0.45 with Rg/√N, −0.55 with `shd`). A free power-law
fit gives N^1.74 over all chains and N^2.45 over N ≥ 100 (±~0.5), consistent with 2.2.
τ_Ree is typically 2-4x τ_Rg (up to 11x). At multiplier 6 every protein got at least 48
independent Rg samples (median 143), Rg errors of 0.4-3%, Ree errors up to 6%, and 1 of 35
drifted. Independent samples = multiplier / (2c), so 8 gives at least 64 (median ~190).
Ree and heights get fewer; replicates are the cheaper way to tighten them.

Run time grows like N^3.2 for a brush (100 × N residues × N^2.2 steps). Rough costs
on one core (1 million residue-steps per second):

| N | single chain | 10 x 10 brush | brush trajectory |
|---|---|---|---|
| 50 | 5 min | 1.9 h | 30 MB |
| 100 | 50 min | 17 h | 60 MB |
| 200 | 8 h | 6.7 days | 120 MB |
| 400 | 3 days | 61 days | 240 MB |

### 3.7 Speed choices

`lj_r_cut` 2.0 nm (the LJ tail beyond 2 nm is worth < 0.005 kT per residue) and
`dt` 0.02 (the time step of the original HPS model; it gives the same temperature and
energies as 0.01). Electrostatics are skipped for pairs that aren't both charged. The
neighbour list is a tree for single chains (mostly empty box) and a cell list for
brushes. One simulation runs on one core; a run gets its speed from many jobs side by
side. (The code would use a GPU if HOOMD were built for one; the group's cluster is
CPU only, and everything is set up and tested for CPUs.)

---

## 4. The pipeline, step by step

### 4.1 `noodle new <name>`

Creates `runs/<name>/run.yaml` from `src/run_template.yaml` and gives the run the next
**number** in `runs/runs.csv`. Names use letters, digits, `. _ -`, and can't be only
digits (those mean run numbers).

### 4.2 `noodle launch <run>`, the first time

1. **Checks** `run.yaml` (unknown or missing settings, with "did you mean"; impossible
   values) and every physics set it names (unknown settings, missing settings,
   invalid choices, grafting density 0-1, square number of chains). Any problem stops
   the launch before anything is written.
2. **Reads the sequences**: the CSV named by `sequences_csv`, from the run folder if it
   is there, else from `sequences/`; the rows in `rows`, counted from 1: `all`, a row
   (`7`), a range (`101-200`), or a list of rows and ranges (`[1-50, 75]`), taken in file
   order. Names are
   cleaned to safe characters and must be unique. Sequences with unknown letters, and
   those longer than `max_residues` (if not 0), are skipped and counted.
3. **Builds the tasks**: every protein × physics set × replicate. Each protein and
   physics set draws its replicates' seeds with `random.sample(range(1, 65536),
   replicates)`, so they are always different. Task ids count from 1, protein by
   protein.
4. **Shows the plan**: sequences, physics sets (with surface potential warnings),
   task count, steps for the shortest and longest protein, rough cost (residues ×
   steps / 1 million per second per core), the number of tasks longer than the time
   limit, trajectory disk space (12 bytes per residue per frame), and the number of
   jobs. Then asks y/N.
5. **Freezes the run**: writes `physics.yaml` (the sets fully written out), `tasks.csv`
   (with every sequence), and `manifest.csv` (one row per protein and physics set:
   identity, conditions, sequence descriptors, section 6). From now on only the cluster
   settings of `run.yaml` (`max_jobs`, `time_limit_days`, `max_running`,
   `extra_sbatch`) still matter; they apply from the next launch.
6. **Groups the tasks into jobs** and starts them (4.3, 4.4).

### 4.3 Grouping tasks into jobs

The grouping uses only what is known exactly before anything runs: each task's
**work**, its residues × the steps it still has to compute (chains × N × steps, with
steps ∝ N^2.2, section 3.6; for a resumed task, the steps left). Speed estimates play
no part, so the same tasks always give the same jobs. Within one physics set, ordering by
work is ordering by chain length.

The run is split into at most `max_jobs` jobs (`run.yaml`, 1000 when not set; 0 = no
limit):

1. With no more tasks than `max_jobs`, every task gets its own job (a small batch runs
   fully in parallel).
2. Otherwise the tasks are sorted from the most work to the least (ties by task id) and
   cut into consecutive groups. Each job is given a share: the work still left divided
   by the jobs still available. It takes the next task, then keeps taking the next one
   while that brings it closer to its share. A task with more work than the share is a
   job of its own; the last job available takes whatever is left.

So the longest proteins run alone, the shortest are pooled many to a job, and every job
holds proteins of similar length and about the same total work. For the 12,000 proteins
of `human_idrs_12k.csv` as single chains, 1000 jobs means: proteins of 331-599 residues
one per job, down to 3942 proteins of 20-37 residues in the last job. The plan shows
this, with a rough warning when a job of several tasks may exceed `time_limit_days` (at
1 million residue-steps per second); it then stops at the limit and continues at the
next launch, like a single long task (4.7). A task longer than the time limit always
needs several launches (4.6). Jobs are numbered after the run's existing ones, from the one with the most work per
task (the longest proteins) to the one with the least, and in each job the longest tasks
come first. Each job gets its file `jobs/job_NNNN.csv`, with a row per task (status
`submitted`), carrying over any earlier progress of the task (attempts, run time, steps).

### 4.4 Starting the jobs

**On a cluster** (when `sbatch` exists), the jobs go to SLURM as job arrays:

```
sbatch --parsable --job-name=<run> --array=<indices>[%max_running] --time=<minutes>
       --chdir=<noodle-md> --open-mode=append --output=runs/<run>/jobs/slurm.out
       <extra_sbatch> src/slurm_job.sh runs/<run> <offset>
```

Array indices must stay below the cluster's `MaxArraySize` (read from `scontrol show
config`), so a run with more jobs is sent as several arrays with an offset each: job
number = offset + array index. Every submitted job is recorded in `jobs/submitted.csv`
(job number, SLURM id like `1234_7`, time), and the command checks the queue.

`src/slurm_job.sh` asks for one node, one task, one CPU, and `--signal=B:USR1@600`
(SLURM sends SIGUSR1 ten minutes before the time limit). It sets one thread for numpy,
sends everything the job prints to `jobs/job_NNNN.out`, and `exec`s the Python job
runner, so Python itself receives SLURM's signals.

New runs get `extra_sbatch: ["--constraint=cascadelake|icelake|sapphirerapids"]` from
`src/run_template.yaml`. It keeps jobs off the 2.1 GHz Skylake and 1.8 GHz nodes. This is
mainly for speed. It also makes per-step times more uniform, which helps the step-time
fit behind `noodle status`. The constraint is passed straight to sbatch with no shell, so
the `|` needs no quoting. The names come from `sinfo -o "%N %f"` on the group's cluster.
The conda-forge HOOMD build is generic x86-64, so it doesn't use AVX-512 itself: the
constraint selects nodes by processor generation and clock, not by instruction set.

**Without SLURM** (your own computer), the jobs run one after another in the terminal,
through the same job runner.

### 4.5 A job (`src/simulate.py`)

The job reads its job file and runs every task not yet completed, in order. For each:

1. **Build**: from the task's checkpoint if there is one (and cut its trajectory back
   to the checkpoint's step), else from the starting configuration (3.5).
2. **Run** in chunks: the first of 1,000 steps, then each sized to take about a minute.
   Frames are written only during production, every `production / frames` steps,
   counted from the start of production. A minute after the start, then every 10
   minutes, the task saves a checkpoint and writes its progress (step, phase, speed, run
   time) into its row of the job file. A checkpoint is one small file per task (a few MB
   even for a brush), replaced each time, so this costs next to nothing.
3. **Stop requests**: SIGUSR1 (time limit), SIGTERM (scancel) and SIGINT (Ctrl-C) only
   note the request. After the current chunk the task saves a checkpoint, records
   `stopped` with the reason, and the job ends (exit code 3), so nothing is lost. A
   second Ctrl-C quits at once (the last checkpoint, at most 10 minutes old, stays).
4. **Done**: the task records `completed` with its frame count and deletes its
   checkpoint.
5. **Analyse** (section 7): the task's averages go into its row of the job file
   (`analysed_at`, then every observable).
   If the analysis fails, the error goes into `analysis_error` and `noodle analyze`
   tries again later; the simulation still counts as completed.

A task that raises an error is recorded as `failed` (short message in `error`, full
traceback in the `.out` file), and the job goes on with its next task.

A resumed task is as valid as an uninterrupted one, but not frame-for-frame identical
(the thermostat's random numbers differ after a restart).

### 4.6 `noodle status [run]`

A task's **current state** is its row in the newest job file that lists it, combined with
the SLURM queue (`squeue --array --name=<run>`, every user's jobs, matched by the SLURM ids
in the run's `jobs/submitted.csv`, so a run launched by someone else in a shared copy is
seen the same way):

| row status | SLURM job in the queue? | shown as |
|---|---|---|
| `completed` | - | done |
| `failed` | - | failed |
| no row | - | not launched |
| `submitted` | yes | waiting |
| `submitted` | no | stopped (its job ended before reaching it) |
| `stopped` | - | stopped (time limit, cancelled, Ctrl-C: `stop_reason`) |
| `running` | yes | running (flagged `stuck?` if silent for over 30 minutes, 3 missed reports) |
| `running` | no | stopped (the job died without saying so) |

Without SLURM, a `running` row counts as running if it reported within 30 minutes, and
`submitted` rows of a job with such a task count as waiting.

The report shows a progress bar, counts (per physics set too), the time left, and the
tasks grouped by state, with failures grouped by error.

**Time left** comes from the speeds measured in this run (unknown until some task has
reported, a minute after it starts). A task that has measured its own speed uses it.
Any other task uses the time per step of the run's measured tasks, which grows with the
number of particles (chains × N): a straight line, overhead + time per particle, is
fitted to them (same physics set if it has two sizes measured) and read at the task's
size. The report then gives two numbers:

- **until done**: how long until the running and waiting tasks are all done. The tasks of
  a job run one after another and jobs run side by side, at most `max_running` at a time
  (one at a time without SLURM): the time of the slowest such line of jobs. Waiting in
  the SLURM queue comes on top.
- **computing in all**: the sum over every unfinished task, i.e. the core-hours still to
  spend, including stopped tasks that need `noodle launch`.

The launch plan gives the same "until done" estimate, at the rough speed of its cost
estimate. It writes `log.csv` (5.3) and refreshes the
run's row in `runs/runs.csv`. Without a run, it refreshes and prints every run.

### 4.7 `noodle launch <run>` again

Sorts the tasks: done (left alone), in the queue (their job's SLURM id is in the
queue: left alone), failed (listed with their errors; restarted only if you agree),
and the rest, to start: resumed from their checkpoint if they have one, else from the
beginning. The tasks to start are grouped into new jobs (4.3) and started (4.4). So
relaunching is always safe, and is the way to continue tasks that hit the time limit.

### 4.8 `noodle cancel <run>`

Finds the run's SLURM ids (from `jobs/submitted.csv`) still in the queue, says how many
are running and waiting, asks, then `scancel`s their arrays (each array belongs to one
run, so this cancels exactly the run's jobs). Running tasks receive SIGTERM and stop
cleanly (4.5); `noodle launch` continues them.

### 4.9 `noodle analyze <run>`

Gathers every analysed task from the job files. Completed tasks without results
(their job's analysis failed, or the job was killed during it) are analysed here, from
their trajectories, into `analysis/reanalysed.csv`, once. Then it writes
`analysis/replicates.csv` (one row per task), the results columns of `manifest.csv`
(replicates pooled: OBSERVABLES.md, section 10) and `analysis/report.json` (5.3), refreshes `runs.csv`, and prints a short report. It takes
seconds, can run at any time, and works the same after the trajectories are deleted
(the results live in `jobs/`).

### 4.10 `noodle export <run> [<run> ...]`

Runs `analyze` on each run, then writes `exports/<name>/`: the runs' `manifest.csv`,
`replicates.csv` stacked with a first column `run`, their rows of
`runs.csv`, and `physics/<run>.yaml`. The folder is written under a temporary name and
renamed at the end, so it is never half-written. It warns about unfinished runs and
about a physics set name that means different settings in different runs.

### 4.11 `noodle clean <run>`

Runs `analyze`, then offers to delete the trajectories (and leftover checkpoints) of
completed, analysed tasks. Before deleting, it adds their SHA-256 checksums to
`trajectories.sha256` (the format `sha256sum -c` reads), and it deletes only after you
type the run's name. Unfinished tasks keep their files.

### 4.12 Run numbers and choosing a run

`runs/runs.csv` has one row per run ever made in this copy of the repository (5.4). A run
folder without a row (copied in, for example) gets one when the file is read. Numbers
are never reused; if a deleted run's name is used again, the new run gets a new number
and the old number says "no longer exists".

Wherever a command needs a run, it accepts, in this order: a run folder's name; a path
to a run folder or anything inside it (relative to where you are); a run number. With
nothing given, it uses the run whose folder you are in. The `noodle` wrapper passes your
current folder to Python in `NOODLE_CWD` before moving to the repository folder.

### 4.13 The wrapper, help and completion

`noodle` (bash) maps each command to its program in `src/` and runs it with the tool's
own Python (`env/bin/python`), with one thread for numpy. `noodle help` prints the
command list; `noodle help <command>` prints the description at the top of that
command's program, so the help can't drift from the code. `noodle-completion.bash`
completes commands, then run names (folders of `runs/` with a `run.yaml`). `setup.sh`
builds `env/` from `src/environment.yml`, runs `noodle check`, and prints the alias and
`source` lines for `~/.bashrc`.

**Shared copies.** When the `noodle-md` folder has the setgid bit (`drwxrwsr-x`, the
usual setup of a group folder), `noodle`, `setup.sh` and `src/slurm_job.sh` set
`umask 002`: every file and folder they and their jobs create is group-writable, and the
setgid bit gives it the folder's group. So any member can edit any run's `run.yaml`,
continue a run someone else launched (a resumed task appends to the existing
trajectory), and update `runs/runs.csv` and `log.csv`. In a personal folder (no setgid
bit) the user's own umask applies, as before.

---

## 5. Files

### 5.1 A run folder, and who writes what

Every file has one writer at any time, and every table is written to a temporary file
that then replaces the old one, so commands on the login node and jobs on the compute
nodes never overwrite each other and never leave half a file.

| file | written by | what |
|---|---|---|
| `run.yaml` | you | the run's settings |
| `physics_sets.yaml` | you (optional) | physics sets only this run uses; may extend the shared ones |
| `tasks.csv` | launch, once | one row per task |
| `physics.yaml` | launch, once | the physics sets used, fully written out |
| `manifest.csv` | launch creates it; analyze rewrites its results columns | one row per protein and physics set |
| `jobs/submitted.csv` | launch, one line per submitted job | job number → SLURM id |
| `jobs/job_NNNN.csv` | launch, just before submitting; then only job NNNN | the job's tasks, their progress and results |
| `jobs/job_NNNN.out` | job NNNN | everything it printed |
| `jobs/slurm.out` | SLURM | SLURM's own messages (normally empty) |
| `log.csv` | status | every task's current state |
| `analysis/` | analyze | pooled results, report, re-analyses |
| `trajectories/task_NNNN_<name>.gsd` | the task's job | production frames (positions, single precision) |
| `checkpoints/task_NNNN.gsd` | the task's job | full state (double precision), until the task is done |
| `trajectories.sha256` | clean | checksums of deleted trajectories |

Tables are plain CSV: numbers with 6 significant digits, empty cells for "none" or
"does not apply".

### 5.2 Task and job tables

**`tasks.csv`**: `task_id`, `sequence_name`, `physics`, `replicate`, `seed`, `length`,
`sequence`.

**`jobs/job_NNNN.csv`**, one row per task of the job:

| column | meaning |
|---|---|
| `task_id`, `sequence_name`, `physics`, `replicate`, `N`, `chains` | which task |
| `status` | `submitted`, `running`, `stopped`, `failed` or `completed` |
| `phase` | `equilibration`, `production` or `done` |
| `step`, `total_steps` | progress |
| `steps_per_second` | last measured speed |
| `run_time_s` | run time so far, over all attempts |
| `attempts` | how many times the task was started (resumes included) |
| `frames` | frames saved (when completed) |
| `host`, `started_at`, `last_update`, `completed_at` | where and when |
| `stop_reason`, `error` | why it stopped or failed |
| `analysed_at`, `analysis_error` | the in-job analysis |
| then every observable and its `_err`, and the quality checks | OBSERVABLES.md |

**`jobs/submitted.csv`**: `job`, `slurm_job`, `submitted_at`.

### 5.3 Status and analysis tables

**`log.csv`**: `task_id`, `sequence_name`, `physics`, `replicate`, `N`, `state` (4.6),
`percent`, `phase`, `job`, `slurm_job`, `steps_per_second`, `run_time_h`, `time_left_h`,
`attempts`, `host`, `started_at`, `last_update`, `completed_at`, `stop_reason`, `error`.

**`manifest.csv`**, one row per protein and physics set, in four groups of columns:

1. identity: `name`, `physics`, `sequence`, `source_csv`, `source_row`
2. conditions: `geometry`, `temperature_K`, `ionic_concentration_M`, `graft_terminus`,
   `wall_charge_density`, `wall_attraction`, `grafting_density` (grafted only),
   `chains`
3. descriptors (OBSERVABLES.md, section 4): `N`, `net_charge`, `ncpr`, `fcr`, `f_plus`, `f_minus`,
   `mean_hydro`, `kappa`, `scd`, `shd`
4. results (OBSERVABLES.md, sections 5–10): `replicates`, each observable with its `_err`,
   `reduced_density`, `tau_rg_tu`, `tau_ree_tu`, `independent_samples`, `drift`

**`analysis/replicates.csv`**: `task_id`, `name`, `physics`, `replicate`, `seed`, `N`,
`chains`, `frames`, `run_time_h`, each observable with its block error, the diagnostics.

**`analysis/reanalysed.csv`**: like the job files'
results, for tasks analysed by `analyze` instead of their job.

**`analysis/report.json`**:

```
run, written_at
tasks:         total, completed, analysed, failed, could_not_analyse
physics_sets:  {set: rows, analysed_rows, observables: {column: {min, median, max}}}
look_at:       drift: ["name (physics)", ...]
               few_independent_samples: [...]          (fewer than 20, all replicates together)
               could_not_analyse: [{task_id, name, physics, error}, ...]
```

### 5.4 Repository-wide files

**`runs/runs.csv`**: `number`, `name`, `description`, `created`, `launched`, `finished`,
`proteins`, `physics_sets`, `replicates`, `tasks`, `jobs`, `done`, `failed`, `analysed`,
`core_hours` (run time so far, all tasks), `trajectories_gb` (on disk now),
`last_checked`. Written only by `noodle` commands (new, launch, status, analyze), never
by jobs.

**`exports/<name>/`**: `manifest.csv`, `replicates.csv` (first column `run`),
`runs.csv`, `physics/<run>.yaml`.

**`physics.yaml`**: named physics sets. A set may `extends: <other set>` and list only
what changes. The sets `baseline` (bulk) and `grafted` list every allowed setting.
**`sequences/`**: the shared sequence sets (`human_idrs_12k.csv`, `mobidb_idps_15k.csv`),
each ordered so that every row is as different as possible from the rows above it.
Next to the sequence and where it comes from (`pick` is the row's place in that order),
each row has the descriptors of section 6, computed with
`src/descriptors.py`, so they equal the manifest's. The engineered features the order was
originally picked with (`scd_pos`, `scd_neg`, `contact_energy`, `patch_neg`, `patch_neu`,
`frac_sticky`, `log_seq_length`, `radius_at_pick`) were removed: they were built to
predict Rg, not to describe sequences. The row order is unchanged.

---

## 6. Sequence descriptors

`src/descriptors.py` computes them at launch, from the sequence alone, with the model's
own charges (histidine +0.5) and hydropathies, and writes them into `manifest.csv`:
`N`, `net_charge`, `ncpr`, `fcr`, `f_plus`, `f_minus`, `mean_hydro`, `kappa`, `scd`,
`shd`. The formulas and column names follow the group's idp-rg
project, so the tables join directly; the values agree with idp-rg's to float precision
on both shared sequence sets. Their definitions: OBSERVABLES.md, section 4.

---

## 7. Observables and statistics

`src/observables.py` measures each finished task, inside its job, on the compute node,
from its trajectory (4.5). Every column it writes is defined exactly in
[OBSERVABLES.md](https://github.com/amgualdron/noodle-md/blob/HEAD/OBSERVABLES.md), with its error bar (section 8 there), the quality
checks (section 9) and how replicates are combined into the manifest (section 10).

How the measurement runs:

- Frames are read 50 at a time, so memory stays small. Each chain is first made whole
  (its periodic images joined back along its bonds).
- Every quantity is computed per frame (and per chain), then averaged over the frames
  (and chains). The per-frame series give the block errors and the quality checks.
- The pairwise part (Rh, ν, contacts) is O(N²) per chain per frame, done for groups of
  chains of up to 10⁶ residue pairs at a time (about 50 MB). A brush's contacts between
  chains use a k-d tree, periodic in x and y.
- Cost: a 300-residue chain with 3000 frames takes seconds; a 10 x 10 brush with 500
  frames about half a minute.

**Only averages are kept.** Curves (density against height, distance against separation
along the chain, values per residue) are not computed: for a run of thousands of tasks
they would be thousands of times larger than the tables, while the trajectories already
hold everything needed to compute them later, for any task, on any computer
(OBSERVABLES.md, section 11 shows how).

---

## 8. How it is checked

**`noodle check`** (`src/check_physics.py`, also run by `setup.sh`) verifies the model
against its formulas and runs the whole path on short simulations:

- residue pairs: energies of 8 pairs at 6 distances against the Ashbaugh-Hatch +
  Debye-Hückel formulas; no interaction beyond cutoffs for neutral pairs;
- brush geometry: box width, chain count growth when the box would be too narrow;
- electrostatics: Bjerrum length (~0.71 nm) and Debye length (~0.79 nm at 150 mM);
- surfaces: repulsive walls give zero force beyond contact, only push away, go
  smoothly to zero, push back a bead that slipped past; attractive walls have the
  right well depth and no force at its bottom, arginine (λ = 0) feels none; the
  charged surface matches its analytic energy and force, and attracts K, repels E;
- simulations: bulk, grafted (attractive and charged, both termini), charged wall and
  two brushes, run as one job: every task completes with finite positions, anchors
  never move, beads stay between surface and ceiling;
- analysis: every task analysed in its job, manifest columns per geometry, brush
  columns filled with sensible values, reduced density.

The observables were also checked against shapes with known answers (straight rod,
ring, ideal random walk with ν ≈ 0.5, upright grafted rod, regular brush),
interchain contacts against a brute-force count, and Rg, κ² and R_ee against MDAnalysis
on simulated trajectories (identical). The descriptors were checked against idp-rg's
values on both shared sequence sets, and κ on the published limiting cases of
Das & Pappu.

---

## 9. The code, its rules, and how to extend it

### 9.1 Map

```
README.md               using it          REFERENCE.md   how it works (this file)
OBSERVABLES.md          every result column, exactly    CHANGELOG.md   what changed
noodle                  the command (bash): maps commands to programs, help
noodle-completion.bash  Tab completion
setup.sh                builds env/, runs the check, prints the ~/.bashrc lines
physics.yaml            the physics sets
src/forcefield.py       the model: residue table, electrostatics, every force
src/simulate.py         the job runner: build, run, checkpoint, stop, analyse
src/observables.py      the analysis of one trajectory
src/descriptors.py      sequence descriptors
src/runs.py             shared helpers: CSV reading/writing, tasks and jobs, run numbers,
                        choosing a run, physics conditions, run lengths, time estimates, SLURM queue,
                        the names of the result columns
src/launch_run.py       noodle launch: checks, plan, tasks, manifest, jobs, submission
src/status.py           noodle status (and log.csv)
src/analyze_run.py      noodle analyze
src/export_runs.py      noodle export
src/cancel_run.py       noodle cancel
src/clean_run.py        noodle clean
src/new_run.py          noodle new        src/run_template.yaml   the run.yaml it copies
src/check_physics.py    noodle check
src/slurm_job.sh        the SLURM job script
src/environment.yml     the packages of env/ (HOOMD-blue CPU, gsd, numpy, scipy, pandas,
                        pyyaml, MDAnalysis for users' own analyses)
```

Start reading at `src/forcefield.py`, whose header explains how HOOMD represents forces.
Each program starts with a description of what it does; `noodle help <command>` prints it.

### 9.2 Rules the code follows

- **No command-line options.** Every setting lives in `run.yaml` or `physics.yaml`, or
  is an automatic default. Commands take only run names or numbers.
- **Plain, readable outputs.** CSV and JSON only, uncompressed, readable without code.
- **Readable code.** Written for readers new to software: short functions, comments
  that say why, no clever tricks.
- **One writer per file**, atomic replacement (5.1).
- **Frozen runs.** A launched run never changes when `physics.yaml` or `sequences/`
  change; it keeps its own copies.
- **The physics is fixed** unless a change is a deliberate, checked model change.
- **Compute only.** Compute nodes never submit jobs; code is never copied into run
  folders; CPU only.
- **Launching stays simple**: `noodle launch` shows the plan and asks once.

### 9.3 How to extend it

**A new physics setting.** Add it with its default to `baseline` (any geometry) or
`grafted` (grafted only) in `physics.yaml`: those two sets define which settings exist,
so launch accepts it and reports misspellings. Use it in `src/forcefield.py` or
`src/simulate.py`. If it is a condition worth a manifest column, add it to
`CONDITION_SETTINGS` in `src/runs.py`. Add a test to `src/check_physics.py`.

**A new observable.** Compute it in `analyse_task` (`src/observables.py`) with the
`add(name, values)` helper (values per frame and chain) and `average(name, ...)`, add its
name to the right list in `src/runs.py` (`CHAIN_OBSERVABLES`, `GRAFTED_OBSERVABLES`,
`BRUSH_OBSERVABLES`: this sets the column order; they live there so that `analyze` and
`export` never load numpy), and define it in OBSERVABLES.md (and in the list at the top
of `src/observables.py`). `analyze` pools it automatically. Tasks analysed before the
change don't have it: to add it, empty their `analysed_at` cells in the job files (of
jobs that have finished) and run `noodle analyze`, which re-analyses them into
`analysis/reanalysed.csv` (their trajectories must still be there).

**A new command.** Write `src/<program>.py` beginning with a description (shown by
`noodle help <command>`), using `run_from_command_line` from `src/runs.py` to get the
run. Add it to `program_of` and the help text in `noodle`, and to `_NOODLE_COMMANDS` in
`noodle-completion.bash`.

**A new geometry** touches `system_geometry`, the start configuration and `add_forces`
in `src/simulate.py`, the checks in `src/launch_run.py`, `brush_size`/`conditions` in
`src/runs.py`, and the height-based parts of `src/observables.py`.

---

## 10. Limitations and open items

- **Brush run lengths are provisional.** `production_multiplier: 3`, `frames: 500` for
  brushes come from relaxation times measured on 4 x 4 brushes, not yet on 10 x 10
  ones. Measuring equilibration and relaxation times of 10 x 10 brushes over lengths
  and densities is the main open task; until then, check `drift` and
  `independent_samples` on every brush row. Dense brushes (~0.3 chains/nm²) likely
  need a longer production.
- **Bulk run length is provisional.** `baseline` uses `production_multiplier: 6`
  (down from 22), not yet measured on free chains: run `length_test_bulk` on varied
  proteins and check `independent_samples` and `drift` (physics.yaml explains how).
- **Model times are not real times** (3.4): compare τ's between tasks only.
- **Rh** from the Kirkwood formula is an approximation for flexible chains; compare
  trends, not absolute values, with experiments.
- **κ above 1** for ~2% of sequences (OBSERVABLES.md, section 4).
- **Linearised surface electrostatics** overestimate the pull above ~1 kT/e surface
  potential; launch warns.
- **`max_jobs` defaults to 1000**, which fits one job array on most clusters
  (`MaxArraySize` is often 1001). With a large `max_jobs`, a job can hold very few
  tasks of the shortest proteins and HOOMD's start-up (seconds) is paid more often; with
  a small one, jobs of many short tasks may need more than one time limit.
- **`max_running` applies per job array.** A run with more jobs than the cluster's
  `MaxArraySize` is sent as several arrays, each with its own limit.
- **Run numbers are per copy of the repository.** Each user's clone numbers its own runs.
- **On your own computer** the jobs run one at a time.
- **Re-analysing tasks** that were already analysed has no command of its own (9.3).
- **To verify on the cluster** (they can only be simulated off it): SLURM's SIGUSR1
  ten minutes before the time limit reaching the task (a run with a very short
  `time_limit_days`), `noodle cancel` on running and pending jobs, and a run with more
  jobs than `MaxArraySize`.

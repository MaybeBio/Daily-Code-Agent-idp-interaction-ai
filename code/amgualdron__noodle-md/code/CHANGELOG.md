# Changelog

What changed between versions of noodle-md, for anyone who pulls an update.
Runs launched with one version can't be continued or analysed by another version
with a different run layout: finish them with the version that launched them
(`git checkout v1`, for example).

## v2 (in progress, branch `release`)

### Added
- `runs/<name>/manifest.csv`, written at the first launch: one row per protein and
  physics set with its name, sequence, where it came from (`source_csv`, `source_row`),
  the physics conditions (`geometry`, `temperature_K`, `ionic_concentration_M`,
  `graft_terminus`, `wall_charge_density`, `wall_attraction`, `grafting_density`,
  `chains`) and the sequence descriptors (`src/descriptors.py`): `N`, `net_charge`,
  `ncpr`, `fcr`, `f_plus`, `f_minus`, `mean_hydro`, `kappa`, `scd`, `scd_pos`,
  `scd_neg`, `shd`. HPS1 charges (histidine +0.5) and idp-rg's formulas and names, so
  the values agree with idp-rg (checked: to float32 precision on both shared sets).
  `kappa` uses Das & Pappu's block-arrangement delta_max (the E/K test sequences give
  0.0009 and 1, as published); it can exceed 1 for ~2% of sequences, mostly ones with
  few charges.

- **Every task is analysed by its own job, right after it finishes** (`src/observables.py`;
  seconds for a chain, about half a minute for a brush). Averages, each with a
  block-average error bar (`<name>_err`), go into the task's row of its job file.
  Lengths in nm, times in model units (tu). Every column is defined in `OBSERVABLES.md`.
  - every chain: `rg_nm`, `rg_sq_nm2`, `rg_sd_nm`, `rh_nm` (Kirkwood), `ree_nm`,
    `ree_sq_nm2`, `asphericity`, `shape_anisotropy` (relative shape anisotropy kappa^2,
    not the sequence kappa), `prolateness`, `transverse_rms_nm`, `nu` and `r0_nm`
    (R_ij = r0 |i-j|^nu for |i-j| >= 5), `bond_cos`, `persistence_nm` (Flory's
    projection), `contacts_per_residue` (|i-j| >= 3, closer than 1.5 x the mean diameter)
  - grafted chains and brushes: `com_height_nm`, `end_height_nm`, `max_height_nm`,
    `height_95_nm`, `rg_perp_nm`, `rg_par_nm`, `end_lateral_rms_nm`, `tilt_cos`,
    `bond_order` (P2 of the bonds with the normal), `adsorbed_fraction`
  - brushes: `brush_height_nm` (2 x mean residue height), `concentration_mg_ml` (below
    height_95), `interchain_contacts_per_residue`, `lateral_heterogeneity`; and, in the
    manifest, `reduced_density` (grafting density x pi <Rg^2> of the same protein as a
    single chain in the same run)
  - diagnostics: `tau_rg_tu`, `tau_ree_tu`, `independent_samples`, `drift`
  Checked on shapes with known answers (rod, ring, random walk, upright grafted rod,
  synthetic brush), interchain contacts against brute force, and Rg, kappa^2 and Ree
  against MDAnalysis on simulated trajectories (identical).
- **Run numbers**: `runs/runs.csv` lists every run with a number (given by `./noodle new`,
  never reused) and run-wide figures: description, created / launched / finished dates,
  proteins, physics sets, replicates, tasks, jobs, done, failed, analysed, core-hours,
  trajectory size (GB), last checked. `new`, `launch`, `status` and `analyze` keep it up
  to date (never the jobs). Runs made before it get a number the first time it is read.
  `./noodle status` without a run prints it.
- **Choosing a run without typing its name**: wherever a run is asked for, give its name,
  its number, a path to its folder (what Tab completes to, e.g. `runs/my_run/`), or
  nothing when you are inside its folder (or a subfolder).
- `./noodle cancel <run>`: counts the run's jobs in the queue, asks, then `scancel`s them.
  Running tasks save a checkpoint first; `./noodle launch` continues the run later.
- `./noodle help` lists the commands; `./noodle help <command>` explains one in detail
  (the description at the top of its program). Unknown commands show the list too.
- Tab completion (`noodle-completion.bash`): commands, then run names. `setup.sh` ends by
  printing two lines to add to `~/.bashrc`: an alias, so `noodle` works from any folder,
  and the line that loads the completion.
- `./noodle export <run> [<run> ...]`: brings each run's analysis up to date, then writes
  `exports/<name>/` with the runs' `manifest.csv` and `replicates.csv`
  merged (first column `run`), their rows of `runs.csv`, and `physics/<run>.yaml`. One
  run: named after it; several: it asks, suggesting the common part of their names. Warns
  about unfinished runs and about a physics set name with different settings in different
  runs. `exports/` is ignored by git.
- scipy is listed in `src/environment.yml` (it was already installed through MDAnalysis).

### Changed
- **One copy can be shared by a research group.** In a folder with the setgid bit (a group
  folder, `drwxrwsr-x`), `noodle`, `setup.sh` and the jobs make everything group-writable
  (`umask 002`), and the SLURM queue is read for every user (`squeue` without `--user`;
  jobs are still matched by the run's own SLURM ids). Before, a run launched by one member
  looked "stopped" to another, whose `noodle launch` would then have started its queued
  tasks a second time. README section 6 says how to set a group up.
- The shared sequence sets (`human_idrs_12k.csv`, `mobidb_idps_15k.csv`) carry the standard
  sequence descriptors `N`, `net_charge`, `ncpr`, `fcr`, `f_plus`, `f_minus`, `mean_hydro`,
  `kappa`, `scd`, `shd`, computed with `src/descriptors.py`, the same code as the
  manifest. The engineered Rg-prediction features are gone (`scd_pos`, `scd_neg`,
  `contact_energy`, `patch_neg`, `patch_neu`, `frac_sticky`, `log_seq_length`,
  `radius_at_pick`). Rows, ids and sequences are unchanged, so `rows:` picks the same
  proteins as before. `scd_pos` and `scd_neg` are also gone from `manifest.csv`.
- Shorter productions, from a measurement of relaxation times (REFERENCE.md, section 3.6):
  grafted chains `production_multiplier` 22 -> 8 and `frames` 3000 -> 2000 (still at
  least 64 independent Rg samples per protein); bulk (`baseline`) 22 -> 6, provisional
  until measured. Brushes unchanged (3, 500). Runs already launched keep their lengths.
- `length_test` (the grafted measurement) is replaced by `length_test_bulk`: `baseline`
  with 6000 frames, to measure the bulk run length the same way.
- New runs use only the cluster's faster nodes:
  `extra_sbatch: ["--constraint=cascadelake|icelake|sapphirerapids"]` is now the default
  in `run.yaml`. It leaves out the Skylake and 1.8 GHz nodes. `extra_sbatch: []` allows
  every node again. Existing runs keep their own `run.yaml` until you add the line.
- Run names can't be only digits (those are run numbers now).
- `rows` also takes a list of rows and ranges, `[1-50, 75, 200-210]`: in `run.yaml` a
  list in brackets always means "several". A range is still written with a dash, `1-50`.
- README rewritten short and plain, for users with no software background: setup,
  a first run in six steps, the commands, the physics sets, the main result columns,
  and fixes for common problems. `src/run_template.yaml` explains the per-job time limit.
- New `REFERENCE.md`: the model, the pipeline step by step, every file and column, how
  it is checked, the code's rules, how to extend it, and the open items.
- New `OBSERVABLES.md`: the exact definition of every column of `manifest.csv` and
  `analysis/replicates.csv` (descriptors, observables, error bars, quality checks, how
  replicates are combined), and how to compute anything else, such as curves, from the
  trajectories, with tested examples.
- **No more profiles (curves).** Density against height, distances against separation
  and per-residue values are no longer computed or saved: an 18-task test already wrote
  ~11,000 rows of them, and a run of thousands of tasks would write millions. The
  trajectories hold everything needed to compute them afterwards, on any computer
  (`OBSERVABLES.md`, section 11). Gone: `jobs/job_NNNN_profiles.csv`,
  `analysis/profiles.csv`, `analysis/reanalysed_profiles.csv`, and `profiles.csv` in
  exports. `noodle analyze` no longer needs pandas.
- `./noodle check` runs its test simulations as one job and checks the new tables
  (in-job analysis, manifest columns per geometry, brush columns, reduced density).
- **`./noodle analyze` pools instead of recomputing**: it gathers the job files into
  `analysis/replicates.csv` (one row per task) and `analysis/report.json` (counts, the range of
  each observable per physics set, rows that drift or have fewer than 20 independent
  samples, tasks that could not be analysed), and fills the results columns of
  `manifest.csv`: `replicates`, then every observable averaged over replicates, with the
  standard error between replicates (or a lone replicate's block error). It takes
  seconds. A task its job could not analyse is analysed by `analyze` instead
  (`analysis/reanalysed.csv`). `analysis/summary.csv`, `analysis/timeseries/` and
  `analysis/profiles/` are gone: the per-task results in `jobs/` replace them, and keep
  working after `./noodle clean` deletes the trajectories.
- **Tasks are grouped into jobs.** A job is one SLURM array element (one core) that runs
  its tasks one after another. A run is split into at most `max_jobs` jobs (new
  `run.yaml` setting, default 1000, the value when missing; 0 = one job per task) by
  chain length: tasks are sorted by their work (chains x N x steps, known exactly before
  running, no speed estimates), and cut into consecutive groups of about equal work. The
  longest proteins get a job each, the shortest share one, and a run with no more tasks
  than `max_jobs` gets one job per task, so small batches run fully in parallel. The
  same tasks always give the same jobs. The jobs with the longest proteins get the lowest
  numbers, so they start first. A failed task doesn't stop its job. Thousands of
  minute-long tasks no longer flood the queue, and HOOMD is loaded once per job.
  `max_jobs` can be changed after launching, like the other cluster settings. `max_running` now counts jobs (one core each). Launching looks the same;
  the plan's Cluster line says how many jobs.
- **New bookkeeping layout** (fewer files): `logs/` (one JSON per task) and `slurm/` (one
  output file per task) are replaced by `jobs/`:
  - `jobs/job_0007.csv`: the tasks of job 7 and how each is going (status, phase, step,
    speed, run time, attempts, host, dates, error), written by launch just before the job
    is submitted, then only by the job itself. A task's current state is its row in the
    newest job file that lists it (a resumed or retried task gets a row in a new job).
  - `jobs/job_0007.out`: everything job 7 printed (all of its tasks).
  - `jobs/submitted.csv`: which SLURM job runs which job.
  - `jobs/slurm.out`: SLURM's own messages (normally empty).
- `./noodle status` writes `log.csv` (one row per task: state, progress, job, SLURM id,
  speed, run time, time left, attempts, host, dates, error) instead of `status.txt`.
  Tasks waiting in a job that is still running other tasks count as waiting.
- A task that is about to hit the cluster time limit stops cleanly: SLURM warns it 10
  minutes ahead (`#SBATCH --signal=B:USR1@600` in `src/slurm_job.sh`), and the task saves
  a checkpoint and records `stopped (time limit)`. Nothing is lost (before: up to an
  hour). `scancel` and Ctrl-C (own computer) stop it the same way; a second Ctrl-C quits
  at once.
- The simulation runs in chunks of about a minute (sized from the measured speed)
  instead of 10,000 steps, so a stop request is always answered within about a minute.

### Fixed
- `noodle analyze`, `export` and `clean` seemed to freeze for a minute or more on the
  cluster, printing nothing: they loaded HOOMD, which is slow from a shared filesystem,
  without needing it. HOOMD is now loaded only where forces are built (`forcefield.py`),
  so these commands start in about a second. `analyze` says what it is doing from its
  first line, and warns before loading HOOMD to re-analyse tasks their job didn't.
- `noodle analyze`, `export` and `clean` no longer load numpy, scipy or gsd either (about
  130 Python modules instead of 580), which took up to half a minute on the cluster's
  shared filesystem. The result column names moved from `observables.py` to `runs.py`,
  and the replicates are pooled with Python's `statistics` module (same numbers). The
  analysis code is loaded only to re-analyse a task its job didn't.
- `./noodle status` no longer crashes on a task in its first minute (no run time yet).
- Progress is reported every 10 minutes instead of every hour (with a checkpoint, which
  costs next to nothing), so `noodle status` percentages move between checks; a task
  silent for 30 minutes is flagged `stuck?`.
- `noodle status` "Time left" is now how long until the running and waiting tasks are
  done, with jobs running side by side (at most `max_running` at a time), instead of the
  total computing time, which was easily read as days of waiting; the total is still
  shown below it. Tasks without a speed of their own are estimated from a fit of time per
  step against particle count over the run's measured tasks, instead of borrowing the
  steps per second of the measured task closest in length. The launch plan also shows how
  long until all are done.
- On your own computer, the tasks waiting in a running job show as waiting, not stopped.
- Replicates of the same protein under the same physics set always get different
  random seeds (before, two could draw the same seed and be the same simulation).

### Removed
- `./noodle calibrate` and `experiments/` (the brush calibration experiment). It was a
  private test harness and is not part of the released tool. It lives on in the tag
  `v1` and the local branch `calibration`.
- A run's own `physics_sets.yaml` is still supported (sets that only one run needs).

## v1 (tag `v1`)

The first version used on the cluster: bulk, grafted and charged-wall single chains,
10 x 10 brushes, `./noodle new | launch | status | analyze | clean | check | calibrate`.

# Coarse-grained simulation of the archive's disordered proteins

The design, written before any run, is [`docs/SIMULATION_DESIGN.md`](https://github.com/Natalija-Stepurko/bioemu-saxs-campaign/blob/HEAD/docs/SIMULATION_DESIGN.md): CALVADOS 3
single-chain simulations (OpenMM, CPU) of the 41 clean disordered profiles with a BioEmu-1 ensemble, under each
entry's recorded temperature, pH and ionic strength; the simulated radius of gyration against the measured one,
paired with BioEmu-1's; a computed SAXS curve per ensemble (PULCHRA back-mapping, Pepsi-SAXS as for the
archive's curves) scored as the archive's curves are; a convergence check; and a circularity check against the
proteins CALVADOS was trained and validated on.

| Script | Environment | Does |
|---|---|---|
| `01_prepare_inputs.py` | main | the input table: proteins, construct sequences, temperature, pH and ionic strength from the SASBDB record by the design's rule (`simlib.conditions`), denaturants, run settings, seeds |
| `02_fetch_sources.py` | simulation | at pinned versions: the CALVADOS 3 residue parameters, PULCHRA 3.06 (compiled with gcc), Pepsi-SAXS 3.0, and the CALVADOS 1, 2 and 3 training and validation protein sets |
| `03_circularity.py` | simulation | each simulated protein against every training and validation protein (local alignment; identity over the shorter sequence; flag at ≥ 90%) |
| `04_plan.py` | simulation | the compute plan from the pilot timings: a power law of speed against chain length, each protein's projected time, and the design's subset rule (§8) |
| `05_run.py` | simulation | the simulations: one single-thread process per protein, at most `--workers` at a time, in chunks of 50 frames with a checkpoint after each, so an interrupted run continues; `--frames N` stops after N frames (the pilot) |
| `06_saxs.py` | simulation | 100 evenly spaced frames per finished run, back-mapped and passed to Pepsi-SAXS with the measured profile; the conformer curves in the archive's format |
| `07_analyse.py` | main | Rg of every production frame, ratios, convergence, χ², NRMSD and curve Rg scored with `bsc.saxs`, the paired comparison with BioEmu-1 over all proteins and the design's subsets, E10 and E11, run metadata, and the page figure |

`simlib.py` holds the functions the scripts share (input table, buffer rule, run settings, radius of gyration,
chain unwrapping, convergence, the Guinier radius of a Gaussian chain); `tests/test_simulations.py` tests them.

Tracked outputs in `results/simulations/`:

| File | Contents |
|---|---|
| `inputs.csv` | per protein: sequence, measured Rg, conditions with their source, parsed buffer species, denaturant, box, save interval, frames, seed |
| `plan.json` | the pilot timings, the speed fit, projected times, frames per protein and the proteins left out |
| `sources.json` | URLs, commits, versions and sha256 of every download |
| `calvados_sets.csv` | the CALVADOS training and validation proteins with sequences, by set |
| `circularity.csv` | per protein: nearest CALVADOS set protein, identity, flag, the sets with a hit at ≥ 90% |
| `per_protein.csv` | per protein: simulated and measured Rg, CALVADOS and BioEmu-1 ratios, χ², NRMSD, curve Rg, conditions, convergence, flags |
| `convergence.csv` | per protein: Rg of five consecutive blocks, relative standard error, difference between halves |
| `summary.json` | the paired comparison over all proteins and each subset, the sensitivity analyses, χ², E10 and E11 |
| `run_metadata.json` | per protein: frames, simulated time, wall time, steps per second; totals |

Trajectories, back-mapped frames and per-frame curves stay in `SIM_WORK` (not tracked).

## Rerun

```bash
uv venv --python 3.12 /path/to/venv-sim
uv pip install --python /path/to/venv-sim --no-deps -r analysis/simulations/requirements.txt
SIM_PY=/path/to/venv-sim/bin/python SIM_WORK=/path/with/10GB/free BSC_ROOT=<where bsc fetch put data/> \
  WORKERS=6 analysis/simulations/run_all.sh
```

Every step resumes. Pepsi-SAXS is distributed for academic use only; `02_fetch_sources.py` downloads it from its
authors and it is not redistributed here. PULCHRA is MIT-licensed; CALVADOS is GPL-3.0 and installed from its
repository at the pinned commit.

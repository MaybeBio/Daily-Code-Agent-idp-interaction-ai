# Installation

## Requirements

| What | Needed for |
|---|---|
| Python ≥ 3.10 (RDKit, NumPy, PyYAML are installed with the package) | everything |
| GROMACS (`gmx`) | running the simulations and `JOB = analyze`; its `oplsaa.ff` for `FF = oplsaa` |
| Packmol (optional) | `PACKER = packmol` |
| AmberTools (optional) | `FF = gaff2` |

## On your computer

The simplest way is one conda environment with everything, including GROMACS. Get conda from
[Miniforge](https://conda-forge.org/download). On Windows, do this inside WSL (Ubuntu): GROMACS
does not run natively on Windows.

1. Install (once), in the Zpoly folder:

   ```
   conda create -n zpoly -c conda-forge python=3.12 rdkit numpy pyyaml gromacs packmol ambertools
   conda activate zpoly
   pip install -e .
   zpoly help
   ```

   - `pip install -e .` links the software to this folder: later changes (new monomers, updates)
     work without reinstalling.
   - `packmol` and `ambertools` are optional: leave them out if you only use `PACKER = emc` and
     `FF = oplsaa`.
   - Without conda: `python -m venv ~/envs/zpoly`, `source ~/envs/zpoly/bin/activate`,
     `pip install -e .`, and install GROMACS separately.

2. Each time:

   ```
   conda activate zpoly
   cd ~/my_simulations                # any working folder
   zpoly template slab > slab.inp     # or: build, bulk, dft, analyze; edit the file
   zpoly slab.inp                     # builds the system folder
   bash DMAPS_n12_KPS_slab/run.sh     # runs GROMACS (threads: NT=8 bash ...)
   ```

   After the run: an input with `JOB = analyze` and `OUTPUT = DMAPS_n12_KPS_slab`, then
   `zpoly analyze.inp`.

3. Check:

   ```
   zpoly list forcefields             # oplsaa: available; gaff2: available with AmberTools
   zpoly examples/inputs/build.inp    # a 12-mer with its force field, in seconds
   ```

   If `oplsaa` is not available, GROMACS is not in the active environment: `conda activate zpoly`,
   or set `OPLSAA_FF = /path/to/oplsaa.ff` in the input file.

## On an HPC cluster (your own space, no admin rights)

1. Copy the package (one file, with all built-in data):

   ```
   scp dist/zpoly-*.whl user@cluster:~/
   ```

   (To make a new one: `pip wheel . --no-deps -w dist`.)

2. Make your own environment:

   ```
   ssh user@cluster
   module load python/3.11                # any Python >= 3.10 (module avail python)
   python -m venv ~/envs/zpoly          # or /project/<you>/envs/zpoly
   source ~/envs/zpoly/bin/activate
   pip install --upgrade pip
   pip install ~/zpoly-*.whl    # also installs RDKit, NumPy, PyYAML
   zpoly help
   ```

   Without a recent Python, use micromamba in your home folder:

   ```
   curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
   ./bin/micromamba create -p ~/envs/zpoly -c conda-forge python=3.12 rdkit numpy pyyaml
   ./bin/micromamba run -p ~/envs/zpoly pip install ~/zpoly-*.whl
   ```

   Install on a login node (it usually has internet); the compute nodes then only use it.

3. Each session (and at the top of `run_slurm.sh`, after the `#SBATCH` lines):

   ```
   source ~/envs/zpoly/bin/activate
   module load gromacs
   ```

4. Update later: `pip install --force-reinstall --no-deps ~/zpoly-<new>.whl`.

## Settings

| Variable / keyword | Use |
|---|---|
| `GMX` | GROMACS command for `run.sh` (default `gmx`; e.g. `export GMX=gmx_mpi`) |
| `NT` | threads for `mdrun` (`run_slurm.sh` sets it from Slurm) |
| `OPLSAA_FF = path` or `ZPOLY_OPLSAA_FF` | your `oplsaa.ff`, if it is not found next to `gmx` |
| `GMXDATA` | GROMACS data folder (for `oplsaa.ff`, water boxes), if not found automatically |

## Prepare locally, simulate on the cluster

1. Build the system on your computer: `zpoly slab.inp` → e.g. `DMAPS_n12_KPS_slab/`.
2. Send the folder (it has only relative paths):
   `rsync -av DMAPS_n12_KPS_slab/ user@cluster:/scratch/user/DMAPS_n12_KPS_slab/`
3. On the cluster: edit the `#SBATCH` lines of `run_slurm.sh`, then `sbatch run_slurm.sh`
   (one temperature per job after equilibration: `bash run_production.sh 300`).
4. Analyse — preferably on the cluster (`zpoly analyze.inp`, small results:
   `phase_diagram.csv`, `density_profile_<T>K.csv`, or `tg_result.yaml`, `tg_density.csv`), or bring
   back only what the analysis reads:

   | System | Files |
   |---|---|
   | slab | `protocol.yaml`, `index.ndx`, `prod_<T>K.tpr`, `prod_<T>K.xtc` |
   | bulk | `protocol.yaml`, `cool_<T>K.edr` |

   ```
   rsync -av --include='protocol.yaml' --include='index.ndx' --include='prod_*K.tpr' \
         --include='prod_*K.xtc' --exclude='*' user@cluster:/scratch/user/DMAPS_n12_KPS_slab/ DMAPS_n12_KPS_slab/
   ```

   The local GROMACS must be the same version as the cluster's, or newer (it reads the `.tpr`).

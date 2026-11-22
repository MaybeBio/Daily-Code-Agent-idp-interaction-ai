## Installation / Environment Setup (ubelix server)

To run the simulations on **ubelix**:

### 1. Load CUDA:

```bash
module load CUDA/11.8.0
```

### 2. Create the OpenMM conda environment:

```bash
conda env create -f utils/environment.yml -n openmm
conda activate openmm
```
Ensure that the openmm is activated properly by submitting the `utils/test.sh` job
```bash
sbatch utils/test.sh
```
The output should look like this:
```
OpenMM Version: 7.7
Git Revision: 130124a3f9277b054ec40927360a6ad20c8f5fa6

There are 4 Platforms available:

1 Reference - Successfully computed forces
2 CPU - Successfully computed forces
3 CUDA - Successfully computed forces
4 OpenCL - Successfully computed forces

Median difference in forces between platforms:

Reference vs. CPU: 6.31347e-06
Reference vs. CUDA: 6.72976e-06
CPU vs. CUDA: 7.17819e-07
Reference vs. OpenCL: 6.74399e-06
CPU vs. OpenCL: 7.50779e-07
CUDA vs. OpenCL: 2.16192e-07

All differences are within tolerance.
```

### 3. Create AmberTools 23 environment (needed for ligand preparation in `build.sh`)

```bash
conda create -n ambertools23 -c conda-forge python=3.9 ambertools=23
conda install -c conda-forge openbabel
conda activate ambertools23
```
Ensure the correct installation of `AmberTools23` by running:
```bash
tleap -v
antechamber -h
```

### 4. VMD 1.9.4a57 Installation on Ubelix
This guide explains how to install VMD 1.9.4a57 on the Ubelix server (Linux RHEL 7+, 64-bit Intel/AMD) without CUDA or SIMD packages.

 1. Download VMD from https://www.ks.uiuc.edu/Research/vmd/
 2. Choose Linux version and copy to a folder on the server, e.g.:
 `LINUX_64 (RHEL 7+) OpenGL, CUDA, OptiX RTX, OSPRay (Linux (RHEL 7+) 64-bit Intel/AMD x86_64 SSE/AVX+ with CUDA 10, OptiX6.5 RTX, OSPRay185.opengl)`
 3. Transfer this file to the ubelix server in the location where you want VMD installed.
 4. Extract the package:
 ```bash
 gunzip vmd-1.9.4a57.bin.LINUXAMD64-CUDA102-OptiX650-OSPRay185.opengl.tar.gz
 tar -xvf vmd-1.9.4a57.bin.LINUXAMD64-CUDA102-OptiX650-OSPRay185.opengl.tar
 cd vmd-1.9.4a57
 ```

 5. Decide on the installation folder (`/storage/homefs/<usrname>/Tools`) and add it inside the `configure` file as shown below:
 
 ``` bash
 # Name of shell script used to start program; this is the name used by users
 $install_name = "vmd";

 # Directory where VMD startup script is installed, should be in users' paths.
 $install_bin_dir="/storage/homefs/<usrname>/Tools";
 
 # Directory where VMD files and executables are installed
 $install_library_dir="/storage/homefs/<usrname>/Tools/lib";
 ```
 
  7. Run the configure script without CUDA or SIMD:
  ```bash
  ./configure LINUXAMD64 TCL IMD COLVARS SILENT
  ```

  8. Build VMD: go to the `src` folder and finalize installation
  ```bash
  cd src
  make
  make install
  cd ../
  ```

  9. Add vmd to path, by adding this line `export PATH=“/storage/homefs/<usrname>/Tools/vmd:$PATH”` at the end of the `~/.bashrc` file. Then source ~/.bashrc.
  10. Check again the installation by typing just `vmd` from any folder. This should output something like:
  ```
  Info) VMD for LINUXAMD64, version 1.9.4a57 (April 27, 2022)
  Info) http://www.ks.uiuc.edu/Research/vmd/                         
  Info) Email questions and bug reports to vmd@ks.uiuc.edu           
  Info) Please include this reference in published work using VMD:   
  Info)    Humphrey, W., Dalke, A. and Schulten, K., `VMD - Visual   
  Info)    Molecular Dynamics', J. Molec. Graphics 1996, 14.1, 33-38.
  Info) -------------------------------------------------------------
  Info) Multithreading available, 128 CPUs.
  Info)   CPU features: SSE2 SSE4.1 AVX AVX2 FMA F16 HT 
  Info) Free system memory: 48GB (38%)
  Info) No CUDA accelerator devices available.
  Info) Dynamically loaded 3 plugins in directory:
  Info) /storage/homefs/cr22r967/vmd-1.9.4a57/plugins/LINUXAMD64/molfile
  vmd > 
  ```
       
> OpenMM v7.7 and AmberTools 23 are required. Using older versions may cause errors.

---

### 5. Faster engine (optional): NVIDIA cuPotential for `run.sh --faster`

A separate conda environment, `openmm_faster`: OpenMM 8.4 + the cuPotential OpenMM plugin, built
from source. The classic `openmm` environment is not touched and keeps working without `--faster`.

**Requirements**
* NVIDIA GPU with compute capability 7.5-12.0 (RTX 4090 = 8.9, A100 = 8.0, H100 = 9.0).
* GPU driver for **CUDA 13** (`nvidia-smi -q | grep -i "cuda"` on a GPU node). A CUDA 12 driver is
  too old: ask the cluster admins.
* conda, git and internet access (GitHub, conda-forge, PyPI). About 10 minutes and ~8 GB.

**Install** (no `module load CUDA` needed: conda provides the toolkit)

```bash
# on a GPU node (the script reads the driver and the GPU type):
bash utils/install_faster.sh

# or on a login node without GPUs: give the driver's CUDA version (from nvidia-smi on a GPU node)
bash utils/install_faster.sh --driver-cuda 13.3
```

Options: `--env <name>` (default `openmm_faster`), `--src <dir>` (cuPotential clone, default
`~/software/cuPotential`), `--arch "80;89;90"` (GPU types to compile for; default: the GPUs of the
machine, or 7.5-9.0 when there are none). Re-running the script updates an existing install.

The script, step by step:
1. clones cuPotential and downloads its prebuilt runtime from the GitHub release, and **reads which
   CUDA toolkit that runtime needs** (release 0.3.5: CUDA 13.4; NVIDIA's blog post still says 13.1);
2. creates the env with OpenMM 8.4, that CUDA toolkit and the compilers (conda-forge);
3. builds the plugin and installs it (`nv_openmm_plugin`);
4. **driver check**: if the driver supports an older CUDA minor than the toolkit (e.g. driver 13.3,
   toolkit 13.4), stock OpenMM fails with `CUDA_ERROR_UNSUPPORTED_PTX_VERSION`, because it compiles
   its kernels at run time (NVRTC) and the driver rejects newer code. The script then installs the
   NVRTC that matches the driver and preloads it on `conda activate openmm_faster` (hook in
   `$CONDA_PREFIX/etc/conda/activate.d/nvrtc_driver.sh`; delete it once the driver is updated, or
   just re-run the script). cuPotential's own kernels are precompiled and unaffected.

**Check compatibility** (on a GPU node; exit code 0 = ready)

```bash
conda activate openmm_faster
python utils/check_faster.py                                         # generated water box
python utils/check_faster.py step3_input.parm7 step3_input.rst7      # a built system
```

It checks OpenMM/plugin versions, driver vs runtime CUDA, the NVRTC version, the GPU type, then
compares the classic and faster engines **against a converged reference** (the faster engine must
be at least as accurate as the classic one) and runs 1000 NPT steps at 4 fs. Expected end:

```
PASS  energy error vs converged reference          faster 6.4e-06, classic 7.8e-06 ...
PASS  force error vs converged reference           faster 4.0e-04, classic 8.2e-04 ...
PASS  faster engine 1000 steps NPT, 4 fs           T = 302.6 K (target 303.15), finite = True
RESULT: READY for run.sh --faster
```

**Test the pipeline and measure the speed-up** on a system built with `build.sh` (GPU node):

```bash
bash utils/test_faster_pipeline.sh /path/to/built/system 0       # full run.sh --faster, shortened
python utils/benchmark_faster.py step3_input.parm7 step5_1.rst    # classic vs faster on this GPU
```

`test_faster_pipeline.sh` works on a temporary copy (the system folder is not changed) and ends
with `RESULT: run.sh --faster pipeline OK`. `benchmark_faster.py` is best given an equilibrated
state (`step4.*_equilibration.rst` or `step5_*.rst`).

**Run** (in `utils/job.sh`: `conda activate openmm_faster` and `bash utils/run.sh '0' --faster`)

```bash
conda activate openmm_faster
bash utils/run.sh 0 --faster
```

**Use one GPU per run** (`run.sh 0 --faster`, `--gres=gpu:1`): two GPUs work but were slower in
our test (996 vs 1,380 ns/day). On a new cluster, or after a driver update, re-run
`check_faster.py` before submitting.

---

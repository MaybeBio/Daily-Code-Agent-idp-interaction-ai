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

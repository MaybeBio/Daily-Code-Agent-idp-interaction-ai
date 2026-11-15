## Repository Structure

The repository is organized into three main directories, corresponding to different stages of the computational investigation of Synapsin I condensates using the CALVADOS 3.0 coarse-grained model.

### 1. [Effective_Models](https://github.com/leandrokrott/Synapsin-simulations/tree/main/Effective_Models)

Contains the LAMMPS input files, simulation configurations, and analysis scripts associated with the effective pair potentials derived from the two-protein Potential of Mean Force (PMF). Different analytical potential models (A and B) are implemented and used to investigate the phase behavior and structural properties of effective particles for Synapsin I systems through larger-scale simulations.

### 2. [PMF](https://github.com/leandrokrott/Synapsin-simulations/tree/main/PMF)

Contains the input files, molecular configurations, and simulation scripts used to perform two-protein simulations with CALVADOS 3.0. Umbrella sampling is employed to calculate the Potential of Mean Force (PMF) as a function of the center-of-mass separation between the ordered domain of two Synapsin I protein chains. The resulting PMF provides the reference effective interaction profile used to parameterize the analytical pair potentials.

### 3. [Synapsin_system](https://github.com/leandrokrott/Synapsin-simulations/tree/main/Synapsin_system)

Contains the input files, molecular configurations, and analysis scripts for collective simulations of a system consisting of 100 Synapsin I protein chains using CALVADOS 3.0. These simulations aim to investigate temperature-dependent phase behavior, liquid–liquid phase separation (LLPS), protein-rich domain formation, and the structural and thermodynamic properties of the system.

## Computational Workflow

The three directories are connected through the following computational workflow:

1. **PMF calculation:** Two-protein CALVADOS 3.0 simulations are performed using umbrella sampling to obtain the effective interaction profile.
2. **Effective potential parameterization:** The PMF is used to construct and parameterize analytical effective pair potentials.
3. **Collective simulations:** The resulting effective models are employed in larger-scale simulations to investigate the phase behavior of 100-protein Synapsin I systems.

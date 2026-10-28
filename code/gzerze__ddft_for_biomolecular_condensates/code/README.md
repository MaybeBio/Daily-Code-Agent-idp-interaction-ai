# DDFT for Biomolecular Condensates

This repository contains a research implementation of an adiabatic dynamical
density functional theory (DDFT) framework for the nonequilibrium evolution of
biomolecular condensates. The code combines a sequence-dependent
random-phase-approximation/Flory–Huggins (RPA+FH) bulk free energy with a
square-gradient interfacial term and evolves the protein density field using
the finite-volume package [FiPy](https://www.ctcms.nist.gov/fipy/).

The current release was developed primarily for two-dimensional, single-protein
calculations of droplet relaxation, interface evolution, and coalescence. It is
research code rather than a general-purpose DDFT package.

## Model and current scope

In reduced units, the implemented dynamics have the form

$$
\frac{\partial \phi}{\partial t}
= \nabla \cdot \left[\phi\nabla\mu\right],
\qquad
\mu = \mu_{\mathrm{RPA+FH}}(\phi)-\kappa\nabla^2\phi,
$$

where $\phi$ is the protein monomer density (volume fraction in the default
mapping), $\mu_{\mathrm{RPA+FH}}$ is the sequence-dependent bulk chemical
potential, and $\kappa$ is the square-gradient or influence parameter. The
main driver uses an adaptive time step and supports zero-flux and Robin-type
boundary conditions.

The released workflow should currently be used with `no_comps,1`. The example
input uses zero added salt. Experimental two-component routines are included,
but they are not yet connected to a complete two-field DDFT driver. Hydrodynamic
interactions and viscous dissipation are not included. Mapping the reduced time
to physical units therefore requires a diffusivity/mobility choice, and absolute
times depend on that mapping.

## Requirements

The code requires Python 3.10 or newer and the following packages:

- NumPy
- SciPy
- Matplotlib
- FiPy
- PETSc and `petsc4py`
- `mpi4py`
- MoviePy 1.x (`moviepy.editor` is imported by the current code)

One possible Conda installation is:

```bash
conda create -n ddft -c conda-forge python=3.10 numpy scipy matplotlib \
    fipy petsc4py mpi4py "moviepy<2"
conda activate ddft
```

Movie creation additionally requires FFmpeg and `gifsicle`. On an HPC system,
the MPI, PETSc, and Python modules should be adjusted to match the local
software stack.

## Quick start

From the repository root:

```bash
mkdir -p Results/example/states

python scripts/DDFT/RPA_ddft.py \
    --i scripts/DDFT/input_params.txt \
    --o Results/example/ \
    --s Results/example/
```

The trailing `/` on the output and state paths is required by the current path
handling. For cluster runs, the files in `scripts/bash/` provide SLURM examples,
but their module names, absolute paths, resources, and account-specific settings
must be edited before use.

If `state.pickle` is present in the state directory and its grid, salt,
temperature, and \(\kappa\) settings match the new input, the calculation resumes
from that checkpoint. Use a new output/state directory to begin an independent
run.

## Input parameters

Simulation settings are read from a comma-separated text file such as
`scripts/DDFT/input_params.txt`. The most important entries are:

| Parameter | Meaning |
| --- | --- |
| `dimension`, `nx`, `ny`, `dx` | Dimensionality and finite-volume grid; the supplied example is 2D. |
| `dt`, `dt_min`, `dt_max` | Initial, minimum, and maximum adaptive time steps. |
| `duration`, `total_steps`, `tolerance` | Stopping limits and nonlinear-sweep tolerance. |
| `seq_name1` | Protein sequence name defined in `Protein_RPA/seq_list.py`. |
| `phis` | Salt volume fraction; the supplied example uses zero. |
| `temperature_inverse` | Reduced inverse-temperature/electrostatic-coupling variable \(u\). |
| `enthalpy1`, `entropy1` | Parameters defining the short-range FH contribution, \(\chi(u)=e_hu+e_s\) in the default mapping. |
| `kappa` | Square-gradient/influence parameter controlling the interfacial penalty. |
| `idense`, `idil` | Initial dense- and dilute-phase protein densities. |
| `initial_density` | Name of an initial-condition method in `density_class.py`; the default creates two nearby circular droplets. |
| `boundary_type` | Boundary condition; the example uses `zero_flux`. |
| `phi_boundary`, `h_mul` | Reservoir density and transfer coefficient used by Robin-type boundaries. |
| `n_capture`, `state_save` | Frequencies for image generation and checkpoint storage. |

To add a protein, define its one-letter amino-acid sequence in the `polymers`
class in `scripts/DDFT/Protein_RPA/seq_list.py` and use that attribute name for
`seq_name1`. The current DDFT driver uses fixed integer charges for charged
residues; pH-dependent charging routines exist in `Protein_RPA` but are not
currently activated by the main input file.

## Output

A 2D calculation produces:

- `P_step_*.png`: two-dimensional protein-density fields;
- `D_step_*.png`: density profiles through the center of the domain;
- `M_step_*.png`: chemical-potential fields;
- `state.pickle`: the latest restart state; and
- `states/state-*.pickle`: saved checkpoint history.

`scripts/DDFT/movie.py` can convert the PNG series to GIF and MP4 files.

## Code organization

- `scripts/DDFT/RPA_ddft.py`: main time-dependent DDFT driver.
- `scripts/DDFT/RPA_cdft.py`: auxiliary cDFT driver.
- `scripts/DDFT/initialize.py`, `BC.py`, and `density_class.py`: initialization,
  boundary conditions, and initial density fields.
- `scripts/DDFT/chempot_table.py`: tabulation and interpolation of the bulk
  chemical potential and its derivative.
- `scripts/DDFT/Protein_RPA/`: sequence-dependent RPA+FH free-energy routines.
- `scripts/bash/`: example HPC/SLURM launch scripts.

## Code Provenance

The RPA-related routines in `scripts/DDFT/Protein_RPA/` were adapted from
[Protein_RPA](https://github.com/laphysique/Protein_RPA), developed by Yi-Hsuan
Lin and Hue Sun Chan under the MIT License. 
Y.-H. Lin, J. D. Forman-Kay, and H. S. Chan, Phys. Rev. Lett. 117, 178101 (2016)

The biomolecular condensate movie-generation workflow was adapted from the methodology provided in [Droplet_analysis](https://github.com/krishna-shrinivas/2020_Henninger_Oksuz_Shrinivas_RNA_feedback/tree/master/Droplet_analysis) , developed by Jonathan E. Henninger  under MIT License.
Henninger, J. E., Oksuz, O., Shrinivas, K., et al., Cell 184, 207–225 (2021)

If you use this DDFT implementation, please also cite the accompanying DDFT
publication. The complete citation will be added here once it is available.
Preprint: https://www.biorxiv.org/content/10.1101/2025.10.02.680159v1.abstract 

## License

The DDFT code in this repository is released under the MIT License; see
`LICENSE`. The adapted Protein_RPA components retain their upstream MIT license
and copyright notice in `scripts/DDFT/Protein_RPA/LICENSE`.


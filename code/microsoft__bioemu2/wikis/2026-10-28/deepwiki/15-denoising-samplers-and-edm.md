# Denoising Samplers and EDM

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [docs/sampling.md](docs/sampling.md)
- [src/bioemu2/biomoa/models/edm_schedule.py](src/bioemu2/biomoa/models/edm_schedule.py)
- [src/bioemu2/biomoa/modules/diffusion.py](src/bioemu2/biomoa/modules/diffusion.py)

</details>



## Purpose and Scope

This page covers the technical implementation of the denoising pipeline in BioEmu2, focusing on the `src/bioemu2/denoising` subpackage and its integration with the Elucidated Diffusion Model (EDM) framework. It details the noise schedules, sampling algorithms, auxiliary utilities (`augmentation`, `rigid_align`, `cp_sync`), and the operational differences between structure prediction and equilibrium sampling modes.

Sources: `src/bioemu2/biomoa/models/edm_schedule.py:1-87](), [docs/sampling.md:1-120]()`

---

## 1. Elucidated Diffusion Model (EDM) Schedule and Scaling

BioEmu2 implements the score-based diffusion formulation introduced by Karras et al., adapted for biomolecular atomic coordinates. The core scaling functions and noise conditioning are defined in `src/bioemu2/biomoa/models/edm_schedule.py` via the `EDMSchedule` class.

### Key Scaling Coefficients

The denoised coordinates ($\mathbf{x}_{\text{denoised}}$) are reconstructed from noisy inputs ($\mathbf{x}_{\text{noisy}}$) and the network prediction ($\mathbf{f}_{\theta}$) using preconditioning coefficients that depend on the noise level $\sigma$:

$$\mathbf{x}_{\text{denoised}} = c_{\text{skip}}(\sigma) \mathbf{x}_{\text{noisy}} + c_{\text{out}}(\sigma) \mathbf{f}_{\theta}(c_{\text{in}}(\sigma) \mathbf{x}_{\text{noisy}}, \sigma)$$

| Coefficient | Function | Mathematical Definition | Role |
| --- | --- | --- | --- |
| `c_skip` | `EDMSchedule.c_skip(sigma)` | $\frac{\sigma_{\text{data}}^2}{\sigma^2 + \sigma_{\text{data}}^2}$ | Retains a proportion of the noisy input signal. |
| `c_out` | `EDMSchedule.c_out(sigma)` | $\frac{\sigma \sigma_{\text{data}}}{\sqrt{\sigma_{\text{data}}^2 + \sigma^2}}$ | Scales the neural network's coordinate update. |
| `c_in` | `EDMSchedule.c_in(sigma)` | $\frac{1}{\sqrt{\sigma^2 + \sigma_{\text{data}}^2}}$ | Normalizes input atomic coordinates to unit variance. |
| `c_noise`| `EDMSchedule.c_noise(sigma)` | $\frac{1}{4} \ln\left(\frac{\sigma}{\sigma_{\text{data}}}\right)$ | Computes noise level conditioning for time/sigma embeddings. |

The default data standard deviation is set to $\sigma_{\text{data}} = 16.0$ Å, with noise bounds ranging from $\sigma_{\min} = 0.0004$ to $\sigma_{\max} = 160.0$.

```mermaid
title "EDMSchedule Data Flow and Coordinate Transformation"

graph TD
    InputCoords["Input Coordinates x_noisy"] --> Cin["EDMSchedule.c_in(sigma)"]
    InputCoords --> Cskip["EDMSchedule.c_skip(sigma)"]
    Sigma["Noise Level sigma"] --> Cin
    Sigma --> Cskip
    Sigma --> Cout["EDMSchedule.c_out(sigma)"]
    Sigma --> Cnoise["EDMSchedule.c_noise(sigma)"]
    
    Cin --> Model["DualStreamDiffusionModule.forward()"]
    Cnoise --> Model
    Model --> NetOut["Network Output f_theta"]
    
    NetOut --> Cout
    Cout --> Combine["Linear Combination"]
    Cskip --> Combine
    InputCoords --> Combine
    
    Combine --> DenoisedCoords["Output Coordinates x_denoised"]

    classDef default fill:#fff,stroke:#333,stroke-width:1px;
```

Sources: `[src/bioemu2/biomoa/models/edm_schedule.py:1-86](), [src/bioemu2/biomoa/modules/diffusion.py:1-135]()`

---

## 2. Samplers and Time Discretization

The sampling engine orchestrates the reverse SDE/ODE integration from pure noise ($\sigma_{\max}$) down to data scale ($\sigma_{\min}$). 

### Sampler Execution Pipeline

Samplers in `src/bioemu2/denoising/samplers.py` manage discretization schedules, stochasticity parameters (such as stochastic Langevin-like noise additions), and step scaling heuristics.

```mermaid
title "Sampler Execution Flow in Code Entity Space"

graph TD
    Session["EmuSession.sample()"] --> SamplerEntry["samplers.sample_loop()"]
    SamplerEntry --> ScheduleGen["Generate Sigma Schedule (t__steps)"]
    
    ScheduleGen --> Loop{For each sigma_i}
    Loop --> PredictScore["DualStreamDiffusionModule()"]
    PredictScore --> ApplyEDM["EDMSchedule (c_skip, c_out, c_in)"]
    ApplyEDM --> EulerStep["Discretization Step (Euler / Heun)"]
    EulerStep --> CheckCP["cp_sync.synchronize_ranks()"]
    CheckCP --> Align["rigid_align.kabsch_align()"]
    Align --> NextSigma{More steps?}
    NextSigma -- Yes --> Loop
    NextSigma -- No --> ReturnSamples["Return Tensor of Conformations"]

    classDef default fill:#fff,stroke:#333,stroke-width:1px;
```

Sources: `[src/bioemu2/biomoa/modules/diffusion.py:1-160](), [docs/sampling.md:72-120]()`

---

## 3. Equilibrium vs. Structure Modes

BioEmu2 distinguishes between two distinct operational modes tailored for different biophysical applications. These modes modify parameters like step scaling (`step_scale`) and sampling noise during execution.

| Mode | CLI Subcommand | Session Configuration | Target Objective | Key Hyperparameters |
| --- | --- | --- | --- | --- |
| **Structure Prediction** | `bioemu2 predict-structure` | `EmuSession(..., mode="structure")` | High-accuracy individual folded structure | Heuristic sharpening with `step_scale = 1.5` |
| **Equilibrium Sampling** | `bioemu2 sample-equilibrium`| `EmuSession(..., mode="equilibrium")`| Broad thermodynamic conformational ensemble | `step_scale = 1.0`, unsharpened trajectory |

### Mode Selection via Python API

When initializing an `EmuSession`, the mode string determines internal diffusion step multipliers and output configurations:

```python
from bioemu2 import EmuSession, System, load_model

system = System.from_pdb_id("1AKE")
model = load_model()

# Structure prediction mode (applies heuristic step scaling)
struct_session = EmuSession(system, model=model, mode="structure")

# Equilibrium sampling mode (standard thermodynamic path)
eq_session = EmuSession(system, model=model, mode="equilibrium")
```

Sources: `[docs/sampling.md:107-120]()`

---

## 4. Supporting Denoiser Utilities

The denoising subsystem relies on several specialized modules for coordinate manipulation, symmetry, and parallel execution:

* **`augmentation`**: Applies random spatial rotations and translations during training or data preprocessing to ensure SO(3) invariance.
* **`rigid_align`**: Implements root-mean-square deviation (RMSD) minimization routines (such as the Kabsch algorithm) to superimpose generated frames onto reference topologies or input templates.
* **`cp_sync`**: Manages tensor synchronization across multiple GPUs when executing context-parallel sampling runs for large macromolecules.

Sources: `[docs/sampling.md:90-96]()`

---
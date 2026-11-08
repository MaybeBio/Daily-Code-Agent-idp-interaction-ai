# Experimental Guidance Potentials

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [idpforge/utils/potential.py](idpforge/utils/potential.py)
- [sample_idp.py](sample_idp.py)
- [sample_ldr.py](sample_ldr.py)

</details>



IDPForge incorporates a gradient-based guidance system that allows experimental data to steer the reverse diffusion process. By defining differentiable "potentials" based on physical observables (Radius of Gyration, FRET, PRE, J-Couplings), the model can generate ensembles that are statistically consistent with experimental measurements while remaining within the learned structural manifold of intrinsically disordered proteins (IDPs).

## Potential Guidance Framework

The guidance system is implemented in `idpforge/utils/potential.py`. It operates by calculating the gradient of a potential function with respect to the atomic coordinates during each step of the reverse diffusion process. These gradients are then used to shift the predicted denoised structure toward states that better satisfy the experimental constraints.

### The Potential Base Class
All guidance terms inherit from the `Potential` interface [idpforge/utils/potential.py:15-18](). This class defines the contract for computing potentials and their corresponding gradients.

*   **`compute(xyz)`**: Calculates a scalar potential value from the current coordinates `[B, L, 5, 3]`. The system is designed to **maximize** this value [idpforge/utils/potential.py:19-31]().
*   **`get_potential_gradients(xyz)`**: Uses PyTorch's `autograd` to compute the gradient of the potential with respect to the input coordinates. It specifically handles zeroing out gradients for atoms not involved in the potential calculation and ensures that gradients are propagated to the backbone atoms (N, Ca, C) [idpforge/utils/potential.py:33-53]().

### Guidance Logic Data Flow

The following diagram illustrates how the `Potential` classes interface with the sampling loop in `sample_idp.py`.

**Potential Integration Overview**
```mermaid
graph TD
    subgraph "Inference Loop (sample_idp.py)"
        START["Start Reverse Diffusion"] --> INIT["Initialize xt, tor_list"]
        INIT --> MODEL_SAMPLE["model.sample()"]
    end

    subgraph "Guidance Application (idpforge/model.py)"
        MODEL_SAMPLE --> POT_CHECK{"Potential Configured?"}
        POT_CHECK -- Yes --> CALC_GRAD["Potential.get_potential_gradients()"]
        CALC_GRAD --> APPLY_STEP["Apply Gradient Step to Denoised Coord"]
        APPLY_STEP --> NEXT_T["Next Diffusion Step"]
    end

    subgraph "Potential Implementations (idpforge/utils/potential.py)"
        CALC_GRAD -.-> RoG["RoG (Radius of Gyration)"]
        CALC_GRAD -.-> Contact["Contact (PRE/NOE)"]
        CALC_GRAD -.-> Efret["Efret (smFRET)"]
        CALC_GRAD -.-> JCoup["JCoup (J-Couplings)"]
    end
```
**Sources:** [sample_idp.py:142-156](), [idpforge/utils/potential.py:15-53]()

## Available Potential Types

IDPForge provides several specialized potential classes to handle different types of experimental data. These are mapped in the `Potentials` dispatch table [idpforge/utils/potential.py:162-168]().

| Potential Key | Class | Description | Implementation Detail |
| :--- | :--- | :--- | :--- |
| `rg` | `RoG` | Radius of Gyration | Minimizes squared error between ensemble mean Rg and target [idpforge/utils/potential.py:56-67](). |
| `contact` | `Contact` | Distance Restraints (PRE/NOE) | Uses a differentiable coordination-style function to enforce distance bounds [idpforge/utils/potential.py:70-98](). |
| `fret` | `Efret` | smFRET Efficiency | Calculates FRET efficiency using $E = 1 / (1 + (r/R_0)^6)$ [idpforge/utils/potential.py:101-123](). |
| `jcoup` | `JCoup` | Scalar J-Couplings | Uses Karplus equations to relate $\phi$ dihedrals to J-values [idpforge/utils/potential.py:126-148](). |
| `multiple` | `Multiple` | Weighted Sum | Combines multiple potential terms into a single objective [idpforge/utils/potential.py:151-159](). |

### Implementation Details

#### Radius of Gyration (`RoG`)
Calculates the distance of Ca atoms from the centroid. The potential returned is the negative squared error, encouraging the model to match the `target` Rg [idpforge/utils/potential.py:61-67]().

#### Contact Potential (`Contact`)
Used for PRE (Paramagnetic Resonance Enhancement) or NOE (Nuclear Overhauser Effect) data. It calculates a differentiable distance gram (`dgram`) and applies a loss if distances fall outside the specified `contact_bounds` [idpforge/utils/potential.py:93-95](). It includes an `exp_mask_p` parameter to randomly mask experimental constraints during sampling, which helps maintain ensemble diversity [idpforge/utils/potential.py:96-97]().

#### J-Couplings (`JCoup`)
Calculates the $\phi$ dihedral angle using `get_dih` [idpforge/utils/potential.py:143]() and applies a standard Karplus equation:
$J(\phi) = 6.51 \cos^2(\phi - 60^\circ) - 1.76 \cos(\phi - 60^\circ) + 1.6$ [idpforge/utils/potential.py:144-145]().

**Sources:** [idpforge/utils/potential.py:56-168]()

## Configuration and Application

In `sample_idp.py`, potentials are initialized based on the `sample.yml` configuration.

**Potential Initialization Logic**
```mermaid
graph TD
    CFG["sample.yml Config"] --> PARSE["Parse settings['potential_cfg']"]
    PARSE --> PRE_CHECK{"Key == 'pre'?"}
    PRE_CHECK -- Yes --> MAP["get_contact_map()"]
    MAP --> INST["Instantiate Contact()"]
    
    PARSE --> RG_CHECK{"Key == 'rg'?"}
    RG_CHECK -- Yes --> INST_RG["Instantiate RoG()"]
    
    INST --> MULTI["Wrap in Multiple()"]
    INST_RG --> MULTI
    
    MULTI --> SAMPLE_CALL["model.sample(..., potential_cfgs=potential_cfg)"]
```
**Sources:** [sample_idp.py:64-87](), [idpforge/utils/potential.py:151-159]()

### Key Parameters in Configuration
*   **`timescale`**: Controls the diffusion steps during which guidance is active [sample_idp.py:66]().
*   **`grad_clip`**: Prevents unstable updates by clipping the magnitude of the potential gradient [sample_idp.py:67]().
*   **`weights`**: A dictionary mapping potential types to their relative importance in the `Multiple` potential sum [idpforge/utils/potential.py:155]().
*   **`exp_mask_p`**: Probability of keeping an experimental constraint in a given step; used to prevent over-fitting to specific restraints [idpforge/utils/potential.py:81]().

### Gradient Step Application
The gradients obtained from `Potential.get_potential_gradients()` are applied to the predicted coordinates. For distance-based potentials (RoG, Contact, Efret), the system ensures gradients are propagated to the Ca atoms [idpforge/utils/potential.py:46-48](). For J-Couplings, gradients are calculated for the N-Ca-C backbone atoms to influence the torsion angles [idpforge/utils/potential.py:36-37]().

**Sources:** [idpforge/utils/potential.py:33-53](), [sample_idp.py:64-87]()

---
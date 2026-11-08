# Back-Calculation of Experimental Observables

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [scoring/calculator.py](scoring/calculator.py)
- [scoring/parser.py](scoring/parser.py)

</details>



This page documents the technical implementation of back-calculators used to bridge the gap between 3D protein conformers and experimental measurements. These utilities, primarily located in `scoring/calculator.py`, allow IDPForge to validate generated ensembles against various biophysical data types, including J-couplings, NOE, PRE, smFRET, chemical shifts, and Radius of Gyration (Rg).

### Overview and Data Flow

The back-calculation process transforms atomic coordinates from PDB files into predicted observables. This is a critical component of the **X-EISD** (Experimental Enhanced Integrative Structural Dissemination) scoring framework.

The system uses a dispatch table `BACK_Calculators` to map experimental data types to their respective calculation functions [scoring/calculator.py:207-214]().

**Sources:**
- [scoring/calculator.py:1-13]()
- [scoring/calculator.py:207-214]()

---

### Implementation of Back-Calculators

#### 1. Scalar Couplings (J-couplings)
The `jc_backcalc` function predicts 3-bond $^3J_{HNH\alpha}$ couplings based on the protein backbone dihedral angle $\phi$. It utilizes a Karplus-like relationship where the predicted value is proportional to the cosine of the angle [scoring/calculator.py:108-123]().

*   **Logic**: Calculates $\phi$ for each residue and applies a phase-shifted cosine transformation: $\cos(\phi - 60^\circ)$ [scoring/calculator.py:122]().
*   **Parameters**: Uses Karplus coefficients (A, B, C) defined in the `Stack` initialization for J-couplings [scoring/parser.py:76-79]().

#### 2. Distance-Based Observables (NOE, PRE, smFRET)
These calculators rely on inter-atomic distances, often requiring specific handling for hydrogen atoms or paramagnetic labels.

| Observable | Function | Implementation Detail |
| :--- | :--- | :--- |
| **smFRET** | `fret_backcalc` | Calculates efficiency $E = 1 / (1 + (d/R_0)^6)$ between CA atoms of specified residues [scoring/calculator.py:126-140](). |
| **PRE** | `pre_backcalc` | Calculates $1/r^6$ distances between a spin-label residue and target protons. Supports virtual atom substitution [scoring/calculator.py:165-182](). |
| **NOE** | `noe_backcalc` | Similar to PRE, but typically involves proton-proton distances within a specific cutoff [scoring/calculator.py:148-162](). |

#### 3. Chemical Shifts (CSpred/UCBShift)
Chemical shifts are calculated via `run_cspred`, which acts as a wrapper for the external **UCBShift** (CSpred) tool.
*   **Process**: It generates a temporary PDB, executes a subprocess running `CSpred.py`, and parses the resulting text output [scoring/calculator.py:143-145]().
*   **Environment**: Requires `CSPRED_PATH` and `CSPRED_PYTHON` environment variables to be set [scoring/calculator.py:14-17]().

#### 4. Radius of Gyration (Rg)
The `calc_rg` function computes the mass-weighted radius of gyration.
*   **Mass Mapping**: Uses a standard `atomic_mass` dictionary (C:12, O:16, N:14, S:32, H:1) [scoring/calculator.py:94]().
*   **Formula**: $R_g = \sqrt{\frac{\sum m_i (r_i - r_{cm})^2}{\sum m_i}}$ where $r_{cm}$ is the center of mass.

**Sources:**
- [scoring/calculator.py:94-95]()
- [scoring/calculator.py:108-140]()
- [scoring/calculator.py:165-182]()
- [scoring/parser.py:76-80]()

---

### Hydrogen Atom Substitution Logic

Since many experimental observables (NOE, PRE, Chemical Shifts) involve hydrogen atoms that may not be present in standard PDB outputs, the system implements a mapping logic to identify the correct heavy-atom "parent" or to handle methyl/methylene group equivalence.

#### Mapping Dictionaries
*   **`hydrogen_abbrev`**: Maps shorthand hydrogen names (e.g., "HB" for Alanine) to specific PDB atom names (e.g., "HB1", "HB2", "HB3") [scoring/calculator.py:44-71]().
*   **`heavy_atom_hydrogen`**: Maps hydrogen names back to their bonded heavy atoms (e.g., Alanine "HA" is bonded to "CA") [scoring/calculator.py:72-93]().

#### Virtual Atom Logic
In functions like `pre_backcalc`, if a hydrogen atom is missing, the calculator can fallback to the coordinates of the parent heavy atom defined in `heavy_atom_hydrogen` to estimate the distance [scoring/calculator.py:178-180]().

**Sources:**
- [scoring/calculator.py:44-93]()
- [scoring/calculator.py:175-181]()

---

### System Entity Map: Back-Calculation Dispatch

The following diagram illustrates how the high-level observable types requested by the user are dispatched to specific code entities and internal data structures.

**Calculator Dispatch Architecture**
```mermaid
graph TD
    subgraph "Natural Language Space"
        J["J-Couplings"]
        CS["Chemical Shifts"]
        FRET["smFRET"]
        Dists["Distances (NOE/PRE)"]
    end

    subgraph "Code Entity Space: scoring/calculator.py"
        Dispatch["BACK_Calculators (dict)"]
        JC_Fn["jc_backcalc()"]
        CS_Fn["run_cspred()"]
        FRET_Fn["fret_backcalc()"]
        PRE_Fn["pre_backcalc()"]
        NOE_Fn["noe_backcalc()"]
        
        Dispatch -->|'jc'| JC_Fn
        Dispatch -->|'cs'| CS_Fn
        Dispatch -->|'fret'| FRET_Fn
        Dispatch -->|'pre'| PRE_Fn
        Dispatch -->|'noe'| NOE_Fn
    end

    subgraph "External Dependencies"
        UCB["UCBShift (Subprocess)"]
        DSSP["mkdssp (Binary)"]
        
        CS_Fn --> UCB
        CS_Fn --> DSSP
    end

    J --- Dispatch
    CS --- Dispatch
    FRET --- Dispatch
    Dists --- Dispatch
```
**Sources:**
- [scoring/calculator.py:14-20]()
- [scoring/calculator.py:207-214]()

---

### Data Registry and Uncertainties

Experimental data is managed via the `EXP_DATA_LIB`, which lazily scans the directory specified by `IDPFORGE_EXP_DATA` [scoring/parser.py:24-72]().

#### Back-Calculation Errors (`BC_ERRORS`)
The system maintains a registry of inherent uncertainties for each back-calculation method, which are used during X-EISD scoring to weight the likelihood of the ensemble [scoring/parser.py:19-22]().

| Data Type | Uncertainty ($\sigma$) | Source/Note |
| :--- | :--- | :--- |
| `pre` | 0.001 | Fixed |
| `noe` | 0.01 | Fixed |
| `fret` | 0.0074 | Fixed |
| `cs (CA)` | 0.97 | UCBShift Error |
| `cs (H)` | 0.38 | UCBShift Error |
| `jc` | Variable | A: 0.14, B: 0.03, C: 0.08 (sqrt) |

**Sources:**
- [scoring/parser.py:19-22]()
- [scoring/parser.py:75-80]()

---

### Ensemble Processing Workflow

The back-calculation typically follows a multi-step workflow when processing large ensembles:

**Ensemble Processing Flow**
```mermaid
flowchart LR
    PDBs["Trajectory (traj.pdb)"] --> Split["split_pdbs()"]
    Split --> Files["Individual .pdb files"]
    Files --> Calc["BACK_Calculators"]
    
    subgraph "Back-Calculator Logic"
        Calc --> Parser["Bio.PDB.PDBParser"]
        Parser --> Loop["Iterate Models & Residues"]
        Loop --> Math["Geometric/External Calculation"]
    end
    
    Math --> Results["Results Matrix (N_models x N_obs)"]
```

**Sources:**
- [scoring/calculator.py:97-106]()
- [scoring/calculator.py:115-123]()
- [scoring/calculator.py:135-140]()

---
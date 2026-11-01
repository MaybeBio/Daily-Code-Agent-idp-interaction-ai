# Data Conversion CLI: starling-converter Tools

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [starling/scripts/starling_converter.py](starling/scripts/starling_converter.py)
- [starling/structure/ensemble.py](starling/structure/ensemble.py)
- [starling/utilities.py](starling/utilities.py)

</details>



The `starling-converter` suite provides a comprehensive set of command-line utilities for transforming STARLING ensemble data between various formats. These tools follow a consistent **Load-Transform-Save** pattern, leveraging the `Ensemble` class for data manipulation and `soursop` for structural trajectory handling.

## Overview of Entrypoints

All conversion tools are defined in `starling/scripts/starling_converter.py` and are typically exposed as CLI commands after installation.

| Command | Input | Output | Description |
| :--- | :--- | :--- | :--- |
| `starling2xtc` | `.starling` | `.xtc` + `.pdb` | Generates 3D coordinates and saves as GROMACS trajectory. |
| `starling2pdb` | `.starling` | `.pdb` | Generates 3D coordinates and saves as a multi-model PDB. |
| `starling2numpy` | `.starling` | `.npy` | Extracts distance maps as a raw NumPy array. |
| `starling2sequence` | `.starling` | `stdout` | Prints the amino acid sequence to the terminal. |
| `starling2info` | `.starling` | `stdout` | Prints metadata (version, date, weights used, mean Rg). |
| `starling2starling` | `.starling` | `.starling` | Performs error checking and filtering on an existing ensemble. |
| `numpy2starling` | `.npy` + seq | `.starling` | Wraps raw distance maps into a STARLING ensemble object. |
| `xtc2starling` | `.xtc` + `.pdb` | `.starling` | Converts a coordinate trajectory back into distance maps. |

**Sources:**
* `starling/scripts/starling_converter.py:14-306`()

---

## Load-Transform-Save Pattern

The converters follow a unified execution flow. The system first parses the output path using `parse_output_path` `[starling/utilities.py:89-123]()`, then loads the data into an `Ensemble` object `[starling/structure/ensemble.py:42-130]()`, performs requested transformations (like 3D reconstruction), and finally serializes the result.

### Code Entity Space: Conversion Logic Flow

The following diagram illustrates how the CLI functions map to the underlying `Ensemble` methods and utility functions.

```mermaid
graph TD
    subgraph "CLI Entrypoints (starling_converter.py)"
        S2X["starling2xtc()"]
        S2P["starling2pdb()"]
        N2S["numpy2starling()"]
        X2S["xtc2starling()"]
    end

    subgraph "Core Entities (ensemble.py)"
        ENS["class Ensemble"]
        LOAD["load_ensemble()"]
        BUILD["Ensemble.build_ensemble_trajectory()"]
        SAVE_T["Ensemble.save_trajectory()"]
        SAVE_E["Ensemble.save()"]
    end

    subgraph "Utilities (utilities.py)"
        PARSE["parse_output_path()"]
        SYM["symmetrize_distance_maps()"]
    end

    S2X --> PARSE
    S2X --> LOAD
    LOAD --> ENS
    S2X --> BUILD
    S2X --> SAVE_T

    N2S --> SYM
    N2S --> ENS
    N2S --> SAVE_E
```

**Sources:**
* `starling/scripts/starling_converter.py:14-63`()
* `starling/scripts/starling_converter.py:228-268`()
* `starling/structure/ensemble.py:536-574`()
* `starling/utilities.py:89-123`()

---

## Structural Reconstruction (starling2xtc / starling2pdb)

When converting to coordinate-based formats (`.xtc` or `.pdb`), the tools check if 3D coordinates are already present. If not, they trigger the MDS-based reconstruction via `build_ensemble_trajectory` `[starling/structure/ensemble.py:536-574]()`.

### Error Checking Workflow
Both `starling2xtc` and `starling2pdb` support the `--remove-errors` flag. This invokes `check_for_errors_trajectory` `[starling/structure/ensemble.py:246-302]()`, which scans the generated coordinates for physical impossibilities (e.g., $C\alpha-C\alpha$ distances $< 2.0$ Å or $> 5.0$ Å for adjacent residues).

**Key Implementation Details:**
1.  **Device Selection:** Users can specify `--device` (cpu, cuda, mps) to accelerate the MDS reconstruction process `[starling/scripts/starling_converter.py:33-37]()`.
2.  **Topology:** A $C\alpha$-only topology is automatically generated using `create_ca_topology_from_coords` `[starling/structure/coordinates.py:36-39]()`.
3.  **Persistence:** Data is saved using `Ensemble.save_trajectory`, which utilizes `soursop.sstrajectory.SSTrajectory` for writing binary formats `[starling/structure/ensemble.py:586-639]()`.

**Sources:**
* `starling/scripts/starling_converter.py:14-112`()
* `starling/structure/ensemble.py:246-302`()
* `starling/structure/ensemble.py:536-574`()

---

## Data Flow: numpy2starling and xtc2starling

These tools allow external data (from MD simulations or other models) to be imported into the STARLING ecosystem for analysis or search indexing.

### numpy2starling
This tool takes a `.npy` file containing distance maps and a sequence string.
-   It ensures maps are symmetrized using `symmetrize_distance_maps` `[starling/utilities.py:10]()`.
-   It validates that the sequence length matches the matrix dimensions `[starling/structure/ensemble.py:157-165]()`.
-   It produces a compressed `.starling` file (using `pickle` and `lzma`) `[starling/structure/ensemble.py:511-534]()`.

### xtc2starling
This tool extracts distance maps from a 3D trajectory.
-   It uses `soursop.sstrajectory.SSTrajectory` to load the `.xtc` and `.pdb` topology `[starling/scripts/starling_converter.py:288-290]()`.
-   It calculates pairwise $C\alpha$ distances for every frame.
-   The resulting distance maps are stored in a new `Ensemble` object.

### System Mapping: Data Ingestion

```mermaid
graph LR
    subgraph "External Data"
        NPY[".npy file"]
        XTC[".xtc / .pdb"]
    end

    subgraph "Ingestion Logic"
        N2S_FN["numpy2starling()"]
        X2S_FN["xtc2starling()"]
        SYM_UTIL["symmetrize_distance_maps()"]
    end

    subgraph "Internal State"
        ENS_OBJ["Ensemble Instance"]
        META["Metadata (Date, Version)"]
    end

    NPY --> N2S_FN
    XTC --> X2S_FN
    N2S_FN --> SYM_UTIL
    SYM_UTIL --> ENS_OBJ
    X2S_FN --> ENS_OBJ
    ENS_OBJ --> META
    META --> SAVE[".starling file"]
```

**Sources:**
* `starling/scripts/starling_converter.py:228-306`()
* `starling/structure/ensemble.py:77-130`()
* `starling/utilities.py:307-331`()

---

## Error Correction: starling2starling

The `starling2starling` utility is specifically designed for cleaning up "badly behaved" sequences. It performs a pass over the distance maps to identify frames that cannot be physically reconstructed into 3D space.

**Workflow:**
1.  Load the source `.starling` file.
2.  Call `Ensemble.check_for_errors(remove_errors=True)` `[starling/structure/ensemble.py:181-244]()`.
3.  This method checks for:
    -   Adjacent $C\alpha$ distances outside the $3.8 \pm 1.2$ Å range.
    -   Non-adjacent distances that are physically impossible (steric clashes).
4.  Save the pruned ensemble to a new file.

**Sources:**
* `starling/scripts/starling_converter.py:208-225`()
* `starling/structure/ensemble.py:181-244`()

---
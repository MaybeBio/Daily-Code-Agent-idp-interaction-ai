# Template Preparation Utilities

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [mk_flex_template.py](mk_flex_template.py)
- [mk_ldr_template.py](mk_ldr_template.py)
- [slice_pdb.py](slice_pdb.py)

</details>



The template preparation utilities provide the necessary pre-processing to generate `.npz` files used by the inference engine (`sample_ldr.py`). These utilities handle the extraction of structural data from PDB files, size-capping for large proteins, and the initialization of disordered region coordinates via geometric sampling or random domain placement.

## Core Utilities Overview

The pipeline utilizes two primary scripts to prepare templates based on the protein's architecture:
1.  **`mk_ldr_template.py`**: Used for Intrinsically Disordered Regions (IDRs) that are flanking or embedded within a static folded context (e.g., N/C-terminal tails or loops).
2.  **`mk_flex_template.py`**: Used for flexible multi-domain proteins where a linker connects two folded domains that can rotate and translate relative to one another.

### System Mapping: Utilities to Code Entities

The following diagram maps the logical preparation steps to the specific functions and classes in the codebase.

**Template Preparation Data Flow**
```mermaid
graph TD
    PDB["Input PDB File"] --> PP["process_pdb()"]
    PP --> GCA["get_chi_angles()"]
    PP --> DSSP["md.compute_dssp()"]
    
    subgraph "mk_ldr_template.py Logic"
        LDR_S["sample_fix_distance()"]
        LDR_T["_truncate_to_size()"]
    end

    subgraph "mk_flex_template.py Logic"
        FLEX_S["sample_fold_orientation()"]
        FLEX_T["_truncate_linker_domains()"]
    end

    GCA --> NPZ[".npz Template File"]
    DSSP --> NPZ
    LDR_T --> NPZ
    FLEX_T --> NPZ
```
**Sources:** `mk_ldr_template.py` [132-143](), `mk_flex_template.py` [181-195]()

---

## Folded Domain Extraction and Feature Encoding

Both utilities begin by extracting atomic coordinates and sequence information using `process_pdb` [mk_ldr_template.py:132](). The scripts then calculate essential structural features required by the IDPForge model:

*   **Torsion Angles**: Side-chain $\chi$ angles are calculated via `get_chi_angles` and converted into $(\sin, \cos)$ pairs [mk_ldr_template.py:133-134]().
*   **Secondary Structure (SS)**: SS is derived using MDTraj's DSSP implementation [mk_ldr_template.py:137](). For residues not classified as Helix (H) or Sheet (E), the `assign_rama` function classifies them based on $(\phi, \psi)$ angles into Ramachandran-based categories (e.g., 'P' for Polyproline II, 'B' for Beta-like) [mk_ldr_template.py:142-143]().

**Sources:** `mk_ldr_template.py` [132-143](), `mk_flex_template.py` [181-195]()

---

## Static IDR Initialization (`mk_ldr_template.py`)

For proteins with static folded domains, the script must provide "seed" coordinates for the disordered regions to initialize the diffusion process.

### Seed-Distance Sampling
The function `sample_fix_distance` generates initial CA coordinates for the first residue of an IDR tail or the midpoint of an IDR loop.
*   It samples points on a sphere with a radius between `SEED_FLOOR` (6.46 Å) and `SEED_CEILING` (9.12 Å) [mk_ldr_template.py:23-24]().
*   To prevent steric clashes in the initial template, it drops any sampled points within `min_fold_dist` of existing folded CA atoms [mk_ldr_template.py:40-56]().

### Size-Capping and Truncation
For proteins exceeding the model's capacity (e.g., `max_residues`), the `_truncate_to_size` function trims the folded domains while preserving the junction interface.
*   **Budget Allocation**: It ensures the total system (IDR + kept fold) fits within `max_residues`.
*   **Junction Preservation**: It maintains a minimum of 10 residues (default) per side of the IDR to ensure the model captures the junction geometry [mk_ldr_template.py:74-81]().
*   **Graft Specification**: The script outputs a `graft_offset` and `graft_fold_range` in the `.npz` file, which is used by the AlphaFlex Step 4 stitching logic to reassemble the full-length protein [mk_ldr_template.py:111-116]().

**Sources:** `mk_ldr_template.py` [40-56](), [59-125]()

---

## Flexible Domain Placement (`mk_flex_template.py`)

When modeling a linker between two domains, the domains are treated as rigid bodies that are randomly oriented to create a diverse starting ensemble.

### Random SO(3) Transformations
The script applies random rotations and translations to the second domain (`fold2`) relative to the first (`fold1`).
*   **Rotation**: Uses `scipy.spatial.transform.Rotation.random()` to sample uniformly from the SO(3) group [mk_flex_template.py:29-31]().
*   **Translation**: The distance $d$ between domains is sampled from a normal distribution based on Flory’s polymer theory: $d = 5.51 \times L^{0.588}$, where $L$ is the linker length [mk_flex_template.py:98-99]().

### Collision Detection
The function `random_nonclashing_transform_scaled` uses a `KDTree` to ensure that the randomly placed domains do not clash (minimum distance > 3.8 Å) [mk_flex_template.py:48-84](). If a non-clashing pose cannot be found after `max_attempts`, it falls back to a potentially clashing pose to allow the diffusion process to resolve the geometry [mk_flex_template.py:127-139]().

**Sources:** `mk_flex_template.py` [29-31](), [48-84](), [98-102]()

---

## The .npz Template Format

The output of these utilities is a compressed NumPy file containing the following keys:

| Key | Description |
| :--- | :--- |
| `coord` | `(L, 37, 3)` array of atomic coordinates (all-atom). |
| `torsion` | `(L, 7, 2)` array of $(\sin, \cos)$ for $\chi$ angles. |
| `mask` | `(L,)` boolean mask (True for folded, False for IDR). |
| `seq` | One-letter amino acid sequence string. |
| `sec` | Secondary structure encoding string (H, E, C, P, B, etc.). |
| `graft_offset` | (Optional) Index offset for re-stitching truncated domains. |

**Sources:** `mk_ldr_template.py` [108-125](), `mk_flex_template.py` [220-230]()

---

## Auxiliary PDB Slicing

The `slice_pdb.py` script is a standalone utility used to extract specific residue ranges from a PDB file. It is frequently used to isolate a single domain from a multi-domain AF2 structure before passing it to the template preparation scripts.

*   **Logic**: It parses `ATOM` and `HETATM` records and filters by the `resSeq` field [slice_pdb.py:15-20]().
*   **Renumbering**: It can optionally renumber the output PDB starting from 1 using the `--renumber` flag [slice_pdb.py:21-22]().

**Sources:** `slice_pdb.py` [7-40]()

---
# Step 4: Kinematic Stitching & Ensemble Assembly

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [AlphaFlex/Step_4_ldr_stitch.py](AlphaFlex/Step_4_ldr_stitch.py)
- [AlphaFlex/utils/graft_back.py](AlphaFlex/utils/graft_back.py)
- [AlphaFlex/utils/stitch.py](AlphaFlex/utils/stitch.py)
- [AlphaFlex/utils/truncate.py](AlphaFlex/utils/truncate.py)

</details>



Step 4 is the final stage of the AlphaFlex pipeline, responsible for integrating independently sampled Intrinsically Disordered Regions (IDRs) into their respective folded domain contexts to produce full-length protein ensembles. This process utilizes kinematic chain assembly to ensure geometric continuity, followed by physics-based relaxation and a multi-tiered structural validation gate.

### Overview of Assembly Workflow

The stitching process in `Step_4_ldr_stitch.py` handles three primary structural categories:
1.  **Category-0 (Full IDPs):** Fast-pass assembly that bypasses kinematic stitching, focusing on ensemble PDB finalization.
2.  **Category-1/2 (Tails and Linkers):** Kinematic assembly of one or more IDRs onto a static folded template.
3.  **Category-3 (Loops):** Complex stitching requiring two-point kinematic anchors to maintain loop closure.

For proteins that were truncated during Step 2 to fit within memory limits, this step performs a "Graft-Back" operation to restore the full-length sequence before relaxation.

Sources: [AlphaFlex/Step_4_ldr_stitch.py:1-12]() | [AlphaFlex/utils/stitch.py:112-128]()

---

### Implementation Architecture

The following diagram maps the high-level logic of `Step_4_ldr_stitch.py` to the underlying code entities and utility functions.

**Kinematic Assembly Logic Flow**
```mermaid
graph TD
    Start["process_protein()"] --> CatCheck{"Check Category"}
    
    CatCheck -- "Category 0" --> FastPass["Fast-pass Ensemble Finalization"]
    CatCheck -- "Category 1-3" --> LocateEnsembles["find_ensemble_dirs()"]
    
    LocateEnsembles --> TruncCheck{"Is Truncated?"}
    TruncCheck -- "Yes" --> Renum["_maybe_renumber_truncated_pool()"]
    Renum --> Graft["graft_conformer_multi()"]
    
    Graft --> StitchLoop["Monte Carlo Stitching Loop"]
    TruncCheck -- "No" --> StitchLoop
    
    subgraph "Stitching Core"
        StitchLoop --> Kinematic["assemble_kinematic_chain()"]
        Kinematic --> RelaxRepair["_relax_repair_validate()"]
    end
    
    RelaxRepair --> Finalize["finalize_model_pdb()"]
    Finalize --> End["mkensemble()"]
```
Sources: [AlphaFlex/Step_4_ldr_stitch.py:23-50](), [AlphaFlex/Step_4_ldr_stitch.py:191-220](), [AlphaFlex/utils/stitch.py:132-165]()

---

### Kinematic Chain Assembly

The core assembly logic resides in `assemble_kinematic_chain` [AlphaFlex/utils/stitch.py:254-325](). It treats the protein as a series of segments (folded domains and disordered regions) and uses BioPython's `Superimposer` to align them at junctions.

#### Segment Mapping and Alignment
The system builds a segment map using `build_segment_map` [AlphaFlex/utils/stitch.py:207-251](), which identifies the sequence ranges for all structural components. For each disordered segment, the assembler:
1.  Identifies the **Alignment Stub**: A small set of residues (defined by `ALIGNMENT_STUB_HALF_SIZE` [AlphaFlex/utils/stitch.py:20]()) at the junction between the folded domain and the IDR.
2.  Performs **Rigid Transformation**: Uses the Kabsch algorithm (via `Superimposer`) to transform the IDR conformer coordinates into the coordinate frame of the folded template.
3.  **Kinematic Propagation**: As each segment is added, subsequent segments are transformed relative to the new cumulative structure, maintaining the kinematic chain.

Sources: [AlphaFlex/utils/stitch.py:254-325](), [AlphaFlex/utils/graft_back.py:47-56]()

---

### Graft-Back for Truncated Pools

For large proteins where IDRs were sampled using truncated templates (Step 2), the `graft_back.py` utility restores the missing folded domains.

| Function | Purpose |
| :--- | :--- |
| `graft_conformer` | Restores distal folded domains for a single-IDR fragment using a 3-step Kabsch-based graft [AlphaFlex/utils/graft_back.py:93-160](). |
| `graft_conformer_multi` | Orchestrates independent grafting of multiple domains for complex linkers [AlphaFlex/utils/graft_back.py:162-180](). |
| `renumber_pool` | Maps truncated residue indices back to native UniProt numbering [AlphaFlex/utils/stitch.py:114](). |

The graft quality is monitored via the `junction_gap` and `seam_gap` metrics, which measure the Euclidean distance between the $C\alpha$ atoms at the splice site [AlphaFlex/utils/graft_back.py:145-159]().

Sources: [AlphaFlex/utils/graft_back.py:93-180](), [AlphaFlex/utils/truncate.py:5-63]()

---

### The _relax_repair_validate Quality Gate

Every stitched conformer must pass a multi-stage quality gate before being accepted into the final ensemble. This is encapsulated in `_relax_repair_validate` [AlphaFlex/Step_4_ldr_stitch.py:256-320]().

**Validation and Repair Pipeline**
```mermaid
graph LR
    subgraph "Relaxation & Repair"
        A["relax_with_established_method()"] --> B["fix_histidine_naming()"]
        B --> C["repair_chirality()"]
    end
    
    subgraph "Validation Gates"
        C --> D{"check_bond_integrity()"}
        D -- "Pass" --> E{"validate_structure_post_relax()"}
        E -- "Pass" --> F{"_check_fold_curvature()"}
    end
    
    F -- "Success" --> Accept["Accept Conformer"]
    D -- "Fail" --> Reject["Reject & Retry"]
```

#### AMBER Relaxation with Restraints
The function `relax_with_established_method` [AlphaFlex/Step_4_ldr_stitch.py:128-189]() performs energy minimization using the AMBER99SB forcefield via OpenMM. Crucially, it applies harmonic restraints to the folded domains while allowing the IDRs to relax freely. This is controlled by the `idr_indices` passed to the `viol_mask` in `relax_protein` [AlphaFlex/Step_4_ldr_stitch.py:163-171]().

#### Structural Integrity Checks
- **Chirality Repair:** `repair_chirality` [AlphaFlex/Step_4_ldr_stitch.py:102]() corrects D-amino acids and improper planarity resulting from the stitching or relaxation process.
- **Bond Integrity:** `check_bond_integrity` [AlphaFlex/Step_4_ldr_stitch.py:103]() ensures no unphysical bond lengths exist at the kinematic junctions.
- **Fold Curvature:** `_check_fold_curvature` [AlphaFlex/Step_4_ldr_stitch.py:91]() validates that the static folded domains have not been distorted during relaxation.

Sources: [AlphaFlex/Step_4_ldr_stitch.py:128-189](), [AlphaFlex/Step_4_ldr_stitch.py:256-320](), [idpforge/utils/relax.py:89]()

---

### Ensemble Finalization

Once the target number of conformers (defined in `config.py`) is reached, the single-model PDBs are aggregated.

1.  **PDB Finalization:** `finalize_model_pdb` [AlphaFlex/Step_4_ldr_stitch.py:90]() adds necessary metadata, TER records, and standardizes atom naming.
2.  **Ensemble Concatenation:** `mkensemble` [AlphaFlex/Step_4_ldr_stitch.py:38-48]() reads the validated PDB files and concatenates them into a single multi-MODEL PDB file, removing redundant headers like `CRYST1` or `MASTER` records.
3.  **Stale Cleanup:** `_clear_stale_ensembles` [AlphaFlex/Step_4_ldr_stitch.py:62-69]() ensures that previous failed or incomplete runs do not contaminate the final output directory.

Sources: [AlphaFlex/Step_4_ldr_stitch.py:38-69](), [idpforge/misc.py:90]()

---
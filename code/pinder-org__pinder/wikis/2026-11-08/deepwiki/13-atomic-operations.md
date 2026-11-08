# Atomic Operations

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/pinder-core/pinder/core/loader/structure.py](src/pinder-core/pinder/core/loader/structure.py)
- [src/pinder-core/pinder/core/structure/atoms.py](src/pinder-core/pinder/core/structure/atoms.py)

</details>



## Purpose and Scope

This page documents the low-level atomic manipulation functions that form the foundation of PINDER's structural operations. These functions operate directly on biotite's `AtomArray` and `AtomArrayStack` objects, providing capabilities for file I/O, filtering, alignment, and property calculations.

For higher-level structure manipulation through the `Structure` class wrapper, see [Structure Class](#3.2). For system-level operations that combine multiple structures, see [PinderSystem](#3.1).

## Overview of Atomic Operations

The atomic operations module provides a comprehensive toolkit for manipulating protein structures at the atomic level. These operations are organized into several functional categories:

```mermaid
graph TB
    subgraph "File I/O"
        READ["atom_array_from_pdb_file()<br/>Read PDB/CIF files"]
        WRITE["write_pdb()<br/>Write PDB files"]
        CONVERT["cif_to_pdb()<br/>Format conversion"]
    end
    
    subgraph "Filtering & Masking"
        FILTER["filter_atoms()<br/>Filter by atom type"]
        BBMASK["backbone_mask()<br/>Backbone atom mask"]
        APPLY["apply_mask()<br/>Apply boolean mask"]
        RESLIST["mask_from_res_list()<br/>Residue-based mask"]
    end
    
    subgraph "Sequence Operations"
        SEQALN["get_seq_alignments()<br/>Pairwise alignment"]
        STRUCSALN["get_seq_aligned_structures()<br/>Align structures by sequence"]
        CHAINMAP["get_per_chain_seq_alignments()<br/>Per-chain mappings"]
        SEQID["get_seq_identity()<br/>Calculate identity"]
    end
    
    subgraph "Structure Processing"
        STACK["stack_filter_intersect()<br/>Create AtomArrayStack"]
        RENUMBER["renumber_res_ids()<br/>Renumber residues"]
        STANDARD["standardize_atom_array()<br/>Standardize order"]
        ORIENT["normalize_orientation()<br/>Align principal axes"]
    end
    
    subgraph "Property Calculations"
        SASA["get_buried_sasa()<br/>Calculate dSASA"]
        VDW["atom_vdw_radius()<br/>Van der Waals radii"]
        RESOLVED["get_resolved_resi()<br/>Resolved residues"]
    end
    
    READ --> FILTER
    FILTER --> APPLY
    BBMASK --> APPLY
    SEQALN --> STRUCSALN
    STRUCSALN --> CHAINMAP
    FILTER --> STACK
    STACK --> SASA
    WRITE -.-> CONVERT
```

**Diagram: Atomic Operations Function Categories**

Sources: [src/pinder-core/pinder/core/structure/atoms.py:1-1010]()

## File I/O Operations

### Reading Structure Files

The primary entry point for reading protein structures is `atom_array_from_pdb_file()`, which supports multiple parsing backends:

| Function | Purpose | Backend Options |
|----------|---------|----------------|
| `atom_array_from_pdb_file()` | Read PDB files into AtomArray | `fastpdb` (default), `biotite` |
| `biotite_pdbfile()` | Get biotite PDB parser | biotite native |
| `rust_pdbfile()` | Get fastpdb parser (Rust-based) | fastpdb with fallback |
| `biotite_pdbxfile()` | Get CIF/mmCIF parser | biotite native |

**Example Usage Pattern:**
```python
# Fast parsing with fastpdb (default)
structure = atom_array_from_pdb_file(pdb_path, backend="fastpdb")

# Parse with additional fields
structure = atom_array_from_pdb_file(
    pdb_path, 
    backend="fastpdb",
    extra_fields=["b_factor"]
)
```

The function automatically handles path resolution and can accept either `Path` objects or `AtomArray` objects (pass-through). When using the `fastpdb` backend, parsing is significantly faster due to Rust implementation, but falls back to biotite if fastpdb is not installed.

Sources: [src/pinder-core/pinder/core/structure/atoms.py:53-71](), [src/pinder-core/pinder/core/structure/atoms.py:27-51]()

### Writing Structure Files

Writing structures back to disk is handled by `write_pdb()` and related functions:

```mermaid
graph LR
    ATOMARRAY["AtomArray<br/>In-memory structure"]
    WRITEPDB["write_pdb()"]
    PDBFILE["PDB File<br/>On disk"]
    CIFTOPDB["cif_to_pdb()"]
    CIFFILE["CIF File"]
    RENAME["rename_chains()"]
    
    ATOMARRAY --> WRITEPDB
    WRITEPDB --> PDBFILE
    CIFFILE --> CIFTOPDB
    CIFTOPDB --> PDBFILE
    PDBFILE --> RENAME
    RENAME --> PDBFILE
```

**Diagram: File Writing Operations Flow**

The `write_pdb()` function automatically creates parent directories if they don't exist and uses biotite's `strucio.save_structure()` internally. For format conversion, `cif_to_pdb()` provides CIF to PDB conversion with optional chain filtering and renaming.

Sources: [src/pinder-core/pinder/core/structure/atoms.py:890-909](), [src/pinder-core/pinder/core/structure/atoms.py:946-982](), [src/pinder-core/pinder/core/structure/atoms.py:911-944]()

## Atom Filtering and Masking

### Creating Boolean Masks

Masks are boolean arrays that specify which atoms to keep or filter. PINDER provides several mask creation functions:

| Function | Description | Returns |
|----------|-------------|---------|
| `backbone_mask()` | Creates mask for backbone atoms | `NDArray[np.bool_]` |
| `mask_from_res_list()` | Creates mask for specific residues | `NDArray[np.bool_]` |
| CA atom selection | Via `atom_array.atom_name == "CA"` | `NDArray[np.bool_]` |
| Heavy atom selection | Via `atom_array.element != "H"` | `NDArray[np.bool_]` |

The `backbone_mask()` function supports two backbone definitions:

```python
# DockQ backbone (C, CA, N, O)
mask = backbone_mask(atoms, backbone_definition="dockq")

# Standard peptide backbone (via biotite)
mask = backbone_mask(atoms, backbone_definition="standard")
```

Sources: [src/pinder-core/pinder/core/structure/atoms.py:87-94](), [src/pinder-core/pinder/core/structure/atoms.py:97-101](), [src/pinder-core/pinder/core/structure/atoms.py:21-21]()

### Applying Masks

The `apply_mask()` function handles both `AtomArray` and `AtomArrayStack` objects:

```mermaid
graph TB
    INPUT["Input: atoms + mask"]
    CHECK{"isinstance<br/>AtomArray?"}
    ARRAY["Return atoms[mask]"]
    STACK["Return atoms[..., mask]"]
    OUTPUT["Filtered structure"]
    
    INPUT --> CHECK
    CHECK -->|Yes| ARRAY
    CHECK -->|No, AtomArrayStack| STACK
    ARRAY --> OUTPUT
    STACK --> OUTPUT
```

**Diagram: Mask Application Logic**

The function correctly handles the different indexing semantics for single structures versus stacks of structures (e.g., multiple decoy models).

Sources: [src/pinder-core/pinder/core/structure/atoms.py:104-124]()

### High-Level Filtering

The `filter_atoms()` function provides a convenient interface for common filtering operations:

```python
# Filter to only C-alpha atoms
ca_only = filter_atoms(atoms, calpha_only=True)

# Filter to backbone atoms (no side chains)
backbone = filter_atoms(
    atoms, 
    backbone_only=True,
    backbone_definition="dockq"
)

# Filter to heavy atoms (no hydrogen)
heavy = filter_atoms(atoms, heavy_only=True)

# Combine multiple filters
ca_heavy = filter_atoms(
    atoms,
    calpha_only=True,
    heavy_only=True
)
```

Sources: [src/pinder-core/pinder/core/structure/atoms.py:127-141]()

## Sequence Alignment Operations

### Pairwise Sequence Alignment

Sequence alignment is performed using biotite's BLOSUM62 matrix with configurable gap penalties:

```mermaid
graph TB
    REF["Reference Sequence<br/>ref_seq: str"]
    SUBJ["Subject Sequence<br/>subject_seq: str"]
    
    ALIGN["get_seq_alignments()<br/>BLOSUM62 matrix<br/>gap_penalty=(-10, -1)"]
    
    OUTPUT["List[Alignment]<br/>Optimal global alignments"]
    
    SEQID["get_seq_identity()<br/>Calculate % identity"]
    IDENTITY["float: identity score"]
    
    ALIGN_SEQ["align_sequences()<br/>With numbering"]
    MAPPED["Aligned sequences<br/>+ numbering"]
    
    REF --> ALIGN
    SUBJ --> ALIGN
    ALIGN --> OUTPUT
    OUTPUT --> SEQID
    SEQID --> IDENTITY
    
    REF --> ALIGN_SEQ
    SUBJ --> ALIGN_SEQ
    ALIGN_SEQ --> MAPPED
```

**Diagram: Sequence Alignment Functions**

The `get_seq_alignments()` function returns optimal global alignments, while `get_seq_identity()` calculates the percentage identity based on the shortest sequence. For aligning with residue numbering, `align_sequences()` preserves both sequence and numbering information.

Sources: [src/pinder-core/pinder/core/structure/atoms.py:301-341](), [src/pinder-core/pinder/core/structure/atoms.py:344-368](), [src/pinder-core/pinder/core/structure/atoms.py:402-469]()

### Structure-Level Sequence Alignment

The `get_seq_aligned_structures()` function aligns two structures based on their sequences and renumbers residues to match:

```python
# Align two structures
ref_aligned, subject_aligned = get_seq_aligned_structures(
    ref=ref_structure,
    subject=subject_structure,
    pdb_engine="fastpdb"
)

# After alignment, residue numbers match for common sequence
assert ref_aligned.res_id[0] == subject_aligned.res_id[0]
```

**Processing Flow:**

```mermaid
graph TB
    REF["Reference Structure"]
    SUBJ["Subject Structure"]
    
    GETCHAINS["_get_paired_structures_and_chains()<br/>Extract chains"]
    
    PERCHAIN["For each chain pair"]
    
    GETSEQ["_get_structure_and_res_info()<br/>Extract sequences + numbering"]
    
    ALIGNMAP["_align_and_map_sequences()<br/>Create residue mapping"]
    
    GETALN["_get_seq_aligned_structures()<br/>Apply mapping to structures"]
    
    COMBINE["Combine aligned chains"]
    
    OUTPUT["ref_aligned, subject_aligned<br/>Matching residue numbering"]
    
    REF --> GETCHAINS
    SUBJ --> GETCHAINS
    GETCHAINS --> PERCHAIN
    PERCHAIN --> GETSEQ
    GETSEQ --> ALIGNMAP
    ALIGNMAP --> GETALN
    GETALN --> COMBINE
    COMBINE --> OUTPUT
```

**Diagram: Structure Alignment Processing Pipeline**

The function processes each chain independently, aligns their sequences, creates residue mappings, and then reconstructs the multi-chain structure with consistent numbering.

Sources: [src/pinder-core/pinder/core/structure/atoms.py:598-668](), [src/pinder-core/pinder/core/structure/atoms.py:538-567](), [src/pinder-core/pinder/core/structure/atoms.py:570-595]()

### Per-Chain Residue Mappings

For cases where you need residue number mappings without modifying structures, use `get_per_chain_seq_alignments()`:

```python
# Get residue mappings for each chain
mappings = get_per_chain_seq_alignments(
    ref=ref_structure,
    subject=subject_structure
)

# mappings = {
#     'A': {7: 8, 8: 9, 9: 10, ...},  # subject res -> ref res
#     'B': {15: 20, 16: 21, ...}
# }
```

The returned dictionary maps chain IDs to residue number mappings, where keys are subject residue IDs and values are reference residue IDs. This is particularly useful for transferring annotations or identifying corresponding residues between structures.

To invert these mappings (ref -> subject instead of subject -> ref), use `invert_chain_seq_map()`:

```python
# Invert mappings
ref_to_subject = invert_chain_seq_map(mappings)
```

Sources: [src/pinder-core/pinder/core/structure/atoms.py:671-731](), [src/pinder-core/pinder/core/structure/atoms.py:734-750]()

## Structure Processing and Standardization

### Creating Atom Array Stacks

When working with multiple models (e.g., docking decoys), `stack_filter_intersect()` creates an `AtomArrayStack` by ensuring all structures have matching atoms:

```mermaid
graph TB
    INPUT["List[AtomArray]<br/>Multiple decoy structures"]
    
    ITER["Iterate through decoys"]
    
    INTERSECT["struc.filter_intersection()<br/>Find common atoms"]
    
    CHECK{"All atoms<br/>matched?"}
    
    RMELEMENTS["Remove element<br/>annotations"]
    
    STANDARDIZE["standardize_atom_array()<br/>Renumber and order"]
    
    STACK["stack()<br/>Create AtomArrayStack"]
    
    OUTPUT["AtomArrayStack<br/>Uniform atom structure"]
    
    INPUT --> ITER
    ITER --> INTERSECT
    INTERSECT --> CHECK
    CHECK -->|No, mismatch| RMELEMENTS
    RMELEMENTS --> INTERSECT
    CHECK -->|Yes| STANDARDIZE
    STANDARDIZE --> STACK
    STACK --> OUTPUT
```

**Diagram: Stack Creation with Intersection Filtering**

The function iteratively applies `filter_intersection` to ensure all structures have the same atoms in the same order. If mismatches occur, it progressively removes problematic annotations (like element types) and can renumber/reorder atoms to achieve consistency.

Sources: [src/pinder-core/pinder/core/structure/atoms.py:201-298]()

### Standardization Operations

| Function | Purpose | Key Operation |
|----------|---------|---------------|
| `standardize_atom_array()` | Standardize atom order and numbering | `struc.info.standardize_order()` + `renumber_res_ids()` |
| `renumber_res_ids()` | Create continuous residue numbering | `struc.create_continuous_res_ids()` |
| `normalize_orientation()` | Align structure to principal axes | `struc.orient_principal_components()` |

These functions ensure structures are in a consistent state for downstream operations:

```python
# Standardize a structure
standardized = standardize_atom_array(arr, standardize=True)

# Renumber residues starting from 1
renumbered = renumber_res_ids(arr)

# Orient structure along principal components
normalized = normalize_orientation(arr)
```

Sources: [src/pinder-core/pinder/core/structure/atoms.py:190-198](), [src/pinder-core/pinder/core/structure/atoms.py:180-187](), [src/pinder-core/pinder/core/structure/atoms.py:984-1010]()

## SASA and Property Calculations

### Buried Surface Area

The `get_buried_sasa()` function calculates the change in solvent-accessible surface area (dSASA) upon complex formation:

**Calculation Formula:**
```
dSASA = SASA(A) + SASA(B) - SASA(A+B)
```

Where:
- `SASA(A)` = surface area of chain A in isolation
- `SASA(B)` = surface area of chain B in isolation  
- `SASA(A+B)` = surface area of the complex

**Implementation Details:**

```mermaid
graph TB
    CHAIN_A["Chain A AtomArray"]
    CHAIN_B["Chain B AtomArray"]
    
    SASA_A["struc.sasa(A)<br/>point_number=200"]
    SASA_B["struc.sasa(B)<br/>point_number=200"]
    SASA_AB["struc.sasa(A+B)<br/>point_number=200"]
    
    SUM_A["Sum SASA A"]
    SUM_B["Sum SASA B"]
    SUM_AB["Sum SASA AB"]
    
    CALC["dSASA = A + B - AB"]
    
    ERROR{"Error with<br/>ProtOr radii?"}
    
    FALLBACK["Retry with<br/>vdw_radii='Single'"]
    
    ROUND["Round to int"]
    OUTPUT["dSASA (int)"]
    
    CHAIN_A --> SASA_A
    CHAIN_B --> SASA_B
    CHAIN_A --> SASA_AB
    CHAIN_B --> SASA_AB
    
    SASA_A --> SUM_A
    SASA_B --> SUM_B
    SASA_AB --> SUM_AB
    
    SUM_A --> CALC
    SUM_B --> CALC
    SUM_AB --> CALC
    
    CALC --> ERROR
    ERROR -->|Yes| FALLBACK
    FALLBACK --> CALC
    ERROR -->|No| ROUND
    ROUND --> OUTPUT
```

**Diagram: Buried SASA Calculation with Error Handling**

The function has built-in error handling for cases where ProtOr van der Waals radii fail, automatically falling back to single-element radii. It also handles NaN values and returns the rounded integer value.

Sources: [src/pinder-core/pinder/core/structure/atoms.py:753-844]()

### Van der Waals Radii

The `atom_vdw_radius()` function retrieves appropriate van der Waals radii for atoms:

```python
# Get VDW radius for an atom
radius = atom_vdw_radius(atom)

# Logic:
# 1. If hetero atom -> use element-based radius
# 2. If standard residue -> use ProtOr residue/atom-specific radius
# 3. If ProtOr fails -> fall back to element-based radius
# 4. If all fails -> use nitrogen radius as default
```

Sources: [src/pinder-core/pinder/core/structure/atoms.py:73-84]()

### Resolved Residues

Two functions provide information about which residues are resolved in a structure:

```python
# From file path
resolved = get_resolved_resi(pdb_path)
# Returns: {'A': [1, 2, 3, ...], 'B': [1, 2, 3, ...]}

# From AtomArray
resolved = get_resolved_resi_from_atom_array(atom_array)
# Returns: same format, dict of chain -> list of residue IDs
```

Both functions check for the presence of C-alpha atoms to determine resolution, returning a per-chain dictionary of resolved residue IDs.

Sources: [src/pinder-core/pinder/core/structure/atoms.py:847-887]()

## Utility Functions

### Sequence Conversion

The `resn2seq()` function converts three-letter residue codes to one-letter amino acid codes:

```python
# Convert residue names to sequence
residue_names = ["ALA", "CYS", "ASP", "GLU"]
sequence = resn2seq(residue_names)
# Returns: "ACDE"
```

The function handles non-canonical amino acids by mapping them to 'X' and includes special handling for selenocysteine (U->C) and pyrrolysine (O->K).

Sources: [src/pinder-core/pinder/core/structure/atoms.py:521-535]()

### Receptor/Ligand Assignment

The `assign_receptor_ligand()` function automatically assigns chains based on size:

```python
chains = ["A", "B", "C"]
receptor_chains, ligand_chains = assign_receptor_ligand(arr, chains)

# Logic: 
# - Receptor = all but smallest chain
# - Ligand = smallest chain
```

This is particularly useful for asymmetric complexes where the binding partner needs to be identified based on size.

Sources: [src/pinder-core/pinder/core/structure/atoms.py:163-177]()

### Chain Manipulation

Two functions handle chain naming operations:

| Function | Purpose | Usage |
|----------|---------|-------|
| `rename_chains()` | Rename chains in PDB file | `rename_chains(pdb_path, {"A": "R", "B": "L"})` |
| `cif_to_pdb()` | Convert CIF with chain filtering | `cif_to_pdb(cif_path, pdb_path, chains={"A": "R"})` |

Both functions can accept either a string (to rename all chains) or a dictionary mapping for selective renaming.

Sources: [src/pinder-core/pinder/core/structure/atoms.py:911-944](), [src/pinder-core/pinder/core/structure/atoms.py:946-982]()

## Integration with Structure Class

The atomic operations described here are heavily used by the `Structure` class ([Structure Class](#3.2)). Key integration points include:

```mermaid
graph LR
    STRUCTURE["Structure Class<br/>High-level wrapper"]
    
    READ["atom_array_from_pdb_file()"]
    FILTER["filter_atoms()"]
    ALIGN["get_seq_aligned_structures()"]
    WRITE["write_pdb()"]
    CHAINMAP["get_per_chain_seq_alignments()"]
    
    STRUCTURE -->|"read_pdb()"| READ
    STRUCTURE -->|"filter()"| FILTER
    STRUCTURE -->|"align_common_sequence()"| ALIGN
    STRUCTURE -->|"to_pdb()"| WRITE
    STRUCTURE -->|"get_per_chain_seq_alignments()"| CHAINMAP
```

**Diagram: Structure Class Using Atomic Operations**

The `Structure` class (defined in [src/pinder-core/pinder/core/loader/structure.py]()) wraps these atomic operations with additional functionality like UniProt mapping, contact calculation, and property access methods.

Sources: [src/pinder-core/pinder/core/loader/structure.py:50-64](), [src/pinder-core/pinder/core/loader/structure.py:122-142](), [src/pinder-core/pinder/core/loader/structure.py:156-209](), [src/pinder-core/pinder/core/loader/structure.py:66-83](), [src/pinder-core/pinder/core/loader/structure.py:147-154]()

## Backend Performance Considerations

### FastPDB vs Biotite

PINDER supports two parsing backends with different performance characteristics:

| Backend | Implementation | Speed | Fallback Behavior |
|---------|---------------|-------|-------------------|
| `fastpdb` | Rust-based (external dependency) | Fast (10-100x faster) | Falls back to biotite if not installed |
| `biotite` | Pure Python | Slower | Always available |

The default is `fastpdb` for performance, but the system gracefully degrades to biotite when fastpdb is unavailable. This is particularly important for large-scale dataset loading operations.

Sources: [src/pinder-core/pinder/core/structure/atoms.py:39-50](), [src/pinder-core/pinder/core/loader/structure.py:48-48]()

## Common Usage Patterns

### Reading and Filtering a Structure

```python
from pinder.core.structure.atoms import (
    atom_array_from_pdb_file,
    filter_atoms,
    write_pdb
)

# Read structure
atoms = atom_array_from_pdb_file(pdb_path, backend="fastpdb")

# Filter to backbone heavy atoms
backbone = filter_atoms(
    atoms,
    backbone_only=True,
    heavy_only=True,
    backbone_definition="dockq"
)

# Write filtered structure
write_pdb(backbone, output_path)
```

### Aligning Two Structures

```python
from pinder.core.structure.atoms import get_seq_aligned_structures

# Align structures based on sequence
ref_aligned, mobile_aligned = get_seq_aligned_structures(
    ref=reference_pdb,
    subject=mobile_pdb
)

# Structures now have matching residue numbering
# for aligned sequence segments
```

### Calculating Interface Properties

```python
from pinder.core.structure.atoms import (
    get_buried_sasa,
    filter_atoms
)

# Get chains
chain_a = atoms[atoms.chain_id == "A"]
chain_b = atoms[atoms.chain_id == "B"]

# Calculate buried surface area
dsasa = get_buried_sasa(chain_a, chain_b, point_number=200)
print(f"Buried SASA: {dsasa} Å²")
```

Sources: [src/pinder-core/pinder/core/structure/atoms.py]()

---
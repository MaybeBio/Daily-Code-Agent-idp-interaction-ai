# Core Data Structures

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/pinder-core/pinder/core/index/system.py](src/pinder-core/pinder/core/index/system.py)
- [src/pinder-core/pinder/core/loader/structure.py](src/pinder-core/pinder/core/loader/structure.py)

</details>



This document describes the fundamental data structures used throughout PINDER to represent protein systems, individual structures, and associated metadata. These classes form the foundation for all data loading, manipulation, and analysis operations in the codebase.

For information about how these structures are used in data loading pipelines, see [Data Loading System](#4). For details on atomic-level manipulation functions, see [Atomic Operations](#3.3).

## Overview

PINDER organizes protein interaction data using a hierarchical set of data structures:

- **`PinderSystem`**: High-level container managing all structural variants of a protein complex
- **`Structure`**: Wrapper around biotite's `AtomArray` with additional functionality
- **`IndexEntry`**: Metadata schema defining available files and properties for each system
- **`MetadataEntry`**: Extended annotation schema with 45+ biological and structural properties

```mermaid
graph TB
    IndexEntry["IndexEntry<br/>(index.parquet)"]
    MetadataEntry["MetadataEntry<br/>(metadata.parquet)"]
    PinderSystem["PinderSystem"]
    Native["Structure<br/>(native)"]
    HoloR["Structure<br/>(holo_receptor)"]
    HoloL["Structure<br/>(holo_ligand)"]
    ApoR["Structure<br/>(apo_receptor)"]
    ApoL["Structure<br/>(apo_ligand)"]
    PredR["Structure<br/>(pred_receptor)"]
    PredL["Structure<br/>(pred_ligand)"]
    AtomArray["biotite.AtomArray"]
    UniprotMap["pd.DataFrame<br/>(uniprot_map)"]
    
    IndexEntry -->|"__init__ param"| PinderSystem
    MetadataEntry -.->|"metadata property"| PinderSystem
    PinderSystem -->|"native"| Native
    PinderSystem -->|"@cached_property"| HoloR
    PinderSystem -->|"@cached_property"| HoloL
    PinderSystem -->|"@cached_property"| ApoR
    PinderSystem -->|"@cached_property"| ApoL
    PinderSystem -->|"@cached_property"| PredR
    PinderSystem -->|"@cached_property"| PredL
    
    Native --> AtomArray
    HoloR --> AtomArray
    HoloL --> AtomArray
    ApoR --> AtomArray
    ApoL --> AtomArray
    PredR --> AtomArray
    PredL --> AtomArray
    
    Native -.-> UniprotMap
    HoloR -.-> UniprotMap
    HoloL -.-> UniprotMap
```

**Sources:** [src/pinder-core/pinder/core/index/system.py:89-162](), [src/pinder-core/pinder/core/loader/structure.py:42-49]()

## Structure Class Hierarchy

The `Structure` class is a pydantic dataclass that wraps biotite's `AtomArray` to provide domain-specific functionality for protein structure manipulation.

### Core Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `filepath` | `Path` | Path to PDB file on disk |
| `atom_array` | `AtomArray` | Biotite structure containing atomic coordinates and annotations |
| `uniprot_map` | `Path \| pd.DataFrame \| None` | UniProt residue mapping for coordinate system conversion |
| `pinder_id` | `str \| None` | Unique identifier within PINDER dataset |
| `pdb_engine` | `str` | Parser engine ("fastpdb" or "biotite") |

**Sources:** [src/pinder-core/pinder/core/loader/structure.py:42-49]()

### Key Properties and Methods

```mermaid
graph LR
    Structure["Structure"]
    
    subgraph "Reading/Writing"
        ReadPDB["read_pdb()<br/>(static)"]
        ToPDB["to_pdb()"]
    end
    
    subgraph "Sequence Operations"
        Sequence["sequence property"]
        ChainSeq["chain_sequence property"]
        Fasta["fasta property"]
        TokenSeq["tokenized_sequence property"]
    end
    
    subgraph "Coordinate Operations"
        Coords["coords property"]
        Chains["chains property"]
        Residues["residues property"]
        Superimpose["superimpose()"]
        AlignSeq["align_common_sequence()"]
    end
    
    subgraph "Filtering"
        Filter["filter()"]
        BackboneMask["backbone_mask property"]
        CAlphaMask["calpha_mask property"]
    end
    
    subgraph "Interface Analysis"
        GetContacts["get_contacts()"]
        GetInterfaceRes["get_interface_residues()"]
        GetInterfaceMask["get_interface_mask()"]
    end
    
    subgraph "UniProt Mapping"
        UniprotMapping["uniprot_mapping property"]
        ResolvedMapping["resolved_mapping property"]
        PDB2Uniprot["resolved_pdb2uniprot property"]
        Uniprot2PDB["resolved_uniprot2pdb property"]
    end
    
    Structure --> ReadPDB
    Structure --> ToPDB
    Structure --> Sequence
    Structure --> ChainSeq
    Structure --> Coords
    Structure --> Chains
    Structure --> Superimpose
    Structure --> AlignSeq
    Structure --> Filter
    Structure --> GetContacts
    Structure --> UniprotMapping
```

**Sources:** [src/pinder-core/pinder/core/loader/structure.py:50-512]()

## PinderSystem Architecture

`PinderSystem` is the primary interface for working with protein complexes in PINDER. It manages all structural variants of a protein-protein interaction and provides methods for complex creation, alignment, and difficulty assessment.

### Initialization and Data Access

```mermaid
graph TB
    Init["PinderSystem.__init__()"]
    Entry["entry: IndexEntry"]
    DatasetPath["dataset_path: Path"]
    
    subgraph "Directory Structure"
        PDBsPath["pdbs_path"]
        TestPDBsPath["test_pdbs_path"]
        MappingsPath["mappings_path"]
    end
    
    subgraph "Eager Loading"
        Native["native: Structure"]
        NativeR["native_R: Structure"]
        NativeL["native_L: Structure"]
    end
    
    subgraph "Lazy Loading (@cached_property)"
        HoloR["holo_receptor: Structure"]
        HoloL["holo_ligand: Structure"]
        AlignedHoloR["aligned_holo_R: Structure"]
        AlignedHoloL["aligned_holo_L: Structure"]
        ApoR["apo_receptor: Structure"]
        ApoL["apo_ligand: Structure"]
        PredR["pred_receptor: Structure"]
        PredL["pred_ligand: Structure"]
    end
    
    Init --> Entry
    Init --> DatasetPath
    Init --> PDBsPath
    Init --> TestPDBsPath
    Init --> MappingsPath
    
    Init --> Native
    Native --> NativeR
    Native --> NativeL
    
    PinderSystem -->|"First access triggers load"| HoloR
    PinderSystem -->|"First access triggers load"| HoloL
    PinderSystem -->|"Computed from holo_receptor"| AlignedHoloR
    PinderSystem -->|"Computed from holo_ligand"| AlignedHoloL
    PinderSystem -->|"First access triggers load"| ApoR
    PinderSystem -->|"First access triggers load"| ApoL
    PinderSystem -->|"First access triggers load"| PredR
    PinderSystem -->|"First access triggers load"| PredL
```

**Sources:** [src/pinder-core/pinder/core/index/system.py:110-162](), [src/pinder-core/pinder/core/index/system.py:163-261]()

### Structure Loading Pattern

The `PinderSystem` uses Python's `@cached_property` decorator to implement lazy loading:

1. **Eager loading** (during `__init__`): Only the native complex is loaded immediately
2. **Lazy loading**: Individual monomers are loaded on first access and cached thereafter
3. **Alignment**: `aligned_holo_R` and `aligned_holo_L` are superimposed to native coordinates

This pattern minimizes memory usage and I/O operations when only a subset of structures is needed.

**Sources:** [src/pinder-core/pinder/core/index/system.py:163-261]()

## Complex Creation Methods

`PinderSystem` provides three methods for creating dimer complexes from monomers:

| Method | Purpose | Key Parameters |
|--------|---------|----------------|
| `create_complex()` | Generic complex assembly from any receptor/ligand pair | `receptor`, `ligand`, `remove_differing_atoms`, `renumber_residues` |
| `create_apo_complex()` | Assemble apo monomers (fallback to holo if unavailable) | `remove_differing_atoms`, `renumber_residues` |
| `create_pred_complex()` | Assemble predicted monomers (fallback to holo if unavailable) | `remove_differing_atoms`, `renumber_residues` |
| `create_masked_bound_unbound_complexes()` | Create equal-shaped holo/apo/pred complexes via pairwise masking | `monomer_types`, `remove_differing_atoms`, `renumber_residues` |

### Complex Creation Workflow

```mermaid
graph TB
    CreateComplex["create_complex()"]
    Receptor["receptor: Structure"]
    Ligand["ligand: Structure"]
    
    AlignSeqR["receptor.align_common_sequence()<br/>(to aligned_holo_R)"]
    AlignSeqL["ligand.align_common_sequence()<br/>(to aligned_holo_L)"]
    
    SuperimposeR["receptor.superimpose()<br/>(to holo_R)"]
    SuperimposeL["ligand.superimpose()<br/>(to holo_L)"]
    
    Combine["R_super + L_super<br/>(Structure.__add__)"]
    Complex["binary: Structure"]
    
    CreateComplex --> Receptor
    CreateComplex --> Ligand
    Receptor --> AlignSeqR
    Ligand --> AlignSeqL
    AlignSeqR --> SuperimposeR
    AlignSeqL --> SuperimposeL
    SuperimposeR --> Combine
    SuperimposeL --> Combine
    Combine --> Complex
```

**Sources:** [src/pinder-core/pinder/core/index/system.py:376-422](), [src/pinder-core/pinder/core/index/system.py:424-488]()

## IndexEntry and MetadataEntry

These dataclasses define the schema for PINDER's index and metadata files stored as Parquet tables.

### IndexEntry Schema

`IndexEntry` contains 58 fields organized into categories:

| Category | Example Fields | Description |
|----------|---------------|-------------|
| Identifiers | `id`, `pdb_id`, `entry_id` | Unique system identifiers |
| File paths | `pinder_pdb`, `holo_R_pdb`, `apo_L_pdb`, `predicted_R_pdb` | Relative paths to structure files |
| Splits | `split`, `cluster_id`, `cluster_id_R`, `cluster_id_L` | Dataset split assignments |
| Availability | `apo_R`, `apo_L`, `predicted_R`, `predicted_L` | Boolean flags for structure availability |
| Quality | `resolution`, `method`, `date` | Crystallographic metadata |

### Property Access Pattern

```mermaid
graph LR
    PS["PinderSystem"]
    Entry["entry: IndexEntry"]
    Metadata["metadata: MetadataEntry"]
    
    IndexPqt["index.parquet<br/>(2.4M+ rows)"]
    MetadataPqt["metadata.parquet<br/>(45+ annotations)"]
    
    GetIndex["get_index()"]
    GetMetadata["get_metadata()"]
    
    PS --> Entry
    PS -.->|"Lazy load via property"| Metadata
    
    GetIndex --> IndexPqt
    GetMetadata --> MetadataPqt
    
    IndexPqt --> Entry
    MetadataPqt --> Metadata
```

**Sources:** [src/pinder-core/pinder/core/index/system.py:12-18](), [src/pinder-core/pinder/core/index/system.py:675-693]()

## Structure Operations

The `Structure` class provides several categories of operations:

### Sequence Alignment and Superposition

```mermaid
graph TB
    AlignCommon["align_common_sequence()"]
    
    subgraph "Sequence Alignment"
        GetPerChain["get_per_chain_seq_alignments()"]
        InvertMap["invert_chain_seq_map()"]
        GetSeqAligned["get_seq_aligned_structures()"]
    end
    
    subgraph "Atom Intersection"
        FixAnnotation["surgery.fix_annotation_mismatch()"]
        FilterIntersection["struc.filter_intersection()"]
    end
    
    subgraph "Residue Renumbering"
        ApplyMapping["Apply res_id mapping"]
    end
    
    Superimpose["superimpose()"]
    SuperimposeCommon["_superimpose_common_atoms()"]
    CalcRMSD["struc.rmsd()"]
    
    AlignCommon --> GetPerChain
    GetPerChain --> InvertMap
    InvertMap --> GetSeqAligned
    GetSeqAligned --> FixAnnotation
    FixAnnotation --> FilterIntersection
    FilterIntersection --> ApplyMapping
    
    Superimpose --> SuperimposeCommon
    SuperimposeCommon --> CalcRMSD
```

**Sources:** [src/pinder-core/pinder/core/loader/structure.py:156-209](), [src/pinder-core/pinder/core/loader/structure.py:296-326](), [src/pinder-core/pinder/core/loader/structure.py:636-674]()

### Interface Analysis

The `Structure` class provides methods to identify and analyze protein-protein interfaces:

| Method | Return Type | Description |
|--------|-------------|-------------|
| `get_contacts()` | `set[tuple[str, str, int, int]]` | Pairs of contacting residues within radius |
| `get_interface_residues()` | `dict[str, list[int]]` | Residue IDs at interface per chain |
| `get_interface_mask()` | `NDArray[np.bool_]` | Boolean mask for interface atoms |

**Sources:** [src/pinder-core/pinder/core/loader/structure.py:211-294]()

## UniProt Coordinate Mapping

PINDER maintains mappings between PDB residue numbering and UniProt canonical numbering to enable cross-structure comparisons.

### Mapping Properties

```mermaid
graph LR
    Structure["Structure"]
    UniprotMap["uniprot_map<br/>(Path or DataFrame)"]
    
    subgraph "Derived Properties"
        UniprotMapping["uniprot_mapping<br/>(DataFrame)"]
        ResolvedMapping["resolved_mapping<br/>(filtered to resolved res)"]
        PDB2Uni["resolved_pdb2uniprot<br/>(dict[int, int])"]
        Uni2PDB["resolved_uniprot2pdb<br/>(dict[int, int])"]
    end
    
    Structure --> UniprotMap
    UniprotMap --> UniprotMapping
    UniprotMapping --> ResolvedMapping
    ResolvedMapping --> PDB2Uni
    PDB2Uni --> Uni2PDB
```

The `mask_common_uniprot()` function uses these mappings to crop structures to their common residue sets:

1. Load UniProt mappings for both structures
2. Find intersection of UniProt residue IDs
3. Map back to PDB residue IDs for each structure
4. Filter `AtomArray` to keep only common residues

**Sources:** [src/pinder-core/pinder/core/loader/structure.py:329-371](), [src/pinder-core/pinder/core/loader/structure.py:532-606]()

## Data Flow Summary

```mermaid
graph TB
    Disk["Disk Storage<br/>(pdbs/, mappings/)"]
    
    subgraph "Index Layer"
        IndexPqt["index.parquet"]
        MetadataPqt["metadata.parquet"]
        IndexEntry["IndexEntry"]
        MetadataEntry["MetadataEntry"]
    end
    
    subgraph "System Layer"
        PS["PinderSystem"]
        Download["download_entry()"]
    end
    
    subgraph "Structure Layer"
        LoadStructure["load_structure()"]
        Structure["Structure"]
        AtomArray["biotite.AtomArray"]
        ReadPDB["read_pdb()"]
    end
    
    subgraph "Operations"
        Filter["filter()"]
        Superimpose["superimpose()"]
        AlignSeq["align_common_sequence()"]
        CreateComplex["create_complex()"]
    end
    
    IndexPqt --> IndexEntry
    MetadataPqt --> MetadataEntry
    IndexEntry --> PS
    MetadataEntry -.-> PS
    
    PS --> Download
    Download --> Disk
    
    PS --> LoadStructure
    LoadStructure --> Disk
    LoadStructure --> Structure
    Structure --> ReadPDB
    ReadPDB --> AtomArray
    
    Structure --> Filter
    Structure --> Superimpose
    Structure --> AlignSeq
    PS --> CreateComplex
```

**Sources:** [src/pinder-core/pinder/core/index/system.py:1-806](), [src/pinder-core/pinder/core/loader/structure.py:1-675]()

---
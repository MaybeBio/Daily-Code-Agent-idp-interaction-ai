# Dataset Structure

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [README.md](README.md)

</details>



This page describes the physical organization of the PINDER dataset on disk, including directory layout, file naming conventions, and storage requirements. For information about how to load and access this data programmatically, see [Loading Your First System](#2.1). For details on the dataset splits and benchmarks, see [Dataset Splits and Benchmarks](#1.4).

## Dataset Location

The PINDER dataset is stored in a version-specific directory determined by the following hierarchy:

1. `PINDER_BASE_DIR` environment variable (if set)
2. `$XDG_DATA_HOME/pinder` (typically `~/.local/share/pinder` on Mac and Linux)
3. Default fallback location

The complete path follows the pattern: `<base_dir>/<release_version>/`, where `<release_version>` is in `YYYY-MM` format (e.g., `2024-02`).

```mermaid
graph TB
    ENV["PINDER_BASE_DIR env var"]
    XDG["$XDG_DATA_HOME/pinder"]
    DEFAULT["Default location"]
    
    CHECK1{"PINDER_BASE_DIR<br/>set?"}
    CHECK2{"XDG_DATA_HOME<br/>set?"}
    
    BASE["Base Directory"]
    VERSION["/<release_version>/"]
    FINAL["Final Dataset Path"]
    
    ENV --> CHECK1
    CHECK1 -->|Yes| BASE
    CHECK1 -->|No| CHECK2
    CHECK2 -->|Yes| XDG
    CHECK2 -->|No| DEFAULT
    XDG --> BASE
    DEFAULT --> BASE
    
    BASE --> VERSION
    VERSION --> FINAL
```

**Diagram: Dataset location resolution logic**

The current location can be queried programmatically using `get_pinder_location()` from [src/pinder-core/pinder/core/utils/paths.py]().

Sources: [README.md:122-145]()

## Root Directory Structure

The dataset root contains five main components organized as follows:

```mermaid
graph TB
    ROOT["Dataset Root<br/>~/.local/share/pinder/2024-02/"]
    
    PDBS["pdbs/<br/>~672GB"]
    TEST_PDBS["test_set_pdbs/<br/>~705MB"]
    MAPPINGS["mappings/<br/>~25GB"]
    INDEX["index.parquet<br/>~50MB"]
    METADATA["metadata.parquet<br/>~100MB"]
    
    ROOT --> PDBS
    ROOT --> TEST_PDBS
    ROOT --> MAPPINGS
    ROOT --> INDEX
    ROOT --> METADATA
    
    PDBS_DESC["Individual monomer and<br/>ground-truth dimer PDBs"]
    TEST_DESC["Pre-transformed test<br/>set structures"]
    MAP_DESC["PDB↔UniProt mappings<br/>Assembly information"]
    IDX_DESC["Master index of<br/>2.4M+ entries"]
    META_DESC["Detailed metadata<br/>45+ annotation fields"]
    
    PDBS -.-> PDBS_DESC
    TEST_PDBS -.-> TEST_DESC
    MAPPINGS -.-> MAP_DESC
    INDEX -.-> IDX_DESC
    METADATA -.-> META_DESC
```

**Diagram: Root directory structure and component sizes**

Sources: [README.md:162-177](), [README.md:198-208]()

## Storage Requirements

| Component | Compressed | Unpacked | Description |
|-----------|-----------|----------|-------------|
| `pdbs.zip` | 144 GB | 672 GB | Structural data files |
| `test_set_pdbs.zip` | 149 MB | 705 MB | Test set structures |
| `mappings.zip` | 6.8 GB | 25 GB | Mapping files |
| `index.parquet` | - | ~50 MB | Master index |
| `metadata.parquet` | - | ~100 MB | Detailed annotations |
| **Total** | **~151 GB** | **~698 GB** | Complete dataset |

Sources: [README.md:198-208]()

## pdbs/ Directory

The `pdbs/` directory contains individual PDB files for monomers and native complexes. Files are organized by two-character PDB ID prefixes:

```mermaid
graph TB
    PDBS["pdbs/"]
    
    subgraph "Two-character prefix directories"
        DIR_1A["1a/"]
        DIR_1B["1b/"]
        DIR_XX["xx/"]
    end
    
    subgraph "Structure types"
        NATIVE["Native complexes<br/>&lt;id&gt;.pdb"]
        HOLO_R["Holo receptors<br/>&lt;id&gt;-R.pdb"]
        HOLO_L["Holo ligands<br/>&lt;id&gt;-L.pdb"]
        APO["Apo monomers<br/>&lt;uniprot&gt;--&lt;pdb&gt;__&lt;chain&gt;.pdb"]
        PRED["Predicted monomers<br/>&lt;uniprot&gt;.pdb"]
    end
    
    PDBS --> DIR_1A
    PDBS --> DIR_1B
    PDBS --> DIR_XX
    
    DIR_1A --> NATIVE
    DIR_1A --> HOLO_R
    DIR_1A --> HOLO_L
    DIR_1A --> APO
    DIR_1A --> PRED
```

**Diagram: pdbs/ directory organization**

### File Naming Conventions

The `pdbs/` directory uses specific naming patterns for different structure types:

| Structure Type | Naming Pattern | Example | Chain ID |
|---------------|----------------|---------|----------|
| Native complex | `<pdb_id>.pdb` | `1a2k.pdb` | Multiple chains |
| Holo receptor | `<pdb_id>-R.pdb` | `1a2k-R.pdb` | R |
| Holo ligand | `<pdb_id>-L.pdb` | `1a2k-L.pdb` | L |
| Apo monomer | `<uniprot>--<pdb_id>__<chain>_<uniprot>.pdb` | `P12345--2xyz__A1_P12345.pdb` | A |
| Predicted monomer | `<uniprot>.pdb` | `P12345.pdb` | A |

**Note:** Apo and predicted monomers are stored with chain ID `A` to reduce duplication when the same structure maps to both receptor and ligand roles across multiple systems. The `PinderSystem` class automatically re-chains these to `R` or `L` as needed.

Sources: [README.md:173-176](), [README.md:404-405]()

## test_set_pdbs/ Directory

The `test_set_pdbs/` directory contains pre-rotated and pre-translated monomer structures for test sets, ready for immediate use in evaluation:

```mermaid
graph TB
    TEST_ROOT["test_set_pdbs/"]
    
    SPLIT_S["pinder_s/"]
    SPLIT_XL["pinder_xl/"]
    SPLIT_AF2["pinder_af2/"]
    
    TEST_ROOT --> SPLIT_S
    TEST_ROOT --> SPLIT_XL
    TEST_ROOT --> SPLIT_AF2
    
    subgraph "Per-system directories"
        SYS_DIR["&lt;system_id&gt;/"]
        HOLO_R["holo_receptor.pdb"]
        HOLO_L["holo_ligand.pdb"]
        APO_R["apo_receptor.pdb"]
        APO_L["apo_ligand.pdb"]
        PRED_R["predicted_receptor.pdb"]
        PRED_L["predicted_ligand.pdb"]
    end
    
    SPLIT_XL --> SYS_DIR
    SYS_DIR --> HOLO_R
    SYS_DIR --> HOLO_L
    SYS_DIR --> APO_R
    SYS_DIR --> APO_L
    SYS_DIR --> PRED_R
    SYS_DIR --> PRED_L
```

**Diagram: test_set_pdbs/ directory structure for benchmark evaluation**

Each test set subdirectory (`pinder_s/`, `pinder_xl/`, `pinder_af2/`) contains system-specific folders with all available structure variants pre-transformed for direct use.

Sources: [README.md:162-177](), [README.md:244-252]()

## mappings/ Directory

The `mappings/` directory stores alignment and mapping information between different coordinate systems and databases:

```mermaid
graph TB
    MAPPINGS["mappings/"]
    
    subgraph "Mapping types"
        UNIPROT_MAP["UniProt mappings<br/>PDB ↔ UniProt coordinates"]
        ASSEMBLY_MAP["Assembly information<br/>Original PDB assemblies"]
        CHAIN_MAP["Chain mappings<br/>Author chains ↔ Entity IDs"]
    end
    
    subgraph "File organization"
        PREFIX_DIR["&lt;two-char-prefix&gt;/"]
        SYSTEM_FILES["Mapping files per system"]
    end
    
    MAPPINGS --> UNIPROT_MAP
    MAPPINGS --> ASSEMBLY_MAP
    MAPPINGS --> CHAIN_MAP
    
    MAPPINGS --> PREFIX_DIR
    PREFIX_DIR --> SYSTEM_FILES
    
    UNIPROT_MAP -.->|"Used by"| RESIDUE_NUMBERING["Residue numbering<br/>alignment"]
    ASSEMBLY_MAP -.->|"Used by"| STRUCTURE_UTILS["Structure utilities<br/>and validation"]
```

**Diagram: mappings/ directory contents and usage**

These mappings enable:
- Conversion between PDB residue numbering and UniProt canonical sequences
- Tracking of biological assembly transformations
- Resolution of author chain IDs to structural entities

Sources: [README.md:173-176]()

## index.parquet File

The `index.parquet` file is the master index containing one row per PINDER entry. It defines which structures are available for each system and their dataset split assignment.

### Index Schema Overview

The index contains 58 fields organized into several categories:

```mermaid
graph LR
    INDEX["index.parquet"]
    
    subgraph "Core identification"
        ID_FIELDS["id<br/>pdb_id<br/>cluster_id"]
    end
    
    subgraph "Split assignment"
        SPLIT_FIELDS["split<br/>pinder_xl<br/>pinder_s<br/>pinder_af2"]
    end
    
    subgraph "Structure availability"
        STRUCT_FIELDS["holo_R/L<br/>apo_R/L<br/>predicted_R/L"]
    end
    
    subgraph "File paths"
        PATH_FIELDS["holo_R_pdb<br/>apo_R_pdbs<br/>predicted_R_pdb"]
    end
    
    subgraph "Quality metrics"
        QUAL_FIELDS["apo_R_quality<br/>apo_L_quality<br/>chain1_neff<br/>chain2_neff"]
    end
    
    subgraph "Annotations"
        ANNOT_FIELDS["contains_antibody<br/>contains_antigen<br/>contains_enzyme"]
    end
    
    INDEX --> ID_FIELDS
    INDEX --> SPLIT_FIELDS
    INDEX --> STRUCT_FIELDS
    INDEX --> PATH_FIELDS
    INDEX --> QUAL_FIELDS
    INDEX --> ANNOT_FIELDS
```

**Diagram: index.parquet field categories**

### Key Index Fields

| Field | Type | Description |
|-------|------|-------------|
| `id` | string | Unique system identifier: `<receptor_spec>--<ligand_spec>` |
| `pdb_id` | string | Source PDB entry identifier |
| `split` | string | Dataset split: `train`, `val`, `test`, `invalid` |
| `cluster_id` | string | Interface cluster identifier for redundancy control |
| `holo_R` / `holo_L` | boolean | Availability of holo receptor/ligand structures |
| `apo_R` / `apo_L` | boolean | Availability of apo receptor/ligand structures |
| `predicted_R` / `predicted_L` | boolean | Availability of AlphaFold2 predicted structures |
| `apo_R_pdbs` / `apo_L_pdbs` | string | Pipe-delimited list of alternate apo structure IDs |
| `apo_R_quality` / `apo_L_quality` | string | Apo pairing quality: `high`, `low`, or empty string |

The index is accessed via the `PinderLoader` class and `get_index()` function from [src/pinder-core/pinder/core/index/utils.py]().

Sources: [README.md:173-176](), [README.md:529-567]()

## metadata.parquet File

The `metadata.parquet` file provides detailed structural and biophysical annotations with 45+ fields per entry:

```mermaid
graph TB
    METADATA["metadata.parquet"]
    
    subgraph "Experimental details"
        EXP["method<br/>resolution<br/>date<br/>release_date"]
    end
    
    subgraph "Structural properties"
        STRUCT["length1/2<br/>length_resolved_1/2<br/>number_of_components_1/2<br/>num_atom_types"]
    end
    
    subgraph "Interface characteristics"
        INTERFACE["buried_sasa<br/>intermolecular_contacts<br/>n_residue_pairs<br/>planarity<br/>link_density"]
    end
    
    subgraph "Contact types"
        CONTACTS["charged_charged_contacts<br/>charged_polar_contacts<br/>polar_polar_contacts<br/>apolar_apolar_contacts"]
    end
    
    subgraph "Quality indicators"
        QUALITY["interface_atom_gaps_4A/8A<br/>missing_interface_residues_4A/8A<br/>max_var_1/2"]
    end
    
    subgraph "Biological relevance"
        BIO["label<br/>probability<br/>complex_type<br/>oligomeric_count"]
    end
    
    subgraph "Domain annotations"
        DOMAIN["ECOD_names_R/L<br/>entity_id_R/L<br/>pdb_strand_id_R/L"]
    end
    
    METADATA --> EXP
    METADATA --> STRUCT
    METADATA --> INTERFACE
    METADATA --> CONTACTS
    METADATA --> QUALITY
    METADATA --> BIO
    METADATA --> DOMAIN
```

**Diagram: metadata.parquet field categories and relationships**

### Critical Metadata Fields

| Field | Type | Description |
|-------|------|-------------|
| `resolution` | float | Experimental resolution in Ångströms |
| `method` | string | Structure determination method (X-RAY, CRYO-EM, NMR) |
| `label` | string | PRODIGY-cryst classification: biological vs crystal contact |
| `probability` | float | Probability of being a true biological complex |
| `buried_sasa` | float | Buried solvent-accessible surface area (Ų) |
| `intermolecular_contacts` | int | Number of inter-chain contacts (≤5Å cutoff) |
| `interface_atom_gaps_4A` | int | Interface atoms within 4Å of structural gaps |
| `num_atom_types` | int | Number of distinct atom types (detects CA-only structures) |
| `planarity` | float | Interface shape complementarity measure |

The metadata is accessed via `system.metadata` on any `PinderSystem` instance or directly through `get_metadata()` from [src/pinder-core/pinder/core/index/utils.py]().

Sources: [README.md:173-177](), [README.md:569-624]()

## Data Flow and Access Patterns

The following diagram illustrates how the different dataset components are accessed during typical usage:

```mermaid
graph TB
    USER["User Code"]
    
    subgraph "High-level APIs"
        LOADER["PinderLoader"]
        SYSTEM["PinderSystem"]
    end
    
    subgraph "Index layer"
        GET_INDEX["get_index()"]
        GET_META["get_metadata()"]
        INDEX_FILE["index.parquet"]
        META_FILE["metadata.parquet"]
    end
    
    subgraph "Structure files"
        PDBS_DIR["pdbs/&lt;prefix&gt;/"]
        TEST_DIR["test_set_pdbs/&lt;split&gt;/"]
    end
    
    subgraph "Mapping files"
        MAP_DIR["mappings/&lt;prefix&gt;/"]
    end
    
    USER --> LOADER
    LOADER --> GET_INDEX
    LOADER --> SYSTEM
    
    SYSTEM --> GET_META
    SYSTEM --> PDBS_DIR
    SYSTEM --> TEST_DIR
    SYSTEM --> MAP_DIR
    
    GET_INDEX --> INDEX_FILE
    GET_META --> META_FILE
```

**Diagram: Data access patterns through the PINDER API**

Sources: [src/pinder-core/pinder/core/loader/loader.py](), [src/pinder-core/pinder/core/index/system.py](), [src/pinder-core/pinder/core/index/utils.py]()

## Download and Synchronization

The dataset can be downloaded using command-line utilities:

| Command | Purpose | Options |
|---------|---------|---------|
| `pinder_download` | Download complete dataset | `--pinder_base_dir`, `--pinder_release`, `--skip_inflation` |
| `pinder_update_index` | Sync index/metadata only | `--pinder_base_dir`, `--pinder_release` |
| `pinder_sync_data` | Sync changed structure files | `--pinder_base_dir`, `--pinder_release` |

These utilities handle downloading from Google Cloud Storage (`gs://pinder` bucket) and automatic decompression.

Sources: [README.md:146-240]()

---
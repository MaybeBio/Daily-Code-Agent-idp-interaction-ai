# mmCIF Processing

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [tests/data/test_get_data.py](tests/data/test_get_data.py)

</details>



## Purpose and Scope

This document describes the mmCIF (macromolecular Crystallographic Information File) processing system within the PINDER data generation pipeline. The mmCIF processing module is responsible for parsing raw structure files downloaded from RCSB PDB, generating biological assemblies, extracting entity information, identifying interacting chains, and producing structured output files for downstream dataset construction.

The mmCIF processing system is located in the `pinder.data.get_data` module and is a critical component of the dataset generation workflow. For information about how mmCIF files are acquired from RCSB, see [Data Acquisition from RCSB](#5.1). For details on task orchestration and parallel processing, see [Task Orchestration](#5.2). For information about GraphQL-based annotation fetching, see [RCSB Annotations](#5.4).

---

## mmCIF Processing Pipeline Overview

The following diagram shows the complete mmCIF processing pipeline, from input files to structured output:

```mermaid
graph TB
    subgraph "Input Files"
        MMCIF["mmCIF File<br/>pdb_{pdb_id}_xyz-enrich.cif.gz"]
    end
    
    subgraph "Core Processing Steps"
        READ["read_mmcif_file()<br/>Parse mmCIF format"]
        BIO_ASM["generate_bio_assembly()<br/>Create biological assembly"]
        ENTITIES["get_entities()<br/>Extract entity info"]
        SEQ_MAP["sequence_mapping()<br/>Map PDB to UniProt"]
        CHAINS["get_structure_chains()<br/>Identify chains"]
        INTERACT["get_interacting_chains()<br/>Contact analysis"]
        META["get_metadata()<br/>Extract metadata"]
    end
    
    subgraph "Processing Functions"
        FILTER["filter_amino_acids()<br/>Protein-only filtering"]
        UNIPROT["infer_uniprot_from_mapping()<br/>Determine UniProt ID"]
        SASA["Buried SASA calculation<br/>Interface analysis"]
    end
    
    subgraph "Output Files"
        ASSEMBLY_CIF["{pdb_id}-assembly.cif<br/>Biological assembly"]
        METADATA_TSV["{pdb_id}-metadata.tsv<br/>Structure metadata"]
        ENTITY_MAP["Entity mapping DataFrame<br/>chain_id → entity_id"]
        INTERACTION_DF["Interacting chains DataFrame<br/>Contact residues & SASA"]
    end
    
    MMCIF --> READ
    READ --> BIO_ASM
    READ --> META
    
    BIO_ASM --> FILTER
    BIO_ASM --> ENTITIES
    BIO_ASM --> ASSEMBLY_CIF
    
    FILTER --> CHAINS
    ENTITIES --> SEQ_MAP
    SEQ_MAP --> UNIPROT
    
    CHAINS --> INTERACT
    UNIPROT --> INTERACT
    
    INTERACT --> SASA
    SASA --> INTERACTION_DF
    
    META --> METADATA_TSV
    BIO_ASM --> ENTITY_MAP
    
    style READ fill:#e1f5ff
    style BIO_ASM fill:#fff4e1
    style INTERACT fill:#ffe1f5
    style ASSEMBLY_CIF fill:#e1ffe1
```

**Sources:** [tests/data/test_get_data.py:217-240](), [tests/data/test_get_data.py:242-267](), [tests/data/test_get_data.py:348-362]()

---

## Reading and Parsing mmCIF Files

The `read_mmcif_file()` function parses compressed mmCIF files downloaded from RCSB PDB. These files follow the naming convention `pdb_{pdb_id}_xyz-enrich.cif.gz` and contain the complete structural information including atomic coordinates, metadata, and biological assembly instructions.

```mermaid
graph LR
    subgraph "Input Format"
        GZ_CIF["pdb_00007cm8_xyz-enrich.cif.gz<br/>Compressed mmCIF"]
    end
    
    subgraph "Parsing Process"
        READ_FUNC["read_mmcif_file()<br/>Decompress & parse"]
        PDBX["CIFFile object<br/>biotite.structure.io.pdbx"]
    end
    
    subgraph "Extracted Data"
        ATOMS["Atomic coordinates<br/>xyz positions"]
        CHAINS["Chain information<br/>asym_id mapping"]
        ASSEMBLY["Assembly instructions<br/>Transformations"]
        METADATA["Structure metadata<br/>Resolution, method, date"]
    end
    
    GZ_CIF --> READ_FUNC
    READ_FUNC --> PDBX
    PDBX --> ATOMS
    PDBX --> CHAINS
    PDBX --> ASSEMBLY
    PDBX --> METADATA
    
    style READ_FUNC fill:#e1f5ff
    style PDBX fill:#fff4e1
```

The `CIFFile` object from biotite provides access to all mmCIF data categories, including `_atom_site`, `_pdbx_struct_assembly`, `_entity_poly`, and others.

**Sources:** [tests/data/test_get_data.py:114](), [tests/data/test_get_data.py:168]()

---

## Biological Assembly Generation

The `generate_bio_assembly()` function transforms the asymmetric unit structure into the biologically relevant assembly using transformation matrices defined in the mmCIF file. This step is critical because crystallographic structures may contain only a portion of the functional biological unit.

```mermaid
graph TB
    subgraph "Input"
        ASYMMETRIC["Asymmetric Unit<br/>Original mmCIF structure"]
    end
    
    subgraph "Assembly Generation"
        PARSE_ASSEMBLY["Parse assembly instructions<br/>_pdbx_struct_assembly"]
        APPLY_TRANSFORM["Apply transformation matrices<br/>Rotation + translation"]
        CHAIN_RENAME["Rename chains<br/>asym_id → chain_id (e.g., A1, A2)"]
    end
    
    subgraph "Entity Mapping"
        BUILD_MAP["Build entity mapping<br/>chain_id → entity_id"]
        ASYM_MAP["Track asym_id<br/>Original chain identifier"]
    end
    
    subgraph "Output"
        BIO_ASSEMBLY["AtomArray<br/>Biological assembly structure"]
        ENTITY_DF["Entity mapping DataFrame<br/>chain_id, asym_id, entity_id"]
    end
    
    ASYMMETRIC --> PARSE_ASSEMBLY
    PARSE_ASSEMBLY --> APPLY_TRANSFORM
    APPLY_TRANSFORM --> CHAIN_RENAME
    CHAIN_RENAME --> BUILD_MAP
    CHAIN_RENAME --> BIO_ASSEMBLY
    BUILD_MAP --> ASYM_MAP
    ASYM_MAP --> ENTITY_DF
    
    style APPLY_TRANSFORM fill:#fff4e1
    style BIO_ASSEMBLY fill:#e1ffe1
    style ENTITY_DF fill:#e1f5ff
```

The generated biological assembly may contain significantly more or fewer chains than the asymmetric unit. For example:
- PDB 7KMX: 14 asym_id chains expand to 840 assembly chains
- PDB 7BDU: 7 asym_id chains remain as 7 assembly chains
- PDB 2A79: 13 asym_id chains expand to 34 assembly chains

**Sources:** [tests/data/test_get_data.py:348-362](), [tests/data/test_get_data.py:116]()

---

## Entity Extraction and Sequence Mapping

Entity extraction identifies distinct molecular entities (proteins) within the structure and maps them to their sequences and UniProt accessions.

```mermaid
graph TB
    subgraph "Entity Identification"
        GET_ENT["get_entities(pdbx_file, entry_id)<br/>Extract from _entity_poly"]
        ENT_DF["Entity DataFrame<br/>entity_id, asym_id, pdb_strand_id"]
    end
    
    subgraph "Sequence Mapping"
        SEQ_MAP["sequence_mapping(pdbx_file, entity_id)<br/>PDB to UniProt mapping"]
        RES_MAP["Residue-level mapping<br/>pdb_num → uniprot_num"]
        RESOLVE["Track resolved residues<br/>Missing vs present"]
    end
    
    subgraph "Chain Association"
        ITER_ASYM["Iterate asym_id values<br/>For each entity"]
        MATCH_CHAIN["Match to assembly chains<br/>asym_id prefix matching"]
        UNIPROT["infer_uniprot_from_mapping()<br/>Determine UniProt accession"]
    end
    
    subgraph "Output"
        CHAIN_MAP["chain_mapping Dictionary<br/>chain → mapping DataFrame"]
        ENRICHED_ENT["Enriched entity DataFrame<br/>+ chain, uniprot, length_resolved"]
    end
    
    GET_ENT --> ENT_DF
    ENT_DF --> SEQ_MAP
    SEQ_MAP --> RES_MAP
    RES_MAP --> RESOLVE
    
    ENT_DF --> ITER_ASYM
    ITER_ASYM --> MATCH_CHAIN
    RES_MAP --> UNIPROT
    
    MATCH_CHAIN --> CHAIN_MAP
    UNIPROT --> ENRICHED_ENT
    RESOLVE --> ENRICHED_ENT
    
    style SEQ_MAP fill:#fff4e1
    style UNIPROT fill:#e1f5ff
    style ENRICHED_ENT fill:#e1ffe1
```

The sequence mapping process creates a residue-by-residue correspondence between PDB numbering and UniProt canonical sequences, tracking which residues are structurally resolved. This is essential for later alignment and superposition operations.

**Sources:** [tests/data/test_get_data.py:118-149]()

---

## Interacting Chain Identification

The `get_interacting_chains()` function performs contact analysis to identify which chain pairs form protein-protein interfaces. This is a core step in determining which structures qualify for inclusion in the PINDER dataset.

### Contact Detection Algorithm

```mermaid
graph TB
    subgraph "Input"
        BIO_ASM["Biological assembly<br/>AtomArray"]
        ENTITIES["Entity DataFrame<br/>Chain metadata"]
        PARAMS["Parameters<br/>contact_threshold (5Å)<br/>backbone_only (False)"]
    end
    
    subgraph "Atom Selection"
        FILTER_TYPE["Filter by atom type<br/>CA only if backbone_only"]
        CHAIN_PAIRS["Generate chain pairs<br/>All unique combinations"]
    end
    
    subgraph "Distance Calculation"
        DIST_CALC["Calculate pairwise distances<br/>Between chain atoms"]
        THRESHOLD["Apply threshold<br/>distances ≤ contact_threshold"]
        COUNT_PAIRS["Count atom pairs<br/>Below threshold"]
    end
    
    subgraph "Residue Analysis"
        MAP_RES["Map atoms to residues<br/>Unique residue IDs"]
        COUNT_RES_PAIRS["Count residue pairs<br/>In contact"]
        TOTAL_RES["Count total residues<br/>In interface"]
        RES_LIST["List contacting residues<br/>Comma-separated"]
    end
    
    subgraph "SASA Calculation"
        CALC_SASA["Calculate buried SASA<br/>Interface area"]
    end
    
    subgraph "Output"
        INTERACT_DF["Interacting chains DataFrame<br/>chain_1, chain_2, asym_id_1, asym_id_2<br/>n_residue_pairs, n_residues<br/>buried_sasa, residue lists"]
    end
    
    BIO_ASM --> FILTER_TYPE
    ENTITIES --> CHAIN_PAIRS
    PARAMS --> FILTER_TYPE
    
    FILTER_TYPE --> DIST_CALC
    CHAIN_PAIRS --> DIST_CALC
    PARAMS --> THRESHOLD
    
    DIST_CALC --> THRESHOLD
    THRESHOLD --> COUNT_PAIRS
    COUNT_PAIRS --> MAP_RES
    
    MAP_RES --> COUNT_RES_PAIRS
    MAP_RES --> TOTAL_RES
    MAP_RES --> RES_LIST
    
    COUNT_RES_PAIRS --> INTERACT_DF
    TOTAL_RES --> INTERACT_DF
    RES_LIST --> INTERACT_DF
    
    BIO_ASM --> CALC_SASA
    CALC_SASA --> INTERACT_DF
    
    style DIST_CALC fill:#fff4e1
    style CALC_SASA fill:#e1f5ff
    style INTERACT_DF fill:#e1ffe1
```

### Contact Analysis Output

The `get_interacting_chains()` function returns a DataFrame with the following fields:

| Field | Description |
|-------|-------------|
| `chain_1` | First chain identifier from assembly |
| `chain_2` | Second chain identifier from assembly |
| `asym_id_1` | Asymmetric unit ID for chain 1 |
| `asym_id_2` | Asymmetric unit ID for chain 2 |
| `pdb_strand_id_1` | Author chain ID for chain 1 |
| `pdb_strand_id_2` | Author chain ID for chain 2 |
| `entry_id` | PDB ID |
| `n_residue_pairs` | Number of residue pairs in contact |
| `n_residues` | Total unique residues in interface |
| `buried_sasa` | Buried solvent accessible surface area (Ų) |
| `chain_1_residues` | Comma-separated list of contacting residue numbers |
| `chain_2_residues` | Comma-separated list of contacting residue numbers |

### Parameter Effects

The contact detection behavior varies based on parameters:

| `contact_threshold` | `backbone_only` | Effect |
|---------------------|-----------------|--------|
| 5.0 Å | False | Default configuration - all heavy atoms, tight contacts |
| 10.0 Å | False | Expanded interface - captures more distant interactions |
| 5.0 Å | True | Backbone only - Cα atoms, stricter interface definition |
| 10.0 Å | True | Relaxed backbone - broader interface region |

**Sources:** [tests/data/test_get_data.py:40-91](), [tests/data/test_get_data.py:93-163]()

---

## Metadata Extraction

The `get_metadata()` function extracts structural and experimental metadata from the mmCIF file that characterizes the quality and type of the structure.

```mermaid
graph LR
    subgraph "mmCIF Categories"
        EXPTL["_exptl<br/>Experimental method"]
        REFINE["_refine<br/>Resolution"]
        DATES["_pdbx_database_status<br/>Dates"]
        ASSEMBLY_META["_pdbx_struct_assembly<br/>Assembly info"]
        OLIGO["_pdbx_struct_assembly_prop<br/>Oligomeric state"]
    end
    
    subgraph "Extracted Fields"
        METHOD["method<br/>X-RAY DIFFRACTION, etc."]
        RES["resolution<br/>Angstroms"]
        DATE["date, release_date<br/>YYYY-MM-DD"]
        ASM_DETAILS["assembly_details<br/>Assembly type"]
        OLIGO_STATE["oligomeric_details<br/>oligomeric_count"]
    end
    
    subgraph "Output"
        META_DF["Metadata DataFrame<br/>Single row per structure"]
    end
    
    EXPTL --> METHOD
    REFINE --> RES
    DATES --> DATE
    ASSEMBLY_META --> ASM_DETAILS
    OLIGO --> OLIGO_STATE
    
    METHOD --> META_DF
    RES --> META_DF
    DATE --> META_DF
    ASM_DETAILS --> META_DF
    OLIGO_STATE --> META_DF
    
    style META_DF fill:#e1ffe1
```

Example metadata record:

| Field | Value |
|-------|-------|
| `entry_id` | 7CM8 |
| `method` | X-RAY DIFFRACTION |
| `date` | 2020-07-25 |
| `release_date` | 2020-08-19 |
| `resolution` | 1.90 |
| `assembly` | 1 |
| `assembly_details` | author_and_software_defined_assembly |
| `oligomeric_details` | dimeric |
| `oligomeric_count` | 2 |

**Sources:** [tests/data/test_get_data.py:165-189](), [tests/data/test_get_data.py:191-216]()

---

## Main Processing Function

The `process_mmcif()` function orchestrates the complete processing pipeline, coordinating all the steps described above.

```mermaid
graph TB
    subgraph "Input"
        MMCIF_PATH["Path to mmCIF file<br/>pdb_{pdb_id}_xyz-enrich.cif.gz"]
    end
    
    subgraph "process_mmcif() Workflow"
        EXTRACT_ID["Extract pdb_id<br/>From filename"]
        READ["read_mmcif_file()<br/>Load mmCIF"]
        GEN_ASM["generate_bio_assembly()<br/>Create assembly & mapping"]
        GET_META["get_metadata()<br/>Extract metadata"]
        FILTER_AA["filter_amino_acids()<br/>Protein-only"]
        GET_ENT["get_entities()<br/>Entity information"]
        GET_CHAINS["get_structure_chains()<br/>Chain identifiers"]
        SEQ_MAPPING["sequence_mapping()<br/>Per entity & asym_id"]
        INFER_UNIPROT["infer_uniprot_from_mapping()<br/>UniProt accession"]
        BUILD_CHAIN_MAP["Build chain_mapping dict<br/>Chain → mapping"]
        GET_INTERACT["get_interacting_chains()<br/>Contact analysis"]
        CLASSIFY["Classify complex_type<br/>homomer vs heteromer"]
        SET_STATUS["Set status<br/>complete/incomplete"]
    end
    
    subgraph "Output Files"
        ASSEMBLY_OUT["{pdb_id}-assembly.cif<br/>Biological assembly"]
        METADATA_OUT["{pdb_id}-metadata.tsv<br/>Structure metadata"]
        CHAINS_OUT["Optional: Chain-level data<br/>For downstream processing"]
    end
    
    MMCIF_PATH --> EXTRACT_ID
    EXTRACT_ID --> READ
    READ --> GEN_ASM
    READ --> GET_META
    GEN_ASM --> FILTER_AA
    GEN_ASM --> ASSEMBLY_OUT
    FILTER_AA --> GET_ENT
    FILTER_AA --> GET_CHAINS
    GET_ENT --> SEQ_MAPPING
    SEQ_MAPPING --> INFER_UNIPROT
    GET_CHAINS --> BUILD_CHAIN_MAP
    INFER_UNIPROT --> BUILD_CHAIN_MAP
    BUILD_CHAIN_MAP --> GET_INTERACT
    GET_INTERACT --> CLASSIFY
    GET_META --> CLASSIFY
    CLASSIFY --> SET_STATUS
    SET_STATUS --> METADATA_OUT
    BUILD_CHAIN_MAP --> CHAINS_OUT
    
    style process_mmcif fill:#fff4e1
    style ASSEMBLY_OUT fill:#e1ffe1
    style METADATA_OUT fill:#e1ffe1
```

The `process_mmcif()` function enriches the metadata with derived fields:

- `complex_type`: Determined by comparing UniProt accessions - "homomer" if all chains have the same UniProt ID, "heteromer" otherwise
- `status`: Set to "complete" if processing succeeds without errors, otherwise may indicate issues

**Sources:** [tests/data/test_get_data.py:217-240]()

---

## Output File Formats

### Assembly CIF File

The `{pdb_id}-assembly.cif` file contains the biological assembly with renamed chain identifiers. Chains are named using a prefix system where the asym_id from the asymmetric unit is preserved with numeric suffixes (e.g., `A1`, `A2` for multiple copies of asym_id `A`).

**Chain naming example for PDB 7CM8:**
```
Asymmetric unit:    A
Biological assembly: A1, A2  (dimer)
```

The assembly CIF can be converted to PDB format with specific chain selections using the `atoms.cif_to_pdb()` function.

**Sources:** [tests/data/test_get_data.py:242-267]()

### Metadata TSV File

The `{pdb_id}-metadata.tsv` file is a tab-separated values file containing a single row with the metadata fields extracted by `get_metadata()` plus derived fields from `process_mmcif()`:

```
entry_id    method    date    release_date    resolution    assembly    ...
7cm8        X-RAY...  2020... 2020-08-19     1.9           1          ...
```

**Sources:** [tests/data/test_get_data.py:227-236]()

---

## Special Cases and Edge Cases

### Large Oligomeric Assemblies

Some structures produce extremely large biological assemblies:

| PDB ID | Asym IDs | Assembly Chains | Entity IDs |
|--------|----------|-----------------|------------|
| 7KMX | 14 | 840 | 2 |
| 6RW4 | 125 | 125 | 42 |
| 2Y26 | 40 | 120 | 2 |

These structures require significant memory and processing time during contact analysis.

**Sources:** [tests/data/test_get_data.py:284-362]()

### Monomer Processing

Structures containing only a single protein chain are processed but may not produce interacting chain data since no protein-protein interface exists. The `process_mmcif()` function handles these cases gracefully.

**Sources:** [tests/data/test_get_data.py:364-368]()

### Invalid or Incomplete Structures

Some mmCIF files may lack required data or contain invalid structures. The `process_mmcif()` function may raise exceptions for these cases, which are caught during task orchestration to skip problematic entries.

**Sources:** [tests/data/test_get_data.py:237-240]()

---

## Integration with Dataset Pipeline

The mmCIF processing outputs feed directly into the dataset construction workflow:

1. **Assembly CIF files** are used to create monomer structures (receptor/ligand) by extracting individual chains
2. **Metadata TSV files** provide filtering criteria (resolution, method, oligomeric state)
3. **Entity mappings** enable pairing of apo structures with holo complexes based on UniProt accessions
4. **Interacting chain data** defines which chain pairs constitute valid protein-protein interfaces for the dataset

The processed data is organized in the file system by two-character PDB codes (e.g., `cm/pdb_00007cm8/`) and subsequently indexed for efficient access by the data loading system (see [PinderLoader](#4.1)).

**Sources:** [tests/data/test_get_data.py:1-368]()

---
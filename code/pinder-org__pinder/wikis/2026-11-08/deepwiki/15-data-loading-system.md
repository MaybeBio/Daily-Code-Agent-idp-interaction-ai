# Data Loading System

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/pinder-core/pinder/core/loader/dataset.py](src/pinder-core/pinder/core/loader/dataset.py)
- [src/pinder-core/pinder/core/loader/loader.py](src/pinder-core/pinder/core/loader/loader.py)
- [src/pinder-core/pinder/core/utils/timer.py](src/pinder-core/pinder/core/utils/timer.py)

</details>



The Data Loading System provides the infrastructure for loading PINDER protein systems and converting them into ML-ready data structures. This system bridges the gap between raw PDB files stored in the PINDER dataset and batched tensors or graphs consumable by PyTorch and PyTorch Geometric models.

The system is organized around three core components: **PinderLoader** for system iteration and filtering, **PinderDataset** for tensor-based workflows, and **PPIDataset** for graph neural network workflows. All three leverage a common filtering and transformation pipeline.

For details on the underlying data structures, see [Core Data Structures](#3). For information on specific filter implementations, see [Filtering System](#4.4). For evaluation workflows, see [Evaluation](#6).

## System Architecture

The data loading pipeline consists of multiple layers that progressively transform raw structural data into ML-ready batches:

```mermaid
graph TB
    subgraph Storage["Storage Layer"]
        INDEX["index.parquet<br/>2.4M+ entries"]
        METADATA["metadata.parquet<br/>Annotations"]
        PDBS["pdbs/<br/>Structural files"]
    end
    
    subgraph Orchestration["Orchestration Layer"]
        LOADER["PinderLoader<br/>System iteration & filtering"]
        MONOMER_SEL["select_monomer()<br/>Monomer selection logic"]
        CREATE_COMPLEX["_create_target_feature_complex()<br/>Complex construction"]
    end
    
    subgraph Filtering["Filtering Pipeline"]
        BASE_FILT["PinderFilterBase<br/>System-level filters"]
        SUB_FILT["PinderFilterSubBase<br/>Monomer-level filters"]
        STRUCT_FILT["StructureFilter<br/>Atomic-level filters"]
    end
    
    subgraph Datasets["Dataset Layer"]
        PINDER_DS["PinderDataset<br/>torch.utils.data.Dataset"]
        PPI_DS["PPIDataset<br/>torch_geometric.data.Dataset"]
    end
    
    subgraph DataLoaders["DataLoader Layer"]
        TORCH_DL["torch.utils.data.DataLoader<br/>collate_batch()"]
        GEO_DL["torch_geometric.loader.DataLoader"]
    end
    
    subgraph Output["Output Format"]
        TENSORS["Dict[str, Tensor]<br/>atom_coordinates<br/>atom_types<br/>residue_types"]
        HETERO["HeteroData<br/>ligand_atom<br/>receptor_atom<br/>edges"]
    end
    
    INDEX --> LOADER
    METADATA --> LOADER
    PDBS --> LOADER
    
    LOADER --> MONOMER_SEL
    MONOMER_SEL --> CREATE_COMPLEX
    
    BASE_FILT --> LOADER
    SUB_FILT --> LOADER
    STRUCT_FILT --> LOADER
    
    LOADER --> PINDER_DS
    LOADER --> PPI_DS
    
    PINDER_DS --> TORCH_DL
    PPI_DS --> GEO_DL
    
    TORCH_DL --> TENSORS
    GEO_DL --> HETERO
```

**Sources:** [src/pinder-core/pinder/core/loader/dataset.py:1-655](), [src/pinder-core/pinder/core/loader/loader.py:1-439]()

## Data Flow Pipeline

The following diagram illustrates how a single system progresses through the loading pipeline, from index entry to batched output:

```mermaid
graph LR
    subgraph Input["Input"]
        ROW["pd.Series<br/>Index row"]
    end
    
    subgraph PinderLoader["PinderLoader.__getitem__()"]
        LOAD_SYS["PinderSystem(row.id)"]
        APPLY_DIMER["apply_dimer_filters()<br/>base_filters + sub_filters"]
        SELECT["select_monomer()<br/>Choose apo/holo/pred"]
        CREATE["_create_target_feature_complex()<br/>Build complexes"]
        APPLY_STRUCT["apply_structure_filters()<br/>structure_filters"]
        TRANSFORM["Apply transforms<br/>structure_transforms_*"]
    end
    
    subgraph Output["Output"]
        TUPLE["(PinderSystem,<br/>Structure feature,<br/>Structure target)"]
    end
    
    ROW --> LOAD_SYS
    LOAD_SYS --> APPLY_DIMER
    APPLY_DIMER --> SELECT
    SELECT --> CREATE
    CREATE --> APPLY_STRUCT
    APPLY_STRUCT --> TRANSFORM
    TRANSFORM --> TUPLE
```

**Sources:** [src/pinder-core/pinder/core/loader/loader.py:357-405]()

### Pipeline Stages

The loading pipeline executes the following stages for each system:

| Stage | Component | Purpose | Can Filter? |
|-------|-----------|---------|-------------|
| **System Loading** | `PinderSystem(row.id)` | Load system from disk | No |
| **Dimer Filtering** | `apply_dimer_filters()` | Apply base and sub filters | Yes |
| **Monomer Selection** | `select_monomer()` | Choose monomer types (apo/holo/pred) | No |
| **Complex Creation** | `_create_target_feature_complex()` | Build target and feature complexes | No |
| **Structure Filtering** | `apply_structure_filters()` | Apply structure-level filters | Yes |
| **Transformation** | `structure_transforms_*` | Apply augmentations | No |

**Sources:** [src/pinder-core/pinder/core/loader/loader.py:357-405]()

## Core Components

### PinderLoader

The `PinderLoader` class is the central orchestrator of the loading pipeline. It implements the iterator protocol and provides indexed access to filtered and transformed systems.

```mermaid
graph TB
    subgraph Initialization["PinderLoader Initialization"]
        PARAMS["Constructor Parameters"]
        SPLIT["split: str<br/>train/val/test"]
        IDS["ids: list[str]<br/>Specific system IDs"]
        FILTERS["Filters<br/>base/sub/structure"]
        MONOMER_PRIORITY["monomer_priority: str<br/>apo/holo/pred/random"]
        QUERIES["index_query<br/>metadata_query"]
    end
    
    subgraph State["Internal State"]
        INDEX_DF["self.index: pd.DataFrame<br/>Filtered index"]
        METADATA_DF["self.metadata: pd.DataFrame<br/>Filtered metadata"]
        BASE_FILTERS["self.base_filters<br/>list[PinderFilterBase]"]
        SUB_FILTERS["self.sub_filters<br/>list[PinderFilterSubBase]"]
        STRUCT_FILTERS["self.structure_filters<br/>list[StructureFilter]"]
        TRANSFORMS["self.structure_transforms_*<br/>list[StructureTransform]"]
    end
    
    subgraph Methods["Key Methods"]
        LEN["__len__()<br/>Returns len(self.index)"]
        GETITEM["__getitem__(idx)<br/>Returns (system, feature, target)"]
        ITER["__iter__()<br/>Yields tuples"]
        APPLY_DIMER_FILTERS["apply_dimer_filters()<br/>Static method"]
        APPLY_STRUCT_FILTERS["apply_structure_filters()<br/>Instance method"]
    end
    
    PARAMS --> INDEX_DF
    PARAMS --> METADATA_DF
    PARAMS --> BASE_FILTERS
    PARAMS --> SUB_FILTERS
    PARAMS --> STRUCT_FILTERS
    PARAMS --> TRANSFORMS
    
    INDEX_DF --> LEN
    INDEX_DF --> GETITEM
    BASE_FILTERS --> APPLY_DIMER_FILTERS
    SUB_FILTERS --> APPLY_DIMER_FILTERS
    STRUCT_FILTERS --> APPLY_STRUCT_FILTERS
    TRANSFORMS --> GETITEM
```

**Key Parameters:**

- `split` (str): Dataset split to load (train/val/test)
- `monomer_priority` (str): Strategy for selecting monomers (holo/apo/pred/random/random_mixed)
- `base_filters` (list[PinderFilterBase]): System-level filters
- `sub_filters` (list[PinderFilterSubBase]): Monomer-level filters  
- `structure_filters` (list[StructureFilter]): Structure-level filters
- `crop_equal_monomer_shapes` (bool): Whether to crop target and feature to same shape
- `fallback_to_holo` (bool): Whether to fall back to holo when preferred monomer unavailable
- `max_load_attempts` (int): Maximum retries when system fails filters (default: 10)

**Sources:** [src/pinder-core/pinder/core/loader/loader.py:192-265]()

### Monomer Selection Logic

The `select_monomer()` function determines which monomer type (holo/apo/pred) to use for receptor and ligand based on availability and priority:

```mermaid
graph TD
    START["select_monomer(row, monomer_priority)"]
    GET_AVAILABLE["get_available_monomers(row)<br/>Returns available types per side"]
    
    PRIORITY_CHECK{"monomer_priority?"}
    RANDOM_MIXED["random_mixed<br/>Random for R and L separately"]
    RANDOM["random<br/>Same random for both"]
    SPECIFIC["apo/holo/pred<br/>Use specified type"]
    
    PAIRED_CHECK{"Type available<br/>for both sides?"}
    FALLBACK_CHECK{"fallback_to_holo?"}
    USE_HOLO["Use holo for both"]
    ERROR["Raise ValueError"]
    
    APO_CHECK{"Selection == apo<br/>and not canonical?"}
    GET_ALT["get_alternate_apo_codes()<br/>Get non-canonical apo"]
    RANDOM_ALT["Randomly select alternate"]
    
    RESULT["Return dict[str, str]<br/>R and L selections"]
    
    START --> GET_AVAILABLE
    GET_AVAILABLE --> PRIORITY_CHECK
    
    PRIORITY_CHECK -->|"random_mixed"| RANDOM_MIXED
    PRIORITY_CHECK -->|"random"| RANDOM
    PRIORITY_CHECK -->|"apo/holo/pred"| SPECIFIC
    
    RANDOM_MIXED --> APO_CHECK
    RANDOM --> APO_CHECK
    
    SPECIFIC --> PAIRED_CHECK
    PAIRED_CHECK -->|"Yes"| APO_CHECK
    PAIRED_CHECK -->|"No"| FALLBACK_CHECK
    FALLBACK_CHECK -->|"True"| USE_HOLO
    FALLBACK_CHECK -->|"False"| ERROR
    USE_HOLO --> RESULT
    
    APO_CHECK -->|"Yes"| GET_ALT
    APO_CHECK -->|"No"| RESULT
    GET_ALT --> RANDOM_ALT
    RANDOM_ALT --> RESULT
```

**Monomer Priority Modes:**

| Mode | Behavior | Use Case |
|------|----------|----------|
| `holo` | Always use bound structures | Baseline performance, fastest loading |
| `apo` | Use unbound structures when available | Realistic docking scenario |
| `pred` | Use AlphaFold predictions when available | Evaluate prediction quality |
| `random` | Randomly select same type for R and L | Diverse training with consistent pairs |
| `random_mixed` | Randomly select types independently | Maximum diversity |

**Sources:** [src/pinder-core/pinder/core/loader/loader.py:81-149](), [src/pinder-core/pinder/core/loader/loader.py:23]()

## Dataset Implementations

### PinderDataset

`PinderDataset` is a PyTorch `Dataset` implementation that yields dictionary-based tensor structures suitable for standard neural networks:

```mermaid
graph TB
    subgraph Initialization["PinderDataset"]
        LOADER_INIT["self.loader = PinderLoader(...)"]
        TRANSFORM["self.transform<br/>Structure -> Tensor"]
        TARGET_TRANSFORM["self.target_transform<br/>Structure -> Tensor"]
    end
    
    subgraph GetItem["__getitem__(idx)"]
        GET_FROM_LOADER["system, feature, target = self.loader[idx]"]
        APPLY_TRANSFORM["transform(feature)<br/>target_transform(target)"]
        BUILD_DICT["Build output dict"]
    end
    
    subgraph Output["Output Dictionary"]
        TARGET_COMPLEX["target_complex:<br/>dict[str, Tensor]"]
        FEATURE_COMPLEX["feature_complex:<br/>dict[str, Tensor]"]
        IDS["id: str<br/>sample_id: str<br/>target_id: str"]
    end
    
    subgraph TensorStructure["Tensor Dictionary Keys"]
        ATOM_COORDS["atom_coordinates<br/>(N_atoms, 3)"]
        ATOM_TYPES["atom_types<br/>(N_atoms,)"]
        ELEMENT_TYPES["element_types<br/>(N_atoms,)"]
        RES_COORDS["residue_coordinates<br/>(N_residues, 3)"]
        RES_TYPES["residue_types<br/>(N_residues,)"]
        RES_IDS["residue_ids<br/>(N_residues,)"]
        CHAIN_IDS["chain_ids<br/>(N_residues,)"]
    end
    
    LOADER_INIT --> GET_FROM_LOADER
    GET_FROM_LOADER --> APPLY_TRANSFORM
    APPLY_TRANSFORM --> BUILD_DICT
    BUILD_DICT --> TARGET_COMPLEX
    BUILD_DICT --> FEATURE_COMPLEX
    BUILD_DICT --> IDS
    
    TARGET_COMPLEX --> ATOM_COORDS
    TARGET_COMPLEX --> ATOM_TYPES
    TARGET_COMPLEX --> ELEMENT_TYPES
    TARGET_COMPLEX --> RES_COORDS
    TARGET_COMPLEX --> RES_TYPES
    TARGET_COMPLEX --> RES_IDS
    TARGET_COMPLEX --> CHAIN_IDS
```

**Default Transform:** The `structure2tensor_transform()` function converts a `Structure` object into a dictionary of tensors using `structure2tensor()` from the geodata module.

**Sources:** [src/pinder-core/pinder/core/loader/dataset.py:246-318](), [src/pinder-core/pinder/core/loader/dataset.py:55-66]()

### PPIDataset

`PPIDataset` is a PyTorch Geometric `Dataset` that creates heterogeneous graph representations (`HeteroData`) suitable for graph neural networks:

```mermaid
graph TB
    subgraph Initialization["PPIDataset"]
        NODE_TYPES["node_types:<br/>set[NodeRepresentation]"]
        MONOMER_SPEC["monomer1: str<br/>monomer2: str"]
        LOADER_INIT["self.loader = PinderLoader"]
        PARALLEL["parallel: bool<br/>max_workers: int"]
    end
    
    subgraph Processing["process() / process_parallel()"]
        ITERATE["Iterate over loader"]
        CREATE_PAIRED["PairedPDB.from_pinder_system()"]
        APPLY_PRE_TRANSFORM["Apply pre_transform"]
        SAVE_PT["torch.save(data, .pt file)"]
        TRACK_FILES["Track filenames"]
    end
    
    subgraph Storage["Processed Files"]
        PROCESSED_DIR["processed/<br/>*.pt files"]
        FILENAMES_TXT["filenames/filenames.txt"]
    end
    
    subgraph GetItem["get(idx)"]
        LOAD_PT["torch.load(filename.pt)"]
        PAIRED_PDB["Returns PairedPDB"]
    end
    
    subgraph HeteroData["HeteroData Structure"]
        LIGAND_ATOM["ligand_atom<br/>Node features"]
        RECEPTOR_ATOM["receptor_atom<br/>Node features"]
        LIGAND_RES["ligand_residue<br/>Node features"]
        RECEPTOR_RES["receptor_residue<br/>Node features"]
        EDGES["Edge indices<br/>k-NN connections"]
    end
    
    NODE_TYPES --> CREATE_PAIRED
    MONOMER_SPEC --> CREATE_PAIRED
    LOADER_INIT --> ITERATE
    
    ITERATE --> CREATE_PAIRED
    CREATE_PAIRED --> APPLY_PRE_TRANSFORM
    APPLY_PRE_TRANSFORM --> SAVE_PT
    SAVE_PT --> PROCESSED_DIR
    TRACK_FILES --> FILENAMES_TXT
    
    PROCESSED_DIR --> LOAD_PT
    LOAD_PT --> PAIRED_PDB
    
    PAIRED_PDB --> LIGAND_ATOM
    PAIRED_PDB --> RECEPTOR_ATOM
    PAIRED_PDB --> LIGAND_RES
    PAIRED_PDB --> RECEPTOR_RES
    PAIRED_PDB --> EDGES
```

**Key Features:**

- **Preprocessing:** Converts all systems to `.pt` files during initialization
- **Caching:** Reuses processed files across runs when `use_cache=True`
- **Parallel Processing:** Supports multiprocessing via `parallel=True`
- **Node Types:** Supports atom-level and residue-level representations
- **Edges:** Automatically adds k-NN edges between nodes when `add_edges=True`

**Sources:** [src/pinder-core/pinder/core/loader/dataset.py:320-623]()

## Filtering Architecture

The filtering system operates at three distinct levels:

```mermaid
graph TB
    subgraph Level1["Level 1: Base Filters"]
        BASE_DESC["System-level properties<br/>Applied to PinderSystem"]
        BASE_EXAMPLES["Examples:<br/>- Cluster size<br/>- Quality metrics<br/>- Missing structures"]
    end
    
    subgraph Level2["Level 2: Sub Filters"]
        SUB_DESC["Monomer-level properties<br/>Modifies PinderSystem"]
        SUB_EXAMPLES["Examples:<br/>- Sequence length<br/>- Resolution<br/>- Completeness"]
        SUB_ACTION["Can set monomers to None"]
    end
    
    subgraph Level3["Level 3: Structure Filters"]
        STRUCT_DESC["Atomic-level properties<br/>Applied to Structure"]
        STRUCT_EXAMPLES["Examples:<br/>- Atom count<br/>- Coordinate validity<br/>- Clashes"]
    end
    
    subgraph Application["Filter Application"]
        LOADER_APPLY["PinderLoader.__getitem__()"]
        STEP1["1. apply_dimer_filters()"]
        STEP2["2. select_monomer()"]
        STEP3["3. apply_structure_filters()"]
    end
    
    BASE_DESC --> BASE_EXAMPLES
    SUB_DESC --> SUB_EXAMPLES
    SUB_EXAMPLES --> SUB_ACTION
    STRUCT_DESC --> STRUCT_EXAMPLES
    
    LOADER_APPLY --> STEP1
    STEP1 --> STEP2
    STEP2 --> STEP3
    
    BASE_EXAMPLES -.->|"Input to"| STEP1
    SUB_EXAMPLES -.->|"Input to"| STEP1
    STRUCT_EXAMPLES -.->|"Input to"| STEP3
```

**Filter Behavior:**

| Filter Type | Input | Output | Modifies System? |
|-------------|-------|--------|------------------|
| `PinderFilterBase` | `PinderSystem` | `bool` | No |
| `PinderFilterSubBase` | `PinderSystem` | `PinderSystem` | Yes (can set monomers to None) |
| `StructureFilter` | `Structure` | `bool` | No |

**Sources:** [src/pinder-core/pinder/core/loader/loader.py:415-435](), [src/pinder-core/pinder/core/loader/loader.py:407-413]()

## Complex Construction

The `_create_target_feature_complex()` function builds the target (ground truth holo) and feature (potentially apo/pred) complexes:

```mermaid
graph TD
    START["_create_target_feature_complex()"]
    INPUTS["Inputs:<br/>system: PinderSystem<br/>selected_monomers: dict<br/>crop_equal_monomer_shapes: bool"]
    
    GET_HOLO["holo_R = system.aligned_holo_R<br/>holo_L = system.aligned_holo_L"]
    
    CHECK_ALT{"Selected monomer<br/>is PDB code?"}
    SET_APO_CODE["system.apo_receptor_pdb_code = code<br/>selected_monomers[side] = 'apo'"]
    
    GET_DECOY["decoy_R = getattr(system, monomer + '_receptor')<br/>decoy_L = getattr(system, monomer + '_ligand')"]
    
    CHECK_NULL{"decoy is None?"}
    FALLBACK{"fallback_to_holo?"}
    SET_HOLO["decoy = holo monomer"]
    RAISE_ERROR["Raise ValueError"]
    
    BOTH_HOLO{"Both selections<br/>are holo?"}
    CROP{"crop_equal_<br/>monomer_shapes?"}
    ALIGN_CROP["_align_monomers_with_mask()<br/>Crop to common atoms"]
    
    BUILD_COMPLEXES["target = holo_R + holo_L<br/>feature = decoy_R + decoy_L"]
    RETURN["Return (target, feature)"]
    
    START --> INPUTS
    INPUTS --> GET_HOLO
    GET_HOLO --> CHECK_ALT
    
    CHECK_ALT -->|"Yes"| SET_APO_CODE
    CHECK_ALT -->|"No"| GET_DECOY
    SET_APO_CODE --> GET_DECOY
    
    GET_DECOY --> CHECK_NULL
    CHECK_NULL -->|"Yes"| FALLBACK
    CHECK_NULL -->|"No"| BOTH_HOLO
    
    FALLBACK -->|"True"| SET_HOLO
    FALLBACK -->|"False"| RAISE_ERROR
    SET_HOLO --> BOTH_HOLO
    
    BOTH_HOLO -->|"No"| CROP
    BOTH_HOLO -->|"Yes"| BUILD_COMPLEXES
    
    CROP -->|"True"| ALIGN_CROP
    CROP -->|"False"| BUILD_COMPLEXES
    
    ALIGN_CROP --> BUILD_COMPLEXES
    BUILD_COMPLEXES --> RETURN
```

**Shape Alignment:** When `crop_equal_monomer_shapes=True`, the function uses `_align_monomers_with_mask()` to crop holo and decoy structures to matching atom sets based on sequence alignment. This ensures the target and feature complexes have identical shapes, which is required by many loss functions.

**Sources:** [src/pinder-core/pinder/core/loader/loader.py:152-189]()

## Batching and Collation

### Padding Strategy

Variable-sized protein structures are padded to a uniform size within each batch using the `pad_and_stack()` and `collate_complex()` functions:

```mermaid
graph TB
    subgraph Input["Batch Input"]
        LIST["list[dict[str, Tensor]]<br/>Variable sizes"]
    end
    
    subgraph CollateComplex["collate_complex()"]
        EXTRACT["Extract tensors:<br/>atom_coordinates<br/>atom_types<br/>residue_coordinates<br/>etc."]
        
        PAD_STACK["For each tensor type:<br/>pad_and_stack()"]
    end
    
    subgraph PadAndStack["pad_and_stack()"]
        FIND_MAX["Find max shape:<br/>shapes.max(dim=0)"]
        PAD_EACH["For each tensor:<br/>pad_to_max_length()"]
        STACK["torch.stack()"]
    end
    
    subgraph PadToMaxLength["pad_to_max_length()"]
        CALC_PAD["Calculate padding:<br/>max_length - current_length"]
        APPLY_PAD["torch.nn.functional.pad()"]
    end
    
    subgraph Output["Batch Output"]
        PADDED["dict[str, Tensor]<br/>Uniform shape (B, N, ...)"]
    end
    
    LIST --> EXTRACT
    EXTRACT --> PAD_STACK
    PAD_STACK --> FIND_MAX
    FIND_MAX --> PAD_EACH
    PAD_EACH --> CALC_PAD
    CALC_PAD --> APPLY_PAD
    APPLY_PAD --> STACK
    STACK --> PADDED
```

**Padding Values:**

| Tensor Type | Padding Value | Purpose |
|-------------|---------------|---------|
| `atom_coordinates` | -100 | Distinguishable from valid coordinates |
| `residue_coordinates` | -100 | Distinguishable from valid coordinates |
| `atom_types` | -1 | Invalid atom type index |
| `element_types` | -1 | Invalid element type index |
| `residue_types` | -1 | Invalid residue type index |
| `residue_ids` | -99 | Invalid residue ID |
| `chain_ids` | -1 | Invalid chain ID |

**Sources:** [src/pinder-core/pinder/core/loader/dataset.py:69-163](), [src/pinder-core/pinder/core/loader/dataset.py:165-202](), [src/pinder-core/pinder/core/loader/dataset.py:205-243]()

### Collation Functions

The `collate_batch()` function merges a list of dataset items into a single batch:

```mermaid
graph LR
    subgraph BatchInput["Batch Input"]
        ITEMS["list[dict]<br/>Each dict has:<br/>- id<br/>- sample_id<br/>- target_id<br/>- target_complex<br/>- feature_complex"]
    end
    
    subgraph Processing["collate_batch()"]
        SEPARATE["Separate into lists:<br/>ids<br/>sample_ids<br/>target_ids<br/>target_structures<br/>feature_structures"]
        
        COLLATE_TARGET["collate_complex(<br/>target_structures)"]
        COLLATE_FEATURE["collate_complex(<br/>feature_structures)"]
    end
    
    subgraph Output["Batch Output"]
        BATCH["dict:<br/>- target_complex: dict[str, Tensor]<br/>- feature_complex: dict[str, Tensor]<br/>- id: list[str]<br/>- sample_id: list[str]<br/>- target_id: list[str]"]
    end
    
    ITEMS --> SEPARATE
    SEPARATE --> COLLATE_TARGET
    SEPARATE --> COLLATE_FEATURE
    COLLATE_TARGET --> BATCH
    COLLATE_FEATURE --> BATCH
```

**Sources:** [src/pinder-core/pinder/core/loader/dataset.py:205-243]()

## DataLoader Construction

### PyTorch DataLoader

The `get_torch_loader()` function creates a standard PyTorch `DataLoader` for `PinderDataset`:

```python
def get_torch_loader(
    dataset: PinderDataset,
    batch_size: int = 2,
    shuffle: bool = True,
    sampler: Sampler | None = None,
    num_workers: int = 1,
    collate_fn: Callable = collate_batch,
    **kwargs
) -> DataLoader[PinderDataset]:
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=num_workers,
        sampler=sampler,
        collate_fn=collate_fn,
        **kwargs,
    )
```

**Key Parameters:**
- `collate_fn`: Defaults to `collate_batch()` which handles variable-sized structures
- `sampler`: Can provide custom sampling logic (e.g., diversity sampling)
- `num_workers`: Number of subprocesses for data loading

**Sources:** [src/pinder-core/pinder/core/loader/dataset.py:637-655]()

### PyTorch Geometric DataLoader

The `get_geo_loader()` function creates a PyTorch Geometric `DataLoader` for `PPIDataset`:

```python
def get_geo_loader(
    dataset: PPIDataset,
    batch_size: int = 2,
    num_workers: int = 1
) -> TorchGeoDataLoader:
    return TorchGeoDataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=True,
        drop_last=True,
        num_workers=num_workers,
    )
```

**Graph Batching:** PyTorch Geometric's `DataLoader` automatically handles batching of `HeteroData` objects by concatenating node features and adjusting edge indices with offsets.

**Sources:** [src/pinder-core/pinder/core/loader/dataset.py:625-634]()

## Usage Patterns

### Basic PinderDataset Example

```python
from pinder.core.loader.dataset import PinderDataset, get_torch_loader

# Create dataset
dataset = PinderDataset(
    split="train",
    monomer_priority="apo",
    fallback_to_holo=True,
)

# Create dataloader
loader = get_torch_loader(
    dataset,
    batch_size=4,
    shuffle=True,
    num_workers=2
)

# Iterate
for batch in loader:
    target = batch["target_complex"]  # Ground truth
    feature = batch["feature_complex"]  # Apo structures
    ids = batch["id"]  # System IDs
    # Train model...
```

### Basic PPIDataset Example

```python
from pinder.core.loader.dataset import PPIDataset, get_geo_loader
from pinder.core.loader.geodata import NodeRepresentation

# Create dataset
dataset = PPIDataset(
    node_types={NodeRepresentation.ATOM},
    split="train",
    monomer1="apo_receptor",
    monomer2="apo_ligand",
    add_edges=True,
    k=10,  # 10-NN edges
)

# Create dataloader
loader = get_geo_loader(dataset, batch_size=4, num_workers=2)

# Iterate
for batch in loader:
    receptor_feats = batch["receptor_atom"]
    ligand_feats = batch["ligand_atom"]
    # Train GNN model...
```

### Advanced: Custom Filtering and Transforms

```python
from pinder.core.loader.dataset import PinderDataset
from pinder.core.loader.filters import PinderFilterBase, StructureFilter
from pinder.core.loader.transforms import SelectAtomTypes

# Define custom filters
class MinAtomFilter(StructureFilter):
    def __init__(self, min_atoms: int):
        self.min_atoms = min_atoms
    
    def __call__(self, structure: Structure) -> bool:
        return len(structure.atom_array) >= self.min_atoms

# Create dataset with filters and transforms
dataset = PinderDataset(
    split="train",
    structure_filters=[MinAtomFilter(min_atoms=100)],
    structure_transforms_feature=[
        SelectAtomTypes(atom_types=["CA"])  # Use only alpha carbons
    ],
    monomer_priority="random_mixed",
    crop_equal_monomer_shapes=True,
)
```

## Retry Logic and Error Handling

When a system fails filters or cannot be loaded, `PinderLoader` automatically retries with a random index:

```mermaid
graph TD
    START["__getitem__(idx)"]
    INIT["valid_idx = False<br/>attempts = 0"]
    
    LOOP_START{"valid_idx or<br/>attempts >= max?"}
    INCREMENT["attempts += 1"]
    LOAD["Load row at idx<br/>Create PinderSystem"]
    
    APPLY_FILTERS["apply_dimer_filters()"]
    CHECK_SYSTEM{"Is valid<br/>PinderSystem?"}
    RANDOM_IDX["idx = random.choice(range(len(self)))"]
    
    BUILD["select_monomer()<br/>_create_target_feature_complex()"]
    APPLY_STRUCT["apply_structure_filters()"]
    CHECK_VALID{"Pass filters?"}
    
    APPLY_TRANSFORMS["Apply structure_transforms_*"]
    RETURN["Return (system, feature, target)"]
    ERROR["Raise IndexError"]
    
    START --> INIT
    INIT --> LOOP_START
    
    LOOP_START -->|"Continue"| INCREMENT
    LOOP_START -->|"Stop - valid"| APPLY_TRANSFORMS
    LOOP_START -->|"Stop - exceeded"| ERROR
    
    INCREMENT --> LOAD
    LOAD --> APPLY_FILTERS
    APPLY_FILTERS --> CHECK_SYSTEM
    
    CHECK_SYSTEM -->|"No"| RANDOM_IDX
    CHECK_SYSTEM -->|"Yes"| BUILD
    
    BUILD --> APPLY_STRUCT
    APPLY_STRUCT --> CHECK_VALID
    
    CHECK_VALID -->|"No"| RANDOM_IDX
    CHECK_VALID -->|"Yes"| LOOP_START
    
    RANDOM_IDX --> LOOP_START
    APPLY_TRANSFORMS --> RETURN
```

**Configurable Retry:** The `max_load_attempts` parameter (default: 10) controls how many times the loader will retry before raising an `IndexError`.

**Sources:** [src/pinder-core/pinder/core/loader/loader.py:357-405]()

## Performance Considerations

### Caching Strategies

| Strategy | Dataset | Benefit | Cost |
|----------|---------|---------|------|
| **Index Filtering** | Both | Pre-filter index to reduce iteration | Higher memory for large filters |
| **Processed Files** | PPIDataset | Avoid repeated graph construction | Disk space (~1-2 MB per system) |
| **Multiprocessing** | PPIDataset | Parallel preprocessing | CPU and memory overhead |
| **DataLoader Workers** | Both | Parallel batch preparation | Memory per worker |

### Memory Optimization

**PinderDataset:**
- Loads structures on-demand (no caching)
- Memory scales with batch size and number of workers
- Use `num_workers=0` for low-memory environments

**PPIDataset:**
- Preprocesses all systems to disk once
- Memory usage during iteration is low
- Set `use_cache=True` to skip reprocessing

### Parallel Processing

`PPIDataset` supports parallel preprocessing via multiprocessing:

```python
dataset = PPIDataset(
    ...,
    parallel=True,
    max_workers=8,  # Number of parallel workers
    force_reload=False,  # Skip existing processed files
)
```

**Sources:** [src/pinder-core/pinder/core/loader/dataset.py:511-541]()

## Integration with Training Loops

The data loading system integrates seamlessly with standard PyTorch training patterns:

```python
# Standard PyTorch training loop
model = MyModel()
optimizer = torch.optim.Adam(model.parameters())

for epoch in range(num_epochs):
    for batch in train_loader:
        target = batch["target_complex"]
        feature = batch["feature_complex"]
        
        # Forward pass
        pred = model(feature)
        loss = criterion(pred, target)
        
        # Backward pass
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
```

**Batch Structure:** Both `PinderDataset` and `PPIDataset` return batches with consistent keys, making it easy to switch between tensor-based and graph-based models with minimal code changes.

**Sources:** [src/pinder-core/pinder/core/loader/dataset.py:1-19]()

---
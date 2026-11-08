# AlphaFold Integration

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [example/output.tar.gz](example/output.tar.gz)
- [example/pae_model_4_multimer_v3_pred_4.json](example/pae_model_4_multimer_v3_pred_4.json)
- [example/test.fasta](example/test.fasta)
- [example/unrelaxed_model_4_multimer_v3_pred_4.pdb](example/unrelaxed_model_4_multimer_v3_pred_4.pdb)
- [run_disobind.py](run_disobind.py)

</details>



## Purpose and Scope

This page explains how to combine Disobind predictions with AlphaFold2 (AF2) and AlphaFold3 (AF3) structural predictions to achieve enhanced accuracy. The integration leverages structural confidence metrics (pLDDT, PAE) to filter high-confidence contacts and combines them with Disobind's sequence-based predictions.

For basic prediction usage without AlphaFold integration, see [Running Predictions](run_disobind.py). For information about processing AlphaFold predictions during analysis and evaluation, see [Processing AlphaFold Predictions](analysis/get_af_prediction.py).

---

## Overview

Disobind can be used standalone or enhanced with pre-computed AlphaFold structures. When AlphaFold predictions are available, the system:

1. Extracts contact maps from AF2/AF3 structures using distance thresholds.
2. Filters contacts using confidence metrics (pLDDT ≥ 70, PAE ≤ 5).
3. Applies coarse-graining to match Disobind resolution (1, 5, or 10 residues per bead).
4. Combines predictions using a max operation: `max(Disobind, AF)`.

This combination approach captures both sequence-based patterns (Disobind) and structural information (AlphaFold), typically improving prediction accuracy on well-folded complexes.

Sources: [run_disobind.py:1-14](), [run_disobind.py:69-76](), [run_disobind.py:620-661]()

---

## Integration Workflow

```mermaid
flowchart TD
    Input["CSV Input with AF Files"] --> CheckAF{"AF Files<br/>Provided?"}
    
    CheckAF -->|No| DisoOnly["Disobind Only<br/>Prediction"]
    CheckAF -->|Yes| ParseAF["Parse AF Structure<br/>AfPrediction class"]
    
    ParseAF --> GetParser["get_parser()<br/>PDBParser/MMCIFParser"]
    GetParser --> LoadStruct["get_structure()<br/>BioPython Structure"]
    
    LoadStruct --> ExtractData["Extract Data"]
    ExtractData --> GetCoords["get_coordinates()<br/>CA coordinates"]
    ExtractData --> GetPAE["get_pae()<br/>PAE matrix"]
    ExtractData --> ExtractPLDDT["extract_perresidue_quantity()<br/>pLDDT values"]
    
    GetCoords --> CreateCMap["create_contact_map()<br/>Distance < 8Å"]
    GetPAE --> FilterConf["Filter by Confidence"]
    ExtractPLDDT --> FilterConf
    
    CreateCMap --> FilterConf
    FilterConf --> ConfCMap["Confident Contact Map<br/>pLDDT ≥ 70, PAE ≤ 5"]
    
    ConfCMap --> ProcessAF["process_af2_pred()<br/>Apply CG & Task Type"]
    
    DisoOnly --> DisoOutput["Disobind Output"]
    ProcessAF --> Combine["Combine Predictions<br/>max(Disobind, AF)"]
    Combine --> Combined["Combined Output<br/>Disobind+AF"]
    
    DisoOutput --> SaveAll["Save All Predictions"]
    Combined --> SaveAll
    
    SaveAll --> OutputCSV["CSV Files"]
    SaveAll --> OutputNPY["NPY Dictionary"]
    
    style CheckAF fill:#f9f9f9
    style FilterConf fill:#f9f9f9
    style Combine fill:#f9f9f9
```

**Diagram: AlphaFold Integration Workflow**

This diagram shows the complete pipeline for integrating AlphaFold predictions within `run_disobind.py`. The key decision point is the presence of AF2/AF3 file paths in the input CSV [run_disobind.py:10]().

Sources: [run_disobind.py:620-661](), [run_disobind.py:831-1053]()

---

## Input Format with AlphaFold Files

To use AlphaFold integration, provide additional columns in the input CSV as defined in the `read_csv_input` logic [run_disobind.py:228-255]():

```csv
Uni_ID1,start_res1,end_res1,Uni_ID2,start_res2,end_res2,af2_struct_file,af2_json_file,chain1,chain2,offset1,offset2
```

### Input Columns

| Column | Description |
|--------|-------------|
| `af2_struct_file` | Path to AF2/AF3 PDB or CIF structure file [run_disobind.py:10]() |
| `af2_json_file` | Path to JSON file containing PAE matrix [run_disobind.py:10]() |
| `chain1` | Chain ID for protein 1 in the structure [run_disobind.py:10]() |
| `chain2` | Chain ID for protein 2 in the structure [run_disobind.py:10]() |
| `offset1` | Residue numbering offset for protein 1 [run_disobind.py:10]() |
| `offset2` | Residue numbering offset for protein 2 [run_disobind.py:10]() |

### Offsets

Offsets account for differences between structure residue numbering and UniProt positions. The `add_offset` method [run_disobind.py:958-963]() handles this mapping:
- **Corrected position** = Structure position + Offset.
- This is critical when using full-length AF predictions for fragment-based Disobind tasks.

Sources: [run_disobind.py:10](), [run_disobind.py:232-253](), [run_disobind.py:958-963]()

---

## The AfPrediction Class

The `AfPrediction` class in `run_disobind.py` handles all AlphaFold structure processing:

```mermaid
classDiagram
    class AfPrediction {
        +str struct_file_path
        +str data_file_path
        +List chains
        +List offsets
        +float dist_threshold = 8
        +int plddt_cutoff = 70
        +int pae_cutoff = 5
        +Structure structure
        +ndarray pae
        
        +get_parser() parser
        +get_structure(parser) Structure
        +get_chains() Generator
        +get_residues(chain) Generator
        +extract_perresidue_quantity(residue, quantity) value
        +get_data_dict() dict
        +get_pae() void
        +add_offset(curr_res, offset) int
        +is_fragment_residue(curr_res, prot_res) bool
        +get_required_residues(prot1_res, prot2_res) Generator
        +get_confident_interactions(prot1_res, prot2_res) ndarray
    }
```

**Diagram: AfPrediction Class Structure**

The `AfPrediction` class [run_disobind.py:831-1053]() is instantiated per protein pair. It encapsulates all logic for parsing structures, extracting coordinates and confidence metrics, and generating confident contact maps.

Sources: [run_disobind.py:831-1053]()

### Key Methods

| Method | Purpose | Source |
|--------|---------|--------|
| `get_parser()` | Selects `PDBParser` or `MMCIFParser` based on file extension | [run_disobind.py:856-869]() |
| `get_structure(parser)` | Loads structure using BioPython | [run_disobind.py:871-876]() |
| `get_pae()` | Loads and symmetrizes PAE matrix: (PAE + PAE.T)/2 | [run_disobind.py:942-956]() |
| `get_required_residues()` | Filters residues belonging to specified fragments using offsets | [run_disobind.py:986-1005]() |
| `get_confident_interactions()` | Main method: generates filtered contact map using pLDDT and PAE | [run_disobind.py:1007-1053]() |

Sources: [run_disobind.py:856-1053]()

---

## Extracting Structure Data

### Coordinate Extraction

The system extracts Cα (CA) atom coordinates and pLDDT values (stored in the B-factor field) for each residue:

```python
# Logic from AfPrediction.extract_perresidue_quantity
if quantity == "coords":
    return residue["CA"].coord
elif quantity == "plddt":
    return residue["CA"].bfactor
```

The `get_required_residues()` generator [run_disobind.py:986-1005]() ensures only residues within the user-defined `[start, end]` range are processed.

Sources: [run_disobind.py:909-923](), [run_disobind.py:986-1005]()

### Contact Map Generation

Contact maps are generated from inter-residue distances within the `get_confident_interactions` method [run_disobind.py:1032-1048]():

```python
# Calculate pairwise Euclidean distances
dist = np.linalg.norm(coords1[:, np.newaxis] - coords2, axis=2)
# Apply threshold (default 8Å)
contact_map = (dist <= self.dist_threshold).astype(int)
```

Sources: [run_disobind.py:1032-1048](), [run_disobind.py:72]()

---

## Confidence Filtering

### Filtering Logic

Disobind applies a strict confidence mask to AlphaFold predictions to avoid false positives from low-confidence structural regions:

```mermaid
flowchart TD
    RawCMap["Raw Contact Map<br/>from Structure"] --> CreateMasks["Create Confidence Masks"]
    
    CreateMasks --> PLDDTMask["pLDDT Mask<br/>pLDDT₁ ≥ 70 AND pLDDT₂ ≥ 70"]
    CreateMasks --> PAEMask["PAE Mask<br/>PAE ≤ 5"]
    
    PLDDTMask --> CombineMasks["Combine Masks<br/>CMap * pLDDT_mask * PAE_mask"]
    PAEMask --> CombineMasks
    
    CombineMasks --> ConfidentCMap["Confident Contact Map<br/>Only high-confidence contacts"]
```

**Diagram: Confidence Filtering Process**

The implementation in `get_confident_interactions` [run_disobind.py:1007-1028]() creates a pLDDT mask by ensuring both residues in a pair meet the threshold, and a PAE mask for the specific pair error.

Sources: [run_disobind.py:1007-1028](), [analysis/get_af_prediction.py:315-361]()

---

## Combining Disobind and AlphaFold Predictions

### Max Operation

The combination strategy uses an element-wise maximum in `process_af2_pred` [run_disobind.py:629-633]():

```python
# Combine Disobind output and AF2 prediction
diso_af2 = np.stack([output.reshape(-1), af2_pred.reshape(-1)], axis=1)
diso_af2 = np.max(diso_af2, axis=1).reshape(m, n)
```

This approach overlays AF's confident contacts (binary 0/1) onto Disobind's probabilities (0-1).

Sources: [run_disobind.py:629-633]()

---

## Task-Specific Processing

### Interface Conversion

If the task is `interface`, the contact map is converted to binary interface vectors [run_disobind.py:682-691]():

```python
# Convert 2D cmap to 1D interface vector
idx = np.where(af2_pred == 1)
p1 = np.zeros((af2_pred.shape[0], 1))
p1[idx[0]] = 1
p2 = np.zeros((af2_pred.shape[1], 1))
p2[idx[1]] = 1
af2_pred = np.concatenate([p1, p2], axis=0)
```

### Coarse-Graining

AlphaFold predictions are coarse-grained to match Disobind's resolution using `nn.MaxPool2d` [run_disobind.py:676-680]():

```python
if self.objective[1] > 1:
    m = nn.MaxPool2d(kernel_size=self.objective[1], stride=self.objective[1])
    af2_pred = m(af2_pred.unsqueeze(0).unsqueeze(0))
```

Sources: [run_disobind.py:676-691]()

---

## Batch Processing for Evaluation

For large-scale evaluation, the `AF2MPredictions` class in `analysis/get_af_prediction.py` processes multiple structures:

```mermaid
classDiagram
    class AF2MPredictions {
        +int version
        +str base_path
        +float dist_threshold = 8
        +int plddt_threshold = 70
        +int pae_threshold = 5
        +float iptm_cutoff = 0.0
        +forward() void
        +get_best_model(model_path, header) tuple
        +get_PAE_matrix(data_file_path) ndarray
        +get_af2m_prediction() void
        +get_af3_prediction() void
    }
```

**Diagram: AF2MPredictions Class for Batch Processing**

This class [analysis/get_af_prediction.py:28-497]() automates model selection (choosing the structure with the highest ipTM) and handles the differing data formats between AF2 (PDB + JSON) and AF3 (mmCIF + JSON).

Sources: [analysis/get_af_prediction.py:28-497](), [analysis/get_af_prediction.py:79-93]()

---
# Step 1: Case Labeling & Subset Filtering

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [AlphaFlex/Data_Inputs/AF2_9606_HUMAN_v4_num_residues.json](AlphaFlex/Data_Inputs/AF2_9606_HUMAN_v4_num_residues.json)
- [AlphaFlex/Data_Inputs/AlphaFlex_database_Nov2025.json](AlphaFlex/Data_Inputs/AlphaFlex_database_Nov2025.json)
- [AlphaFlex/Data_Inputs/Test_Structures/O14653.pdb](AlphaFlex/Data_Inputs/Test_Structures/O14653.pdb)
- [AlphaFlex/Step_1B_subset_label.py](AlphaFlex/Step_1B_subset_label.py)
- [AlphaFlex/Step_1_case_label.py](AlphaFlex/Step_1_case_label.py)
- [AlphaFlex/utils/file_ops.py](AlphaFlex/utils/file_ops.py)
- [AlphaFlex/utils/smart_scoring.py](AlphaFlex/utils/smart_scoring.py)

</details>



This stage of the AlphaFlex pipeline is responsible for the topological classification of proteins based on their Intrinsically Disordered Regions (IDRs) and the subsequent filtering of these proteins into specific subsets for ensemble generation. It processes the raw AlphaFold2-derived database to assign hierarchical categories (IDP, Tails, Linkers, Loops) and enables users to define custom datasets based on length and IDR counts.

## IDR Classification Logic

The script `Step_1_case_label.py` implements a hierarchical classification system. It evaluates the position of IDRs relative to folded domains and uses Predicted Aligned Error (PAE) data to distinguish between internal disordered segments.

### Hierarchical Categories
Proteins are assigned exactly one category based on the most complex IDR type present:

| Category | Name | Description | Implementation |
| :--- | :--- | :--- | :--- |
| **0** | **IDP** | Fully disordered protein; no folded domains detected. | `is_full_idp` [AlphaFlex/Step_1_case_label.py:111-123]() |
| **1** | **Tails** | Folded domain(s) with only N-terminal or C-terminal IDRs. | `has_tail` [AlphaFlex/Step_1_case_label.py:191-193]() |
| **2** | **Linkers** | Internal IDRs between non-interacting folded domains. | `has_linker` [AlphaFlex/Step_1_case_label.py:188-190]() |
| **3** | **Loops** | Internal IDRs between interacting domains (or within one). | `has_loop` [AlphaFlex/Step_1_case_label.py:185-187]() |

### Implementation Flow
The classifier iterates through the `master_db` and uses residue ranges from the `num_residues_db` to determine terminal status. For internal IDRs, it checks the `interaction_set` (derived from `mean_pae` thresholds in the input JSON) to decide between "Loop" and "Linker" labels.

### Data Entity Mapping: Classification
The following diagram bridges the logical classification categories to the specific code variables and data structures used in `Step_1_case_label.py`.

**Classification Logic Flow**
```mermaid
graph TD
    DB["AlphaFlex_database_Nov2025.json"] -- "prot_id, data" --> Loop["Main Classification Loop"]
    Loop -- "check range [1, num_res]" --> C0["Category 0: is_full_idp"]
    Loop -- "check is_n_term | is_c_term" --> C1["Category 1: has_tail"]
    Loop -- "check interaction_set" --> C3["Category 3: has_loop"]
    Loop -- "else" --> C2["Category 2: has_linker"]
    
    subgraph "Code Entities"
    C0 --> L1["labeled_idrs.append(type='IDP')"]
    C1 --> L2["labeled_idrs.append(type='Tail IDR')"]
    C2 --> L3["labeled_idrs.append(type='Linker IDR')"]
    C3 --> L4["labeled_idrs.append(type='Loop IDR')"]
    end
```
Sources: [AlphaFlex/Step_1_case_label.py:83-200](), [AlphaFlex/Data_Inputs/AlphaFlex_database_Nov2025.json:1-55]()

## Subset Filtering & Generation

`Step_1B_subset_label.py` acts as a query engine for the labeled database. It allows researchers to extract specific UniProt IDs that meet structural criteria for downstream simulation.

### Filter Criteria
1.  **Protein Length**: Filters based on total residue count using `min_len` and `max_len` [AlphaFlex/Step_1B_subset_label.py:71-77]().
2.  **IDR Count**: Enforces specific numbers of Tails, Linkers, or Loops. Can be `exact` (e.g., exactly 2 tails) or `minimum` (at least 2 tails) [AlphaFlex/Step_1B_subset_label.py:32-48]().
3.  **IDR Length**: Filters individual IDR segments by length (e.g., only proteins where all IDRs are between 20-50 residues) [AlphaFlex/Step_1B_subset_label.py:126-148]().
4.  **Sampling**: Randomly selects a maximum number of IDs via `max_samples` [AlphaFlex/Step_1B_subset_label.py:170-171]().

### Data Entity Mapping: Filtering
This diagram shows how CLI arguments in `Step_1B_subset_label.py` map to the internal processing of the labeled JSON database.

**Subset Filtering Entity Map**
```mermaid
graph LR
    subgraph "CLI Arguments (args)"
    A1["--tail_count"]
    A2["--idr_min_len"]
    A3["--exact"]
    end

    subgraph "Internal Functions"
    F1["count_idr_types()"]
    F2["check_type_counts()"]
    F3["idr_details loop"]
    end

    subgraph "Outputs"
    O1["custom_subsets/.txt"]
    O2["advanced_info/_report.txt"]
    end

    A1 --> F1
    A3 --> F2
    F1 --> F2
    A2 --> F3
    F2 --> O1
    F3 --> O2
```
Sources: [AlphaFlex/Step_1B_subset_label.py:22-63](), [AlphaFlex/Step_1B_subset_label.py:107-167]()

## Data Schemas

### AlphaFlex_database_Nov2025.json
The input database contains structural metadata for each UniProt entry.
*   `idrs`: List of `[start, end]` residue ranges [AlphaFlex/Data_Inputs/AlphaFlex_database_Nov2025.json:3-8]().
*   `mean_pae`: A dictionary mapping domain pairs (e.g., "F1-F2") to their average PAE value. Low PAE indicates interaction [AlphaFlex/Data_Inputs/AlphaFlex_database_Nov2025.json:41-48]().
*   `interactions`: A list of domain pairs (e.g., `["F1", "F2"]`) pre-determined to be interacting [AlphaFlex/Data_Inputs/AlphaFlex_database_Nov2025.json:49-54]().

### Output ID List Format
The primary output of Step 1B is a plain text file containing one UniProt ID per line, located in `output_root/custom_subsets/`. This file is used as the `--split_index` or input list for Step 2 and Step 3 [AlphaFlex/Step_1B_subset_label.py:178-192]().

## Summary Reports
For every subset generated, the system produces a filtration report (`_report.txt`) including:
*   Applied length constraints.
*   Breakdown of IDR types found.
*   Detailed list of every IDR range and length for the included proteins.
*   A statistical summary of categories (0-3) represented in the subset.

Sources: [AlphaFlex/Step_1B_subset_label.py:194-210](), [AlphaFlex/Step_1_case_label.py:215-220]()

---
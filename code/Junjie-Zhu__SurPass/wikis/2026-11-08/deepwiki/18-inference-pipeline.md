# Inference Pipeline

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [configs/inference.yaml](configs/inference.yaml)
- [src/inference.py](src/inference.py)

</details>



The inference pipeline executes trained `ResOnly` models on protein-protein interaction (PPI) test sets [src/inference.py:17](). It ingests structural protein sequences from FASTA files [src/inference.py:78](), pre-computed protein Language Model (pLM) embeddings [src/inference.py:114](), and test metadata TSVs [src/inference.py:167](). It coordinates data batching [src/inference.py:151](), long-sequence windowing via `plan_ppi_windows` [src/inference.py:30](), multi-head prediction scoring across protein, residue, and pair levels [src/inference.py:33](), and output generation.

Sources: [src/inference.py:1-47](), [src/inference.py:78-155](), [src/inference.py:157-192]()

---

## Configuration System (`configs/inference.yaml`)

The inference configuration governs model architecture hyperparameters, input data directories, threshold settings, and execution flags. It maps directly to model construction and dataset instantiation routines.

| Configuration Block | Key Parameters | Code Mapping |
| :--- | :--- | :--- |
| `data` | `test_set_tsv`, `fasta_path`, `embedding_dir`, `max_full_length`, `contact_threshold` | `InferencePairDataset` [src/inference.py:157-192]() |
| `model` | `n_layers`, `dim_token`, `dim_pair`, `n_heads`, `plm_in_dim`, `plm_out_dim` | `ResOnly` [src/inference.py:17]() |
| `metrics` | `positive_weight` | Evaluation weighting [src/inference.py:22-24]() |
| `performance` | `checkpoint_pair_blocks` | Gradient checkpointing flags [src/inference.py:46-48]() |

```yaml
# configs/inference.yaml excerpt
task_prefix: infer
seed: 42
ckpt_dir: null
data:
  test_set_tsv: ${data.root_dir}/test_set.tsv
  fasta_path: ${data.root_dir}/test_set.fasta
  embedding_dir: ${data.root_dir}/embeddings
model:
  n_layers: 8
  dim_token: 768
  dim_pair: 512
```

Sources: [configs/inference.yaml:1-52]()

---

## Data Loading and Preprocessing (`InferencePairDataset`)

The inference dataset infrastructure parses raw FASTA strings, loads tensorized pLM features, encodes residue types into vocabulary indices, and constructs batch representations suitable for dual-protein complex evaluation.

```mermaid
graph TD
    A["test_set.tsv"] --> B["InferencePairDataset"]
    C["test_set.fasta"] --> D["load_fasta_sequences"]
    D --> B
    E["embedding_dir"] --> F["load_plm_embedding"]
    F --> B
    B --> G["build_residue_features"]
    G --> H["encode_residue_types"]
    H --> I["inference_collate_fn"]
    I --> J["collate_concatenated_pairs"]

    sub-graph Code_Entities
        D
        F
        G
        H
        I
        J
        B
    end

    style B fill:#fff,stroke:#000,stroke-width:2px
    style D fill:#fff,stroke:#000,stroke-width:2px
    style F fill:#fff,stroke:#000,stroke-width:2px
    style G fill:#fff,stroke:#000,stroke-width:2px
    style H fill:#fff,stroke:#000,stroke-width:2px
    style I fill:#fff,stroke:#000,stroke-width:2px
    style J fill:#fff,stroke:#000,stroke-width:2px
```
*Figure 5.1.1: Data loading flow from raw disk inputs to collated token tensors.*

Sources: [src/inference.py:78-192]()

### FASTA and Embedding Parsers
- **`load_fasta_sequences`**: Reads FASTA files sequentially, validates identifiers and non-empty sequences, and outputs a dictionary mapping protein IDs to uppercase sequence strings [src/inference.py:78-111]().
- **`load_plm_embedding`**: Loads pre-computed embedding files supporting raw `torch.Tensor` objects or nested dictionary wrappers containing the `plm_emb` key, coercing types to `torch.float32` [src/inference.py:114-120]().
- **`encode_residue_types`**: Maps residue characters to integer vocabularies using `restype_order` with fallback to `unk_restype_index` [src/inference.py:12-14](), [src/inference.py:123-127]().
- **`build_residue_features`**: Combines token embeddings, residue codes, positions, and masks into a standard feature dictionary [src/inference.py:130-148]().

Sources: [src/inference.py:78-148]()

### Dataset and Collation
`InferencePairDataset` iterates through the test set TSV metadata, resolves protein pairs via `parse_protein_pair`, checks sequence presence, and caches valid query items [src/inference.py:157-192](). The helper `inference_collate_fn` unpacks batch tuples and delegates pairing concatenation to `collate_concatenated_pairs` [src/inference.py:151-155]().

Sources: [src/inference.py:151-192]()

---

## Inference Execution & Prediction Heads

During execution, complexes exceeding length thresholds are processed via long-sequence windowing (`plan_ppi_windows`) [src/inference.py:30](). The `ResOnly` model computes predictions across defined PPI heads (`PPI_HEADS = ("protein", "residue", "pair")`) [src/inference.py:17](), [src/inference.py:33]().

```mermaid
graph TD
    A["InferencePairDataset"] --> B["inference_collate_fn"]
    B --> C["plan_ppi_windows"]
    C --> D["ResOnly"]
    D --> E["pair_bind_scores"]
    D --> F["_binary_classification_metrics"]
    E --> G["PREDICTION_COLUMNS"]
    F --> G

    sub-graph Code_Execution
        B
        C
        D
        E
        F
    end

    style A fill:#fff,stroke:#000,stroke-width:2px
    style B fill:#fff,stroke:#000,stroke-width:2px
    style C fill:#fff,stroke:#000,stroke-width:2px
    style D fill:#fff,stroke:#000,stroke-width:2px
    style E fill:#fff,stroke:#000,stroke-width:2px
    style F fill:#fff,stroke:#000,stroke-width:2px
    style G fill:#fff,stroke:#000,stroke-width:2px
```
*Figure 5.1.2: Execution pipeline mapping data batches to model evaluation heads and metrics.*

Sources: [src/inference.py:17-33](), [src/inference.py:151-155]()

Sources: [src/inference.py:17-33](), [src/inference.py:151-155]()

---

## Metrics and Output Generation

Results are compiled into structured dataframes containing prediction probabilities, binary classification decisions, contact counts, and forward pass statistics. Output columns are explicitly defined by `PREDICTION_COLUMNS` [src/inference.py:34-46](). Serialized outputs are written to disk as TSV, CSV, or PyTorch tensor (`.pt`) files.

```python
PREDICTION_COLUMNS = [
    "Protein pairs",
    "Category",
    "p_protein",
    "pred_protein",
    "p_residue",
    "pred_residue",
    "n_pos_residues",
    "p_pair",
    "pred_pair",
    "n_inter_contacts",
    "n_forwards",
]
```

Sources: [src/inference.py:34-46]()

---
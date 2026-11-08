# Model Architecture

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/model/__init__.py](src/model/__init__.py)
- [src/model/surpass.py](src/model/surpass.py)

</details>



## Purpose and Scope

This page provides a high-level overview of the `ResOnly` network architecture [src/model/surpass.py:28-46] implemented in SurPass. The model processes protein sequences and protein Language Model (pLM) embeddings through dual concurrent streams: a 1D residue stream (`dim_token`) and a 2D pair stream (`dim_pair`). By iteratively exchanging information between token representations and pairwise distances/contacts via attention pair biases, outer product means, and triangle updates, SurPass predicts inter-chain residue contacts, residue-level binding scores, and overall protein-protein interaction (PPI) logits.

This is a parent overview page. Detailed technical implementations, tensor manipulations, and layer-by-layer mathematics are documented in the child pages linked below:
- **[ResOnly Model and SurpassOutput](#3.1)** — Forward pass orchestration, layer stacking, activation checkpointing, and prediction output heads.
- **[Embedding Layer](#3.2)** — `ResidueEmbedder`: pLM feature extraction, residue type encodings, relative position encoding (`r_max`), and initial pair construction.
- **[Transformer and Attention Primitives](#3.3)** — `ResidueTransformer`, multi-head attention with pair bias, QK LayerNorm, and attention primitives.
- **[Pair Stack: Outer Product Mean and Triangle Updates](#3.4)** — `EfficientOuterProductMean`, `TriangleMultiplicationOutgoing`, `TriangleMultiplicationIncoming`, and `PairTransition`.

Sources: [src/model/surpass.py:20-46]()

---

## High-Level Architecture Overview

The SurPass architecture bridges sequence-level embeddings with 2D structural relationship maps. The network is structured around a stack of $N$ alternating blocks (`n_layers`) [src/model/surpass.py:31-88] where the residue stream and pair stream co-evolve.

```mermaid
graph TD
    sub5 ["Input Processing"]
    A["ResidueInput / pLM"] --> B["ResidueEmbedder"]
    
    subgraph Stream ["Dual-Stream Architecture"]
        B --> C["SingleRepresentation (dim_token)"]
        B --> D["PairRepresentation (dim_pair)"]
        
        C --> E["ResidueTransformer"]
        D --> F["EfficientOuterProductMean"]
        D --> G["TriangleMultiplicationOutgoing/Incoming"]
        D --> H["PairTransition"]
        
        E <-->|"Pair Bias & OPM"| F
        F --> G --> H
    end
    
    subgraph Outputs ["SurpassOutput"]
        E --> I["residue_bind_linear"]
        H --> J["pair_out_linear"]
        E --> K["protein_head"]
    end

    classDef default fill:#fff,stroke:#000,stroke-width:1px;
```
*Figure 1: High-level data flow through the `ResOnly` dual-stream architecture from `ResidueEmbedder` to `SurpassOutput` [src/model/surpass.py:143-192].*

Sources: [src/model/surpass.py:143-192]()

---

## 3.1 ResOnly Model and SurpassOutput

The core orchestrator of the architecture is the `ResOnly` module [src/model/surpass.py:28-46], which manages the iterative update loops across `n_layers`. During the forward pass [src/model/surpass.py:143-192], input dictionaries are passed to the embedder, generating initial single and pair representations. Within each layer, `ResidueTransformer` updates the token stream, while `EfficientOuterProductMean`, triangle multiplication blocks, and `PairTransition` update the pair stream.

The model encapsulates its predictions in a structured `SurpassOutput` dataclass [src/model/surpass.py:20-27] containing `pair_logits`, `residue_logits`, `protein_logit`, and `pair_mask`. 

For deep technical details on the forward pass, layer sequencing, activation checkpointing (`checkpoint_pair_blocks`), and output head formulations, see **[ResOnly Model and SurpassOutput](#3.1)**.

Sources: [src/model/surpass.py:20-192]()

---

## 3.2 Embedding Layer

The embedding subsystem (`ResidueEmbedder`) serves as the bridge between raw inputs (such as precomputed protein Language Model representations and sequence identities) and the internal vector spaces of `ResOnly`. It projects input features into a token dimension (`dim_token`) and initializes a 2D pairwise representation (`dim_pair`) incorporating relative position encodings governed by maximum distance parameters (`r_max`).

For detailed information on embedding combinations, feature dimensions, and relative positional encoding schemes, see **[Embedding Layer](#3.2)**.

Sources: [src/model/surpass.py:68-72]()

---

## 3.3 Transformer and Attention Primitives

At the heart of the 1D feature processing stream is `ResidueTransformer` [src/model/surpass.py:76-88]. It implements multi-head attention mechanisms enhanced by 2D pair biases, allowing pairwise structural context to directly modulate residue-to-residue attention weights. It incorporates QK LayerNorm (`use_qkln`) for training stability and residue-level transition expansions.

For comprehensive architectural specifics regarding attention mechanisms, pair bias integration, and stable layer normalization, see **[Transformer and Attention Primitives](#3.3)**.

Sources: [src/model/surpass.py:75-88]()

---

## 3.4 Pair Stack: Outer Product Mean and Triangle Updates

To allow 1D representations to inform 2D pairwise features, SurPass uses `EfficientOuterProductMean` [src/model/surpass.py:90-96]. Additionally, geometric constraints are propagated across the pair stream using `TriangleMultiplicationOutgoing` [src/model/surpass.py:98-103], `TriangleMultiplicationIncoming` [src/model/surpass.py:105-110], and `PairTransition` [src/model/surpass.py:112-118]. These modules ensure triangle inequality constraints and relational consistency across interacting residues.

For detailed equations, hidden dimension configurations (`dim_opm_inner`, `dim_triangle_hidden`), and tensor contraction details, see **[Pair Stack: Outer Product Mean and Triangle Updates](#3.4)**.

Sources: [src/model/surpass.py:90-118]()

---
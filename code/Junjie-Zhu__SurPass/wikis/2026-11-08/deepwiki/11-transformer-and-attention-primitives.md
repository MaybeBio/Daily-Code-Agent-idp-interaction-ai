# Transformer and Attention Primitives

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/model/components/calvados_functions.py](src/model/components/calvados_functions.py)
- [src/model/components/calvados_residues.pickle](src/model/components/calvados_residues.pickle)
- [src/model/components/primitives.py](src/model/components/primitives.py)
- [src/model/components/transformer.py](src/model/components/transformer.py)
- [src/model/components/triangle_attention.py](src/model/components/triangle_attention.py)

</details>



## Purpose and Scope

This page details the foundational building blocks, attention mechanisms, update primitives, and physical energy functions residing in `src/model/components/`. These primitives form the structural backbone of SurPass's residue and pair streams, enabling residue-level feature updates, pair-conditioned multi-head attention, geometric triangle attention, and integration with the CALVADOS coarse-grained force field.

---

## 1. Primitive Modules and Activation Functions

SurPass implements custom neural network primitives tailored for structural biology models, featuring adaptive normalizations, gated activation functions, and transition blocks.

### Adaptive Layer Normalization (`AdaptiveLayerNorm` and `AdaptiveLayerNormOutputScale`)
`AdaptiveLayerNorm` applies a standard LayerNorm with element-wise affine parameters ($\gamma$ and $\beta$) dynamically generated from conditioning variables (e.g., global or context embeddings) via linear projections and a sigmoid activation on $\gamma$ [src/model/components/primitives.py:11-43](). `AdaptiveLayerNormOutputScale` applies zero-initialized gating to scale representations conditioned on auxiliary signals [src/model/components/primitives.py:45-73]().

### Gated Linear Units and Transitions (`SwiGLU`, `Transition`, `PairTransition`)
* **`SwiGLU`**: Splits the input tensor along the final dimension into two halves, applying the SiLU (Swish) activation function to the gate half before element-wise multiplication [src/model/components/primitives.py:75-88]().
* **`Transition` / `PairTransition`**: Feed-forward blocks that expand hidden dimensions by an `expansion_factor` (default 4), apply `SwiGLU`, and project back to the model dimension, with an optional initial LayerNorm [src/model/components/primitives.py:90-151]().

```mermaid
graph TD
    subcode["Code Entity Space"]
        IN["Input Tensor (x)"] --> LN["torch.nn.LayerNorm"]
        LN --> SWI["SwiGLU (src/model/components/primitives.py:75)"]
        SWI --> LIN["Linear Out (src/model/components/primitives.py:106)"]
        LIN --> MASK["Mask Multiplication"]
    end
    
    style IN fill:none,stroke:#000,stroke-width:1px
    style LN fill:none,stroke:#000,stroke-width:1px
    style SWI fill:none,stroke:#000,stroke-width:1px
    style LIN fill:none,stroke:#000,stroke-width:1px
    style MASK fill:none,stroke:#000,stroke-width:1px
```
*Sources: `[src/model/components/primitives.py:11-151]()`*

---

## 2. Attention Mechanisms and QK LayerNorm

Multi-head attention in SurPass incorporates query-key layer normalization (`qk_layernorm`), rotary position embeddings (`RotaryEmbedding`), and pair-wise bias terms.

### `GatedAttention`
`GatedAttention` computes standard self-attention over residue sequences with integrated rotary position embeddings and an output gating mechanism [src/model/components/transformer.py:40-138](). 

* **QKV Projection**: Projects input features `d_model` into queries, keys, and values of inner dimension `n_heads * d_head` [src/model/components/transformer.py:59]() alongside a gate projection `g_proj` [src/model/components/transformer.py:61]().
* **QK LayerNorm**: Applies `nn.LayerNorm` independently to queries and keys prior to rotary embedding application to stabilize attention logits [src/model/components/transformer.py:64-70]().
* **Rotary Embedding**: Applies `RotaryEmbedding` to the unflattened head dimensions of $Q$ and $K$ [src/model/components/transformer.py:71-80]().
* **Scaled Dot-Product Attention**: Leverages `torch.nn.functional.scaled_dot_product_attention` (SDPA) with key padding masks [src/model/components/transformer.py:114-121]().
* **Output Gating**: Multiplies attention outputs by the sigmoid of the gate projection `g` before the final linear projection [src/model/components/transformer.py:135]().

### `PairBiasAttention`
`PairBiasAttention` integrates pair-wise representations (`pair_feats`) directly into the attention score matrix as additive bias terms [src/model/components/transformer.py:141-194]().

* **Pair Bias Projection**: When `pair_dim` is provided, pair features are normalized and projected via `nn.Linear(pair_dim, heads, bias=False)` to yield per-head spatial biases [src/model/components/transformer.py:171-175]().
* **Logit Modification**: The projected pair bias is added directly to the $QK^T$ attention matrix before softmax normalization.

```mermaid
graph TD
    subcode["Attention Architecture & Code Entities"]
        InputNode["node_feats (b,n,d_s)"] --> QKV["to_qkv (Transformer.py:165)"]
        InputPair["pair_feats (b,n,n,d_e)"] --> PairNorm["pair_norm (Transformer.py:173)"]
        PairNorm --> ToBias["to_bias (Transformer.py:172)"]
        QKV --> QKLN["q_layer_norm / k_layer_norm (Transformer.py:169)"]
        QKLN --> QKMatMul["Scaled Dot Product + Pair Bias"]
        ToBias --> QKMatMul
        QKMatMul --> Gate["to_g Output Gate (Transformer.py:166)"]
        Gate --> OutNode["to_out_node (Transformer.py:167)"]
    end

    style InputNode fill:none,stroke:#000,stroke-width:1px
    style InputPair fill:none,stroke:#000,stroke-width:1px
    style QKV fill:none,stroke:#000,stroke-width:1px
    style PairNorm fill:none,stroke:#000,stroke-width:1px
    style ToBias fill:none,stroke:#000,stroke-width:1px
    style QKLN fill:none,stroke:#000,stroke-width:1px
    style QKMatMul fill:none,stroke:#000,stroke-width:1px
    style Gate fill:none,stroke:#000,stroke-width:1px
    style OutNode fill:none,stroke:#000,stroke-width:1px
```
*Sources: `[src/model/components/transformer.py:40-194]()`*

---

## 3. Triangle Attention Primitives

To update pair representations along geometric triangle inequalities (AlphaFold2/OpenFold style), SurPass implements specialized triangle attention modules.

### `BiasGatedAttention`
`BiasGatedAttention` executes multi-head attention on pair tensors `[*, I, J, C]` using multiple additive bias tensors (e.g., masking and triangle bias) and an output gating projection [src/model/components/triangle_attention.py:12-56]().

### `TriangleAttention` (`TriangleAttentionStartingNode` & `TriangleAttentionEndingNode`)
`TriangleAttention` operates on pair representations `[*, I, J, C_in]` by treating either starting nodes or ending nodes as the sequence dimension [src/model/components/triangle_attention.py:58-164]().

* **Orientation**: 
  * `TriangleAttentionStartingNode` sets `starting=True` (Algorithm 13 equivalent) [src/model/components/triangle_attention.py:151-156]().
  * `TriangleAttentionEndingNode` sets `starting=False` (Algorithm 14 equivalent), transposing dimensions `-2` and `-3` prior to processing [src/model/components/triangle_attention.py:159-164]().
* **Triangle Bias**: Computes a linear projection of LayerNorm-normalized pair features, permutes dimensions, and adds it as an attention bias alongside sequence masks [src/model/components/triangle_attention.py:138-139]().
* **Chunking**: Supports optional sequence chunking (`_chunk`) to reduce memory consumption during large-scale pair updates [src/model/components/triangle_attention.py:96-115]().

*Sources: `[src/model/components/triangle_attention.py:1-164]()`*

---

## 4. CALVADOS Force Field and Energy Functions

SurPass includes support for the CALVADOS coarse-grained polymer physics model (`calvados_functions.py`), enabling physical energy calculations based on residue-specific interaction parameters [src/model/components/calvados_functions.py:1-121]().

### `CalvadosModel` Dataclass
Encapsulates physical interaction maps and parameters [src/model/components/calvados_functions.py:9-19]():
* `sigmamap`: DataFrame of pairwise Lennard-Jones distance parameters $\sigma_{ij}$ [src/model/components/calvados_functions.py:11]().
* `lambdamap`: DataFrame of hydrophobicity scaling parameters $\lambda_{ij}$ [src/model/components/calvados_functions.py:12]().
* `yukawa_eps`: DataFrame of electrostatic coupling amplitudes [src/model/components/calvados_functions.py:13]().
* `yukawa_kappa`: Debye-Huckel screening parameter derived from salt concentration [src/model/components/calvados_functions.py:14]().

### Key Functions
* `load_calvados_model`: Loads residue pickle tables (`calvados_residues.pickle`), adjusts histidine charge states based on pH, computes Bjerrum length, and constructs pairwise parameter matrices [src/model/components/calvados_functions.py:34-80]().
* `compute_ashbaugh_hatch`: Calculates the Ashbaugh-Hatch functional form for hydrophobicity-scaled Lennard-Jones interactions with switching functions [src/model/components/calvados_functions.py:100-110]().
* `compute_yukawa`: Computes screened electrostatic interactions via Yukawa potential [src/model/components/calvados_functions.py:112-116]().
* `compute_calvados_energy`: Combines Ashbaugh-Hatch and Yukawa terms to yield total coarse-grained interaction energies [src/model/components/calvados_functions.py:118-121]().

```mermaid
graph TD
    subcode["CALVADOS Physics Integration"]
        Pickle["calvados_residues.pickle"] --> Load["load_calvados_model (calvados_functions.py:34)"]
        Load --> Model["CalvadosModel Dataclass (calvados_functions.py:9)"]
        Model --> AH["compute_ashbaugh_hatch (calvados_functions.py:100)"]
        Model --> YU["compute_yukawa (calvados_functions.py:112)"]
        AH --> Energy["compute_calvados_energy (calvados_functions.py:118)"]
        YU --> Energy
    end

    style Pickle fill:none,stroke:#000,stroke-width:1px
    style Load fill:none,stroke:#000,stroke-width:1px
    style Model fill:none,stroke:#000,stroke-width:1px
    style AH fill:none,stroke:#000,stroke-width:1px
    style YU fill:none,stroke:#000,stroke-width:1px
    style Energy fill:none,stroke:#000,stroke-width:1px
```
*Sources: `[src/model/components/calvados_functions.py:1-121](), [src/model/components/calvados_residues.pickle:1-8]()`*

---
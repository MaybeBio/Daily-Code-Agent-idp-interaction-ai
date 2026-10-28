# Model Architecture (biomoa)

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/bioemu2/biomoa/__init__.py](src/bioemu2/biomoa/__init__.py)
- [src/bioemu2/biomoa/models/base.py](src/bioemu2/biomoa/models/base.py)
- [src/bioemu2/biomoa/models/score_model.py](src/bioemu2/biomoa/models/score_model.py)

</details>



This page provides a high-level overview of the BioMoA score model stack within `bioemu2`, covering sequence embedding trunks, dual-stream diffusion modules, core mathematical operators, and runtime precision and tiling utilities. BioMoA serves as the foundational biomolecular emulator scoring network, combining sequence-level representations with atom- and token-level spatial diffusion.

As a parent page, this section outlines the architectural components and directs you to dedicated child pages for deep technical dives into individual modules.

---

## 3.1. Score Model and Diffusion Schedule

The core architecture centers around `BioMoaCModel`, defined in `src/bioemu2/biomoa/models/score_model.py`, which integrates a sequence embedding trunk with a dual-stream diffusion module. Input features are processed into unified single and pair representations, while noise conditioning is governed by the empirical EDM (Elucidating the Design Space of Diffusion-Based Generative Models) schedule.

Outputs from the score model are wrapped in `ScoreModelOutput` (`src/bioemu2/biomoa/models/base.py`), containing the unscaled score predictions, optional distograms over distance bins, and denoised coordinates. Intermediate trunk embeddings can also be serialized or loaded via `TrunkEmbeddings` to avoid redundant sequence processing during repeated sampling runs.

For deep technical details on base model abstractions, schedule parameters, and preconditioning, see [Score Model and Diffusion Schedule](#3.1).

```mermaid
graph TD
    subcode["Code Entity Space"]
        BMC["BioMoaCModel\nsrc/bioemu2/biomoa/models/score_model.py"] --> SMO["ScoreModelOutput\nsrc/bioemu2/biomoa/models/base.py"]
        BMC --> TE["TrunkEmbeddings\nsrc/bioemu2/biomoa/models/score_model.py"]
        EDM["EDMSchedule\nsrc/bioemu2/biomoa/models/score_model.py"] --> BMC
    end
    
    sub_nlp["Natural Language Space"] --> subcode
```
*Sources: `src/bioemu2/biomoa/models/base.py:1-13`, `src/bioemu2/biomoa/models/score_model.py:29-159`*

---

## 3.2. Network Modules

BioMoA utilizes a modular network topology comprising input embedders, pairwise relative position encoders, sequence pairmixers, and dual-stream diffusion transformers. The sequence trunk processes tokenized inputs, MSA features, and pocket classes, generating rich single and pair representations that condition the spatial denoising process.

The diffusion module alternates between atom-level and token-level attention layers, incorporating time embeddings and structural conditioning. Atom features are embedded directly within the dual-stream architecture, allowing seamless integration across hierarchical resolution scales.

For comprehensive details on embedders, encoders, conditioning transformers, pairmixers, and positional embeddings, see [Network Modules](#3.2).

```mermaid
graph TD
    subcode2["Code Entity Space"]
        IE["InputEmbedder\nsrc/bioemu2/biomoa/modules/embedders.py"] --> PM["PairMixerSequenceEmbedder\nsrc/bioemu2/biomoa/modules/pairmixer.py"]
        PM --> RP["FusedRelativePositionEncoder\nsrc/bioemu2/biomoa/modules/relative_position_encoder.py"]
        RP --> DSD["DualStreamDiffusionModule\nsrc/bioemu2/biomoa/modules/diffusion.py"]
        DM["DistogramModule\nsrc/bioemu2/biomoa/modules/embedders.py"] --> DSD
    end

    sub_nlp2["Natural Language Space"] --> subcode2
```
*Sources: `src/bioemu2/biomoa/models/score_model.py:10-18`*

---

## 3.3. Core Operators and Triton Kernels

The layers within BioMoA rely on specialized mathematical operators designed for biomolecular modeling. These include attention with pair biases (`attention_pair_bias`), triangular update operations (multiplication and outgoing/incoming edges), outer product means (OPM), pair-weighted averaging, and transition blocks.

Performance-critical sections—such as the outer product mean—are accelerated using custom Triton kernels to optimize memory bandwidth and reduce intermediate tensor allocations during training and inference.

For an in-depth breakdown of individual operators, dropouts, atom/token operations, and Triton kernel implementations, see [Core Operators and Triton Kernels](#3.3).

```mermaid
graph TD
    subcode3["Code Entity Space"]
        APB["attention_pair_bias"] --> OP["Core Operators & Triton Kernels"]
        TM["triangle_multiplication"] --> OP
        OPM["outer_product_mean\n(Triton Kernel)"] --> OP
        PWA["pair_weighted_averaging"] --> OP
    end

    sub_nlp3["Natural Language Space"] --> subcode3
```
*Sources: `src/bioemu2/biomoa/models/score_model.py:95-159`*

---

## 3.4. Model Utilities: Precision, Tiling, Strategy

To handle large macromolecular complexes without exceeding GPU memory limits, BioMoA incorporates execution utilities such as dense batch-to-dictionary feature converters (`convert_densebatch_to_dictionary`), precision managers (supporting mixed `bf16` execution via `matmul_precision_ctx`), and memory tiling resolvers (`resolve_tiling_cfg`).

These utilities dynamically adjust execution strategies based on system dimensions, device capabilities, and configuration parameters, issuing warnings or fallback paths when memory thresholds are approached.

For detailed documentation on feature conversion, precision contexts, tiling configurations, and execution strategies, see [Model Utilities: Precision, Tiling, Strategy](#3.4).

```mermaid
graph TD
    subcode4["Code Entity Space"]
        FC["convert_densebatch_to_dictionary\nsrc/bioemu2/biomoa/utils/feature_converter.py"] --> RT["Model Runtime Utilities"]
        RT_CFG["resolve_tiling_cfg\nsrc/bioemu2/biomoa/utils/tiling.py"] --> RT
        MPC["matmul_precision_ctx\nsrc/bioemu2/precision.py"] --> RT
    end

    sub_nlp4["Natural Language Space"] --> subcode4
```
*Sources: `src/bioemu2/biomoa/models/score_model.py:19-22`*

---
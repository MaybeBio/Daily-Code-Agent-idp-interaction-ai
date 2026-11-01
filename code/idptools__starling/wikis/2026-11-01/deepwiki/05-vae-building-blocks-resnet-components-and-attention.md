# VAE Building Blocks: ResNet Components and Attention

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [starling/data/data_wrangler.py](starling/data/data_wrangler.py)
- [starling/models/attention.py](starling/models/attention.py)
- [starling/models/blocks.py](starling/models/blocks.py)
- [starling/models/normalization.py](starling/models/normalization.py)
- [starling/models/transformer.py](starling/models/transformer.py)
- [starling/models/vae_components.py](starling/models/vae_components.py)

</details>



This page provides a technical reference for the primitive neural building blocks used within the STARLING Variational Autoencoder (VAE). These components handle spatial processing of distance maps, normalization, conditioning, and attention-based refinement.

## ResNet Components

The VAE architecture relies on specialized Residual Blocks for both the encoder and decoder. These blocks are designed to handle 2D spatial data (distance maps) and can optionally incorporate conditioning information.

### ResBlockEncBasic
The `ResBlockEncBasic` is the fundamental building block for the `ResNet_Encoder` [starling/models/vae_components.py:13-21](). It follows the standard ResNet-18/34 design, consisting of two convolutional layers with a residual connection [starling/models/blocks.py:131-166]().

*   **Structure**: Two $3 \times 3$ (default) convolutions. The first convolution can perform spatial downsampling via the `stride` parameter [starling/models/blocks.py:180-186]().
*   **Conditioning**: Supports FiLM-style conditioning via a `time_mlp` if a `timestep` dimension is provided [starling/models/blocks.py:195-199]().

### ResBlockDecBasic
The `ResBlockDecBasic` is used in the `ResNet_Decoder` [starling/models/vae_components.py:108-117](). It is designed to mirror the encoder's structure but uses `ResizeConv2d` for upsampling to avoid checkerboard artifacts [starling/models/blocks.py:246-281]().

*   **Upsampling**: If `stride > 1`, the block utilizes `ResizeConv2d` to increase spatial dimensions [starling/models/blocks.py:302-311]().
*   **Contraction**: The block typically reduces channel depth (contraction) as it increases spatial resolution [starling/models/blocks.py:247-248]().

### ResizeConv2d
A replacement for `ConvTranspose2d`, this module uses `F.interpolate` followed by a standard `Conv2d` [starling/models/blocks.py:65-104](). This approach is mathematically preferred for generating smooth distance maps without grid-like artifacts [starling/models/blocks.py:79-82]().

**ResNet Block Data Flow**

```mermaid
graph TD
    Input["Input Tensor (B, C, H, W)"] --> Conv1["Conv2d (3x3, Stride=S)"]
    
    subgraph "Conditioning (Optional)"
        Time["Timestep/Label Embedding"] --> MLP["time_mlp (SiLU + Linear)"]
        MLP --> Scale["Scale (gamma)"]
        MLP --> Shift["Shift (beta)"]
    end

    Conv1 --> FiLM["FiLM: (x * gamma) + beta"]
    FiLM --> Norm1["Normalization (BN/IN/LN/GN/RMS)"]
    Norm1 --> Act1["ReLU"]
    Act1 --> Conv2["Conv2d (3x3, Stride=1)"]
    Conv2 --> Norm2["Normalization"]
    
    Input --> Skip["Shortcut Connection (Identity or Conv)"]
    Norm2 --> Add["Add (Residual)"]
    Skip --> Add
    Add --> Act2["ReLU (Output)"]

    style FiLM stroke-dasharray: 5 5
```
**Sources:** [starling/models/blocks.py:131-244](), [starling/models/blocks.py:65-128]()

---

## Normalization Strategies

STARLING implements a flexible normalization interface supporting multiple strategies, selectable via string configuration in the VAE factory functions [starling/models/vae_components.py:26-31]().

| Type | Class/Function | Description |
| :--- | :--- | :--- |
| **BatchNorm** | `nn.BatchNorm2d` | Standard batch-wise normalization. |
| **InstanceNorm** | `nn.InstanceNorm2d` | Normalization per sample and per channel. |
| **LayerNorm** | `LayerNorm` | Custom implementation supporting `channels_first` format (B, C, H, W) by permuting dimensions [starling/models/blocks.py:9-41](). |
| **RMSNorm** | `RMSNorm` | Root Mean Square Layer Normalization, scaling by the square root of the channel dimension [starling/models/normalization.py:6-12](). |
| **GroupNorm** | `nn.GroupNorm` | Normalizes groups of channels; STARLING defaults to 32 groups [starling/models/blocks.py:190](). |
| **AdaLayerNorm** | `AdaLayerNorm` | Adaptive Layer Norm used in Transformers, where scale and shift are predicted from a conditioning vector [starling/models/transformer.py:103-137](). |

**Sources:** [starling/models/blocks.py:9-41](), [starling/models/normalization.py:6-12](), [starling/models/transformer.py:103-137]()

---

## Attention Mechanisms

The VAE and Diffusion backbones utilize several attention variants for spatial and sequence-based conditioning.

### MultiHeadAttention
A generic implementation supporting both self-attention and cross-attention [starling/models/attention.py:11-38](). It uses `torch.nn.functional.scaled_dot_product_attention` for optimized computation (FlashAttention where available) [starling/models/attention.py:74]().

### SelfAttentionConv
A specialized attention block for 2D feature maps. It treats spatial locations as tokens and computes attention across the height and width dimensions [starling/models/attention.py:162-185]().
*   **Input**: (B, C, H, W)
*   **Mechanism**: Flattens spatial dimensions to (B, H*W, C), applies Multi-Head Attention, and reshapes back to (B, C, H, W) [starling/models/attention.py:221-236]().

### CrossAttention
Used primarily in the Diffusion model's `SequenceEncoder` to condition the 2D latents on 1D amino acid sequences [starling/models/attention.py:82-120](). It performs normalization on both query and context inputs before projection [starling/models/attention.py:129-130]().

**Attention Entity Mapping**

```mermaid
graph LR
    subgraph "Code Entity Space"
        MHA["MultiHeadAttention"]
        SA["SelfAttention"]
        SAC["SelfAttentionConv"]
        CA["CrossAttention"]
    end

    subgraph "Logical Function"
        D1["Global Spatial Context"] --> SAC
        D2["Sequence Conditioning"] --> CA
        D3["Transformer Blocks"] --> MHA
    end

    SAC -- "calls" --> SA
    CA -- "uses" --> F_SDPA["F.scaled_dot_product_attention"]
    MHA -- "uses" --> F_SDPA
```
**Sources:** [starling/models/attention.py:11-160](), [starling/models/attention.py:162-240]()

---

## Conditioning and Encodings

### FiLM (Feature-wise Linear Modulation)
Integrated into `ResBlockEncBasic` and `ResBlockDecBasic`. It modulates feature maps by applying a learned affine transformation (scale and shift) derived from an external vector (e.g., timestep or class label) [starling/models/blocks.py:218-220]().

### SinusoidalPosEmb
Generates sinusoidal embeddings for timestep encoding in the diffusion process [starling/models/transformer.py:15-30](). It maps a scalar time value to a high-dimensional vector using varying frequencies of sine and cosine functions [starling/models/transformer.py:51-57]().

**Sources:** [starling/models/blocks.py:218-220](), [starling/models/transformer.py:15-58]()

---
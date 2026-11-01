# Diffusion Model (DDPM)

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [starling/configs/configs.yaml](starling/configs/configs.yaml)
- [starling/configs/sequence_encoder/sequence_encoder.yaml](starling/configs/sequence_encoder/sequence_encoder.yaml)
- [starling/data/data_wrangler.py](starling/data/data_wrangler.py)
- [starling/models/attention.py](starling/models/attention.py)
- [starling/models/blocks.py](starling/models/blocks.py)
- [starling/models/diffusion.py](starling/models/diffusion.py)
- [starling/models/transformer.py](starling/models/transformer.py)
- [starling/models/vae_components.py](starling/models/vae_components.py)
- [starling/models/vit.py](starling/models/vit.py)
- [starling/training/diffusion_train.py](starling/training/diffusion_train.py)

</details>



The Diffusion Model in STARLING is a discrete-time Denoising Diffusion Probabilistic Model (DDPM) that operates within the latent space defined by the [Variational Autoencoder (VAE)](#2.1). It is responsible for generating compressed distance map representations conditioned on protein sequences and environmental factors like ionic strength.

## Overview of the Diffusion Process

The `DiffusionModel` class, implemented in [starling/models/diffusion.py:55-188](), manages the transition between noise and structured latent representations. It utilizes a forward process to add noise to VAE latents and a learned reverse process to denoise them, guided by a Vision Transformer (ViT) backbone.

### Forward and Reverse Mathematics
The model supports multiple beta schedules to control the noise injection rate: `linear`, `cosine`, and `sigmoid` [starling/models/diffusion.py:65-69](). These schedules populate buffers such as `alphas_cumprod` and `sqrt_one_minus_alphas_cumprod` used during training and sampling [starling/models/diffusion.py:162-184]().

### Latent Space Scaling
To ensure the latent space has unit variance (as recommended in Rombach et al., 2021), the model initializes a `latent_space_scaling_factor` [starling/models/diffusion.py:158-160](). During the first training step, this factor is calculated by computing the standard deviation of the initial batch of latents across all distributed ranks [starling/models/diffusion.py:228-243]().

## Architecture and Conditioning

The diffusion process is driven by a `ViT` backbone that processes the latent "images" and is conditioned on sequence information via a `SequenceEncoder`.

### DDPM System Components
The following diagram illustrates the relationship between the core diffusion classes and the data flow.

**DDPM Architecture Flow**
```mermaid
graph TD
    subgraph "Conditioning Space"
        A["Sequence (Amino Acids)"] --> B["SequenceEncoder"]
        IS["Ionic Strength"] --> B
    end

    subgraph "Latent Diffusion Space"
        Noise["Gaussian Noise (z_t)"] --> ViT["ViT (Backbone)"]
        B -- "Cross-Attention" --> ViT
        T["Timestep (t)"] -- "SinusoidalPosEmb" --> ViT
        ViT --> Pred["Predicted Noise (ε_θ)"]
    end

    subgraph "Code Entities"
        B["starling.models.transformer.SequenceEncoder"]
        ViT["starling.models.vit.ViT"]
        DDPM["starling.models.diffusion.DiffusionModel"]
    end

    DDPM -- "wraps" --> ViT
    DDPM -- "wraps" --> B
```
Sources: [starling/models/diffusion.py:135-137](), [starling/models/vit.py:31-95](), [starling/models/transformer.py:236-324]()

### Vision Transformer (ViT) Backbone
The `ViT` backbone treats the VAE latent as a single-channel image. It uses `PatchEmbed` to divide the latent into patches [starling/models/vit.py:9-30]() and processes them through a series of `DiTBlock` (Diffusion Transformer Blocks) [starling/models/vit.py:76-79]().
For details, see [Vision Transformer (ViT) Backbone and Transformer Blocks](#2.2.1).

### Sequence and Ionic Strength Conditioning
The `SequenceEncoder` processes protein sequences using a transformer architecture [starling/models/transformer.py:236-242]().
- **Ionic Strength**: This is treated as a continuous conditioning variable. It is concatenated with the sequence embeddings or processed through dropout during training to enable classifier-free guidance [starling/models/transformer.py:315-320]().
- **Cross-Attention**: The `ViT` backbone interacts with the sequence embeddings through cross-attention layers within each transformer block [starling/models/vit.py:113-114]().

## Training Objectives

The model is trained to minimize the difference between the added noise and the predicted noise.

### Min-SNR Loss Weighting
To improve training stability and convergence, STARLING implements Min-SNR loss weighting [starling/models/diffusion.py:255-271](). This strategy weights the loss at different timesteps based on the Signal-to-Noise Ratio (SNR), capped by a gamma hyperparameter (default 5.0) [starling/models/diffusion.py:80]().

### Training Logic Map
This diagram bridges the mathematical objectives to the specific implementation in the training loop.

**DDPM Training Loop Logic**
```mermaid
graph LR
    Data["Latent (z_0)"] --> Scaling["Apply latent_space_scaling_factor"]
    Scaling --> NoiseInject["q_sample: Add Noise"]
    NoiseInject --> LossCalc["Training Step"]
    
    subgraph "Loss Computation"
        direction TB
        Pred["model.forward()"]
        SNR["Min-SNR Weighting"]
        MSE["F.mse_loss"]
    end
    
    LossCalc --> Pred
    Pred --> MSE
    SNR --> MSE
    MSE --> Opt["Optimizer Update"]

    classDef code font-family:monospace
    Scaling:::code
    NoiseInject:::code
    Pred:::code
    SNR:::code
```
Sources: [starling/models/diffusion.py:202-225](), [starling/models/diffusion.py:255-276]()

## Configuration
The DDPM is configured via `starling/configs/diffusion/diffusion.yaml`. Key parameters include:
- `timesteps`: Total diffusion steps (default 1000) [starling/models/diffusion.py:77]().
- `beta_scheduler`: The noise schedule type [starling/models/diffusion.py:76]().
- `min_snr_loss`: Boolean flag to enable SNR weighting [starling/models/diffusion.py:79]().

For detailed information on the internal transformer blocks and adaptive normalization, see the child page: [Vision Transformer (ViT) Backbone and Transformer Blocks](#2.2.1).

Sources: [starling/models/diffusion.py:55-280](), [starling/models/vit.py:31-122](), [starling/models/transformer.py:236-324](), [starling/training/diffusion_train.py:81-112]()

---
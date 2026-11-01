# STARLING Overview

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [README.md](README.md)
- [docs/index.rst](docs/index.rst)
- [docs/usage/ensemble.rst](docs/usage/ensemble.rst)
- [docs/usage/sequence_encoder.rst](docs/usage/sequence_encoder.rst)
- [pyproject.toml](pyproject.toml)
- [starling/__init__.py](starling/__init__.py)
- [starling/configs.py](starling/configs.py)
- [starling_logo-1.png](starling_logo-1.png)

</details>



STARLING (**conSTruction of intrinsically disordered proteins ensembles efficiently via multi-dimensional Generative models**) is a high-performance machine learning framework designed to predict conformational ensembles of intrinsically disordered proteins (IDPs) and regions (IDRs) directly from their amino acid sequences.

By utilizing a two-stage generative architecture—combining a Variational Autoencoder (VAE) with a Latent Diffusion Model (DDPM)—STARLING can generate thousands of physically plausible protein conformations in seconds [README.md:19-27]().

## System Architecture

The STARLING pipeline operates by encoding amino acid sequences into a latent space, generating diverse latent representations via diffusion, and decoding those representations into pairwise distance maps. These maps are then reconstructed into 3D Cartesian coordinates.

### Logic Flow: Sequence to Ensemble

The following diagram illustrates the transformation of data from a raw sequence string to a complete `Ensemble` object.

**Diagram: Data Transformation Pipeline**
```mermaid
graph TD
    subgraph "Input Space"
        Input["Amino Acid Sequence"]
    end

    subgraph "Encoding (Code: sequence_encoder_backend)"
        Tokenizer["StarlingTokenizer"]
        SeqEnc["SequenceEncoder (ViT)"]
        Ionic["Ionic Strength Conditioning"]
    end

    subgraph "Generative Loop (Code: generate_backend)"
        Noise["Gaussian Noise"]
        DDPM["DiffusionModel (ViT Backbone)"]
        Sampler["DDIM / PLMS Sampler"]
        VAE_Dec["ResNet Decoder"]
    end

    subgraph "Reconstruction (Code: torch_mds / sklearn MDS)"
        DistMap["Distance Map (N x N)"]
        MDS["Multidimensional Scaling"]
        Traj["3D Trajectory (PDB/XTC)"]
    end

    Input --> Tokenizer
    Tokenizer --> SeqEnc
    Ionic --> SeqEnc
    SeqEnc --> DDPM
    Noise --> Sampler
    DDPM --> Sampler
    Sampler --> VAE_Dec
    VAE_Dec --> DistMap
    DistMap --> MDS
    MDS --> Traj
```
**Sources:** [starling/frontend/ensemble_generation.py:8-9](), [starling/inference/model_loading.py:20-40](), [starling/structure/ensemble.py:4-20]()

---

## Key Concepts

### 1. Latent Diffusion for IDPs
Unlike traditional molecular dynamics, STARLING does not simulate physics over time. Instead, it uses a Diffusion Model to denoise random Gaussian noise into latent representations that correspond to valid protein topologies. These latents are decoded by a VAE into distance maps, ensuring the generated ensembles capture the heterogeneous nature of disordered proteins.

### 2. Sequence Conditioning and Ionic Strength
The generation process is conditioned on the protein sequence using a Vision Transformer (ViT)-based `SequenceEncoder`. Furthermore, STARLING explicitly models environmental conditions by accepting **Ionic Strength** (e.g., 20, 150, or 300 mM) as a conditioning variable, allowing the model to predict how ensembles expand or collapse in different salt concentrations [starling/configs.py:23-24]().

### 3. The Ensemble Data Structure
The `Ensemble` class is the primary container for STARLING outputs. It stores distance maps and metadata lazily, only calculating 3D coordinates or biophysical properties (like Radius of Gyration or Hydrodynamic Radius) when requested by the user [docs/usage/ensemble.rst:4-30]().

---

## Major Subsystems

### Inference & Generation
The frontend provides a simple `generate()` function that abstracts away the complexities of model loading and batching. It manages the `ModelManager` singleton to handle lazy weight downloading and `torch.compile` optimizations [starling/__init__.py:8-11](), [starling/configs.py:27-35]().
* For details, see [Quick Start: Generating Ensembles](#1.2).

### 3. Structural Reconstruction (MDS)
Because the generative model outputs distance maps, STARLING includes a high-performance 3D reconstruction engine. It supports both CPU-based `sklearn` MDS and GPU-accelerated `torch_mds` for converting these maps into PDB/XTC trajectories [starling/configs.py:19-20]().

### 4. Search and Indexing
STARLING includes a FAISS-based search engine that allows users to query large databases of protein sequences to find those with similar ensemble-aware embeddings [starling/configs.py:128-156]().

---

## Code Entity Map

This diagram maps high-level system components to their specific implementations in the codebase.

**Diagram: System Component Mapping**
```mermaid
graph LR
    subgraph "User Interface"
        CLI["starling CLI"] -- calls --> API["generate()"]
        API -- uses --> MM["ModelManager"]
    end

    subgraph "Model Components"
        MM -- loads --> VAE["VAE (starling/models/vae.py)"]
        MM -- loads --> DIFF["DiffusionModel (starling/models/diffusion.py)"]
        MM -- loads --> ENC["SequenceEncoder (starling/models/vit.py)"]
    end

    subgraph "Data & Analysis"
        API -- returns --> ENS["Ensemble (starling/structure/ensemble.py)"]
        ENS -- calculates --> BME["BME (starling/structure/bme.py)"]
        ENS -- reconstructs --> MDS["MDS (starling/structure/mds_batch.py)"]
    end
```
**Sources:** [pyproject.toml:53-64](), [starling/inference/model_loading.py:10-70](), [starling/structure/ensemble.py:1-100]()

---

## Getting Started

To begin using STARLING, you can install it via PyPI or clone the repository for development.

* **Installation:** `pip install idptools-starling` [README.md:54]()
* **Basic Usage:** `starling "ACDEF..." -c 400 -r` [README.md:73-77]()

For detailed setup instructions, including environment configuration and weight management, see **[Getting Started: Installation and Configuration](#1.1)**.

For a guide on generating your first ensemble via the Python API, see **[Quick Start: Generating Ensembles](#1.2)**.

**Sources:** [README.md:40-86](), [pyproject.toml:7-43]()

---
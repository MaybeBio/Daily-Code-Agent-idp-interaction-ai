# Glossary

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [pyproject.toml](pyproject.toml)
- [starling/configs.py](starling/configs.py)
- [starling/configs/vae_model/model.yaml](starling/configs/vae_model/model.yaml)
- [starling/inference/constraints.py](starling/inference/constraints.py)
- [starling/models/diffusion.py](starling/models/diffusion.py)
- [starling/models/vae.py](starling/models/vae.py)
- [starling/samplers/ddim_sampler.py](starling/samplers/ddim_sampler.py)
- [starling/samplers/plms_sampler.py](starling/samplers/plms_sampler.py)
- [starling/scripts/starling_pretokenize.py](starling/scripts/starling_pretokenize.py)
- [starling/search/__init__.py](starling/search/__init__.py)
- [starling/search/builder.py](starling/search/builder.py)
- [starling/search/search_engine.py](starling/search/search_engine.py)
- [starling/search/search_utils.py](starling/search/search_utils.py)
- [starling/search/similarity_search.py](starling/search/similarity_search.py)
- [starling/search/store.py](starling/search/store.py)
- [starling/structure/ensemble.py](starling/structure/ensemble.py)
- [starling/training/diffusion_train.py](starling/training/diffusion_train.py)

</details>



This page provides technical definitions and implementation details for core concepts, architectural components, and domain-specific terminology used within the STARLING codebase. It serves as a reference for engineers to map high-level protein physics and generative modeling concepts to specific classes and functions.

## Architectural Components

### VAE (Variational Autoencoder)
The first stage of the STARLING pipeline. It compresses high-dimensional protein distance maps (typically $384 \times 384$) into a compact, continuous latent space (typically $24 \times 24 \times 1$).
*   **Implementation**: Defined in `VAE` class `[starling/models/vae.py:86-104]()`.
*   **Encoder/Decoder**: Uses ResNet-based architectures (e.g., `Resnet18_Encoder`) `[starling/models/vae.py:157-166]()`.
*   **Latent Distribution**: Parameterized via `DiagonalGaussianDistribution` `[starling/models/vae.py:15-15]()`.

### Diffusion Model (DDPM)
The second stage of the pipeline. A Denoising Diffusion Probabilistic Model that learns to generate new samples within the VAE's latent space, conditioned on amino acid sequences and ionic strength.
*   **Implementation**: `DiffusionModel` class `[starling/models/diffusion.py:55-70]()`.
*   **Backbone**: Uses a Vision Transformer (`ViT`) to predict noise `[starling/training/diffusion_train.py:88-88]()`.
*   **Conditioning**: Uses a `SequenceEncoder` to process protein sequences into embeddings for the `ViT` cross-attention layers `[starling/models/diffusion.py:136-136]()`.

### Sampler
The algorithm responsible for the reverse diffusion process (denoising).
*   **DDIM (Denoising Diffusion Implicit Models)**: A deterministic/stochastic non-Markovian sampler that allows for faster generation by skipping timesteps `[starling/samplers/ddim_sampler.py:19-58]()`.
*   **PLMS (Pseudo Linear Multi-Step)**: An accelerated sampler using Adams-Bashforth methods `[starling/samplers/plms_sampler.py:62-72]()`.

## Domain Concepts

### Distance Map
A 2D matrix representing pairwise distances between all $C\alpha$ atoms in a protein chain.
*   **Symmetrization**: Distance maps are enforced to be symmetric ($D_{ij} = D_{ji}$) and have a zero diagonal ($D_{ii} = 0$) using `symmetrize_distance_maps` `[starling/inference/constraints.py:12-38]()`.

### Latent Space Scaling Factor
A normalization constant used to scale VAE latents to unit variance before diffusion training, ensuring the diffusion model operates on a standardized distribution.
*   **Buffer**: Stored as `latent_space_scaling_factor` in the `DiffusionModel` `[starling/models/diffusion.py:158-160]()`.

### Ionic Strength
A conditioning variable representing the salt concentration (in mM) of the environment, which influences the electrostatic interactions and compactness of IDPs.
*   **Default**: 150 mM `[starling/configs.py:23-23]()`.
*   **Processing**: Injected into the `SequenceEncoder` to modulate sequence embeddings `[starling/samplers/ddim_sampler.py:68-70]()`.

## Data Structures & Search

### Ensemble
The primary data container in STARLING, holding a collection of distance maps and the corresponding protein sequence.
*   **Implementation**: `Ensemble` class `[starling/structure/ensemble.py:42-75]()`.
*   **Lazy Computation**: Biophysical properties like Radius of Gyration (`rg_vals`) are computed only when requested `[starling/structure/ensemble.py:106-107]()`.

### Search Engine
A FAISS-backed system for finding similar protein ensembles based on sequence embeddings.
*   **Implementation**: `SearchEngine` class `[starling/search/search_engine.py:93-134]()`.
*   **Candidate**: A single search result containing a global ID (`gid`), score, and metadata `[starling/search/search_utils.py:220-228]()`.
*   **SequenceStore**: An SQLite database mapping GIDs to raw sequences and metadata `[starling/search/store.py:1-10]()`.

## Mapping: Natural Language to Code Space

### Pipeline Data Flow
This diagram maps the conceptual flow of "Sequence to Ensemble" to the specific classes and methods involved.

Title: STARLING Inference Data Flow
```mermaid
graph TD
    Input["Protein Sequence (str)"] --> Tokenizer["StarlingTokenizer.encode()"]
    Tokenizer --> SeqEnc["SequenceEncoder.forward()"]
    SeqEnc --> Diffusion["DiffusionModel (DDPM/DDIM Sampler)"]
    
    subgraph "Latent Space Operations"
    Diffusion -- "Predicts" --> Noise["Latent Noise (torch.Tensor)"]
    Noise -- "Denoise Loop" --> Latent["Sampled Latent (24x24)"]
    end
    
    Latent -- "VAE.decode()" --> DistMap["Distance Map (384x384)"]
    DistMap -- "symmetrize_distance_maps()" --> SymMap["Symmetric DistMap"]
    SymMap -- "Ensemble()" --> Final["Ensemble Object"]
    
    Final -- "generate_3d_coordinates_from_distances()" --> PDB["3D Coordinates (MDS)"]
```
**Sources**: `[starling/samplers/ddim_sampler.py:173-214]()`, `[starling/structure/ensemble.py:77-104]()`, `[starling/inference/constraints.py:12-38]()`.

### Constraint System Architecture
This diagram associates the concept of "Guided Diffusion" with the internal `Constraint` class hierarchy and the sampling loop.

Title: Constraint Application Lifecycle
```mermaid
graph LR
    Sampler["DDIMSampler.sample()"] --> Apply["Constraint.apply()"]
    
    subgraph "Constraint.apply() Logic"
    Apply --> Decode["VAE.decode(latents)"]
    Decode --> Loss["Constraint.compute_loss(dist_maps)"]
    Loss --> Grad["torch.autograd.grad(loss, latents)"]
    Grad --> Update["latents = latents - grad * weight"]
    end
    
    subgraph "Concrete Constraints"
    Loss --- Rg["RgConstraint"]
    Loss --- Hel["HelicityConstraint"]
    Loss --- Dist["DistanceConstraint"]
    end
```
**Sources**: `[starling/inference/constraints.py:41-94]()`, `[starling/inference/constraints.py:203-230]()`, `[starling/samplers/ddim_sampler.py:218-223]()`.

## Key Terms Summary Table

| Term | Code Pointer | Description |
| :--- | :--- | :--- |
| **KLD Weight** | `[starling/models/vae.py:21-33]()` | The weight of the Kullback-Leibler Divergence in the VAE ELBO loss. |
| **MDS** | `[starling/structure/coordinates.py:38-39]()` | Multidimensional Scaling; used to reconstruct 3D coordinates from distance maps. |
| **BME** | `[starling/structure/bme.py:35-35]()` | Bayesian Maximum Entropy; used to reweight ensembles based on experimental data. |
| **Gid** | `[starling/search/search_utils.py:224-224]()` | Global Identifier; a unique integer index for a sequence in the `SequenceStore`. |
| **nprobe** | `[starling/search/search_engine.py:25-25]()` | FAISS parameter defining how many clusters to search during ANN retrieval. |
| **FiLM** | `[starling/models/vae.py:146-148]()` | Feature-wise Linear Modulation; a conditioning technique used in ResNet blocks. |

**Sources**: `[starling/models/vae.py:21-151]()`, `[starling/structure/ensemble.py:1-40]()`, `[starling/search/search_engine.py:1-154]()`.
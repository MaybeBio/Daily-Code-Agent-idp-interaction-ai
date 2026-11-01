# Model Loading and Inference Infrastructure

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [hubconf.py](hubconf.py)
- [starling/configs.py](starling/configs.py)
- [starling/data/ddpm_loader_tar.py](starling/data/ddpm_loader_tar.py)
- [starling/inference/model_loading.py](starling/inference/model_loading.py)

</details>



This page details the infrastructure responsible for managing model lifecycles, weight acquisition, and inference-time optimizations. The system is built around a centralized manager that handles lazy loading, remote weight retrieval from GitHub/Zenodo, and optional `torch.compile` integration for high-performance sampling.

## ModelManager Singleton

The `ModelManager` class [starling/inference/model_loading.py:16-17]() acts as the central orchestrator for the two primary neural components: the **VAE** (Encoder/Decoder) and the **Diffusion Model** (DDPM). It ensures that models are only instantiated when needed and provides a unified interface for loading weights from diverse sources.

### Lazy Loading and Initialization
The manager implements a lazy-loading pattern via `get_models()` [starling/inference/model_loading.py:63-100](). When models are requested, the manager checks if they are already resident in memory; if not, it triggers the full loading and compilation sequence.

**Data Flow for Model Initialization:**

1.  **Path Resolution**: The manager resolves local paths or URLs for both the Encoder and DDPM weights using defaults from `starling.configs` [starling/inference/model_loading.py:36-40]().
2.  **Weight Acquisition**: If a path starts with `http`, the `load_from_path_or_url` helper downloads the file to the local `torch.hub` cache [starling/inference/model_loading.py:24-33]().
3.  **Component Assembly**:
    *   Instantiates a `SequenceEncoder` [starling/inference/model_loading.py:49]().
    *   Loads the `DiffusionModel` using `load_from_checkpoint`, injecting a `ViT` backbone and the `SequenceEncoder` [starling/inference/model_loading.py:50-56]().
    *   Loads the `VAE` using `load_from_checkpoint` [starling/inference/model_loading.py:57-60]().
4.  **Compilation**: If enabled in configuration, the models are passed through the `compile()` method [starling/inference/model_loading.py:95-97]().

### Sources:
- `starling/inference/model_loading.py:16-131`()
- `starling/configs.py:74-91`()

## Weight Management and URL Downloads

STARLING supports automatic weight downloading to facilitate ease of use. The system defaults to specific GitHub release versions but can be overridden via environment variables.

| Configuration Key | Default Value / URL | Purpose |
| :--- | :--- | :--- |
| `DEFAULT_ENCODE_WEIGHTS` | `STARLING_v2.0.0_ViT_VAE_2025_10_14.ckpt` | Default VAE checkpoint name [starling/configs.py:14]() |
| `GITHUB_ENCODER_URL` | `.../v2.0.0/STARLING_v2.0.0_ViT_VAE_2025_10_14.ckpt` | Remote VAE weight source [starling/configs.py:82-84]() |
| `GITHUB_DDPM_URL` | `.../v2.0.0/STARLING_v2.0.0_ViT_DDPM_2025_10_14.ckpt` | Remote DDPM weight source [starling/configs.py:85]() |
| `STARLING_ENCODER_PATH` | Environment Variable | Override for local/remote VAE weights [starling/configs.py:88-90]() |

### Implementation Detail: `load_from_path_or_url`
This internal function in `ModelManager` leverages `torch.hub.download_url_to_file` to manage the local cache directory (typically `~/.cache/torch/hub/checkpoints/`) [starling/inference/model_loading.py:24-33]().

### Sources:
- `starling/configs.py:8-25`()
- `starling/configs.py:81-91`()
- `starling/inference/model_loading.py:24-33`()

## Torch Compile Integration

To accelerate inference, STARLING integrates with `torch.compile`. This is particularly beneficial for the diffusion sampling loop, which requires multiple iterations of the ViT backbone.

### Compilation Logic
The `compile()` method [starling/inference/model_loading.py:102-130]() applies `torch.compile` to specific sub-modules:
- `diffusion_model.model` (the ViT backbone) [starling/inference/model_loading.py:109-111]().
- `encoder_model.decoder` (used to transform latents back to distance maps) [starling/inference/model_loading.py:112-114]().

### Configuration Settings
Compilation behavior is controlled via the `TORCH_COMPILATION` dictionary in `configs.py` [starling/configs.py:27-35]():
- `enabled`: Boolean flag (default `False`).
- `options`: Dictionary passed to `torch.compile`, including `mode` (e.g., `reduce-overhead`), `fullgraph`, and `backend`.

### Sources:
- `starling/inference/model_loading.py:102-131`()
- `starling/configs.py:26-35`()

## PyTorch Hub Entrypoint

The `hubconf.py` file allows users to load STARLING models directly using the standard PyTorch Hub API without manually cloning the repository.

### Function: `starling_model`
This entrypoint [hubconf.py:6-18]() provides a simplified interface:
1.  It instantiates a `ModelManager`.
2.  It calls `model_manager.get_models(device=device)`.
3.  It returns the tuple `(encoder_model, diffusion_model)`.

**Usage Example:**
```python
import torch
encoder, diffusion = torch.hub.load('idptools/starling', 'starling_model', device='cuda')
```

### Sources:
- `hubconf.py:1-19`()

## Infrastructure Architecture Diagrams

### System Entity Mapping: Inference Setup
This diagram bridges the conceptual "Loading" phase to the specific classes and files involved.

```mermaid
graph TD
    subgraph "Natural Language Space"
        User["User / API Call"]
        Weights["Model Weights (CKPT)"]
        Compiled["Optimized Models"]
    end

    subgraph "Code Entity Space"
        Hub["hubconf.py:starling_model"]
        MM["ModelManager [starling/inference/model_loading.py]"]
        Cfg["configs.py [starling/configs]"]
        VAE["VAE [starling/models/vae.py]"]
        DDPM["DiffusionModel [starling/models/diffusion.py]"]
        TC["torch.compile"]
    end

    User --> Hub
    Hub --> MM
    MM -- "reads defaults" --> Cfg
    MM -- "downloads/loads" --> Weights
    Weights --> VAE
    Weights --> DDPM
    MM -- "applies" --> TC
    TC --> Compiled
```
**Sources:** `hubconf.py:6-18`(), `starling/inference/model_loading.py:16-100`(), `starling/configs.py:74-91`()

### Data Flow: Weight Acquisition and Compilation
This diagram illustrates the logic flow inside `ModelManager.get_models`.

```mermaid
flowchart TD
    Start["Call get_models()"] --> CheckLoaded{"Models in memory?"}
    CheckLoaded -- "Yes" --> Return["Return (VAE, DDPM)"]
    CheckLoaded -- "No" --> ResolvePaths["Resolve Paths/URLs from configs.py"]
    
    ResolvePaths --> IsURL{"Is path a URL?"}
    IsURL -- "Yes" --> Download["torch.hub.download_url_to_file"]
    IsURL -- "No" --> LoadCKPT["Load via .load_from_checkpoint()"]
    Download --> LoadCKPT
    
    LoadCKPT --> CheckCompile{"TORCH_COMPILATION['enabled']?"}
    CheckCompile -- "Yes" --> CompileSubmodules["torch.compile(ViT) & torch.compile(Decoder)"]
    CheckCompile -- "No" --> Return
    CompileSubmodules --> Return
```
**Sources:** `starling/inference/model_loading.py:21-61`(), `starling/inference/model_loading.py:63-100`(), `starling/inference/model_loading.py:102-130`()

---
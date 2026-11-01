# VAE Training

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [starling/configs/dataloader/vae_dataloader.yaml](starling/configs/dataloader/vae_dataloader.yaml)
- [starling/configs/trainer/vae_trainer.yaml](starling/configs/trainer/vae_trainer.yaml)
- [starling/configs/vae_configs.yaml](starling/configs/vae_configs.yaml)
- [starling/configs/vae_model/model.yaml](starling/configs/vae_model/model.yaml)
- [starling/data/VAE_loader_tar.py](starling/data/VAE_loader_tar.py)
- [starling/data/argument_parser.py](starling/data/argument_parser.py)
- [starling/models/vae.py](starling/models/vae.py)
- [starling/training/config.yaml](starling/training/config.yaml)
- [starling/training/vae_train.py](starling/training/vae_train.py)

</details>



The Variational Autoencoder (VAE) training pipeline is responsible for learning a compressed latent representation of protein distance maps. This stage is critical for the subsequent Diffusion model training, as it defines the latent space where diffusion occurs. The pipeline is built on PyTorch Lightning and utilizes Hydra for configuration management, supporting large-scale distributed training on protein conformational data.

### Training Entrypoint and Data Flow

The primary entrypoint for VAE training is the `train_vae` function located in `starling/training/vae_train.py` [starling/training/vae_train.py:106-191](). This script orchestrates the instantiation of the model, the data module, and the PyTorch Lightning `Trainer`.

**System Architecture and Data Flow**
The following diagram illustrates the relationship between the configuration files, the training script, and the core VAE classes.

"VAE Training Architecture"
```mermaid
graph TD
    subgraph "Configuration Space (Hydra)"
        VC["vae_configs.yaml"] --> DL["dataloader/vae_dataloader.yaml"]
        VC --> TR["trainer/vae_trainer.yaml"]
        VC --> VM["vae_model/model.yaml"]
    end

    subgraph "Code Entity Space (Execution)"
        TV["train_vae() in vae_train.py"]
        SVM["setup_vae_model()"]
        SDM["setup_data_module()"]
        VAE_CLASS["class VAE in vae.py"]
        VLD["class VAEdataloader in VAE_loader_tar.py"]
    end

    VC -- "Parsed by @hydra.main" --> TV
    VM -- "instantiate(cfg.vae_model)" --> SVM
    SVM -- "returns" --> VAE_CLASS
    DL -- "instantiate(cfg.dataloader)" --> SDM
    SDM -- "returns" --> VLD
    TV -- "trainer.fit(vae, dataset)" --> VAE_CLASS
    VLD -- "yields distance_maps" --> VAE_CLASS
```
Sources: [starling/training/vae_train.py:101-191](), [starling/configs/vae_configs.yaml:1-5](), [starling/models/vae.py:86-104]()

### Configuration Hierarchy

Training is driven by a hierarchical Hydra configuration defined in `starling/configs/vae_configs.yaml` [starling/configs/vae_configs.yaml:1-5]().

| Config Component | File Path | Key Responsibilities |
| :--- | :--- | :--- |
| **Model** | `vae_model/model.yaml` | Architecture (ResNet type), latent dimensions, loss types (`mse`/`nll`), and KLD scheduling [starling/configs/vae_model/model.yaml:1-17](). |
| **Trainer** | `trainer/vae_trainer.yaml` | GPU/Node counts, precision (`bf16-mixed`), checkpoint paths, and fine-tuning flags [starling/configs/trainer/vae_trainer.yaml:1-12](). |
| **Dataloader** | `dataloader/vae_dataloader.yaml` | Dataset paths, batch size, and worker configuration for `tar` or `h5` backends [starling/configs/dataloader/vae_dataloader.yaml:1-17](). |

### VAE Model Initialization and Fine-Tuning

The `setup_vae_model` function in `vae_train.py` handles the logic for creating a new model or loading weights from an existing checkpoint [starling/training/vae_train.py:79-98]().

*   **Fresh Training**: If `fine_tune` is `false` and no checkpoint is provided, the model is instantiated from the Hydra config [starling/training/vae_train.py:96]().
*   **Fine-Tuning**: If `fine_tune` is `true`, the model is first instantiated with the current configuration, and then weights are loaded from the checkpoint via `load_state_dict` [starling/training/vae_train.py:83-93]().
*   **Resuming**: If `fine_tune` is `false` but a checkpoint exists, the `Trainer.fit` method handles resuming the full training state (optimizer, epoch count) [starling/training/vae_train.py:141-148](), [starling/training/vae_train.py:187-191]().

Sources: [starling/training/vae_train.py:79-98](), [starling/training/vae_train.py:187-191]()

### Distributed Training and Batch Size

Effective batch size is calculated to ensure consistency across distributed environments. For the `tar` (WebDataset) loader, the effective batch size is determined by:
`effective_batch_size = cfg.trainer.cuda * cfg.trainer.num_nodes * cfg.dataloader.tar.batch_size` [starling/training/vae_train.py:152-154]().

The `VAEdataloader` uses this value to calculate the number of training and validation batches per epoch, which is required for `webdataset` to provide a consistent epoch length [starling/data/VAE_loader_tar.py:47-49]().

Sources: [starling/training/vae_train.py:152-156](), [starling/data/VAE_loader_tar.py:47-49]()

### KLD Weight Scheduling

The VAE employs a `KLDWeightScheduler` to manage the trade-off between reconstruction accuracy and latent space regularity (the KL divergence term in the ELBO loss) [starling/models/vae.py:21-84]().

*   **Cyclical Annealing**: The default scheduler type is `cyclical` [starling/models/vae.py:55](). It divides the total training steps into cycles (defaulting to 5 cycles, each 20% of total steps) [starling/models/vae.py:57]().
*   **Warmup**: Within each cycle, the weight ramps linearly from 0 to `max_weight` during the `warmup_fraction` of the cycle [starling/models/vae.py:60-67]().

"VAE Training Logic and Loss Calculation"
```mermaid
sequenceDiagram
    participant D as VAEdataloader
    participant V as VAE (LightningModule)
    participant S as KLDWeightScheduler
    participant W as WandB

    D->>V: batch (distance_maps)
    V->>V: encode(distance_maps) -> posterior
    V->>V: decode(latent_samples) -> reconstruction
    V->>S: get_weight(current_step)
    S-->>V: current_kld_weight
    V->>V: compute ELBO (recon_loss + kld_weight * kld_loss)
    V->>W: log(epoch_val_loss, current_kld_weight)
```
Sources: [starling/models/vae.py:21-84](), [starling/models/vae.py:236-267](), [starling/training/vae_train.py:168-169]()

### Checkpoint and Experiment Management

The pipeline integrates with **Weights & Biases (WandB)** for experiment tracking.
*   **Initialization**: `wandb_init` is called on `rank_zero` to prevent duplicate runs in distributed settings [starling/training/vae_train.py:19-21]().
*   **Callbacks**: The `ModelCheckpoint` callback is configured to monitor `epoch_val_loss` and save the best model as well as a `last.ckpt` for resiliency [starling/training/vae_train.py:128-138]().
*   **Logging**: The `LearningRateMonitor` tracks scheduler behavior across steps [starling/training/vae_train.py:180]().

Sources: [starling/training/vae_train.py:19-21](), [starling/training/vae_train.py:128-138](), [starling/training/vae_train.py:168-184]()

---
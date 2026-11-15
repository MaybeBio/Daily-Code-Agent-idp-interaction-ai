#!/usr/bin/env python3
"""
Train a sparse autoencoder on pre-computed embeddings.

Usage:
    python scripts/03_train_sae.py paths=colab model=esm2_8m sae=l1_32x
    python scripts/03_train_sae.py paths=colab model=esm2_650m sae=topk_8x model.layer=33
"""

import logging
import sys
from pathlib import Path

import hydra
import numpy as np
import torch
from omegaconf import DictConfig, OmegaConf

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from idpinterp.models.memmap_store import load_mmap_config
from idpinterp.sae.architecture import build_sae
from idpinterp.sae.dataset import BufferedShuffleDataset
from idpinterp.sae.trainer import SAETrainer

logger = logging.getLogger(__name__)

try:
    import wandb
    HAS_WANDB = True
except ImportError:
    HAS_WANDB = False


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
    )

    # Determine which layer to train on
    # If model.layer is specified, use it. Otherwise, default to the last layer.
    layer = cfg.model.get("layer")
    if layer is None:
        layer = cfg.model.layers_to_extract[-1]

    _tag = cfg.sae.get("tag", None)
    _tag_suffix = f"_{_tag}" if _tag else ""
    run_name = f"{cfg.model.name}_{cfg.sae.type}_{cfg.sae.expansion_factor}x{_tag_suffix}_layer{layer}"
    logger.info(f"=== SAE training [{run_name}] ===")

    # Idempotent sweep: if this config already finished, skip (a partially-trained
    # config with only step checkpoints will still resume below).
    final_ckpt = Path(cfg.paths.checkpoint_dir) / f"{run_name}_final.pt"
    if final_ckpt.exists() and not cfg.sae.get("force_retrain", False):
        logger.info(
            f"Final checkpoint exists: {final_ckpt}. Skipping training "
            f"(set sae.force_retrain=true to override)."
        )
        return

    # Set seeds
    torch.manual_seed(cfg.seed)
    np.random.seed(cfg.seed)
    if torch.cuda.is_available():
        torch.backends.cudnn.deterministic = True

    # Determine device
    device = "cuda" if torch.cuda.is_available() else "cpu"

    # Load memmap config
    embed_dir = Path(cfg.paths.embedding_dir) / cfg.model.name / cfg.dataset.name
    mmap_path = embed_dir / f"layer{layer}.mmap"
    mmap_cfg = load_mmap_config(mmap_path)
    total_residues = mmap_cfg["total_residues"]
    hidden_dim = mmap_cfg["hidden_dim"]

    logger.info(f"Embedding file: {mmap_path}")
    logger.info(f"  Residues: {total_residues:,}, Hidden dim: {hidden_dim}")

    # Initialize dataset
    dataset = BufferedShuffleDataset(
        mmap_path, 
        total_residues, 
        hidden_dim, 
        batch_size=cfg.sae.batch_size
    )

    # Build SAE
    sae = build_sae(cfg.sae, cfg.model)
    n_params = sum(p.numel() for p in sae.parameters())
    logger.info(f"SAE: {cfg.sae.type}, {cfg.sae.expansion_factor}x, {n_params:,} params")

    # Initialize W&B
    if cfg.wandb.get("enabled", True) and HAS_WANDB:
        wandb.init(
            project=cfg.wandb.project,
            entity=cfg.wandb.get("entity"),
            name=run_name,
            config=OmegaConf.to_container(cfg, resolve=True),
        )

    # Check for existing checkpoint to resume from
    ckpt_dir = Path(cfg.paths.checkpoint_dir)
    resume_from = None
    existing_ckpts = sorted(ckpt_dir.glob(f"{run_name}_step*.pt"))
    if existing_ckpts:
        resume_from = existing_ckpts[-1]
        logger.info(f"Found existing checkpoint: {resume_from}")

    # Train
    trainer = SAETrainer(
        batch_size=cfg.sae.batch_size,
        lr=cfg.sae.lr,
        total_steps=cfg.sae.total_steps,
        checkpoint_every=cfg.sae.checkpoint_every,
        log_every=cfg.sae.log_every,
        dead_feature_window=cfg.sae.dead_feature_window,
        lr_warmup_frac=cfg.sae.get("lr_warmup_frac", 0.05),
        checkpoint_dir=ckpt_dir,
        use_wandb=cfg.wandb.get("enabled", True),
        device=device,
        resample_dead=cfg.sae.get("resample_dead", False),
        resample_every=cfg.sae.get("resample_every", 25_000),
    )

    trained_sae = trainer.train(
        sae=sae,
        dataset=dataset,
        run_name=run_name,
        resume_from=resume_from,
    )

    if HAS_WANDB and wandb.run is not None:
        wandb.finish()

    logger.info(f"=== Training complete: {run_name} ===")


if __name__ == "__main__":
    main()

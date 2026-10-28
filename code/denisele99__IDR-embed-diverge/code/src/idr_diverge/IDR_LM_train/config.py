from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import yaml


#Config utilities for IDR_LM pretraining script.
#Dataclass structures and config loading logic were generated with assistance from ChatGPT (OpenAI) and modified by Denise Le.


@dataclass
class TrainingConfig:
    # Optimization
    batch_size: int
    test_size: float
    num_train_epochs: int
    warmup_proportion: float
    learning_rate: float
    no_cuda: bool
    seed: int
    gradient_accumulation_steps: int

    # Tokenizer / MLM settings
    max_seq_length: int
    masked_lm_prob: float
    max_predictions_per_seq: int

    # Logging
    wandb_project: str
    wandb_run_name: str

    # Paths
    train_data_path: Path
    model_config_path: Path
    resume_checkpoint_path: Optional[Path]
    checkpoint_output_dir: Path
    final_model_dir: Path


def load_config(config_path: str | Path) -> TrainingConfig:
    config_path = Path(config_path)

    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    return TrainingConfig(
        # Optimization
        batch_size=cfg.get("batch_size", 64),
        test_size=cfg.get("test_size", 0.15),
        num_train_epochs=cfg.get("num_train_epochs", 30),
        warmup_proportion=cfg.get("warmup_proportion", 0.1),
        learning_rate=float(cfg.get("learning_rate", 1e-4)),
        no_cuda=cfg.get("no_cuda", False),
        seed=cfg.get("seed", 42),
        gradient_accumulation_steps=cfg.get("gradient_accumulation_steps", 1),

        # Tokenizer / MLM settings
        max_seq_length=cfg.get("max_seq_length", 512),
        masked_lm_prob=cfg.get("masked_lm_prob", 0.05),
        max_predictions_per_seq=cfg.get("max_predictions_per_seq", 50),

        # Logging
        wandb_project=cfg.get("wandb_project", "idr_lm"),
        wandb_run_name=cfg.get("wandb_run_name", "default_run"),

        # Paths
        train_data_path=Path(cfg["train_data_path"]),
        model_config_path=Path(cfg["model_config_path"]),
        resume_checkpoint_path=(
            Path(cfg["resume_checkpoint_path"])
            if cfg.get("resume_checkpoint_path")
            else None
        ),
        checkpoint_output_dir=Path(cfg["checkpoint_output_dir"]),
        final_model_dir=Path(cfg["final_model_dir"]),
    )
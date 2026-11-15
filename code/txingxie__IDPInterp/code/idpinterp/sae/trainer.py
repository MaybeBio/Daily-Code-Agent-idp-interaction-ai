import logging
import os
import time
from pathlib import Path
from typing import Optional

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from tqdm import tqdm

try:
    import wandb
    HAS_WANDB = True
except ImportError:
    HAS_WANDB = False

from idpinterp.sae.architecture import BaseSAE

logger = logging.getLogger(__name__)


class SAETrainer:
    """Trainer for Sparse Autoencoders."""

    def __init__(
        self,
        batch_size: int = 2048,
        lr: float = 1e-4,
        total_steps: int = 500_000,
        checkpoint_every: int = 10_000,
        log_every: int = 100,
        dead_feature_window: int = 5000,
        lr_warmup_frac: float = 0.05,
        checkpoint_dir: str | Path = "outputs/checkpoints",
        use_wandb: bool = True,
        device: str = "cuda",
        resample_dead: bool = False,
        resample_every: int = 25_000,
    ):
        self.batch_size = batch_size
        self.lr = lr
        self.total_steps = total_steps
        self.checkpoint_every = checkpoint_every
        self.log_every = log_every
        self.dead_feature_window = dead_feature_window
        self.lr_warmup_frac = lr_warmup_frac
        # Opt-in dead-feature resampling (Anthropic-style). Off by default; the
        # TopK sweep usually keeps far more features alive, so enable only if a
        # config still shows >~10% dead.
        self.resample_dead = resample_dead
        self.resample_every = resample_every
        self.checkpoint_dir = Path(checkpoint_dir)
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        self.use_wandb = use_wandb and HAS_WANDB
        self.device = device

    def _get_lr_multiplier(self, step: int) -> float:
        """Linear warmup for learning rate."""
        warmup_steps = int(self.total_steps * self.lr_warmup_frac)
        if step < warmup_steps:
            return step / max(1, warmup_steps)
        return 1.0

    def _get_l1_multiplier(self, step: int) -> float:
        """Linear warmup for L1 penalty to allow features to initialize."""
        warmup_steps = int(self.total_steps * self.lr_warmup_frac)
        if step < warmup_steps:
            return step / max(1, warmup_steps)
        return 1.0

    @torch.no_grad()
    def _resample_dead_features(self, sae, optimizer, batch, steps_since_active) -> int:
        """
        Reinitialise features that have been dead longer than the window toward
        directions the SAE currently reconstructs poorly (high residual norm),
        and reset their Adam moments. Returns the number resampled.

        Weight shapes (see architecture.py): W_enc (d_in, d_latent),
        W_dec (d_latent, d_in), b_enc (d_latent,).
        """
        dead = (steps_since_active > self.dead_feature_window).nonzero(as_tuple=True)[0]
        if dead.numel() == 0:
            return 0

        x_hat, _ = sae(batch)
        resid = batch - x_hat
        err = (resid ** 2).sum(dim=1)
        if float(err.sum()) <= 0:
            return 0

        # Sample candidate input directions with prob ∝ squared reconstruction error
        probs = err / err.sum()
        n = min(int(dead.numel()), batch.shape[0])
        pick = torch.multinomial(probs, num_samples=n, replacement=True)
        dirs = resid[pick]
        dirs = dirs / (dirs.norm(dim=1, keepdim=True) + 1e-8)         # (n, d_in) unit
        dead = dead[:n]

        alive = steps_since_active <= self.dead_feature_window
        enc_norm = (
            sae.W_enc.data[:, alive].norm(dim=0).mean().item() if alive.any() else 1.0
        )

        sae.W_dec.data[dead] = dirs                                    # (n, d_in)
        sae.W_enc.data[:, dead] = (dirs * (0.2 * enc_norm)).T          # (d_in, n)
        sae.b_enc.data[dead] = 0.0

        # Reset Adam moments for the touched slices so they get a fresh start
        for p, sl in (
            (sae.W_enc, (slice(None), dead)),
            (sae.W_dec, (dead,)),
            (sae.b_enc, (dead,)),
        ):
            st = optimizer.state.get(p)
            if st:
                for key in ("exp_avg", "exp_avg_sq"):
                    if key in st:
                        st[key][sl] = 0.0

        steps_since_active[dead] = 0
        return int(n)

    def train(
        self,
        sae: BaseSAE,
        dataset: torch.utils.data.Dataset,
        run_name: str,
        resume_from: Optional[str | Path] = None,
    ) -> BaseSAE:
        
        sae = sae.to(self.device)
        optimizer = optim.Adam(sae.parameters(), lr=self.lr)
        
        # Load checkpoint if resuming
        start_step = 0
        if resume_from:
            logger.info(f"Resuming from checkpoint: {resume_from}")
            ckpt = torch.load(resume_from, map_location=self.device)
            sae.load_state_dict(ckpt["model_state_dict"])
            optimizer.load_state_dict(ckpt["optimizer_state_dict"])
            start_step = ckpt["step"]
            logger.info(f"Resumed at step {start_step}")

        # Infinite data loader
        is_iterable = isinstance(dataset, torch.utils.data.IterableDataset)
        if is_iterable:
            loader = DataLoader(
                dataset,
                batch_size=None,  # Dataset already yields batches
                num_workers=2,
                pin_memory=True,
            )
        else:
            loader = DataLoader(
                dataset, 
                batch_size=self.batch_size, 
                shuffle=True, 
                num_workers=2,
                drop_last=True,
                pin_memory=True,
            )
        def cycle(iterable):
            while True:
                for x in iterable:
                    yield x
        iterator = cycle(loader)

        # Track dead features (uses d_hidden property which returns d_latent)
        steps_since_active = torch.zeros(sae.d_hidden, device=self.device)
        
        # Use Automatic Mixed Precision (AMP) for faster training on GPUs
        use_amp = self.device == "cuda"
        scaler = torch.amp.GradScaler(enabled=use_amp)

        sae.train()
        pbar = tqdm(total=self.total_steps, initial=start_step, desc="Training SAE")

        for step in range(start_step, self.total_steps):
            t0 = time.perf_counter()
            batch = next(iterator).to(self.device)
            t1 = time.perf_counter()
            data_ms = (t1 - t0) * 1000.0

            # Schedulers
            lr_mult = self._get_lr_multiplier(step)
            for pg in optimizer.param_groups:
                pg["lr"] = self.lr * lr_mult
                
            l1_mult = self._get_l1_multiplier(step)

            optimizer.zero_grad()

            t2 = time.perf_counter()

            # Dead-feature mask from the PRIOR steps' activity (current-step
            # deadness isn't known yet). TopKSAE uses it for the AuxK loss.
            dead_mask = steps_since_active > self.dead_feature_window

            with torch.amp.autocast(device_type="cuda", enabled=use_amp):
                # sae.loss() is defined in BaseSAE — returns (loss_dict, activations)
                loss_dict, acts = sae.loss(batch, l1_mult=l1_mult, dead_mask=dead_mask)
                loss = loss_dict["loss"]

            scaler.scale(loss).backward()
            
            # Constrain decoder weights to unit norm before the optimizer step
            sae.normalize_decoder()

            scaler.step(optimizer)
            scaler.update()
            
            t3 = time.perf_counter()
            compute_ms = (t3 - t2) * 1000.0

            # Track dead features; GPU-only ops, no host sync
            with torch.no_grad():
                fired = acts.sum(dim=0) > 1e-5
                steps_since_active[fired] = 0
                steps_since_active[~fired] += 1

            # Opt-in dead-feature resampling
            if (
                self.resample_dead
                and step > 0
                and step % self.resample_every == 0
            ):
                n_res = self._resample_dead_features(
                    sae, optimizer, batch, steps_since_active
                )
                if n_res:
                    logger.info(f"Step {step}: resampled {n_res} dead features.")

            pbar.update(1)

            # Logging (compute the host-sync'd metrics only here, not every step)
            if step % self.log_every == 0:
                l0 = (acts > 1e-5).float().sum(dim=1).mean().item()
                dead_count = (steps_since_active > self.dead_feature_window).sum().item()
                fve = loss_dict.get("fve", 0.0)
                
                metrics = {
                    "step": step,
                    "loss/total": loss.item(),
                    "loss/mse": loss_dict["mse_loss"],
                    "loss/l1": loss_dict.get("l1_loss", 0.0),
                    "loss/aux": loss_dict.get("aux_loss", 0.0),
                    "metrics/l0": l0,
                    "metrics/dead_features": dead_count,
                    "metrics/fve": fve,
                    "hyperparams/lr_mult": lr_mult,
                    "hyperparams/l1_mult": l1_mult,
                    "profiling/data_ms": data_ms,
                    "profiling/compute_ms": compute_ms,
                }
                
                if self.use_wandb:
                    wandb.log(metrics, step=step)
                    
                pbar.set_postfix(
                    mse=f"{loss_dict['mse_loss']:.3f}",
                    l0=f"{l0:.1f}",
                    fve=f"{fve:.3f}",
                    dead=dead_count,
                    aux=f"{loss_dict.get('aux_loss', 0.0):.3f}",
                    data_ms=f"{data_ms:.1f}",
                    comp_ms=f"{compute_ms:.1f}"
                )

            # Checkpointing
            if step > 0 and step % self.checkpoint_every == 0:
                ckpt_path = self.checkpoint_dir / f"{run_name}_step{step}.pt"
                torch.save({
                    "step": step,
                    "model_state_dict": sae.state_dict(),
                    "optimizer_state_dict": optimizer.state_dict(),
                }, ckpt_path)
                
                # Delete older checkpoints to save space (keep only last 2)
                old_ckpts = sorted(self.checkpoint_dir.glob(f"{run_name}_step*.pt"))
                if len(old_ckpts) > 2:
                    for old in old_ckpts[:-2]:
                        old.unlink()

        # Save final model
        final_path = self.checkpoint_dir / f"{run_name}_final.pt"
        torch.save({
            "step": self.total_steps,
            "model_state_dict": sae.state_dict(),
            "optimizer_state_dict": optimizer.state_dict(),
        }, final_path)
        logger.info(f"Final model saved to {final_path}")

        pbar.close()
        return sae

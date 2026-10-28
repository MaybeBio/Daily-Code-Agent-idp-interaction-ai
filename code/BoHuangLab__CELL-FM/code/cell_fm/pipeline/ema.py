# -*- coding: utf-8 -*-
import json
import os

import torch
from torch.nn.parallel import DistributedDataParallel

from diffusers.training_utils import EMAModel
from transformers import TrainerCallback

from cell_fm.logging import logger

EMA_WEIGHTS_NAME = "ema_pytorch_model.bin"
EMA_STATE_NAME = "ema_state.json"


def unwrap_model(model):
    while isinstance(model, DistributedDataParallel):
        model = model.module
    return model


def add_ema_callback(trainer, callback):
    """Register an `EMACallback` ahead of the reporting callbacks.

    `Trainer` appends user callbacks after the wandb/reporting ones, which would read `logs` in
    `on_log` before the callback adds `ema_decay` to it.
    """
    trainer.add_callback(callback)
    callbacks = trainer.callback_handler.callbacks
    callbacks.remove(callback)
    callbacks.insert(0, callback)
    return callback


class EMACallback(TrainerCallback):
    """Track an exponential moving average of the trainable weights.

    The averaged weights are written next to each checkpoint as `ema_pytorch_model.bin`, a full
    state dict interchangeable with `pytorch_model.bin`, so inference loads them through the
    ordinary `--infer --loadcheck_path` path with no code change.

    Only parameters with `requires_grad=True` are averaged; frozen submodules (e.g. a pretrained
    VAE or ESM trunk) are copied through from the live model.
    """

    def __init__(
        self,
        decay: float = 0.9999,
        min_decay: float = 0.0,
        update_after_step: int = 0,
        use_ema_warmup: bool = True,
        inv_gamma: float = 1.0,
        power: float = 2 / 3,
        foreach: bool = True,
    ):
        self.decay = decay
        self.min_decay = min_decay
        self.update_after_step = update_after_step
        self.use_ema_warmup = use_ema_warmup
        self.inv_gamma = inv_gamma
        self.power = power
        self.foreach = foreach

        self.ema = None
        self.param_names = None
        self.params = None

    def _checkpoint_dir(self, args, state) -> str:
        return os.path.join(args.output_dir, f"checkpoint-{state.global_step}")

    def on_train_begin(self, args, state, control, model=None, **kwargs):
        model = unwrap_model(model)

        named = [(n, p) for n, p in model.named_parameters() if p.requires_grad]
        self.param_names = [n for n, _ in named]
        self.params = [p for _, p in named]

        self.ema = EMAModel(
            self.params,
            decay=self.decay,
            min_decay=self.min_decay,
            update_after_step=self.update_after_step,
            use_ema_warmup=self.use_ema_warmup,
            inv_gamma=self.inv_gamma,
            power=self.power,
            foreach=self.foreach,
        )

        num_params = sum(p.numel() for p in self.params)
        logger.info(
            f"EMA enabled | decay={self.decay} warmup={self.use_ema_warmup} "
            f"inv_gamma={self.inv_gamma} power={self.power:.4f} "
            f"update_after_step={self.update_after_step}"
        )
        logger.info(
            f"EMA tracking {len(self.params)} tensors / {num_params:,} trainable params"
        )

        # A non-zero global_step means the Trainer restored trainer_state.json, i.e. this is a
        # resume rather than a fresh run or a --ft start.
        if state.global_step > 0:
            self._restore(self._checkpoint_dir(args, state))

    def _restore(self, checkpoint_dir: str):
        weights_path = os.path.join(checkpoint_dir, EMA_WEIGHTS_NAME)
        state_path = os.path.join(checkpoint_dir, EMA_STATE_NAME)

        if not (os.path.isfile(weights_path) and os.path.isfile(state_path)):
            logger.warning(
                f"EMA resume: {EMA_WEIGHTS_NAME}/{EMA_STATE_NAME} not found in {checkpoint_dir}; "
                "starting EMA fresh from the current weights"
            )
            return

        with open(state_path) as f:
            self.ema.load_state_dict(json.load(f))

        checkpoint_state = torch.load(weights_path, map_location="cpu")
        missing = [n for n in self.param_names if n not in checkpoint_state]
        if missing:
            logger.warning(
                f"EMA resume: {len(missing)} tracked params absent from {weights_path} "
                f"(e.g. {missing[:3]}); starting EMA fresh from the current weights"
            )
            return

        self.ema.shadow_params = [
            checkpoint_state[n].to(device=p.device, dtype=p.dtype)
            for n, p in zip(self.param_names, self.params)
        ]
        logger.info(
            f"EMA resumed from {checkpoint_dir} at optimization_step="
            f"{self.ema.optimization_step}"
        )

    def on_step_end(self, args, state, control, **kwargs):
        self.ema.step(self.params)

    def on_log(self, args, state, control, logs=None, **kwargs):
        if logs is not None and self.ema is not None and self.ema.cur_decay_value is not None:
            logs["ema_decay"] = self.ema.cur_decay_value

    def on_save(self, args, state, control, model=None, **kwargs):
        if not args.should_save:
            return

        model = unwrap_model(model)
        checkpoint_dir = self._checkpoint_dir(args, state)
        os.makedirs(checkpoint_dir, exist_ok=True)

        # Keys come from the live model, so they land in its own namespace and stay loadable by
        # the name-filtered `load_pretrained_weights`. Frozen params/buffers pass through as-is.
        checkpoint_state = {k: v.detach().cpu() for k, v in model.state_dict().items()}
        for name, shadow in zip(self.param_names, self.ema.shadow_params):
            checkpoint_state[name] = shadow.detach().cpu()
        torch.save(checkpoint_state, os.path.join(checkpoint_dir, EMA_WEIGHTS_NAME))

        ema_state = self.ema.state_dict()
        ema_state.pop("shadow_params")
        with open(os.path.join(checkpoint_dir, EMA_STATE_NAME), "w") as f:
            json.dump(ema_state, f, indent=2)

        logger.info(
            f"Saved EMA weights to {os.path.join(checkpoint_dir, EMA_WEIGHTS_NAME)} "
            f"(decay={self.ema.cur_decay_value})"
        )

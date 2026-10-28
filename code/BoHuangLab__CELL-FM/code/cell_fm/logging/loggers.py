# -*- coding: utf-8 -*-
import sys
from dataclasses import dataclass, fields, is_dataclass
from typing import Dict, Union

import torch
from loguru import logger

from cell_fm.utils.dist_utils import is_master_node
from transformers import TrainerCallback

# wandb is only used by the training callback below, so it is optional rather than a
# required install: the Colab notebooks import cell_fm through this module and Colab has
# no wandb, which would otherwise make every notebook fail at its first import.
try:
    import wandb  # isort:skip
except ImportError:  # pragma: no cover
    wandb = None

handlers = {}


def get_logger():
    if not handlers:
        logger.remove()  # remove default handler
        handlers["console"] = logger.add(
            sys.stdout,
            format="[<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green>][<cyan>{level}</cyan>]: {message}",
            colorize=True,
            filter=console_log_filter,
            enqueue=True,
        )

    return logger


# Custom function to handle tensor attributes
def dataclass_to_dict(dataclass_obj: Union[dataclass, Dict]) -> Dict:
    if isinstance(dataclass_obj, dict):
        return dataclass_obj
    result = {}
    for field in fields(dataclass_obj):
        value = getattr(dataclass_obj, field.name)
        if isinstance(value, torch.Tensor) and not value.is_leaf:
            result[field.name] = value.clone().detach()
        else:
            result[field.name] = value
    return result


class MetricLogger(object):
    def log(self, metrics, prefix=""):
        if not is_master_node():
            return

        log_data = {}
        if type(metrics) is dict:
            log_data = metrics
        elif is_dataclass(metrics):
            log_data = dataclass_to_dict(metrics)

            if "extra_output" in log_data:
                extra_output = log_data["extra_output"]
                if extra_output is not None:
                    for k, v in extra_output.items():
                        log_data[k] = v
                del log_data["extra_output"]

        for k in log_data:
            if isinstance(log_data[k], torch.Tensor):
                log_data[k] = log_data[k].detach().item()

        logger.info(" | ".join([f"{k}={v:.4g}" for k, v in log_data.items()]))
        if wandb is not None and wandb.run is not None:
            # Add prefix
            if prefix:
                log_data = {f"{prefix}/{k}": v for k, v in log_data.items()}
            wandb.log(log_data)


def console_log_filter(record):
    # For message with level INFO, we only log it on master node
    # For others, we log it on all nodes
    if record["level"].name != "INFO":
        return True

    return is_master_node()


class StableDiffusionLoggingCallback(TrainerCallback):
    def __init__(self):
        self.logger = logger
    
    def on_log(self, args, state, control, logs=None, **kwargs):
        if logs is None:
            return
        
        formatted_log = (
            f"loss={logs.get('loss', 0):.4f} | grad_norm={logs.get('grad_norm', 0):.4f} | "
            f"lr={logs.get('learning_rate', 0):.7f} | epoch={logs.get('epoch', 0):.4f} | "
            f"global_step={state.global_step}"
        )
        
        self.logger.info(formatted_log)


class CELLFMLoggingCallback(TrainerCallback):
    def __init__(self):
        self.logger = logger

    def on_log(self, args, state, control, logs=None, **kwargs):
        if not logs:
            return

        # key -> format string
        field_fmt = {
            "loss": "{:.4f}",
            "sequence_loss": "{:.4f}",
            "image_diff_loss": "{:.4f}",
            "image_recon_loss": "{:.4f}",
            "grad_norm": "{:.4f}",
            "learning_rate": "{:.7f}",
            "epoch": "{:.4f}",
            "ema_decay": "{:.6f}",
        }

        parts = []
        for k, fmt in field_fmt.items():
            if k in logs and logs[k] is not None:
                v = logs[k]
                parts.append(f"{k}={fmt.format(float(v))}")

        parts.append(f"global_step={state.global_step}")

        formatted_log = " | ".join(parts)
        self.logger.info(formatted_log)


class CELLFMClsLoggingCallback(TrainerCallback):
    def __init__(self):
        self.logger = logger
    
    def on_log(self, args, state, control, logs=None, **kwargs):
        if logs is None:
            return

        loss = logs.get("loss", 0.0)
        image_class_loss = logs.get("image_class_loss", 0.0)
        grad_norm = logs.get("grad_norm", 0.0)
        learning_rate = logs.get("learning_rate", 0.0)
        epoch = logs.get("epoch", 0.0)

        formatted_log = (
            f"loss={loss:.4f} | image_class_loss={image_class_loss:.4f} | "
            f"grad_norm={grad_norm:.4f} | lr={learning_rate:.7f} | "
            f"epoch={epoch:.4f} | global_step={state.global_step}"
        )

        self.logger.info(formatted_log)


class VAELoggingCallback(TrainerCallback):
    def __init__(self):
        self.logger = logger
    
    def on_log(self, args, state, control, logs=None, **kwargs):
        if logs is None:
            return
        
        loss = logs.get("loss", 0.0)
        recon_loss = logs.get("recon_loss", 0.0)
        kl_loss = logs.get("kl_loss", 0.0)
        grad_norm = logs.get("grad_norm", 0.0)
        learning_rate = logs.get("learning_rate", 0.0)
        epoch = logs.get("epoch", 0.0)

        formatted_log = (
            f"loss={loss:.4f} | recon_loss={recon_loss:.4f} | "
            f"kl_loss={kl_loss:.4f} | grad_norm={grad_norm:.4f} | "
            f"lr={learning_rate:.7f} | epoch={epoch:.4f} | global_step={state.global_step}"
        )

        self.logger.info(formatted_log)


class VanillaLoggingCallback(TrainerCallback):
    def __init__(self):
        self.logger = logger
    
    def on_log(self, args, state, control, logs=None, **kwargs):
        if logs is None:
            return

        loss = logs.get("loss", 0.0)
        grad_norm = logs.get("grad_norm", 0.0)
        learning_rate = logs.get("learning_rate", 0.0)
        epoch = logs.get("epoch", 0.0)

        formatted_log = (
            f"loss={loss:.4f} | "
            f"grad_norm={grad_norm:.4f} | lr={learning_rate:.7f} | "
            f"epoch={epoch:.4f} | global_step={state.global_step}"
            f"grad_norm={grad_norm:.4f} | lr={learning_rate:.7f} | "
            f"epoch={epoch:.4f} | global_step={state.global_step}"
        )

        self.logger.info(formatted_log)

import torch
import torch.nn.functional as F
import torchvision.transforms.functional as TF
import random
from typing import Sequence
import numpy as np

class RandomRotation:
    def __init__(self, angles: Sequence[int]):
        self.angles = angles

    def __call__(self, x):
        angle = random.choice(self.angles)
        return TF.rotate(x, angle)

def normalize(img, global_min, global_max):
    img = np.clip(img, global_min, global_max)
    return (img - global_min) / (global_max - global_min)


class CropOrPad:
    def __init__(
        self,
        target_height: int,
        target_width: int,
        crop_mode: str = "center",   # or 'random'
        pad_mode: str = "constant",  # or 'reflect', 'replicate', etc.
        pad_value: float = 0.0       # only used if pad_mode == 'constant'
    ):
        assert crop_mode in {"center", "random"}, "crop_mode must be 'center' or 'random'"
        self.target_height = target_height
        self.target_width = target_width
        self.crop_mode = crop_mode
        self.pad_mode = pad_mode
        self.pad_value = pad_value

    def __call__(self, x: torch.Tensor) -> torch.Tensor:
        """
        x: Tensor of shape (..., H, W)
        """
        original_shape = x.shape[:-2]  # e.g. (B, C)
        H, W = x.shape[-2:]
        th, tw = self.target_height, self.target_width

        x = x.reshape(-1, H, W)  # Flatten leading dims

        # -------------------- Pad -------------------- #
        pad_top = max((th - H) // 2, 0)
        pad_bottom = max(th - H - pad_top, 0)
        pad_left = max((tw - W) // 2, 0)
        pad_right = max(tw - W - pad_left, 0)

        if pad_top > 0 or pad_bottom > 0 or pad_left > 0 or pad_right > 0:
            padding = [pad_left, pad_right, pad_top, pad_bottom]  # (left, right, top, bottom)
            x = F.pad(x, padding, mode=self.pad_mode, value=self.pad_value)

        # -------------------- Crop -------------------- #
        H_new, W_new = x.shape[-2:]
        if H_new > th or W_new > tw:
            if self.crop_mode == "center":
                top = (H_new - th) // 2
                left = (W_new - tw) // 2
            else:  # random
                top = random.randint(0, H_new - th)
                left = random.randint(0, W_new - tw)
            x = TF.crop(x, top, left, th, tw)

        # Reshape back to original leading dimensions
        return x.reshape(*original_shape, th, tw)
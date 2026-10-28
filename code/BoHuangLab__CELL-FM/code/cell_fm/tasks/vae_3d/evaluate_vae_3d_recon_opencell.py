# -*- coding: utf-8 -*-
import os
import sys

import numpy as np
import torch
import tifffile

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([".", ".."])

from pathlib import Path

from cell_fm.data.opencell_3d_crop_data.dataset import OpenCell3DCropImageOnlyDataset
from cell_fm.models.vae_3d.vae_3d_config import VAE3DConfig
from cell_fm.models.vae_3d.vae_3d_model import VAE3DModel
from cell_fm.utils.cli_utils import cli


def to_uint16(vol: np.ndarray) -> np.ndarray:
    """Normalize a float volume to uint16 [0, 65535]."""
    vol = vol - vol.min()
    denom = vol.max()
    if denom > 1e-8:
        vol = vol / denom
    return (vol * 65535).clip(0, 65535).astype(np.uint16)


@cli(VAE3DConfig)
def main(args) -> None:
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    # gene_name = 'TUBB4B'  # Microtubules
    gene_name = 'AAMP'  # Microtubules
    valset = OpenCell3DCropImageOnlyDataset(args, split_key="all")
    # Filter to target gene using the flat img_paths / gene_name_labels lists
    mask = [g == gene_name for g in valset.gene_name_labels]
    valset.img_paths = [p for p, m in zip(valset.img_paths, mask) if m]
    valset.gene_name_labels = [gene_name] * len(valset.img_paths)

    model = VAE3DModel(config=VAE3DConfig(**vars(args)))
    model.to(device)
    model.eval()

    output_dir = Path(args.output_dir)

    with torch.no_grad():
        for i, data in enumerate(valset):
            print(f"Processing {i+1}/{len(valset)}")
            # protein_img: (1, D, H, W)
            protein_img = data['protein_img'].unsqueeze(0).to(device)  # (1, 1, D, H, W)

            latent_dist = model.encode(protein_img)
            latents = latent_dist.mean                                  # (1, C_lat, D_lat, H_lat, W_lat)

            print(f"Latent shape: {latents.shape}, min: {latents.min().item():.4f}, max: {latents.max().item():.4f}")

            protein_img_recon = model.decode(latent_dist.sample()).sample  # (1, 1, D, H, W)

            save_dir = output_dir / '{:04d}'.format(i + 1)
            save_dir.mkdir(parents=True, exist_ok=True)

            # Save full 3D volumes as TIF: (D, H, W) float32 → uint16
            orig_np  = protein_img[0, 0].cpu().numpy()       # (D, H, W)
            recon_np = protein_img_recon[0, 0].cpu().numpy()  # (D, H, W)

            tifffile.imwrite(save_dir / 'protein_img.tif',       to_uint16(orig_np))
            tifffile.imwrite(save_dir / 'protein_img_recon.tif', to_uint16(recon_np))

if __name__ == "__main__":
    main()

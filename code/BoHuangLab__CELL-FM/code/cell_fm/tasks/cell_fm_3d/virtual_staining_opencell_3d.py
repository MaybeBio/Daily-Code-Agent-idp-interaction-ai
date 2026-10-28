# -*- coding: utf-8 -*-
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([".", ".."])

import torch
import numpy as np
import tifffile as tiff
from pathlib import Path
from tqdm import tqdm

from cell_fm.data.opencell_3d_crop_data.dataset import OpenCell3DCropAllImageDataset
from cell_fm.models.cell_fm_3d.cell_fm_3d_config import CELLFM3DConfig
from cell_fm.models.cell_fm_3d.cell_fm_3d_model import CELLFM3DModel
from cell_fm.utils.cli_utils import cli


def save_tif(image, output_path):
    """Save tensor (C, D, H, W) uint16 TIFF. Expects values in [0, 1]."""
    tensor_np = image.cpu().numpy()
    tensor_np = tensor_np.clip(0, 1)
    tensor_np = np.round(tensor_np * 65535).astype(np.uint16)
    tensor_np = np.transpose(tensor_np, (1, 0, 2, 3))  # (C, D, H, W) -> (D, C, H, W)
    tiff.imwrite(output_path, tensor_np, imagej=True)


@cli(CELLFM3DConfig)
def main(args) -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"

    config = CELLFM3DConfig(**vars(args))
    valset = OpenCell3DCropAllImageDataset(config, split_key=config.split_key)
    model = CELLFM3DModel(config=config)

    model.to(device)
    model.eval()

    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Use one fixed nucleus image (ATG7, index 1) for all generations
    ref_gene = 'ATG7'
    ref_index = 2
    ref_valset = OpenCell3DCropAllImageDataset(config, split_key=config.split_key)
    ref_valset.meta_data = ref_valset.meta_data[
        ref_valset.meta_data['gene_name'] == ref_gene
    ].reset_index(drop=True)
    chosen_nucleus_img = ref_valset[0]['nucleus_imgs'][ref_index].unsqueeze(0).to(device)

    batch_size = 2

    for i, data in enumerate(tqdm(valset, desc="Generating 3D protein volumes")):
        protein_seq = data['protein_seq'].unsqueeze(0).to(device)

        if protein_seq.shape[1] > args.max_protein_sequence_len + 2:
            continue

        # num_samples = len(data['protein_imgs'])
        num_samples = 64

        save_dir = output_dir / f"{data['gene_name']}_{i + 1:04d}"
        save_dir.mkdir(parents=True, exist_ok=True)

        for batch_start in range(0, num_samples, batch_size):
            batch_end = min(batch_start + batch_size, num_samples)
            n = batch_end - batch_start
            protein_seq_batch = protein_seq.repeat(n, 1)

            # Generate 3D protein volume: (B, 1, D, H, W) in [-1, 1]
            nucleus_batch = chosen_nucleus_img.repeat(n, 1, 1, 1, 1)  # (B, 1, D, H, W)
            sample = model.sequence_to_image(protein_seq_batch, nucleus_batch, num_steps=args.num_steps)

            # Stack fixed nucleus + generated protein: (B, 2, D, H, W)
            cat_img = torch.cat([nucleus_batch, sample], dim=1)        # (B, 2, D, H, W)

            for j in range(n):
                # Convert from [-1,1] to [0,1]
                vol = (cat_img[j] + 1.0) / 2.0
                save_tif(vol, save_dir / f"{batch_start + j + 1:04d}.tif")


if __name__ == "__main__":
    main()

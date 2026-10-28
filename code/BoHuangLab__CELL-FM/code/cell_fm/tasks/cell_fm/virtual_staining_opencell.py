# -*- coding: utf-8 -*-
import os
import sys

import torch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([".", ".."])

from pathlib import Path

from cell_fm.criterions.cell_fm.unidiffuser import UniDiffCriterions
from cell_fm.data.opencell_crop_data.dataset import OpenCellCropAllImageDataset
from cell_fm.models.cell_fm.cell_fm_config import CELLFMConfig
from cell_fm.models.cell_fm.cell_fm_model import CELLFMModel
from cell_fm.utils.cli_utils import cli

import numpy as np
import tifffile as tiff
from tqdm import tqdm


def save_tif(image, output_path):
    tensor_np = image.cpu().numpy()
    tensor_np = tensor_np.clip(0, 1)
    tensor_np = np.round(tensor_np * 65535)

    tensor_np = tensor_np.astype(np.uint16)    
    tiff.imwrite(output_path, tensor_np, imagej=True)


@cli(CELLFMConfig)
def main(args) -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"

    config = CELLFMConfig(**vars(args))
    valset = OpenCellCropAllImageDataset(config, split_key=config.split_key)
    model = CELLFMModel(config=config, loss_fn=UniDiffCriterions)

    # chosen images
    gene_name = 'ATG7'
    index = 1
    valset.meta_data = valset.meta_data[valset.meta_data['gene_name'] == gene_name]
    chosen_data = valset.__getitem__(0)
    chosen_nucleus_img = chosen_data['nucleus_imgs'][index].unsqueeze(0).to(device)

    model.to(device)
    model.eval()

    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # batch size
    batch_size = 32
    valset = OpenCellCropAllImageDataset(config, split_key=config.split_key)

    for i, data in enumerate(tqdm(valset, desc="Processing valset")):
        protein_seq = data['protein_seq'].unsqueeze(0).to(device)

        if protein_seq.shape[1] > args.max_protein_sequence_len + 2:
            continue

        num_samples = len(data['nucleus_imgs'])

        for batch_start in range(0, num_samples, batch_size):
            batch_end = min(batch_start + batch_size, num_samples)
            # cell_imgs_batch = nucleus_imgs[batch_start:batch_end]
            cell_imgs_batch = chosen_nucleus_img.repeat(batch_end - batch_start, 1, 1, 1)
            protein_seq_batch = protein_seq.repeat(batch_end - batch_start, 1)

            # Generate images
            sample = model.sequence_to_image(
                protein_seq_batch, 
                cell_imgs_batch, 
                num_steps=args.num_steps, 
            )

            cat_img = torch.cat([cell_imgs_batch, sample], dim=1)

            save_file = output_dir / (data['gene_name']+f'_{i + 1:04d}')
            save_file.mkdir(parents=True, exist_ok=True)

            # Save sample images
            for j in range(batch_start, batch_end):
                # Save the generated image
                cat_img_j = cat_img[j - batch_start]
                # convert to [0, 1]
                cat_img_j = (cat_img_j + 1) / 2
                cat_img_j = cat_img_j.clamp(0, 1)

                save_tif(cat_img_j, save_file / ('{:04d}.tif'.format(j + 1)))


if __name__ == "__main__":
    main()
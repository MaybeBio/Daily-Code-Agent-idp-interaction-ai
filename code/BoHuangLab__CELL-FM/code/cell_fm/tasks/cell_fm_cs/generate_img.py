# -*- coding: utf-8 -*-
import os
import sys
import torch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([".", ".."])

from pathlib import Path

from cell_fm.criterions.cell_fm.unidiffuser import UniDiffCriterions
from cell_fm.data.condenseq_data.dataset import CondenSeqAllImageDataset
from cell_fm.models.cell_fm.cell_fm_cs_config import CELLFMCSConfig
from cell_fm.models.cell_fm.cell_fm_cs_model import CELLFMCSModel
from cell_fm.utils.cli_utils import cli
from cell_fm.logging import logger

from esm.utils import decoding
import tifffile as tiff
import numpy as np
import matplotlib.pyplot as plt


def save_tif(image, output_path):
    image = image * 65535

    image = image.astype(np.uint16)
    tiff.imwrite(output_path, image, imagej=True)

@cli(CELLFMCSConfig)
def main(args) -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"

    config = CELLFMCSConfig(**vars(args))
    model = CELLFMCSModel(config=config, loss_fn=UniDiffCriterions)

    model.to(device)
    model.eval()

    output_dir = Path(config.output_dir)
    output_dir = output_dir / config.split_key

    valset = CondenSeqAllImageDataset(args, split_key=args.split_key)
    vocab = valset.vocab

    batch_size = 32

    for i, data in enumerate(valset):
        protein_seq = data['protein_seq'].unsqueeze(0).to(device)

        if protein_seq.shape[1] > config.max_protein_sequence_len + 2:
            continue

        protein_seq = data['protein_seq'].unsqueeze(0).to(device)

        protein_imgs = torch.stack(data['protein_imgs']).to(device)
        nucleus_imgs = torch.stack(data['nucleus_imgs']).to(device)
        mask_imgs = torch.stack(data['mask_imgs']).to(device).round().long()

        protein_intensity_levels = torch.tensor(data['protein_intensity_levels']).float()
        protein_intensity_levels = protein_intensity_levels.unsqueeze(-1).to(device)

        logger.info(data['index'])

        save_file = output_dir / ('{:04d}_'.format(i+1) + str(data['index']))
        save_file.mkdir(parents=True, exist_ok=True)

        real_max_minus_medians = []
        real_intensity_levels = []

        protein_imgs = ((protein_imgs + 1) / 2).clamp(0, 1)

        for i in range(protein_imgs.shape[0]):
            masked_img = protein_imgs[i][mask_imgs[i] == 1]

            max_val = masked_img.max().item()
            median_val = masked_img.median().item()
            diff = max_val - median_val

            real_max_minus_medians.append(diff)
            real_intensity_levels.append(protein_intensity_levels[i].item())

        real_protein_seq = decoding.decode_sequence(protein_seq.squeeze(), vocab)
        logger.info(real_protein_seq)

        generated_imgs = []
        max_minus_medians_all = []

        for batch_start in range(0, protein_imgs.shape[0], batch_size):
            batch_end = min(batch_start + batch_size, protein_imgs.shape[0])

            chosen_cell_img_batch = nucleus_imgs[batch_start:batch_end]
            chosen_mask_img = mask_imgs[batch_start:batch_end]

            protein_intensity_level_batch = protein_intensity_levels[batch_start:batch_end]

            protein_seq_batch = protein_seq.repeat(protein_intensity_level_batch.shape[0], 1)

            sample = model.sequence_to_image(
                protein_seq_batch, 
                chosen_cell_img_batch, 
                protein_intensity_level_batch, 
                num_steps=args.num_steps, 
            )
            generated_imgs.append(sample.squeeze(1))

            sample = ((sample + 1) / 2).clamp(0, 1)
            mask = chosen_mask_img.round().long()

            for i in range(sample.shape[0]):
                masked_img = sample[i][mask[i] == 1]
                max_val = masked_img.max().item()
                median_val = masked_img.median().item()
                diff = max_val - median_val
                max_minus_medians_all.append(diff)

        # Concatenate all images
        generated_imgs = torch.cat(generated_imgs, dim=0)

        # rescale images to [0, 1]
        generated_imgs = (generated_imgs + 1) / 2
        generated_imgs = generated_imgs.clamp(0, 1)

        # save as tif file
        save_tif(generated_imgs.detach().cpu().numpy(), save_file / 'imgs.tif')

        # save real protein images
        protein_imgs = protein_imgs.clamp(0, 1)
        save_tif(protein_imgs.detach().cpu().numpy(), save_file / 'real_protein_imgs.tif')

        # Step 0: prepare figure and axis
        fig, ax = plt.subplots(figsize=(12, 6))

        # Step 1: plot scatter for generated data
        ax.scatter(protein_intensity_levels.squeeze().detach().cpu().numpy(), np.array(max_minus_medians_all), color='orange', alpha=0.6, s=15, label='Generated Data')

        # Step 2: plot scatter on same axis
        ax.scatter(real_intensity_levels, real_max_minus_medians, color='darkblue', alpha=0.6, s=15, label='Real Data')

        # Step 3: formatting
        ax.set_xlabel('Protein Intensity Level')
        ax.set_ylabel('Max - Median')
        ax.set_title('Max - Median vs Protein Intensity Level')
        ax.grid(axis='y')
        ax.legend()
        plt.tight_layout()
        plt.savefig(save_file / 'combined_max_minus_median_normalized_img.png')
        plt.close()

    logger.info("Done!")

if __name__ == "__main__":
    main()
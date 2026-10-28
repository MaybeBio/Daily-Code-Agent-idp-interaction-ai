# -*- coding: utf-8 -*-
import os
import sys
import torch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([".", ".."])

from pathlib import Path

from cell_fm.criterions.cell_fm.unidiffuser import UniDiffCriterions
from cell_fm.data.condenseq_data.dataset import CondenSeqAllImageDataset, CondenSeqDataset
from cell_fm.models.cell_fm.cell_fm_cs_config import CELLFMConfig
from cell_fm.models.cell_fm.cell_fm_cs_model import CELLFMCSModel
from cell_fm.utils.cli_utils import cli
from cell_fm.logging import logger

from esm.utils import decoding, encoding
import tifffile as tiff
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass


def save_tif(image, output_path):
    image = image * 65535

    image = image.astype(np.uint16)
    tiff.imwrite(output_path, image, imagej=True)

def mutate_sequence(protein_seq: str, mutation_type: str) -> str:
    if mutation_type == 'WT':
        mutated_sequence = protein_seq
    else:
        type1, type2 = mutation_type.split('2')
        mutated_sequence = protein_seq.replace(type1, type2)
    return mutated_sequence

@dataclass
class ExtendedConfig:
    mutation_type: str = 'WT'


@cli(CELLFMConfig, ExtendedConfig)
def main(args) -> None:

    device = "cuda" if torch.cuda.is_available() else "cpu"

    config = CELLFMConfig(**vars(args))
    model = CELLFMCSModel(config=config, loss_fn=UniDiffCriterions)

    model.to(device)
    model.eval()

    output_dir = Path(config.output_dir)
    output_dir = output_dir / config.split_key

    mutation_type = args.mutation_type
    output_dir = output_dir / mutation_type

    # chosen images
    protein_index = 12626
    index = 0

    valset = CondenSeqAllImageDataset(args, split_key='all')
    valset.meta_data = valset.meta_data[valset.meta_data['index'] == protein_index]
    chosen_data = valset.__getitem__(0)

    chosen_nucleus_img = chosen_data['nucleus_imgs'][index].unsqueeze(0).to(device)
    chosen_mask_img = chosen_data['mask_imgs'][index].unsqueeze(0).to(device).round().long()

    chosen_cell_img = chosen_nucleus_img

    batch_size = 512
    num_protein_intensity_levels = 512

    min_protien_intensity_level = 20 # to avoid log(0) = -inf
    max_protein_intensity_level = 2048

    valset = CondenSeqDataset(args, split_key=args.split_key)
    vocab = valset.vocab

    for i, data in enumerate(valset):

        from cell_fm.pipeline.accelerator.trainer import seed_everything
        seed_everything(args.seed)

        protein_seq = data['protein_seq'].unsqueeze(0).to(device)

        if protein_seq.shape[1] > config.max_protein_sequence_len + 2:
            continue

        protein_seq_txt = decoding.decode_sequence(protein_seq.squeeze(), vocab)
        mutated_protein_seq_txt = mutate_sequence(protein_seq_txt, mutation_type)

        protein_seq = encoding.tokenize_sequence(mutated_protein_seq_txt, vocab, True)
        protein_seq = protein_seq.unsqueeze(0).to(device)

        logger.info(data['index'])

        save_file = output_dir / ('{:04d}_'.format(i+1) + str(data['index']))
        save_file.mkdir(parents=True, exist_ok=True)

        real_protein_seq = decoding.decode_sequence(protein_seq.squeeze(), vocab)
        logger.info(real_protein_seq)

        generated_imgs = []
        max_minus_medians_all = []

        # create logspace of protein intensity levels from min_protien_intensity_level to max_protein_intensity_level
        protein_intensity_levels_logspace = torch.logspace(
            start=np.log10(min_protien_intensity_level),
            end=np.log10(max_protein_intensity_level),
            steps=num_protein_intensity_levels,
        ).to(device)

        for batch_start in range(0, num_protein_intensity_levels, batch_size):
            batch_end = min(batch_start + batch_size, num_protein_intensity_levels)

            protein_intensity_level_batch = protein_intensity_levels_logspace[batch_start:batch_end].unsqueeze(-1)

            protein_seq_batch = protein_seq.repeat(protein_intensity_level_batch.shape[0], 1)
            chosen_cell_img_batch = chosen_cell_img.repeat(protein_intensity_level_batch.shape[0], 1, 1, 1)

            sample = model.sequence_to_image(
                protein_seq_batch, 
                chosen_cell_img_batch, 
                protein_intensity_level_batch, 
                num_steps=config.num_steps, 
            )
            generated_imgs.append(sample.squeeze(1))

            sample = ((sample + 1) / 2).clamp(0, 1).squeeze(1)
            mask = chosen_mask_img.squeeze().round().long()

            masked_samples = sample[:, mask == 1] # shape (B, N), where N = number of masked voxels

            max_vals = masked_samples.max(dim=1).values
            median_vals = masked_samples.median(dim=1).values

            max_minus_medians = (max_vals - median_vals).tolist()
            max_minus_medians_all.append(max_minus_medians)

        # Concatenate all images
        generated_imgs = torch.cat(generated_imgs, dim=0)

        # rescale images to [0, 1]
        generated_imgs = (generated_imgs + 1) / 2
        generated_imgs = generated_imgs.clamp(0, 1)

        # save as tif file
        save_tif(generated_imgs.detach().cpu().numpy(), save_file / 'imgs.tif')
        # save protein_intensity_levels_logspace as npy
        np.save(save_file / 'protein_intensity_levels.npy', protein_intensity_levels_logspace.detach().cpu().numpy())
        # save real_protein_seq as txt
        with open(save_file / 'protein_seq.txt', 'w') as f:
            f.write(real_protein_seq)

        # Step 0: prepare figure and axis
        fig, ax = plt.subplots(figsize=(12, 6))

        # Step 1: plot generated data scatter
        x_vals = protein_intensity_levels_logspace.detach().cpu().tolist()
        y_vals = [v for batch in max_minus_medians_all for v in batch]
        ax.scatter(x_vals, y_vals, color='orange', alpha=0.6, s=15, label='Generated Data')

        # Step 2: formatting
        ax.set_xscale('log')
        ax.set_xlim(min_protien_intensity_level, max_protein_intensity_level)
        ax.set_xlabel('Protein Intensity Level (log scale)')
        ax.set_ylabel('Max - Median')
        ax.set_title('Max - Median vs Protein Intensity Level (log scale)')
        ax.grid(axis='y')
        ax.legend()
        plt.tight_layout()
        plt.savefig(save_file / 'Max_minus_median_normalized_img.png')
        plt.close()

    logger.info("Done!")

if __name__ == "__main__":
    main()
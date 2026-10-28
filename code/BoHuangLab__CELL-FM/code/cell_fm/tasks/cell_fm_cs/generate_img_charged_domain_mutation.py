# -*- coding: utf-8 -*-
import os
import sys
import torch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([".", ".."])

from pathlib import Path

from cell_fm.criterions.cell_fm.unidiffuser import UniDiffCriterions
from cell_fm.data.condenseq_data.dataset import CondenSeqAllImageDataset
from cell_fm.models.cell_fm.cell_fm_cs_config import CELLFMConfig
from cell_fm.models.cell_fm.cell_fm_cs_model import CELLFMCSModel
from cell_fm.utils.cli_utils import cli
from cell_fm.logging import logger

import tifffile as tiff
import numpy as np

from tqdm import tqdm
from esm.utils import encoding
import pandas as pd


def save_tif(image, output_path):
    # image = np.transpose(image, (1, 0, 2, 3))
    image = image * 65535

    image = image.astype(np.uint16)
    tiff.imwrite(output_path, image, imagej=True)


@cli(CELLFMConfig)
def main(args) -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"

    config = CELLFMConfig(**vars(args))
    model = CELLFMCSModel(config=config, loss_fn=UniDiffCriterions)

    model.to(device)
    model.eval()

    meta_data = pd.read_csv("cell_fm/tasks/cell_fm_cs/data/random_KRDE/panel.csv")

    all_seq_dir = {}

    for index, row in tqdm(meta_data.iterrows(), total=len(meta_data)):
        protein_index = row['index']
        protein_seq = row['protein_seq']

        all_seq_dir[f'{index+1:04d}_{protein_index}'] = protein_seq

    # chosen images
    protein_index = 12626
    index = 0

    valset = CondenSeqAllImageDataset(args, split_key='all')
    valset.meta_data = valset.meta_data[valset.meta_data['index'] == protein_index]
    chosen_data = valset.__getitem__(0)

    vocab = valset.vocab

    chosen_nucleus_img = chosen_data['nucleus_imgs'][index].unsqueeze(0).to(device)

    chosen_cell_img = chosen_nucleus_img

    min_protein_intensity_level = 20  # to avoid log(0) = -inf
    max_protein_intensity_level = 2048
    num_protein_intensity_levels = 512
    batch_size = 512

    protein_intensity_levels_logspace = torch.logspace(
        start=np.log10(min_protein_intensity_level),
        end=np.log10(max_protein_intensity_level),
        steps=num_protein_intensity_levels,
    ).to(device)

    for key, protein_seq in tqdm(all_seq_dir.items()):
        output_folder_name = key

        sub_seq_dir = {}

        sub_seq_dir["WT"] = protein_seq

        from generate_mutation_kappa import make_alternating_pattern, make_blocky_pattern, make_mutant_from_pattern

        # blocky (highly segregated)
        blocky_pattern = make_blocky_pattern(protein_seq)
        seq_blocky = make_mutant_from_pattern(protein_seq, blocky_pattern)

        sub_seq_dir[f"blocky"] = seq_blocky

        # alternating (well-mixed)
        alt_pattern = make_alternating_pattern(protein_seq)
        seq_alt = make_mutant_from_pattern(protein_seq, alt_pattern)

        sub_seq_dir[f"alternating"] = seq_alt

        for mut_key, seq in sub_seq_dir.items():

            output_dir = Path(config.output_dir)
            output_dir = output_dir / output_folder_name / mut_key
            output_dir.mkdir(parents=True, exist_ok=True)

            protein_seq = encoding.tokenize_sequence(seq, vocab, True)
            protein_seq = protein_seq.unsqueeze(0).to(device)

            generated_imgs = []

            for batch_start in range(0, num_protein_intensity_levels, batch_size):
                batch_end = min(batch_start + batch_size, num_protein_intensity_levels)
                current_batch_size = batch_end - batch_start

                protein_intensity_level_batch = protein_intensity_levels_logspace[batch_start:batch_end]
                protein_intensity_level_batch = protein_intensity_level_batch.unsqueeze(1).to(device)
                protein_seq_batch = protein_seq.repeat(current_batch_size, 1).to(device)
                chosen_cell_img_batch = chosen_cell_img.repeat(current_batch_size, 1, 1, 1).to(device)

                # Generate images
                sample = model.sequence_to_image(
                    protein_seq_batch, 
                    chosen_cell_img_batch, 
                    protein_intensity_level_batch, 
                    num_steps=args.num_steps, 
                )
                generated_imgs.append(sample)

            generated_imgs = torch.cat(generated_imgs, dim=0)

            # rescale images to [0, 1]
            combined_generated_imgs = torch.cat([chosen_nucleus_img.repeat(generated_imgs.size(0), 1, 1, 1), generated_imgs], dim=1)  # [N, 2, H, W]
            combined_generated_imgs = (combined_generated_imgs + 1) / 2.0  # shift to [0,1]
            combined_generated_imgs = combined_generated_imgs.clamp(0.0, 1.0)

            save_tif(
                combined_generated_imgs.cpu().numpy(),
                output_dir / "combined_images.tif", 
            )

            # save protein_intensity_levels_logspace as npy
            np.save(
                output_dir / "protein_intensity_levels.npy",
                protein_intensity_levels_logspace.cpu().numpy(),
            )

    logger.info("Done!")

if __name__ == "__main__":
    main()
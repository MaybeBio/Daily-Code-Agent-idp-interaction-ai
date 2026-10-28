# -*- coding: utf-8 -*-
import os
import sys
import torch

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.extend([".", ".."])

from pathlib import Path

from cell_fm.criterions.cell_fm.unidiffuser import UniDiffCriterions
from cell_fm.data.hpa_data.dataset import HPAAllImageDataset
from cell_fm.models.cell_fm.cell_fm_config import CELLFMConfig
from cell_fm.models.cell_fm.cell_fm_model import CELLFMModel
from cell_fm.utils.cli_utils import cli
from cell_fm.logging import logger

from torchvision.utils import save_image
from esm.utils import decoding

def colorize_image(tensor, color):
    # Create a zero tensor with the same size as the input tensor but with three channels
    colored_image = torch.zeros((3, tensor.size(0), tensor.size(1)), dtype=tensor.dtype, device=tensor.device)

    if color == 'blue':
        colored_image[2] = tensor  # Set blue channel
    elif color == 'red':
        colored_image[0] = tensor  # Set red channel
    elif color == 'green':
        colored_image[1] = tensor  # Set green channel
    elif color == 'yellow':
        colored_image[0] = tensor  # Set red and green channels to get yellow
        colored_image[1] = tensor
    return colored_image

def save_colored_image(tensor, filename, color):
    colored_tensor = colorize_image(tensor, color)
    save_image(colored_tensor, filename, normalize=True, value_range=(0, 1))


@cli(CELLFMConfig)
def main(args) -> None:
    if not torch.cuda.is_available():
        device = "cpu"
    else:
        device = "cuda"

    config = CELLFMConfig(**vars(args))

    valset = HPAAllImageDataset(args, split_key=config.split_key)
    model = CELLFMModel(config=config, loss_fn=UniDiffCriterions)

    vocab = valset.vocab

    model.to(device)
    model.eval()

    output_dir = Path(config.output_dir)
    output_dir = output_dir / config.split_key

    batch_size = 2

    for i, data in enumerate(valset):
        protein_seq = data['protein_seq'].unsqueeze(0).to(device)

        if protein_seq.shape[1] > config.max_protein_sequence_len + 2:
            continue

        protein_seq = data['protein_seq'].unsqueeze(0).to(device)

        protein_imgs = torch.stack(data['protein_imgs']).to(device)
        nucleus_imgs = torch.stack(data['nucleus_imgs']).to(device)
        microtubules_imgs = torch.stack(data['microtubules_imgs']).to(device)
        ER_imgs = torch.stack(data['ER_imgs']).to(device)

        cell_imgs = torch.cat([nucleus_imgs, ER_imgs, microtubules_imgs], dim=1)
        
        logger.info(data['gene_name'])

        save_file = output_dir / ('{:04d}_'.format(i+1) + data['gene_name'])
        save_file.mkdir(parents=True, exist_ok=True)

        real_protein_seq = decoding.decode_sequence(protein_seq.squeeze(), vocab)
        logger.info(real_protein_seq)

        protein_seq = protein_seq.repeat(cell_imgs.shape[0], 1)  # (B, L)

        for batch_start in range(0, cell_imgs.shape[0], batch_size):
            batch_end = min(batch_start + batch_size, cell_imgs.shape[0])
            batch_cell_imgs = cell_imgs[batch_start:batch_end].to(device)
            batch_protein_seq = protein_seq[batch_start:batch_end].to(device)
            batch_prot_imgs = protein_imgs[batch_start:batch_end].to(device)

            # Generate images
            sample = model.sequence_to_image(
                batch_protein_seq, 
                batch_cell_imgs, 
                num_steps=config.num_steps, 
            )

            generated_save_dir = save_file / 'generated'
            real_save_dir = save_file / 'real'

            generated_save_dir.mkdir(parents=True, exist_ok=True)
            real_save_dir.mkdir(parents=True, exist_ok=True)

            # save images
            for j in range(batch_start, batch_end):
                # save samples
                save_image(sample[j - batch_start], generated_save_dir / f'img_{j}.png', normalize=True, value_range=(-1, 1))
                # save protein images
                save_image(batch_prot_imgs[j - batch_start], real_save_dir / f'img_{j}.png', normalize=True, value_range=(-1, 1))

    logger.info("Done!")

if __name__ == "__main__":
    main()

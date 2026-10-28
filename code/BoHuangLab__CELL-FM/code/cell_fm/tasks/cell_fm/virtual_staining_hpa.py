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

from torchvision.utils import save_image
from esm.utils import encoding, decoding


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
    valset = HPAAllImageDataset(config, split_key=config.split_key)
    model = CELLFMModel(config=config, loss_fn=UniDiffCriterions)

    vocab = valset.vocab

    model.to(device)
    model.eval()

    output_dir = Path(config.output_dir)
    output_dir = output_dir / config.split_key

    chosen_data = valset.__getitem__(1)
    index = 0

    chosen_nucleus_img = chosen_data['nucleus_imgs'][index].unsqueeze(0).to(device)
    chosen_microtubules_img = chosen_data['microtubules_imgs'][index].unsqueeze(0).to(device)
    chosen_ER_img = chosen_data['ER_imgs'][index].unsqueeze(0).to(device)

    num_samples = 5

    with torch.no_grad():
        if config.cell_image == 'nucl':
            chosen_cell_img = chosen_nucleus_img
        elif config.cell_image == 'nucl,er':
            if config.test_cell_image == 'nucl':
                chosen_ER_img = torch.full_like(chosen_nucleus_img, fill_value=-2)
            chosen_cell_img = torch.cat([chosen_nucleus_img, chosen_ER_img], dim=1)
        elif config.cell_image == 'nucl,mt':
            if config.test_cell_image == 'nucl':
                chosen_microtubules_img = torch.full_like(chosen_nucleus_img, fill_value=-2)
            chosen_cell_img = torch.cat([chosen_nucleus_img, chosen_microtubules_img], dim=1)
        elif config.cell_image == 'nucl,er,mt':
            if config.test_cell_image == 'nucl':
                chosen_ER_img = torch.full_like(chosen_nucleus_img, fill_value=-2)
                chosen_microtubules_img = torch.full_like(chosen_nucleus_img, fill_value=-2)
            elif config.test_cell_image == 'nucl,er':
                chosen_microtubules_img = torch.full_like(chosen_nucleus_img, fill_value=-2)
            elif config.test_cell_image == 'nucl,mt':
                chosen_ER_img = torch.full_like(chosen_nucleus_img, fill_value=-2)
            chosen_cell_img = torch.cat([chosen_nucleus_img, chosen_ER_img, chosen_microtubules_img], dim=1)
        else:
            raise ValueError(f"Cell image type: {config.cell_image} is not supported")

    for i, data in enumerate(valset):
        protein_seq = data['protein_seq'].unsqueeze(0).to(device)

        if protein_seq.shape[1] > args.max_protein_sequence_len:
            continue

        protein_img = data['protein_imgs'][index].unsqueeze(0).to(device)
        nucleus_img = data['nucleus_imgs'][index].unsqueeze(0).to(device)
        microtubules_img = data['microtubules_imgs'][index].unsqueeze(0).to(device)
        ER_img = data['ER_imgs'][index].unsqueeze(0).to(device)

        print(data['gene_name'])

        save_file = output_dir / ('{:04d}_'.format(i+1) + data['gene_name'])
        save_file.mkdir(parents=True, exist_ok=True)

        real_protein_seq = decoding.decode_sequence(protein_seq.squeeze(), vocab)
        print(real_protein_seq)
        
        real_img = torch.cat([torch.full_like(protein_img, -1), protein_img, nucleus_img], dim=1)
        save_image(real_img, save_file / 'real_img.png', normalize=True, value_range=(-1, 1))
        save_image(nucleus_img, save_file / 'real_nucleus_img.png', normalize=True, value_range=(-1, 1))
        save_image(protein_img, save_file / 'real_protein_img.png', normalize=True, value_range=(-1, 1))

        for j in range(num_samples):
            sample = model.sequence_to_image(
                protein_seq, 
                chosen_cell_img, 
                num_steps=config.num_steps, 
            )
            pred_img = torch.cat([torch.full_like(nucleus_img, -1), sample, chosen_nucleus_img], dim=1)

            save_image(pred_img, save_file / f'generated_img_{j}.png', normalize=True, value_range=(-1, 1))
            save_image(sample, save_file / f'generated_protein_img_{j}.png', normalize=True, value_range=(-1, 1))
            save_colored_image((sample.squeeze() + 1) / 2, save_file / f'generated_protein_img_green_{j}.png', 'green')

        save_image(chosen_nucleus_img, save_file / 'chosen_real_nucleus_img.png', normalize=True, value_range=(-1, 1))

        save_colored_image((chosen_nucleus_img.squeeze() + 1) / 2, save_file / 'chosen_nucleus_img_blue.png', 'blue')
        save_colored_image((chosen_microtubules_img.squeeze() + 1) / 2, save_file / 'chosen_microtubules_img_red.png', 'red')
        save_colored_image((chosen_ER_img.squeeze() + 1) / 2, save_file / 'chosen_ER_img_yellow.png', 'yellow')

        save_colored_image((nucleus_img.squeeze() + 1) / 2, save_file / 'real_nucleus_img_blue.png', 'blue')
        save_colored_image((microtubules_img.squeeze() + 1) / 2, save_file / 'real_microtubules_img_red.png', 'red')
        save_colored_image((protein_img.squeeze() + 1) / 2, save_file / 'real_protein_img_green.png', 'green')
        save_colored_image((ER_img.squeeze() + 1) / 2, save_file / 'real_ER_img_yellow.png', 'yellow')


if __name__ == "__main__":
    main()




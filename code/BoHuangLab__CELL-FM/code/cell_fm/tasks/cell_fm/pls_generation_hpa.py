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
from dataclasses import dataclass, field

from torchvision.utils import save_image
import csv
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


@dataclass
class Args:
    pls_type : str = 'nls'  # 'nls', 'nes'
    num_aas: list = field(default_factory=lambda: [10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25])
    num_gen: int = 20
    gen_pos: str = 'end'  # 'start', 'end'


@cli(CELLFMConfig, Args)
def main(args) -> None:
    device = 'cuda' if torch.cuda.is_available() else 'cpu'

    config = CELLFMConfig(**vars(args))

    valset = HPAAllImageDataset(args, split_key=args.split_key)
    model = CELLFMModel(config=config, loss_fn=UniDiffCriterions)

    vocab = valset.vocab

    model.to(device)
    model.eval()

    output_dir = Path(config.output_dir)
    
    if config.pls_type == 'nls':
        chosen_data = HPAAllImageDataset(args, split_key='all').__getitem__(8259)
        index = 10
    elif config.pls_type == 'nes':
        chosen_data = HPAAllImageDataset(args, split_key='test').__getitem__(24)
        index = 0

    output_dir = Path(output_dir)

    chosen_nucleus_img = chosen_data['nucleus_imgs'][index].unsqueeze(0).to(device)
    chosen_ER_img = chosen_data['ER_imgs'][index].unsqueeze(0).to(device)
    chosen_microtubules_img = chosen_data['microtubules_imgs'][index].unsqueeze(0).to(device)
    chosen_protein_img = chosen_data['protein_imgs'][index].unsqueeze(0).to(device)

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
    
    output_dir.mkdir(parents=True, exist_ok=True)

    save_image(chosen_nucleus_img, output_dir / 'chosen_nucleus_img.png', normalize=True, value_range=(-1, 1))
    save_image(chosen_protein_img, output_dir / 'chosen_protein_img.png', normalize=True, value_range=(-1, 1))
    save_image(chosen_ER_img, output_dir / 'chosen_ER_img.png', normalize=True, value_range=(-1, 1))
    save_image(chosen_microtubules_img, output_dir / 'chosen_microtubules_img.png', normalize=True, value_range=(-1, 1))

    cat_img = torch.cat([torch.full_like(chosen_nucleus_img, -1), chosen_protein_img, chosen_nucleus_img], dim=1)
    save_image(cat_img, output_dir / 'cat_img.png', normalize=True, value_range=(-1, 1))

    num_aas = args.num_aas

    num_gen = args.num_gen
    num_val = 0
    gen_pos = args.gen_pos

    gen_signals = []

    # Prompt sequence for generation
    # seq = 'MPRQGSLGAAPPKVAPDSSETVVG'
    seq = 'MPSQGSLGAAPPEVAPDSSETEEG'

    for num_aa in num_aas:
        if gen_pos == 'start':
            masked_seq = "M" + "<mask>" * num_aa + seq[1:]
        elif gen_pos == 'end':
            masked_seq = seq + "<mask>" * num_aa

        print(masked_seq)

        protein_seq = encoding.tokenize_sequence(masked_seq, vocab, True)
        protein_seq = protein_seq.unsqueeze(0).to(device)
        protein_seq_mask = (protein_seq == vocab.mask_token_id)

        for i in range(num_gen):
            sample = model.image_to_sequence(
                protein_seq, protein_seq_mask, chosen_protein_img, 
                chosen_cell_img, order='random', temperature=1.0, 
                progress=False, 
            )
            gen_seq = decoding.decode_sequence(sample.squeeze(), vocab)

            if gen_pos == 'start':
                gen_sig = gen_seq[1:num_aa+1]
            elif gen_pos == 'end':
                gen_sig = gen_seq[-num_aa:]

            if gen_sig not in gen_signals:
                gen_signals.append(gen_sig)

                seq_save = [gen_seq, gen_sig]
                save_dir = output_dir / '{:04d}_{:02d}'.format(num_aa, i+1)
                save_dir.mkdir(parents=True, exist_ok=True)

                with open(save_dir / 'seq.csv', 'w', newline='') as file:
                    writer = csv.writer(file)
                    for item in seq_save:
                        writer.writerow([item])

                for j in range(num_val):
                    protein_image_gen = model.sequence_to_image(
                        sample, chosen_cell_img, config.num_steps, 
                    )

                    save_image(
                        protein_image_gen, 
                        save_dir / 'protein_image_gen_{:02d}.png'.format(j+1), 
                        normalize=True, 
                        value_range=(-1, 1)
                    )
    # save gen_signals
    with open(output_dir / 'gen_signals.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        for signal in gen_signals:
            writer.writerow([signal])

if __name__ == "__main__":
    main()
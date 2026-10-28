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
from esm.utils import encoding
from cell_fm.pipeline.accelerator.trainer import seed_everything


def colorize_image(tensor, color):
    colored_image = torch.zeros((3, tensor.size(0), tensor.size(1)), dtype=tensor.dtype, device=tensor.device)
    if color == 'blue':
        colored_image[2] = tensor
    elif color == 'red':
        colored_image[0] = tensor
    elif color == 'green':
        colored_image[1] = tensor
    elif color == 'yellow':
        colored_image[0] = tensor
        colored_image[1] = tensor
    return colored_image

def save_colored_image(tensor, filename, color):
    colored_tensor = colorize_image(tensor, color)
    save_image(colored_tensor, filename, normalize=True, value_range=(0, 1))


@cli(CELLFMConfig)
def main(args) -> None:
    device = "cuda" if torch.cuda.is_available() else "cpu"

    config = CELLFMConfig(**vars(args))
    dataset = HPAAllImageDataset(config, split_key='all')
    vocab = dataset.vocab

    # ── anchor cell image ─────────────────────────────────────────────────────
    dataset.meta_data = dataset.meta_data[dataset.meta_data['gene_name'] == 'H3C13']
    chosen_data = dataset.__getitem__(0)
    index = 33

    model = CELLFMModel(config=config, loss_fn=UniDiffCriterions)
    model.to(device)
    model.eval()

    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    chosen_nucleus_img      = chosen_data['nucleus_imgs'][index].unsqueeze(0).to(device)
    chosen_microtubules_img = chosen_data['microtubules_imgs'][index].unsqueeze(0).to(device)
    chosen_ER_img           = chosen_data['ER_imgs'][index].unsqueeze(0).to(device)

    with torch.no_grad():
        if config.cell_image == 'nucl':
            chosen_cell_img = chosen_nucleus_img
        elif config.cell_image == 'nucl,er':
            if config.test_cell_image == 'nucl':
                chosen_ER_img = torch.full_like(chosen_nucleus_img, fill_value=-2.0)
            chosen_cell_img = torch.cat([chosen_nucleus_img, chosen_ER_img], dim=1)
        elif config.cell_image == 'nucl,mt':
            if config.test_cell_image == 'nucl':
                chosen_microtubules_img = torch.full_like(chosen_nucleus_img, fill_value=-2.0)
            chosen_cell_img = torch.cat([chosen_nucleus_img, chosen_microtubules_img], dim=1)
        elif config.cell_image == 'nucl,er,mt':
            if config.test_cell_image == 'nucl':
                chosen_ER_img           = torch.full_like(chosen_nucleus_img, fill_value=-2.0)
                chosen_microtubules_img = torch.full_like(chosen_nucleus_img, fill_value=-2.0)
            elif config.test_cell_image == 'nucl,er':
                chosen_microtubules_img = torch.full_like(chosen_nucleus_img, fill_value=-2.0)
            elif config.test_cell_image == 'nucl,mt':
                chosen_ER_img = torch.full_like(chosen_nucleus_img, fill_value=-2.0)
            chosen_cell_img = torch.cat([chosen_nucleus_img, chosen_ER_img, chosen_microtubules_img], dim=1)
        else:
            raise ValueError(f"Unsupported cell_image: {config.cell_image}")

    num_samples = 64
    batch_size  = 64

    # ── sequence dict: sliding-window section on PRRSV ──────────────────────
    wt_name = 'PRRSV'
    wt_seq  = 'MPNNNGKQQKRKKGDGQPVNQLCQMLGKIIAQQNQSRGKGPGKKNKKKNPEKPHFPLATEDDVRHHFTPSERQLCLSSIQTAFNQGAGTCTLSDSGRISYTVEFSLPTHHTVRLIRVTASPSA'
    window  = 25

    seq_dict = {}
    for i in range(1, len(wt_seq) - window + 1):
        mut_seq  = "M" + wt_seq[i:i+window]
        mut_name = f"{i+1}-{i+window}"
        seq_dict[mut_name] = mut_seq

    # ── save anchor cell images once ─────────────────────────────────────────
    anchor_dir = output_dir / wt_name
    anchor_dir.mkdir(parents=True, exist_ok=True)
    save_colored_image((chosen_nucleus_img.squeeze()      + 1) / 2, anchor_dir / 'chosen_nucleus.png',      color='blue')
    save_colored_image((chosen_microtubules_img.squeeze() + 1) / 2, anchor_dir / 'chosen_microtubules.png', color='red')
    save_colored_image((chosen_ER_img.squeeze()           + 1) / 2, anchor_dir / 'chosen_ER.png',           color='green')

    # ── generate ──────────────────────────────────────────────────────────────
    for prot_name, prot_seq in seq_dict.items():
        print(f"prot_seq: {prot_seq}")
        
        protein_seq = encoding.tokenize_sequence(prot_seq, vocab, True).unsqueeze(0).to(device)

        seed_everything(args.seed) # Fix random seed for each protein sequence to reduce randomness in generation.

        save_file = output_dir / wt_name / prot_name
        save_file.mkdir(parents=True, exist_ok=True)

        print(f"[{prot_name}]  len={len(prot_seq)}  tokens={protein_seq.shape[1]}")

        for batch_start in range(0, num_samples, batch_size):
            batch_end          = min(batch_start + batch_size, num_samples)
            current_batch_size = batch_end - batch_start

            samples = model.sequence_to_image(
                protein_seq.repeat(current_batch_size, 1),
                chosen_cell_img.repeat(current_batch_size, 1, 1, 1),
                num_steps=config.num_steps,
            )
            for j in range(current_batch_size):
                save_image(samples[j], save_file / f'generated_img_{batch_start + j}.png',
                           normalize=True, value_range=(-1, 1))


if __name__ == "__main__":
    main()

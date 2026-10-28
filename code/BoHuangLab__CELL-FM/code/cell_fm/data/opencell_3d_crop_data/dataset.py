# -*- coding: utf-8 -*-
import os
import random
from typing import List

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset

import tifffile as tiff
from esm.utils import encoding
from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer

from cell_fm.data.utils.seq_utils import OAARDM_sequence_masking
from cell_fm.data.opencell_crop_data.locations import all_locations
from .collater import collate_fn


def _normalize_volume(vol: np.ndarray, p_lo: float = 0.1, p_hi: float = 99.9) -> np.ndarray:
    """Percentile clip (p_lo, p_hi) then min-max normalize to [0, 1]."""
    lo, hi = np.percentile(vol, p_lo), np.percentile(vol, p_hi)
    vol = np.clip(vol, lo, hi)
    denom = hi - lo
    if denom < 1e-8:
        return np.zeros_like(vol)
    return (vol - lo) / denom


class OpenCell3DCropDataset(Dataset):
    """Dataset for 3D OpenCell volumetric crops.

    Each TIF has shape (48, 2, 192, 192): 48 z-slices, channel 0 = nucleus, channel 1 = protein.
    Returns protein_img and nucleus_img as (1, 48, 192, 192) tensors in [-1, 1].
    """

    def __init__(self, args, split_key: str) -> None:
        super().__init__()
        self.args = args
        self.split_key = split_key
        self.data_path = args.data_path

        self.vocab = EsmSequenceTokenizer()
        # Single metadata CSV covers all proteins (no per-split files)
        self.meta_data = pd.read_csv(
            os.path.join(self.data_path, f'{split_key}_3d_192_crop_merged_meta_data.csv')
        )

        self.data_aug = args.data_aug
        self.seq_zero_mask_ratio = args.seq_zero_mask_ratio
        self.max_protein_sequence_len = args.max_protein_sequence_len
        self.pre_pad_seq = getattr(args, 'pre_pad_seq', False)
        self.seq_length_control = getattr(args, 'seq_length_control', 'filter')
        self.phase = args.phase

        self.gene_names = sorted(self.meta_data['gene_name'].unique().tolist())

    def get_img(self, image_path: str):
        """Load a 3D TIF and return (protein_img, nucleus_img), each (1, D, H, W) in [-1, 1]."""
        img = tiff.imread(os.path.join(self.data_path, image_path)).astype(np.float32)
        # img: (D, 2, H, W)  where D=48, H=W=192
        nucleus_vol = img[:, 0, :, :]   # (D, H, W)
        protein_vol = img[:, 1, :, :]   # (D, H, W)

        nucleus_vol = _normalize_volume(nucleus_vol)
        protein_vol = _normalize_volume(protein_vol)

        # (D, H, W) → (1, D, H, W)
        nucleus_img = torch.from_numpy(nucleus_vol).unsqueeze(0)
        protein_img = torch.from_numpy(protein_vol).unsqueeze(0)

        if self.data_aug:
            stacked = torch.cat([protein_img, nucleus_img], dim=0)  # (2, D, H, W)
            if random.random() < 0.5:
                stacked = stacked.flip(-1)   # horizontal flip
            if random.random() < 0.5:
                stacked = stacked.flip(-2)   # vertical flip
            if random.random() < 0.5:
                stacked = stacked.flip(1)    # z-axis flip
            protein_img, nucleus_img = stacked[0:1], stacked[1:2]

        # [0, 1] → [-1, 1]
        protein_img = protein_img * 2.0 - 1.0
        nucleus_img = nucleus_img * 2.0 - 1.0

        return protein_img, nucleus_img

    def __getitem__(self, index: int) -> dict:
        # In 'crop' modes the collater shortens over-long sequences instead of dropping them,
        # so resampling here would keep them out of training entirely.
        if self.phase == 'train' and self.seq_length_control == 'filter':
            while True:
                meta_data = self.meta_data.iloc[index]
                if len(meta_data['sequence']) <= self.max_protein_sequence_len:
                    break
                index = random.randint(0, len(self) - 1)
        else:
            meta_data = self.meta_data.iloc[index]

        item = {}

        image_path = random.choice(meta_data['image_paths'].split(',')).strip()
        item['protein_img'], item['nucleus_img'] = self.get_img(image_path)

        protein_seq = meta_data['sequence']
        protein_seq_masked, protein_seq_mask, zm_label = OAARDM_sequence_masking(protein_seq, self.seq_zero_mask_ratio)

        if self.pre_pad_seq:
            sequence_len = len(protein_seq)
            protein_seq = '<cls>' + protein_seq + '<eos>' + '<pad>' * (self.max_protein_sequence_len - sequence_len)
            protein_seq_masked = '<cls>' + protein_seq_masked + '<mask>' + '<mask>' * (self.max_protein_sequence_len - sequence_len)
            protein_seq_token = encoding.tokenize_sequence(protein_seq, self.vocab, False)
            protein_seq_masked_token = encoding.tokenize_sequence(protein_seq_masked, self.vocab, False)
            protein_seq_mask = np.pad(protein_seq_mask, (1, 1 + len(protein_seq_token) - 2 - len(protein_seq_mask)), mode='constant', constant_values=(False, True))
        else:
            protein_seq_token = encoding.tokenize_sequence(protein_seq, self.vocab, True)
            protein_seq_masked_token = encoding.tokenize_sequence(protein_seq_masked, self.vocab, True)
            protein_seq_mask = np.insert(protein_seq_mask, 0, False)
            protein_seq_mask = np.append(protein_seq_mask, False)

        item['zm_label'] = torch.Tensor([zm_label]).bool()
        item['protein_seq'] = protein_seq_token
        item['protein_seq_masked'] = protein_seq_masked_token
        item['protein_seq_mask'] = torch.from_numpy(protein_seq_mask).bool()
        item['gene_name'] = meta_data['gene_name']

        locations = meta_data['locations']
        if pd.isna(locations):
            locations = "NA"
        location_list = [loc.strip() for loc in locations.split(',')]
        item['location_label'] = torch.Tensor([loc in location_list for loc in all_locations]).bool()
        item['locations'] = location_list
        item['protein_idx'] = torch.Tensor([self.gene_names.index(meta_data['gene_name'])]).long()

        return item

    def __len__(self) -> int:
        return len(self.meta_data)

    def collate(self, samples: List[dict]) -> dict:
        return collate_fn(
            samples, self.vocab, self.max_protein_sequence_len, 0,
            seq_length_control=self.seq_length_control,
        )


class OpenCell3DCropImageOnlyDataset(Dataset):
    """Image-only dataset for 3D VAE training.

    Flattens all image paths across all proteins into a flat list.
    Returns protein_img and nucleus_img as (1, D, H, W) tensors in [-1, 1],
    plus protein_idx for gene identity.  No sequence handling.
    """

    def __init__(self, args, split_key: str) -> None:
        super().__init__()
        self.data_path = args.data_path
        self.data_aug = args.data_aug

        meta_data = pd.read_csv(
            os.path.join(self.data_path, f'{split_key}_3d_192_crop_merged_meta_data.csv')
        )
        self.gene_names = sorted(meta_data['gene_name'].unique().tolist())

        self.img_paths: List[str] = []
        self.gene_name_labels: List[str] = []
        for row in meta_data.itertuples(index=False):
            for image_path in row.image_paths.split(','):
                self.img_paths.append(image_path.strip())
                self.gene_name_labels.append(row.gene_name)

    def get_img(self, image_path: str):
        img = tiff.imread(os.path.join(self.data_path, image_path)).astype(np.float32)
        nucleus_vol = _normalize_volume(img[:, 0, :, :])
        protein_vol = _normalize_volume(img[:, 1, :, :])

        nucleus_img = torch.from_numpy(nucleus_vol).unsqueeze(0)
        protein_img = torch.from_numpy(protein_vol).unsqueeze(0)

        if self.data_aug:
            stacked = torch.cat([protein_img, nucleus_img], dim=0)
            if random.random() < 0.5:
                stacked = stacked.flip(-1)
            if random.random() < 0.5:
                stacked = stacked.flip(-2)
            if random.random() < 0.5:
                stacked = stacked.flip(1)
            protein_img, nucleus_img = stacked[0:1], stacked[1:2]

        protein_img = protein_img * 2.0 - 1.0
        nucleus_img = nucleus_img * 2.0 - 1.0
        return protein_img, nucleus_img

    def __len__(self) -> int:
        return len(self.img_paths)

    def __getitem__(self, idx: int) -> dict:
        protein_img, nucleus_img = self.get_img(self.img_paths[idx])
        protein_idx = torch.Tensor([self.gene_names.index(self.gene_name_labels[idx])]).long()
        return {"protein_img": protein_img, "nucleus_img": nucleus_img, "protein_idx": protein_idx}

    @staticmethod
    def collate(samples: List[dict]) -> dict:
        protein_img = torch.stack([s["protein_img"] for s in samples])
        nucleus_img = torch.stack([s["nucleus_img"] for s in samples])
        protein_idx = torch.stack([s["protein_idx"] for s in samples])
        return {"batched_data": {"protein_img": protein_img, "nucleus_img": nucleus_img, "protein_idx": protein_idx}}


class OpenCell3DCropAllImageDataset(OpenCell3DCropDataset):
    """Returns all images for each protein as lists (instead of sampling one)."""

    def __getitem__(self, index: int) -> dict:
        meta_data = self.meta_data.iloc[index]
        item = {}

        item['protein_imgs'] = []
        item['nucleus_imgs'] = []
        for image_path in meta_data['image_paths'].split(','):
            protein_img, nucleus_img = self.get_img(image_path.strip())
            item['protein_imgs'].append(protein_img)
            item['nucleus_imgs'].append(nucleus_img)

        protein_seq = meta_data['sequence']
        protein_seq_masked, protein_seq_mask, zm_label = OAARDM_sequence_masking(protein_seq, self.seq_zero_mask_ratio)

        protein_seq_token = encoding.tokenize_sequence(protein_seq, self.vocab, True)
        protein_seq_masked_token = encoding.tokenize_sequence(protein_seq_masked, self.vocab, True)
        protein_seq_mask = np.insert(protein_seq_mask, 0, False)
        protein_seq_mask = np.append(protein_seq_mask, False)

        item['zm_label'] = torch.Tensor([zm_label]).bool()
        item['protein_seq'] = protein_seq_token
        item['protein_seq_masked'] = protein_seq_masked_token
        item['protein_seq_mask'] = torch.from_numpy(protein_seq_mask).bool()
        item['gene_name'] = meta_data['gene_name']

        locations = meta_data['locations']
        if pd.isna(locations):
            locations = "NA"
        location_list = [loc.strip() for loc in locations.split(',')]
        item['location_label'] = torch.Tensor([loc in location_list for loc in all_locations]).bool()
        item['locations'] = location_list
        item['protein_idx'] = torch.Tensor([self.gene_names.index(meta_data['gene_name'])]).long()

        return item

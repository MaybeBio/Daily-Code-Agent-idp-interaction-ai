# -*- coding: utf-8 -*-
import os
import random
from typing import List, Sequence, Union

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


def _parse_spatial_size(value: Union[str, Sequence[int]]) -> List[int]:
    """Accept '48,192,192' (raw argparse Namespace) or [48, 192, 192] (parsed config)."""
    if isinstance(value, str):
        return [int(s) for s in value.split(',')]
    return [int(s) for s in value]


def _normalize_volume(vol: torch.Tensor) -> torch.Tensor:
    """Clip to the 1% - 99% range through torch.quantile, then min-max normalize to [0, 1].

    Matches the raw-FOV normalization in cell_fm/data/opencell_data/dataset.py.
    A 48x192x192 crop is 1.77M voxels, well under torch.quantile's 2**24 input limit.
    """
    vol = torch.clamp(vol, torch.quantile(vol, 0.01), torch.quantile(vol, 0.99))
    vol = (vol - vol.min()) / max((vol.max() - vol.min()), 1e-8)
    return vol.clamp(0, 1)


def _crop_or_pad_axis(length: int, target: int, random_crop: bool):
    """Return (offset, pad_before, pad_after) to turn `length` into `target` along one axis."""
    if length >= target:
        offset = random.randint(0, length - target) if random_crop else (length - target) // 2
        return offset, 0, 0

    pad_before = (target - length) // 2
    return 0, pad_before, target - length - pad_before


def _crop_zyx(vol: np.ndarray, out_d: int, out_h: int, out_w: int, random_crop: bool) -> np.ndarray:
    """Crop (or edge-pad) a (D, C, H, W) volume to (out_d, C, out_h, out_w).

    Nucleus and protein share a single window, so the two channels stay registered.
    """
    d, _, h, w = vol.shape
    d_off, d_pre, d_post = _crop_or_pad_axis(d, out_d, random_crop)
    h_off, h_pre, h_post = _crop_or_pad_axis(h, out_h, random_crop)
    w_off, w_pre, w_post = _crop_or_pad_axis(w, out_w, random_crop)

    vol = vol[d_off: d_off + min(d, out_d),
              :,
              h_off: h_off + min(h, out_h),
              w_off: w_off + min(w, out_w)]

    if d_pre or d_post or h_pre or h_post or w_pre or w_post:
        vol = np.pad(vol, ((d_pre, d_post), (0, 0), (h_pre, h_post), (w_pre, w_post)), mode='edge')

    return vol


class _RawVolumeLoaderMixin:
    """Shared raw-FOV loading: partial z-read, joint crop, per-crop normalization, flips."""

    def _init_img_loader(self, args) -> None:
        self.data_path = args.data_path
        self.data_aug = args.data_aug
        self.input_spatial_size = _parse_spatial_size(args.input_spatial_size)

    def get_img(self, image_path: str):
        """Load a raw FOV and return (protein_img, nucleus_img), each (1, D, H, W) in [-1, 1]."""
        out_d, out_h, out_w = self.input_spatial_size
        random_crop = self.data_aug

        with tiff.TiffFile(os.path.join(self.data_path, image_path)) as tf:
            d, c, _, _ = tf.series[0].shape

            # Read only the z-slices we keep: pages are flattened over (D, C).
            z0, _, _ = _crop_or_pad_axis(d, out_d, random_crop)
            n_z = min(d, out_d)
            vol = tf.asarray(key=range(z0 * c, (z0 + n_z) * c))

        vol = vol.reshape(n_z, c, *vol.shape[-2:]).astype(np.float32)
        # z is already selected above; _crop_zyx handles H/W (and pads z only if a stack were short)
        vol = _crop_zyx(vol, out_d, out_h, out_w, random_crop)

        # (D, H, W) → (1, D, H, W)
        nucleus_img = torch.from_numpy(np.ascontiguousarray(vol[:, 0])).unsqueeze(0)
        protein_img = torch.from_numpy(np.ascontiguousarray(vol[:, 1])).unsqueeze(0)

        nucleus_img = _normalize_volume(nucleus_img)
        protein_img = _normalize_volume(protein_img)

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


class OpenCell3DDataset(_RawVolumeLoaderMixin, Dataset):
    """Dataset for raw 3D OpenCell fields of view.

    Each TIF has shape (D, 2, 600, 600) with ragged D (51-131): channel 0 = nucleus,
    channel 1 = protein.  A 48x192x192 window is cropped on the fly — randomly when
    data_aug is on, centred otherwise.  Returns protein_img and nucleus_img as
    (1, 48, 192, 192) tensors in [-1, 1].
    """

    def __init__(self, args, split_key: str) -> None:
        super().__init__()
        self.args = args
        self.split_key = split_key

        self.vocab = EsmSequenceTokenizer()
        # Single metadata CSV covers all proteins (no per-split files)
        self.meta_data = pd.read_csv(
            os.path.join(args.data_path, f'{split_key}_3d_merged_meta_data.csv')
        )

        self._init_img_loader(args)
        self.seq_zero_mask_ratio = args.seq_zero_mask_ratio
        self.max_protein_sequence_len = args.max_protein_sequence_len
        self.pre_pad_seq = getattr(args, 'pre_pad_seq', False)
        self.phase = args.phase

        self.gene_names = sorted(self.meta_data['gene_name'].unique().tolist())

    def __getitem__(self, index: int) -> dict:
        if self.phase == 'train':
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
        return collate_fn(samples, self.vocab, self.max_protein_sequence_len, 0)


class OpenCell3DImageOnlyDataset(_RawVolumeLoaderMixin, Dataset):
    """Image-only dataset for 3D VAE training on raw fields of view.

    Flattens all image paths across all proteins into a flat list.
    Returns protein_img and nucleus_img as (1, D, H, W) tensors in [-1, 1],
    plus protein_idx for gene identity.  No sequence handling.
    """

    def __init__(self, args, split_key: str) -> None:
        super().__init__()
        self._init_img_loader(args)

        meta_data = pd.read_csv(
            os.path.join(args.data_path, f'{split_key}_3d_merged_meta_data.csv')
        )
        self.gene_names = sorted(meta_data['gene_name'].unique().tolist())

        self.img_paths: List[str] = []
        self.gene_name_labels: List[str] = []
        for row in meta_data.itertuples(index=False):
            for image_path in row.image_paths.split(','):
                self.img_paths.append(image_path.strip())
                self.gene_name_labels.append(row.gene_name)

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


class OpenCell3DAllImageDataset(OpenCell3DDataset):
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

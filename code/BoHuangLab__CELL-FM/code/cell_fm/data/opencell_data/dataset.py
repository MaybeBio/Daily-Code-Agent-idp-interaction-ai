# -*- coding: utf-8 -*-
import os
from typing import List
import pandas as pd
import random

import numpy as np
import torch

from torch.utils.data import Dataset

from .collater import collate_fn
from esm.utils import encoding
from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer

from cell_fm.data.utils.seq_utils import OAARDM_sequence_masking
from cell_fm.data.utils.img_utils import RandomRotation

from torchvision import transforms
from torchvision.transforms.functional import to_tensor
import tifffile as tiff
from .locations import all_locations
from pathlib import Path


class OpenCellDataset(Dataset):
    def __init__(self, args, split_key) -> None:
        super().__init__()
        self.args = args

        self.split_key = split_key
        self.data_path = self.args.data_path

        self.vocab = EsmSequenceTokenizer()
        self.meta_data = pd.read_csv(os.path.join(self.data_path, f'{self.split_key}_merged_meta_data.csv'))

        self.data_aug = self.args.data_aug
        self.img_resize = self.args.img_resize
        self.img_crop_size = self.args.img_crop_size
        self.seq_zero_mask_ratio = self.args.seq_zero_mask_ratio
        self.max_protein_sequence_len = self.args.max_protein_sequence_len
        self.pre_pad_seq = self.args.pre_pad_seq
        self.phase = self.args.phase

        self.gene_names = sorted(self.meta_data['gene_name'].unique().tolist())

        self.transform = self._build_transform()

    def _build_transform(self) -> transforms.Compose:

        # To maintain the resolution.
        # First crop then resize.

        t_forms = []

        if self.data_aug:
            t_forms.append(transforms.RandomCrop(self.img_crop_size))
            t_forms.append(transforms.Resize(self.img_resize, antialias=None))
            t_forms.append(transforms.RandomHorizontalFlip(p=0.5))
            t_forms.append(RandomRotation([0, 90, 180, 270]))
        else:
            t_forms.append(transforms.CenterCrop(self.img_crop_size))
            t_forms.append(transforms.Resize(self.img_resize, antialias=None))

        t_forms = transforms.Compose(t_forms)

        return t_forms

    def __getitem__(self, index: int) -> dict:
        if self.phase == 'train':
            while True:
                meta_data = self.meta_data.iloc[index]
                protein_seq = meta_data['sequence']
                if len(protein_seq) <= self.max_protein_sequence_len:
                    break
                index = random.randint(0, self.__len__() - 1)
        else:
            meta_data = self.meta_data.iloc[index]

        item = {}

        image_paths = meta_data['image_paths'].split(',')
        image_path = random.choice(image_paths).strip()
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
            """
            - convert string sequence to int index
            """
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
        location_label = [loc in location_list for loc in all_locations]
        item['location_label'] = torch.Tensor(location_label).bool()
        item['locations'] = location_list

        item['protein_idx'] = torch.Tensor([self.gene_names.index(meta_data['gene_name'])]).long()

        return item

    def get_img(self, image_path):
        img = tiff.imread(os.path.join(self.data_path, image_path)).astype(np.float32)        
        
        nucleus_img = img[0]
        protein_img = img[1]

        nucleus_img = to_tensor(nucleus_img)
        protein_img = to_tensor(protein_img)

        img = torch.stack([protein_img, nucleus_img], dim=0)
        protein_img, nucleus_img = self.transform(img)

        # clip the images to 1% - 99% range through torch.quantile
        nucleus_img = torch.clamp(nucleus_img, torch.quantile(nucleus_img, 0.01), torch.quantile(nucleus_img, 0.99))
        protein_img = torch.clamp(protein_img, torch.quantile(protein_img, 0.01), torch.quantile(protein_img, 0.99))

        # normalize the images to [0, 1]
        nucleus_img = (nucleus_img - nucleus_img.min()) / max((nucleus_img.max() - nucleus_img.min()), 1e-8)
        protein_img = (protein_img - protein_img.min()) / max((protein_img.max() - protein_img.min()), 1e-8)

        nucleus_img = nucleus_img.clamp(0, 1)
        protein_img = protein_img.clamp(0, 1)

        nucleus_img = nucleus_img * 2 - 1
        protein_img = protein_img * 2 - 1

        return protein_img, nucleus_img

    def __len__(self) -> int:
        return self.meta_data.shape[0]

    def collate(self, samples: List[dict]) -> dict:
        return collate_fn(samples, self.vocab, self.max_protein_sequence_len, 0)


class OpenCellImageOnlyDataset(Dataset):
    def __init__(self, args, split_key: str):
        super().__init__()
        self.args = args
        self.split_key = split_key
        self.data_path = Path(args.data_path)
        self.img_resize = args.img_resize
        self.img_crop_size = self.args.img_crop_size
        self.data_aug = args.data_aug

        self.meta_data = pd.read_csv(os.path.join(self.data_path, f'{self.split_key}_merged_meta_data.csv'))
        self.gene_names = sorted(self.meta_data['gene_name'].unique().tolist())

        # gather image paths
        self.img_paths: List[Path] = []
        self.gene_name_labels = []
        for row in self.meta_data.itertuples(index=False):
            image_paths = row.image_paths.split(',')
            for image_path in image_paths:
                self.img_paths.append(self.data_path / image_path.strip())
                self.gene_name_labels.append(row.gene_name)

        # build transform once
        self.transform = self._build_transform()

    def _build_transform(self) -> transforms.Compose:

        # To maintain the resolution.
        # First crop then resize.

        t_forms = []

        if self.data_aug:
            t_forms.append(transforms.RandomCrop(self.img_crop_size))
            t_forms.append(transforms.Resize(self.img_resize, antialias=None))
            t_forms.append(transforms.RandomHorizontalFlip(p=0.5))
            t_forms.append(RandomRotation([0, 90, 180, 270]))
        else:
            t_forms.append(transforms.CenterCrop(self.img_crop_size))
            t_forms.append(transforms.Resize(self.img_resize, antialias=None))

        t_forms = transforms.Compose(t_forms)

        return t_forms
    
    def get_img(self, image_path):
        img = tiff.imread(os.path.join(self.data_path, image_path)).astype(np.float32)        
        
        nucleus_img = img[0]
        protein_img = img[1]

        nucleus_img = to_tensor(nucleus_img)
        protein_img = to_tensor(protein_img)

        img = torch.stack([protein_img, nucleus_img], dim=0)
        protein_img, nucleus_img = self.transform(img)

        # clip the images to 1% - 99% range through torch.quantile
        nucleus_img = torch.clamp(nucleus_img, torch.quantile(nucleus_img, 0.01), torch.quantile(nucleus_img, 0.99))
        protein_img = torch.clamp(protein_img, torch.quantile(protein_img, 0.01), torch.quantile(protein_img, 0.99))

        # normalize the images to [0, 1]
        nucleus_img = (nucleus_img - nucleus_img.min()) / max((nucleus_img.max() - nucleus_img.min()), 1e-8)
        protein_img = (protein_img - protein_img.min()) / max((protein_img.max() - protein_img.min()), 1e-8)

        nucleus_img = nucleus_img.clamp(0, 1)
        protein_img = protein_img.clamp(0, 1)

        nucleus_img = nucleus_img * 2 - 1
        protein_img = protein_img * 2 - 1

        return protein_img, nucleus_img

    def __len__(self) -> int:
        return len(self.img_paths)

    def __getitem__(self, idx: int) -> dict:
        item = {}
        item['protein_img'], item['nucleus_img'] = self.get_img(self.img_paths[idx])
        item['protein_idx'] = torch.Tensor([self.gene_names.index(self.gene_name_labels[idx])]).long()

        return item

    @staticmethod
    def collate(samples):
        protein_img = torch.stack([s["protein_img"] for s in samples], 0)
        nucleus_img = torch.stack([s["nucleus_img"] for s in samples], 0)
        protein_idx = torch.stack([s["protein_idx"] for s in samples], 0)

        return {"batched_data": {"protein_img": protein_img, "nucleus_img": nucleus_img, "protein_idx": protein_idx}}
# -*- coding: utf-8 -*-
import os
from typing import List
from pathlib import Path
import pandas as pd
import random

import numpy as np
import torch

from torch.utils.data import Dataset

from .collater import collate_fn
from esm.utils import encoding
from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer

from cell_fm.data.utils.seq_utils import OAARDM_sequence_masking
from cell_fm.data.utils.img_utils import RandomRotation, CropOrPad

from torchvision import transforms
from torchvision.transforms.functional import to_tensor
import tifffile as tiff
from glob import glob


class CondenSeqDataset(Dataset):
    def __init__(self, args, split_key) -> None:
        super().__init__()
        self.args = args

        self.split_key = split_key
        self.data_path = self.args.data_path

        self.vocab = EsmSequenceTokenizer()
        self.meta_data = pd.read_csv(os.path.join(self.data_path, f'{self.split_key}_meta_data_large_pool_updated.csv'))

        self.data_aug = self.args.data_aug
        self.img_resize = self.args.img_resize
        self.seq_zero_mask_ratio = self.args.seq_zero_mask_ratio
        self.max_protein_sequence_len = self.args.max_protein_sequence_len
        self.img_type = self.args.img_type
        self.phase = self.args.phase

        self.background_value_gfp = 133.3
        self.background_value_snap = 124.3

        self.index_to_path = {}

    def __getitem__(self, index: int) -> dict:

        if self.phase == 'train':
            while True:
                meta_data = self.meta_data.iloc[index]
                protein_seq = meta_data['protein_seq']
                if len(protein_seq) <= self.max_protein_sequence_len:
                    break
                index = random.randint(0, self.__len__() - 1)
        else:
            meta_data = self.meta_data.iloc[index]

        item = {}

        if index not in self.index_to_path:
            if self.img_type == 'GFP':
                barcode_paths = meta_data['GFP_paths'].split(',')
            elif self.img_type == 'SNAP':
                barcode_paths = meta_data['SNAP_paths'].split(',')

            all_gfp_paths = []
            all_dapi_paths = []
            all_mask_paths = []

            for barcode_path in barcode_paths:
                gfp_paths = sorted(glob(os.path.join(self.data_path, barcode_path.strip(), '*_GFP.tif')))
                dapi_paths = sorted(glob(os.path.join(self.data_path, barcode_path.strip(), '*_DAPI.tif')))
                mask_paths = sorted(glob(os.path.join(self.data_path, barcode_path.strip(), '*_mask.npy')))

                all_gfp_paths.extend(gfp_paths)
                all_dapi_paths.extend(dapi_paths)
                all_mask_paths.extend(mask_paths)

            image_paths = list(zip(all_gfp_paths, all_dapi_paths, all_mask_paths))
            self.index_to_path[index] = image_paths
        else:
            image_paths = self.index_to_path[index]

        assert len(image_paths) > 0, f"No images found for {self.split_key} split with img_type {self.img_type} with index {meta_data['index']}."

        # randomly choose one image_path
        image_path = random.choice(image_paths)

        # get images
        item['protein_img'], item['nucleus_img'], item['mask_img'], item['protein_intensity_level'] = self.get_img(image_path)

        protein_seq = meta_data['protein_seq']
        protein_seq_masked, protein_seq_mask, zm_label = OAARDM_sequence_masking(protein_seq, self.seq_zero_mask_ratio)

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
        item['protein_seq_mask'] = torch.from_numpy(protein_seq_mask).long()
        item['protein_intensity_level'] = torch.Tensor([item['protein_intensity_level']]).float()
        item['index'] = meta_data['index']

        return item

    def get_img(self, image_path):
        gfp_path, dapi_path, mask_path = image_path

        gfp_img = tiff.imread(gfp_path)
        dapi_img = tiff.imread(dapi_path)
        mask_img = np.load(mask_path)

        gfp_img = gfp_img.astype(np.float32)
        dapi_img = dapi_img.astype(np.float32)

        if self.img_type == "GFP":
            protein_intensity_level = np.mean(gfp_img[mask_img]) - self.background_value_gfp
        elif self.img_type == "SNAP":
            protein_intensity_level = np.mean(gfp_img[mask_img]) - self.background_value_snap

        protein_img = gfp_img
        nucleus_img = dapi_img

        # normalize the images to [0, 1]
        nucleus_img = (nucleus_img - nucleus_img.min()) / (nucleus_img.max() - nucleus_img.min())
        protein_img = (protein_img - protein_img.min()) / (protein_img.max() - protein_img.min())

        nucleus_img = to_tensor(nucleus_img)
        protein_img = to_tensor(protein_img)
        mask_img = torch.from_numpy(mask_img).float().unsqueeze(0)

        t_forms = []

        t_forms.append(CropOrPad(
            target_height=self.img_resize,
            target_width=self.img_resize,
            crop_mode='center',
            pad_mode='constant',
            pad_value=0.0,
        ))

        if self.data_aug:
            t_forms.append(transforms.RandomHorizontalFlip(p=0.5))
            t_forms.append(RandomRotation([0, 90, 180, 270]))

        t_forms = transforms.Compose(t_forms)

        img = torch.stack([protein_img, nucleus_img, mask_img], dim=0)
        protein_img, nucleus_img, mask_img = t_forms(img)

        # normalize the images to [-1, 1]
        protein_img = 2 * protein_img - 1
        nucleus_img = 2 * nucleus_img - 1

        return protein_img, nucleus_img, mask_img, protein_intensity_level

    def __len__(self) -> int:
        return self.meta_data.shape[0]

    def collate(self, samples: List[dict]) -> dict:
        return collate_fn(samples, self.vocab, self.max_protein_sequence_len, 0)


class CondenSeqAllImageDataset(CondenSeqDataset):
    def __getitem__(self, index: int) -> dict:
        meta_data = self.meta_data.iloc[index]

        item = {}

        if self.img_type == 'GFP':
            barcode_paths = meta_data['GFP_paths'].split(',')
        elif self.img_type == 'SNAP':
            barcode_paths = meta_data['SNAP_paths'].split(',')

        all_gfp_paths = []
        all_dapi_paths = []
        all_mask_paths = []

        for barcode_path in barcode_paths:
            gfp_paths = sorted(glob(os.path.join(self.data_path, barcode_path.strip(), '*_GFP.tif')))
            dapi_paths = sorted(glob(os.path.join(self.data_path, barcode_path.strip(), '*_DAPI.tif')))
            mask_paths = sorted(glob(os.path.join(self.data_path, barcode_path.strip(), '*_mask.npy')))

            all_gfp_paths.extend(gfp_paths)
            all_dapi_paths.extend(dapi_paths)
            all_mask_paths.extend(mask_paths)

        image_paths = list(zip(all_gfp_paths, all_dapi_paths, all_mask_paths))
        assert len(image_paths) > 0, f"No images found for {self.split_key} split with img_type {self.img_type} with index {meta_data['index']}."

        item['protein_imgs'] = []
        item['nucleus_imgs'] = []
        item['mask_imgs'] = []
        item['protein_intensity_levels'] = []

        for image_path in image_paths:
            protein_img, nucleus_img, mask_img, protein_intensity_level = self.get_img(image_path)
            item['protein_imgs'].append(protein_img)
            item['nucleus_imgs'].append(nucleus_img)
            item['mask_imgs'].append(mask_img)
            item['protein_intensity_levels'].append(protein_intensity_level)

        protein_seq = meta_data['protein_seq']
        protein_seq_masked, protein_seq_mask, zm_label = OAARDM_sequence_masking(protein_seq, self.seq_zero_mask_ratio)

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
        item['protein_seq_mask'] = torch.from_numpy(protein_seq_mask).long()
        item['index'] = meta_data['index']

        return item


class CondenSeqProteinImageOnlyDataset(Dataset):
    def __init__(self, args, split_key: str):
        super().__init__()
        self.args = args
        self.split_key = split_key
        self.data_path = Path(args.data_path)
        self.img_resize = args.img_resize
        self.data_aug = args.data_aug

        # read metadata
        meta_csv = self.data_path / f"{split_key}_meta_data_large_pool.csv"
        self.meta_data = pd.read_csv(meta_csv)

        # # choose the first few rows for testing
        # self.meta_data = self.meta_data.head(100)

        # gather image paths
        self.img_paths: List[Path] = []
        for row in self.meta_data.itertuples(index=False):
            for paths in (row.GFP_paths, row.SNAP_paths):
                for p in paths.split(","):
                    self.img_paths.extend(
                        sorted((self.data_path / p.strip()).glob("*_GFP.tif"))
                    )

        # build transform once
        self.transform = self._build_transform()

    # ------------------------------------------------------------------ #
    #                         helper functions                           #
    # ------------------------------------------------------------------ #
    def _build_transform(self) -> transforms.Compose:
        ops = [CropOrPad(self.img_resize, self.img_resize, crop_mode="center")]
        if self.data_aug:
            ops += [
                transforms.RandomHorizontalFlip(0.5),
                RandomRotation([0, 90, 180, 270]),
            ]
        ops += [transforms.Normalize([0.5], [0.5])]
        return transforms.Compose(ops)

    @staticmethod
    def _load_tiff(fpath: Path) -> torch.Tensor:
        z = tiff.imread(fpath).astype(np.float32)
        z = (z - z.min()) / (z.ptp() + 1e-6)
        return torch.from_numpy(z).unsqueeze(0)  # (1, H, W)

    # ------------------------------------------------------------------ #
    #                       Dataset interface                             #
    # ------------------------------------------------------------------ #
    def __len__(self) -> int:
        return len(self.img_paths)

    def __getitem__(self, idx: int) -> dict:
        img = self._load_tiff(self.img_paths[idx])
        return {"protein_img": self.transform(img)}

    @staticmethod
    def collate(samples):
        imgs = torch.stack([s["protein_img"] for s in samples], 0)
        return {"batched_data": {"protein_img": imgs}}


class CondenSeqIDRImageDataset(CondenSeqDataset):
    """
    Single-cell images labeled by IDR id, for the ViT IDR classifier (vit_cls).
    Each item is one random cell of one IDR; the label is the IDR's row position in the split csv.
    """
    def get_image_paths(self, meta_data):
        if self.img_type == 'GFP':
            barcode_paths = meta_data['GFP_paths'].split(',')
        elif self.img_type == 'SNAP':
            barcode_paths = meta_data['SNAP_paths'].split(',')

        image_paths = []
        for barcode_path in barcode_paths:
            for gfp_path in sorted(glob(os.path.join(self.data_path, barcode_path.strip(), '*_GFP.tif'))):
                stem = gfp_path[:-len('_GFP.tif')]
                image_paths.append((gfp_path, stem + '_DAPI.tif', stem + '_mask.npy'))

        assert len(image_paths) > 0, f"No images found for {self.split_key} split with img_type {self.img_type} with index {meta_data['index']}."

        return image_paths

    def __getitem__(self, index: int) -> dict:
        image_paths = self.get_image_paths(self.meta_data.iloc[index])

        item = {}
        item['protein_img'], item['nucleus_img'], _, _ = self.get_img(random.choice(image_paths))
        item['protein_idx'] = torch.Tensor([index]).long()

        return item

    @staticmethod
    def collate(samples):
        protein_img = torch.stack([s["protein_img"] for s in samples], 0)
        nucleus_img = torch.stack([s["nucleus_img"] for s in samples], 0)
        protein_idx = torch.stack([s["protein_idx"] for s in samples], 0)

        return {"batched_data": {"protein_img": protein_img, "nucleus_img": nucleus_img, "protein_idx": protein_idx}}


class CondenSeqIDRImageEmbedDataset(CondenSeqIDRImageDataset):
    """
    All cells (or up to max_cells_per_idr random cells) of every IDR in the split, for embedding.
    """
    def __init__(self, args, split_key, max_cells_per_idr=-1, seed=0) -> None:
        super().__init__(args, split_key)

        rng = random.Random(seed)

        self.cells = []
        for index in range(len(self.meta_data)):
            image_paths = self.get_image_paths(self.meta_data.iloc[index])
            if 0 < max_cells_per_idr < len(image_paths):
                image_paths = rng.sample(image_paths, max_cells_per_idr)
            self.cells.extend((image_path, index) for image_path in image_paths)

    def __getitem__(self, idx: int) -> dict:
        image_path, index = self.cells[idx]

        item = {}
        item['protein_img'], item['nucleus_img'], _, _ = self.get_img(image_path)
        item['protein_idx'] = torch.Tensor([index]).long()

        return item

    def __len__(self) -> int:
        return len(self.cells)

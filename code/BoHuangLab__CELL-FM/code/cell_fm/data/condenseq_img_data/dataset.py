# -*- coding: utf-8 -*-
import os
from typing import List
import pandas as pd

import numpy as np
import torch

from torch.utils.data import Dataset
from cell_fm.data.utils.img_utils import RandomRotation, CropOrPad

from torchvision import transforms
from torchvision.transforms.functional import to_tensor
import tifffile as tiff

class CondenSeqImageDataset(Dataset):
    def __init__(self, args, split_key) -> None:
        super().__init__()
        self.args = args

        self.split_key = split_key
        self.data_path = self.args.data_path

        self.meta_data = pd.read_csv(os.path.join(self.data_path, f'{self.split_key}_image_labeled.csv'))

        self.data_aug = self.args.data_aug
        self.img_resize = self.args.img_resize

    def __getitem__(self, index: int) -> dict:

        meta_data = self.meta_data.iloc[index]

        item = {}
        
        gfp_path = meta_data['gfp_path']
        dapi_path = meta_data['dapi_path']
        class_label = int(meta_data['label_idx'])

        if class_label == 2:
            label = 1 # condensate
        else:
            label = 0 # non-condensate

        item['protein_img'], item['nucleus_img'] = self.get_img(gfp_path, dapi_path)
        item['label'] = torch.tensor(label, dtype=torch.long)

        return item

    def get_img(self, gfp_path, dapi_path):

        gfp_img = tiff.imread(os.path.join(self.data_path, gfp_path)).astype(np.float32)
        dapi_img = tiff.imread(os.path.join(self.data_path, dapi_path)).astype(np.float32)

        protein_img = gfp_img
        nucleus_img = dapi_img

        # normalize the images to [0, 1]
        nucleus_img = (nucleus_img - nucleus_img.min()) / (nucleus_img.max() - nucleus_img.min())
        protein_img = (protein_img - protein_img.min()) / (protein_img.max() - protein_img.min())

        nucleus_img = to_tensor(nucleus_img)
        protein_img = to_tensor(protein_img)

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

        img = torch.stack([protein_img, nucleus_img], dim=0)
        protein_img, nucleus_img = t_forms(img)

        protein_img = protein_img.clamp(0.0, 1.0)
        nucleus_img = nucleus_img.clamp(0.0, 1.0)

        # normalize the images to [-1, 1]
        protein_img = 2 * protein_img - 1
        nucleus_img = 2 * nucleus_img - 1

        return protein_img, nucleus_img

    def __len__(self) -> int:
        return self.meta_data.shape[0]

    def collate(self, samples: List[dict]) -> dict:
        batch = dict()

        batch["protein_img"] = torch.cat(
            [s["protein_img"].unsqueeze(0) for s in samples]
        )

        batch["nucleus_img"] = torch.cat(
            [s["nucleus_img"].unsqueeze(0) for s in samples]
        )

        batch["label"] = torch.cat(
            [s["label"].unsqueeze(0) for s in samples]
        )

        return {'batched_data': batch}



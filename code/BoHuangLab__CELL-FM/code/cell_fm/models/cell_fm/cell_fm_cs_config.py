# -*- coding: utf-8 -*-
from dataclasses import dataclass
from cell_fm.models.cell_fm.cell_fm_config import CELLFMConfig


@dataclass
class CELLFMCSConfig(CELLFMConfig):
    img_type: str = "GFP"  # "GFP", "SNAP"

    def __init__(self, **kwargs):
        self.img_type = kwargs.pop("img_type", "GFP")
        super().__init__(**kwargs)

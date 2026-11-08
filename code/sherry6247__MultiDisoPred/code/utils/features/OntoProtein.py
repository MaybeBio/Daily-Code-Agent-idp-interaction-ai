'''
Author: Sicen Liu
Date: 2024-12-05 11:46:30
LastEditTime: 2025-06-18 20:08:06
FilePath: /liusicen/shared_files/MOFE/utils/features/OntoProtein.py
Description: 

Copyright (c) 2024 by ${liusicen_cs@outlook.com}, All Rights Reserved. 
'''
import numpy as np
from transformers import AutoTokenizer, AutoModelForMaskedLM, FeatureExtractionPipeline
import torch

# tokenizer_onto = AutoTokenizer.from_pretrained("zjukg/OntoProtein")
# model_onto = AutoModelForMaskedLM.from_pretrained("zjukg/OntoProtein")

print("loading OntoProtein model...")
path = "/home/suhanghua/MoSE/Datasets/preTrainModels/OntoProtein"
tokenizer_onto = AutoTokenizer.from_pretrained(path)
model_onto = AutoModelForMaskedLM.from_pretrained(path)
bert_model_onto = FeatureExtractionPipeline(model=model_onto, tokenizer=tokenizer_onto)


def encode_by_onto(seq):
    # maxsize = 500
    seq = " ".join(list(seq))
    res = bert_model_onto(seq)[0][1:-1]
    res = np.squeeze(res)#[:maxsize]
    # res = np.pad(res, ((0, 500 - len(res)), (0, 0)), constant_values=0)
    return res


"""
seq = "ARCYCS" * 100
x = encode_by_onto(seq)
print(x.shape)
"""

# seq = "ATCYCRTGRCATRESLSGVCRISGRLYRLCCR"
# x = encode_by_onto(seq)
# print(x.shape)


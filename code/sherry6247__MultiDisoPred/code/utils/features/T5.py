'''
Author: Sicen Liu
Date: 2024-12-05 11:46:30
LastEditTime: 2025-06-18 20:08:00
FilePath: /liusicen/shared_files/MOFE/utils/features/T5.py
Description: 

Copyright (c) 2024 by ${liusicen_cs@outlook.com}, All Rights Reserved. 
'''
import re
import torch
import numpy as np
import torch.nn as nn

from transformers import T5Tokenizer, T5Model

print("loading T5...")
# tokenizer_t5 = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_uniref50', do_lower_case=False)
# model_t5 = T5Model.from_pretrained("Rostlab/prot_t5_xl_uniref50")

path = "/home/suhanghua/MoSE/Datasets/preTrainModels/protTrans/prot_t5_xl_uniref50"
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
tokenizer_t5 = T5Tokenizer.from_pretrained(path, do_lower_case=False)
model_t5 = T5Model.from_pretrained(path)#.to(device)
# if torch.cuda.device_count() > 1:
#     model_t5 = nn.DataParallel(model_t5)  # 使用 DataParallel

def encode_by_t5(seq):
    sequences_Example = [" ".join(list(seq.upper()))]
    sequences_Example = [re.sub(r"[UZOB]", "X", sequence) for sequence in sequences_Example]
    ids = tokenizer_t5.batch_encode_plus(sequences_Example,
                                         add_special_tokens=True,
                                         padding=False)

    input_ids = torch.tensor(ids['input_ids'])#.to(device)
    attention_mask = torch.tensor(ids['attention_mask'])#.to(device)

    with torch.no_grad():
        embedding = model_t5(input_ids=input_ids,
                             attention_mask=attention_mask,
                             decoder_input_ids=input_ids,
                             )

    # maxsize = 500
    # embedding_size = 600

    encoder_embedding = embedding.last_hidden_state[0,:-1].detach().cpu().numpy()
    # return encoder_embedding[:maxsize, :embedding_size]
    return encoder_embedding


"""
seq = "ARC"*200
x = encode_by_t5(seq)
print(x.shape)
"""

# seq = "ARC"*200
# x = encode_by_t5(seq)
# print(x.shape)
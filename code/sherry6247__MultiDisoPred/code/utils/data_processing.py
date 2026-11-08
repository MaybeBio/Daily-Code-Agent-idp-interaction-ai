'''
Author: Sicen Liu
Date: 2024-12-10 09:12:48
LastEditTime: 2026-03-03 00:00:52
FilePath: /liusicen/methods/MoSE/utils/data_processing.py
Description: Data processing code

Copyright (c) 2024 by ${liusicen_cs@outlook.com}, All Rights Reserved. 
'''
import sys, getopt
import torch
import numpy as np
import re
from transformers import AutoModel
from transformers import AutoTokenizer
import os
import textwrap
from tqdm import tqdm
import gc

from utils.features.traditionalFeature import  encode_by_pssm
from torch.utils.data import Dataset

from sklearn.model_selection import KFold
ESM2_PATH = "/home/suhanghua/MoSE/Datasets/preTrainModels/esm-main"

class fastaExample(object):
    def __init__(self, id, name, seq, label):
        self.id = id
        self.name = name
        self.seq = seq
        self.label = label

class CAID3_fastaExample(object):
    def __init__(self, id, name, seq, label, unexcluded_label_idx):
        self.id = id
        self.name = name
        self.seq = seq
        self.label = label
        self.unexcluded_label_idx = unexcluded_label_idx

def parseFasta_with_label(input_file):
    samples = []
    load_f = []
    with open(input_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line=line.strip('\n')
            if line:
                load_f.append(line)
    i = 0
    while i < len(load_f):
        if load_f[i].startswith('>'):
            id = load_f[i].strip('>').strip('\n\r')
            seq = load_f[i+1] if i+1 < len(load_f) else ''
            if i+2 < len(load_f) and not load_f[i+2].startswith('>'):
                label = load_f[i+2]
                i += 3
            else:
                label = '0' * len(seq)
                i += 2
            samples.append(fastaExample(id, id, seq, label))
        else:
            i += 1
    return samples 

def parseCAID3_Fasta_with_label(input_file):
    samples = []
    load_f = []
    with open(input_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            line=line.strip('\n')
            load_f.append(line)
    for i in range(len(load_f)):
        if i % 3 == 0:
            id = load_f[i].strip('>').strip('\n\r')
            seq = load_f[i+1]
            label = load_f[i+2]
            
            unexcluded_label_idx = []
            temp_label = []
            for i, (s, l) in enumerate(zip(seq, label)):
                if l == "-" or l== '2':
                    temp_label.append('0')
                    continue
                else:
                    temp_label.append(l)
                    unexcluded_label_idx.append(i)
            temp_label = ''.join(temp_label)
            samples.append(CAID3_fastaExample(id, id, seq, temp_label, unexcluded_label_idx))
    return samples 


from typing import List, Tuple, Any, Optional
import torch
def _contiguous_spans_from_labels(labels) -> List[Tuple[int, int]]:
    """把逐残基 0/1 标签转换为连续的 [start, end) 片段（1 段）。"""
    spans = []
    start = None
    # 加一个哨兵 0 以便在末尾 flush
    for i, v in enumerate(list(map(int, labels)) + [0]):
        if v == 1 and start is None:
            start = i
        elif v != 1 and start is not None:
            spans.append((start, i))
            start = None
    return spans
def collect_span_embeddings_from_dataset(
    dataset: List[Any],
    short_long_threshold: int = 30,
    device: Optional[str] = None,
) -> Tuple[torch.Tensor, torch.LongTensor]:
    """
    从 dataset（元素含属性 id, name, seq, features, label）提取：
      - span_embeddings: (N_spans, D)   每个无序片段的 mean 表示
      - span_labels:     (N_spans,)     短=0 (<threshold), 长=1 (>=threshold)

    要求：
      - features: (L, D) 的 per-residue embedding（Tensor/NDArray/可转 Tensor）
      - label:    长度为 L 的 0/1 序列，1 表示无序
    """
    span_vecs = []
    span_labs = []

    for rec in dataset:
        seq = getattr(rec, "seq", None)
        feats = getattr(rec, "features", None)   # 期望 (L, D)
        labels = getattr(rec, "label", None)     # 期望 长度 L 的 0/1

        if seq is None or feats is None or labels is None:
            raise ValueError("Each record must have attributes: seq, features, label")

        # to tensor
        E = torch.as_tensor(feats, device=device).float()
        if E.ndim == 3 and E.size(0) == 1:  # 兼容 (1,L,D)
            E = E.squeeze(0)
        if E.ndim != 2:
            raise ValueError(f"features must be (L,D), got {tuple(E.shape)}")

        L, D = E.shape
        # 对齐长度，避免不一致
        L_aligned = min(L, len(labels), len(seq))
        if L_aligned != L:
            E = E[:L_aligned]
        labels = list(labels)[:L_aligned]

        # 取出所有无序片段
        spans = _contiguous_spans_from_labels(labels)
        for s, e in spans:
            if not (0 <= s < e <= L_aligned):
                continue
            vec = E[s:e].mean(dim=0)                   # (D,)
            length = e - s
            lab = 0 if length < short_long_threshold else 1
            span_vecs.append(vec.detach().cpu())
            span_labs.append(lab)

    if not span_vecs:
        # 没有无序片段时给出明确提示
        raise RuntimeError("No disorder spans (label=1) found in the dataset.")

    span_embeddings = torch.stack(span_vecs, dim=0)          # (N_spans, D)
    span_labels = torch.tensor(span_labs, dtype=torch.long)  # (N_spans,)
    return span_embeddings, span_labels



"""--------------------------------------------------------"""

device = "cuda" if torch.cuda.is_available() else "cpu"
# ESM-2 model will be loaded only when needed
model_esm2 = None
alphabet_esm2 = None
"""--------------------------------------------------------"""
# DR-BERT model will be loaded only when needed
tokenizer_drbert = None
model_drbert = None
"""--------------------------------------------------------"""
def encode_by_drbert(seq):
    global tokenizer_drbert, model_drbert
    # Load DR-BERT model if not already loaded
    if tokenizer_drbert is None or model_drbert is None:
        print("loading DR-BERT...")
        # model_name_drbert = "lib/DR-BERT"
        model_name_drbert = "/home/suhanghua/MoSE/Datasets/preTrainModels/DR-BERT"
        tokenizer_drbert = AutoTokenizer.from_pretrained(model_name_drbert)
        model_drbert = AutoModel.from_pretrained(model_name_drbert).to(device).eval()
        gc.collect()

    def get_hidden_states(encoded, model):
        with torch.no_grad():
            output = model(**encoded)
        # Get last hidden state
        return output.last_hidden_state

    # if len(seq) > 1022, cut it    
    seq = seq.upper()
    if len(seq) > 1022:
        embedding_list = []
        seq_list = textwrap.wrap(seq, 1022)
        for seq in seq_list:
            encoded = tokenizer_drbert(seq, return_tensors="pt").to(device)
            embedding = get_hidden_states(encoded, model_drbert).detach().cpu().numpy()#[0, 1:-1, :embedding_dim]
            embedding_list.append(embedding[0, 1:-1])
        embedding = np.concatenate(embedding_list, axis=0)
    else:
        encoded = tokenizer_drbert(seq, return_tensors="pt").to(device)
        embedding = get_hidden_states(encoded, model_drbert).detach().cpu().numpy()[0, 1:-1]#[0, 1:-1, :embedding_dim]

    # return embedding[:maxsize]
    return embedding

def encode_by_esm2(seq):
    global model_esm2, alphabet_esm2
    # Load ESM-2 model if not already loaded
    if model_esm2 is None or alphabet_esm2 is None:
        print("loading ESM-2...")
        model_esm2, alphabet_esm2 = torch.hub.load(ESM2_PATH, "esm2_t33_650M_UR50D", source="local")
        gc.collect()
    
    batch_converter = alphabet_esm2.get_batch_converter()
    model_esm2.eval()  # disables dropout for deterministic results
    # Prepare data (first 2 sequences from ESMStructuralSplitDataset superfamily / 4)
    data = [
        ("tmpDeepDRP", seq.upper()),
    ]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    # if len(seq) > 4000:#len(seq)>2000: len(seq)>3088: 
    #     cut_num = 2000 # 1000
    # else:
    #     cut_num = 4000 # 2000
    # # if seq len > 2000, cut it
    # if len(seq) > cut_num:
    #     # 定义切分间隔
    #     interval = cut_num
    #     # 切分张量
    #     tensors = [batch_tokens[:, i:i + interval] for i in range(0, batch_tokens.size(1), interval)]
    #     results = []
    #     for tensor in tensors:
    #         with torch.no_grad():
    #             result = model_esm2(tensor, repr_layers=[33], return_contacts=True)
    #             gc.collect()
    #         token_representation = result["representations"][33] 
    #         results.append(token_representation)
    #     token_representations = torch.cat(results, dim=1)
    #     token_representations = token_representations.detach().cpu().numpy()[0][1:-1,:]
    # else:
    #     with torch.no_grad():
    #         results = model_esm2(batch_tokens, repr_layers=[33], return_contacts=True)
    #     token_representations = results["representations"][33]
    #     token_representations = token_representations.detach().cpu().numpy()[0][1:-1,:]
    
    with torch.no_grad():
        results = model_esm2(batch_tokens, repr_layers=[33], return_contacts=True)
    token_representations = results["representations"][33]
    token_representations = token_representations.detach().cpu().numpy()[0][1:-1,:]

    return token_representations


class InputFastaExample(object):
    """A single example feature of input data"""    
    def __init__(self, id, name, seq, feature1, feature2, feature3, label):
        self.id = id
        self.name = name
        self.seq = seq
        self.feature1 = feature1
        self.feature2 = feature2
        self.feature3 = feature3
        self.label = label

def convert_example_to_features(args, examples, data_type='train', max_esm2_dim=None, max_drbert_dim=None):
    """Loads a data file into a list of InputFastaExample."""
    datasets = []
    skipped_samples = []
    if args.debug:
        if data_type in ['SL329', 'DISORDER723_test', 'CASP']:
            examples = examples
        else:
            examples = examples[:50]
    if data_type in ["Train_723_bc25", "Train_494_bc25", "Train_casp_bc25", "Train_329_bc25", "Train_504_bc25", "Train_CAID3_nox_bc25", "Train_CAID3_pdb_bc25", "Train_CAID2_nox_bc25", "Train_CAID2_pdb_bc25", "DM4845_nr25_from_CAID3"]:
        data_type = "DM4845_training"
    elif data_type in ["DM4845_training_nr25_SL329", "DM4845_training_nr25_MXD494", "DM4845_training_nr25_DISORDER723_test", "DM4845_training_nr25_CAID2_disorder_pdb", "DM4845_training_nr25_CAID3_disorder_pdb"]:
        data_type = "DM4845_training"
    elif data_type == "lipid_TR552_processed_merged":
        data_type = "lipid_TR552_processed_merged"
    elif data_type == "TR552_idr_binding":
        data_type = "lipid_TR552_processed_merged"
    elif data_type == "TR552":
        data_type = "TR552_processed"
    
    # 默认截断维度
    if max_esm2_dim is None:
        max_esm2_dim = 640
    if max_drbert_dim is None:
        max_drbert_dim = 384
    for i, example in tqdm(enumerate(examples), desc="Example:"):
        try:
            # Perform preprocessing on the sequence
            id = example.id
            seq = example.seq
            label = example.label
            features = []
            
            if args.pssm:
                pssm_feature = encode_by_pssm(args, data_type, id)
                if pssm_feature.size > 0:
                    features.append(pssm_feature)
            else:
                pssm_feature = np.zeros((len(seq), 0))
            
            # 由于计算预训练的 feature 时间过长，需要对计算过程进行单独保持
            # 按照氨基酸的 id 进行预训练特征的保存
            if args.esm2:
                # 首先尝试加载 .pt 格式的特征文件（新格式）
                preT_path_pt = os.path.join(args.preT_feature_path, data_type, "{}.pt".format(id))
                # 其次尝试 {data_type}_esm2 子目录下的 .pt 文件
                preT_path_pt_esm2 = os.path.join(args.preT_feature_path, data_type + "_esm2", "{}.pt".format(id))
                # 然后尝试加载 .npy 格式的特征文件（旧格式）
                preT_path_npy = os.path.join(args.preT_feature_path, data_type, "esm2", "{}.npy".format(id))
                
                if os.path.exists(preT_path_pt):
                    # 加载 .pt 格式特征
                    preT_esm2_feature = torch.load(preT_path_pt, map_location='cpu').numpy()
                    # 截断ESM2特征维度，减少内存使用
                    if preT_esm2_feature.shape[1] > max_esm2_dim:
                        preT_esm2_feature = preT_esm2_feature[:, :max_esm2_dim]
                elif os.path.exists(preT_path_pt_esm2):
                    # 加载 {data_type}_esm2 子目录下的 .pt 格式特征
                    preT_esm2_feature = torch.load(preT_path_pt_esm2, map_location='cpu').numpy()
                    # 截断ESM2特征维度，减少内存使用
                    if preT_esm2_feature.shape[1] > max_esm2_dim:
                        preT_esm2_feature = preT_esm2_feature[:, :max_esm2_dim]
                elif os.path.exists(preT_path_npy):
                    # 加载 .npy 格式特征
                    preT_esm2_feature = np.load(preT_path_npy)
                    # 截断ESM2特征维度，减少内存使用
                    if preT_esm2_feature.shape[1] > max_esm2_dim:
                        preT_esm2_feature = preT_esm2_feature[:, :max_esm2_dim]
                else:
                    # 如果预训练特征不存在，则动态计算
                    dir_name = os.path.dirname(preT_path_npy)
                    if not os.path.exists(dir_name):
                        os.makedirs(dir_name)
                    preT_esm2_feature = encode_by_esm2(seq)
                    np.save(preT_path_npy, preT_esm2_feature) # save as numpy datatype
                    # 截断ESM2特征维度，减少内存使用
                    if preT_esm2_feature.shape[1] > max_esm2_dim:
                        np.save(preT_path_npy, preT_esm2_feature)
                        preT_esm2_feature = preT_esm2_feature[:, :max_esm2_dim]
                if preT_esm2_feature.size > 0:
                    features.append(preT_esm2_feature)
            else:
                preT_esm2_feature = np.zeros((len(seq), 0))

            if args.drBERT:
                # 首先尝试加载 .pt 格式的特征文件（新格式）
                preT_path_pt = os.path.join(args.preT_feature_path, data_type + "_drbert", "{}.pt".format(id))
                # 然后尝试加载 .npy 格式的特征文件（旧格式）
                preT_path_npy = os.path.join(args.preT_feature_path, data_type, "drBERT", "{}.npy".format(id))
                
                if os.path.exists(preT_path_pt):
                    # 加载 .pt 格式特征
                    preT_drbert_feature = torch.load(preT_path_pt, map_location='cpu').numpy()
                    # 截断DR-BERT特征维度，减少内存使用
                    if preT_drbert_feature.shape[1] > max_drbert_dim:
                        preT_drbert_feature = preT_drbert_feature[:, :max_drbert_dim]
                elif os.path.exists(preT_path_npy):
                    # 加载 .npy 格式特征
                    preT_drbert_feature = np.load(preT_path_npy)
                    # 截断DR-BERT特征维度，减少内存使用
                    if preT_drbert_feature.shape[1] > max_drbert_dim:
                        preT_drbert_feature = preT_drbert_feature[:, :max_drbert_dim]
                else:
                    # 如果预训练特征不存在，则动态计算
                    dir_name = os.path.dirname(preT_path_npy)
                    if not os.path.exists(dir_name):
                        os.makedirs(dir_name)
                    preT_drbert_feature = encode_by_drbert(seq)
                    np.save(preT_path_npy, preT_drbert_feature) # save as numpy datatype
                    # 截断DR-BERT特征维度，减少内存使用
                    if preT_drbert_feature.shape[1] > max_drbert_dim:
                        np.save(preT_path_npy, preT_drbert_feature)
                        preT_drbert_feature = preT_drbert_feature[:, :max_drbert_dim]
                if preT_drbert_feature.size > 0:
                    features.append(preT_drbert_feature)
            else:
                preT_drbert_feature = np.zeros((len(seq), 0))
            
            # 确保至少有一个特征
            if len(features) == 0:
                # 如果没有任何特征，使用 aaindex 特征作为默认
                aaindex_feature = encode_by_aaindex(seq)
                features.append(aaindex_feature)
            
            features = np.concatenate(features, axis=-1)
            datasets.append(InputFastaExample(
                id = id,
                name = example.name,
                seq = seq,
                feature1 = pssm_feature,
                feature2 = preT_esm2_feature,
                feature3 = preT_drbert_feature,
                label = label
            ))
            
            # 及时释放内存
            del features, pssm_feature, preT_esm2_feature, preT_drbert_feature
            if i % 10 == 0:
                gc.collect()
                
        except Exception as e:
            print(f"\nWarning: 跳过样本 {example.id} (第{i+1}个): {e}")
            skipped_samples.append(example.id)
            continue
    
    if skipped_samples:
        print(f"\n总共跳过了 {len(skipped_samples)} 个样本: {skipped_samples}")
    
    return datasets    


class CAID3_InputFastaExample(object):
    """A single example feature of input data"""    
    def __init__(self, id, name, seq, feature1, feature2, feature3, label, unexcluded_label_idx):
        self.id = id
        self.name = name
        self.seq = seq
        self.feature1 = feature1
        self.feature2 = feature2
        self.feature3 = feature3
        self.label = label
        self.unexcluded_label_idx = unexcluded_label_idx
        
def ensure_numpy_array(obj):
    if not isinstance(obj, np.ndarray):
        obj = np.array(obj)
    return obj

def CAID3_convert_example_to_features(args, examples, data_type='train', max_esm2_dim=None, max_drbert_dim=None):
    """Loads a data file into a list of InputFastaExample."""
    datasets = []
    skipped_samples = []
    if args.debug:
        if data_type in ['SL329', 'DISORDER723_test', 'CASP']:
            examples = examples
        else:
            examples = examples[:50]
    # 根据数据集类型选择特征路径
    if data_type == 'SL329':
        preT_base = "/home/suhanghua/MoSE/Data/PreT_features"
    elif data_type == 'DM4845_training_nr25_SL329':
        preT_base = "/home/suhanghua/MoSE/Datasets/PreT_features"
    else:
        preT_base = args.preT_feature_path
    
    # 特征目录名映射（数据集名可能与特征文件夹名不同）
    feature_dir_map = {
        'DM4845_training_nr25_SL329': 'DM4845_training',
    }
    feature_data_type = feature_dir_map.get(data_type, data_type)
    
    # SL329不截断，其他数据集默认截断
    if max_esm2_dim is None:
        max_esm2_dim = 1280 if data_type == 'SL329' else 640
    if max_drbert_dim is None:
        max_drbert_dim = 768 if data_type == 'SL329' else 384
    
    for i, example in tqdm(enumerate(examples), desc="Example:"):
        try:
            # Perform preprocessing on the sequence
            id = example.id
            seq = example.seq
            label = example.label
            unexcluded_label_idx = example.unexcluded_label_idx
                
            if args.pssm:
                pssm_feature = encode_by_pssm(args, feature_data_type, id)
            else:
                pssm_feature = np.zeros((len(seq), 0))
            
            # 由于计算预训练的feature时间过长，需要对计算过程进行单独保持
            # 按照氨基酸的id进行预训练特征的保存
            if args.esm2:
                preT_path = os.path.join(preT_base, feature_data_type, "esm2", "{}.npy".format(id))
                if os.path.exists(preT_path):
                    preT_esm2_feature = np.load(preT_path)
                    # 截断ESM2特征维度，减少内存使用
                    if preT_esm2_feature.shape[1] > max_esm2_dim:
                        preT_esm2_feature = preT_esm2_feature[:, :max_esm2_dim]
                else:
                    dir_name = os.path.dirname(preT_path)
                    if not os.path.exists(dir_name):
                        os.makedirs(dir_name)
                    preT_esm2_feature = encode_by_esm2(seq)
                    np.save(preT_path, preT_esm2_feature) # save as numpy datatype
                    # 截断ESM2特征维度，减少内存使用
                    if preT_esm2_feature.shape[1] > max_esm2_dim:
                        np.save(preT_path, preT_esm2_feature)
                        preT_esm2_feature = preT_esm2_feature[:, :max_esm2_dim]
            else:
                preT_esm2_feature = np.zeros((len(seq), 0))

            if args.drBERT:
                preT_path = os.path.join(preT_base, feature_data_type, "drBERT", "{}.npy".format(id))
                if os.path.exists(preT_path):
                    preT_drbert_feature = np.load(preT_path)
                    # 截断DR-BERT特征维度，减少内存使用
                    if preT_drbert_feature.shape[1] > max_drbert_dim:
                        preT_drbert_feature = preT_drbert_feature[:, :max_drbert_dim]
                else:
                    dir_name = os.path.dirname(preT_path)
                    if not os.path.exists(dir_name):
                        os.makedirs(dir_name)
                    preT_drbert_feature = encode_by_drbert(seq)
                    np.save(preT_path, preT_drbert_feature) # save as numpy datatype
                    # 截断DR-BERT特征维度，减少内存使用
                    if preT_drbert_feature.shape[1] > max_drbert_dim:
                        np.save(preT_path, preT_drbert_feature)
                        preT_drbert_feature = preT_drbert_feature[:, :max_drbert_dim]
            else:
                preT_drbert_feature = np.zeros((len(seq), 0))
            
            datasets.append(CAID3_InputFastaExample(
                id = id,
                name = example.name,
                seq = seq,
                feature1 = pssm_feature,
                feature2 = preT_esm2_feature,
                feature3 = preT_drbert_feature,
                label = label,
                unexcluded_label_idx = unexcluded_label_idx
            ))
            
            # 及时释放内存
            del pssm_feature, preT_esm2_feature, preT_drbert_feature
            if i % 10 == 0:
                gc.collect()
                
        except Exception as e:
            print(f"\nWarning: 跳过样本 {example.id} (第{i+1}个): {e}")
            skipped_samples.append(example.id)
            continue
    
    if skipped_samples:
        print(f"\n总共跳过了 {len(skipped_samples)} 个样本: {skipped_samples}")
    
    return datasets 

def load_dataset(args, fn, dataset_type='train', label=1, max_esm2_dim=None, max_drbert_dim=None):
    '''
    Author: Sicen Liu
    Description: process fasta file
    param {*} args : args
    param {*} fn : source file path of fasta format with label
    param {*} label : default label for disordered residues
    param {*} max_esm2_dim : ESM2特征维度，None表示不截断（使用完整1280维）
    param {*} max_drbert_dim : DR-BERT特征维度，None表示不截断（使用完整768维）
    return {
        ids: name list
        seqs: residue list
        labels: label list
    }
    '''
    # 从args获取特征维度，如果没有设置则使用默认值
    if max_esm2_dim is None:
        max_esm2_dim = getattr(args, 'max_esm2_dim', 1280)
    if max_drbert_dim is None:
        max_drbert_dim = getattr(args, 'max_drbert_dim', 768)
    
    if dataset_type in ['CAID3_disorder_pdb', 'CAID3_disorder_nox','CASP','SL329','CAID2_disorder_pdb', 'CAID2_disorder_nox']:
        prot_list = parseCAID3_Fasta_with_label(fn)
        dataset = CAID3_convert_example_to_features(args, prot_list, data_type=dataset_type, max_esm2_dim=max_esm2_dim, max_drbert_dim=max_drbert_dim)
    elif dataset_type == 'binding_idr':
        # binding_idr使用CAID3格式解析（标签包含-和2需要mask）
        prot_list = parseCAID3_Fasta_with_label(fn)
        dataset = CAID3_convert_example_to_features(args, prot_list, data_type=dataset_type, max_esm2_dim=max_esm2_dim, max_drbert_dim=max_drbert_dim)
    else:
        prot_list = parseFasta_with_label(fn)
        # 使用指定的特征维度截断参数
        dataset = convert_example_to_features(args, prot_list, data_type=dataset_type, max_esm2_dim=max_esm2_dim, max_drbert_dim=max_drbert_dim)
    return dataset

class proteinDataset(Dataset):
    def __init__(self, samples):
        self.samples = samples

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, index):
        return self.samples[index]

    def get_id(self, index):
        return self.samples[index].id
    def get_seq(self, index):
        return self.samples[index].seq
    
    
def kfold_split(dataset, k=5, seed=42):
    train_datasets = []
    val_datasets = []
    kf = KFold(n_splits=k, shuffle=True, random_state=seed)
    for train_index, test_index in kf.split(dataset):
        train_data = [dataset[i] for i in train_index]
        val_data = [dataset[i] for i in test_index]
        train_datasets.append(train_data)
        val_datasets.append(val_data)
    return train_datasets, val_datasets
        
        
def collate_fn_batch(batchs, args):
    ids = []
    seqs = []
    ori_labels = []
    labels = []
    features1 = []
    features2 = []
    batch_size = len(batchs)
    max_seq_len = max([len(seq.label) for seq in batchs])
    labels = torch.full((batch_size, max_seq_len), 0.)
    input_mask = torch.full((batch_size, max_seq_len), 0.)

    feature_dim1 = batchs[0].feature1.shape[-1]
    feature_dim2 = batchs[0].feature2.shape[-1]
    feature_dim3 = batchs[0].feature3.shape[-1]
    input_feature1 = torch.full((batch_size, max_seq_len, feature_dim1), 0.)
    input_feature2 = torch.full((batch_size, max_seq_len, feature_dim2), 0.)
    input_feature3 = torch.full((batch_size, max_seq_len, feature_dim3), 0.)

    for b, batch in enumerate(batchs):
        id = batch.id
        name = batch.name
        seq = batch.seq
        feature1 = batch.feature1
        feature2 = batch.feature2
        feature3 = batch.feature3
        label = batch.label
        for s, (fea1, fea2, fea3) in enumerate(zip(feature1, feature2, feature3)):
            input_feature1[b, s, :] = torch.from_numpy(fea1)
            input_feature2[b, s, :] = torch.from_numpy(fea2)
            input_feature3[b, s, :] = torch.from_numpy(fea3)
        
        for s in range(len(label)):
            input_mask[b, s] = 1.
            # 标签处理：1→1.0, 0→0.0, 2→mask(不参与计算)
            if label[s] == '1':
                labels[b, s] = 1.
            elif label[s] == '0':
                labels[b, s] = 0.
            elif label[s] == '2':
                # 标签为 2，设置为 mask，不参与 loss 计算
                input_mask[b, s] = 0.
                labels[b, s] = 0.  # 任意值，反正不会被使用
            else:
                labels[b, s] = 0.0
            
        ids.append(id)
        seqs.append(id)
        ori_labels.append(label)       
       
        # 及时释放内存
        del feature1, feature2, feature3, label
    
    return  id, seq, ori_labels, input_feature1, input_feature2, input_feature3, input_mask, labels

def CAID3_collate_fn_batch(batchs, args):
    ids = []
    seqs = []
    ori_labels = []
    labels = []
    batch_size = len(batchs)
    
    max_seq_len = max([len(seq.label) for seq in batchs])
    feature_dim1 = batchs[0].feature1.shape[-1]
    feature_dim2 = batchs[0].feature2.shape[-1]
    feature_dim3 = batchs[0].feature3.shape[-1]
    input_feature1 = torch.full((batch_size, max_seq_len, feature_dim1), 0.)
    input_feature2 = torch.full((batch_size, max_seq_len, feature_dim2), 0.)
    input_feature3 = torch.full((batch_size, max_seq_len, feature_dim3), 0.)
    
    labels = torch.full((batch_size, max_seq_len), 0.)
    input_mask = torch.full((batch_size, max_seq_len), 0.)
    for b, batch in enumerate(batchs):
        id = batch.id
        name = batch.name
        seq = batch.seq
        feature1 = batch.feature1
        feature2 = batch.feature2
        feature3 = batch.feature3
        label = batch.label
        unexcluded_label_idx = batch.unexcluded_label_idx
        for s, (fea1, fea2, fea3) in enumerate(zip(feature1, feature2, feature3)):
            input_feature1[b, s, :] = torch.from_numpy(fea1)
            input_feature2[b, s, :] = torch.from_numpy(fea2)
            input_feature3[b, s, :] = torch.from_numpy(fea3)
        
        for s in range(len(label)):
            if s in unexcluded_label_idx:
                input_mask[b, s] = 1.
            # 标签处理：1→1.0, 0→0.0, 2→mask(不参与计算)
            if label[s] == '1':
                labels[b, s] = 1.
            elif label[s] == '0':
                labels[b, s] = 0.
            elif label[s] == '2':
                # 标签为 2，设置为 mask，不参与 loss 计算
                input_mask[b, s] = 0.
                labels[b, s] = 0.  # 任意值，反正不会被使用
            else:
                labels[b, s] = 0.
            
        ids.append(id)
        seqs.append(id)
        ori_labels.append(label)        
       
    return  id, seq, ori_labels, input_feature1, input_feature2, input_feature3, input_mask, labels
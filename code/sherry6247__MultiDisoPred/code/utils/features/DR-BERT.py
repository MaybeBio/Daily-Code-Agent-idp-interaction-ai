import numpy as np
from transformers import AutoModel, AutoTokenizer
import torch
import os

# 使用环境变量指定的 GPU（如果没有设置则使用 cuda:0）
gpu_id = os.environ.get('CUDA_VISIBLE_DEVICES', '0')
if ',' in gpu_id:
    # 如果指定了多个 GPU，使用第一个
    gpu_id = gpu_id.split(',')[0]
device = f"cuda:{gpu_id}" if torch.cuda.is_available() else "cpu"

# 使用绝对路径
model_name = "/home/suhanghua/MoSE/Datasets/preTrainModels/DR-BERT"

tokenizer_drbert = AutoTokenizer.from_pretrained(model_name)
model_drbert = AutoModel.from_pretrained(model_name).eval().to(device)


def encode_by_drbert(seq, use_sliding_window=True, window_size=400, stride=200):
    """
    提取 DR-BERT 特征
    
    Args:
        seq: 蛋白质序列
        use_sliding_window: 是否使用滑动窗口处理长序列
        window_size: 窗口大小（默认 400，留有余量给特殊 token）
        stride: 滑动步长（默认 200）
    
    Returns:
        特征矩阵 (seq_len, 600)
    """
    def get_hidden_states(encoded, model):
        with torch.no_grad():
            output = model(**encoded)
        # Get last hidden state
        return output.last_hidden_state

    embedding_dim = 600
    max_input_length = 510  # DR-BERT 最大输入长度（留有余量）
    
    seq = seq.upper()
    
    # 如果序列较短，直接处理
    if len(seq) <= window_size:
        encoded = tokenizer_drbert.encode_plus(seq, return_tensors="pt").to(device)
        embedding = get_hidden_states(encoded, model_drbert).detach().cpu().numpy()[0, 1:-1, :embedding_dim]
        return embedding
    
    # 如果启用滑动窗口，处理长序列
    if use_sliding_window:
        chunks = []
        positions = []
        
        # 滑动窗口分割
        start = 0
        while start < len(seq):
            end = min(start + window_size, len(seq))
            chunk = seq[start:end]
            chunks.append(chunk)
            positions.append((start, end))
            start += stride
            
            # 如果已经到达序列末尾，退出
            if end >= len(seq):
                break
        
        # 对每个窗口提取特征
        all_embeddings = []
        for i, (chunk, (start_pos, end_pos)) in enumerate(zip(chunks, positions)):
            encoded = tokenizer_drbert.encode_plus(chunk, return_tensors="pt").to(device)
            embedding = get_hidden_states(encoded, model_drbert).detach().cpu().numpy()[0, 1:-1, :embedding_dim]
            all_embeddings.append((embedding, start_pos, end_pos))
        
        # 合并重叠区域（取平均）
        full_length = len(seq)
        merged = np.zeros((full_length, embedding_dim))
        counts = np.zeros(full_length)
        
        for embedding, start_pos, end_pos in all_embeddings:
            actual_len = min(len(embedding), end_pos - start_pos)
            merged[start_pos:start_pos + actual_len] += embedding[:actual_len]
            counts[start_pos:start_pos + actual_len] += 1
        
        # 归一化（除以重叠次数）
        merged[counts > 0] /= counts[counts > 0][:, np.newaxis]
        
        return merged
    
    # 不使用滑动窗口，直接截断
    encoded = tokenizer_drbert.encode_plus(seq[:max_input_length], return_tensors="pt").to(device)
    embedding = get_hidden_states(encoded, model_drbert).detach().cpu().numpy()[0, 1:-1, :embedding_dim]
    return embedding

"""
seq = 'arck' * 150
tmpDeepDRP = encode_by_drbert(seq)
print(tmpDeepDRP.shape)
"""
seq = 'arck' * 150
tmpDeepDRP = encode_by_drbert(seq)
print(tmpDeepDRP.shape)

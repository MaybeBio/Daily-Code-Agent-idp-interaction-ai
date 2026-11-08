import torch
import os

# 使用环境变量指定的 GPU
gpu_id = os.environ.get('CUDA_VISIBLE_DEVICES', '0')
if ',' in gpu_id:
    # 如果指定了多个 GPU，使用第一个
    gpu_id = gpu_id.split(',')[0]

# 在加载模型前设置 CUDA 设备
if torch.cuda.is_available():
    with torch.cuda.device(int(gpu_id)):
        model_esm2, alphabet_esm2 = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")
else:
    model_esm2, alphabet_esm2 = torch.hub.load("facebookresearch/esm:main", "esm2_t33_650M_UR50D")


def encode_by_esm2(seq):
    batch_converter = alphabet_esm2.get_batch_converter()
    model_esm2.eval()  # disables dropout for deterministic results
    # Prepare data (first 2 sequences from ESMStructuralSplitDataset superfamily / 4)
    data = [
        ("tmpDeepDRP", seq.upper()),
    ]
    batch_labels, batch_strs, batch_tokens = batch_converter(data)

    with torch.no_grad():
        results = model_esm2(batch_tokens, repr_layers=[33], return_contacts=True)
    token_representations = results["representations"][33]
    token_representations = token_representations.detach().cpu().numpy()[0][1:-1,:]

    embedding = 600
    maxsize = 500

    return token_representations[:maxsize, :embedding]


"""
seq = "arcKM" * 120
x = encode_by_esm2(seq)
print(x.shape)
"""

seq = "arcKM" * 120
x = encode_by_esm2(seq)
print(x.shape)
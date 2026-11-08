'''
Author: Sicen Liu
Date: 2025-05-16 14:37:10
LastEditTime: 2025-06-19 15:45:36
FilePath: /liusicen/shared_files/MOFE/src/feature_level/model.py
Description: model part

Copyright (c) 2025 by ${liusicen_cs@outlook.com}, All Rights Reserved. 
'''

import torch 
import torch.nn as nn

import torch.nn.functional as F
from collections import OrderedDict


from src.feature_level.modules import TransformerEncoder
class Encoder(nn.Module):
    def __init__(self,feature1_dim, feature2_dim, hidden_dim, heads, n_layers=2, dropout=0.1):
        super(Encoder, self).__init__()
        self.feature1_dim = feature1_dim
        self.feature2_dim = feature2_dim
        self.hidden_dim = hidden_dim
        self.trans = nn.Linear(feature1_dim+feature2_dim, hidden_dim)
        self.encoder = TransformerEncoder(hidden_dim, heads=heads, depth=n_layers, dropout=dropout)
        self.output_layer = nn.Linear(hidden_dim, 1)
        
        
    def forward(self, evo_input, preT_input):
        if evo_input is not None:
            inputs = evo_input
        elif preT_input is not None:
            inputs = preT_input
        elif evo_input is not None and preT_input is not None:
            inputs = torch.cat([evo_input, preT_input], dim=-1)

        batch_size, seq_len, _ = inputs.size()

        inputs = self.trans(inputs)
        output = self.encoder(inputs)
        output = self.output_layer(output)
        output = output.squeeze(-1)
        # output = torch.sigmoid(output)
        return output     
     
def DiffSoftmax(logits, tau=1.0, hard=False, dim=-1):
    # y_soft = (logits / tau).softmax(dim)
    gumbel_noise = -torch.empty_like(logits).exponential_().log()
    noisy_logits = (logits + gumbel_noise) / tau
    y_soft = nn.functional.softmax(noisy_logits, dim=-1)
    if hard:
        # Straight through.
        index = y_soft.max(dim, keepdim=True)[1]
        y_hard = torch.zeros_like(logits, memory_format=torch.legacy_contiguous_format).scatter_(dim, index, 1.0)
        ret = y_hard - y_soft.detach() + y_soft
    else:
        # Reparametrization trick.
        ret = y_soft
    return ret



class BiLSTM(nn.Module):
    def __init__(self, feature1_dim, feature2_dim, hidden_dim, n_layers=1, dropout=0.1):
        super(BiLSTM, self).__init__()
        self.feature1_dim = feature1_dim
        self.feature2_dim = feature2_dim
        self.hidden_dim = hidden_dim
        self.bilstm = nn.LSTM(
            input_size=feature1_dim + feature2_dim,
            hidden_size=hidden_dim,
            num_layers=n_layers,
            batch_first=True,
            bidirectional=True,
            dropout=dropout if n_layers > 1 else 0
        )
        self.output_layer = nn.Linear(hidden_dim * 2, hidden_dim)
        
    def forward(self, input1, input2):
        if input1 is not None and input2 is not None:
            inputs = torch.cat([input1, input2], dim=-1)
        elif input1 is not None:
            inputs = input1
        elif input2 is not None:
            inputs = input2
        
        batch_size, seq_len, _ = inputs.size()

        output, _ = self.bilstm(inputs)
        output = self.output_layer(output)
        output = output.squeeze(-1)
        return output


class SkipGateESANetBiLSTM(nn.Module):
    def __init__(self,
                 num_classes=1,
                 feature1_dim=1280,  
                 feature2_dim=20,   
                 hidden_dim=512,     # BiLSTM隐藏层维度
                 fusion_dim=512,      # 融合后维度
                 temp=1,
                 block_rule=[1,1,1,1]):
        
        super(SkipGateESANetBiLSTM, self).__init__()
        
        self.block_rule = block_rule if block_rule else [1, 1, 1, 1]
        self.temp = temp
        
        # 特征处理层,统一维度
        self.feature1_proj = nn.Linear(feature1_dim, hidden_dim)
        self.feature2_proj = nn.Linear(feature2_dim, hidden_dim)
        
        # 编码器 进一步编码捕捉上下文信息
        self.encoder_feature1 = nn.LSTM(
            input_size=hidden_dim,
            hidden_size=hidden_dim // 2,
            num_layers=1,
            batch_first=True,
            bidirectional=True
        )
        
        self.encoder_feature2 = nn.LSTM(
            input_size=hidden_dim,
            hidden_size=hidden_dim // 2,
            num_layers=1,
            batch_first=True,
            bidirectional=True
        )
        
        # 4个block中的融合模块
        self.fuse_layer1 = BiLSTM(hidden_dim, hidden_dim, fusion_dim)
        self.fuse_layer2 = BiLSTM(hidden_dim, hidden_dim, fusion_dim)
        self.fuse_layer3 = BiLSTM(hidden_dim, hidden_dim, fusion_dim)
        self.fuse_layer4 = BiLSTM(hidden_dim, hidden_dim, fusion_dim)
        
        # 门控网络
        self.gate_layer = GlobalGateSeq(branch_num=5, input_dim=hidden_dim * 2)
        
        # 输出层
        self.output_layer = nn.Sequential(
            nn.Linear(fusion_dim, fusion_dim // 2),
            nn.ReLU(),
            nn.Linear(fusion_dim // 2, num_classes)
        )
        
        # 初始化
        self.baseline = False
        self.ini_stage = False
        self.hard_gate = False
        self.save_weight_info = False
        self.weight_list = torch.Tensor()

    def forward(self, feature1, feature2, test=False, return_weight=False):
        # 特征投影统一维度
        feature1 = self.feature1_proj(feature1)
        feature2 = self.feature2_proj(feature2)
        
        # 初始编码
        feature1, _ = self.encoder_feature1(feature1)
        feature2, _ = self.encoder_feature2(feature2)
        
        # 门控权重生成 1 基线模式 2 初始化阶段随机生成 3 默认情况使用序列的均值作为门控输入
        if self.baseline:
            weight = torch.zeros(feature1.size(0), 5).to(feature1.device)
            weight[:, 4] = 1
        elif self.ini_stage:
            weight = torch.zeros(feature1.size(0), 5).to(feature1.device)
            idx = torch.randint(0, 5, (feature1.size(0),))
            weight[range(feature1.size(0)), idx] = 1
        else:
            gate_input = torch.mean(torch.cat([feature1, feature2], dim=-1), dim=1)
            weight = self.gate_layer(gate_input, self.temp, self.hard_gate)
        
        if self.save_weight_info:
            self.weight_list = torch.cat((self.weight_list, weight.cpu()))
        

        # Block 1
        feature1 = feature1
        fused1 = self.fuse_layer1(feature1, feature2)
        w = weight[:, 0].view(-1, 1, 1)
        output1 = w * feature1 + (1 - w) * fused1
        
        # Block 2
        feature1 = output1
        fused2 = self.fuse_layer2(feature1, feature2)
        w = (weight[:, 0] + weight[:, 1]).view(-1, 1, 1)
        output2 = w * feature1 + (1 - w) * fused2
        
        # Block 3
        feature1 = output2
        fused3 = self.fuse_layer3(feature1, feature2)
        w = (weight[:, 0] + weight[:, 1] + weight[:, 2]).view(-1, 1, 1)
        output3 = w * feature1 + (1 - w) * fused3
        
        # Block 4
        feature1 = output3
        fused4 = self.fuse_layer4(feature1, feature2)
        w = weight[:, 4].view(-1, 1, 1)
        output4 = (1 - w) * feature1 + w * fused4
        
        # 最终输出
        out = self.output_layer(output4)
        out = out.squeeze(-1)
        
        if test and return_weight:
            return out, weight
        else:
            return out





#门控网络，适用于序列数据
class GlobalGateSeq(nn.Module):
    def __init__(self, branch_num, input_dim=512, hidden_dim=64):
        super(GlobalGateSeq, self).__init__()
        self.bnum = branch_num
        self.fc = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, branch_num),
            nn.Softmax(dim=1)
        )
        
    def forward(self, x, temp=1.0, hard=False):
        # x: [batch_size, input_dim]
        logits = self.fc(x)
        if hard:
            # 硬选择
            indices = torch.argmax(logits, dim=1)
            weights = F.one_hot(indices, num_classes=self.bnum).float()
        else:
            # Gumbel softmax采样
            weights = F.gumbel_softmax(logits, tau=temp, hard=hard)
        return weights
    
def main():
    model = SkipGateESANetBiLSTM(
        num_classes=1,         # 输出类别数
        feature1_dim=1280,      # ESM-2特征维度
        feature2_dim=20,       # PSSM特征维度
        hidden_dim=512,         # BiLSTM隐藏层维度
        fusion_dim=512          # 融合后维度
    )

# 假设输入特征
    esm2_features = torch.randn(16, 500, 1280)  # [batch_size, seq_len, feature_dim]
    pssm_features = torch.randn(16, 500, 20)

# 前向传播
    output, _ = model(esm2_features, pssm_features)
    print(output.shape)  # torch.Size([16, 500, 20])

if __name__ == "__main__":
    main()

"""     # 初始化模型
model = SkipGateESANetBiLSTM(
    num_classes=20,         # 输出类别数
    feature1_dim=1280,      # ESM-2特征维度
    feature2_dim=20,       # PSSM特征维度
    hidden_dim=512,         # BiLSTM隐藏层维度
    fusion_dim=256          # 融合后维度
)

# 假设输入特征
esm2_features = torch.randn(16, 500, 1280)  # [batch_size, seq_len, feature_dim]
pssm_features = torch.randn(16, 500, 20)

# 前向传播
output, loss = model(esm2_features, pssm_features)
print(output.shape)  # torch.Size([16, 500, 20]) """
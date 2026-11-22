import os
import json
import math
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader

# ============================================================
# CONFIG
# ============================================================

SEQ_LEN = 501
CHUNK_SIZE = 25
STRIDE = 15

DIM = 512
BATCH_SIZE = 128
EPOCHS = 2000
LR = 5e-4

DEVICE = "cuda:6" if torch.cuda.is_available() else "cpu"

SAVE_DIR = "hierarchical_protein_model_cs=50_stride=25"
os.makedirs(SAVE_DIR, exist_ok=True)

# ============================================================
# VOCAB
# ============================================================

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
EXTRA = ["X"]
PAD_TOKEN = "[PAD]"

vocab = AMINO_ACIDS + EXTRA + [PAD_TOKEN]

stoi = {ch: i for i, ch in enumerate(vocab)}
PAD_ID = stoi[PAD_TOKEN]
VOCAB_SIZE = len(vocab)

with open(os.path.join(SAVE_DIR, "vocab.json"), "w") as f:
    json.dump({"stoi": stoi, "PAD_ID": PAD_ID}, f)

# ============================================================
# DATASET
# ============================================================

class ProteinDataset(Dataset):
    def __init__(self, file):
        self.seqs = open(file).read().strip().split("\n")

    def encode(self, seq):
        x = [stoi.get(a, stoi["X"]) for a in seq.upper()]

        if len(x) > SEQ_LEN:
            x = x[:SEQ_LEN]
        else:
            x += [PAD_ID] * (SEQ_LEN - len(x))

        return torch.tensor(x, dtype=torch.long)

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, idx):
        return self.encode(self.seqs[idx])

# ============================================================
# CHUNKING
# ============================================================

def create_chunks(x):
    B, L = x.shape
    chunks = []

    for i in range(0, L - CHUNK_SIZE + 1, STRIDE):
        chunks.append(x[:, i:i+CHUNK_SIZE])

    return torch.stack(chunks, dim=1)

# ============================================================
# POSITIONAL ENCODING
# ============================================================

class PositionalEncoding(nn.Module):
    def __init__(self, dim, max_len=512):
        super().__init__()

        pe = torch.zeros(max_len, dim)

        pos = torch.arange(max_len).unsqueeze(1)

        div = torch.exp(torch.arange(0, dim, 2) * (-math.log(10000.0) / dim))

        pe[:, 0::2] = torch.sin(pos * div)
        pe[:, 1::2] = torch.cos(pos * div)

        self.register_buffer("pe",pe.unsqueeze(0))

    def forward(self, x):
        return x + self.pe[:, :x.size(1)]

# ============================================================
# LOCAL ENCODER (MLM + CLS)
# ============================================================

class LocalEncoder(nn.Module):

    def __init__(self):
        super().__init__()

        self.embed = nn.Embedding(VOCAB_SIZE,DIM,padding_idx=PAD_ID)

        self.pos = PositionalEncoding(DIM,max_len=CHUNK_SIZE + 1)

        self.cls = nn.Parameter(torch.randn(1, 1, DIM))

        enc = nn.TransformerEncoderLayer(d_model=DIM,nhead=4, batch_first=True)

        self.encoder = nn.TransformerEncoder(enc,num_layers=2)

        self.norm = nn.LayerNorm(DIM)

        self.lm_head = nn.Linear(DIM,VOCAB_SIZE)

    def forward(self, x):

        B, L = x.shape

        h = self.embed(x)

        cls = self.cls.expand(B, -1, -1)

        h = torch.cat([cls, h],dim=1)

        h = self.pos(h)

        pad_mask = (x == PAD_ID)

        pad_mask = torch.cat([torch.zeros(B,1,dtype=torch.bool,device=x.device),pad_mask],dim=1)

        h = self.encoder(h,src_key_padding_mask=pad_mask)

        h = self.norm(h)

        cls_token = h[:, 0]

        token_repr = h[:, 1:]

        logits = self.lm_head(token_repr)

        return cls_token, logits
    


# ============================================================
# GLOBAL ENCODER
# ============================================================

class GlobalEncoder(nn.Module):

    def __init__(self):
        super().__init__()

        self.global_cls = nn.Parameter(torch.randn(1, 1, DIM))

        self.pos = PositionalEncoding(DIM,max_len=512)

        enc = nn.TransformerEncoderLayer(d_model=DIM,nhead=4,batch_first=True)

        self.encoder = nn.TransformerEncoder(enc, num_layers=2)

        self.norm = nn.LayerNorm(DIM)

    def forward(self,chunk_embeddings,chunk_valid):

        B, N, D = chunk_embeddings.shape

        cls = self.global_cls.expand(B,-1,-1)

        h = torch.cat([cls, chunk_embeddings],dim=1)

        h = self.pos(h)

        global_pad_mask = torch.cat([torch.zeros(B,1,dtype=torch.bool,device=h.device),~chunk_valid],dim=1)

        h = self.encoder(h,src_key_padding_mask=global_pad_mask)

        h = self.norm(h)

        return h[:, 0]

# ============================================================
# FULL MODEL
# ============================================================

class HierarchicalModel(nn.Module):

    def __init__(self):
        super().__init__()

        self.local = LocalEncoder()

        self.global_enc = GlobalEncoder()

    def forward(self, chunks):

        B, N, C = chunks.shape

        cls_tokens = []

        logits_list = []

        chunk_valid = ((chunks != PAD_ID).sum(dim=2)>= 4)

        for i in range(N):

            cls_token, logits = self.local(chunks[:, i])

            cls_tokens.append(cls_token)

            logits_list.append(logits)

        cls_tokens = torch.stack(cls_tokens,dim=1)

        z = self.global_enc(cls_tokens,chunk_valid)

        return z, logits_list

# ============================================================
# AUGMENTATION
# ============================================================
def augment(x, p=0.2):

    valid = (x != PAD_ID)

    mask = (torch.rand_like(x.float()) < p) & valid

    x2 = x.clone()

    noise = torch.randint(low=0,high=PAD_ID,size=x.shape,device=x.device)

    x2[mask] = noise[mask]

    return x2, mask

# ============================================================
# LOSSES
# ============================================================

def contrastive_loss(z1, z2, temp=0.1):
    z1 = F.normalize(z1, dim=-1)
    z2 = F.normalize(z2, dim=-1)

    logits = (z1 @ z2.T) / temp

    labels = torch.arange(z1.size(0), device=z1.device)

    return (F.cross_entropy(logits, labels) +
            F.cross_entropy(logits.T, labels)) / 2
            
def lm_loss(logits_list, target_chunks):

    losses = []

    for i in range(len(logits_list)):

        targets = target_chunks[:, i]

        valid = (targets != PAD_ID)

        if valid.sum() == 0:
            continue

        ce = F.cross_entropy(logits_list[i].reshape(-1, VOCAB_SIZE),targets.reshape(-1),ignore_index=PAD_ID)

        losses.append(ce)

    if len(losses) == 0:
        return torch.tensor( 0.0,device=target_chunks.device)

    return torch.stack(losses).mean()
# ============================================================
# DATA
# ============================================================


dataset = ProteinDataset("sequences_len_le_500.txt")
loader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True)

model = HierarchicalModel().to(DEVICE)
opt = torch.optim.AdamW(model.parameters(), lr=LR)

loss_log = []

# ============================================================
# TRAIN LOOP
# ============================================================
for epoch in range(EPOCHS):

    total_loss = 0.0
    total_mlm = 0.0
    total_con = 0.0

    model.train()

    for x in loader:

        x = x.to(DEVICE)

        # --------------------------------------------------
        # Two augmented views
        # --------------------------------------------------

        x1, _ = augment(x)
        x2, _ = augment(x)

        # --------------------------------------------------
        # Original targets
        # --------------------------------------------------

        target_chunks = create_chunks(x)

        # --------------------------------------------------
        # Chunk augmented sequences
        # --------------------------------------------------

        c1 = create_chunks(x1)
        c2 = create_chunks(x2)

        # --------------------------------------------------
        # Forward
        # --------------------------------------------------

        z1, logits1 = model(c1)
        z2, logits2 = model(c2)

        # --------------------------------------------------
        # MLM losses
        # --------------------------------------------------

        loss_mlm_1 = lm_loss(logits1,target_chunks)

        loss_mlm_2 = lm_loss(logits2,target_chunks)

        loss_mlm = (loss_mlm_1 +loss_mlm_2) / 2

        # --------------------------------------------------
        # Contrastive loss
        # --------------------------------------------------

        loss_con = contrastive_loss(z1,z2)

        # --------------------------------------------------
        # Total loss
        # --------------------------------------------------

        loss = loss_mlm + 0.5 * loss_con

        # --------------------------------------------------
        # Backprop
        # --------------------------------------------------

        opt.zero_grad()

        loss.backward()

        torch.nn.utils.clip_grad_norm_(model.parameters(),1.0)

        opt.step()

        total_loss += loss.item()
        total_mlm += loss_mlm.item()
        total_con += loss_con.item()

    avg_loss = total_loss / len(loader)
    avg_mlm = total_mlm / len(loader)
    avg_con = total_con / len(loader)

    loss_log.append(avg_loss)

    print(f"Epoch {epoch+1:4d} | MLM {avg_mlm:.4f} | CON {avg_con:.4f} | TOTAL {avg_loss:.4f}", flush=True)
# ============================================================
# SAVE
# ============================================================

torch.save(model.state_dict(), os.path.join(SAVE_DIR, "model.pt"))

with open(os.path.join(SAVE_DIR, "loss.json"), "w") as f:
    json.dump(loss_log, f)

print("Training complete. Saved model + logs.")

import os
import json
import math
import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F


# ============================================================
# CONFIG
# ============================================================

SEQ_LEN = 501

CHUNK_SIZE = 25
STRIDE = 15

DIM = 512

DEVICE = "cuda:7" if torch.cuda.is_available() else "cpu"

MODEL_DIR = "hierarchical_protein_model_cs=50_stride=25"

INPUT_FILE = ("sequences_len_le_500.txt")

SAVE_DIR = ("TokenEmbeddings_seq_500")

os.makedirs(SAVE_DIR, exist_ok=True)

# Directory for individual token representations
TOKEN_SAVE_DIR = os.path.join(SAVE_DIR,"token_embeddings")

os.makedirs(TOKEN_SAVE_DIR,exist_ok=True)


# ============================================================
# LOAD VOCAB
# ============================================================

with open(os.path.join(MODEL_DIR,"vocab.json")) as f:

    vocab_data = json.load(f)


stoi = vocab_data["stoi"]

PAD_ID = vocab_data["PAD_ID"]

itos = { v: k for k, v in stoi.items()}

VOCAB_SIZE = len(stoi)

print("Vocab loaded")
print("Vocabulary size:", VOCAB_SIZE)
print("PAD ID:", PAD_ID)


# ============================================================
# CHUNKING
# ============================================================

def create_chunks(x):

    starts = list(range(0,SEQ_LEN - CHUNK_SIZE + 1,STRIDE))

    # Make sure the final region is included
    if starts[-1] != SEQ_LEN - CHUNK_SIZE:

        starts.append(SEQ_LEN - CHUNK_SIZE)

    chunks = [x[:, s:s + CHUNK_SIZE] for s in starts]

    return torch.stack(chunks,dim=1)


# ============================================================
# CHUNK START POSITIONS
# ============================================================

def get_chunk_starts():

    starts = list(range(0,SEQ_LEN - CHUNK_SIZE + 1,STRIDE))
    
    if starts[-1] != SEQ_LEN - CHUNK_SIZE:

        starts.append(SEQ_LEN - CHUNK_SIZE)

    return starts


CHUNK_STARTS = get_chunk_starts()

print("Chunk starts:",CHUNK_STARTS)

print("Number of chunks:",len(CHUNK_STARTS))


# ============================================================
# POSITIONAL ENCODING
# ============================================================

class PositionalEncoding(nn.Module):

    def __init__(self,dim,max_len=512):

        super().__init__()

        pe = torch.zeros(max_len,dim)

        pos = torch.arange(max_len).unsqueeze(1)

        div = torch.exp(torch.arange(0,dim,2)*(-math.log(10000.0)/ dim))
        pe[:, 0::2] = torch.sin(pos * div)

        pe[:, 1::2] = torch.cos(pos * div)
        

        self.register_buffer("pe",pe.unsqueeze(0))


    def forward(self, x):

        return ( x + self.pe[:, :x.size(1)])


# ============================================================
# LOCAL ENCODER
# ============================================================

class LocalEncoder(nn.Module):

    def __init__(self):

        super().__init__()

        self.embed = nn.Embedding(VOCAB_SIZE,DIM,padding_idx=PAD_ID)

        self.pos = PositionalEncoding( DIM,CHUNK_SIZE + 1)

        self.cls = nn.Parameter(torch.randn(1,1,DIM))

        enc = nn.TransformerEncoderLayer(d_model=DIM,nhead=4,batch_first=True)

        self.encoder = nn.TransformerEncoder(enc,num_layers=2)

        self.norm = nn.LayerNorm(DIM)

        self.lm_head = nn.Linear(DIM,VOCAB_SIZE)


    def forward(self, x):

        B, L = x.shape

        # ----------------------------------------------------
        # Token embedding
        # ----------------------------------------------------

        h = self.embed(x)

        # ----------------------------------------------------
        # CLS token
        # ----------------------------------------------------

        cls = self.cls.expand(B,-1,-1)

        h = torch.cat([cls,h],dim=1)

        # ----------------------------------------------------
        # Positional encoding
        # ----------------------------------------------------

        h = self.pos(h)

        # ----------------------------------------------------
        # Padding mask
        # ----------------------------------------------------

        pad_mask = (x == PAD_ID)

        pad_mask = torch.cat([torch.zeros(B,1,dtype=torch.bool,device=x.device),pad_mask],dim=1)

        # ----------------------------------------------------
        # Transformer
        # ----------------------------------------------------

        h = self.encoder(h,src_key_padding_mask=pad_mask)

        # ----------------------------------------------------
        # Layer normalization
        # ----------------------------------------------------

        h = self.norm(h)

        # ----------------------------------------------------
        # CLS representation
        # ----------------------------------------------------

        cls_token = h[:, 0]

        # ----------------------------------------------------
        # Individual token representations
        # ----------------------------------------------------

        token_repr = h[:, 1:]

        # ----------------------------------------------------
        # Language model prediction
        # ----------------------------------------------------

        logits = self.lm_head(token_repr)

        return (cls_token,token_repr,logits)


# ============================================================
# GLOBAL ENCODER
# ============================================================

class GlobalEncoder(nn.Module):

    def __init__(self):

        super().__init__()

        self.global_cls = nn.Parameter(torch.randn(1,1,DIM))

        self.pos = PositionalEncoding( DIM,512)

        enc = nn.TransformerEncoderLayer(d_model=DIM,nhead=4,batch_first=True)

        self.encoder = nn.TransformerEncoder(enc,num_layers=2)

        self.norm = nn.LayerNorm(DIM)


    def forward(self,chunk_embeddings,chunk_valid):

        B, N, D = chunk_embeddings.shape

        # ----------------------------------------------------
        # Global CLS
        # ----------------------------------------------------

        cls = self.global_cls.expand(B,-1,-1)
        h = torch.cat([cls,chunk_embeddings],dim=1)

        # ----------------------------------------------------
        # Positional encoding
        # ----------------------------------------------------

        h = self.pos(h)

        # ----------------------------------------------------
        # Padding mask
        # ----------------------------------------------------

        global_pad_mask = torch.cat([torch.zeros(B,1,dtype=torch.bool, device=h.device),~chunk_valid],dim=1)

        # ----------------------------------------------------
        # Transformer
        # ----------------------------------------------------

        h = self.encoder(h,src_key_padding_mask=global_pad_mask)

        # ----------------------------------------------------
        # Layer normalization
        # ----------------------------------------------------

        h = self.norm(h)

        # ----------------------------------------------------
        # Global sequence representation
        # ----------------------------------------------------

        return h[:, 0]


# ============================================================
# FULL HIERARCHICAL MODEL
# ============================================================

class HierarchicalModel(nn.Module):

    def __init__(self):

        super().__init__()

        self.local = LocalEncoder()

        self.global_enc = GlobalEncoder()


    def forward(self, chunks):

        B, N, C = chunks.shape

        cls_tokens = []

        token_reprs = []

        logits_list = []

        # ----------------------------------------------------
        # Valid chunks
        # ----------------------------------------------------

        chunk_valid = ((chunks != PAD_ID).sum(dim=2)>= 4)

        # ----------------------------------------------------
        # Local encoding
        # ----------------------------------------------------

        for i in range(N):

            (cls_token,token_repr,logits) = self.local(chunks[:, i])

            cls_tokens.append( cls_token)

            token_reprs.append(token_repr)

            logits_list.append(logits)

        # ----------------------------------------------------
        # Stack CLS representations
        # ----------------------------------------------------

        cls_tokens = torch.stack(cls_tokens,dim=1)

        # Shape:
        #
        # [B, N_chunks, 512]

        # ----------------------------------------------------
        # Stack token representations
        # ----------------------------------------------------

        token_reprs = torch.stack(token_reprs, dim=1)

        # Shape:
        #
        # [B, N_chunks, 50, 512]

        # ----------------------------------------------------
        # Global encoder
        # ----------------------------------------------------

        z = self.global_enc(cls_tokens,chunk_valid)

        return (z,cls_tokens,token_reprs,logits_list)


# ============================================================
# LOAD MODEL
# ============================================================

model = HierarchicalModel().to(DEVICE)

model.load_state_dict(torch.load(os.path.join(MODEL_DIR,"model.pt"),map_location=DEVICE))

model.eval()

print("Model loaded")


# ============================================================
# ENCODE SEQUENCE
# ============================================================

def encode_sequence(seq):

    seq = seq.strip().upper()

    tokens = [stoi.get(a,stoi["X"]) for a in seq]

    length = min(len(tokens), SEQ_LEN)

    # --------------------------------------------------------
    # Truncate
    # --------------------------------------------------------

    if len(tokens) > SEQ_LEN:

        tokens = tokens[:SEQ_LEN]

    # --------------------------------------------------------
    # Padding
    #
    # Padding is ONLY used internally so that the model can
    # process sequences of different lengths.
    #
    # We DO NOT save these padded representations.
    # --------------------------------------------------------

    else:

        tokens += [PAD_ID] * (SEQ_LEN - len(tokens))

    return (torch.tensor(tokens),length)


# ============================================================
# MASK SEQUENCE
# ============================================================

def mask_sequence(x,p=0.2):

    x = x.clone()

    valid = (x != PAD_ID)

    mask = (torch.rand_like( x.float() ) < p) & valid

    # Replace masked residues with X
    x[mask] = stoi["X"]

    return x


# ============================================================
# LOAD SEQUENCES
# ============================================================

sequences = []

with open(INPUT_FILE) as f:

    for line in f:

        line = line.strip()

        if line:

            sequences.append(line)


print("Number of sequences:",len(sequences))


# ============================================================
# STORAGE
# ============================================================

# Mean chunk-level CLS representations
local_embeddings = []

# Global sequence representations
global_embeddings = []

# Sequence lengths
sequence_lengths = []

# Reconstruction comparison
comparison = []

overall_correct = 0
overall_total = 0


# ============================================================
# INFERENCE
# ============================================================

for seq_idx, seq in enumerate(sequences):

    # --------------------------------------------------------
    # Encode sequence
    # --------------------------------------------------------

    x, length = encode_sequence(seq)

    sequence_lengths.append(length)

    x = x.unsqueeze( 0).to(DEVICE)

    # --------------------------------------------------------
    # Mask sequence
    # --------------------------------------------------------

    x = mask_sequence(x,p=0.2)

    # --------------------------------------------------------
    # Create overlapping chunks
    # --------------------------------------------------------

    chunks = create_chunks(x)

    # --------------------------------------------------------
    # Model inference
    # --------------------------------------------------------

    with torch.no_grad():

        (z,cls_tokens,token_reprs,logits_list) = model(chunks)


    # ========================================================
    # GLOBAL REPRESENTATION
    # ========================================================

    global_emb = (z.squeeze(0).cpu().numpy())

    global_embeddings.append(global_emb.tolist())


    # ========================================================
    # LOCAL REPRESENTATION
    # ========================================================

    local_emb = (cls_tokens.mean(dim=1).squeeze(0).cpu().numpy())

    local_embeddings.append(local_emb.tolist())


    # ========================================================
    # TOKEN REPRESENTATIONS
    # ========================================================

    # Current shape:
    #
    # [1, N_chunks, 50, 512]
    #
    # Remove batch dimension.

    token_reprs = (token_reprs.squeeze(0).cpu().numpy())


    # --------------------------------------------------------
    # Accumulate representations for original residues
    # --------------------------------------------------------

    residue_repr = np.zeros((length,DIM),dtype=np.float32)

    residue_count = np.zeros(length,dtype=np.float32)


    # --------------------------------------------------------
    # Map overlapping chunks back to original sequence
    # --------------------------------------------------------

    for chunk_idx, start in enumerate(CHUNK_STARTS):

        chunk_repr = token_reprs[chunk_idx]

        for j in range(CHUNK_SIZE):

            pos = start + j

            # Don't include padded residues
            if pos >= length:
                continue

            residue_repr[pos] += ( chunk_repr[j])

            residue_count[pos] += 1.0


    # --------------------------------------------------------
    # Average overlapping representations
    # --------------------------------------------------------

    valid_positions = (residue_count > 0)

    residue_repr[valid_positions] /= (residue_count[valid_positions][:, None])


    # ========================================================
    # SAVE TOKEN REPRESENTATION
    # ========================================================

    # IMPORTANT:
    #
    # No padding is saved.
    #
    # Example:
    #
    # sequence length = 187
    #
    # saved shape = (187, 512)
    #
    # sequence length = 251
    #
    # saved shape = (251, 512)

    token_file = os.path.join(TOKEN_SAVE_DIR,f"sample_{seq_idx + 1}.npy")

    np.save(token_file,residue_repr)


    # ========================================================
    # PREDICTIONS
    # ========================================================

    pred_chunks = []

    for logits in logits_list:

        pred = logits.argmax(dim=-1)

        pred_chunks.append(pred)

    pred_chunks = torch.stack(pred_chunks,dim=1)


    # ========================================================
    # ACCURACY
    # ========================================================

    target_chunks = chunks

    correct = (pred_chunks == target_chunks)

    valid = (target_chunks != PAD_ID)

    overall_correct += (correct & valid).sum().item()

    overall_total += (valid).sum().item()


    # ========================================================
    # RECONSTRUCTION
    # ========================================================

    recon_tokens = np.full(SEQ_LEN,PAD_ID,dtype=np.int64)


    for chunk_idx, start in enumerate(CHUNK_STARTS):

        pred_chunk = (pred_chunks[0,chunk_idx].cpu().numpy())

        for j in range(CHUNK_SIZE):

            pos = start + j

            if pos >= length:
                continue

            if recon_tokens[pos] == PAD_ID:

                recon_tokens[pos] = (pred_chunk[j])


    # --------------------------------------------------------
    # Convert reconstructed tokens
    # --------------------------------------------------------

    recon_seq = "".join([itos[int(tok)] for tok in recon_tokens[:length]])



    comparison.append(f"ORIG:  {seq}")

    comparison.append(f"RECON: {recon_seq}")

    comparison.append("")


    # ========================================================
    # PROGRESS
    # ========================================================

    if (seq_idx + 1) % 100 == 0:

        print(f"Processed "
            f"{seq_idx + 1}/"
            f"{len(sequences)}" , flush=True)


# ============================================================
# RECONSTRUCTION ACCURACY
# ============================================================

accuracy = (overall_correct/overall_total)

print("\nReconstruction accuracy:",f"{accuracy:.4f}")


# ============================================================
# SAVE LOCAL EMBEDDINGS
# ============================================================

with open(os.path.join(SAVE_DIR,"local_embeddings.json"),"w") as f:

    json.dump(local_embeddings,f)


# ============================================================
# SAVE GLOBAL EMBEDDINGS
# ============================================================

with open(os.path.join(SAVE_DIR,"global_embeddings.json"),"w") as f:

    json.dump(global_embeddings,f)


# ============================================================
# SAVE SEQUENCE LENGTHS
# ============================================================

np.save(os.path.join(SAVE_DIR,"sequence_lengths.npy"),np.array(sequence_lengths,dtype=np.int32))


# ============================================================
# SAVE COMPARISON
# ============================================================

with open(os.path.join(SAVE_DIR,"comparison.txt"),"w") as f:

    f.write("\n".join(comparison))


# ============================================================
# SUMMARY
# ============================================================

print("\n" + "=" * 60)
print("DONE")
print("=" * 60)

print("Number of sequences:",len(sequences))

print("Local embeddings:",len(local_embeddings))

print("Global embeddings:",len(global_embeddings))

print("Token embeddings saved in:",TOKEN_SAVE_DIR)

print("Sequence lengths saved in:",os.path.join(SAVE_DIR,"sequence_lengths.npy"))

print("Reconstruction accuracy:",f"{accuracy:.4f}")

print("\nExample token embedding files:")

for i in range(min(5, len(sequences))):

    file = os.path.join(TOKEN_SAVE_DIR,f"sample_{i + 1}.npy")

    emb = np.load(file)

    print(f"sample_{i + 1}.npy : ", f"{emb.shape}")

print("=" * 60)

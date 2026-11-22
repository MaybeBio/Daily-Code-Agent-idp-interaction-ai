# PyTorch model class for downstream transformer encoder
import torch
import torch.nn as nn
import math


class SequenceTransformer(nn.Module):
    """
    Model class for downstream transformer encoder, to be used on foundation model embeddings.
    Inputs sequences into multi-head self-attention layers and outputs per-token embeddings.
    """
    def __init__(self, embedding_dim: int, num_layers: int = 4, 
                 num_heads: int = 8, dim_feedforward: int = 320, 
                 dropout: float = 0.1):
        """
        :param embedding_dim: dimension of each embedding vector, 
            the higher embedding_dim the larger the model
        :type embedding_dim: int
        :param num_layers: number of attention layers
        :type num_layers: int
        :param num_heads: number of heads in multi-head attention layers
        :type num_heads: int
        :param dim_feedforward: dimension of output for the feedforward layer
            after each attention
        :type dim_feedforward: int
        :param drpoout: dropout rate in positional embedding and attention layers
        :type dropout: float
        """
        super().__init__()
        self.embedding_dim = embedding_dim
        # learnable [CLS] token
        self.cls_token = nn.Parameter(torch.randn(1, 1, embedding_dim))
        # positional embedding
        self.pos_embed = PositionalEmbedding(
            embedding_dim=embedding_dim, dropout=dropout, max_tokens=6000
        )
        # transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=embedding_dim,
            nhead=num_heads,
            dim_feedforward=dim_feedforward,
            dropout=dropout,
            batch_first=True
        )
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        # store attention weights
        self.encoder_weights = []

    def forward(self, x, attention_mask, output_attentions: bool = False):
        """
        :param x: padded input embedding; 
            shape (B, L, D), B: batch size, L: maximum length, D: embedding dimension
        :type x: torch.Tensor
        :param attention_mask: padding mask; shape (B, L), 1 for real tokens, 0 for padding
        :type attention_mask: torch.Tensor
        :param output_attentions: whether to output attention weight matrices
        :type output_attentions: bool
        :return: embeddings or Tuple of embeddings and attention weight matrices
        :rtype: torch.Tensor | tuple
        """
        B, L, D = x.shape
        # # expand and prepend [CLS] token
        # cls_token = self.cls_token.expand(B, -1, -1)  # (B, 1, D)
        # x = torch.cat([cls_token, x], dim=1)          # (B, L+1, D)
        # apply positional enembedding
        x = self.pos_embed(x)   # (B, L+1, D)
        # # update attention mask to account for [CLS] token
        # cls_mask = torch.ones(B, 1, device=x.device)
        # attention_mask = torch.cat([cls_mask, attention_mask], dim=1)  # (B, L+1)
        # generate transformer key padding mask (False = keep, True = mask out)
        key_padding_mask = ~attention_mask.bool()  # invert for TransformerEncoder
        # transformer encoding
        self.transformer_encoder.use_nested_tensor = True if self.training else False
        if output_attentions:
            self.encoder_weights.clear()
            x = self.transformer_encoder(x, src_key_padding_mask=key_padding_mask)
            return x, self.encoder_weights
        x = self.transformer_encoder(x, src_key_padding_mask=key_padding_mask)
        return x


class PositionalEmbedding(nn.Module):
    """
    Model class for positional embedding to be used in the 
    """
    def __init__(self, embedding_dim: int, dropout: float = 0.1, max_tokens: int = 6000):
        """
        :param embedding_dim: dimension of each embedding vector
        :type embedding_dim: int
        :param drpoout: dropout rate in positional embedding and attention layers
        :type dropout: float
        :param max_tokens: maximum length of sequence allowed
        :type max_tokens: int
        """
        super().__init__()
        self.dropout = nn.Dropout(dropout)
        position = torch.arange(max_tokens).unsqueeze(1)   # (max_tokens, 1)
        div_term = torch.exp(
            torch.arange(0, embedding_dim, 2) * (-math.log(10000.0) / embedding_dim)
        )
        pe = torch.zeros(max_tokens, embedding_dim)        # (max_tokens, D)
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        pe = pe.unsqueeze(0)  # (1, max_tokens, D) so it can be broadcasted across batch
        self.register_buffer('pe', pe)
    
    def forward(self, x):
        """
        :param x: padded input embedding; 
            shape (B, L, D), B: batch size, L: maximum length, D: embedding dimension
        :type x: torch.Tensor
        :return: input with positional embedding added shape;
            shape (B, L, D), B: batch size, L: maximum length, D: embedding dimension
        :rtype: torch.Tensor
        """
        x = x + self.pe[:, :x.size(1), :]
        return self.dropout(x)


# Example usage
if __name__ == '__main__':
    # input config
    B = 3       # batch size
    L = 680     # sequence length
    D = 320     # embedding dimension
    SQN_LENGTHS = [573, 215, 680]
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # random input
    embeddings = torch.randn(B, L, D, device=DEVICE)
    attention_masks = torch.zeros(B, L, device=DEVICE)
    for i, l in enumerate(SQN_LENGTHS):
        attention_masks[i, :l] = 1

    # initialize model
    encoder = SequenceTransformer(
        embedding_dim=D, 
        num_layers=1, 
        num_heads=4, 
        dim_feedforward=D, 
        dropout=0.1
    )
    encoder = encoder.to(DEVICE)

    # forward pass
    outputs = encoder(embeddings, attention_masks)

    print('Model: Multi-head self attention')
    print('Batch encoder output shape:', outputs.shape)
    print('Batch size:', len(outputs))
    
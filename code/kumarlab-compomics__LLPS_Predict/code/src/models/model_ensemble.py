# PyTorch ensemble model classes
# ESM2 + SaProt
# IntermapCNN + ESM2
# IntermapCNN + ESM2 + SaProt
import torch
import torch.nn as nn
from src.models.embed_seq import get_sqn_embed_model
from src.models.embed_intermap import IntermapEmbed
from src.models.attention import SequenceTransformer


class EnsembleClassifier(nn.Module):
    """
    Model class that uses ensemble of sequence-based embeddings and structure-based embeddings for prediction.
    e.g., ESM2 + SaProt
    The sequence-based embeddings are cross attended to the structure-based embeddings.
    """
    def __init__(self, embed_model1: str, embed_model2: str,
                 embedding_dim: int, num_layers: int, num_heads: int, 
                 dim_feedforward: int, dropout: float = 0.1, 
                 local_files_only: bool = False):
        """
        :param embed_model1: name of the first embedding model
        :type embed_model1: str
        :param embed_model2: name of the second embedding model
        :type embed_model2: str
        :param embedding_dim: dimension of the embedding vector
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
        :param local_files_only: True if load model checkpoints from local,
            False if from online hub (e.g., HuggingFaceHub)
        :type local_files_only: bool
        """
        super().__init__()
        self.embed_model1 = get_sqn_embed_model(embed_model1, local_files_only=local_files_only)
        self.embed_model2 = get_sqn_embed_model(embed_model2, local_files_only=local_files_only)
        self.encoder1 = SequenceTransformer(
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            dim_feedforward=dim_feedforward,
            dropout=dropout
        )
        self.encoder2 = SequenceTransformer(
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            dim_feedforward=dim_feedforward,
            dropout=dropout
        )
        # cross attention
        self.cross_attention = nn.MultiheadAttention(
            embed_dim=embedding_dim, num_heads=num_heads, batch_first=True
        )
        # binary classifier (uses [CLS] token)
        self.classifier = nn.Linear(embedding_dim, 1)
    
    def tokenize(self, sequence_list):
        self.embed_model1.eval()
        self.embed_model2.eval()
        with torch.no_grad():
            tokens1 = self.embed_model1.tokenize(sequence_list)
            tokens2 = self.embed_model2.tokenize(sequence_list)
            return tokens1, tokens2
    
    def forward(self, sequence_list1, sequence_list2):
        # embedding model
        self.embed_model1.eval()
        self.embed_model2.eval()
        with torch.no_grad():
            embeddings1, attention_masks1 = self.embed_model1(sequence_list1)
            embeddings2, attention_masks2 = self.embed_model2(sequence_list2)
        embeddings1 = embeddings1.detach()
        embeddings2 = embeddings2.detach()
        attention_masks1 = attention_masks1.detach()
        attention_masks2 = attention_masks2.detach()
        # transformer encoder
        x1 = self.encoder1(embeddings1, attention_masks1)
        x2 = self.encoder2(embeddings2, attention_masks2)
        # cross attention
        attn_output, _ = self.cross_attention(query=x2, key=x1, value=x1)  # enc1 attends to enc2
        x1 = x1 + attn_output  # residual fusion
        # use the [CLS] token output
        cls_output = x1[:, 0, :]  # (B, D)
        # final binary prediction
        logits = self.classifier(cls_output).squeeze(-1)  # (B,)
        return logits


class EnsembleIntermapClassifier(nn.Module):
    """
    Model class that uses ensemble of intermap CNN and sequence-based embeddings for prediction.
    e.g., FINCHES + ESM2
    The intermap CNN embeddings are added to the sequence-based embeddings.
    """
    def __init__(self, embed_model: str, 
                 embedding_dim: int, num_layers: int, num_heads: int, 
                 dim_feedforward: int, dropout: float = 0.1, 
                 local_files_only: bool = False):
        """
        :param embed_model: name of the embedding model
        :type embed_model: str
        :param embedding_dim: dimension of the embedding vector
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
        :param local_files_only: True if load model checkpoints from local,
            False if from online hub (e.g., HuggingFaceHub)
        :type local_files_only: bool
        """
        super().__init__()
        self.embed_model1 = get_sqn_embed_model(embed_model, local_files_only=local_files_only)
        self.embed_model2 = IntermapEmbed(in_channels=1, output_dim=embedding_dim, dropout=dropout)
        self.encoder1 = SequenceTransformer(
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            dim_feedforward=dim_feedforward,
            dropout=dropout
        )
        # binary classifier (uses [CLS] token)
        self.classifier = nn.Linear(embedding_dim, 1)
    
    def tokenize(self, sequence_list):
        self.embed_model1.eval()
        with torch.no_grad():
            tokens1 = self.embed_model1.tokenize(sequence_list)
            return tokens1
    
    def forward(self, sequence_list1, sequence_list2):
        # embedding model
        self.embed_model1.eval()
        self.embed_model2.train()
        with torch.no_grad():
            embeddings1, attention_masks1 = self.embed_model1(sequence_list1)
        embeddings2 = self.embed_model2(sequence_list2)
        embeddings1 = embeddings1.detach()
        # embeddings2 = embeddings2.detach()
        attention_masks1 = attention_masks1.detach()
        # transformer encoder
        x1 = self.encoder1(embeddings1, attention_masks1)
        # use the [CLS] token output
        cls_output = x1[:, 0, :]  # (B, D)
        # concatenate intermap embedding
        cls_output += embeddings2
        # final binary prediction
        logits = self.classifier(cls_output).squeeze(-1)  # (B,)
        return logits


class EnsembleIntermapClassifier2(nn.Module):
    """
    Model class that uses ensemble of intermap CNN and sequence-based embeddings for prediction.
    e.g., FINCHES + ESM2 + SaProt
    The intermap CNN embeddings are added to the sequence-based embeddings.
    """
    def __init__(self, embed_model1: str, embed_model2: str,
                 embedding_dim: int, num_layers: int, num_heads: int, 
                 dim_feedforward: int, dropout: float = 0.1, 
                 local_files_only: bool = False):
        """
        :param embed_model: name of the embedding model
        :type embed_model: str
        :param embedding_dim: dimension of the embedding vector
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
        :param local_files_only: True if load model checkpoints from local,
            False if from online hub (e.g., HuggingFaceHub)
        :type local_files_only: bool
        """
        super().__init__()
        self.embed_model1 = get_sqn_embed_model(embed_model1, local_files_only=local_files_only)
        self.embed_model2 = get_sqn_embed_model(embed_model2, local_files_only=local_files_only)
        self.embed_model3 = IntermapEmbed(in_channels=1, output_dim=embedding_dim, dropout=dropout)
        self.encoder1 = SequenceTransformer(
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            dim_feedforward=dim_feedforward,
            dropout=dropout
        )
        self.encoder2 = SequenceTransformer(
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            dim_feedforward=dim_feedforward,
            dropout=dropout
        )
        # cross attention
        self.cross_attention = nn.MultiheadAttention(
            embed_dim=embedding_dim, num_heads=num_heads, batch_first=True
        )
        # binary classifier (uses [CLS] token)
        self.classifier2 = nn.Linear(2*embedding_dim, 1)
    
    def tokenize(self, sequence_list):
        self.embed_model1.eval()
        with torch.no_grad():
            tokens1 = self.embed_model1.tokenize(sequence_list)
            return tokens1
    
    def forward(self, sequence_list1, sequence_list2, sequence_list3):
        # embedding model
        self.embed_model1.eval()
        self.embed_model2.train()
        with torch.no_grad():
            embeddings1, attention_masks1 = self.embed_model1(sequence_list1)
            embeddings2, attention_masks2 = self.embed_model2(sequence_list2)
        embeddings1 = embeddings1.detach()
        embeddings2 = embeddings2.detach()
        embeddings3 = self.embed_model3(sequence_list3)
        # embeddings3 = embeddings3.detach()
        attention_masks1 = attention_masks1.detach()
        attention_masks2 = attention_masks2.detach()
        # transformer encoder
        x1 = self.encoder1(embeddings1, attention_masks1)
        x2 = self.encoder2(embeddings2, attention_masks1)
        # cross attention
        attn_output, _ = self.cross_attention(query=x2, key=x1, value=x1)  # enc1 attends to enc2
        x1 = x1 + attn_output  # residual fusion
        # use the [CLS] token output
        cls_output = x1[:, 0, :]  # (B, D)
        # concatenate intermap embedding
        cls_output = torch.cat((cls_output.detach(), embeddings3), dim=-1)
        # final binary prediction
        logits = self.classifier2(cls_output).squeeze(-1)  # (B,)
        return logits


# Example usage
if __name__ == "__main__":
    # command
    # cd LLPS_Predict
    # python -m src.models.model_ensemble

    # input config
    SQN_EMBED_MODEL1 = 'esm2_650m'
    SQN_EMBED_MODEL2 = 'saprot_650m'
    SQN_EMBED_LOAD_LOCAL = True
    SQN_EMBED_DIM = 1280
    DROPOUT = 0.6
    N_LAYER = 1
    N_HEAD = 4
    FEED_FORWARD_DIM = 1280
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # sample input
    sequence_list = ('MEVQLVQYK', 'MEVQLVQYK')
    saprot_sequence_list = ('M#EvVpQpL#VyQdYaKv', 'M#EvVpQpL#VyQdYaKv')

    # model
    classifier = EnsembleClassifier(
        SQN_EMBED_MODEL1,
        SQN_EMBED_MODEL2,
        embedding_dim=SQN_EMBED_DIM,
        num_layers=N_LAYER,
        num_heads=N_HEAD,
        dim_feedforward=FEED_FORWARD_DIM,
        dropout=DROPOUT,
        local_files_only=SQN_EMBED_LOAD_LOCAL
    )
    classifier = classifier.to(DEVICE)

    # forward pass
    x = classifier(sequence_list, saprot_sequence_list)

    print('Model: Ensemble ESM2 + SaProt')
    print('Logits shape:', x.shape)     # Expected: (2,)
    print('Logits:', x)

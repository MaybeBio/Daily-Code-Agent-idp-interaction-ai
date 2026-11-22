# PyTorch sequence-based model classes
import torch
import torch.nn as nn
from src.models.embed_seq import get_sqn_embed_model
from src.models.attention import SequenceTransformer


class TransformerClassifier(nn.Module):
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
        self.embed_model = get_sqn_embed_model(embed_model, local_files_only=local_files_only)
        self.encoder = SequenceTransformer(
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            dim_feedforward=dim_feedforward,
            dropout=dropout
        )
        # binary classifier (uses [CLS] token)
        self.classifier = nn.Linear(embedding_dim, 1)
    
    def tokenize(self, sequence_list: list):
        """
        Tokenize list of sequences.

        :param sequence_list: list of B protein sequences, maximum length L
        :type sequence_list: List[str]
        :return: tuple of the following:
            * token ids, padded with zeros; 
            shape (B, L), B: number of proteins, L: maximum length
            * binary masks indicating which tokens are padded; 
            shape (B, L), B: number of proteins, L: maximum length
        :rtype: tuple[torch.Tensor, torch.Tensor]
        """
        self.embed_model.eval()
        with torch.no_grad():
            return self.embed_model.tokenize(sequence_list)
    
    def forward(self, sequence_list: list):
        """
        Classify list of sequences.

        :param sequence_list: list of B protein sequences, maximum length L
        :type sequence_list: List[str]
        :return: tuple of the following:
            * logits; shape (B,), B: number of proteins
            * list of list of attention weights, 
            list of (L, L) matrices from first to last attention layers for each sequence
        :rtype: torch.Tensor | tuple[torch.Tensor, torch.Tensor]
        """
        # embedding model
        self.embed_model.eval()
        with torch.no_grad():
            embeddings, attention_masks = self.embed_model(sequence_list)
        embeddings = embeddings.detach()
        attention_masks = attention_masks.detach()
        # transformer encoder
        x = self.encoder(embeddings, attention_masks)            
        # use the [CLS] token output
        cls_output = x[:, 0, :]  # (B, D)
        # final binary prediction
        logits = self.classifier(cls_output).squeeze(-1)  # (B,)
        return logits


class DenseClassifier(nn.Module):
    def __init__(self, embed_model, embedding_dim, dropout=0.1, local_files_only=False):
        super().__init__()
        self.embed_model = get_sqn_embed_model(embed_model, local_files_only=local_files_only)
        # binary classifier (uses [CLS] token)
        self.classifier = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(embedding_dim, embedding_dim),
            nn.Tanh(),
            nn.Dropout(dropout),
            nn.Linear(embedding_dim, 1)
        )
    
    def forward(self, sequence_list):
        # embedding model
        self.embed_model.eval()
        with torch.no_grad():
            embeddings, attention_masks = self.embed_model(sequence_list)
        embeddings = embeddings.detach()
        attention_masks = attention_masks.detach()
        # use the [CLS] token output
        cls_output = embeddings[:, 0, :]  # (B, D)
        # final binary prediction
        logits = self.classifier(cls_output).squeeze(-1)  # (B,)
        return logits
    

class MaxPoolClassifier(nn.Module):
    def __init__(self, embed_model,
                 embedding_dim, num_layers, num_heads, dim_feedforward, dropout=0.1, 
                 local_files_only=False):
        super().__init__()
        self.embed_model = get_sqn_embed_model(embed_model, local_files_only=local_files_only)
        self.encoder = SequenceTransformer(
            embedding_dim=embedding_dim,
            num_layers=num_layers,
            num_heads=num_heads,
            dim_feedforward=dim_feedforward,
            dropout=dropout
        )
        # binary classifier (uses [CLS] token)
        self.pool = nn.AdaptiveMaxPool1d(1)
        self.classifier = nn.Sequential(
            # nn.Dropout(dropout),
            nn.Linear(embedding_dim, 48),
            nn.ReLU(),
            nn.Linear(48, 1)
        )
    
    def forward(self, sequence_list):
        # embedding model
        self.embed_model.eval()
        with torch.no_grad():
            embeddings, attention_masks = self.embed_model(sequence_list)
        embeddings = embeddings.detach()
        attention_masks = attention_masks.detach()
        # transformer encoder
        x = self.encoder(embeddings, attention_masks)
        # use max pooling over sequence tokens
        cls_output = x[:, 1:, :]  # (B, L, D)
        # max pooling
        cls_output = cls_output.permute(0, 2, 1)
        cls_output = self.pool(cls_output)
        cls_output = cls_output.permute(0, 2, 1).squeeze()  # (B, D)
        # final binary prediction
        logits = self.classifier(cls_output).squeeze(-1)  # (B,)
        return logits


# Example usage
if __name__ == '__main__':
    # command
    # cd LLPS_Predict
    # python -m src.models.model

    # SQN_EMBED_MODEL options:
    # 'esmc_300m': ESMC 300M, embedding shape (sqn_length+2, 360)
    # 'esmc_600m': ESMC 600M, embedding shape (sqn_length+2, 1152)
    # 'esm3_sm': ESM3 SM, embedding shape (sqn_length+2, 1536)
    # 'esm2_8m': ESM2 6 Layer 8M, embedding shape (sqn_length+2, 320)
    # 'esm2_35m': ESM2 12 Layer 35M, embedding shape (sqn_length+2, 480)
    # 'esm2_150m': ESM2 30 Layer 150M, embedding shape (sqn_length+2, 640)
    # 'esm2_650m': ESM2 33 Layer 650M, embedding shape (sqn_length+2, 1280)
    # 'saprot_35m': SaProt 35M, embedding shape (sqn_length+2, 480)
    # 'saprot_650m': SaProt 650M, embedding shape (sqn_length+2, 1280)
    SQN_EMBED_MODEL = 'esm2_8m'
    # SQN_EMBED_LOAD_LOCAL: whether to load model weights from local
    # True if model weights have been downloaded
    # False if load online model weights, need Internet connection
    SQN_EMBED_LOAD_LOCAL = True
    SQN_EMBED_DIM = 320
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # sequence data
    sequence_list = (
        'AGFQPQSQGMSLNDFQKQQKQAAPKPKKTLKLVSSSGIKLANATKKVGTKPAESDKKEEEKSAETKEPTKEPTKVEEPVKKEEKPVQTEEKTEEKSELPKVEDLKISESTHNTNNANVTSADALIKEQEEEVDDEVVNDMFGGKDHVSLIFMGHVDAGKSTMGGNLLYLTGSVDKRTIEKYEREAKDAGRQGWYLSWVMDTNKEERNDGKTIEVGKAYFETEKRRYTILDAPGHKMYVSEMIGGASQADVGVLVISARKGEYETGFERGGQTREHALLAKTQGVNKMVVVVNKMDDPTVNWSKERYDQCVSNVSNFLRAIGYNIKTDVVFMPVSGYSGANLKDHVDPKECPWYTGPTLLEYLDTMNHVDRHINAPFMLPIAAKMKDLGTIVEGKIESGHIKKGQSTLLMPNKTAVEIQNIYNETENEVDMAMCGEQVKLRIKGVEEEDISPGFVLTSPKNPIKSVTKFVAQIAIVELKSIIAAGFSCVMHVHTAIEEVHIVKLLHKLEKGTNRKSKKPPAFAKKGMKVIAVLETEAPVCVETYQDYPQLGRFTLRDQGTTIAIGKIVKIAE',
        'ASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRG',
        'CLRKKRKPQAEKVDVIAGSSKMKGFSSSESESSSESSSSDSEDSETEMAPKSKKKGHPGREQKKHHHHHHQQMQQAPAPVPQQPPPPPQQPPPPPPPQQQQQPPPPPPPPSMPQQAAPAMKSSPPPFIATQVPVLEPQLPGSVFDPIGHFTQPILHLPQPELPPHLPQPPEHSTPPHLNQHAVVSPPALHNALPQQPSRPSNRAAALPPKPARPPAVSPALTQTPLLPQPPMAQPPQVLLEDEEPPAPPLTSMQMQLYLQQLQKVQPPTPLLPSVKVQSQPPPPLPPPPHPSVQQQLQQQPPPPPPPQPQPPPQQQHQPPPRPVHLQPMQFSTHIQQPPPPQGQQPPHPPPGQQPPPPQPAKPQQVIQHHHSPRHHKSDPYSTGHLREAPSPLMIHSPQMSQFQSLTHQSPPQQNVQPKKQELRAASVVQPQPLVVVKEEKIHSPIIRSEPFSPSLRPEPPKHPESIKAPVHLPQRPEMKPVDVGRPVIRPPEQNAPPPGAPDKDKQKQEPKTPVAPKKDLKIKNMGSWASLVQKHPTTPSSTAKSSSDSFEQFRRAAREKEEREKALKAQAEHAEKEKERLRQERMRSREDEDALEQARRAHEEARRRQEQQQQQRQEQQQQQQQQAAAVAAAATPQAQSSQPQSMLDQQRELARKREQERRRREAMAATIDMNFQS'
    )

    # initialize model
    classifier = TransformerClassifier(
        embed_model=SQN_EMBED_MODEL, local_files_only=SQN_EMBED_LOAD_LOCAL,
        embedding_dim=SQN_EMBED_DIM,
        num_layers=1,
        num_heads=4,
        dim_feedforward=SQN_EMBED_DIM,
        dropout=0.1
    )
    classifier = classifier.to(DEVICE)

    # forward pass
    logits = classifier(sequence_list)

    print('Model: Embedding Model -> Attention -> Linear')
    print('Logits shape:', logits.shape)
    print('Logits:', logits)

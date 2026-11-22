# PyTorch style model classes for protein language foundation models
import os
import torch
import torch.nn as nn
from torch.nn.utils.rnn import pad_sequence
# ESMC & ESM3: evolutionaryscale/esm
from esm.models.esm3 import ESM3
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein, LogitsConfig, SamplingConfig
# ESM2: hugging face model facebook/esm2
from transformers import AutoTokenizer, AutoModel


# save pretrained foundation models
CHECKPOINT_DIR = 'checkpoints/foundation_models/'


class ESMC_Embed(nn.Module):
    """
    ESM Cambrian embedding model class.
    """
    def __init__(self, model: str, local_files_only: bool = False):
        """
        :param model: path to model checkpoint
        :type model: str
        :param local_files_only: True if load model checkpoints from local,
            False if from online hub (e.g., HuggingFaceHub)
        :type local_files_only: bool
        """
        super().__init__()
        if local_files_only:
            self.model = torch.load(model, weights_only=False)
        else:
            self.model = ESMC.from_pretrained(model)
            checkpoint_path = os.path.join(CHECKPOINT_DIR, model+'.pth')
            torch.save(self.model, checkpoint_path)
    
    def tokenize(self, sequence_list: list):
        """
        Tokenize list of sequences.

        :param sequence_list: list of B protein sequences, maximum length L
        :type sequence_list: List[str]
        :return: tuple of the following:
            * token ids, padded with zeros; shape (B, L), 
            B: number of proteins, L: maximum length
            * binary masks indicating which tokens are padded; shape (B, L),
            B: number of proteins, L: maximum length
        :rtype: tuple[torch.Tensor, torch.Tensor]
        """
        tokens = [
            self.model.encode(ESMProtein(sequence=sequence)).sequence
            for sequence in sequence_list
        ]
        # pad sequences of tokens to the same length
        padded_tokens = pad_sequence(tokens, batch_first=True)
        # create attention masks: 1 for real tokens, 0 for padding
        attention_masks = torch.tensor([
            [1]*token.size(0) + [0]*(padded_tokens.size(1) - token.size(0))
            for token in tokens
        ])
        attention_masks = attention_masks.to(padded_tokens.device)
        return padded_tokens, attention_masks
    
    def embed_sequence(self, sequence: str):
        """
        Embed a single sequence into vectors.

        :param sequence: a single protein sequence, length L
        :type sequence: str
        :return: sequence embedding; 
            shape (L, D), L: sequence length, D: embedding dimension
        :rtype: torch.Tensor
        """
        EMBEDDING_CONFIG = LogitsConfig(
            sequence=True, return_embeddings=True, return_hidden_states=True
        )
        protein = ESMProtein(sequence=sequence)
        protein_tensor = self.model.encode(protein)
        output = self.model.logits(protein_tensor, EMBEDDING_CONFIG)
        return output.embeddings[0]
    
    def forward(self, sequence_list: list):
        """
        Embed list of sequences into vectors.

        :param sequence_list: list of B protein sequences, maximum length L
        :type sequence_list: List[str]
        :return: tuple of the following:
            * sequence embeddings, padded with zeros; shape (B, L, D), 
            B: number of proteins, L: maximum length, D: embedding dimension
            * binary masks indicating which tokens are padded; shape (B, L),
            B: number of proteins, L: maximum length
        :rtype: tuple[torch.Tensor, torch.Tensor]
        """
        embeddings = [self.embed_sequence(seq) for seq in sequence_list]
        # pad sequences to the same length
        padded_embeddings = pad_sequence(embeddings, batch_first=True)
        # create attention masks: 1 for real tokens, 0 for padding
        attention_masks = torch.tensor([
            [1]*emb.size(0) + [0]*(padded_embeddings.size(1) - emb.size(0))
            for emb in embeddings
        ])
        attention_masks = attention_masks.to(padded_embeddings.device)
        return padded_embeddings, attention_masks


class ESM3_Embed(nn.Module):
    """
    ESM3 embedding model class.
    """
    def __init__(self, model: str, local_files_only: bool = False):
        """
        :param model: path to model checkpoint
        :type model: str
        :param local_files_only: True if load model checkpoints from local,
            False if from online hub (e.g., HuggingFaceHub)
        :type local_files_only: bool
        """
        super().__init__()
        if local_files_only:
            self.model = torch.load(model, weights_only=False)
        else:
            self.model = ESM3.from_pretrained(model)
            checkpoint_path = os.path.join(CHECKPOINT_DIR, model+'.pth')
            torch.save(self.model, checkpoint_path)
    
    def tokenize(self, sequence_list: list):
        """
        Tokenize list of sequences.

        :param sequence_list: list of B protein sequences, maximum length L
        :type sequence_list: List[str]
        :return: tuple of the following:
            * token ids, padded with zeros; shape (B, L), 
            B: number of proteins, L: maximum length
            * binary masks indicating which tokens are padded; shape (B, L),
            B: number of proteins, L: maximum length
        :rtype: tuple[torch.Tensor, torch.Tensor]
        """
        tokens = [
            self.model.encode(ESMProtein(sequence=sequence)).sequence
            for sequence in sequence_list
        ]
        # pad sequences of tokens to the same length
        padded_tokens = pad_sequence(tokens, batch_first=True)
        # create attention masks: 1 for real tokens, 0 for padding
        attention_masks = torch.tensor([
            [1]*token.size(0) + [0]*(padded_tokens.size(1) - token.size(0))
            for token in tokens
        ])
        attention_masks = attention_masks.to(padded_tokens.device)
        return padded_tokens, attention_masks
    
    def embed_sequence(self, sequence: str):
        """
        Embed a single sequence into vectors.

        :param sequence: a single protein sequence, length L
        :type sequence: str
        :return: sequence embedding; 
            shape (L, D), L: sequence length, D: embedding dimension
        :rtype: torch.Tensor
        """
        EMBEDDING_CONFIG = SamplingConfig(return_per_residue_embeddings=True)
        protein = ESMProtein(sequence=sequence)
        protein_tensor = self.model.encode(protein)
        output = self.model.forward_and_sample(protein_tensor, EMBEDDING_CONFIG)
        return output.per_residue_embedding
    
    def forward(self, sequence_list: list):
        """
        Embed list of sequences into vectors.

        :param sequence_list: list of B protein sequences, maximum length L
        :type sequence_list: List[str]
        :return: tuple of the following:
            * sequence embeddings, padded with zeros; shape (B, L, D), 
            B: number of proteins, L: maximum length, D: embedding dimension
            * binary masks indicating which tokens are padded; shape (B, L),
            B: number of proteins, L: maximum length
        :rtype: tuple[torch.Tensor, torch.Tensor]
        """
        embeddings = [self.embed_sequence(seq) for seq in sequence_list]
        # pad sequences to the same length
        padded_embeddings = pad_sequence(embeddings, batch_first=True)
        # create attention masks: 1 for real tokens, 0 for padding
        attention_masks = torch.tensor([
            [1]*emb.size(0) + [0]*(padded_embeddings.size(1) - emb.size(0))
            for emb in embeddings
        ])
        attention_masks = attention_masks.to(padded_embeddings.device)
        return padded_embeddings, attention_masks


class ESM2_Embed(nn.Module):
    """
    ESM2 embedding model class.
    """
    def __init__(self, model: str, local_files_only: bool = False):
        """
        :param model: path to model checkpoint
        :type model: str
        :param local_files_only: True if load model checkpoints from local,
            False if from online hub (e.g., HuggingFaceHub)
        :type local_files_only: bool
        """
        super().__init__()
        self.tokenizer = AutoTokenizer.from_pretrained(model, local_files_only=local_files_only)
        self.model = AutoModel.from_pretrained(model, local_files_only=local_files_only, add_pooling_layer=False)
    
    def tokenize(self, sequence_list: list):
        """
        Tokenize list of sequences.

        :param sequence_list: list of B protein sequences, maximum length L
        :type sequence_list: List[str]
        :return: tuple of the following:
            * token ids, padded with zeros; shape (B, L), 
            B: number of proteins, L: maximum length
            * binary masks indicating which tokens are padded; shape (B, L),
            B: number of proteins, L: maximum length
        :rtype: tuple[torch.Tensor, torch.Tensor]
        """
        protein_tensors = self.tokenizer(sequence_list, return_tensors="pt", padding=True)
        protein_tensors = protein_tensors.to(self.model.device)
        return protein_tensors['input_ids'], protein_tensors['attention_mask']
    
    def embed_sequence(self, sequence: str):
        """
        Embed a single sequence into vectors.

        :param sequence: a single protein sequence, length L
        :type sequence: str
        :return: tuple of the following:
            * token ids, padded with zeros; shape (L,), L: sequence length
            * binary masks indicating which tokens are padded; shape (L,), L: sequence length
        :rtype: tuple[torch.Tensor, torch.Tensor]
        """
        protein_tensors = self.tokenizer(sequence, return_tensors="pt")
        protein_tensors = protein_tensors.to(self.model.device)
        output = self.model(**protein_tensors)
        return output.last_hidden_state[0]

    def forward(self, sequence_list: list):
        """
        Embed list of sequences into vectors.

        :param sequence_list: list of B protein sequences, maximum length L
        :type sequence_list: List[str]
        :return: tuple of the following:
            * sequence embeddings, padded with zeros; shape (B, L, D), 
            B: number of proteins, L: maximum length, D: embedding dimension
            * binary masks indicating which tokens are padded; shape (B, L),
            B: number of proteins, L: maximum length
            * list of attention weights, one (L, L) matrix for each sequence
        :rtype: tuple[torch.Tensor, torch.Tensor]
        """
        # get embeddings and attention weights (optional) for all sequences
        embeddings = [self.embed_sequence(seq) for seq in sequence_list]
        # pad sequences to the same length
        padded_embeddings = pad_sequence(embeddings, batch_first=True)
        # create attention masks: 1 for real tokens, 0 for padding
        attention_masks = torch.tensor([
            [1]*emb.size(0) + [0]*(padded_embeddings.size(1) - emb.size(0))
            for emb in embeddings
        ])
        attention_masks = attention_masks.to(padded_embeddings.device)
        return padded_embeddings, attention_masks


def get_sqn_embed_model(model_name: str, local_files_only: bool = False):
    """
    Model factory function that returns model instance given model name.

    :param model_name: model name; list of supported model names:
        esmc_300m, esmc_600m, esm3_sm, esm2_8m, esm2_35m, esm2_150m, esm2_650m, 
        saprot_35m, and saprot_650m
    :type model_name: str
    :param local_files_only: True if load model checkpoints from local,
        False if from online hub (e.g., HuggingFaceHub)
    :type local_files_only: bool
    :return: model instance
    :rtype: ESMC_Embed or ESM3_Embed or ESM2_Embed
    """
    model_paths = {
        'local': {
            'esmc_300m': os.path.join(CHECKPOINT_DIR, 'esmc_300m.pth'),
            'esmc_600m': os.path.join(CHECKPOINT_DIR, 'esmc_600m.pth'),
            'esm3_sm': os.path.join(CHECKPOINT_DIR, 'esm3_sm_open_v1.pth'),
            'esm2_8m': 'facebook/esm2_t6_8M_UR50D',
            'esm2_35m': 'facebook/esm2_t12_35M_UR50D',
            'esm2_150m': 'facebook/esm2_t30_150M_UR50D',
            'esm2_650m': 'facebook/esm2_t33_650M_UR50D',
            'saprot_35m': 'westlake-repl/SaProt_35M_AF2',
            'saprot_650m': 'westlake-repl/SaProt_650M_AF2'
        },
        'online': {
            'esmc_300m': 'esmc_300m',
            'esmc_600m': 'esmc_600m',
            'esm3_sm': 'esm3_sm_open_v1',
            'esm2_8m': 'facebook/esm2_t6_8M_UR50D',
            'esm2_35m': 'facebook/esm2_t12_35M_UR50D',
            'esm2_150m': 'facebook/esm2_t30_150M_UR50D',
            'esm2_650m': 'facebook/esm2_t33_650M_UR50D',
            'saprot_35m': 'westlake-repl/SaProt_35M_AF2',
            'saprot_650m': 'westlake-repl/SaProt_650M_AF2'
        }
    }

    available_models = model_paths['local'].keys()
    if model_name not in available_models:
        raise NotImplementedError(f"Model {model_name} is not available.\nAvailable models: {list(available_models)}")

    # get model path based on flag
    model_path = model_paths['local' if local_files_only else 'online'][model_name]

    # map model prefixes to classes
    model_class_map = {
        'esmc': ESMC_Embed,
        'esm3': ESM3_Embed,
        'esm2': ESM2_Embed,
        'saprot': ESM2_Embed
    }

    # determine appropriate class
    for prefix, cls in model_class_map.items():
        if model_name.startswith(prefix):
            return cls(model_path, local_files_only=local_files_only)

    raise NotImplementedError(f"Model {model_name} is not supported.")


# Example usage
if __name__ == '__main__':
    # command
    # cd LLPS_Predict
    # python -m src.models.embed_seq

    # SQN_EMBED_MODEL options:
    # 'esmc_300m': ESMC 300M, embedding shape (sqn_length+2, 960)
    # 'esmc_600m': ESMC 600M, embedding shape (sqn_length+2, 1152)
    # 'esm3_sm': ESM3 SM, embedding shape (sqn_length+2, 1536)
    # 'esm2_8m': ESM2 6 Layer 8M, embedding shape (sqn_length+2, 320)
    # 'esm2_35m': ESM2 12 Layer 35M, embedding shape (sqn_length+2, 480)
    # 'esm2_150m': ESM2 30 Layer 150M, embedding shape (sqn_length+2, 640)
    # 'esm2_650m': ESM2 33 Layer 650M, embedding shape (sqn_length+2, 1280)
    # 'saprot_35m': SaProt 35M, embedding shape (sqn_length+2, 480)
    # 'saprot_650m': SaProt 650M, embedding shape (sqn_length+2, 1280)
    SQN_EMBED_MODEL = 'esm2_650m'
    # SQN_EMBED_LOAD_LOCAL: whether to load model weights from local
    # True if model weights have been downloaded
    # False if load online model weights, need Internet connection
    SQN_EMBED_LOAD_LOCAL = True
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # sample input
    if SQN_EMBED_MODEL.startswith('esm'):
        # ESMC, ESM3, ESM2 sample input
        sequence_list = (
            'AGFQPQSQGMSLNDFQKQQKQAAPKPKKTLKLVSSSGIKLANATKKVGTKPAESDKKEEEKSAETKEPTKEPTKVEEPVKKEEKPVQTEEKTEEKSELPKVEDLKISESTHNTNNANVTSADALIKEQEEEVDDEVVNDMFGGKDHVSLIFMGHVDAGKSTMGGNLLYLTGSVDKRTIEKYEREAKDAGRQGWYLSWVMDTNKEERNDGKTIEVGKAYFETEKRRYTILDAPGHKMYVSEMIGGASQADVGVLVISARKGEYETGFERGGQTREHALLAKTQGVNKMVVVVNKMDDPTVNWSKERYDQCVSNVSNFLRAIGYNIKTDVVFMPVSGYSGANLKDHVDPKECPWYTGPTLLEYLDTMNHVDRHINAPFMLPIAAKMKDLGTIVEGKIESGHIKKGQSTLLMPNKTAVEIQNIYNETENEVDMAMCGEQVKLRIKGVEEEDISPGFVLTSPKNPIKSVTKFVAQIAIVELKSIIAAGFSCVMHVHTAIEEVHIVKLLHKLEKGTNRKSKKPPAFAKKGMKVIAVLETEAPVCVETYQDYPQLGRFTLRDQGTTIAIGKIVKIAE',
            'ASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNSSSGGGGGGGGGGNYGQDQSSMSSGGGSGGGYGNQDQSGGGGSGGYGQQDRG',
            'CLRKKRKPQAEKVDVIAGSSKMKGFSSSESESSSESSSSDSEDSETEMAPKSKKKGHPGREQKKHHHHHHQQMQQAPAPVPQQPPPPPQQPPPPPPPQQQQQPPPPPPPPSMPQQAAPAMKSSPPPFIATQVPVLEPQLPGSVFDPIGHFTQPILHLPQPELPPHLPQPPEHSTPPHLNQHAVVSPPALHNALPQQPSRPSNRAAALPPKPARPPAVSPALTQTPLLPQPPMAQPPQVLLEDEEPPAPPLTSMQMQLYLQQLQKVQPPTPLLPSVKVQSQPPPPLPPPPHPSVQQQLQQQPPPPPPPQPQPPPQQQHQPPPRPVHLQPMQFSTHIQQPPPPQGQQPPHPPPGQQPPPPQPAKPQQVIQHHHSPRHHKSDPYSTGHLREAPSPLMIHSPQMSQFQSLTHQSPPQQNVQPKKQELRAASVVQPQPLVVVKEEKIHSPIIRSEPFSPSLRPEPPKHPESIKAPVHLPQRPEMKPVDVGRPVIRPPEQNAPPPGAPDKDKQKQEPKTPVAPKKDLKIKNMGSWASLVQKHPTTPSSTAKSSSDSFEQFRRAAREKEEREKALKAQAEHAEKEKERLRQERMRSREDEDALEQARRAHEEARRRQEQQQQQRQEQQQQQQQQAAAVAAAATPQAQSSQPQSMLDQQRELARKREQERRRREAMAATIDMNFQS'
        )
    else:
        # SaProt sample input
        # sequence tokens intertwined with structure tokens
        # use placeholder '#' for low-confidence structure
        sequence_list = ('M#EvVpQpL#VyQdYaKv',)

    # initialize model
    model = get_sqn_embed_model(
        SQN_EMBED_MODEL, local_files_only=SQN_EMBED_LOAD_LOCAL
    )
    model = model.to(DEVICE)

    # embed each protein sequence
    model.eval()
    with torch.no_grad():
        outputs, attention_masks = model(sequence_list)

    print('Model:', SQN_EMBED_MODEL)
    print('Batch embedding shape:', outputs.shape)
    print('Batch attention mask shape:', attention_masks.shape)
    print('Number of proteins:', len(outputs))
    for i, mask in enumerate(attention_masks):
        print(f"Sequence [{i}] length: {int(sum(mask==1))}")

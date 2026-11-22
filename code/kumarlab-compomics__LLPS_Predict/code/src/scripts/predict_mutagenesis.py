##########
# Script for performing mutagenesis analysis for a given protein sequence.
# For each position, substitute the amino acid with every other amino acid, 
# and compute the delta LLPS mutation score.
#
# Usage: 
# cd LLPS_Predict
# python -m src.scripts.predict_mutagenesis [--config c] [--data_path d] [--score_path s]
#   [--dropout p]
#
# Override config parameters with command line arguments if provided
##########

# torch
import torch
from src.models.model_ensemble import EnsembleClassifier
# moving window helper function
from src.utils.window import VOCABULARY, compute_mutagenesis
from tqdm import tqdm
# config parameters
import yaml
import argparse
# results
import pandas as pd


# command line parameters
parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, required=True)
parser.add_argument('--uniprotid', type=str)
parser.add_argument('--genesymbol', type=str)
parser.add_argument('--data_path', type=str)
parser.add_argument('--score_path', type=str)
parser.add_argument('--dropout', type=float)
args = parser.parse_args()
# parse yaml config file
with open(args.config, 'r') as f:
    config = yaml.safe_load(f)
# override with command line arguments if provided
args_dict = {
    k: v for k, v in vars(args).items()
    if v is not None and k != "config"
}
config.update(args_dict)

# datasets
TEST_DATA_PATH = config['data_path']
SCORE_PATH = config['score_path']
SQN_COLUMN1 = config['sqn_column1']
SQN_COLUMN2 = config['sqn_column2']

# model
SQN_EMBED_MODEL1 = config['sqn_embed_model1']
SQN_EMBED_MODEL2 = config['sqn_embed_model2']
SQN_EMBED_LOAD_LOCAL = config['sqn_embed_load_local']
SQN_EMBED_DIM = config['sqn_embed_dim']
N_LAYER = config['n_layer']
N_HEAD = config['n_head']
FEEDFORWARD_DIM = config['feedforward_dim']
CHECKPOINT_PATH = config['checkpoint_path']

# configuration
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
DROPOUT = config['dropout']

# protein info
UNIPROTID = config['uniprotid']
UNIPROTID_COLUMN = config['uniprotid_column']
GENESYMBOL = config['genesymbol']
GENESYMBOL_COLUMN = config['genesymbol_column']

# input
df = pd.read_csv(TEST_DATA_PATH)
mask = df[UNIPROTID_COLUMN] == UNIPROTID
sequence = df[SQN_COLUMN1][mask].iloc[0]
saprot_sequence = df[SQN_COLUMN2][mask].iloc[0]

# model
classifier = EnsembleClassifier(
    embed_model1=SQN_EMBED_MODEL1,
    embed_model2=SQN_EMBED_MODEL2,
    embedding_dim=SQN_EMBED_DIM,
    num_layers=N_LAYER,
    num_heads=N_HEAD,
    dim_feedforward=FEEDFORWARD_DIM,
    dropout=DROPOUT,
    local_files_only=SQN_EMBED_LOAD_LOCAL
)
classifier = classifier.to(DEVICE)

# checkpoint path
checkpoint = torch.load(CHECKPOINT_PATH, weights_only=False)
classifier.load_state_dict(checkpoint['model'], strict=False)

# compute delta LLPS scores by changing the residue at each position
print(f"Running mutagesis analysis for gene: {GENESYMBOL}")
all_logits = compute_mutagenesis(classifier, sequence, saprot_sequence)
df_pred = pd.DataFrame(all_logits, columns=VOCABULARY, index=range(len(sequence)))
df_pred.to_csv(SCORE_PATH, index=False)
print(f"Finished mutagenesis analysis for gene: {GENESYMBOL}")

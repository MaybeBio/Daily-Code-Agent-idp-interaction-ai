##########
# Script for calculating delta LLPS mutation score
#
# Usage: 
# cd LLPS_Predict
# python -m src.scripts.predict_mutation [--config c] [--data_path d] [--score_path s]
#   [--dropout p]
#
# Override config parameters with command line arguments if provided
##########

# torch
import torch
from src.models.model_ensemble import EnsembleClassifier
# moving window helper function
from src.utils.window import compute_llps, compute_llps_window
from tqdm import tqdm
# config parameters
import yaml
import argparse
# results
import pandas as pd


# command line parameters
parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, required=True)
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
REF_COLUMN1 = config['ref_column1']
REF_COLUMN2 = config['ref_column2']
ALT_COLUMN1 = config['alt_column1']
ALT_COLUMN2 = config['alt_column2']
POSITION_COLUMN = config['position_column']

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
MODE = config['mode']


# input
df = pd.read_csv(TEST_DATA_PATH)

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

# compute LLPS scores for ClinVar wild types vs. disease-related mutations
for i in tqdm(range(len(df))):
    ref_seq = df.loc[i, REF_COLUMN1]
    ref_saprot_seq = df.loc[i, REF_COLUMN2]
    alt_seq = df.loc[i, ALT_COLUMN1]
    alt_saprot_seq = df.loc[i, ALT_COLUMN2]
    pos_aa = int(df.loc[i, POSITION_COLUMN])
    if MODE == "protein_level":
        ref_prob = compute_llps(classifier, ref_seq, ref_saprot_seq)
        alt_prob = compute_llps(classifier, alt_seq, alt_saprot_seq)
    else:
        ref_prob = compute_llps_window(classifier, ref_seq, ref_saprot_seq, pos_aa-1)
        alt_prob = compute_llps_window(classifier, alt_seq, alt_saprot_seq, pos_aa-1)
    # mutation score = LLPS(mutation) - LLPS(wild type)
    score = alt_prob - ref_prob
    df.loc[i, 'MutationScore'] = score
    df.loc[i, 'RefScore'] = ref_prob
    df.loc[i, 'AltScore'] = alt_prob

df2 = df[['MutationScore', 'RefScore', 'AltScore']]
df2.to_csv(SCORE_PATH, index=False)

# torch
import torch
from src.models.model_ensemble import EnsembleClassifier
# moving window helper function
from src.utils.window import compute_critical_region
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

# compute LLPS score at residue level
for i in tqdm(range(len(df))):
    sequence = df.loc[i, SQN_COLUMN1]
    saprot_sequence = df.loc[i, SQN_COLUMN2]
    prob_list = compute_critical_region(classifier, sequence, saprot_sequence)
    llps = ';'.join([str(round(p, 2)) for p in prob_list])
    df.loc[i, 'Prediction'] = llps

df2 = df[['Prediction']]
df2.to_csv(SCORE_PATH, index=False)

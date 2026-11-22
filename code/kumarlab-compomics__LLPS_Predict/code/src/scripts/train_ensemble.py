##########
# Script for training sequence-structure model (ESM2 & SaProt)
#
# Usage: 
# cd LLPS_Predict
# python -m src.scripts.train_ensemble [--config c] [--data_path d] [--batch_size n] 
#   [--dropout p] [--epochs e] [--lr l] [--weight_decay wd]
#
# Override config parameters with command line arguments if provided
##########

# torch
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from src.loaders.dataset_pair import LLPSDatasetPair
from src.models.model_ensemble import EnsembleClassifier
# config parameters
import argparse
import yaml
# training loop
from src.utils.trainer import train


# command line parameters
parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, required=True)
parser.add_argument('--data_path', type=str)
parser.add_argument('--batch_size', type=int)
parser.add_argument('--dropout', type=float)
parser.add_argument('--epochs', type=int)
parser.add_argument('--lr', type=float)
parser.add_argument('--weight_decay', type=float)
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
TRAIN_DATA_PATH = config['data_path']
SQN_COLUMN1 = config['sqn_column1']
SQN_COLUMN2 = config['sqn_column2']
LABEL_COLUMN = config['label_column']

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
RESUME_TRAINING = config['resume_training']
SEED = config['seed']
BATCH_SIZE = config['batch_size']
DROPOUT = config['dropout']
EPOCHS = config['epochs']
LR = config['lr']
WEIGHT_DECAY = config['weight_decay']

print('Data path:', TRAIN_DATA_PATH)
print('Sequence embedding model1:', SQN_EMBED_MODEL1)
print('Sequence embedding model2:', SQN_EMBED_MODEL2)
print('Seed:', SEED)
print('Batch size:', BATCH_SIZE)
print('Dropout:', DROPOUT)
print('Epochs:', EPOCHS)
print('Learning rate:', LR)
print('Weight decay:', WEIGHT_DECAY)


# Dataloader
dataset_train = LLPSDatasetPair(
    TRAIN_DATA_PATH, SQN_COLUMN1, SQN_COLUMN2, LABEL_COLUMN
)
dataloader_train = DataLoader(
    dataset_train, batch_size=BATCH_SIZE, shuffle=True, drop_last=True
)

# Training
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

# model checkpoint
if RESUME_TRAINING:
    checkpoint = torch.load(CHECKPOINT_PATH, weights_only=False)
    classifier.load_state_dict(checkpoint['model'], strict=False)

# optimizer
optimizer = torch.optim.Adam(classifier.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
# if RESUME_TRAINING:
#     optimizer.load_state_dict(checkpoint['optimizer'])

# loss function
criterion = nn.BCEWithLogitsLoss()

# train model
torch.manual_seed(SEED)
losses_train = train(
    model=classifier,
    dataloader=dataloader_train,
    optimizer=optimizer,
    criterion=criterion,
    epochs=EPOCHS,
    checkpoint_path=CHECKPOINT_PATH,
    device=DEVICE,
    verbose=True
)

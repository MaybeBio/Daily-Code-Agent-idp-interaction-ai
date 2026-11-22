##########
# Script for training intermap model, 
# CNN model that predicts LLPS propensity from force field intermaps (FINCHES)
#
# Usage:
# cd LLPS_Predict
# python -m src.scripts.train_intermap [--config c] [--data_path d] [--batch_size n] 
#   [--dropout p] [--epochs e] [--lr l] [--weight_decay wd]
#
# Override config parameters with command line arguments if provided
##########

# torch
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from src.loaders.dataset import LLPSDataset
from src.models.model_intermap import IntermapClassifier
# FINCHES force field intermap
from src.utils.force_field import force_field_transform_factory
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
SQN_COLUMN = config['sqn_column']
LABEL_COLUMN = config['label_column']

# model
IMG_SIZE = config['img_size']
N_LAYER = config['n_layer']
N_HEAD = config['n_head']
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
print('Model: Force Field Intermap CNN')
print('Seed:', SEED)
print('Batch size:', BATCH_SIZE)
print('Dropout:', DROPOUT)
print('Epochs:', EPOCHS)
print('Learning rate:', LR)
print('Weight decay:', WEIGHT_DECAY)


# Dataloader
seq_transform = force_field_transform_factory(img_size=IMG_SIZE)
dataset_train = LLPSDataset(
    TRAIN_DATA_PATH, SQN_COLUMN, LABEL_COLUMN,
    seq_transform=seq_transform
)
dataloader_train = DataLoader(
    dataset_train, batch_size=BATCH_SIZE, shuffle=True, drop_last=True
)

# Training
# model
classifier = IntermapClassifier(in_channels=1, dropout=DROPOUT)
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

##########
# Script for testing intermap model, 
# CNN model that predicts LLPS propensity from force field intermaps (FINCHES)
#
# Usage: 
# cd LLPS_Predict
# python -m src.scripts.test_intermap [--config c] [--data_path d] [--score_path s]
#   [--batch_size n] [--dropout p] [threshold t]
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
# evaluation metrics
from src.utils.metrics import classification_metrics
# config parameters
import yaml
import argparse
# predict loop
from src.utils.predictor import predict
# results
import pandas as pd


# command line parameters
parser = argparse.ArgumentParser()
parser.add_argument('--config', type=str, required=True)
parser.add_argument('--data_path', type=str)
parser.add_argument('--score_path', type=str)
parser.add_argument('--batch_size', type=int)
parser.add_argument('--dropout', type=float)
parser.add_argument('--threshold', type=float)
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
SQN_COLUMN = config['sqn_column']
LABEL_COLUMN = config['label_column']
SCORE_PATH = config['score_path']

# model
IMG_SIZE = config['img_size']
N_LAYER = config['n_layer']
N_HEAD = config['n_head']
CHECKPOINT_PATH = config['checkpoint_path']
THRESHOLD = config['threshold']

# configuration
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
BATCH_SIZE = config['batch_size']
DROPOUT = config['dropout']

print('Data path:', TEST_DATA_PATH)

# Dataloader
seq_transform = force_field_transform_factory(img_size=IMG_SIZE)
dataset_test = LLPSDataset(
    TEST_DATA_PATH, SQN_COLUMN, LABEL_COLUMN,
    seq_transform=seq_transform
)
dataloader_test = DataLoader(
    dataset_test, batch_size=BATCH_SIZE, shuffle=False
)

# Testing
# model
classifier = IntermapClassifier(in_channels=1, dropout=DROPOUT)
classifier = classifier.to(DEVICE)

# checkpoint path
checkpoint = torch.load(CHECKPOINT_PATH, weights_only=False)
classifier.load_state_dict(checkpoint['model'], strict=False)

# loss function
criterion = nn.BCEWithLogitsLoss()

# compute test loss
# get predictions
pred = predict(
    model=classifier, 
    dataloader=dataloader_test,
    criterion=criterion,
    device=DEVICE,
    verbose=True
)
sequence = dataset_test.sequences
target = dataset_test.labels
df_pred = pd.DataFrame({
    'Sequence': sequence, 'Label': target, 'Prediction': pred
})
df_pred.to_csv(SCORE_PATH, index=False)


# Compute Metrics
classification_metrics(target, pred, threshold=THRESHOLD)

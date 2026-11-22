# PyTorch style dataset class for LLPS prediction
import torch
import pandas as pd


class LLPSDataset(torch.utils.data.Dataset):
    """
    Dataset Class for sequential data (input) and labels (ground truth).
    e.g., (sequence, label)
    """
    def __init__(self, df_path: str, sqn_column: str, label_column: str, seq_transform=None):
        """
        :param df_path: path to dataframe
        :type df_path: str
        :param sqn_column: column name for input sequences, e.g., 'Sequence', 'SaProtSeq'
        :type sqn_column: str
        :param label_column: column name for ground truth labels, e.g., 'Label'
        :type label_column: str
        :param seq_transform: transform function for input sequence
        :type seq_transform: Callable[str, Any]
        """
        df = pd.read_csv(df_path)
        self.sequences = df[sqn_column]
        self.labels = df[label_column]
        self.seq_transform = seq_transform

    def __len__(self):
        return len(self.sequences)
    
    def __getitem__(self, index):
        seq = self.sequences[index]
        label = self.labels[index]
        if self.seq_transform is not None:
            seq = self.seq_transform(seq)
        return seq, label


if __name__ == '__main__':
    # command:
    # cd LLPS_Predict
    # python -m src.loaders.dataset

    # input config
    LLPS_DATA_PATH = 'data/training/llps_vs_non_llps_train.csv'
    SQN_COLUMN = 'Sequence'
    LABEL_COLUMN = 'Label'

    dataset = LLPSDataset(LLPS_DATA_PATH, sqn_column=SQN_COLUMN, label_column=LABEL_COLUMN)
    
    print('Number of entries:', len(dataset))
    print('Sequence [0]:', dataset[0][0])
    print('Label [0]:', dataset[0][1])

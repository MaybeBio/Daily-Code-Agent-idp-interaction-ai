# PyTorch style dataset classes for LLPS prediction
import torch
import pandas as pd


class LLPSDatasetPair(torch.utils.data.Dataset):
    """
    Dataset class for pairwise sequential data (input) and labels (ground truth).
    e.g., (sequence, SaProt sequence, label)
    """
    def __init__(self, df_path: str, sqn_column1: str, sqn_column2: str, label_column: str, 
                 seq1_transform=None, seq2_transform=None):
        """
        :param df_path: path to dataframe
        :type df_path: str
        :param sqn_column1: column name for first type of input sequences, e.g., 'Sequence'
        :type sqn_column1: str
        :param sqn_column2: column name for second type of input sequences, e.g., 'SaProtSeq'
        :type sqn_column2: str
        :param label_column: column name for ground truth labels, e.g., 'Label'
        :type label_column: str
        :param seq1_transform: transform function for first type of input sequence
        :type seq1_transform: Callable[str, Any]
        :param seq2_transform: transform function for second type of input sequence
        :type seq2_transform: Callable[str, Any]
        """
        df = pd.read_csv(df_path)
        self.sequences1 = df[sqn_column1]
        self.sequences2 = df[sqn_column2]
        self.labels = df[label_column]
        self.seq1_transform = seq1_transform
        self.seq2_transform = seq2_transform

    def __len__(self):
        return len(self.sequences1)
    
    def __getitem__(self, index):
        seq1 = self.sequences1[index]
        seq2 = self.sequences2[index]
        label = self.labels[index]
        if self.seq1_transform is not None:
            seq1 = self.seq1_transform(seq1)
        if self.seq2_transform is not None:
            seq2 = self.seq2_transform(seq2)
        return seq1, seq2, label


if __name__ == '__main__':
    # command:
    # cd LLPS_Predict
    # python -m src.loaders.dataset_pair

    # input config
    LLPS_DATA_PATH = 'data/training/llps_vs_non_llps_train.csv'
    SQN_COLUMN1 = 'Sequence'
    SQN_COLUMN2 = 'SaProtSeq'
    LABEL_COLUMN = 'Label'

    dataset = LLPSDatasetPair(
        LLPS_DATA_PATH, 
        sqn_column1=SQN_COLUMN1, 
        sqn_column2=SQN_COLUMN2, 
        label_column=LABEL_COLUMN
    )
    
    print('Number of entries:', len(dataset))
    print('Sequence [0]:', dataset[0][0])
    print('SaProt sequence [0]:', dataset[0][1])
    print('Label [0]:', dataset[0][2])

# PyTorch style dataset classes for LLPS prediction
import torch
import pandas as pd


class LLPSDatasetTriplet(torch.utils.data.Dataset):
    """
    Dataset class for triplet sequential data (input) and labels (ground truth).
    e.g., (sequence, SaProt sequence, intermapCNN, label)
    """
    def __init__(self, df_path: str, sqn_column1: str, sqn_column2: str, sqn_column3: str,
                 label_column: str, seq1_transform=None, seq2_transform=None, seq3_transform=None):
        """
        :param df_path: path to dataframe
        :type df_path: str
        :param sqn_column1: column name for first type of input sequences, e.g., 'Sequence'
        :type sqn_column1: str
        :param sqn_column2: column name for second type of input sequences, e.g., 'SaProtSeq'
        :type sqn_column2: str
        :param sqn_column2
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
        self.sequences3 = df[sqn_column3]
        self.labels = df[label_column]
        self.seq1_transform = seq1_transform
        self.seq2_transform = seq2_transform
        self.seq3_transform = seq3_transform

    def __len__(self):
        return len(self.sequences1)
    
    def __getitem__(self, index):
        seq1 = self.sequences1[index]
        seq2 = self.sequences2[index]
        seq3 = self.sequences3[index]
        label = self.labels[index]
        if self.seq1_transform is not None:
            seq1 = self.seq1_transform(seq1)
        if self.seq2_transform is not None:
            seq2 = self.seq2_transform(seq2)
        if self.seq3_transform is not None:
            seq3 = self.seq3_transform(seq3)
        return seq1, seq2, seq3, label


if __name__ == '__main__':
    # command:
    # cd LLPS_Predict
    # python -m src.loaders.dataset_triplet

    # transform image
    from src.utils.force_field import force_field_transform_factory

    # input config
    LLPS_DATA_PATH = 'data/training/llps_vs_non_llps_train.csv'
    INTERMAP_SIZE = 320
    SQN_COLUMN1 = 'Sequence'
    SQN_COLUMN2 = 'Sequence'
    SQN_COLUMN3 = 'SaProtSeq'
    LABEL_COLUMN = 'Label'

    # transform from sequence to force field intermap
    force_field_transform = force_field_transform_factory(img_size=INTERMAP_SIZE)

    dataset = LLPSDatasetTriplet(
        LLPS_DATA_PATH, 
        sqn_column1=SQN_COLUMN1,
        sqn_column2=SQN_COLUMN2,
        sqn_column3=SQN_COLUMN3,
        label_column=LABEL_COLUMN,
        seq1_transform=force_field_transform
    )
    
    print('Number of entries:', len(dataset))
    print('Force field intermap [0] shape:', dataset[0][0].shape)
    print('Sequence [0]:', dataset[0][1])
    print('SaProt sequence [0]:', dataset[0][2])
    print('Label [0]:', dataset[0][3])

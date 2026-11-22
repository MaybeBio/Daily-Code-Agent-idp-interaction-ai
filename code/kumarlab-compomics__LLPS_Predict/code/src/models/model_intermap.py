# PyTorch style intermap CNN model class
import torch
import torch.nn as nn
from src.models.embed_intermap import IntermapEmbed


class IntermapClassifier(IntermapEmbed):
    """
    CNN model class that uses FINCHES force field intermap for prediction.
    """
    def __init__(self, in_channels: int = 1, dropout: float = 0.2):
        """
        :param in_channels: number of input channels, 1 for force field intermaps
        :type in_channels: int
        :param dropout: dropout rate in final MLP
        :type dropout: float
        """
        super(IntermapClassifier, self).__init__(
            in_channels, output_dim=1, dropout=dropout
        )

    def forward(self, x):
        """
        :param x: input intermaps; 
            shape (B, L, L), B: batch size, L: sequence length
        :type x: torch.Tensor
        :return: embedding vector; 
            shape (B,), B: batch size
        :rtype: torch.Tensor
        """
        x = self.conv_block(x)
        x = self.global_pool(x)
        x = self.fc(x)
        return x.squeeze(-1)    # return shape (B,)


# Example usage
if __name__ == "__main__":
    # command
    # cd LLPS_Predict
    # python -m src.models.embed_intermap
    
    import torch

    # input config
    IMG_SIZE = 320
    DROPOUT = 0.2
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # random input as force field intermap
    x = torch.randn(3, 1, IMG_SIZE, IMG_SIZE).to(DEVICE)

    # model
    classifier = IntermapClassifier(in_channels=1, dropout=DROPOUT)
    classifier = classifier.to(DEVICE)

    # forward pass
    x = classifier(x)

    print('Model: IntermapCNN Clasiifier')
    print('Logits shape:', x.shape)    # Expected: (3,)
    print('Logits:', x)

# PyTorch style model class for embedding force field intermap
import torch.nn as nn


class IntermapEmbed(nn.Module):
    """
    Intermap CNN model class that embeds intermaps into vectors.
    """
    def __init__(self, in_channels: int = 1, output_dim: int = 1280, dropout: float = 0.2):
        """
        :param in_channels: number of input channels, default 1
        :type in_channels: int
        :param output_dim: dimension of each output vector
        :type output_dim: int
        :param dropout: dropout rate for final MLP
        :type dropout: float
        """
        super(IntermapEmbed, self).__init__()

        self.conv_block = nn.Sequential(
            nn.Conv2d(in_channels, 16, kernel_size=3, padding=1),  # (B, 16, L, L)
            nn.ReLU(),
            nn.BatchNorm2d(16),

            nn.Conv2d(16, 32, kernel_size=3, padding=1),           # (B, 64, L, L)
            nn.ReLU(),
            nn.BatchNorm2d(32),

            nn.Conv2d(32, 64, kernel_size=3, padding=1),           # (B, 256, L, L)
            nn.ReLU(),
            nn.BatchNorm2d(64)
        )

        self.global_pool = nn.AdaptiveAvgPool2d((1, 1))  # (B, 256, 1, 1)
        
        self.fc = nn.Sequential(
            nn.Flatten(),                 # (B, 256)
            nn.Linear(64, 256),           # intermediate
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(256, output_dim) # final embedding
        )

    def forward(self, x):
        """
        :param x: input intermaps; 
            shape (B, L, L), B: batch size, L: sequence length
        :type x: torch.Tensor
        :return: embedding vector; 
            shape (B, D), B: batch size, D: embedding dimension
        :rtype: torch.Tensor
        """
        x = self.conv_block(x)
        x = self.global_pool(x)
        x = self.fc(x)
        return x # return shape (B, D)


if __name__ == '__main__':
    # command
    # cd LLPS_Predict
    # python -m src.models.embed_intermap
    
    import torch

    # input config
    IMG_SIZE = 320
    EMBED_DIM = 1280
    DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # random input as force field intermap
    x = torch.randn(3, 1, IMG_SIZE, IMG_SIZE).to(DEVICE)

    # model
    cnn = IntermapEmbed(in_channels=1, output_dim=EMBED_DIM)
    cnn = cnn.to(DEVICE)

    # forward pass
    x = cnn(x)

    print('Model: IntermapCNN')
    print('Batch embedding shape:', x.shape)    # Expected: (3, 1280)
    print('Batch size:', len(x))        # Expected: 3

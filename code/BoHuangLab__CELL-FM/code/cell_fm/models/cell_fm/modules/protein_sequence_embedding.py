import torch
import torch.nn as nn
from esm.models.esmc import ESMC

class ESMEmbed(nn.Module):
    def __init__(self, model_name, out_dims, fixed_embedding=True):
        super(ESMEmbed, self).__init__()

        if model_name == 'esmc_300m':
            self.model = ESMC.from_pretrained("esmc_300m")
            model_dim = 960
        elif model_name == 'esmc_600m':
            self.model = ESMC.from_pretrained("esmc_600m")
            model_dim = 1152

        if out_dims != model_dim:
            self.scale_layer = nn.Linear(model_dim, out_dims)
        else:
            self.scale_layer = nn.Identity()

        # Determine whether to freeze the model's parameters
        self.fixed_embedding = fixed_embedding
        if self.fixed_embedding:
            for param in self.model.parameters():
                param.requires_grad = False
            self.model = self.model.eval()

    def forward(self, x, **kwargs):
        # If the model's parameters are fixed, use torch.no_grad()
        if self.fixed_embedding:
            with torch.no_grad():
                x = self.model(x).embeddings
        else:
            x = self.model(x).embeddings
        
        x = x.to(torch.float32)
        x = self.scale_layer(x)

        return x
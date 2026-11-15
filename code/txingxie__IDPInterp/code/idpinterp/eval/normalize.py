import logging

import torch
from torch.utils.data import DataLoader
from tqdm import tqdm

logger = logging.getLogger(__name__)


class FeatureNormalizer:
    """Normalizes SAE feature activations to the [0, 1] range based on empirical max."""

    @staticmethod
    def compute_max_activations(
        sae: torch.nn.Module,
        dataset: torch.utils.data.Dataset,
        batch_size: int = 4096,
        device: str = "cuda",
    ) -> torch.Tensor:
        """
        Sweep over the dataset to find the maximum activation for each feature.
        """
        loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=2)
        sae = sae.eval().to(device)

        max_acts = torch.zeros(sae.d_hidden, device=device)

        with torch.no_grad():
            for batch in tqdm(loader, desc="Finding max activations"):
                batch = batch.to(device)
                acts = sae.encode(batch)
                batch_max = acts.max(dim=0).values
                max_acts = torch.maximum(max_acts, batch_max)

        return max_acts

    @staticmethod
    def normalize_sae(sae: torch.nn.Module, max_acts: torch.Tensor) -> None:
        """
        Adjust the SAE decoder weights to output normalized activations.
        
        We divide the decoder weights by the max activation, and scale the encoder 
        to match, so the relative scale of features is preserved but the max is 1.0.
        """
        # Avoid division by zero for dead features
        scale = torch.where(max_acts > 1e-6, max_acts, torch.ones_like(max_acts))
        
        with torch.no_grad():
            # Adjust encoder weights and bias
            sae.W_enc.data /= scale.unsqueeze(0)
            sae.b_enc.data /= scale
            
            # Adjust decoder weights
            sae.W_dec.data *= scale.unsqueeze(1)
            
        logger.info("SAE features normalized to [0, 1] range.")

from src.models.model_ensemble import EnsembleClassifier
from src.utils.window import compute_mutagenesis


def get_mutagenesis(model, sequence: str, saprot_sequence: str):
    if not isinstance(model, EnsembleClassifier):
        raise NotImplementedError(f"Only support ensemble model for mutagenesis analysis.")
    return compute_mutagenesis(model, sequence, saprot_sequence)


if __name__ == '__main__':
    # LLPS models
    import torch
    from src.lib.llps_model import get_llps_model
    # input
    from src.utils.saprotseq import get_saprotseq
    # output
    import numpy as np

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # model
    model = get_llps_model('ensemble', local_files_only=True)
    model = model.to(device)
    # input
    # sequence: str
    # saprotseq: str
    sequence = 'MEELSADEIRRRRLARLAGGQTS'
    saprotseq = get_saprotseq(sequence)
    # inference
    scores = get_mutagenesis(model, sequence, saprotseq)
    scores = np.array(scores)
    print("Sequence length:", len(sequence))
    print("Mutagenesis heatmap shape:", scores.shape)   # Expected: (L, 20)
    print("Mutagenesis heatmap:", scores)

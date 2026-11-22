# FINCHES force field intermap helper function

# FINCHES force field intermap
from finches import Mpipi_frontend
# transform image
from torchvision.transforms import Resize, ToTensor, Compose


def force_field_transform_factory(img_size: int):
    # FINCHES force field intermap
    mf = Mpipi_frontend()
    # reduce intermap dimension to given size for ease of computation
    resize = Compose([ToTensor(), Resize(img_size)])
    # transform function to be used in the Dataset class
    def transform(seq):
        seq = seq.replace('X', 'G')   # replace with dummy amino acid G
        idr, _, _ = mf.intermolecular_idr_matrix(seq, seq, window_size=1)[0]
        idr = resize(idr)
        return idr.float()
    return transform


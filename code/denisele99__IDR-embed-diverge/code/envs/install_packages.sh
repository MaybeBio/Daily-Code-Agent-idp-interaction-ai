
mamba install numpy pandas scipy scikit-learn 

-itertools
pygosemsim
dataclasses
pathlib
typing
yaml
#distances
scipy
math
collections
json
sklearn
h5py
glob
multiprocessing
torch
Biopython
sys
re

import obonet
from goatools

#Go enrichments
mamba install numpy pandas scikit-learn=1.4.0 h5py pyyaml=6.0.1
pip install goatools==1.4.12 biopython==1.81
pip install git+https://github.com/mojaie/pygosemsim.git



y
mamba install scipy=1.12.0 scikit-learn=1.4.0 h5py pyyaml=6.0.1
pip install pygosemsim==0.1.0 goatools==1.4.12 biopython==1.81 torch==2.1.0
pip install git+https://github.com/mojaie/pygosemsim.git
dependenceise

    - biopython=1.81=py312h98912ed_1
  - h5py=3.14.0=nompi_py312h3faca00_100
  - hdf5=1.14.6=nompi_h2d575fe_101
  - pytorch=2.1.0=cuda120py312hfe5e8c6_301

pip
- goatools==1.4.12
  - multiprocess==0.70.16
- pygosemsim==0.1.0
 - pyyaml==6.0.1
- scikit-learn==1.4.0
      - scipy==1.12.0


#go enrichments



Python 3.12.7

#visualization packages
scipy
plotly
umap

mamba install matplotlib seaborn



- matplotlib==3.8.2

#generating embeddings

biopython
python Python 3.12.7
pickle
transformers
yaml
h5py

mamba install h5py pyyaml=6.0.1
pip install biopython==1.81 torch==2.1.0 transformers==4.37.2

#training and/or generating embeddings
wandb
os 
pathlib
numpy
pandas
time
datasets
typing
argparse
transformers
yaml
h5py
pytorch (which one?)
collections
pickle
re


mamba install numpy pandas h5py pyyaml 

pip install torch==2.3.0 transformers datasets wandb
conda env export --from-history > environment.yml


rsync -avz --exclude='/home/moseslab/denise/IDR_LM/data' \
  moseslab@soma.csb.utoronto.ca:/home/moseslab/denise/IDR_LM \
  moseslab@lindt.csb.utoronto.ca:/lindt/denise/soma_data
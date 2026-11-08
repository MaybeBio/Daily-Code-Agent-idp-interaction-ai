#!/bin/bash

# Unzip ProtTrans dir.
tar -xzvf ./ProtTrans.tar.gz

# Export disobind path.
if [ -z "$PYTHONPATH" ]; then
    echo "export PYTHONPATH=`pwd`:$PYTHONPATH" >> ~/.bash_profile
else
    echo "export PYTHONPATH=$PYTHONPATH:`pwd`" >> ~/.bash_profile
fi
source ~/.bash_profile

# Create conda environment for Disobind.
conda create --name diso python=3.9

# Activate the conda environment.
conda activate diso

# Install dependencies.
pip install -r ./requirements.txt

# Deactivate conda environment.
conda deactivate

echo "You are now ready to roll... Don't forget to activate the conda anvironment using 'conda activate diso' before using Disobind."
echo "May the Force serve you well..."
#!/bin/sh
export PATH=/home/xie/anaconda3/bin:${PATH}
export KKTOOLS=/home/xie/local/kktools
export OMEGATK=/home/xie/local/omegagene/toolkit
export PYTHONPATH=$KKTOOLS:$OMEGATK:$PYTHONPATH

python3 ${KKTOOLS}/mda_pca_traj_selectatoms.py \
  --trr ../P_nonP_ca.trr \
  --ref ../run_ca.pdb  \
  -o pca_traj.txt  --dim 5 \
  --o-eig pca_eig.dat \
  --coordinate-type distance \
  --atom-pairs 0-70:71-88 \
  --min-prob 0.00 \
  --i-prob ../prob_mix.txt --non-log-prob  --col-prob 0 \
    > log.txt

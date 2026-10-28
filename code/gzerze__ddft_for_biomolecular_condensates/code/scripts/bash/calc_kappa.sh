#!/bin/bash
#SBATCH -J kappa_calc
#SBATCH -o kappa.o%j
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks-per-node=1 -N 1
#SBATCH -t 96:00:00
#SBATCH --mem-per-cpu=32GB

#module --ignore_cache load "OpenMPI/.4.1.4-GCC-11.3.0"
export OMP_NUM_THREADS=1

export main_dir="/project/zerze/morganmurray/software/RPA_ddft3_2"
#export data_dir=${main_dar}
#export PYTHONPATH=${main_dar}
epoch_time=$(date +%s)
unsortedjob="unsorted_$epoch_time"


jobname=${1:-$unsortedjob}
jobname=Ddx4_u_40_3D
module purge
module load Miniforge3/py3.10
#echo "${main_dir}/env"
source activate ${main_dir}/env

#ml Trilinos/13.4.1-foss-2022a

cd $main_dir
mkdir -p ${main_dir}/Results/$jobname
mkdir -p ${main_dir}/Results/$jobname/states


python  scripts/DDFT/Protein_RPA/calc_kappa.py Ddx4_N1 3.5 0.15 -0.3 100
#python  scripts/DDFT/movie.py --i ${main_dir}/Results/$jobname --o ${main_dir}/Results/$jobname
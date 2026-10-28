#!/bin/bash
#SBATCH -J ddft_logplot
#SBATCH -o ddft_logplot.o%j
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks-per-node=1 -N 1
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=32GB

#module --ignore_cache load "OpenMPI/.4.1.4-GCC-11.3.0"
export OMP_NUM_THREADS=1

export main_dir="/project/zerze/morganmurray/software/RPA_ddft3_2"
#export data_dir=${main_dar}
#export PYTHONPATH=${main_dar}
epoch_time=$(date +%s)
unsortedjob="unsorted_$epoch_time"


jobname=${1:-$unsortedjob}
jobname=EKSCD0667_realKappaZeroFLux
module purge
module load Miniforge3/py3.10
#echo "${main_dir}/env"
source activate ${main_dir}/env

#ml Trilinos/13.4.1-foss-2022a

cd $main_dir
mkdir -p ${main_dir}/Results/$jobname
mkdir -p ${main_dir}/Results/$jobname/states


mpirun -np 1 python  scripts/DDFT/plot_state.py --i ${main_dir}/scripts/DDFT/input_params.txt --o ${main_dir}/Results/$jobname/  --s ${main_dir}/Results/$jobname/states
#python  scripts/DDFT/movie.py --i ${data_dir}/Results/$jobname --o ${data_dir}/Results/$jobname

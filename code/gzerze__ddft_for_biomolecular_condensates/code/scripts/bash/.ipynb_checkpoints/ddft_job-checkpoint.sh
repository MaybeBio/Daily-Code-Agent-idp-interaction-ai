#!/bin/bash
#SBATCH -J ddft_job
#SBATCH -o ddft_job.o%j
#SBATCH -t 48:60:20
#SBATCH -n 16
#SBATCH --mem-per-cpu=32GB 

#module --ignore_cache load "OpenMPI/.4.1.4-GCC-11.3.0"
export OMP_NUM_THREADS=1

export main_dir="/project/zerze/morganmurray/software/RPA_ddft3_2"
#export data_dir=${main_dar}
#export PYTHONPATH=${main_dar}

module load Miniforge3/py3.10
#echo "${main_dir}/env"
source activate ${main_dir}/env

#ml Trilinos/13.4.1-foss-2022a

cd $main_dir
mkdir -p ${main_dir}/Results/states

mpirun -np 8 python  scripts/DDFT/RPA_ddft.py --i ${main_dir}/scripts/DDFT/input_params.txt --o ${main_dir}/Results/  --s ${main_dir}/Results/
#mpirun -np 1 python  scripts/DDFT/movie.py --i ${data_dir}/Results/ --o ${data_dir}/Results/

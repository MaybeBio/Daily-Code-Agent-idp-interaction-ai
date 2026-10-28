#!/bin/bash
#SBATCH -J ddft_job
#SBATCH -o ddft_job.o%j
#SBATCH -t 48:60:20
#SBATCH -n 8
#SBATCH --mem-per-cpu=16GB 

#module --ignore_cache load "OpenMPI/.4.1.4-GCC-11.3.0"
export OMP_NUM_THREADS=1

export main_dir="/project/zerze/morganmurray/software/RPA_ddft3_2"
export data_dir=${main_dar}
export PYTHONPATH=${main_dar}

module 

#ml Trilinos/13.4.1-foss-2022a

cd $main_dir
mkdir -p ${data_dir}/output_8/states

mpirun -np 8 python  scripts/DDFT/RPA_ddft.py --i ${data_dir}/scripts/DDFT/input_params.txt --o ${data_dir}/output_{i}/output_8/  --s ${data_dir}/output_8

#!/bin/bash
#SBATCH -J ddft_movie
#SBATCH -o movie_output.o%j
#SBATCH -t 00:60:20
#SBATCH --mem-per-cpu=32GB 


jobname=${1:-$unsortedjob}
jobname='EKSCD0200_zeroflux_3_0'
export main_dir="/project/zerze/morganmurray/software/RPA_ddft3_2"

#export data_dir=${main_dar}
#export PYTHONPATH=${main_dar}

module load Miniforge3/py3.10
#echo "${main_dir}/env"
source activate ${main_dir}/env

#ml Trilinos/13.4.1-foss-2022a

cd $main_dir
mkdir -p ${main_dir}/Results//${jobname}/movie


python  scripts/DDFT/movie.py --i ${main_dir}/Results/${jobname}/ --o ${main_dir}/Results/${jobname}/movie/
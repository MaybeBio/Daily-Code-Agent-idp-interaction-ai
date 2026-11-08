#!/bin/bash
#$ -S /bin/bash
#$ -l h_rt=02:00:00
#$ -l node_o=1
#$ -cwd

module load  openmpi/5.0.2-gcc
module load cuda/12.1.0
export OMP_NUM_THREADS=24
export LD_LIBRARY_PATH=/gs/bs/hp170025/application/plumed/lib:${LD_LIBRARY_PATH}
GMX=/gs/bs/hp170025/application/gromacs-2022.6-plumed/bin
############################################################

# [1] Continue from the previous iteration
${GMX}/gmx_mpi grompp  -f ../run.mdp -c  ../../../npt_1000_Parrinello-Rahman/n1/posre.gro  -p ../../../top_1000/md.top  -o run.tpr

# [2] Continue from the previous iteration
##${GMX}/gmx_mpi grompp  -f ../run.mdp -c  ../../../init/3/n1/run.gro  -p ../../../sys/complex_P_02.top  -o run.tpr

# [3] Re-setting the initial configurations from the entire ensemble
##${GMX}/gmx_mpi grompp  -f ../run.mdp -c  ../../0/next_conf/n1/run.gro  -p ../../../sys/complex_P_02.top  -o run.tpr


${GMX}/gmx_mpi mdrun  -deffnm run  -plumed ../plumed.dat

#!/bin/bash
#$ -S /bin/bash
#$ -l h_rt=02:00:00
#$ -l node_o=1
#$ -cwd

module load  openmpi/5.0.2-gcc
# New plumd procedure.
export OMP_NUM_THREADS=24
export LD_LIBRARY_PATH=/gs/bs/hp170025/application/plumed/lib:${LD_LIBRARY_PATH}
GMX=/gs/bs/hp170025/application/gromacs-2022.6-plumed/bin
############################################################

## Minimization
${GMX}/gmx_mpi grompp -f ../../basic/minim.mdp -c complex_P_md_n1.gro -p complex_P_before_md1.top -o minim.tpr
${GMX}/gmx_mpi mdrun -deffnm minim 

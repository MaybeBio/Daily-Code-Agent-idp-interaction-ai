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

## Posre
${GMX}/gmx_mpi grompp -f ../../basic/posre_Parrinello-Rahman.mdp -c ../../npt_1000_Berendsen/n1/posre.gro -p ../../top_1000/complex_P_before_md1.top -r ../../npt_1000_Berendsen/n1/posre.gro -o posre.tpr -maxwarn 1
${GMX}/gmx_mpi mdrun -deffnm posre 

#!/bin/bash

# [3] Iteration 3 (Re-setting the initial configurations from the entire ensemble)
${GMX}/gmx_mpi grompp  -f run.mdp -c  iteration2/next_conf/run.gro  -p ../parameters/md.top  -o run.tpr

${GMX}/gmx_mpi mdrun  -deffnm run  -plumed plumed.dat

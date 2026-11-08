#!/bin/bash

# [2] Iteration 2 (Continue from the previous iteration)
${GMX}/gmx_mpi grompp  -f run.mdp -c  ../iteration1/run.gro  -p ../parameters/md.top  -o run.tpr

${GMX}/gmx_mpi mdrun  -deffnm run  -plumed plumed.dat

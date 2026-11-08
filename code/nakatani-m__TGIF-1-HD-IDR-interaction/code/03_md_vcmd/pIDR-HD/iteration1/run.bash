#!/bin/bash

# [1] NVT ensemble
${GMX}/gmx_mpi grompp  -f run.mdp -c  posre.gro  -p ../parameters/md.top  -o run.tpr

${GMX}/gmx_mpi mdrun  -deffnm run  -plumed plumed.dat  



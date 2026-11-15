#!/bin/bash
#$ -cwd
#$ -j y
#$ -v

SaveFolder=$HOME/demoZ
dt=dt001
N=N512
L=L4
SAVE=s200


if [ -d "$SaveFolder" ]; then
	echo 'saving file in $SaveFolder';
else
	echo 'Folder' $SaveFolder 'does not exist';
	echo 'Creating a New Folder First';
	mkdir $SaveFolder
fi

bsub -q gpu -o ${SaveFolder}/MYGPUJOB.o%J ./LLPS ${SaveFolder} -s ${dt},$N,$L,${SAVE} >> ${SaveFolder}/out.txt

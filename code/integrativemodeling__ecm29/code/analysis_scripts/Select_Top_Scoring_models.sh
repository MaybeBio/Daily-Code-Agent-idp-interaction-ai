#!/bin/bash
#$ -l arch=linux-x64 
#$ -l h_rt=18:0:0
#$ -l mem_free=2G
#$ -S /bin/bash 
#$ -N pref
#$ -cwd 
#$ -j y 
#$ -R y
#$ -t 1-500


i=$(($SGE_TASK_ID - 1))

source ~/.bash_profile
#python prefilter.py -prefilter 70 -nmods 10000 -dir /scrapp/ilan/ECM29/19S/Running_Directory_4/modeling${i}

process_output.py -f /scrapp/ilan/ECM29/19S/Running_Directory_4/modeling${i}/output/stat.0.out -s Total_Score CrossLinkingMassSpectrometryRestraint_Data_Score_Lan ExcludedVolumeSphere_all ExcludedVolumeSphere_ecm29 --nframe | awk '{print $1, $3, $4, $5, $6}' > Scores_${i}.txt
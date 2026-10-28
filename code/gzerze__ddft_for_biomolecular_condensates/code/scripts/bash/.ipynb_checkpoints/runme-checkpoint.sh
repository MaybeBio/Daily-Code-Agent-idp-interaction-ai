#!/bin/bash

#rm ddft*

#scancel -u mgmurray
#rm ddft*


for i in 8;do
mkdir -p ./Results/output_${i}
cat ./scripts/bash/template.sh | sed -e"s/nproc/${i}/g"|sed -e"s/Results/output_${i}/g" > ./scripts/bash/run_prog_${i}.sh
echo "hello"
sbatch ./scripts/bash/run_prog_${i}.sh
done


tail -f ddft*


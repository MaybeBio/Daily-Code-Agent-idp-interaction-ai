#!/bin/bash

submit () {
   sbr="$(/usr/bin/sbatch "$@")"
   wait
   if [[ "$sbr" =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
       echo "${BASH_REMATCH[1]}"
   else
       echo "sbatch failed"
       exit 1 
   fi
}


if [ $# == 3 ]; then
   jid=$3
   echo "dependency: $jid"
   jid=`submit --dependency=afterany:$jid $1`
else
   echo "no dependency"
   jid=`submit $1`
fi
wait
echo "SUBMITTED JOBID: $jid"

for ((i=1; i<=$2; i++))
do
    jid="$(submit --dependency=afterany:$jid $1)"
    wait
    echo "SUBMITTED JOBID: $jid"
done

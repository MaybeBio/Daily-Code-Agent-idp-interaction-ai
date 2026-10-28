#!/bin/bash

for i in `seq 10000 10000 100000`;
do

    for j in `seq 1 50`; 
    do
	shuf $1  | awk '{print $1}' | shuf | shuf | head -${i} | sort -n  | head -1 >> tmp_1;
    done

    avg1=$(awk '{sum+=$1; b++}END{print sum/b}' tmp_1)
    std1=$(awk '{sum+=$1; sum2+=$1*$1; b++}END{print sqrt(sum2/b-(sum/b)*(sum/b))}' tmp_1)

    echo $i $avg1 $std1
    rm tmp_1

done

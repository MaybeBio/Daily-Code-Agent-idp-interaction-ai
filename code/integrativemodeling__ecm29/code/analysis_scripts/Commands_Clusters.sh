#!/usr/bin/bash

mkdir cluster.0 cluster.1
while read p; do grep " $p 00" Identities2.txt; done < cluster.0.all.txt | awk '{print $1}' > cluster.0/Models.txt
while read p; do grep " $p 00" Identities2.txt; done < cluster.1.all.txt | awk '{print $1}' > cluster.1/Models.txt

for i in `seq 0 1`; do while read p; do cp $p cluster.$i/ ; done < cluster.$i/Models.txt; done



#!/bin/bash

DATE=20261014
EXP=run3/s

for S in 300; do
    DIR="${EXP}${S}"
    echo "Pulling $DIR ..."
    mkdir -p "data/${DATE}/${DIR}"
    scp -r "tohme@thanos.phys.ens.fr:/home/tohme/.julia/dev/Droplet/data/${DATE}/${DIR}/*" \
        "data/${DATE}/${DIR}"
done

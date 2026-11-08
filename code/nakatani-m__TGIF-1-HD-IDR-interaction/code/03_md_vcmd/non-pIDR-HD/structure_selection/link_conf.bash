#!/bin/bash

#. /etc/profile.d/modules.sh
#module load python/3.8.3

python link_conf_dim3.py -s $1 -n $2 -inp inp  > st${1}_next_conf.bash

bash st${1}_next_conf.bash 

#mv st${1}_next_conf.bash  ./bash_files/

rm st${1}_next_conf.bash 


#!/bin/sh
#06 packmol
for i in `seq 1 1 1000`
do


Random=`LC_CTYPE=C tr -dc 0-9 < /dev/urandom | fold -w 6 | head -1`

echo "${i} ${Random}"
cat <<- EOF > input/packmol_n${i}.in
seed ${Random}
tolerance 2.0
filetype pdb
output PDB/6fqp_P_in_water_n${i}.pdb

#discale=1.5
#maxit=20
#avoid_overlap=yes

structure 6fqp_after_tleap.pdb
  number 1
  center
  fixed 51.8 51.8 51.8 40. 40. 40.
end structure

structure IDP_P_Min.pdb
  number 1
  center
  #inside box 0. 0. 0. 103.6/ 103.6/ 103.6/
  outside sphere 51.8 51.8 51.8 30.0
  inside  sphere 51.8 51.8 51.8 51.8
end structure

structure spc.pdb
  number 40000
  inside box 0. 0. 0. 103.6/ 103.6/ 103.6/
end structure
EOF

done 

count=1
for i in `seq 1 1 1000`
do
if [ $count -eq 8 ]; then
	count=0
	wait
fi
count=$((count + 1))

echo "${i} ${count}"
#sleep 10 & 
packmol <  input/packmol_n${i}.in > LOG/log_n${i}.txt &


done




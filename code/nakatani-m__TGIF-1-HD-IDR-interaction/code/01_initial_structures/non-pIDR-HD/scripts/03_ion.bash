#!/bin/bash
for i in {1..1000}; do
#01 make box
gmx editconf -f PDBS/6fqp_nonP_in_water_n${i}.pdb -o GROS/complex_nonP_in_water_n${i}.gro -box 10.36 10.36 10.36

#02 input water file
cat complex_nonP_edit_01.top |
awk '{
       if ( substr($0,1,16) == "[ moleculetype ]" && inserted == 0 ){
               printf "; Include atomtypes for water and ions\n"
               printf "\#include \"./water_and_ions_nonbonded.itp\"\n"
               printf "\n"
               print $0
	       inserted = 1
       }else if(substr($0,1,10) == "[ system ]" ){
               printf "; Include topology for water\n"
               printf "\#include \"./spce.itp\"\n"
               printf "; Include topology for ions\n"
               printf "\#include \"./ions.itp\"\n"
               printf "\n"
               print $0
       }else{
               print $0
       }
}END{
       printf "SOL       40000\n"
}' > top/tmp_complex_nonP_before_md${i}.top

cp top/tmp_complex_nonP_before_md${i}.top top/complex_nonP_before_md${i}.top

#03 input ions
gmx grompp -f ions.mdp -c GROS/complex_nonP_in_water_n${i}.gro -p top/complex_nonP_before_md${i}.top -o ions${i}.tpr
gmx genion -s ions${i}.tpr -o GROS2/complex_nonP_md_n${i}.gro -p top/complex_nonP_before_md${i}.top -conc 0.15 -neutral yes << EOF
13
EOF

done


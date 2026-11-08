#!/bin/sh
##HD
#01 get pdb 
cat 6fqp.pdb | grep "ATOM  " | grep -v "REMARK" | grep " A " > 6fqp_edit_01.pdb

#02 modeller
cat rcsb_pdb_6FQP.fasta | sed -n 4p |
awk 'BEGIN{
        printf ">P1;6fqp_model\n"
                printf "sequence:6fqp_model::::::::\n"
}{
        gsub("GPMDIPLDLSSSAGSGKRRR","KRRR",$0)
        gsub("NQFTISRRGAKISE","",$0)
        printf "%s*\n", $0
}' > unknownStructure.pir

echo " " >> unknownStructure.pir

cat rcsb_pdb_6FQP.fasta | sed -n 4p |
awk 'BEGIN{
        printf ">P1;REFERENCE\n"
        printf "structure:6fqp_modeller_inp:5:A:71:A::::\n"
}{
        gsub("GPMDIPLDLSSSAGSGKRRR","----",$0)
        gsub("NQFTISRRGAKISE","",$0)
        printf "%s*\n", $0
}' >> unknownStructure.pir



cat 6fqp_edit_01.pdb |
awk 'BEGIN{
        nai_start = 4
        bef_num=0
}{
        num=substr($0,23,4)
        if(num != bef_num ){nai_start+=1}
        #print bef_num
        #print num
        #print $0
        printf "%s%4d%s\n",substr($0,1,22),nai_start,substr($0,27,54)
        bef_num=num

}' > 6fqp_modeller_inp.pdb




cat <<- EOF > use_modeller.py
#!/usr/bin/python

from modeller import *
from modeller.automodel import *

class MyModel(automodel):
        def special_restraints(self, aln):
                rsr = self.restraints
                at = self.atoms
                rsr.rigid_bodies.append(rigid_body(
                selection(a.residue_range('5:A', '71:A')) ))


env = environ()

a = MyModel(env,
              alnfile        = 'unknownStructure.pir',
              knowns         = 'REFERENCE',
              sequence       = '6fqp_model',
              assess_methods = (assess.DOPE, assess.GA341)
             )

a.starting_model = 1
a.ending_model   = 10
a.make()
EOF

mod10.5 use_modeller.py


selectd_model=`cat use_modeller.log | grep pdb | sort -k 3n | head -n 1 | awk '{print $1}'`
echo ${selectd_model}

cat ${selectd_model} | grep "ATOM  " | grep -v "OXT" |
awk 'BEGIN{
        no_structure = 163
        std_num=0 
}{
        residue=substr($0,23,4)
        if(residue != std_num ){no_structure+=1}
        #print std_num
        #print residue
        #print $0
        printf "%s%4d%s\n",substr($0,1,22),no_structure,substr($0,27,54)
        std_num = residue

}' > 6fqp_modeller_oup.pdb

cat <<- EOF > cap.py

from pymol import editor, cmd

cmd.load('6fqp_modeller_oup.pdb')

cmd.select('N_term',  '(chain A and resi 164) and name N ')
editor.attach_amino_acid('N_term', 'ACE')

cmd.select('C_term',  '(chain A and resi 234) and name C')
editor.attach_amino_acid('C_term', 'NME')


cmd.save('cappedchain.pdb')

EOF

/Applications/PyMOL.app/Contents/MacOS/python cap.py

cat cappedchain.pdb | grep -v "TER" > cappedchain_02.pdb

#03 tleap
cat <<- EOF > tleap.in
source /Users/madoka/Apps/AMBER/amber22/dat/leap/cmd/oldff/leaprc.ff99SB
chainA = loadpdb  cappedchain_02.pdb
saveamberparm chainA 6fqp_after_tleap.prmtop 6fqp_after_tleap.inpcrd 
quit
EOF

tleap -f tleap.in

#04 output file
cat <<- EOF > cpptraj.in
parm    6fqp_after_tleap.prmtop
trajin  6fqp_after_tleap.inpcrd 
trajout 6fqp_after_tleap.pdb pdb
go
EOF

cpptraj -i cpptraj.in


##IDR
#05 fasta
cat  <<- EOF > IDP_nonP.fasta
TLGRPLSPKPSSPGSVLARPS
EOF

#06 rename
SEQ=`cat IDP_nonP.fasta | fold -w 1 |
awk 'BEGIN{printf "ACE "}{
 if(substr($0,1,1) == "A") {gsub("A","ALA",$0)}
 else if(substr($0,1,1) == "T") {gsub("T","THR",$0)}
 else if(substr($0,1,1) == "L") {gsub("L","LEU",$0)}
 else if(substr($0,1,1) == "G") {gsub("G","GLY",$0)}
 else if(substr($0,1,1) == "R") {gsub("R","ARG",$0)}
 else if(substr($0,1,1) == "P") {gsub("P","PRO",$0)}
 else if(substr($0,1,1) == "S") {gsub("S","SER",$0)}
 else if(substr($0,1,1) == "K") {gsub("K","LYS",$0)}
 else if(substr($0,1,1) == "V") {gsub("V","VAL",$0)}
printf "%s ",$0
#printf "%d %s \n",NR,$0
}END{printf "NME\n"}'`
echo ${SEQ}

#07 tleap
cat <<- EOF > tleap.in
source ${AMBERHOME}/dat/leap/cmd/leaprc.protein.ff14SB
source ${AMBERHOME}/dat/leap/cmd/leaprc.phosaa14SB
chainA = sequence { ${SEQ} }
saveamberparm chainA IDP_nonP_after_tleap.prmtop IDP_nonP_after_tleap.inpcrd
savepdb chainA IDP_nonP_after_tleap.pdb
quit
EOF

tleap -f tleap.in

#08 minimize
cat <<- EOF > Min.in
Minimize
 &cntrl
  imin=1,
  ntx=1,
  irest=0,
  maxcyc=2000,
  ncyc=1000,
  ntpr=100,
  ntwx=0,
  igb=6
  cut=999.9,
 /
EOF

sander -O -i Min.in -o Min.out -p IDP_nonP_after_tleap.prmtop -c IDP_nonP_after_tleap.inpcrd -r IDP_nonP_Min.rst -inf Min.mdinfo

#09 binary to pdb
cat <<- EOF > cpptraj.in
parm    IDP_nonP_after_tleap.prmtop
trajin  IDP_nonP_Min.rst
trajout IDP_nonP_Min.pdb pdb
go
EOF

cpptraj -i cpptraj.in


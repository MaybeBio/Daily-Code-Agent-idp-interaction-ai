mol new step3_input.pdb

# Move lower corner to 0,0,0
lassign [measure minmax [atomselect top all]] v1 v2
[atomselect top all] moveby [vecinv $v1]  ;# <-- Fixed closing bracket
[atomselect top all] writepdb step3_input.pdb
[atomselect top all] writerst7 step3_input.rst7

# Write the restraints file
set fp [open "restraints/prot_pos.txt" w]  ;# <-- Fixed file path format (removed leading /)
set bb [[atomselect top "noh and backbone"] get index]
foreach b $bb {
    puts $fp "$b BB"
}
set sc [[atomselect top "sidechain and noh"] get index]
foreach s $sc {
    puts $fp "$s SC"
}
close $fp

quit

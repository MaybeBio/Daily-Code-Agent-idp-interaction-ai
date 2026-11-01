#################### Inputs ####################
# 0 --> psf file of the system to remove hydrogens and extract only the protein
# 1 --> pdb location to be saved
# 2 --> ligand's resname

# Check if argv has at least 3 elements
if {[llength $argv] > 2} {
    set resname [lindex $argv 2]
    set seltext "protein and noh and not resname $resname"
} else {
    # No argv2 provided → use alternate selection, there is no ligand to unselect
    set seltext "protein and noh"
}

mol new [lindex $argv 0] type pdb
[atomselect top $seltext] writepdb [lindex $argv 1]
quit

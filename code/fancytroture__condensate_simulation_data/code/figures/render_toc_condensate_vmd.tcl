# Render a tagged chain within its real local condensate environment.
# Usage: vmd -dispdev text -e render_toc_condensate_vmd.tcl -args input.pdb scene.dat color_id

set input_pdb [lindex $argv 0]
set output_scene [lindex $argv 1]
set color_id [lindex $argv 2]

display projection Orthographic
display depthcue off
axes location Off
color Display Background white
color change rgb 7 0.184 0.498 0.451
color change rgb 1 0.710 0.294 0.373
color change rgb 6 0.68 0.71 0.73

mol new $input_pdb type pdb waitfor all
mol delrep 0 top

# Actual neighboring chains: translucent bead chains in the condensate.
mol representation Bonds 0.13 16
mol color ColorID 6
mol material Transparent
mol selection "not segname TARG"
mol addrep top
mol representation VDW 0.44 16
mol color ColorID 6
mol material Transparent
mol selection "not segname TARG"
mol addrep top

# Tagged chain: saturated pearl-necklace representation.
mol representation Bonds 0.30 20
mol color ColorID $color_id
mol material AOShiny
mol selection "segname TARG"
mol addrep top
mol representation VDW 1.05 20
mol color ColorID $color_id
mol material AOShiny
mol selection "segname TARG"
mol addrep top

display resetview
rotate x by -8
rotate y by 10
scale by 1.08
if {[catch {render Tachyon $output_scene} message]} {
    puts "RENDER ERROR: $message"
}
quit

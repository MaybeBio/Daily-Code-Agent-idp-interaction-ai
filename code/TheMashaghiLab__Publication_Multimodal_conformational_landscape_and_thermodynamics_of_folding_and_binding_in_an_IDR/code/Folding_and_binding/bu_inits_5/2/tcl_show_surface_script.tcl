# === Clear previous drawings ===
graphics top delete all

# === PARAMETERS ===
set a 0.01            ;# Profile: r(z) = a * z^2

set h 18.0
set f 3.0
set b 0.25
set sc 30.0

set zmin 5.0
set zmax 41.0
set nz 150
set ntheta 60

# === SELECT TWO ATOMS TO DEFINE THE AXIS ===
set atom1 [atomselect top "serial 117"]

set comsel [atomselect top "serial 326 707"]
set atom2 [atomselect top "serial 326"]
set com [measure center $comsel weight mass]

foreach {x1 y1 z1} [lindex [$atom1 get {x y z}] 0] break
foreach {x2 y2 z2} [lindex [$atom2 get {x y z}] 0] break
foreach {x2 y2 z2} $com break

# === Compute axis vector and rotation matrix ===
set vx [expr {$x2 - $x1}]
set vy [expr {$y2 - $y1}]
set vz [expr {$z2 - $z1}]
set vlen [expr {sqrt($vx*$vx + $vy*$vy + $vz*$vz)}]

# Normalize
set ux [expr {$vx / $vlen}]
set uy [expr {$vy / $vlen}]
set uz [expr {$vz / $vlen}]

# Now create rotation matrix to rotate z-axis to this vector
# This is a rotation from vector (0,0,1) → (ux, uy, uz)
# Based on Rodrigues’ rotation formula

# Compute cross product (z-axis × u)
set rx [expr {0.0 * $uz - 1.0 * $uy}]
set ry [expr {1.0 * $ux - 0.0 * $uz}]
set rz [expr {0.0 * $uy - 0.0 * $ux}]
set rlen [expr {sqrt($rx*$rx + $ry*$ry + $rz*$rz)}]

if {$rlen > 1e-6} {
    # Normalize rotation axis
    set rx [expr {$rx / $rlen}]
    set ry [expr {$ry / $rlen}]
    set rz [expr {$rz / $rlen}]
    
    # Angle between vectors
    set dot [expr {$uz}]
    set angle [expr {acos($dot)}]

    # Precompute Rodrigues terms
    set c [expr {cos($angle)}]
    set s [expr {sin($angle)}]
    set t [expr {1.0 - $c}]
    
    # Rotation matrix R (3x3)
    set R [list \
        [expr {$t*$rx*$rx + $c}]       [expr {$t*$rx*$ry - $s*$rz}]   [expr {$t*$rx*$rz + $s*$ry}] \
        [expr {$t*$ry*$rx + $s*$rz}]   [expr {$t*$ry*$ry + $c}]       [expr {$t*$ry*$rz - $s*$rx}] \
        [expr {$t*$rz*$rx - $s*$ry}]   [expr {$t*$rz*$ry + $s*$rx}]   [expr {$t*$rz*$rz + $c}]     \
    ]
} else {
    # Vectors are parallel; identity matrix
    set R {1 0 0 0 1 0 0 0 1}
}

# === FUNCTION: apply rotation matrix to vector ===
proc rotate {R vec} {
    set x [lindex $vec 0]
    set y [lindex $vec 1]
    set z [lindex $vec 2]
    return [list \
        [expr {[lindex $R 0]*$x + [lindex $R 1]*$y + [lindex $R 2]*$z}] \
        [expr {[lindex $R 3]*$x + [lindex $R 4]*$y + [lindex $R 5]*$z}] \
        [expr {[lindex $R 6]*$x + [lindex $R 7]*$y + [lindex $R 8]*$z}] \
    ]
}

# === FUNCTION: translate vector ===
proc translate {vec dx dy dz} {
    return [list \
        [expr {[lindex $vec 0] + $dx}] \
        [expr {[lindex $vec 1] + $dy}] \
        [expr {[lindex $vec 2] + $dz}] \
    ]
}

# === BUILD THE SURFACE IN Z-AXIS COORDINATES ===
set dz [expr {($zmax - $zmin) / double($nz)}]
set dtheta [expr {2.0 * acos(-1.0) / double($ntheta)}]

for {set iz 0} {$iz < $nz} {incr iz} {
    set z0 [expr {$zmin + $dz * $iz}]
    set z1_loc [expr {$zmin + $dz * ($iz+1)}]
    set r0 [expr {$h*(1./(1.+exp($b*($z0-$sc))))+$f}]
    set r1 [expr {$h*(1./(1.+exp($b*($z1_loc-$sc))))+$f}]

    for {set itheta 0} {$itheta < $ntheta} {incr itheta} {
        set t0 [expr {$itheta * $dtheta}]
        set t1 [expr {($itheta + 1) * $dtheta}]
        
        # Points in local Z-axis frame
        set p00 [list [expr {$r0 * cos($t0)}] [expr {$r0 * sin($t0)}] $z0]
        set p01 [list [expr {$r0 * cos($t1)}] [expr {$r0 * sin($t1)}] $z0]
        set p10 [list [expr {$r1 * cos($t0)}] [expr {$r1 * sin($t0)}] $z1_loc]
        set p11 [list [expr {$r1 * cos($t1)}] [expr {$r1 * sin($t1)}] $z1_loc]

        # Rotate all points to align with axis
        set p00r [rotate $R $p00]
        set p01r [rotate $R $p01]
        set p10r [rotate $R $p10]
        set p11r [rotate $R $p11]

        # Translate to atom1 position
        set p00t [translate $p00r $x1 $y1 $z1]
        set p01t [translate $p01r $x1 $y1 $z1]
        set p10t [translate $p10r $x1 $y1 $z1]
        set p11t [translate $p11r $x1 $y1 $z1]

        # Draw triangles
        graphics top triangle $p00t $p01t $p10t
        graphics top triangle $p01t $p11t $p10t
    }
}

graphics top material Transparent

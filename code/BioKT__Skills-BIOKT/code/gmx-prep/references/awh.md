# AWH Free Energy (COM distance)

Build on NPT-equilibrated system. Requires index groups for the two molecules/domains.

**Additional mdp sections** (add to npt.mdp base):

```ini
; Pull code — COM distance as AWH external-potential coordinate
pull                     = yes
pull-ngroups             = 2
pull-group1-name         = MOL_A       ; must exist in index.ndx
pull-group2-name         = MOL_B
pull-ncoords             = 1
pull-coord1-type         = external-potential
pull-coord1-potential-provider = AWH
pull-coord1-geometry     = distance
pull-coord1-groups       = 1 2
pull-coord1-start        = no
pull-coord1-dim          = Y Y Y       ; 3-D scalar distance

; AWH
awh                      = yes
awh-potential            = convolved
awh-share-multisim       = no
awh-seed                 = -1
awh-nstout               = 5000        ; write bias every 10 ps (multiple of nstenergy)

; AWH bias 1
awh-nbias                = 1
awh1-error-init          = 10.0        ; kJ/mol
awh1-growth              = exp-linear
awh1-equilibrate-histogram = yes
awh1-target              = constant
awh1-ndim                = 1

; Dimension 1
awh1-dim1-coord-index    = 1
awh1-dim1-start          = START_NM    ; nm — closest approach
awh1-dim1-end            = END_NM      ; nm — fully separated
awh1-dim1-force-constant = 4000        ; kJ/mol/nm²
awh1-dim1-diffusion      = 5e-5        ; nm²/ps (initial estimate)
```

**Index group requirement:**
```bash
gmx make_ndx -f prep/system.gro -o prep/index.ndx << 'EOF'
r FIRST_RESID-LAST_RESID_A
name N MOL_A
r FIRST_RESID-LAST_RESID_B
name N+1 MOL_B
q
EOF
```

**AWH run command:**
```bash
gmx grompp -f mdp/awh.mdp -c runs/npt/npt.gro -t runs/npt/npt.cpt \
    -p prep/topol.top -n prep/index.ndx -o runs/awh/awh.tpr -maxwarn 2
gmx mdrun -v -deffnm runs/awh/awh -ntmpi 1 -ntomp 4 -pme gpu
```

**Extract PMF:**
```bash
gmx awh -s runs/awh/awh.tpr -f runs/awh/awh.edr -o analysis/pmf.xvg -more
```

# H-REMD (Hamiltonian Replica Exchange)

Uses free-energy perturbation with vdw_lambdas to scale water–solute interactions.

**Rep-0 mdp template** (`mdp/hremd_rep0.mdp`):

```ini
integrator               = sd
dt                       = 0.002
nsteps                   = 5000000      ; 10 ns per replica
nstlog                   = 1000
nstenergy                = 1000
nstxout-compressed       = 1000
compressed-x-grps        = non-Water
continuation             = yes
constraint_algorithm     = lincs
constraints              = h-bonds
lincs-iter               = 1
lincs-order              = 4
cutoff-scheme            = Verlet
nstlist                  = 10
coulombtype              = PME
rcoulomb                 = 1.0
rvdw                     = 1.0
DispCorr                 = EnerPres
tc-grps                  = Protein Non-Protein
tau_t                    = 1.0    1.0
ref_t                    = 300    300
pcoupl                   = no
gen_vel                  = no

; H-REMD settings
free-energy              = yes
init-lambda-state        = 0            ; replica index (0-based)
nstdhdl                  = 1000
vdw-lambdas              = 0.00 0.15 0.30 0.45 0.60 0.75 0.90 1.00
```

**Generate per-replica mdp files:**
```bash
NREP=8
for i in $(seq 1 $((NREP-1))); do
    awk -v i=$i '/init-lambda-state/{print "init-lambda-state        = " i; next}{print}' \
        mdp/hremd_rep0.mdp > mdp/hremd_rep${i}.mdp
done
```

**Prepare tpr files (one per replica):**
```bash
for r in $(seq 0 $((NREP-1))); do
    gmx grompp -f mdp/hremd_rep${r}.mdp \
        -c runs/npt/npt.gro -t runs/npt/npt.cpt \
        -p prep/topol.top -n prep/index.ndx \
        -o runs/hremd/hremd_rep${r}.tpr -maxwarn 1
done
```

**Run (MPI required):**
```bash
mpirun -np ${NREP} gmx_mpi mdrun \
    -s runs/hremd/hremd_rep \
    -deffnm runs/hremd/hremd_rep \
    -multi ${NREP} \
    -replex 1000 \
    -ntomp 4
```

**Demux trajectories** (after run):
```bash
# requires demux.pl from GROMACS tools or gmx trjcat
demux.pl runs/hremd/hremd_rep0.log
```

# Umbrella Sampling (GROMACS pull code)

For PLUMED-based umbrella sampling (RESTRAINT inputs, WHAM with error bars), use the `plumed-us` skill instead.

Three phases: pulling → window extraction → WHAM.

**Phase 1 — Steered MD (generate windows):**

```ini
; pull_steer.mdp — append to npt.mdp base
pull                     = yes
pull-ngroups             = 2
pull-group1-name         = REFERENCE
pull-group2-name         = LIGAND
pull-ncoords             = 1
pull-coord1-type         = umbrella
pull-coord1-geometry     = distance
pull-coord1-groups       = 1 2
pull-coord1-rate         = 0.001        ; nm/ps (slow pull)
pull-coord1-k            = 1000         ; kJ/mol/nm²
pull-coord1-start        = yes
```

**Phase 2 — Window simulations:**
```bash
# Extract one frame per window using gmx trjconv
# Run each window with harmonic restraint at target COM distance
# pull-coord1-rate = 0 ; pull-coord1-init = TARGET_NM
```

**Phase 3 — WHAM:**
```bash
gmx wham -it tpr-files.dat -if pullf-files.dat -o analysis/pmf_wham.xvg \
    -hist analysis/histograms.xvg -unit kJ
```

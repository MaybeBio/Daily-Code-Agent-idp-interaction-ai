# Condensate-only analysis

This folder contains analysis scripts for the condensate-only simulations, corresponding to Extended Data Fig. 2 in the manuscript.

Panels a-c were analyzed using GROMACS built-in tools:

**a) Monomer RMSD**

```bash
gmx_mpi rms -f protein_fit.xtc -s protein.tpr -n index_md.ndx -o rmsd.xvg
```

**b) Maximum cluster size in multichain simulations**

```bash
gmx_mpi clustsize -f protein_only.xtc -s protein.tpr -n index_md.ndx -cut 0.528 -mc maxcluster.xvg -nc ncluster.xvg -mcn maxcluster.ndx
```

With BB beads selected, the largest-cluster size in `maxcluster.xvg` is a BB-bead count. Divide this value by the number of residues per chain to express the cluster size in chain equivalents.

**c) SASA in multichain simulations**

```bash
gmx_mpi sasa -f protein_only.xtc -s protein.tpr -probe 0.191 -o sasa.xvg
```

**d-f) Interchain contact analysis**

Panels d-f are based on interchain contact statistics. Keep `vdw_def.py` in the same directory as `Inter Contact-map.py`.

Running analysis code need MDAnalysis.(version 2.10.0 we used)
For MDAnalysis installation instructions, see the [MDAnalysis installation guide](https://www.mdanalysis.org/pages/installation_quick_start/).

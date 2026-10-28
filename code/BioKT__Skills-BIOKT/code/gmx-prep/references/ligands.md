# Ligands / Custom Molecules (ACPYPE → GROMACS)

For parameterising the ligand itself (charges, atom types, validation), use the `amber-parameterize` skill; this file covers bringing the result into a GROMACS system.

```bash
# 1. Parametrise with ACPYPE (GAFF2 force field)
acpype -i ligand.mol2 -c bcc -n 0 -a gaff2

# 2. Copy generated files
cp ligand.acpype/ligand_GMX.itp prep/ligand.itp
cp ligand.acpype/ligand_GMX.gro prep/ligand.gro

# 3. In topol.top — add before [ system ]
#include "ligand.itp"

# 4. In [ molecules ] section, add ligand line
# LIG   1

# 5. Combine coordinates
gmx insert-molecules -f prep/protein_processed.gro \
    -ci prep/ligand.gro -nmol 1 -o prep/complex.gro -radius 0.3

# 6. Continue with editconf → solvate → genion as normal
```

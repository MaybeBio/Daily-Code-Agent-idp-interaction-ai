## The Bug: Pre-Existing deltaHam! Cross-Species Double-Counting for Same-Polymer PolyMonomers
### Root Cause

There are two intertwined pre-existing bugs in deltaHam! for heterogeneous polymers (polymers with multiple monomer types). They both stem from the same fundamental design gap: the energy accounting for cross-species interactions between PolyMonomers on the same polymer was never correctly implemented.

### Bug 1: Wrong nnp dispatch
In deltaHam!, cross-species energy for same-polymer PolyMonomers uses:

This dispatches to nnp(AbstractSpecies, AbstractSpecies, pos) which counts all occupied neighbor sites (latS > 0), regardless of monomer type. But ham! uses:

This dispatches to nnp(AbstractPolyMonomer, AbstractPolyMonomer, pos) which counts only sites matching target.intIdx (latS == target.intIdx).

**Effect:** When a type-1 bead has both type-1 and type-2 neighbors, nnp(p, sp, pos) counts ALL of them, while ham! via nnp(m, sp, pos) counts only the target type. This makes deltaHam! overcount cross-species neighbors.

### Bug 2: Missing bidirectional cross-species energy
ham! counts cross-species interactions without dividing by 2:

This means each cross bond between type-1 bead A and type-2 bead B is counted from both sides:

From A: intMat[1,2] × 1
From B: intMat[2,1] × 1
Total per bond: 2 × intMat[1,2] (since the matrix is symmetric)
But deltaHam! for an addBond only computes from the new beads' perspective:

It doesn't account for the existing cross-type neighbor beads also gaining a new neighbor in the ham! sum. So deltaHam! computes 1 × intMat[1,2] per new cross bond, but ham! change is 2 × intMat[1,2].

### Why This Was Never Caught

The existing hamiltonian test for heterogeneous polymers (test_hamiltonians.jl line 208) sets intPot!(m1, -0.3) and intPot!(m2, -0.5) but never sets intPot!(m1, m2, value). So intMat[1,2] = 0.0, and both bugs produce 0 × anything = 0. The test passes by accident.

Even with a default RNG, the test ran 5000 frames without hitting a configuration where the bugs matter. But with a fixed seed (42) or just running longer, the bugs are easily triggered on main.

### Why It Surfaced Now

Phase 3 added boundary rejection logic in move! that changes which moves get accepted/rejected. This shifts the RNG sequence, causing previously-untriggered lattice configurations (where different-type beads become spatial neighbors, not just chain neighbors) to occur within the first few frames.

### The Correct Fix

The fix requires changing how deltaHam! accounts for same-polymer cross-species energy — this is the same code path that _addSwapEnergy! and _removeSwapEnergy! also modify.

Both bugs need to be fixed together (changing nnp dispatch without the factor of 2, or vice versa, produces wrong results).
The swap energy functions (_addSwapEnergy!/_removeSwapEnergy!) also have cross-species terms that need the same treatment.
The fix is:

In deltaHam!: change nnp(p, sp, ...) → nnp(m, sp, ...) for same-polymer PolyMonomers AND multiply the cross-species energy by 2
In _addSwapEnergy!/_removeSwapEnergy!: multiply all crossInteraction terms by 2
This applies to addBond, kinkBond, and removeBond sections
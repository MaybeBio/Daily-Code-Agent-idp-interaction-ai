# Heterogeneous Polymer Implementation Notes

## Overview
This document captures key lessons learned during the implementation of heterogeneous polymer support (multiple monomer types with distinct interaction potentials).

---

## Architecture

### Key Data Structures
- `Polymer.monomers::Dict{Int32, PolyMonomer}` - Maps `intIdx` to monomer type
- `PolyMonomer.intIdx` - Unique species index for this monomer type
- `PolyMonomer.poly` - Back-reference to parent Polymer
- `Polymer.deltaHP` - Accumulated partner energy change (shared by all monomers of same polymer)

### Species Indexing
- Each `PolyMonomer` type gets its own `intIdx` in the global species registry
- `intSpecies` lists contain `(species, intIdx)` tuples for interaction partners
- Bulk species interacting with heterogeneous polymers will have **multiple** PolyMonomer entries in their `intSpecies`

---

## Critical Bugs & Fixes

### 1. `setDeltaHP` Must Accumulate, Not Overwrite
**Problem:** When bulk moves near a heterogeneous polymer, `setDeltaHP` is called for each monomer type (m1, m2). Original implementation used `=` which meant only the last call's value was kept.

**Fix:** 
```julia
# Wrong
setDeltaHP(m::AbstractPolyMonomer, val) = m.poly.deltaHP = val

# Correct
setDeltaHP(m::AbstractPolyMonomer, val) = m.poly.deltaHP += val
```

**Implication:** Must reset `polymer.deltaHP = 0.0` before accumulation loop.

### 2. Same-Polymer PolyMonomer Check in deltaHam!
**Problem:** When polymer moves and iterates over `intSpecies`, it may encounter other `PolyMonomer` types of the **same** polymer. Calling `setDeltaHP` on these incorrectly modifies `p.deltaHP`.

**Fix:** Check if the interaction partner is a PolyMonomer of the same polymer:
```julia
if sp isa PolyMonomer && sp.poly === p
    p.deltaH += energy_change  # Add directly, don't use setDeltaHP
else
    setDeltaHP(sp, energy_change)
    p.deltaH += getDeltaHP(sp)
end
```

### 3. `nnp` Must Count Specific Monomer Type, Not All Polymer Sites
**Problem:** Generic `nnp(bulk, polymonomer, pos)` counted ALL polymer sites (`latS > 0`), not just sites of the target monomer type.

**Fix:** Added specialized method:
```julia
function nnp(s::AbstractBulkMonomer, target::AbstractPolyMonomer, pos)
    # Count only neighbors matching target.intIdx
    n += (pLatS[...] == target.intIdx)  # Not just > 0
end
```

### 4. Double-Counting in Bulk deltaHam!
**Problem:** In bulk's `deltaHam!` loop over `intSpecies`:
```julia
for (sp, idx) in b.intSpecies
    setDeltaHP(sp, energy)      # Accumulates to polymer.deltaHP
    b.deltaH += getDeltaHP(sp)  # Returns accumulated value!
end
```
After m2, `getDeltaHP` returns `energy_m1 + energy_m2`, causing double-count.

**Fix:** Use the calculated energy directly:
```julia
for (sp, idx) in b.intSpecies
    energy_change = ...
    setDeltaHP(sp, energy_change)
    b.deltaH += energy_change  # Not getDeltaHP(sp)
end
```

---

## Testing Strategy

### Time Evolution Consistency Test
The gold standard test: run many moves, track both:
1. `ham!()` - fresh energy calculation
2. Accumulated `deltaH`/`deltaHP` increments

These must match at every step. Discrepancies reveal bugs in incremental energy calculations.

### Debugging Tips
1. **Print the diff**: `hamss[j][end] - accDeltaHs[j][end]` often equals an interaction potential, revealing which interaction is miscounted
2. **Check `@which`**: Verify the correct specialized method is dispatched
3. **Track all species**: When bulk moves, polymer energy changes too (via `deltaHP`)

---

## Invariants to Maintain

1. **Symmetry**: If A-B interaction changes A's energy by ε, it changes B's energy by ε too
2. **Reset before accumulate**: Always reset `deltaHP = 0.0` before loops that call `setDeltaHP`
3. **Method dispatch**: More specific methods (`AbstractBulkMonomer, AbstractPolyMonomer`) must exist for heterogeneous cases

---

## Files Modified
- `src/utils.jl` - `setDeltaHP` accumulation, specialized `nnp` method
- `src/hamiltonian.jl` - Reset logic, same-polymer checks, double-counting fix
- `test/test_hamiltonians.jl` - Time evolution tests for heterogeneous polymers

# Refactor Hamiltonian to Support PolyMonomer Interactions

## Overview
Modify `hamiltonian.jl` to support energy calculations for polymers with bead-dependent interactions via PolyMonomers. The key change is accessing interaction potentials (`intPot`) and chemical potentials (`chemPot`) from the appropriate PolyMonomer for each bead, rather than from a single species-level value.

Once this refactoring is complete, remove the now-redundant `chemPot` and `intPot` fields from the `Polymer` struct, along with any initialization code that sets them.

## Objectives

### 1. Create a Monomer Potential Accessor Helper
Add a new helper function in `utils.jl` to retrieve potentials for any species:

```julia
function getMonomerPotentials(species::AbstractSpecies, beadKey::Int32) 
    → (chemPot::Float64, intPot::Float64)
```

- **For `Polymer`**: Use `getPolyMonomer(poly, beadKey)` to find the PolyMonomer, then return its `chemPot` and `intPot`
- **For `AbstractBulkMonomer` and other `AbstractSpecies`**: Return the species' own `chemPot` and `intPot`

This follows the same dispatch pattern as existing helpers like `getLatS()` and `getTempVec()`.

### 2. Refactor All Hamiltonian Functions
In `hamiltonian.jl`, replace all direct field access (`s.chemPot`, `s.intPot`) with calls to the new accessor:
- `ham!(s::AbstractSpecies)` and related functions
- `deltaH()` and all its variants and overloads
- Any other energy calculation functions that read these potentials

Use the helper like:
```julia
chemPot, intPot = getMonomerPotentials(species, beadKey)
```

### 3. Remove Obsolete Fields from Polymer
Once all hamiltonian functions use the new accessor:
- Remove `chemPot::Union{Nothing,Float64}` field from the `Polymer` struct
- Remove `intPot::Union{Nothing,Float64}` field from the `Polymer` struct

### 4. Clean Up Initialization Code
In `init.jl`:
- Remove all lines that assign `chemPot` or `intPot` to a `Polymer` object
- The `initMonomers!()` function already initializes these fields on individual `PolyMonomer` objects, so no changes needed there

### 5. Circle back on the code and check that everything looks good
- Do a round of feedback to check that everything is as clean and polished as it could be

### 6. Test the changes
- Write a comprehensive test suite for the new features of the hamiltonian functions with PolyMonomers, ensuring backward compatibility

## Implementation Guidelines
- **Minimize changes**: Update only the necessary functions and lines
- **Follow existing patterns**: Mirror the design of `getLatS()`, `getTempVec()`, and `getBoxlength()` in `utils.jl`
- **Avoid duplication**: Don't create separate functions for Polymer vs. other species; use dispatch
- **Complete the refactor**: Ensure all places that read `chemPot`/`intPot` are updated before removing the fields
- **Verify**: Confirm that both Polymer and Bulk species continue to work correctly with energy calculations

## Files to Modify (Tentative)
1. **utils.jl** - Add `getMonomerPotentials()` helper function
2. **hamiltonian.jl** - Update all energy calculation functions to use the new accessor
3. **types.jl** - Remove `chemPot` and `intPot` fields from `Polymer` struct definition
4. **init.jl** - Remove initialization of `chemPot`/`intPot` for Polymer objects

## Success Criteria
✓ All `Polymer` energy calculations correctly read potentials from PolyMonomer objects  
✓ All `BulkMonomer` and `FastBulkMonomer` energy calculations continue to work as before  
✓ No code duplication between Polymer and other species paths  
✓ `Polymer` struct is simpler (fewer fields, clearer intent)  
✓ Tests pass for both Polymer and Bulk species with various interaction configurations

# Multicomponent Polymer Implementation - Phase 1: Initialization

## Background

The codebase has a `PolyMonomer` type that represents constituents of a `Polymer` structure. These monomers:
- Handle interaction potential and chemical potential for polymers
- Enable multicomponent interactions with bead-dependent interaction energies
- Do NOT store their own position (`latS` and `pos`)
- Do NOT have runtime/initialization containers for position (unlike `BulkMonomer` and `FastBulkMonomer`)
- Only behave like other monomer types for interaction matrix and `intSpecies` machinery (see `interaction.jl`)

## Objective

Enable `PolyMonomer` initialization to support the following API:

```julia
l = Lattice(;
    boxlength = 4,
    nspecies = 3
)

tf = FastBulkMonomer(;
    length = 12,
    initShape = "uniform"
)

poly = Polymer(;
    length = l.boxlength,
    maxlength = 8,
    nMonomerTypes = 3,
    chemPot = 0.0,
    intPot = 0.0
)

initMonomers!(poly)

dna = poly.monomers[1]
promoter = poly.monomers[2]

dna.name = "dna"
promoter.name = "promoter"

intPot!(tf, dna, 0.3)
intPot!(tf, promoter, 2.0)
species = [tf, dna, promoter]

init!(l, species)
```

## Technical Requirements

### 1. Position/Lattice Access Pattern

The hamiltonian and related initialization functions must handle two cases:

**Case 1: Bulk Monomers**
- Use the monomer's own `pos` and `latS` fields

**Case 2: PolyMonomers**
- Delegate to parent `Polymer`'s `pos` and `latS` fields

### 2. Implementation Scope

**In Scope for Phase 1:**
- Complete the initialization pathway for `PolyMonomer` species
- Modify hamiltonian functions ONLY as required for `init!` to work
- Ensure proper delegation of position/lattice access
- Comprehensive test suite for initialization and related functions

**Out of Scope (Future Phases):**
- Move functionality modifications
- Full hamiltonian refactoring beyond initialization needs
- Runtime dynamics for multicomponent polymers

## Deliverables

1. **Code Changes:**
   - Review and complete position/lattice delegation mechanism
   - Ensure `initMonomers!(poly)` works correctly
   - Enable `intPot!` calls between different monomer types
   - Make `init!(l, species)` work with mixed species arrays

2. **Test Suite:**
   - Write comprehensive tests for initialization functions
   - Test multicomponent polymer functionalities
   - Verify interaction potential setup between different monomer types
   - Run all tests and ensure they pass

## Implementation Notes

- The foundation has been started but needs completion
- Focus on minimal changes to achieve initialization functionality
- Maintain separation between bulk and poly monomer behavior
- Preserve existing interaction.jl machinery

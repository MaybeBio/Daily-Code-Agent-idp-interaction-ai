# Three-Way Interaction Implementation Plan

## Overview

This document outlines the implementation of three-way interactions for `BulkMonomer` species. A three-way interaction contributes to the Hamiltonian of a BulkMonomer species A (self) based on nearest-neighbor interactions with species B (partner), when A occupies the same lattice site as species C (context).

**Key Constraints:**
- Only `BulkMonomer` can be the "self" species (index `i`)
- Partners and context (`j`, `k`) can be any species (BulkMonomer or PolyMonomer)
- `i`, `j`, `k` must all be distinct species
- Coupling is symmetric in `(j, k)`: `coupling[i,j,k] == coupling[i,k,j]`
- Self-coupling terms are zero (enforced, not stored)
- No reciprocal energy updates to partners
- `FastBulkMonomer` does not support three-way interactions
- `Polymer` Hamiltonian is unchanged (PolyMonomers don't see these couplings)

---

## Section 1: Fix `nnpq` for PolyMonomer Type Matching

### Description

The current `nnpq` function uses `qLatS[...] > 0` to check for context species overlap. For heterogeneous polymers, different bead types have different `intIdx` values stored in the lattice. We need to check for the specific `intIdx` when the context species is a `PolyMonomer`.

### Changes

**File:** `src/utils.jl`

1. Create a specialized `nnpq` method for when `q` is a `PolyMonomer`:
   ```julia
   function nnpq(s::AbstractSpecies, 
                 p::AbstractSpecies, 
                 q::AbstractPolyMonomer,  # specialized for PolyMonomer context
                 pos::Union{SubArray{Int32}, Vector{Int32}})
   ```
   - Change the overlap check from `qLatS[...] > 0` to `qLatS[...] == q.intIdx`

2. Keep the existing generic `nnpq` for non-PolyMonomer context species (uses `> 0` check).

3. Similarly, add specialized version when `p` is a `PolyMonomer` to check `pLatS[...] == p.intIdx` instead of `> 0`.

### Testing Checklist

- [ ] `nnpq` with BulkMonomer context returns same result as before (regression)
- [ ] `nnpq` with PolyMonomer context only counts when specific `intIdx` matches
- [ ] `nnpq` with heterogeneous polymer correctly distinguishes between different PolyMonomer types at neighboring sites
- [ ] `nnpq` with PolyMonomer partner only counts neighbors matching that PolyMonomer's `intIdx`
- [ ] `nnpq` returns 0 when context species is not present at position
- [ ] `nnpq` returns 0 when no neighbors of partner type exist
- [ ] `nnpq` correctly handles periodic boundary conditions

---

## Section 2: Add Coupling Data Structures to `BulkMonomer`

### Description

Add fields to store three-way coupling potentials in the `BulkMonomer` struct.

### Changes

**File:** `src/types.jl`

1. Add new fields to `BulkMonomer` struct:
   ```julia
   # Three-way coupling fields
   couplingMat::Union{Nothing, Array{Float64, 3}} = nothing  # m × n × n matrix
   couplingIdx::Union{Nothing, Int32} = nothing  # this species' index in couplingMat (1st dimension)
   couplingSpecies::Union{Nothing, Vector{Tuple{AbstractMonomer, Int32}}} = nothing  # shared species list for 2nd/3rd dims
   ```

2. The `couplingMat` has dimensions:
   - 1st dimension (`m`): BulkMonomer species only (owners of couplings)
   - 2nd/3rd dimensions (`n`): All interacting species (same as `intSpecies`)

3. Note: `couplingSpecies` will reference the same list as `intSpecies` for the 2nd/3rd dimensions, ensuring consistent indexing.

### Testing Checklist

- [ ] `BulkMonomer` can be constructed with default `nothing` values for coupling fields
- [ ] Coupling fields do not affect existing functionality when `nothing`
- [ ] `FastBulkMonomer` struct is unchanged (no coupling fields)
- [ ] `PolyMonomer` struct is unchanged (no coupling fields)

---

## Section 3: Implement `couplingPot!` Initialization Function

### Description

Create the interface function `couplingPot!(i, j, k, value)` to set three-way coupling potentials, similar to `intPot!`.

### Changes

**File:** `src/interactions.jl`

1. Add helper function `_ensureCouplingStructures!(b::BulkMonomer)`:
   - Initialize `couplingMat` if `nothing`
   - Set `couplingIdx` for this BulkMonomer
   - Link `couplingSpecies` to `intSpecies` (they share the same species list)
   - Ensure `intSpecies` is initialized first (coupling depends on interaction network)

2. Add helper function `_syncCouplingMatrices!(a::BulkMonomer, b::BulkMonomer)`:
   - When two BulkMonomers need to share coupling information
   - Merge and expand coupling matrices similar to `_syncMatrices!`

3. Implement main function:
   ```julia
   function couplingPot!(i::BulkMonomer, j::AbstractMonomer, k::AbstractMonomer, pot::Float64)
   ```
   - **Validation:**
     - Assert `i !== j && j !== k && i !== k` (all distinct species)
     - Assert `i` is a `BulkMonomer` (not `FastBulkMonomer`)
   - **Initialization:**
     - Call `_ensureIntStructures!` on all three species (needed for indexing)
     - Call `_ensureCouplingStructures!` on `i`
     - Sync matrices if species don't share the same `intMat`
   - **Assignment:**
     - Set `couplingMat[i.couplingIdx, j_idx, k_idx] = pot`
     - Set `couplingMat[i.couplingIdx, k_idx, j_idx] = pot` (symmetry in j,k)

4. Add export for `couplingPot!`

### Design Notes

- The 1st dimension of `couplingMat` uses a separate index (`couplingIdx`) that only counts BulkMonomers
- The 2nd/3rd dimensions use the standard `intIdx` from the interaction matrix
- When expanding matrices due to new species, preserve existing coupling values

### Testing Checklist

- [ ] `couplingPot!(i, j, k, val)` correctly sets the coupling value
- [ ] Symmetry enforced: after setting, `couplingMat[i, j, k] == couplingMat[i, k, j]`
- [ ] Error thrown when `i === j`
- [ ] Error thrown when `j === k`
- [ ] Error thrown when `i === k`
- [ ] Error thrown when `i` is not a `BulkMonomer` (e.g., `FastBulkMonomer`)
- [ ] `intSpecies` is auto-initialized if not already present
- [ ] Multiple calls with different `(j, k)` pairs work correctly
- [ ] Multiple BulkMonomers can have independent couplings
- [ ] Adding new species after initial coupling setup expands matrices correctly
- [ ] Coupling values survive matrix expansion/sync operations
- [ ] `couplingPot!` with PolyMonomer as `j` or `k` works correctly
- [ ] Zero coupling (explicit) can be set and retrieved

---

## Section 4: Update `ham!` for `BulkMonomer`

### Description

Add three-way coupling contributions to the full Hamiltonian calculation for `BulkMonomer`.

### Changes

**File:** `src/hamiltonian.jl`

1. Modify `ham!(b::AbstractBulkMonomer)` function:
   - After the existing interaction loop, add coupling contribution:
   ```julia
   # Three-way coupling contribution
   if !isnothing(b.couplingMat)
       for (sp_j, idx_j) in b.intSpecies
           for (sp_k, idx_k) in b.intSpecies
               # Skip if same species (no self-coupling)
               if idx_j >= idx_k
                   continue  # Avoid double-counting due to symmetry
               end
               coupling_val = b.couplingMat[b.couplingIdx, idx_j, idx_k]
               if coupling_val != 0.0
                   for i in 1:b.length
                       pos = b.pos[i][1:3]
                       b.ham += coupling_val * nnpq(b, sp_j, sp_k, pos)
                   end
               end
           end
       end
   end
   ```

2. Important: The loop structure ensures each `(j, k)` pair is counted once (not twice due to symmetry).

### Testing Checklist

- [ ] `ham!` returns same result when `couplingMat` is `nothing` (regression)
- [ ] `ham!` correctly adds coupling contribution when present
- [ ] Coupling contribution is zero when no species overlap at any position
- [ ] Coupling contribution correctly accounts for overlapping species
- [ ] No double-counting of symmetric `(j, k)` pairs
- [ ] `ham!` with multiple coupling terms sums correctly
- [ ] `ham!` handles PolyMonomer as partner/context correctly
- [ ] Performance: `ham!` completes in reasonable time with couplings

---

## Section 5: Update `deltaHam!` for `BulkMonomer`

### Description

Add three-way coupling contributions to the incremental Hamiltonian calculation when a BulkMonomer bead moves.

### Changes

**File:** `src/hamiltonian.jl`

1. Modify `deltaHam!(b::BulkMonomer)` function:
   - After existing interaction calculations, add:
   ```julia
   # Three-way coupling contribution
   if !isnothing(b.couplingMat)
       for (sp_j, idx_j) in b.intSpecies
           for (sp_k, idx_k) in b.intSpecies
               # Skip if same species or if j > k (symmetry)
               if idx_j >= idx_k
                   continue
               end
               coupling_val = b.couplingMat[b.couplingIdx, idx_j, idx_k]
               if coupling_val != 0.0
                   delta_coupling = coupling_val * (
                       nnpq(b, sp_j, sp_k, b.newPos) - nnpq(b, sp_j, sp_k, b.bi)
                   )
                   b.deltaH += delta_coupling
               end
           end
       end
   end
   ```

2. Note: No `deltaHP` updates for partners (coupling is owned by self only).

### Testing Checklist

- [ ] `deltaHam!` returns same result when `couplingMat` is `nothing` (regression)
- [ ] `deltaHam!` correctly computes coupling delta for move into overlapping site
- [ ] `deltaHam!` correctly computes coupling delta for move out of overlapping site
- [ ] `deltaHam!` returns zero coupling delta when moving between non-overlapping sites
- [ ] `deltaHam!` handles case where old position had overlap, new doesn't
- [ ] `deltaHam!` handles case where new position has overlap, old didn't
- [ ] `deltaHam!` handles case where both positions have different overlaps
- [ ] Consistency: `ham!` before + `deltaHam!` == `ham!` after (for any move)
- [ ] No reciprocal updates to partner `deltaHP` from couplings
- [ ] Multiple coupling terms sum correctly in `deltaH`

---

## Section 6: Ensure `FastBulkMonomer` Unchanged

### Description

Verify that `FastBulkMonomer` continues to work without three-way interactions and explicitly document the exclusion.

### Changes

**File:** `src/hamiltonian.jl`

1. Add a comment in `deltaHam!(b::FastBulkMonomer)` clarifying that three-way couplings are not supported:
   ```julia
   # Note: FastBulkMonomer does not support three-way couplings.
   # Use BulkMonomer if three-way interactions are needed.
   ```

2. No functional changes to `FastBulkMonomer`.

**File:** `src/interactions.jl`

1. In `couplingPot!`, add explicit type check:
   ```julia
   @assert !(i isa FastBulkMonomer) "FastBulkMonomer does not support three-way couplings"
   ```

### Testing Checklist

- [ ] `FastBulkMonomer` `deltaHam!` unchanged (regression)
- [ ] `FastBulkMonomer` `ham!` unchanged (regression, if applicable)
- [ ] `couplingPot!` with `FastBulkMonomer` as first argument throws error
- [ ] Existing spec files using `FastBulkMonomer` still work

---

## Section 7: Integration Tests

### Description

Create comprehensive integration tests that verify the entire three-way interaction system works correctly in realistic simulation scenarios.

### Test Scenarios

1. **Basic Three-Species Test:**
   - Create three BulkMonomer species A, B, C
   - Set two-way interactions: A-B, A-C, B-C
   - Set three-way coupling: `couplingPot!(A, B, C, -1.0)`
   - Manually place particles to create known overlaps
   - Verify `ham!` matches hand calculation
   - Perform move, verify `deltaHam!` matches

2. **BulkMonomer + PolyMonomer Test:**
   - Create BulkMonomer A, B and PolyMonomer P
   - Set three-way coupling: `couplingPot!(A, B, P, -0.5)`
   - Place polymer with specific pattern
   - Place BulkMonomers overlapping polymer
   - Verify coupling only activates at polymer sites

3. **Heterogeneous Polymer Test:**
   - Create BulkMonomer A, B
   - Create Polymer with two PolyMonomer types P1, P2
   - Set `couplingPot!(A, B, P1, -1.0)` (only P1, not P2)
   - Verify coupling only activates at P1 sites, not P2

4. **Monte Carlo Consistency Test:**
   - Run short simulation with couplings
   - Periodically recalculate full `ham!`
   - Verify accumulated `deltaHam!` matches

5. **Detailed Balance Test:**
   - Verify that forward and reverse moves have correct energy changes
   - Important for ensuring correct MC acceptance probabilities

### Testing Checklist

- [ ] Three BulkMonomer species with three-way coupling: energy correct
- [ ] BulkMonomer + PolyMonomer coupling: energy correct at polymer sites
- [ ] Heterogeneous polymer: coupling distinguishes PolyMonomer types
- [ ] MC consistency: accumulated deltas match recalculated Hamiltonian
- [ ] Detailed balance: forward/reverse move energies are symmetric
- [ ] No memory leaks or allocations in hot path
- [ ] Performance acceptable for typical system sizes

---

## Section 8: Documentation and Examples

### Description

Add documentation and example usage for the new three-way interaction feature.

### Changes

1. **File:** `context/domain_knowledge.md` (or appropriate doc file)
   - Add section explaining three-way interactions
   - Document physical interpretation
   - Document API usage

2. **Example in spec file or test:**
   ```julia
   # Example: Cooperative binding where TF-DNA interaction is enhanced 
   # when another TF is already bound nearby
   
   tf1 = BulkMonomer(name="TF1", length=100)
   tf2 = BulkMonomer(name="TF2", length=50)
   dna = poly.monomers[1]  # PolyMonomer
   
   # Two-way interactions
   intPot!(tf1, -1.0)
   intPot!(tf2, -1.0)
   intPot!(tf1, dna, -0.5)
   intPot!(tf2, dna, -0.5)
   
   # Three-way coupling: tf1-dna interaction enhanced when tf2 is also at dna site
   couplingPot!(tf1, tf2, dna, -0.3)  # Extra stabilization
   ```

### Testing Checklist

- [ ] Documentation is clear and accurate
- [ ] Example code runs without errors
- [ ] Example produces expected physical behavior

---

## Implementation Order

1. **Section 1** - Fix `nnpq` (prerequisite for correct coupling calculation)
2. **Section 2** - Add data structures (prerequisite for initialization)
3. **Section 3** - Implement `couplingPot!` (prerequisite for energy functions)
4. **Section 4** - Update `ham!` (needed for validation)
5. **Section 5** - Update `deltaHam!` (core functionality)
6. **Section 6** - Ensure `FastBulkMonomer` excluded
7. **Section 7** - Integration tests
8. **Section 8** - Documentation

Each section should be completed and tested before moving to the next.

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|------------|--------|------------|
| Performance degradation in `deltaHam!` | Medium | High | Early profiling, skip zero couplings |
| Matrix indexing errors | Medium | High | Careful testing, clear documentation of index conventions |
| Inconsistency between `ham!` and `deltaHam!` | Low | Critical | Consistency tests in every section |
| Breaking existing functionality | Low | High | Regression tests before each change |
| Memory overhead from 3D matrix | Low | Low | Matrix is typically small (< 1KB) |

---

## Success Criteria

1. All existing tests pass (no regressions)
2. New tests for three-way interactions pass
3. `ham!` and `deltaHam!` are consistent for all test cases
4. `couplingPot!` interface is intuitive and validates inputs
5. Performance impact is negligible (< 5% slowdown for typical cases)
6. Documentation is complete and accurate

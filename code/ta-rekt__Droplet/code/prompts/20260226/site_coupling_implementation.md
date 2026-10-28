# Site-Coupling Implementation Plan

**Date:** 2026-02-26  
**Goal:** Replace the neighbor-based coupling (`nnpq`) with a purely on-site 
co-location coupling (`cplpq`). Also add two-way on-site coupling via 
`couplingMat[i,j,j]`, and enforce full symmetry across the first two indices 
of `couplingMat` so that A,B,C at the same site contribute the same energy 
to each species.

---

## Summary of Changes

| Old (nnpq) | New (cplpq) |
|---|---|
| "partner = nearest neighbor, context = same site" | "all participants co-located at same site" |
| Returns integer 0–6 (neighbor count) | Returns Boolean 0 or 1 |
| Only three-way (`i ≠ j ≠ k`) | Three-way (`i ≠ j ≠ k, j < k`) + two-way (`j == k, j ≠ i`) |
| Asymmetric in j,k roles | Symmetric in all indices |
| `couplingPot!(self, partner, context, pot)` | Three-way: `couplingPot!(self, j, k, pot)` (j ≠ k ≠ self) + Two-way: `couplingPot!(self, partner, pot)` |
| Owner-only: `couplingMat[i,j,k]` independent per owner | Symmetric: `couplingMat[i,j,k]` enforced equal for all owners sharing the triplet |

### Key Design Decisions

1. **`cplp` is removed.** All coupling uses `cplpq(j, k, pos)`, including 
   the two-way case where `j == k` → `cplpq(j, j, pos)` which degenerates 
   to a single occupancy check, but the uniform interface avoids branching 
   in `ham!`/`deltaHam!`.

2. **`nnpq` is fully replaced** by `cplpq`. The `nnpq` functions (all 4 
   dispatch variants) and `test_nnpq.jl` are deleted.

3. **`couplingMat[i,j,j]`** is the two-way on-site coupling between species 
   `i` and `j`. It means "energy contributed when i and j co-occupy a site."  
   Constraint: `i ≠ j` (no self-coupling). This is distinct from 
   `intMat[i,j]` which is the nearest-neighbor contact energy.

4. **Full first-two-index symmetry in `couplingMat`:** When 
   `couplingPot!(a, b, c, pot)` is called, the implementation sets:
   - `couplingMat[a_idx, b_intIdx, c_intIdx] = pot`
   - `couplingMat[a_idx, c_intIdx, b_intIdx] = pot`  *(j,k symmetry)*
   - `couplingMat[b_idx, a_intIdx, c_intIdx] = pot`  *(first-two symmetry)*
   - `couplingMat[b_idx, c_intIdx, a_intIdx] = pot`  *(j,k symmetry)*
   - `couplingMat[c_idx, a_intIdx, b_intIdx] = pot`  *(first-two symmetry)*
   - `couplingMat[c_idx, b_intIdx, a_intIdx] = pot`  *(j,k symmetry)*

   This means **all three species** involved in a coupling get 
   `_ensureCouplingStructures!` called, and all three share the coupling 
   matrix (via sync). The same coupling value appears from any species' 
   perspective.

   For the two-way case `couplingPot!(a, b, pot)`:
   - `couplingMat[a_idx, b_intIdx, b_intIdx] = pot`
   - `couplingMat[b_idx, a_intIdx, a_intIdx] = pot`  *(first-two symmetry)*

5. **Loop structure in `ham!`:** Two passes over `intSpecies` pairs:
   - **Pass 1 (two-way, `j == k`):** Loop over `(sp_j, idx_j)`, read 
     `couplingMat[couplingIdx, idx_j, idx_j]`. If nonzero, sum 
     `cplpq(sp_j, sp_j, pos)` over all beads of self.
   - **Pass 2 (three-way, `j < k`):** Loop over `(sp_j, idx_j)` and 
     `(sp_k, idx_k)` with `idx_j < idx_k`. Read 
     `couplingMat[couplingIdx, idx_j, idx_k]`. If nonzero, sum 
     `cplpq(sp_j, sp_k, pos)` over all beads of self.

6. **No `selfFactor`** needed in the coupling loop because:
   - Two-way: `j == k` is a single term per pair, not double-counted.
   - Three-way: `j < k` ensures each pair is counted once, and 
     `cplpq(j, k, pos)` is inherently symmetric.

---

## Section 1: `cplpq` — Dispatch Variants & Delete `nnpq`/`cplp`

### 1.1 Add PolyMonomer dispatches for `cplpq`

The generic `cplpq(p, q, pos)` checks `pLatS[pos] > 0` and `qLatS[pos] > 0`.  
For PolyMonomer types, the lattice stores a monomer-type index (`intIdx`), so 
`> 0` is insufficient — it must check `== p.intIdx` to discriminate monomer 
types in a heterogeneous polymer.

**Four dispatch variants (mirrors the `nnpq` pattern):**

```
cplpq(p::AbstractSpecies,      q::AbstractSpecies,      pos)  # already exists
cplpq(p::AbstractPolyMonomer,  q::AbstractSpecies,      pos)  # check p.intIdx
cplpq(p::AbstractSpecies,      q::AbstractPolyMonomer,  pos)  # check q.intIdx
cplpq(p::AbstractPolyMonomer,  q::AbstractPolyMonomer,  pos)  # check both
```

**Implementation detail:** Each variant follows the same pattern:
```julia
function cplpq(p::AbstractPolyMonomer, q::AbstractSpecies, pos)
    pLatS = getLatS(p)
    qLatS = getLatS(q)
    return (pLatS[pos[1], pos[2], pos[3]] == p.intIdx) * (qLatS[pos[1], pos[2], pos[3]] > 0)
end
```

### 1.2 Delete `cplp`

Remove the `cplp` function from `src/utils.jl` (lines ~325–332). It is 
unused in production code and superseded by `cplpq(j, j, pos)`.

### 1.3 Delete `nnpq`

Remove all four `nnpq` dispatch variants from `src/utils.jl` (lines ~356–470+).

### 1.4 Delete `test_nnpq.jl`

Remove `test/test_nnpq.jl` from disk and remove the `include("test_nnpq.jl")` 
line from `test/runtests.jl`.

### Tests for Section 1

**New file: `test/test_cplpq.jl`** (replaces `test_nnpq.jl`)

| Test | Description |
|---|---|
| `cplpq - generic: both present` | Place two BulkMonomers at same site → returns 1 |
| `cplpq - generic: one absent` | Place only one at site → returns 0 |
| `cplpq - generic: both absent` | Neither at site → returns 0 |
| `cplpq - PolyMonomer p: correct intIdx` | p is PolyMonomer at site with matching intIdx → 1 |
| `cplpq - PolyMonomer p: wrong intIdx` | p is heterogeneous, site has other monomer type → 0 |
| `cplpq - PolyMonomer q: correct intIdx` | Symmetric of above |
| `cplpq - PolyMonomer q: wrong intIdx` | Symmetric of above |
| `cplpq - both PolyMonomer` | Both p and q are PolyMonomers from different polymers |
| `cplpq - both PolyMonomer, heterogeneous` | Two PolyMonomers from same heterogeneous polymer |
| `cplpq - self-check (j == j)` | `cplpq(j, j, pos)` when j occupies pos → 1, else 0 |
| `cplpq - periodic boundary` | Site at box edge, verify boundary wrapping (uses `boxmod!` or not? — note: `cplpq` uses raw position indices, no wrapping needed since `pos` is already a valid lattice coordinate) |

**Add `include("test_cplpq.jl")` to `runtests.jl`** in place of `test_nnpq.jl`.

---

## Section 2: `couplingPot!` — API Changes

### 2.1 New two-argument signature

```julia
couplingPot!(self::BulkMonomer, partner::AbstractMonomer, pot::Float64)
```

**Validation:**
- `self !== partner` (no self-coupling)
- `self` must be `BulkMonomer`, not `FastBulkMonomer`

**Behavior (with first-two-index symmetry):**
1. `_ensureIntStructures!` for both self and partner
2. Sync intMat between self and partner
3. `_ensureCouplingStructures!` for **both** self and partner (since both 
   will own entries in `couplingMat`)
4. `_syncCouplingStructures!` between self and partner (so both share 
   the same `couplingMat`)
5. `_expandCouplingMat!` if needed
6. Set:
   - `couplingMat[self.couplingIdx, partner.intIdx, partner.intIdx] = pot`
   - `couplingMat[partner.couplingIdx, self.intIdx, self.intIdx] = pot`

### 2.2 Modify three-argument signature

```julia
couplingPot!(self::BulkMonomer, partner::AbstractMonomer, context::AbstractMonomer, pot::Float64)
```

**Validation changes:**
- Remove `partner === context` check (was: error if j == k; now: 
  `j == k` is valid and delegates to the two-argument form internally, 
  OR we keep both paths separate and just allow `j == k` to set the 
  same entries as the two-argument form)
- Keep `self === partner` and `self === context` checks

**DECISION:** If `partner === context`, internally delegate to the 
two-argument behavior. This avoids duplicating logic.

**Behavior (with first-two-index symmetry):**
1. `_ensureIntStructures!` for self, partner, context
2. Sync intMat across all three
3. `_ensureCouplingStructures!` for **all three** species
4. `_syncCouplingStructures!` pairwise until all three share the same 
   `couplingMat` and `couplingSpecies`
5. `_expandCouplingMat!` if needed
6. Set all 6 entries (see Design Decision 4 above)

### 2.3 Modify `_ensureCouplingStructures!`

Currently only called for the "self" BulkMonomer. Now it must be callable 
for any `BulkMonomer` that participates in a coupling, including partners 
and contexts that happen to be BulkMonomers.

**Key change:** If partner or context is a BulkMonomer, it also gets a 
`couplingIdx` and appears in `couplingSpecies`. If partner/context is a 
PolyMonomer (or other non-BulkMonomer), it does NOT get coupling structures 
— it only has `intIdx` in the shared `intMat`, and the coupling entry is 
written from the BulkMonomer owners' perspective.

**Wait — symmetry implication:** If species C is a PolyMonomer, we cannot 
write `couplingMat[c_idx, ...]` because PolyMonomer doesn't have 
`couplingMat`. In that case, the first-two-index symmetry only applies to 
the BulkMonomer participants. Concretely:

- `couplingPot!(a::BulkMonomer, b::BulkMonomer, pm::PolyMonomer, pot)`:
  - Set `couplingMat[a_idx, b_intIdx, pm_intIdx] = pot` and j,k-symmetric
  - Set `couplingMat[b_idx, a_intIdx, pm_intIdx] = pot` and j,k-symmetric
  - PolyMonomer `pm` cannot own coupling entries → skip

- `couplingPot!(a::BulkMonomer, pm::PolyMonomer, c::BulkMonomer, pot)`:
  - Set `couplingMat[a_idx, pm_intIdx, c_intIdx] = pot` and j,k-symmetric
  - Set `couplingMat[c_idx, pm_intIdx, a_intIdx] = pot` and j,k-symmetric
  - PolyMonomer `pm` cannot own coupling entries → skip

- `couplingPot!(a::BulkMonomer, b::BulkMonomer, c::BulkMonomer, pot)`:
  - Full 6-entry symmetry as described above

This is consistent: each BulkMonomer that participates sees the same coupling 
value from its own `couplingIdx`.

### 2.4 Modify `_syncCouplingStructures!`

Currently syncs between two BulkMonomers. With the new symmetry requirement, 
we need to sync across all BulkMonomer participants. The existing pairwise 
approach works: call `_syncCouplingStructures!(a, b)` then 
`_syncCouplingStructures!(a, c)` (the second call will see that a and b 
already share structures and just adds c).

### 2.5 FastBulkMonomer guard

Keep the existing `couplingPot!(::FastBulkMonomer, ..., pot)` error.  
Add: `couplingPot!(::FastBulkMonomer, ::AbstractMonomer, pot)` error for 
the two-argument form.

### 2.6 `getCouplingPot` updates

Add a two-argument form:
```julia
getCouplingPot(self::BulkMonomer, partner::AbstractMonomer)
```
Returns `couplingMat[self.couplingIdx, partner.intIdx, partner.intIdx]`.

### Tests for Section 2

**Rewrite `test/test_couplingPot.jl`** — the structure changes significantly.

| Part | Test | Description |
|---|---|---|
| **1: Basic** | `couplingPot! 3-arg - basic set/get` | Set (a,b,c,-0.5), verify getCouplingPot |
| | `couplingPot! 2-arg - basic set/get` | Set (a,b,-0.5), verify getCouplingPot two-arg |
| | `couplingPot! 3-arg - j,k symmetry` | `getCouplingPot(a,b,c) == getCouplingPot(a,c,b)` |
| **2: First-two-index symmetry** | `3-arg symmetry - all BulkMonomer` | `couplingPot!(a,b,c,pot)` → verify `getCouplingPot(a,b,c) == getCouplingPot(b,a,c) == getCouplingPot(c,a,b)` |
| | `3-arg symmetry - one PolyMonomer` | `couplingPot!(a,b,pm,pot)` → verify `getCouplingPot(a,b,pm) == getCouplingPot(b,a,pm)`, pm has no coupling |
| | `2-arg symmetry` | `couplingPot!(a,b,pot)` → verify `getCouplingPot(a,b) == getCouplingPot(b,a)` (where getCouplingPot(b,a) reads couplingMat[b_idx, a_intIdx, a_intIdx]) |
| **3: Validation** | `error: self === partner (3-arg)` | `couplingPot!(a,a,b,pot)` errors |
| | `error: self === context (3-arg)` | `couplingPot!(a,b,a,pot)` errors |
| | `j === k delegates to 2-arg` | `couplingPot!(a,b,b,pot)` works and equivalent to `couplingPot!(a,b,pot)` |
| | `error: self === partner (2-arg)` | `couplingPot!(a,a,pot)` errors |
| | `error: FastBulkMonomer (3-arg)` | existing test |
| | `error: FastBulkMonomer (2-arg)` | `couplingPot!(f,b,pot)` errors |
| **4: Auto-init** | `auto-init intSpecies` | Same as existing |
| | `auto-init coupling structures for all BulkMonomer participants` | Verify that partner/context BulkMonomers also get `couplingIdx` |
| **5: Multiple couplings** | `multiple (j,k) pairs` | Same as existing |
| | `multiple BulkMonomers - shared coupling values` | Verify that after `couplingPot!(a,b,c,pot)`, all three BulkMonomers see same pot |
| | `matrix expansion on new species` | Same as existing |
| | `coupling values survive sync` | Same as existing |
| **6: PolyMonomer** | `PolyMonomer as partner (3-arg)` | Same as existing |
| | `PolyMonomer as context (3-arg)` | Same as existing |
| | `PolyMonomer as partner (2-arg)` | `couplingPot!(a, pm, pot)` sets on-site coupling |
| **7: Edge cases** | `zero coupling` | Same |
| | `overwrite coupling` | Same |
| | `getCouplingPot returns nothing` | Same |
| **8: Matrix dimensions** | `correct dimensions after 3-arg` | Adapted: now more BulkMonomers in couplingSpecies |
| | `correct dimensions after 2-arg` | New |

---

## Section 3: `ham!` for `AbstractBulkMonomer` — Coupling Loop

### 3.1 Replace `nnpq` with `cplpq` and restructure loop

Current coupling block in `ham!(b::AbstractBulkMonomer)`:
```julia
if hasfield(typeof(b), :couplingMat) && !isnothing(b.couplingMat)
    for (sp_j, idx_j) in b.intSpecies
        for (sp_k, idx_k) in b.intSpecies
            if idx_j >= idx_k
                continue
            end
            coupling_val = b.couplingMat[b.couplingIdx, idx_j, idx_k]
            if coupling_val != 0.0
                for i in 1:b.length
                    pos = b.pos[i][range(1, 3)]
                    b.ham += coupling_val * (nnpq(b, sp_j, sp_k, pos) + nnpq(b, sp_k, sp_j, pos))
                end
            end
        end
    end
end
```

**New coupling block:**
```julia
if hasfield(typeof(b), :couplingMat) && !isnothing(b.couplingMat)
    # Pass 1: two-way on-site coupling (j == k, j ≠ self)
    for (sp_j, idx_j) in b.intSpecies
        if idx_j == b.intIdx
            continue  # skip self-coupling (couplingMat[i,i,i] disallowed)
        end
        coupling_val = b.couplingMat[b.couplingIdx, idx_j, idx_j]
        if coupling_val != 0.0
            for i in 1:b.length
                pos = b.pos[i][range(1, 3)]
                b.ham += coupling_val * cplpq(sp_j, sp_j, pos)
            end
        end
    end

    # Pass 2: three-way on-site coupling (j < k)
    for (sp_j, idx_j) in b.intSpecies
        for (sp_k, idx_k) in b.intSpecies
            if idx_j >= idx_k
                continue
            end
            coupling_val = b.couplingMat[b.couplingIdx, idx_j, idx_k]
            if coupling_val != 0.0
                for i in 1:b.length
                    pos = b.pos[i][range(1, 3)]
                    b.ham += coupling_val * cplpq(sp_j, sp_k, pos)
                end
            end
        end
    end
end
```

**Key differences from old code:**
- `nnpq(b, sp_j, sp_k, pos) + nnpq(b, sp_k, sp_j, pos)` → `cplpq(sp_j, sp_k, pos)`
  - `cplpq` is inherently symmetric, so no need to sum both orderings
  - No reference to `b` (self) in the call — `cplpq` only checks the two 
    partner species' lattices
- New Pass 1 for two-way terms
- Skip `idx_j == b.intIdx` in Pass 1 to avoid self-coupling

### Tests for Section 3

**Rewrite `test/test_ham_coupling.jl`**

| Part | Test | Description |
|---|---|---|
| **1: Regression** | `ham! unchanged when couplingMat is nothing` | Same as existing |
| | `ham! with zero coupling = no effect` | Same as existing |
| **2: Three-way basic** | `coupling active when all three co-located` | Place A,B,C at same site → ham includes coupling |
| | `coupling zero when only two co-located` | A and B at same site, C far away → coupling 0 |
| | `coupling zero when species separated` | All three at different sites → coupling 0 |
| **3: Two-way basic** | `two-way coupling active` | `couplingPot!(a,b,-1.0)`, place A,B at same site → ham includes -1.0 |
| | `two-way coupling zero when apart` | A,B at different sites → 0 |
| | `two-way distinct from intMat` | Set both `intPot!(a,b,-0.3)` and `couplingPot!(a,b,-0.5)`, verify both contribute independently (intMat via neighbors, coupling via co-location) |
| **4: Manual calculation** | `matches hand calculation (3-way)` | Multi-bead species, some beads co-located, some not. Compute expected by hand. |
| | `matches hand calculation (2-way)` | Similar for two-way. |
| | `matches hand calculation (mixed 2+3-way)` | Both two-way and three-way couplings active simultaneously. |
| **5: Symmetry** | `first-two-index symmetry in energy` | After `couplingPot!(a,b,c,pot)`, verify `ham!(a)`, `ham!(b)`, `ham!(c)` each see the coupling term (from the same configuration). Specifically: each owner's coupling contribution should be `pot * cplpq(j,k,pos)` for the beads of that owner. |
| | `no double counting` | With `couplingPot!(a,b,c,-1.0)`, verify coupling contributes `-1.0` once per co-located site per owner (not 2× or 6×). |
| **6: Multiple couplings** | `multiple coupling terms sum correctly` | Same structure as existing, adapted values. |
| **7: PolyMonomer** | `coupling with PolyMonomer` | Same structure as existing, adapted for co-location instead of neighbor+overlap. |
| | `heterogeneous polymer discriminates types` | `cplpq` with PolyMonomer uses `intIdx` — coupling active at pm1 sites but not pm2 sites. |
| **8: Performance** | `performance with couplings` | Same as existing. |

---

## Section 4: `deltaHam!` for `BulkMonomer` — Coupling Delta

### 4.1 Replace `nnpq` with `cplpq` and add two-way pass

Current coupling block in `deltaHam!(b::BulkMonomer)`:
```julia
if !isnothing(b.couplingMat)
    for (sp_j, idx_j) in b.intSpecies
        for (sp_k, idx_k) in b.intSpecies
            if idx_j >= idx_k
                continue
            end
            coupling_val = b.couplingMat[b.couplingIdx, idx_j, idx_k]
            if coupling_val != 0.0
                delta_coupling = coupling_val * (
                    (nnpq(b, sp_j, sp_k, b.newPos) + nnpq(b, sp_k, sp_j, b.newPos)) -
                    (nnpq(b, sp_j, sp_k, b.bi) + nnpq(b, sp_k, sp_j, b.bi))
                )
                b.deltaH += delta_coupling
            end
        end
    end
end
```

**New coupling block:**
```julia
if !isnothing(b.couplingMat)
    # Pass 1: two-way on-site coupling delta (j == k, j ≠ self)
    for (sp_j, idx_j) in b.intSpecies
        if idx_j == b.intIdx
            continue
        end
        coupling_val = b.couplingMat[b.couplingIdx, idx_j, idx_j]
        if coupling_val != 0.0
            b.deltaH += coupling_val * (cplpq(sp_j, sp_j, b.newPos) - cplpq(sp_j, sp_j, b.bi))
        end
    end

    # Pass 2: three-way on-site coupling delta (j < k)
    for (sp_j, idx_j) in b.intSpecies
        for (sp_k, idx_k) in b.intSpecies
            if idx_j >= idx_k
                continue
            end
            coupling_val = b.couplingMat[b.couplingIdx, idx_j, idx_k]
            if coupling_val != 0.0
                b.deltaH += coupling_val * (cplpq(sp_j, sp_k, b.newPos) - cplpq(sp_j, sp_k, b.bi))
            end
        end
    end
end
```

**Key differences from old code:**
- `nnpq(..., newPos) + nnpq(..., newPos)` → `cplpq(sp_j, sp_k, newPos)` (single call, symmetric)
- Same for `bi` terms
- New Pass 1 for two-way delta

### 4.2 `deltaHam!(::FastBulkMonomer)` — no changes

Already has a comment that coupling is not supported. No code changes needed.

### Tests for Section 4

**Rewrite `test/test_deltaHam_coupling.jl`**

| Part | Test | Description |
|---|---|---|
| **1: Regression** | `deltaHam! unchanged when couplingMat is nothing` | Same |
| **2: Consistency** | `ham_before + delta == ham_after (3-way)` | 100 random moves, verify consistency |
| | `ham_before + delta == ham_after (2-way)` | Same but with two-way coupling |
| | `ham_before + delta == ham_after (mixed)` | Both two-way and three-way active |
| | `ham_before + delta == ham_after (multiple terms)` | 4 species, multiple couplings |
| **3: Controlled** | `move into co-located site (3-way)` | A moves to site where B and C are → delta = -coupling_val |
| | `move out of co-located site (3-way)` | A moves away from B,C site → delta = +coupling_val |
| | `move between non-coupled sites` | No coupling active at either end → delta = 0 |
| | `move into co-located site (2-way)` | A moves to site where B is → delta = -coupling_val |
| | `move out of co-located site (2-way)` | A moves away from B site → delta = +coupling_val |
| **4: No reciprocal deltaHP** | `partner deltaHP unchanged by coupling` | Same as existing |
| **5: Zero coupling** | `zero coupling has no effect` | Same structure |
| **6: Symmetry** | `delta consistent across owners` | After `couplingPot!(a,b,c,pot)`, move a to co-located site and move b to co-located site in separate configs — both should see same magnitude delta |
| **7: PolyMonomer** | `consistency with PolyMonomer` | Same structure, adapted |
| **8: Simulation (self-moves)** | `accumulated deltaH matches ham!` | Same structure |
| **9: Full simulation** | `all species move, single-step` | Same structure |
| **10: Performance** | `performance` | Same |

---

## Section 5: `test_integration_coupling.jl` — Full Integration

**Rewrite `test/test_integration_coupling.jl`** to reflect the new semantics.

### Structure

| Test | Description |
|---|---|
| **FastBulkMonomer exclusion** | `couplingPot!` error, ham!/deltaHam! regression, no couplingMat field — same as existing Section 6 tests |
| **Three BulkMonomer - 3-way** | Place A,B,C at same site, verify ham! = coupling_val. Move A away, verify deltaH = +coupling_val. Adapted from existing Test 1. |
| **Three BulkMonomer - 2-way** | `couplingPot!(a,b,-1.0)`, place A,B at same site, verify ham!/deltaH. New. |
| **BulkMonomer + PolyMonomer** | Coupling with PolyMonomer context, consistency check. Adapted. |
| **Heterogeneous polymer** | Coupling discriminates monomer types via `cplpq` PolyMonomer dispatch. Adapted — no neighbor concept, so use co-location only. |
| **MC consistency (self-moves)** | Accumulated getDeltaHam matches ham! over 500 moves. Same pattern. |
| **MC consistency (all species, single-step)** | ham_before + deltaH == ham_after for each species in round-robin. Same pattern. |
| **Detailed balance** | Forward/reverse energy symmetry. Adapted. |
| **Multiple coupling owners** | All BulkMonomers see same coupling value (first-two-index symmetry). Both A and B consistency checks. Adapted. |
| **Mixed FastBulkMonomer + BulkMonomer** | FastBulkMonomer unaffected by coupling. Same pattern. |
| **BulkMonomer + Polymer full simulation** | With polymer present. Same pattern. |
| **Mixed 2-way + 3-way simulation** | Both types of coupling active simultaneously. New. |
| **Performance** | 1000 move! calls timing. Same. |

---

## Section 6: Cleanup & Documentation

### 6.1 Delete deprecated code

- Remove `nnpq` functions from `src/utils.jl` (all 4 dispatches)
- Remove `cplp` function from `src/utils.jl`
- Remove `test/test_nnpq.jl`

### 6.2 Update `runtests.jl`

- Replace `include("test_nnpq.jl")` with `include("test_cplpq.jl")`
- Other includes remain the same

### 6.3 Update exports

- Check if `nnpq` is exported; if so, remove and export `cplpq` instead
- Export `getCouplingPot` two-argument form if needed

### 6.4 Update documentation

- Update `context/domain_knowledge.md` Section 2.4 to reflect co-location 
  semantics instead of neighbor+overlap
- Update `couplingPot!` docstring
- Update `ham!` inline comments

### Tests for Section 6

Run full test suite (`runtests.jl`) to verify no regressions. All existing 
non-coupling tests should pass unchanged.

---

## Execution Order

1. **Section 1** — `cplpq` dispatches + delete `nnpq`/`cplp` + write `test_cplpq.jl`
2. **Section 2** — `couplingPot!` API changes + rewrite `test_couplingPot.jl`
3. **Section 3** — `ham!` coupling loop + rewrite `test_ham_coupling.jl`
4. **Section 4** — `deltaHam!` coupling delta + rewrite `test_deltaHam_coupling.jl`
5. **Section 5** — Rewrite `test_integration_coupling.jl`
6. **Section 6** — Cleanup, docs, full test suite

Each section should be committed after its tests pass. The full test suite 
(including all non-coupling tests) must pass at the end of each section 
to ensure no regressions.

---

## Risk Areas

1. **`_syncCouplingStructures!` with three+ BulkMonomers:** Pairwise syncing 
   must correctly merge all coupling values without overwriting. The existing 
   "skip if already nonzero" logic in `_syncCouplingStructures!` could mask 
   bugs. Need to verify that after sync, all BulkMonomers point to the same 
   `couplingMat` object.

2. **`cplpq(j, j, pos)` when j is a PolyMonomer:** This degenerates to 
   `(pLatS[pos] == j.intIdx) * (pLatS[pos] == j.intIdx)` = 
   `(pLatS[pos] == j.intIdx)`, which is correct (single occupancy check). 
   But the two-argument `couplingPot!` with a PolyMonomer partner needs to 
   store entries at `[self.couplingIdx, pm.intIdx, pm.intIdx]`, and the 
   PolyMonomer side has no coupling structures. This is fine because the 
   symmetry only applies to BulkMonomer participants.

3. **Interaction between `intMat` and `couplingMat`:** Two-way on-site 
   coupling (`couplingMat[i,j,j]`) is physically distinct from nearest-neighbor 
   interaction (`intMat[i,j]`). Tests must verify they contribute independently 
   and don't interfere.

4. **Existing non-coupling tests:** The deletion of `nnpq` must not break 
   any existing code that calls it. Verify `nnpq` is not used outside of 
   coupling-related code with a `grep_search` before deleting.

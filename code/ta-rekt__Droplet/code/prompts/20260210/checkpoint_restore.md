# Implementation Plan: Checkpoint & Restore from Simulation State

## Overview

This plan implements a checkpoint/restore mechanism that saves the mutable state of all species at the end of a simulation, and restores it at the start of a new simulation. The spec file is always re-run first (re-creating structs, setting interactions), and then `restore!` overwrites positions and dynamic state on top of the freshly-initialized objects.

**Key invariant**: `checkpoint!` followed by `restore!` produces a physically equivalent system — same bead positions, same chain topology, same lattice occupancy, same Hamiltonians (recomputed). Scratch buffers and RNG state are **not** preserved.

**Use case**: Run simulation A → checkpoint → modify `intMat` entries and/or `chemPot` values in the spec → `restore!` → run simulation B from where A left off.

**Format**: Julia `Serialization` (already a dependency). No cross-version or cross-language compatibility needed.

---

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Format | `Serialization.serialize` / `deserialize` | Already a dep, binary, fast, Julia-native round-trip |
| Lattice | Reconstructed from species `pos` | Avoids storing the O(L³) array |
| `intMat` / `intIdx` / `intSpecies` | Saved per-species, re-linked on restore | Simple; same code version assumed |
| `KEYCOUNT` | Recomputed as `max(all keys) + 1` | More robust than saving global state |
| `SPECIES_ID_COUNTER` | Recomputed as `max(all ids) + 1` | Same reasoning |
| `FastBulkMonomer.movableIdxs` | Recomputed from `pos` + `latS` on restore | O(N) rebuild, smaller checkpoint |
| `ham` | Recomputed via `ham!(s)` on restore | Cheap, avoids stale values |
| Scratch buffers | Re-allocated on restore (same as `initSpecies!`) | Transient, never need saving |
| Parameters that may change | `intMat` entries, `chemPot` values only | No concentration, length, or pattern changes between runs |

---

## What Gets Saved vs. Recomputed

### Saved (per Polymer)
| Field | Why |
|-------|-----|
| `pos` (OrderedDict) | Core simulation state — bead positions + chain topology + intIdx per bead |
| `head` | Which key is the chain head |
| `length` | Current polymer length (may differ from initial) |
| `keys_vec` | Key vector (tracks all live keys) |
| `pattern` | Current pattern symbol (`:uniform`, `:custom`, etc.) |
| `customPattern` | The pattern vector if `:custom` |
| `patternTails` | intIdx → tail key mapping |
| `baseIntIdx` | Which intIdx is the base type |
| `id` | To match checkpoint species to live species |

### Saved (per BulkMonomer / FastBulkMonomer)
| Field | Why |
|-------|-----|
| `pos` (OrderedDict) | Core simulation state — bead positions |
| `length` | Number of beads (shouldn't change, but saved for validation) |
| `id` | To match checkpoint species to live species |

### Recomputed on Restore
| Item | How |
|------|-----|
| `latS` (lattice views) | Rebuilt by `init!` or a dedicated `_rebuildLattice!` |
| `lat` (full 4D lattice) | Rebuilt from all species' `pos` |
| `ham` | `ham!(s)` after lattice is populated |
| `KEYCOUNT` | `max(all keys across all pos dicts) + 1` |
| `SPECIES_ID_COUNTER` | `max(all species ids) + 1` |
| `movableIdxs` (FastBulkMonomer) | `initMovableSet!(s)` after lattice is populated |
| Scratch buffers | `tempVec`, `bondVec`, `moveVec`, `newPosi`, `newPosj`, `bi`, `bj`, `moveTypes`, `moveTypeProbs`, `nTypes`, `hastingsFactor`, `deltaH`, `deltaHP`, `freeNNs`, `nnIdxs`, `newPos`, `nnsOld`, `nnsNew` — all re-allocated to zeros/defaults |
| `ghost` species | Re-created (empty Species with latS view on layer 1) |
| `speciesDir` | Rebuilt during lattice reconstruction |
| `name` | Re-assigned (e.g., `"poly1_a"`) |

### NOT Saved, NOT Recomputed (comes from spec re-run)
| Item | Why |
|------|-----|
| `intMat`, `intIdx`, `intSpecies` | Set by `intPot!` calls in the spec — user may have changed them |
| `chemPot` | Set in spec — user may have changed it |
| `intPot` | Set in spec |
| `boxlength`, `maxlength`, `nMonomerTypes` | Structural params, unchanged between runs |
| `monomers` dict (Polymer) | Created by `initMonomers!` + `updateMonomers!` in spec |

---

## Checkpoint File Structure

A single `.jls` file containing a `Dict{Symbol, Any}`:

```julia
Dict(
    :version => 1,                          # schema version (for future use)
    :species => [                            # Vector of per-species dicts
        Dict(
            :type => :Polymer,               # or :BulkMonomer, :FastBulkMonomer
            :id => Int32(3),
            :pos => OrderedDict{...},
            :head => Int32(1),
            :length => Int32(150),
            :keys_vec => [...],
            :pattern => :custom,
            :customPattern => [...],
            :patternTails => Dict{...},
            :baseIntIdx => Int32(1),
        ),
        Dict(
            :type => :BulkMonomer,
            :id => Int32(5),
            :pos => OrderedDict{...},
            :length => Int32(1110),
        ),
        Dict(
            :type => :FastBulkMonomer,
            :id => Int32(7),
            :pos => OrderedDict{...},
            :length => Int32(1480),
        ),
    ]
)
```

---

## API

```julia
# Save checkpoint at end of simulation
checkpoint!("path/to/checkpoint", [poly, tf, rnap])

# Restore at start of new simulation (after spec re-run up to init!)
# replaces init!(l, slist) in the spec
restore!("path/to/checkpoint", l, [poly, tf, rnap])
```

### `checkpoint!(filename, slist)`
1. Build the checkpoint dict from each species' mutable state
2. Serialize to `filename.jls`

### `restore!(filename, l, slist)`
1. Deserialize checkpoint dict
2. Match checkpoint species to live species by `id`
3. Overwrite `pos`, `head`, `length`, `keys_vec`, pattern fields
4. Allocate scratch buffers (same as `initSpecies!` does)
5. Rebuild lattice (`lat`, `latS` views, `speciesDir`)
6. Populate lattice from all species' `pos`
7. Rebuild `movableIdxs` for `FastBulkMonomer`
8. Recompute `KEYCOUNT` and `SPECIES_ID_COUNTER`
9. Recompute Hamiltonians via `ham!(s)`

---

## Phase 1: `checkpoint!` — Save Mutable State

### 1.1 Add `checkpoint!` function

**File**: `src/port.jl`

```julia
function checkpoint!(filename::String, slist::Vector{<:AbstractSpecies})
```

Iterates over `slist`, extracts the mutable state fields listed above into a `Dict`, and serializes to `"$filename.jls"`.

**Polymer fields to save**: `id`, `pos`, `head`, `length`, `keys_vec`, `pattern`, `customPattern`, `patternTails`, `baseIntIdx`.

**BulkMonomer / FastBulkMonomer fields to save**: `id`, `pos`, `length`.

**Type tag**: Store `:Polymer`, `:BulkMonomer`, or `:FastBulkMonomer` so restore can validate type matching.

### 1.2 Export `checkpoint!`

**File**: `src/port.jl`

Add to existing export statement.

### 1.3 Tests for Phase 1

```julia
@testset "checkpoint! saves mutable state" begin

    @testset "Checkpoint creates file" begin
        # Set up a simple system
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=40, nMonomerTypes=2,
                       pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")

        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pattern = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pattern, 0.0); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        # Checkpoint
        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        @test isfile("$tmpfile.jls")
        rm("$tmpfile.jls")
    end

    @testset "Checkpoint contains correct species count" begin
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=40, nMonomerTypes=2,
                       pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")

        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pattern = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pattern, 0.0); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        data = Serialization.deserialize("$tmpfile.jls")
        @test length(data[:species]) == 2
        rm("$tmpfile.jls")
    end

    @testset "Checkpoint stores Polymer fields correctly" begin
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2,
                       pattern=:uniform, initShape="line")

        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pat, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        # Apply a custom pattern to have non-trivial patternTails
        newPattern = vcat(fill(base.intIdx, 6), fill(pat.intIdx, 8), fill(base.intIdx, 6))
        updatePattern!(l, poly, newPattern)

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly])

        data = Serialization.deserialize("$tmpfile.jls")
        sp = data[:species][1]

        @test sp[:type] == :Polymer
        @test sp[:id] == poly.id
        @test sp[:head] == poly.head
        @test sp[:length] == poly.length
        @test sp[:pos] == poly.pos
        @test sp[:keys_vec] == poly.keys_vec
        @test sp[:pattern] == poly.pattern
        @test sp[:customPattern] == poly.customPattern
        @test sp[:patternTails] == poly.patternTails
        @test sp[:baseIntIdx] == poly.baseIntIdx

        rm("$tmpfile.jls")
    end

    @testset "Checkpoint stores BulkMonomer fields correctly" begin
        l = Lattice(boxlength=32, nspecies=1)
        bulk = BulkMonomer(length=100, initShape="uniform")
        intPot!(bulk, -0.5)
        init!(l, [bulk])

        tmpfile = tempname()
        checkpoint!(tmpfile, [bulk])

        data = Serialization.deserialize("$tmpfile.jls")
        sp = data[:species][1]

        @test sp[:type] == :BulkMonomer
        @test sp[:id] == bulk.id
        @test sp[:pos] == bulk.pos
        @test sp[:length] == bulk.length

        rm("$tmpfile.jls")
    end

    @testset "Checkpoint stores FastBulkMonomer fields correctly" begin
        l = Lattice(boxlength=32, nspecies=1)
        fbulk = FastBulkMonomer(length=100, initShape="uniform")
        intPot!(fbulk, -0.5)
        init!(l, [fbulk])

        tmpfile = tempname()
        checkpoint!(tmpfile, [fbulk])

        data = Serialization.deserialize("$tmpfile.jls")
        sp = data[:species][1]

        @test sp[:type] == :FastBulkMonomer
        @test sp[:id] == fbulk.id
        @test sp[:pos] == fbulk.pos
        @test sp[:length] == fbulk.length

        rm("$tmpfile.jls")
    end

    @testset "Checkpoint after simulation evolution captures evolved state" begin
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=10, maxlength=30, nMonomerTypes=1,
                       pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")

        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, -0.5); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        # Run some moves so state evolves
        for _ in 1:5000
            move!(poly)
            move!(bulk)
        end

        evolved_pos = deepcopy(poly.pos)
        evolved_length = poly.length

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        data = Serialization.deserialize("$tmpfile.jls")
        sp = data[:species][1]

        @test sp[:pos] == evolved_pos
        @test sp[:length] == evolved_length

        rm("$tmpfile.jls")
    end
end
```

---

## Phase 2: `_restoreSpecies!` — Overwrite Mutable State onto Live Objects

### 2.1 Add `_restoreSpecies!` helper functions

**File**: `src/port.jl`

Three internal methods dispatched on species type:

```julia
function _restoreSpecies!(p::Polymer, saved::Dict{Symbol, Any})
```
Overwrites: `pos`, `head`, `length`, `keys_vec`, `pattern`, `customPattern`, `patternTails`, `baseIntIdx`.
Allocates scratch buffers: `tempVec`, `bondVec`, `moveVec`, `newPosi`, `newPosj`, `moveTypes`, `moveTypeProbs`, `nTypes`, `hastingsFactor`, `deltaH`, `deltaHP`.
Does **not** touch: `monomers`, `intMat`, `chemPot`, `maxlength`, `nMonomerTypes`, `boxlength`.

```julia
function _restoreSpecies!(b::BulkMonomer, saved::Dict{Symbol, Any})
```
Overwrites: `pos`, `length` (validated unchanged).
Allocates scratch buffers: `freeNNs`, `hastingsFactor`, `nnsOld`, `nnsNew`, `newPos`, `tempVec`, `nnIdxs`, `deltaH`, `deltaHP`.

```julia
function _restoreSpecies!(b::FastBulkMonomer, saved::Dict{Symbol, Any})
```
Same as BulkMonomer plus: `movableIdxs` (empty Set, rebuilt later), `bi`, `moveVec`, `moveIdxNew`.

### 2.2 ID matching and type validation

`_restoreSpecies!` should:
- Verify `saved[:id] == species.id`
- Verify `saved[:type]` matches the Julia type of the live object
- Error with a clear message if either fails

### 2.3 Tests for Phase 2

```julia
@testset "_restoreSpecies! overwrites mutable state" begin

    @testset "Polymer pos and chain topology restored" begin
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2,
                       pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pat, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        # Run some moves
        for _ in 1:2000; move!(poly); end

        # Capture evolved state
        saved = Dict{Symbol, Any}(
            :type => :Polymer,
            :id => poly.id,
            :pos => deepcopy(poly.pos),
            :head => poly.head,
            :length => poly.length,
            :keys_vec => copy(poly.keys_vec),
            :pattern => poly.pattern,
            :customPattern => poly.customPattern === nothing ? nothing : copy(poly.customPattern),
            :patternTails => poly.patternTails === nothing ? nothing : copy(poly.patternTails),
            :baseIntIdx => poly.baseIntIdx,
        )

        # Re-create a fresh polymer (simulating spec re-run)
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2,
                        pattern=:uniform, initShape="line")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]; pat2 = poly2.monomers[2]
        intPot!(base2, 0.0); intPot!(pat2, 0.0)
        updateMonomers!(poly2, base2.intIdx)

        # Hack: give poly2 the same id as poly for matching
        # (In practice, spec re-run produces same ids if same creation order)

        Droplet._restoreSpecies!(poly2, saved)

        @test poly2.pos == saved[:pos]
        @test poly2.head == saved[:head]
        @test poly2.length == saved[:length]
        @test poly2.keys_vec == saved[:keys_vec]
    end

    @testset "Polymer scratch buffers re-allocated" begin
        # After _restoreSpecies!, scratch buffers should be freshly allocated
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1,
                       pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        saved = Dict{Symbol, Any}(
            :type => :Polymer, :id => poly.id,
            :pos => deepcopy(poly.pos), :head => poly.head,
            :length => poly.length, :keys_vec => copy(poly.keys_vec),
            :pattern => poly.pattern, :customPattern => nothing,
            :patternTails => nothing, :baseIntIdx => poly.baseIntIdx,
        )

        Droplet._restoreSpecies!(poly, saved)

        @test poly.tempVec !== nothing && length(poly.tempVec) == 3
        @test poly.bondVec !== nothing && length(poly.bondVec) == 3
        @test poly.moveVec !== nothing && length(poly.moveVec) == 3
        @test poly.newPosi !== nothing && length(poly.newPosi) == 3
        @test poly.newPosj !== nothing && length(poly.newPosj) == 3
        @test poly.moveTypes !== nothing
        @test poly.moveTypeProbs !== nothing
        @test poly.deltaHP == 0.0
    end

    @testset "BulkMonomer pos restored" begin
        l = Lattice(boxlength=32, nspecies=1)
        bulk = BulkMonomer(length=50, initShape="uniform")
        intPot!(bulk, -0.5)
        init!(l, [bulk])

        # Run moves
        for _ in 1:1000; move!(bulk); end

        saved = Dict{Symbol, Any}(
            :type => :BulkMonomer, :id => bulk.id,
            :pos => deepcopy(bulk.pos), :length => bulk.length,
        )

        bulk2 = BulkMonomer(length=50, initShape="uniform")
        intPot!(bulk2, -0.5)

        Droplet._restoreSpecies!(bulk2, saved)

        @test bulk2.pos == saved[:pos]
        @test bulk2.deltaH == 0.0
        @test bulk2.deltaHP == 0.0
    end

    @testset "FastBulkMonomer pos restored" begin
        l = Lattice(boxlength=32, nspecies=1)
        fbulk = FastBulkMonomer(length=50, initShape="uniform")
        intPot!(fbulk, -0.5)
        init!(l, [fbulk])

        for _ in 1:1000; move!(fbulk); end

        saved = Dict{Symbol, Any}(
            :type => :FastBulkMonomer, :id => fbulk.id,
            :pos => deepcopy(fbulk.pos), :length => fbulk.length,
        )

        fbulk2 = FastBulkMonomer(length=50, initShape="uniform")
        intPot!(fbulk2, -0.5)

        Droplet._restoreSpecies!(fbulk2, saved)

        @test fbulk2.pos == saved[:pos]
        @test fbulk2.movableIdxs !== nothing  # allocated (empty, rebuilt later)
    end

    @testset "Type mismatch raises error" begin
        saved = Dict{Symbol, Any}(:type => :BulkMonomer, :id => Int32(99),
                                   :pos => OrderedDict{Int32,Vector{Int32}}(),
                                   :length => Int32(10))

        poly = Polymer(length=10, maxlength=20, nMonomerTypes=1,
                       pattern=:uniform, initShape="line")

        @test_throws ErrorException Droplet._restoreSpecies!(poly, saved)
    end

    @testset "ID mismatch raises error" begin
        poly = Polymer(length=10, maxlength=20, nMonomerTypes=1,
                       pattern=:uniform, initShape="line")

        saved = Dict{Symbol, Any}(:type => :Polymer, :id => Int32(poly.id + 100),
                                   :pos => OrderedDict{Int32,Vector{Int32}}(),
                                   :head => Int32(1), :length => Int32(10),
                                   :keys_vec => Int32[], :pattern => :uniform,
                                   :customPattern => nothing,
                                   :patternTails => nothing, :baseIntIdx => nothing)

        @test_throws ErrorException Droplet._restoreSpecies!(poly, saved)
    end

    @testset "chemPot and intMat are NOT overwritten" begin
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1,
                       pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, -0.5)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        saved = Dict{Symbol, Any}(
            :type => :Polymer, :id => poly.id,
            :pos => deepcopy(poly.pos), :head => poly.head,
            :length => poly.length, :keys_vec => copy(poly.keys_vec),
            :pattern => poly.pattern, :customPattern => nothing,
            :patternTails => nothing, :baseIntIdx => poly.baseIntIdx,
        )

        # Change chemPot AFTER checkpoint (simulating spec parameter change)
        base.chemPot = 2.0

        Droplet._restoreSpecies!(poly, saved)

        # chemPot should still be the NEW value (2.0), not the old one
        @test base.chemPot == 2.0
    end
end
```

---

## Phase 3: `_rebuildLattice!` — Reconstruct Lattice from Species State

### 3.1 Add `_rebuildLattice!` function

**File**: `src/port.jl` (or `src/init.jl`)

```julia
function _rebuildLattice!(l::Lattice, slist::Vector{<:AbstractSpecies})
```

This function does what `init!` does for lattice setup, but **without** calling `initSpecies!` (positions are already restored). Specifically:

1. Allocate `l.lat = zeros(...)` (fresh 4D array)
2. Build `l.speciesDir`
3. Create ghost species with `latS` view on layer 1
4. For each species:
   - Assign `s.latS = view(l.lat, i+1, :, :, :)`
   - Assign `s.ghost = gs`
   - Assign `s.boxlength = l.boxlength`
   - Assign `s.name` (letter suffix)
   - Populate `s.latS` from `s.pos` (same logic as `init!`)
5. For `FastBulkMonomer`: call `initMovableSet!(s)`
6. Recompute `KEYCOUNT` = `max(all keys across all species) + 1`
7. For each species with a `ham` field: call `ham!(s)` to recompute

### 3.2 Handle PolyMonomer in slist

Same filtering as `init!` — if slist contains `PolyMonomer`s, extract their parent `Polymer`s and use those as layered species.

### 3.3 Tests for Phase 3

```julia
@testset "_rebuildLattice! reconstructs lattice from species pos" begin

    @testset "Lattice populated correctly for Polymer" begin
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2,
                       pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pat, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        # Run moves to evolve state
        for _ in 1:2000; move!(poly); end

        # Save current pos
        saved_pos = deepcopy(poly.pos)

        # Wipe lattice
        l.lat .= 0

        # Rebuild
        Droplet._rebuildLattice!(l, [poly])

        # Verify every bead is on the lattice
        key = poly.head
        while key != 0
            x, y, z = poly.pos[key][1], poly.pos[key][2], poly.pos[key][3]
            @test poly.latS[x, y, z] == poly.pos[key][6]  # intIdx
            key = poly.pos[key][5]
        end
    end

    @testset "Lattice populated correctly for FastBulkMonomer" begin
        l = Lattice(boxlength=32, nspecies=1)
        fbulk = FastBulkMonomer(length=100, initShape="uniform")
        intPot!(fbulk, -0.5)
        init!(l, [fbulk])

        for _ in 1:1000; move!(fbulk); end

        l.lat .= 0
        Droplet._rebuildLattice!(l, [fbulk])

        # Verify every bead
        for j in 1:fbulk.length
            x, y, z = fbulk.pos[j][1], fbulk.pos[j][2], fbulk.pos[j][3]
            @test fbulk.latS[x, y, z] == j
        end
    end

    @testset "movableIdxs rebuilt for FastBulkMonomer" begin
        l = Lattice(boxlength=32, nspecies=1)
        fbulk = FastBulkMonomer(length=100, initShape="uniform")
        intPot!(fbulk, -0.5)
        init!(l, [fbulk])

        for _ in 1:1000; move!(fbulk); end
        expected_movable = copy(fbulk.movableIdxs)

        # Wipe and rebuild
        fbulk.movableIdxs = Set{Int32}()
        l.lat .= 0
        Droplet._rebuildLattice!(l, [fbulk])

        @test fbulk.movableIdxs == expected_movable
    end

    @testset "KEYCOUNT recomputed correctly" begin
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=10, maxlength=30, nMonomerTypes=1,
                       pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        # Grow polymer (adds keys)
        growPolymer!(poly)

        maxKey = maximum(keys(poly.pos))

        KEYCOUNT[] = 0  # Wipe global
        Droplet._rebuildLattice!(l, [poly])

        @test KEYCOUNT[] >= maxKey + 1
    end

    @testset "Hamiltonians recomputed after rebuild" begin
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1,
                       pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=100, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, -0.5); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:2000
            move!(poly); move!(bulk)
        end

        # Capture expected ham
        ham!(poly); ham!(bulk)
        expected_ham_poly = poly.ham
        expected_ham_bulk = bulk.ham

        # Wipe hams and lattice
        poly.ham = -999.0; bulk.ham = -999.0
        l.lat .= 0

        Droplet._rebuildLattice!(l, [poly, bulk])

        @test isapprox(poly.ham, expected_ham_poly, atol=1e-10)
        @test isapprox(bulk.ham, expected_ham_bulk, atol=1e-10)
    end

    @testset "Mixed species: Polymer + BulkMonomer + FastBulkMonomer" begin
        l = Lattice(boxlength=32, nspecies=3)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1,
                       pattern=:uniform, initShape="line")
        bm = BulkMonomer(length=30, initShape="uniform")
        fbm = FastBulkMonomer(length=50, initShape="uniform")

        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, 0.0); intPot!(bm, -0.3); intPot!(fbm, -0.5)
        intPot!(base, bm, -0.1); intPot!(base, fbm, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bm, fbm])

        for _ in 1:2000
            move!(poly); move!(bm); move!(fbm)
        end

        l.lat .= 0
        Droplet._rebuildLattice!(l, [poly, bm, fbm])

        # Spot-check: no overlaps
        key = poly.head
        while key != 0
            x, y, z = poly.pos[key][1], poly.pos[key][2], poly.pos[key][3]
            @test poly.latS[x, y, z] != 0
            key = poly.pos[key][5]
        end
    end
end
```

---

## Phase 4: `restore!` — Full End-to-End Restore

### 4.1 Add `restore!` function

**File**: `src/port.jl`

```julia
function restore!(filename::String, l::Lattice, slist::Vector{<:AbstractSpecies})
```

Orchestrates the full restore:
1. Deserialize checkpoint from `"$filename.jls"`
2. Validate species count matches
3. For each species in `slist`, find matching entry in checkpoint by `id`
4. Call `_restoreSpecies!(species, saved_dict)` for each
5. Call `_rebuildLattice!(l, slist)`

### 4.2 Export `restore!`

**File**: `src/port.jl`

Add to existing export statement (alongside existing `save!`, `load!`).

### 4.3 Tests for Phase 4

```julia
@testset "restore! full round-trip" begin

    @testset "checkpoint → restore produces identical positions" begin
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2,
                       pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=100, initShape="uniform")

        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, -0.5); intPot!(pat, -0.5); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        # Evolve
        for _ in 1:5000; move!(poly); move!(bulk); end

        # Snapshot evolved state
        evolved_poly_pos = deepcopy(poly.pos)
        evolved_poly_length = poly.length
        evolved_poly_head = poly.head
        evolved_bulk_pos = deepcopy(bulk.pos)

        # Checkpoint
        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        # Simulate spec re-run: create NEW objects with SAME structure
        # (In practice, user re-runs spec. Here we re-create manually
        #  with the same id by resetting SPECIES_ID_COUNTER)
        SPECIES_ID_COUNTER[] = 0
        KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2,
                        pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=100, initShape="uniform")

        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]; pat2 = poly2.monomers[2]
        intPot!(base2, -0.5); intPot!(pat2, -0.5); intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.1)
        updateMonomers!(poly2, base2.intIdx)

        # Restore INSTEAD of init!
        restore!(tmpfile, l2, [poly2, bulk2])

        # Verify positions match evolved state
        @test poly2.pos == evolved_poly_pos
        @test poly2.length == evolved_poly_length
        @test poly2.head == evolved_poly_head
        @test bulk2.pos == evolved_bulk_pos

        rm("$tmpfile.jls")
    end

    @testset "Restored system can continue running moves" begin
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=30, nMonomerTypes=1,
                       pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")

        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, -0.5); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:2000; move!(poly); move!(bulk); end

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        # Fresh system
        SPECIES_ID_COUNTER[] = 0; KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=30, nMonomerTypes=1,
                        pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]
        intPot!(base2, -0.5); intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.1)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        # Should be able to run more moves without error
        nmoves = 0
        for _ in 1:5000
            if move!(poly2); nmoves += 1; end
            if move!(bulk2); nmoves += 1; end
        end
        @test nmoves > 0

        rm("$tmpfile.jls")
    end

    @testset "Restore with changed intPot uses new potentials" begin
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2,
                       pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")

        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, -0.5); intPot!(pat, -0.3); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:2000; move!(poly); move!(bulk); end
        ham!(poly)
        ham_before = poly.ham

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        # Re-run spec with DIFFERENT intPot
        SPECIES_ID_COUNTER[] = 0; KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2,
                        pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]; pat2 = poly2.monomers[2]

        # Changed potentials
        intPot!(base2, -2.0)  # was -0.5
        intPot!(pat2, -1.5)   # was -0.3
        intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.1)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        # Ham should be different because potentials changed
        @test !isapprox(poly2.ham, ham_before, atol=1e-10)

        # But positions should be the same
        @test poly2.pos == poly.pos

        rm("$tmpfile.jls")
    end

    @testset "Restore with changed chemPot uses new chemPot" begin
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2,
                       pattern=:uniform, initShape="line")

        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pat, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        for _ in 1:1000; move!(poly); end
        ham!(poly)
        ham_before = poly.ham

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly])

        # Re-run spec with different chemPot
        SPECIES_ID_COUNTER[] = 0; KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=1)
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2,
                        pattern=:uniform, initShape="line")
        initMonomers!(poly2, defaultChemPot=5.0)  # was 1.0
        base2 = poly2.monomers[1]; pat2 = poly2.monomers[2]
        intPot!(base2, 0.0); intPot!(pat2, 0.0)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2])

        # Ham uses new chemPot
        @test !isapprox(poly2.ham, ham_before, atol=1e-10)

        rm("$tmpfile.jls")
    end

    @testset "Heterogeneous polymer: pattern state survives round-trip" begin
        l = Lattice(boxlength=64, nspecies=1)
        poly = Polymer(length=10, maxlength=30, nMonomerTypes=2,
                       pattern=:uniform, initShape="line")

        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pat, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        # Grow and assign pattern
        growPolymer!(poly)
        patternVec = vcat(fill(base.intIdx, 10), fill(pat.intIdx, 10), fill(base.intIdx, 10))
        updatePattern!(l, poly, patternVec)

        @test poly.patternTails !== nothing
        saved_tails = copy(poly.patternTails)
        saved_pattern = copy(poly.customPattern)

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly])

        # Fresh system
        SPECIES_ID_COUNTER[] = 0; KEYCOUNT[] = 0
        l2 = Lattice(boxlength=64, nspecies=1)
        poly2 = Polymer(length=10, maxlength=30, nMonomerTypes=2,
                        pattern=:uniform, initShape="line")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]; pat2 = poly2.monomers[2]
        intPot!(base2, 0.0); intPot!(pat2, 0.0)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2])

        @test poly2.patternTails == saved_tails
        @test poly2.customPattern == saved_pattern
        @test poly2.pattern == :custom
        @test poly2.baseIntIdx == base.intIdx

        rm("$tmpfile.jls")
    end
end
```

---

## Phase 5: Integration Test — Full Simulation Round-Trip

### 5.1 End-to-end test mimicking the spec workflow

This phase has no new code — only tests that exercise the full `spec → run → checkpoint → modify spec → restore → run` workflow, validating physical consistency.

### 5.2 Tests for Phase 5

```julia
@testset "Full simulation round-trip integration" begin

    @testset "Hamiltonian continuity across checkpoint" begin
        # Run simulation A for N moves, checkpoint at the end.
        # Restore into simulation B with SAME parameters.
        # Ham immediately after restore should match ham at end of A.

        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2,
                       pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=100, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, -0.5); intPot!(pat, -0.3); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1); intPot!(pat, bulk, -0.2)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:10000; move!(poly); move!(bulk); end

        ham!(poly); ham!(bulk)
        ham_poly_end = poly.ham
        ham_bulk_end = bulk.ham

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        # Restore with SAME parameters
        SPECIES_ID_COUNTER[] = 0; KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2,
                        pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=100, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]; pat2 = poly2.monomers[2]
        intPot!(base2, -0.5); intPot!(pat2, -0.3); intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.1); intPot!(pat2, bulk2, -0.2)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        @test isapprox(poly2.ham, ham_poly_end, atol=1e-10)
        @test isapprox(bulk2.ham, ham_bulk_end, atol=1e-10)

        rm("$tmpfile.jls")
    end

    @testset "deltaHam consistency after restore" begin
        # After restoring, deltaHam should still be consistent with
        # full Hamiltonian recalculation across multiple moves.

        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1,
                       pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, -0.5); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.2)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:3000; move!(poly); move!(bulk); end

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        SPECIES_ID_COUNTER[] = 0; KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=1,
                        pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]
        intPot!(base2, -0.5); intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.2)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        # Run more moves and verify deltaHam consistency
        for _ in 1:500
            ham!(poly2)
            hamBefore = poly2.ham
            if move!(poly2)
                ham!(poly2)
                # deltaH from move should match actual ham difference
                @test isapprox(hamBefore + poly2.deltaH, poly2.ham, atol=1e-10)
            end
        end

        rm("$tmpfile.jls")
    end

    @testset "Lattice self-consistency after restore + moves" begin
        # After restoring and running more moves, every bead's
        # lattice entry should still match its pos.

        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=30, nMonomerTypes=2,
                       pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, -0.5); intPot!(pat, -0.3); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:3000; move!(poly); move!(bulk); end

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        SPECIES_ID_COUNTER[] = 0; KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=30, nMonomerTypes=2,
                        pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]; pat2 = poly2.monomers[2]
        intPot!(base2, -0.5); intPot!(pat2, -0.3); intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.1)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        # Run more moves
        for _ in 1:5000; move!(poly2); move!(bulk2); end

        # Verify polymer lattice consistency
        key = poly2.head
        while key != 0
            x, y, z = poly2.pos[key][1], poly2.pos[key][2], poly2.pos[key][3]
            @test poly2.latS[x, y, z] == poly2.pos[key][6]
            key = poly2.pos[key][5]
        end

        # Verify bulk lattice consistency
        for j in 1:bulk2.length
            x, y, z = bulk2.pos[j][1], bulk2.pos[j][2], bulk2.pos[j][3]
            @test bulk2.latS[x, y, z] == j
        end

        rm("$tmpfile.jls")
    end

    @testset "Multiple checkpoint-restore cycles" begin
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1,
                       pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, -0.5); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        tmpfile = tempname()

        for cycle in 1:3
            # Run some moves
            for _ in 1:2000; move!(poly); move!(bulk); end

            # Checkpoint
            checkpoint!(tmpfile, [poly, bulk])

            # Restore
            SPECIES_ID_COUNTER[] = 0; KEYCOUNT[] = 0
            l = Lattice(boxlength=32, nspecies=2)
            poly = Polymer(length=20, maxlength=20, nMonomerTypes=1,
                           pattern=:uniform, initShape="line")
            bulk = FastBulkMonomer(length=50, initShape="uniform")
            initMonomers!(poly, defaultChemPot=1.0)
            base = poly.monomers[1]
            intPot!(base, -0.5); intPot!(bulk, -0.5)
            intPot!(base, bulk, -0.1)
            updateMonomers!(poly, base.intIdx)

            restore!(tmpfile, l, [poly, bulk])

            # Verify can still run
            nmoves = 0
            for _ in 1:1000
                if move!(poly); nmoves += 1; end
                if move!(bulk); nmoves += 1; end
            end
            @test nmoves > 0
        end

        rm("$tmpfile.jls")
    end
end
```

---

## Summary

| Phase | Description | File(s) Modified | New Functions |
|-------|-------------|------------------|---------------|
| 1 | Save mutable state to file | `src/port.jl` | `checkpoint!` |
| 2 | Overwrite mutable state onto live objects | `src/port.jl` | `_restoreSpecies!` (×3 dispatches) |
| 3 | Reconstruct lattice from species positions | `src/port.jl` | `_rebuildLattice!` |
| 4 | Full restore orchestrator | `src/port.jl` | `restore!` |
| 5 | Integration tests (no new code) | `test/` only | — |

### Files Modified

| File | Changes |
|------|---------|
| `src/port.jl` | Add `checkpoint!`, `_restoreSpecies!` (3 methods), `_rebuildLattice!`, `restore!`; update exports |
| `test/test_checkpoint.jl` | New test file with all phases |
| `test/runtests.jl` | Add `include("test_checkpoint.jl")` |

### Unchanged Files

| File | Why Unchanged |
|------|---------------|
| `src/types.jl` | No new struct fields needed |
| `src/init.jl` | `init!` logic reused internally by `_rebuildLattice!`, not modified |
| `src/interactions.jl` | Interaction setup comes from spec re-run |
| `src/hamiltonian.jl` | `ham!` called as-is after restore |
| `src/move.jl` | No changes needed |

### Usage Pattern in Spec Files

```julia
# ===== Simulation A =====
# ... (create species, set intPot!, init!) ...
# ... (run simulation loop) ...
checkpoint!("data/checkpoint_001", [poly, tf, rnap])

# ===== Simulation B (new Julia session or same) =====
# ... (create species, set NEW intPot!/chemPot values) ...
# REPLACE: init!(l, [poly, tf, rnap])
# WITH:
restore!("data/checkpoint_001", l, [poly, tf, rnap])
# ... (continue simulation loop) ...
```

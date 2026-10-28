
export buildOffsets!, buildMoves!, computeRate!, deltaH_kmc, decodeSlot, isOccupied,
       fenwick_build!, fenwick_update!, fenwick_prefix, fenwick_total, fenwick_find,
       sampleMove, applyMove!, refreshAffected!, rebuildQ!,
       kmc_step!, kmc_run!,
       nslots

# direction encoding for KMC hops (axis ∈ 1..3, sign ∈ {+1,-1})
const KMC_DIRS = ((Int32(1), Int32(+1)), (Int32(1), Int32(-1)),
                  (Int32(2), Int32(+1)), (Int32(2), Int32(-1)),
                  (Int32(3), Int32(+1)), (Int32(3), Int32(-1)))

# slot count per species:
#   BulkMonomer: 6 hops × length
#   Polymer:     4 dirs × maxlength (one slot batch per stable bondIdx, with
#                bondIdx ∈ 0..maxlength-1 reserved one-per-bead so add/remove
#                never resize. The chain tail's 4 slots are always inactive,
#                and any unallocated bondIdx contributes 4 zero-rate slots.)
@inline nslots(b::BulkMonomer) = 6 * Int(b.length)
@inline nslots(p::Polymer)     = 4 * Int(p.maxlength)


############# Phase 1: rate primitives (BulkMonomer-only) #############
# MARK: -Phase1

# fills species_offset and nslots from each species' nslots(s)
function buildOffsets!(sys::KMCSystem)
    n = length(sys.slist)
    resize!(sys.species_offset, n)
    cum = 0
    for (i, s) in enumerate(sys.slist)
        sys.species_offset[i] = cum
        cum += nslots(s)
    end
    sys.nslots = cum
    return sys
end

# slot id → (speciesIdx, slot-within-species sub-index, direction).
# For BulkMonomer the second element is `particleId ∈ 1..length` and dir ∈ 1..6.
# For Polymer  the second element is `bondIdx ∈ 0..maxlength-2`     and dir ∈ 1..4.
@inline function decodeSlot(sys::KMCSystem, slot::Int)
    sIdx = 1
    @inbounds for i in length(sys.species_offset):-1:1
        if slot > sys.species_offset[i]
            sIdx = i
            break
        end
    end
    rel = slot - sys.species_offset[sIdx] - 1
    s = sys.slist[sIdx]
    if s isa Polymer
        return sIdx, div(rel, 4), mod(rel, 4) + 1
    else
        return sIdx, div(rel, 6) + 1, mod(rel, 6) + 1
    end
end

# pure ΔH for a unit hop of BulkMonomer b from src to dst (no struct mutation
# beyond nnp's transient pos arithmetic, which restores src/dst on exit)
function deltaH_kmc(b::BulkMonomer, src::Vector{Int32}, dst::Vector{Int32})
    @assert b.couplingMat === nothing "deltaH_kmc: couplings not supported in v1"
    dh = 0.0
    for (sp, idx) in b.intSpecies
        dn = nnp(b, sp, dst) - nnp(b, sp, src)
        # the moving particle still sits at src, so it is counted as a self
        # neighbor of dst (src ∈ NN(dst) for any unit hop). Subtract 1.
        if b.intIdx == idx
            dn -= 1
        end
        dh += b.intMat[b.intIdx, idx] * dn
    end
    return dh
end

# is `pos` occupied by any species in slist?
@inline function isOccupied(sys::KMCSystem, pos::Vector{Int32})
    @inbounds for s in sys.slist
        if s.latS[pos[1], pos[2], pos[3]] > 0
            return true
        end
    end
    return false
end

# compute rate for a single slot (no side effects on sys.rates / sys.Q)
function computeRate!(sys::KMCSystem, slot::Int)
    sIdx, x, dir = decodeSlot(sys, slot)
    s = sys.slist[sIdx]
    if s isa BulkMonomer
        return _computeRate_bulk!(sys, s, Int(x), Int(dir))
    elseif s isa Polymer
        return _computeRate_polymer!(sys, s, Int(x), Int(dir))
    else
        error("KMC v2 supports BulkMonomer and uniform Polymer only")
    end
end

function _computeRate_bulk!(sys::KMCSystem, s::BulkMonomer, particleId::Int, dir::Int)
    pos_arr = s.pos[Int32(particleId)]
    sys.src[1] = pos_arr[1]; sys.src[2] = pos_arr[2]; sys.src[3] = pos_arr[3]

    axis, sgn = KMC_DIRS[dir]
    sys.dst[1] = sys.src[1]; sys.dst[2] = sys.src[2]; sys.dst[3] = sys.src[3]
    sys.dst[axis] += sgn
    boxmod!(sys.dst, sys.dst, Int32(1), s.boxlength)

    # self-only blocking: per the existing Metropolis convention, only the
    # mover species' own lattice excludes (cross-species particles can co-occupy).
    if s.latS[sys.dst[1], sys.dst[2], sys.dst[3]] > 0
        return 0.0
    end

    dh = deltaH_kmc(s, sys.src, sys.dst)
    return s.k0 * min(1.0, exp(-sys.beta * dh))
end

function _computeRate_polymer!(sys::KMCSystem, p::Polymer, bondIdx::Int, dir::Int)
    @assert p.nMonomerTypes == 1 "KMC v2: uniform polymers only (got nMonomerTypes=$(p.nMonomerTypes))"

    state = sys.polymer_states[p]
    if bondIdx + 1 > length(state.key_at_bondIdx)
        return 0.0
    end
    ki = state.key_at_bondIdx[bondIdx + 1]
    if ki == 0
        return 0.0   # bondIdx not currently assigned to any bead
    end

    identifyBond!(p; ki=ki)
    if p.kj == 0
        return 0.0   # bead is the chain tail, no bond
    end

    moveType = _polyMoveVec!(p, dir)

    # Validity (intra-polymer overlap; cross-species co-occupation is allowed).
    if moveType == 1                                     # add
        if p.length >= p.maxlength - 1
            return 0.0
        end
        if _polymer_blocked(p, p.newPosi, p.bprev) ||
           _polymer_blocked(p, p.newPosj, p.bnext)
            return 0.0
        end
    elseif moveType == 2                                 # kink
        moving = _kink_moving_pos(p)
        merge_  = _kink_merge_pos(p)
        if _polymer_blocked(p, moving, merge_)
            return 0.0
        end
    elseif moveType == 3                                 # remove
        if p.length <= 2
            return 0.0
        end
    else
        return 0.0
    end

    dh = deltaH_kmc(p, ki, moveType)

    k0 = moveType == 1 ? p.k0_add :
         moveType == 2 ? p.k0_kink :
                          p.k0_remove
    return k0 * min(1.0, exp(-sys.beta * dh))
end

# A polymer hop is blocked iff `dst` is occupied by ANOTHER bead of the same
# polymer (mover.latS > 0) AND `dst` does not coincide with the merge target
# (bprev/bnext) — the merge case is the kink/remove geometry where the new
# position legitimately overlays an adjacent chain neighbor.
@inline function _polymer_blocked(p::Polymer, dst::Vector{Int32}, mergePos::Vector{Int32})
    return p.latS[dst[1], dst[2], dst[3]] > 0 && !coords_equal(dst, mergePos)
end

# For a kink, identify the new position of the moving bead and the chain
# neighbor it merges past. Caller has already populated newPosi/newPosj/bprev/bnext.
@inline function _kink_moving_pos(p::Polymer)
    return p.newPosi == p.bprev ? p.newPosj : p.newPosi
end
@inline function _kink_merge_pos(p::Polymer)
    return p.newPosi == p.bprev ? p.bprev : p.bnext
end

# Sets p.moveVec, p.newPosi, p.newPosj for the given direction (1..4) and
# returns the moveType (1 add, 2 kink, 3 remove). Mirrors moveTypes! /
# proposeMove! but for one specific direction without random sampling.
function _polyMoveVec!(p::Polymer, dir::Int)
    dims = findall(iszero, p.bondVec)
    @assert length(dims) == 2 "bondVec must be aligned with x, y, or z axis"

    idx = COORD_MOVES[dir]           # 1, 1, 2, 2
    sgn = MOVE_ORIENTATIONS[dir]     # +1, -1, +1, -1

    p.moveVec[1] = 0; p.moveVec[2] = 0; p.moveVec[3] = 0
    p.moveVec[dims[idx]] = sgn

    # newPosi = bi + moveVec
    for j in 1:3
        p.tempVec[j] = p.bi[j] + p.moveVec[j]
    end
    boxmod!(p.tempVec, p.newPosi, Int32(1), p.boxlength)

    # newPosj = bj + moveVec
    for j in 1:3
        p.tempVec[j] = p.bj[j] + p.moveVec[j]
    end
    boxmod!(p.tempVec, p.newPosj, Int32(1), p.boxlength)

    return Int32(p.newPosi == p.bprev) + Int32(p.newPosj == p.bnext) + 1
end

# Pure ΔH for a polymer move. Caller has already populated p.bi/bj/bprev/bnext/
# newPosi/newPosj via identifyBond! + _polyMoveVec! and verified validity.
# Mirrors deltaHam!(p::Polymer) for the uniform-polymer branches only — no
# heterogeneous swap, no internal cross monomers, no struct mutation of energy
# fields (p.deltaH/p.ham/partner.deltaHP).
function deltaH_kmc(p::Polymer, ki::Int32, moveType::Int)
    @assert p.nMonomerTypes == 1 "deltaH_kmc(Polymer): uniform polymers only"
    m = p.monomers[p.pos[ki][6]]
    @assert !isnothing(m.intMat) "deltaH_kmc(Polymer): intMat unset"

    dh = 0.0
    if moveType == 1                                          # bond addition
        for (sp, idx) in m.intSpecies
            if m.intIdx == idx
                dh += -2 * m.chemPot
                dh += m.intPot * (nnp(m, m, p.newPosi) +
                                  nnp(m, m, p.newPosj) + 1)
            else
                dh += m.intMat[m.intIdx, idx] *
                      (nnp(p, sp, p.newPosi) + nnp(p, sp, p.newPosj))
            end
        end
    elseif moveType == 2                                       # kink
        for (sp, idx) in m.intSpecies
            if p.newPosi == p.bprev
                # bead ki moves from bi to newPosj
                if m.intIdx == idx
                    dh += m.intPot *
                          (nnp(m, m, p.newPosj) - nnp(m, m, p.bi))
                else
                    dh += m.intMat[m.intIdx, idx] *
                          (nnp(p, sp, p.newPosj) - nnp(p, sp, p.bi))
                end
            elseif p.newPosj == p.bnext
                # bead kj moves from bj to newPosi
                if m.intIdx == idx
                    dh += m.intPot *
                          (nnp(m, m, p.newPosi) - nnp(m, m, p.bj))
                else
                    dh += m.intMat[m.intIdx, idx] *
                          (nnp(p, sp, p.newPosi) - nnp(p, sp, p.bj))
                end
            end
        end
    elseif moveType == 3                                       # bond removal
        for (sp, idx) in m.intSpecies
            if m.intIdx == idx
                dh += 2 * m.chemPot
                dh -= m.intPot * (nnp(m, m, p.bi) + nnp(m, m, p.bj) - 1)
            else
                dh += -(m.intMat[m.intIdx, idx] *
                        (nnp(p, sp, p.bi) + nnp(p, sp, p.bj)))
            end
        end
    end
    return dh
end

# Initialize per-polymer KMC state: assign bondIdx 0..L-1 to chain beads in
# chain order; remaining slots in the free list (LIFO so they're consumed
# in increasing order on subsequent addBond! events).
function _initPolymerState!(sys::KMCSystem, p::Polymer)
    M = Int(p.maxlength)
    state = KMCPolymerState(Dict{Int32, Int}(), zeros(Int32, M), Int[])

    key = p.head
    bondIdx = 0
    while key != 0 && bondIdx < M
        state.bondIdx_of_key[key] = bondIdx
        state.key_at_bondIdx[bondIdx + 1] = key
        bondIdx += 1
        key = p.pos[key][5]
    end

    for b in (M-1):-1:bondIdx
        push!(state.free_bondIdxs, b)
    end

    sys.polymer_states[p] = state
    return state
end

# recompute every rate from scratch and rebuild the Fenwick tree.
# Used by buildMoves! at startup and as a periodic drift-correction pass.
function rebuildQ!(sys::KMCSystem)
    @assert sys.nslots > 0 "rebuildQ!: call buildMoves! first"
    for slot in 1:sys.nslots
        sys.rates[slot] = computeRate!(sys, slot)
    end
    fenwick_build!(sys.fenwick, sys.rates)
    sys.Q = fenwick_total(sys.fenwick)
    return sys
end

# allocate slot arrays from slist sizes, init per-polymer state, fill rates + Fenwick.
function buildMoves!(sys::KMCSystem)
    buildOffsets!(sys)
    for s in sys.slist
        if s isa Polymer
            _initPolymerState!(sys, s)
        end
    end
    resize!(sys.rates, sys.nslots)
    fill!(sys.rates, 0.0)
    rebuildQ!(sys)
    return sys
end


############# Phase 2: Fenwick tree + sampling #############
# MARK: -Phase2

# Build Fenwick (BIT) tree in-place from a vector of values. 1-indexed.
function fenwick_build!(tree::Vector{Float64}, vals::Vector{Float64})
    n = length(vals)
    resize!(tree, n)
    copyto!(tree, vals)
    @inbounds for i in 1:n
        j = i + (i & -i)
        if j <= n
            tree[j] += tree[i]
        end
    end
    return tree
end

# Add `delta` to slot i.
@inline function fenwick_update!(tree::Vector{Float64}, i::Int, delta::Float64)
    n = length(tree)
    @inbounds while i <= n
        tree[i] += delta
        i += i & -i
    end
    return tree
end

# Sum of values 1..i.
@inline function fenwick_prefix(tree::Vector{Float64}, i::Int)
    s = 0.0
    @inbounds while i > 0
        s += tree[i]
        i -= i & -i
    end
    return s
end

@inline fenwick_total(tree::Vector{Float64}) = fenwick_prefix(tree, length(tree))

# Smallest i such that prefix(i) > u. With u ∈ [0, total), this skips zero-rate
# slots: any slot with rate 0 leaves the prefix unchanged, so prefix(i) > u
# selects only slots that contribute strictly positive weight.
function fenwick_find(tree::Vector{Float64}, u::Float64)
    n = length(tree)
    bit = 1
    while bit <= n
        bit <<= 1
    end
    bit >>= 1
    i = 0
    @inbounds while bit > 0
        nxt = i + bit
        if nxt <= n && tree[nxt] <= u
            i = nxt
            u -= tree[nxt]
        end
        bit >>= 1
    end
    return i + 1
end

# Sample a slot proportional to its rate. Caller ensures sys.Q > 0.
@inline function sampleMove(sys::KMCSystem)
    @assert sys.Q > 0.0 "sampleMove: total rate is zero"
    return fenwick_find(sys.fenwick, rand() * sys.Q)
end


############# Phase 3: applyMove! + incremental updates #############
# MARK: -Phase3

# inverse direction (KMC_DIRS[d] points the opposite way to KMC_DIRS[reverse_dir(d)])
@inline reverse_dir(d::Integer) = ((d - 1) ⊻ 1) + 1

# refresh every slot affected by a move with centers c1, c2 (the lattice sites
# whose occupancy or chain-pointer pattern changed), then reconcile sys.Q and
# the Fenwick tree by the per-slot deltas.
#
# Bulk slots: NN1(c1) ∪ NN1(c2) seeds source/dest enumeration (existing logic).
# Polymer slots: dispatched by sys.refresh_mode —
#   :fine_polymer   — beads with K.pos within L1 ≤ 2 of c1 or c2 contribute
#                     K's outgoing slot AND K's predecessor's slot.
#   :simple_polymer — every polymer slot is refreshed (correctness reference).
#   :rebuild_all    — full rebuildQ! (handled at the call site, not here).
#
# IMPORTANT: c1/c2 must not alias sys.src or sys.dst once the refresh loop
# starts, since computeRate! overwrites those buffers. The collection phase
# (which reads c1/c2) finishes before any computeRate! call.
function refreshAffected!(sys::KMCSystem, c1::Vector{Int32}, c2::Vector{Int32})
    empty!(sys.refresh_set)
    _collectBulkAround!(sys, c1)
    _collectBulkAround!(sys, c2)
    if !isempty(sys.polymer_states)
        if sys.refresh_mode == :simple_polymer
            _markAllPolymerSlots!(sys)
        else  # :fine_polymer (or any future incremental polymer mode)
            _markFinePolymerSlots!(sys, c1, c2)
        end
    end
    _applyRefresh!(sys)
    return sys
end

# Recompute every slot in sys.refresh_set in place, reconciling sys.Q and the
# Fenwick tree by the per-slot delta. computeRate! has scratch-buffer side
# effects on the species (it reuses identifyBond!'s p.bi/p.bj/...); callers
# must finish reading any aliased buffers before invoking this.
function _applyRefresh!(sys::KMCSystem)
    for slot in sys.refresh_set
        old_rate = sys.rates[slot]
        new_rate = computeRate!(sys, slot)
        if new_rate != old_rate
            sys.rates[slot] = new_rate
            delta = new_rate - old_rate
            fenwick_update!(sys.fenwick, slot, delta)
            sys.Q += delta
        end
    end
    return sys
end

# Enumerate NN1(center) (center plus its 6 face neighbors) and gather the
# bulk slots whose source or destination lands at any of those sites.
# Polymer slot indexing (4 dirs, bondIdx-based) does not match this iteration
# scheme, so polymer entries are skipped here and handled by
# _markFinePolymerSlots! / _markAllPolymerSlots!.
function _collectBulkAround!(sys::KMCSystem, center::Vector{Int32})
    boxlength = first(sys.slist).boxlength
    _addBulkSiteAffected!(sys, center)

    nb = sys.nbuf
    for d in 1:6
        axis, sgn = KMC_DIRS[d]
        nb[1] = center[1]; nb[2] = center[2]; nb[3] = center[3]
        nb[axis] += sgn
        boxmod!(nb, nb, Int32(1), boxlength)
        _addBulkSiteAffected!(sys, nb)
    end
    return sys
end

# at site x: any bulk particle at x has its 6 source-moves added; any bulk
# particle at a face-neighbor of x has its move-toward-x added (dest-at-x).
# Polymer species are skipped — their slots have a different indexing scheme.
function _addBulkSiteAffected!(sys::KMCSystem, x::Vector{Int32})
    boxlength = first(sys.slist).boxlength

    # source-at-x
    for (sIdx, s) in enumerate(sys.slist)
        s isa BulkMonomer || continue
        pid = s.latS[x[1], x[2], x[3]]
        if pid > 0
            base = sys.species_offset[sIdx] + (Int(pid) - 1) * 6
            for d in 1:6
                push!(sys.refresh_set, base + d)
            end
        end
    end

    # dest-at-x: walk NN(x); particles there have a move toward x with dir = reverse_dir(d)
    yp = sys.ybuf
    for d in 1:6
        axis, sgn = KMC_DIRS[d]
        yp[1] = x[1]; yp[2] = x[2]; yp[3] = x[3]
        yp[axis] += sgn
        boxmod!(yp, yp, Int32(1), boxlength)
        rev = reverse_dir(d)
        for (sIdx, s) in enumerate(sys.slist)
            s isa BulkMonomer || continue
            pid = s.latS[yp[1], yp[2], yp[3]]
            if pid > 0
                push!(sys.refresh_set,
                      sys.species_offset[sIdx] + (Int(pid) - 1) * 6 + rev)
            end
        end
    end
    return sys
end

# Add ALL slots of every polymer in slist (correctness fallback — the v2.3
# "simple" refresh mode). Cheap to maintain; useful as a reference oracle.
function _markAllPolymerSlots!(sys::KMCSystem)
    for (sIdx, s) in enumerate(sys.slist)
        s isa Polymer || continue
        base = sys.species_offset[sIdx]
        n = nslots(s)
        for i in 1:n
            push!(sys.refresh_set, base + i)
        end
    end
    return sys
end

# Fine-grained polymer enumeration: a slot keyed by bead K (with successor
# K_succ) depends on lattice occupancy in NN1({K, K_succ, K_newPosi, K_newPosj});
# since the slot's newPos is at L1=1 from K and from K_succ, the dependency on
# any dirty site x reduces to "K or K_succ within L1 ≤ 2 of x". We iterate
# beads, find those within L1 ≤ 2 of c1 or c2 (the move centers), and add
# their outgoing-bond slot AND their predecessor's outgoing-bond slot (which
# has K as its bj). Released-bondIdx slots from a remove are NOT covered here
# (the deleted beads aren't in p.pos anymore) — the caller adds them directly.
function _markFinePolymerSlots!(sys::KMCSystem,
                                c1::Vector{Int32}, c2::Vector{Int32})
    for (sIdx, s) in enumerate(sys.slist)
        s isa Polymer || continue
        state = sys.polymer_states[s]
        base = sys.species_offset[sIdx]
        boxL = Int32(s.boxlength)
        for (k, pos) in s.pos
            if !_pos_within_dist2(pos, c1, c2, boxL)
                continue
            end
            # K's outgoing-bond slot
            if haskey(state.bondIdx_of_key, k)
                b = state.bondIdx_of_key[k]
                base_k = base + b * 4
                for d in 1:4
                    push!(sys.refresh_set, base_k + d)
                end
            end
            # K_pred's outgoing-bond slot (K is its bj)
            kpred = pos[4]
            if kpred != 0 && haskey(state.bondIdx_of_key, kpred)
                bp = state.bondIdx_of_key[kpred]
                base_p = base + bp * 4
                for d in 1:4
                    push!(sys.refresh_set, base_p + d)
                end
            end
        end
    end
    return sys
end

# L1 distance with periodic boundaries on a cubic lattice of side L.
@inline function _pbc_l1(pos::Vector{Int32}, c::Vector{Int32}, L::Int32)
    dx = abs(pos[1] - c[1]); dx = min(dx, L - dx)
    dy = abs(pos[2] - c[2]); dy = min(dy, L - dy)
    dz = abs(pos[3] - c[3]); dz = min(dz, L - dz)
    return dx + dy + dz
end

@inline function _pos_within_dist2(pos::Vector{Int32},
                                   c1::Vector{Int32}, c2::Vector{Int32},
                                   L::Int32)
    return _pbc_l1(pos, c1, L) <= 2 || _pbc_l1(pos, c2, L) <= 2
end

# Apply the hop encoded by `slot`: advance simulation time using the OLD Q,
# then dispatch to the species-specific commit + rate refresh.
function applyMove!(sys::KMCSystem, slot::Int)
    sIdx, x, dir = decodeSlot(sys, slot)
    s = sys.slist[sIdx]

    # Gillespie time advance using the rate distribution that produced this draw
    Q_old = sys.Q
    @assert Q_old > 0.0 "applyMove!: Q == 0"
    sys.t_sim += -log(rand()) / Q_old
    sys.n_events += 1

    if s isa BulkMonomer
        _applyMove_bulk!(sys, sIdx, s, Int(x), Int(dir))
    elseif s isa Polymer
        _applyMove_polymer!(sys, sIdx, s, Int(x), Int(dir))
    else
        error("KMC v2: unsupported species $(typeof(s))")
    end
    return sys
end

function _applyMove_bulk!(sys::KMCSystem, _sIdx::Int, s::BulkMonomer,
                          particleId::Int, dir::Int)
    # source and destination — held in sys.src/sys.dst, then read by refreshAffected!
    # during its collection pass (before any computeRate! call clobbers them).
    pos_arr = s.pos[Int32(particleId)]
    sys.src[1] = pos_arr[1]; sys.src[2] = pos_arr[2]; sys.src[3] = pos_arr[3]

    axis, sgn = KMC_DIRS[dir]
    sys.dst[1] = sys.src[1]; sys.dst[2] = sys.src[2]; sys.dst[3] = sys.src[3]
    sys.dst[axis] += sgn
    boxmod!(sys.dst, sys.dst, Int32(1), s.boxlength)

    # commit the move on lattice + pos
    s.latS[sys.src[1], sys.src[2], sys.src[3]] = Int32(0)
    s.latS[sys.dst[1], sys.dst[2], sys.dst[3]] = Int32(particleId)
    pos_arr[1] = sys.dst[1]; pos_arr[2] = sys.dst[2]; pos_arr[3] = sys.dst[3]

    if sys.refresh_mode == :rebuild_all && !isempty(sys.polymer_states)
        rebuildQ!(sys)
    else
        refreshAffected!(sys, sys.src, sys.dst)
    end
    return sys
end

# Apply a polymer bond move (add / kink / remove) sampled by `bondIdx, dir`.
# Mirrors the Metropolis flow: identifyBond! → _polyMoveVec! → makeMove!. The
# stable bondIdx mapping in sys.polymer_states[p] is updated alongside the
# chain (free list pop on add, push back on remove; kink is mapping-neutral).
# In v2.3 the rate refresh runs incrementally per sys.refresh_mode; the
# :rebuild_all mode preserves the v2.2 behavior as a correctness oracle.
function _applyMove_polymer!(sys::KMCSystem, sIdx::Int, p::Polymer,
                             bondIdx::Int, dir::Int)
    @assert p.nMonomerTypes == 1 "KMC v2: uniform polymers only"

    state = sys.polymer_states[p]
    ki_key = state.key_at_bondIdx[bondIdx + 1]
    @assert ki_key != 0 "applyMove!: sampled an unallocated polymer bondIdx slot"

    identifyBond!(p; ki=ki_key)
    moveType = _polyMoveVec!(p, dir)
    @assert moveType in (1, 2, 3) "applyMove!: invalid moveType $(moveType) from sampled slot"
    p.moveType = Int32(moveType)

    # Snapshot the keys whose bondIdxs we'll release on a remove. removeBond!
    # deletes p.pos[p.ki] and p.pos[p.kj], but the integer keys themselves
    # remain valid handles into the bondIdx_of_key dict.
    old_ki = p.ki
    old_kj = p.kj

    # Stash the move centers (sites whose occupancy / chain pointer pattern
    # changes) into sys.cbuf{1,2}. We must do this BEFORE makeMove! mutates
    # p.bi/p.newPosi/etc. The polymer scratch buffers will be clobbered again
    # by computeRate! during the eventual refresh, so we rely on cbuf{1,2}
    # owned by the system through the refresh.
    c1 = sys.cbuf1; c2 = sys.cbuf2
    if moveType == 1                           # add: newPosi, newPosj newly occupied
        c1[1] = p.newPosi[1]; c1[2] = p.newPosi[2]; c1[3] = p.newPosi[3]
        c2[1] = p.newPosj[1]; c2[2] = p.newPosj[2]; c2[3] = p.newPosj[3]
    elseif moveType == 2                       # kink: moving bead's old & new positions
        if p.newPosi == p.bprev                # ki moves: bi → newPosj
            c1[1] = p.bi[1];      c1[2] = p.bi[2];      c1[3] = p.bi[3]
            c2[1] = p.newPosj[1]; c2[2] = p.newPosj[2]; c2[3] = p.newPosj[3]
        else                                   # kj moves: bj → newPosi
            c1[1] = p.bj[1];      c1[2] = p.bj[2];      c1[3] = p.bj[3]
            c2[1] = p.newPosi[1]; c2[2] = p.newPosi[2]; c2[3] = p.newPosi[3]
        end
    else                                       # moveType == 3 (remove): bi, bj newly empty
        c1[1] = p.bi[1]; c1[2] = p.bi[2]; c1[3] = p.bi[3]
        c2[1] = p.bj[1]; c2[2] = p.bj[2]; c2[3] = p.bj[3]
    end

    # Snapshot the bondIdxs we'll release on remove BEFORE the mapping update
    # — the dict entry is gone afterward.
    released_b1 = -1; released_b2 = -1
    if moveType == 3
        released_b1 = state.bondIdx_of_key[old_ki]
        released_b2 = state.bondIdx_of_key[old_kj]
    end

    makeMove!(p)
    p.length = length(keys(p.pos))

    if moveType == 1
        # addBond! pushes the two new bead keys to the END of p.keys_vec, in
        # chain order (next_key1 then next_key2).
        new_key1 = p.keys_vec[end-1]
        new_key2 = p.keys_vec[end]
        @assert length(state.free_bondIdxs) >= 2 "applyMove!: free bondIdx list exhausted on add"
        b1 = pop!(state.free_bondIdxs)
        b2 = pop!(state.free_bondIdxs)
        state.bondIdx_of_key[new_key1] = b1
        state.bondIdx_of_key[new_key2] = b2
        state.key_at_bondIdx[b1 + 1] = new_key1
        state.key_at_bondIdx[b2 + 1] = new_key2
    elseif moveType == 3
        delete!(state.bondIdx_of_key, old_ki)
        delete!(state.bondIdx_of_key, old_kj)
        state.key_at_bondIdx[released_b1 + 1] = Int32(0)
        state.key_at_bondIdx[released_b2 + 1] = Int32(0)
        push!(state.free_bondIdxs, released_b1)
        push!(state.free_bondIdxs, released_b2)
    end
    # moveType == 2 (kink): bondIdx mapping is unaffected. Bead positions
    # change in-place, keys are stable, chain topology is preserved.

    if sys.refresh_mode == :rebuild_all
        rebuildQ!(sys)
    else
        # Refresh strategy by moveType:
        # - kink (2): p.length unchanged → fine enumeration captures every
        #   slot whose rate can flip (purely local).
        # - add (1) / remove (3): p.length changes, which flips the GLOBAL
        #   validity gates in _computeRate_polymer! ("length >= maxlength - 1"
        #   for add, "length <= 2" for remove). These can flip far-away
        #   slots' rates without any local lattice change, so fine
        #   enumeration alone misses them. We refresh every polymer slot
        #   (across all polymers in slist) for these events. :simple_polymer
        #   always does this regardless of moveType.
        empty!(sys.refresh_set)
        if moveType == 3
            # Explicit: the released bondIdxs' slots' beads are no longer in
            # p.pos, so an enumeration over p.pos would miss them. (Redundant
            # with _markAllPolymerSlots! below, but kept for clarity and so
            # the remove path stays correct if the polymer-marking strategy
            # ever changes again.)
            base = sys.species_offset[sIdx]
            base_b1 = base + released_b1 * 4
            base_b2 = base + released_b2 * 4
            for d in 1:4
                push!(sys.refresh_set, base_b1 + d)
                push!(sys.refresh_set, base_b2 + d)
            end
        end
        _collectBulkAround!(sys, c1)
        _collectBulkAround!(sys, c2)
        if sys.refresh_mode == :simple_polymer || moveType != 2
            _markAllPolymerSlots!(sys)
        else
            _markFinePolymerSlots!(sys, c1, c2)
        end
        _applyRefresh!(sys)
    end
    return sys
end


############# Phase 4: main loop + sim-time sampling #############
# MARK: -Phase4

# one Gillespie event: sample → applyMove! (which advances t_sim and refreshes)
@inline function kmc_step!(sys::KMCSystem)
    applyMove!(sys, sampleMove(sys))
    return sys
end

# Run KMC until any of the caps fire. Sampling is on a fixed grid of simulation
# time t_burn_in, t_burn_in + dt_sample, t_burn_in + 2·dt_sample, ...; whenever
# sys.t_sim crosses the next grid point, sample_callback(sys) is invoked once.
# t_max_wall is in minutes (matches target!() convention).
function kmc_run!(sys::KMCSystem;
                  t_max_wall::Float64=Inf,
                  t_max_sim::Float64=Inf,
                  t_burn_in::Float64=0.0,
                  dt_sample::Float64=Inf,
                  sample_callback=nothing,
                  rebuild_every::Int=1_000_000)
    t_start = time_ns()
    t_next_sample = t_burn_in
    next_rebuild = sys.n_events + rebuild_every
    n_samples = 0

    while sys.t_sim < t_max_sim
        if (time_ns() - t_start) / 6e10 >= t_max_wall
            break
        end
        if sys.Q <= 0.0
            break
        end

        kmc_step!(sys)

        if sys.n_events >= next_rebuild
            rebuildQ!(sys)
            next_rebuild = sys.n_events + rebuild_every
        end

        if sample_callback !== nothing && sys.t_sim >= t_next_sample
            sample_callback(sys)
            n_samples += 1
            t_next_sample += dt_sample
        end
    end

    return n_samples
end

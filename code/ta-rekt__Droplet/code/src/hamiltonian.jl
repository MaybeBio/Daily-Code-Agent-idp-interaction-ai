
export ham!, deltaHam!

############# Polymer Energy functions #############
# MARK: -Polymer

# Hamiltonian of a polymer
function ham!(p::Polymer)
    p.ham = 0
    totalChemEnergy = 0.0
    totalIntEnergy = 0.0

    key = p.head
    
    # Get any monomer to check intMat exists (all share same matrix)
    refMonomer = first(values(p.monomers))
    @assert !isnothing(refMonomer.intMat) "No interaction matrix!"

    for i in 1:p.length
        pos = p.pos[key][1:3]
        m = p.monomers[p.pos[key][6]]  # Dictionary lookup by intIdx
        
        # Chemical potential contribution
        totalChemEnergy -= m.chemPot
        
        for (sp, idx) in m.intSpecies    
            # Interaction energy: count neighbors
            nnp_ = nnp(m, sp, pos)
            
            # For self-interactions, divide by 2 to avoid double-counting bonds
            # For cross-interactions, count all neighbors
            if m.intIdx == idx
                totalIntEnergy += m.intMat[m.intIdx, idx] * nnp_ / 2.0
            else
                totalIntEnergy += m.intMat[m.intIdx, idx] * nnp_
            end
        end

        key = p.pos[key][5]
    end

    @assert key == 0 "Tail is not N'th bead"

    p.ham = totalChemEnergy + totalIntEnergy

    return p
end

# Compute energy change from swapping tail beads from pattern type to base type
# Called during addBond inside a pattern region
function _addSwapEnergy!(p::Polymer, mPattern::PolyMonomer)
    mBase = p.monomers[p.baseIntIdx]
    
    # Get the two tail beads that will be swapped
    tailKey = p.patternTails[mPattern.intIdx]
    prevTailKey = p.pos[tailKey][4]
    tailPos = view(p.pos[tailKey], 1:3)
    prevTailPos = view(p.pos[prevTailKey], 1:3)
    
    # Get positions of new beads being added
    newPosiView = view(p.newPosi, 1:3)
    newPosjView = view(p.newPosj, 1:3)
    
    # Count neighbors at each position BEFORE the swap
    nnpTailPattern = nnp(mPattern, mPattern, tailPos)     # pattern neighbors of tail
    nnpPrevTailPattern = nnp(mPattern, mPattern, prevTailPos)  # pattern neighbors of prev_tail
    nnpTailBase = nnp(mBase, mBase, tailPos)              # base neighbors of tail
    nnpPrevTailBase = nnp(mBase, mBase, prevTailPos)      # base neighbors of prev_tail
    
    # Energy change from self-interactions:
    #
    # BEFORE swap:
    # - Pattern-pattern bonds: (nnpTailPattern + nnpPrevTailPattern - 1) 
    #   (subtract 1 to avoid double-counting the tail-prevTail bond)
    # - Base-base bonds involving tail/prevTail: 0 (they're pattern type)
    #
    # AFTER swap:
    # - The tail-prevTail bond becomes base-base: 1 bond
    # - Tail gains nnpTailBase base-base bonds (existing base neighbors)
    # - PrevTail gains nnpPrevTailBase base-base bonds (existing base neighbors)
    # - Bonds to remaining pattern neighbors become base-pattern CROSS bonds
    #
    # Let's count remaining pattern neighbors (excluding the other swapped bead):
    # - Tail: nnpTailPattern - 1 (exclude prev_tail) = remaining pattern neighbors
    # - PrevTail: nnpPrevTailPattern - 1 (exclude tail) = remaining pattern neighbors
    
    # Pattern-pattern bonds removed (all of them)
    patternBondsRemoved = nnpTailPattern + nnpPrevTailPattern - 1  # -1 for tail-prevTail counted twice
    p.deltaH -= mPattern.intPot * patternBondsRemoved
    
    # Base-base bonds added
    # - 1 for the tail-prevTail bond (now base-base)
    # - nnpTailBase for tail's existing base neighbors
    # - nnpPrevTailBase for prevTail's existing base neighbors
    baseBondsAdded = 1 + nnpTailBase + nnpPrevTailBase
    p.deltaH += mBase.intPot * baseBondsAdded
    
    # Base-pattern cross bonds created
    # - Tail had (nnpTailPattern - 1) pattern neighbors (excluding prev_tail)
    # - Prev_tail had (nnpPrevTailPattern - 1) pattern neighbors (excluding tail)
    # These become base-pattern cross bonds
    crossInteraction = mPattern.intMat[mPattern.intIdx, p.baseIntIdx]
    crossBondsCreated = (nnpTailPattern - 1) + (nnpPrevTailPattern - 1)
    p.deltaH += crossInteraction * crossBondsCreated
    
    # New cross bonds created: remaining pattern neighbors of swapped beads
    # - Tail had (nnpTailPattern - 1) pattern neighbors (excluding prev_tail)
    # - Prev_tail had (nnpPrevTailPattern - 1) pattern neighbors (excluding tail)
    # These become base-pattern cross bonds (counted bidirectionally → ×2)
    crossBondsCreated = (nnpTailPattern - 1) + (nnpPrevTailPattern - 1)
    p.deltaH += 2 * crossInteraction * crossBondsCreated

    # Correction for bonds between new beads (newPosi/newPosj) and swapped tail beads
    # The main deltaHam! already added pattern-pattern energy for these bonds,
    # but after the swap, these bonds become pattern-base cross bonds.
    # We need to correct: subtract pattern-pattern, add cross energy
    bondsToCorrect = 0
    if isneighbor(newPosiView, tailPos)
        bondsToCorrect += 1
    end
    if isneighbor(newPosiView, prevTailPos)
        bondsToCorrect += 1
    end
    if isneighbor(newPosjView, tailPos)
        bondsToCorrect += 1
    end
    if isneighbor(newPosjView, prevTailPos)
        bondsToCorrect += 1
    end
    # Subtract the pattern-pattern energy that was over-counted, add cross energy

    p.deltaH -= mPattern.intPot * bondsToCorrect
    p.deltaH += crossInteraction * bondsToCorrect
    
    # Bulk interactions: swapped beads change their bulk interaction type
    for (sp, idx) in mPattern.intSpecies
        # Skip self-interactions (already handled above)
        if idx == mPattern.intIdx || idx == p.baseIntIdx
            continue
        end
        # This is a bulk species
        nnpTailBulk = nnp(mPattern, sp, tailPos)
        nnpPrevTailBulk = nnp(mPattern, sp, prevTailPos)
        # Remove old pattern-bulk, add new base-bulk
        energy_change = -mPattern.intMat[mPattern.intIdx, idx] * (nnpTailBulk + nnpPrevTailBulk) +
                         mBase.intMat[mBase.intIdx, idx] * (nnpTailBulk + nnpPrevTailBulk)
        p.deltaH += energy_change
        # Notify bulk species of energy change
        setDeltaHP(sp, getDeltaHP(sp) + energy_change)
    end
    
    # Chemical potential change: 2 beads go from pattern to base type
    p.deltaH -= 2 * (mBase.chemPot - mPattern.chemPot)
end

# Compute energy change from swapping beads after tail from base type to pattern type
# Called during removeBond inside a pattern region
function _removeSwapEnergy!(p::Polymer, mPattern::PolyMonomer)
    mBase = p.monomers[p.baseIntIdx]
    
    # Get the two beads after tail that will be swapped to pattern type
    # These beads are currently base type and will become pattern type
    tailKey = p.patternTails[mPattern.intIdx]
    nextAfterTail = p.pos[tailKey][5]
    nextNextAfterTail = p.pos[nextAfterTail][5]
    afterTailPos = view(p.pos[nextAfterTail], 1:3)
    afterAfterTailPos = view(p.pos[nextNextAfterTail], 1:3)
    
    # Get positions of beads being removed (bi, bj)
    biPos = view(p.bi, 1:3)
    bjPos = view(p.bj, 1:3)
    
    # Count neighbors at each position BEFORE the swap
    nnpAfterBase = nnp(mBase, mBase, afterTailPos)           # base neighbors
    nnpAfterAfterBase = nnp(mBase, mBase, afterAfterTailPos) # base neighbors
    nnpAfterPattern = nnp(mPattern, mPattern, afterTailPos)      # pattern neighbors (including tail)
    nnpAfterAfterPattern = nnp(mPattern, mPattern, afterAfterTailPos)
    
    # Exclude bi and bj from pattern neighbor counts (they're being removed!)
    # Check if bi or bj are neighbors of afterTail
    if isneighbor(afterTailPos, biPos)
        nnpAfterPattern -= 1
    end
    if isneighbor(afterTailPos, bjPos)
        nnpAfterPattern -= 1
    end
    # Check if bi or bj are neighbors of afterAfterTail
    if isneighbor(afterAfterTailPos, biPos)
        nnpAfterAfterPattern -= 1
    end
    if isneighbor(afterAfterTailPos, bjPos)
        nnpAfterAfterPattern -= 1
    end
    
    # Energy change from self-interactions:
    #
    # BEFORE swap (beads are base type):
    # - Base-base bonds: (nnpAfterBase + nnpAfterAfterBase - 1)
    #   (subtract 1 to avoid double-counting the after-afterAfter bond)
    # - Pattern-pattern bonds: 0 (they're base type)
    #
    # AFTER swap (beads become pattern type):
    # - The after-afterAfter bond becomes pattern-pattern: 1 bond
    # - after gains nnpAfterPattern pattern-pattern bonds (including current tail)
    # - afterAfter gains nnpAfterAfterPattern pattern-pattern bonds
    # - Bonds to remaining base neighbors become base-pattern CROSS bonds
    
    # Base-base bonds removed (all of them)
    baseBondsRemoved = nnpAfterBase + nnpAfterAfterBase - 1  # -1 for after-afterAfter counted twice
    p.deltaH -= mBase.intPot * baseBondsRemoved
    
    # Pattern-pattern bonds added
    # - 1 for the after-afterAfter bond (now pattern-pattern)
    # - nnpAfterPattern for after's existing pattern neighbors (includes current tail)
    # - nnpAfterAfterPattern for afterAfter's existing pattern neighbors
    patternBondsAdded = 1 + nnpAfterPattern + nnpAfterAfterPattern
    p.deltaH += mPattern.intPot * patternBondsAdded
    
    # Base-pattern cross bonds created
    # - after had (nnpAfterBase - 1) base neighbors (excluding afterAfter)
    # - afterAfter had (nnpAfterAfterBase - 1) base neighbors (excluding after)
    # These become base-pattern cross bonds
    crossInteraction = mPattern.intMat[mPattern.intIdx, p.baseIntIdx]
    crossBondsCreated = (nnpAfterBase - 1) + (nnpAfterAfterBase - 1)
    p.deltaH += crossInteraction * crossBondsCreated
    
    # Base-pattern cross bonds created
    # - after had (nnpAfterBase - 1) base neighbors (excluding afterAfter)
    # - afterAfter had (nnpAfterAfterBase - 1) base neighbors (excluding after)
    # These become base-pattern cross bonds
    # Multiply by 2 because ham! counts cross bonds from BOTH sides (no /2)
    crossBondsCreated = (nnpAfterBase - 1) + (nnpAfterAfterBase - 1)
    p.deltaH += 2 * crossInteraction * crossBondsCreated

    # Bulk interactions: swapped beads change their bulk interaction type
    for (sp, idx) in mPattern.intSpecies
        # Skip self-interactions (already handled above)
        if idx == mPattern.intIdx || idx == p.baseIntIdx
            continue
        end
        # This is a bulk species
        nnpAfterBulk = nnp(mBase, sp, afterTailPos)
        nnpAfterAfterBulk = nnp(mBase, sp, afterAfterTailPos)
        # Remove old base-bulk, add new pattern-bulk
        energy_change = -mBase.intMat[mBase.intIdx, idx] * (nnpAfterBulk + nnpAfterAfterBulk) +
                         mPattern.intMat[mPattern.intIdx, idx] * (nnpAfterBulk + nnpAfterAfterBulk)
        p.deltaH += energy_change
        # Notify bulk species of energy change
        setDeltaHP(sp, getDeltaHP(sp) + energy_change)
    end
    
    # Chemical potential change: 2 beads go from base to pattern type
    p.deltaH -= 2 * (mPattern.chemPot - mBase.chemPot)
end

# Energy increment of a proposed move. Requires proposeMove! to be called
# which sets moveType and newPosi/j. 
function deltaHam!(p::Polymer)

    @assert p.moveType == 1 || p.moveType == 2 || p.moveType == 3 begin
        "p.moveType must be either 1, 2, or 3"
    end
    
    p.deltaH = 0
    p.deltaHP = 0
    
    # Get monomer for the bead being moved (ki)
    m = p.monomers[p.pos[p.ki][6]]
    
    @assert !isnothing(m.intMat) "No interaction matrix!"

    # bond addition
    if p.moveType == 1
        for (sp, idx) in m.intSpecies
            
            # cross-species interactions (different intIdx)
            if m.intIdx != idx  
                
                energy_change = m.intMat[m.intIdx, idx] * nnp(p, sp, p.newPosi) + 
                                m.intMat[m.intIdx, idx] * nnp(p, sp, p.newPosj)

                # For PolyMonomers of the same polymer, use nnp(m, sp, ...) to count
                # only matching intIdx, and multiply by 2 because ham! counts each
                # cross bond from BOTH sides (no /2 for cross terms).
                if sp isa PolyMonomer && sp.poly === p
                    energy_change = 2 * m.intMat[m.intIdx, idx] * (nnp(m, sp, p.newPosi) + 
                                                                    nnp(m, sp, p.newPosj))
                    p.deltaH += energy_change
                else
                    # For external species (bulk monomers), use nnp(p, sp, ...)
                    energy_change = m.intMat[m.intIdx, idx] * (nnp(p, sp, p.newPosi) + 
                                                               nnp(p, sp, p.newPosj))
                    setDeltaHP(sp, energy_change)
                    p.deltaH += getDeltaHP(sp)
                end
            
            # self-interactions
            else
                # Chemical potential: adding 2 beads
                p.deltaH += - 2 * m.chemPot
                
                # Interaction energy: count bonds being added
                # New beads at newPosi and newPosj will be connected (1 bond)
                # New bead at newPosi will connect to nnp(m,m,newPosi) same-type neighbors
                # New bead at newPosj will connect to nnp(m,m,newPosj) same-type neighbors  
                # Total bonds added: 1 + nnp(newPosi) + nnp(newPosj)
                p.deltaH += m.intPot * (nnp(m, m, p.newPosi) + nnp(m, m, p.newPosj) + 1)
            end
        end
        
        # Swap energy for add inside pattern region (only if in contiguous pattern region)
        if m.intIdx != p.baseIntIdx && haskey(p.patternTails, m.intIdx) && isInPatternRegion(p, p.ki, m.intIdx)
            _addSwapEnergy!(p, m)
        end

    # kink
    elseif p.moveType == 2
        for (sp, idx) in m.intSpecies

            # bead ki is moving from bi to newPosj
            if p.newPosi == p.bprev

                # cross-species interactions
                if m.intIdx != idx  

                    # For PolyMonomers of the same polymer, use nnp(m, sp, ...) 
                    # and multiply by 2 for bidirectional counting
                    if sp isa PolyMonomer && sp.poly === p
                        energy_change = 2 * m.intMat[m.intIdx, idx] * (nnp(m, sp, p.newPosj) - nnp(m, sp, p.bi))
                        p.deltaH += energy_change
                    else
                        energy_change = m.intMat[m.intIdx, idx] * (nnp(p, sp, p.newPosj) - nnp(p, sp, p.bi))
                        setDeltaHP(sp, energy_change)
                        p.deltaH += getDeltaHP(sp)
                    end

                # self-interactions
                else
                    # Bead ki moves from bi to newPosj
                    # Energy change = (interactions at newPosj) - (interactions at bi)
                    # Use nnp(m,m,pos) to count only same-type neighbors
                    p.deltaH += m.intPot * (nnp(m, m, p.newPosj) - nnp(m, m, p.bi))
                end

            # bead kj is moving from bj to newPosi
            elseif p.newPosj == p.bnext

                # cross-species interactions
                if m.intIdx != idx  

                    # For PolyMonomers of the same polymer, use nnp(m, sp, ...)
                    # and multiply by 2 for bidirectional counting
                    if sp isa PolyMonomer && sp.poly === p
                        energy_change = 2 * m.intMat[m.intIdx, idx] * (nnp(m, sp, p.newPosi) - nnp(m, sp, p.bj))
                        p.deltaH += energy_change
                    else
                        energy_change = m.intMat[m.intIdx, idx] * (nnp(p, sp, p.newPosi) - nnp(p, sp, p.bj))
                        setDeltaHP(sp, energy_change)
                        p.deltaH += getDeltaHP(sp)
                    end

                # self-interactions
                else
                    # Bead kj moves from bj to newPosi
                    # Energy change = (interactions at newPosi) - (interactions at bj)
                    # Use nnp(m,m,pos) to count only same-type neighbors
                    p.deltaH += m.intPot * (nnp(m, m, p.newPosi) - nnp(m, m, p.bj))
                end

            else
                error("Move is not a kink")
            end
        end

    # bond removal  
    elseif p.moveType == 3
        for (sp, idx) in m.intSpecies
            
            # cross-species interactions
            if m.intIdx != idx  

                energy_change = -(m.intMat[m.intIdx, idx] * nnp(p, sp, p.bi) + 
                                  m.intMat[m.intIdx, idx] * nnp(p, sp, p.bj))

                # For PolyMonomers of the same polymer, use nnp(m, sp, ...)
                # and multiply by 2 for bidirectional counting
                if sp isa PolyMonomer && sp.poly === p
                    energy_change = -2 * (m.intMat[m.intIdx, idx] * (nnp(m, sp, p.bi) + 
                                                                      nnp(m, sp, p.bj)))
                    p.deltaH += energy_change
                else
                    energy_change = -(m.intMat[m.intIdx, idx] * (nnp(p, sp, p.bi) + 
                                                                  nnp(p, sp, p.bj)))
                    setDeltaHP(sp, energy_change)
                    p.deltaH += getDeltaHP(sp)
                end
            
            # self-interactions
            else
                # Chemical potential: removing 2 beads
                p.deltaH += 2 * m.chemPot
                
                # Interaction energy: count bonds being removed
                # Beads bi and bj are connected (1 bond between them)
                # Bead bi has nnp(m,m,bi) - 1 other same-type neighbors (excluding bj)
                # Bead bj has nnp(m,m,bj) - 1 other same-type neighbors (excluding bi)
                # Total bonds removed: 1 + (nnp(bi)-1) + (nnp(bj)-1) = nnp(bi) + nnp(bj) - 1
                p.deltaH -= m.intPot * (nnp(m, m, p.bi) + nnp(m, m, p.bj) - 1)
            end
        end
        
        # Swap energy for remove inside pattern region (only if in contiguous pattern region)
        if m.intIdx != p.baseIntIdx && haskey(p.patternTails, m.intIdx) && isInPatternRegion(p, p.ki, m.intIdx)
            _removeSwapEnergy!(p, m)
        end
    end

    p.ham += p.deltaH

    return p
end

ham!(m::PolyMonomer) = ham!(m.poly)
deltaHam!(m::PolyMonomer) = deltaHam!(m.poly)

############# BulkMonomer energy functions #############
# MARK: -BulkMonomer

# Hamiltonian of a dilute species in the bulkMonomer
function ham!(b::AbstractBulkMonomer)
    
    # chemical potential
    b.ham = - b.chemPot * b.length

    # interaction potentials (including self)
    for (sp, idx) in b.intSpecies
        nnp_ = 0.0   # number of nearest partners
        selfFactor = b.intIdx == idx ? 0.5 : 1.0

        for i in 1:b.length
            pos = b.pos[i][range(1, 3)]
            nnp_ += nnp(b, sp, pos)
        end
        
        b.ham += b.intMat[b.intIdx, idx] * nnp_ * selfFactor
    end

    # On-site coupling contribution (BulkMonomer only, not FastBulkMonomer)
    if hasfield(typeof(b), :couplingMat) && !isnothing(b.couplingMat)
        # Pass 1: two-way on-site coupling (j == k, j ≠ self)
        for (sp_j, idx_j) in b.intSpecies
            if idx_j == b.intIdx
                continue
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

    return b
end

# Energy increment for the bulkMonomer system
function deltaHam!(b::BulkMonomer)

    # energy increment from self-interactions
    b.deltaH = 0.0

    # Reset partner deltaHPs before accumulating (needed for heterogeneous polymers)
    if !isnothing(b.intMat)
        seen_polys = Set{Polymer}()
        for (sp, idx) in b.intSpecies
            if sp isa PolyMonomer && !(sp.poly in seen_polys)
                sp.poly.deltaHP = 0.0
                push!(seen_polys, sp.poly)
            elseif !(sp isa PolyMonomer)
                sp.deltaHP = 0.0
            end
        end
    end

    # energy increment from cross-interactions with each partner
    if !isnothing(b.intMat)
        for (sp, idx) in b.intSpecies
            energy_change = (b.intMat[b.intIdx, idx] 
                       * ((nnp(b, sp, b.newPos) - nnp(b, sp, b.bi))
                       - Int(b.intIdx == idx))) # subtract 1 in self-interaction
            setDeltaHP(sp, energy_change)
            b.deltaH += energy_change
        end
    end

    # On-site coupling contribution (BulkMonomer only, not FastBulkMonomer)
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

    b.ham += b.deltaH

    return b
end


############# FastBulkMonomer energy functions #############
# MARK: -FastBulkMonomer

# Energy increment for the bulkMonomer system
# Note: FastBulkMonomer does not support three-way couplings.
# Use BulkMonomer if three-way interactions are needed.
function deltaHam!(b::FastBulkMonomer)

    # energy increment from self-interactions
    b.deltaH = 0.0

    if !xor(b.moveIdx == 0, b.moveIdxNew == 0)
        return b
    else
        sign = Int(b.moveIdxNew == 0) * 2 - 1
    end

    # Reset partner deltaHPs before accumulating (needed for heterogeneous polymers)
    if !isnothing(b.intMat)
        seen_polys = Set{Polymer}()
        for (sp, idx) in b.intSpecies
            if sp isa PolyMonomer && !(sp.poly in seen_polys)
                sp.poly.deltaHP = 0.0
                push!(seen_polys, sp.poly)
            elseif !(sp isa PolyMonomer)
                sp.deltaHP = 0.0
            end
        end
    end

    # energy increment from cross-interactions with each partner
    if !isnothing(b.intMat)
        for (sp, idx) in b.intSpecies
            energy_change = (b.intMat[b.intIdx, idx] 
                       * (nnp(b, sp, b.newPos) - nnp(b, sp, b.bi))) * sign
            setDeltaHP(sp, energy_change)
            b.deltaH += energy_change
        end
    end

    b.ham += b.deltaH

    return b
end


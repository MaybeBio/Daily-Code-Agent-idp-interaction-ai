using Droplet
using Test

# Test suite for Polymer Hamiltonian functions (non-interacting case)
# Tests ham! and deltaHam! for polymers without external bulk species

# Helper function to create and initialize a polymer
function create_polymer(length::Int32, maxlength::Int32, boxsize::Int32; 
                       chemPot::Float64=1.0, intPot::Float64=0.5)
    p = Polymer(
        length = length,
        maxlength = maxlength,
        nMonomerTypes = Int32(1)
    )
    
    initMonomers!(p, defaultChemPot=chemPot, defaultIntPot=intPot)
    l = Lattice(boxlength=boxsize, nspecies=Int32(1))
    init!(l, [p])
    
    return p
end

@testset "Polymer Hamiltonian Tests (Non-interacting)" begin
    
    # MARK: - Independent Validation Tests
    @testset "Independent Validation: ham!" begin
        
        @testset "Single bond polymer" begin
            # Create a minimal polymer with 2 beads (1 bond)
            p = create_polymer(Int32(2), Int32(10), Int32(16), chemPot=2.0, intPot=1.5)
            
            # Compute Hamiltonian
            ham!(p)
            
            # Manual calculation:
            # - 2 beads × chemPot = -2 × 2.0 = -4.0
            # - Bead 1: neighbors = 1 (bead 2)
            # - Bead 2: neighbors = 1 (bead 1)
            # - Total neighbors counted: 1 + 1 = 2
            # - With division by 2 to avoid double-counting: 2/2 × 1.5 = 1.5
            # - Total = -4.0 + 1.5 = -2.5
            expected_ham = -2 * 2.0 + (1 + 1) / 2.0 * 1.5
            
            @test p.ham ≈ expected_ham atol=1e-10
        end
        
        @testset "Three bond polymer (line)" begin
            # Create polymer with 4 beads (3 bonds) in a line
            p = create_polymer(Int32(4), Int32(10), Int32(16), chemPot=1.0, intPot=0.5)
            
            ham!(p)
            
            # Manual calculation for linear polymer:
            # - 4 beads × chemPot = -4 × 1.0 = -4.0
            # - Bead 1 (end): 1 neighbor (bead 2)
            # - Bead 2: 2 neighbors (beads 1,3)
            # - Bead 3: 2 neighbors (beads 2,4)
            # - Bead 4 (end): 1 neighbor (bead 3)
            # - Total neighbors: 1 + 2 + 2 + 1 = 6
            # - With division by 2: 6/2 × 0.5 = 1.5
            # - Total = -4.0 + 1.5 = -2.5
            expected_ham = -4 * 1.0 + (1 + 2 + 2 + 1) / 2.0 * 0.5
            
            @test p.ham ≈ expected_ham atol=1e-10
        end
        
        @testset "Polymer with zero potentials" begin
            # Test edge case: all potentials are zero
            p = create_polymer(Int32(3), Int32(10), Int32(16), chemPot=0.0, intPot=0.0)
            
            ham!(p)
            
            # Expected: all terms are zero
            @test p.ham ≈ 0.0 atol=1e-10
        end
        
        @testset "Larger polymer (7 beads)" begin
            # Test with a slightly larger polymer
            p = create_polymer(Int32(7), Int32(10), Int32(16), chemPot=0.8, intPot=0.3)
            
            ham!(p)
            
            # Manual calculation for 7-bead linear polymer:
            # - Chemical: -7 × 0.8 = -5.6
            # - Interaction: 2 ends with 1 neighbor each, 5 middle beads with 2 neighbors each
            # - Total neighbors: 2×1 + 5×2 = 2 + 10 = 12
            # - With division by 2: 12/2 × 0.3 = 1.8
            # - Total = -5.6 + 1.8 = -3.8
            expected_ham = -7 * 0.8 + (2 * 1 + 5 * 2) / 2.0 * 0.3
            
            @test p.ham ≈ expected_ham atol=1e-10
        end
    end
    
    # MARK: - Consistency Tests for Bond Addition
    @testset "Consistency: deltaHam! for addBond!" begin
        
        @testset "Add bond at middle of polymer" begin
            p = create_polymer(Int32(3), Int32(10), Int32(16), chemPot=1.0, intPot=0.5)
            
            # Compute initial Hamiltonian
            ham!(p)
            initial_ham = p.ham
            
            # Select middle bond (index 2)
            p.moveIdx = 2
            identifyBond!(p)
            
            # Manually set up a bond addition move
            p.moveType = 1  # bond addition
            
            # Calculate new positions for the bond
            # Bond vector from bi to bj
            p.newPosi = copy(p.bi)
            p.newPosj = copy(p.bj)
            
            # Move perpendicular to the bond (assuming bond is along x-axis)
            if p.bondVec[1] != 0  # bond along x
                p.newPosi[2] += 1
                p.newPosj[2] += 1
            elseif p.bondVec[2] != 0  # bond along y
                p.newPosi[3] += 1
                p.newPosj[3] += 1
            else  # bond along z
                p.newPosi[1] += 1
                p.newPosj[1] += 1
            end
            
            # Apply periodic boundaries
            boxmod!(p.newPosi, p.newPosi, Int32(1), p.boxlength)
            boxmod!(p.newPosj, p.newPosj, Int32(1), p.boxlength)
            
            # Check that new positions don't overlap
            if !overlap(p)
                # Compute deltaHam
                deltaHam!(p)
                delta_ham = p.deltaH
                
                # Apply the move
                addBond!(p)
                p.length += 2
                
                # Recompute final Hamiltonian
                ham!(p)
                final_ham = p.ham
                
                # Check consistency
                @test final_ham ≈ initial_ham + delta_ham atol=1e-8
            end
        end
        
        @testset "Add bond at start of polymer" begin
            p = create_polymer(Int32(2), Int32(10), Int32(16), chemPot=2.0, intPot=1.0)
            
            ham!(p)
            initial_ham = p.ham
            
            # Select first bond
            p.moveIdx = 1
            identifyBond!(p)
            p.moveType = 1
            
            # Set new positions perpendicular to bond
            p.newPosi = copy(p.bi)
            p.newPosj = copy(p.bj)
            
            if p.bondVec[1] != 0
                p.newPosi[2] += 1
                p.newPosj[2] += 1
            elseif p.bondVec[2] != 0
                p.newPosi[3] += 1
                p.newPosj[3] += 1
            else
                p.newPosi[1] += 1
                p.newPosj[1] += 1
            end
            
            boxmod!(p.newPosi, p.newPosi, Int32(1), p.boxlength)
            boxmod!(p.newPosj, p.newPosj, Int32(1), p.boxlength)
            
            if !overlap(p)
                deltaHam!(p)
                delta_ham = p.deltaH
                
                addBond!(p)
                p.length += 2
                
                ham!(p)
                final_ham = p.ham
                
                @test final_ham ≈ initial_ham + delta_ham atol=1e-8
            end
        end
    end
    
    # MARK: - Consistency Tests for Bond Removal
    @testset "Consistency: deltaHam! for removeBond!" begin
        
        @testset "Remove bond from middle of polymer" begin
            p = create_polymer(Int32(5), Int32(10), Int32(16), chemPot=1.5, intPot=0.8)
            
            ham!(p)
            initial_ham = p.ham
            
            # Select middle bond (index 3)
            p.moveIdx = 3
            identifyBond!(p)
            p.moveType = 3  # bond removal
            
            # For removal, newPosi = bprev and newPosj = bnext
            p.newPosi = copy(p.bprev)
            p.newPosj = copy(p.bnext)
            
            deltaHam!(p)
            delta_ham = p.deltaH
            
            removeBond!(p)
            p.length -= 2
            
            ham!(p)
            final_ham = p.ham
            
            @test final_ham ≈ initial_ham + delta_ham atol=1e-8
        end
        
        @testset "Remove bond from start of polymer" begin
            p = create_polymer(Int32(4), Int32(10), Int32(16), chemPot=1.0, intPot=0.5)
            
            ham!(p)
            initial_ham = p.ham
            
            # Select first bond
            p.moveIdx = 1
            identifyBond!(p)
            p.moveType = 3
            
            p.newPosi = copy(p.bprev)
            p.newPosj = copy(p.bnext)
            
            deltaHam!(p)
            delta_ham = p.deltaH
            
            removeBond!(p)
            p.length -= 2
            
            ham!(p)
            final_ham = p.ham
            
            @test final_ham ≈ initial_ham + delta_ham atol=1e-8
        end
        
        @testset "Remove bond from end of polymer" begin
            p = create_polymer(Int32(4), Int32(10), Int32(16), chemPot=0.5, intPot=0.3)
            
            ham!(p)
            initial_ham = p.ham
            
            # Select second-to-last bond (length-2)
            p.moveIdx = p.length - 2
            identifyBond!(p)
            p.moveType = 3
            
            p.newPosi = copy(p.bprev)
            p.newPosj = copy(p.bnext)
            
            deltaHam!(p)
            delta_ham = p.deltaH
            
            removeBond!(p)
            p.length -= 2
            
            ham!(p)
            final_ham = p.ham
            
            @test final_ham ≈ initial_ham + delta_ham atol=1e-8
        end
    end
    
    # MARK: - Consistency Tests for Kinking
    @testset "Consistency: deltaHam! for kinkBond!" begin
        
        @testset "Kink bond: move ki" begin
            p = create_polymer(Int32(4), Int32(10), Int32(16), chemPot=1.2, intPot=0.6)
            
            ham!(p)
            initial_ham = p.ham
            
            # Select a bond in the middle
            p.moveIdx = 2
            identifyBond!(p)
            p.moveType = 2  # kink
            
            # Kink by moving ki: newPosi = bprev, newPosj = new position
            p.newPosi = copy(p.bprev)
            p.newPosj = copy(p.bj)
            
            # Move newPosj perpendicular to bond
            if p.bondVec[1] != 0
                p.newPosj[2] += 1
            elseif p.bondVec[2] != 0
                p.newPosj[3] += 1
            else
                p.newPosj[1] += 1
            end
            boxmod!(p.newPosj, p.newPosj, Int32(1), p.boxlength)
            
            # Check it's a valid kink and doesn't overlap
            if !overlap(p) && p.newPosi == p.bprev
                deltaHam!(p)
                delta_ham = p.deltaH
                
                kinkBond!(p)
                
                ham!(p)
                final_ham = p.ham
                
                @test final_ham ≈ initial_ham + delta_ham atol=1e-8
            end
        end
        
        @testset "Kink bond: move kj" begin
            p = create_polymer(Int32(5), Int32(10), Int32(16), chemPot=0.9, intPot=0.4)
            
            ham!(p)
            initial_ham = p.ham
            
            # Select a bond in the middle
            p.moveIdx = 2
            identifyBond!(p)
            p.moveType = 2  # kink
            
            # Kink by moving kj: newPosi = new position, newPosj = bnext
            p.newPosi = copy(p.bi)
            p.newPosj = copy(p.bnext)
            
            # Move newPosi perpendicular to bond
            if p.bondVec[1] != 0
                p.newPosi[2] += 1
            elseif p.bondVec[2] != 0
                p.newPosi[3] += 1
            else
                p.newPosi[1] += 1
            end
            boxmod!(p.newPosi, p.newPosi, Int32(1), p.boxlength)
            
            # Check it's a valid kink and doesn't overlap
            if !overlap(p) && p.newPosj == p.bnext
                deltaHam!(p)
                delta_ham = p.deltaH
                
                kinkBond!(p)
                
                ham!(p)
                final_ham = p.ham
                
                @test final_ham ≈ initial_ham + delta_ham atol=1e-8
            end
        end
        
        @testset "Multiple kink operations" begin
            # Test consistency over multiple kink operations
            p = create_polymer(Int32(6), Int32(10), Int32(20), chemPot=1.0, intPot=0.5)
            
            ham!(p)
            cumulative_ham = p.ham
            
            # Try several kink moves
            for trial in 1:3
                p.moveIdx = 2 + (trial - 1)
                if p.moveIdx > p.length - 1
                    break
                end
                
                identifyBond!(p)
                p.moveType = 2
                
                # Alternate between moving ki and kj
                if trial % 2 == 1
                    p.newPosi = copy(p.bprev)
                    p.newPosj = copy(p.bj)
                    p.newPosj[2] += trial
                else
                    p.newPosi = copy(p.bi)
                    p.newPosj = copy(p.bnext)
                    p.newPosi[2] += trial
                end
                
                boxmod!(p.newPosi, p.newPosi, Int32(1), Int32(20))
                boxmod!(p.newPosj, p.newPosj, Int32(1), Int32(20))
                
                if !overlap(p)
                    deltaHam!(p)
                    delta_ham = p.deltaH
                    
                    kinkBond!(p)
                    
                    ham!(p)
                    
                    @test p.ham ≈ cumulative_ham + delta_ham atol=1e-8
                    cumulative_ham = p.ham
                end
            end
        end
    end
    
    # MARK: - Edge Cases
    @testset "Edge Cases" begin
        
        @testset "Minimal polymer (length=2)" begin
            p = create_polymer(Int32(2), Int32(10), Int32(16), chemPot=1.0, intPot=1.0)
            
            ham!(p)
            
            # Should compute without error
            @test p.ham !== nothing
            @test isfinite(p.ham)
        end
        
        @testset "Negative potentials" begin
            p = create_polymer(Int32(3), Int32(10), Int32(16), chemPot=-2.0, intPot=-1.0)
            
            ham!(p)
            
            # Manual calculation:
            # - Chemical: -3 × (-2.0) = 6.0
            # - Interaction: 2 ends with 1 neighbor, 1 middle with 2 neighbors
            # - Total neighbors: 2×1 + 1×2 = 4
            # - With division by 2: 4/2 × (-1.0) = -2.0
            # - Total: 6.0 - 2.0 = 4.0
            expected_ham = 3 * 2.0 + (2 * 1 + 1 * 2) / 2.0 * (-1.0)
            
            @test p.ham ≈ expected_ham atol=1e-10
        end
    end
end

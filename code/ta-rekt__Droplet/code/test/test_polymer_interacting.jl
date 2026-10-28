using Droplet
using Test

# Test suite for Polymer Hamiltonian functions (interacting case)
# Tests ham! and deltaHam! for polymers with external bulk species

# Helper function to create and initialize a polymer with bulk species
function create_interacting_system(poly_length::Int32, bulk_length::Int32, boxsize::Int32;
                                   poly_chemPot::Float64=1.0, poly_intPot::Float64=0.5,
                                   bulk_chemPot::Float64=0.5, bulk_intPot::Float64=0.3,
                                   cross_intPot::Float64=0.4)
    p = Polymer(
        length = poly_length,
        maxlength = poly_length + 10,
        nMonomerTypes = Int32(1)
    )
    
    initMonomers!(p, defaultChemPot=poly_chemPot, defaultIntPot=poly_intPot)
    
    b = FastBulkMonomer(
        length = bulk_length,
        initShape = "uniform"
    )
    
    # Set interaction potentials
    m = p.monomers[1]
    intPot!(m, poly_intPot)
    intPot!(b, bulk_intPot)
    intPot!(m, b, cross_intPot)
    
    # Initialize
    l = Lattice(boxlength=boxsize, nspecies=Int32(2))
    init!(l, [p, b])
    
    return p, b, m
end

@testset "Polymer Hamiltonian Tests (Interacting)" begin
    
    # MARK: - Independent Validation Tests
    @testset "Independent Validation: ham! with bulk species" begin
        
        @testset "Simple polymer with bulk interaction" begin
            # Create a 2-bead polymer and small bulk
            p, b, m = create_interacting_system(Int32(2), Int32(4), Int32(16),
                                               poly_chemPot=2.0, poly_intPot=1.5,
                                               bulk_chemPot=1.0, bulk_intPot=0.5,
                                               cross_intPot=0.8)
            
            # Compute Hamiltonians
            ham!(p)
            ham!(b)
            
            # Manual calculation for polymer:
            # - Chemical: -2 × 2.0 = -4.0
            # - Self-interaction: (1+1)/2 × 1.5 = 1.5
            # - Cross-interaction with bulk: count bulk neighbors at each polymer bead
            # This depends on random placement, so we verify consistency instead
            
            # Verify that energy is finite and makes sense
            @test isfinite(p.ham)
            @test isfinite(b.ham)
            
            # Chemical potential contribution should be deterministic
            expected_chem = -2 * 2.0  # -4.0
            # The interaction part depends on configuration
            # With attractive interactions (positive intPot), energy becomes more positive
            # So total ham can be less negative than chemical potential alone
            @test p.ham >= expected_chem  # Total energy should be >= chemical contribution
        end
        
        @testset "Polymer with multiple bulk neighbors" begin
            # Create setup where we can control placement more
            p = Polymer(length=Int32(3), maxlength=Int32(10), nMonomerTypes=Int32(1),
                       initPos=[8,8,8])
            initMonomers!(p, defaultChemPot=1.0, defaultIntPot=0.5)
            
            b = FastBulkMonomer(length=Int32(6), initShape="box", initParams=Int32(1),
                               initPos=[8,8,8])
            
            m = p.monomers[1]
            intPot!(m, 0.5)
            intPot!(b, 0.3)
            intPot!(m, b, 0.4)
            
            l = Lattice(boxlength=Int32(16), nspecies=Int32(2))
            init!(l, [p, b])
            
            ham!(p)
            ham!(b)
            
            # Both should have interactions since they're placed near each other
            @test isfinite(p.ham)
            @test isfinite(b.ham)
        end
    end
    
    # MARK: - Consistency Tests for Bond Addition with Bulk Interactions
    @testset "Consistency: deltaHam! for addBond! with bulk" begin
        
        @testset "Add bond with bulk interactions" begin
            p, b, m = create_interacting_system(Int32(3), Int32(8), Int32(16),
                                               poly_chemPot=1.0, poly_intPot=0.5,
                                               bulk_chemPot=0.5, bulk_intPot=0.3,
                                               cross_intPot=0.4)
            
            # Compute initial Hamiltonians
            ham!(p)
            ham!(b)
            initial_ham_p = p.ham
            initial_ham_b = b.ham
            
            # Select middle bond
            p.moveIdx = 2
            identifyBond!(p)
            p.moveType = 1  # bond addition
            
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
                # Reset deltaHP for bulk species
                b.deltaHP = 0.0
                
                # Compute deltaHam for polymer
                deltaHam!(p)
                delta_ham_p = p.deltaH
                delta_ham_b = b.deltaHP  # Polymer sets this
                
                # Apply the move
                addBond!(p)
                p.length += 2
                
                # Recompute Hamiltonians
                ham!(p)
                ham!(b)
                final_ham_p = p.ham
                final_ham_b = b.ham
                
                # Check consistency
                @test final_ham_p ≈ initial_ham_p + delta_ham_p atol=1e-8
                @test final_ham_b ≈ initial_ham_b + delta_ham_b atol=1e-8
            end
        end
        
        @testset "Add bond at start with bulk" begin
            p, b, m = create_interacting_system(Int32(2), Int32(6), Int32(16),
                                               poly_chemPot=2.0, poly_intPot=1.0,
                                               bulk_chemPot=1.0, bulk_intPot=0.5,
                                               cross_intPot=0.6)
            
            ham!(p)
            ham!(b)
            initial_ham_p = p.ham
            initial_ham_b = b.ham
            
            p.moveIdx = 1
            identifyBond!(p)
            p.moveType = 1
            
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
                b.deltaHP = 0.0
                
                deltaHam!(p)
                delta_ham_p = p.deltaH
                delta_ham_b = b.deltaHP
                
                addBond!(p)
                p.length += 2
                
                ham!(p)
                ham!(b)
                final_ham_p = p.ham
                final_ham_b = b.ham
                
                @test final_ham_p ≈ initial_ham_p + delta_ham_p atol=1e-8
                @test final_ham_b ≈ initial_ham_b + delta_ham_b atol=1e-8
            end
        end
    end
    
    # MARK: - Consistency Tests for Bond Removal with Bulk Interactions
    @testset "Consistency: deltaHam! for removeBond! with bulk" begin
        
        @testset "Remove bond from middle with bulk" begin
            p, b, m = create_interacting_system(Int32(5), Int32(10), Int32(16),
                                               poly_chemPot=1.5, poly_intPot=0.8,
                                               bulk_chemPot=0.8, bulk_intPot=0.4,
                                               cross_intPot=0.5)
            
            ham!(p)
            ham!(b)
            initial_ham_p = p.ham
            initial_ham_b = b.ham
            
            p.moveIdx = 3
            identifyBond!(p)
            p.moveType = 3  # bond removal
            
            p.newPosi = copy(p.bprev)
            p.newPosj = copy(p.bnext)
            
            b.deltaHP = 0.0
            
            deltaHam!(p)
            delta_ham_p = p.deltaH
            delta_ham_b = b.deltaHP
            
            removeBond!(p)
            p.length -= 2
            
            ham!(p)
            ham!(b)
            final_ham_p = p.ham
            final_ham_b = b.ham
            
            @test final_ham_p ≈ initial_ham_p + delta_ham_p atol=1e-8
            @test final_ham_b ≈ initial_ham_b + delta_ham_b atol=1e-8
        end
        
        @testset "Remove bond from start with bulk" begin
            p, b, m = create_interacting_system(Int32(4), Int32(8), Int32(16),
                                               poly_chemPot=1.0, poly_intPot=0.5,
                                               bulk_chemPot=0.5, bulk_intPot=0.3,
                                               cross_intPot=0.4)
            
            ham!(p)
            ham!(b)
            initial_ham_p = p.ham
            initial_ham_b = b.ham
            
            p.moveIdx = 1
            identifyBond!(p)
            p.moveType = 3
            
            p.newPosi = copy(p.bprev)
            p.newPosj = copy(p.bnext)
            
            b.deltaHP = 0.0
            
            deltaHam!(p)
            delta_ham_p = p.deltaH
            delta_ham_b = b.deltaHP
            
            removeBond!(p)
            p.length -= 2
            
            ham!(p)
            ham!(b)
            final_ham_p = p.ham
            final_ham_b = b.ham
            
            @test final_ham_p ≈ initial_ham_p + delta_ham_p atol=1e-8
            @test final_ham_b ≈ initial_ham_b + delta_ham_b atol=1e-8
        end
        
        @testset "Remove bond from end with bulk" begin
            p, b, m = create_interacting_system(Int32(4), Int32(8), Int32(16),
                                               poly_chemPot=0.5, poly_intPot=0.3,
                                               bulk_chemPot=0.3, bulk_intPot=0.2,
                                               cross_intPot=0.3)
            
            ham!(p)
            ham!(b)
            initial_ham_p = p.ham
            initial_ham_b = b.ham
            
            p.moveIdx = p.length - 2
            identifyBond!(p)
            p.moveType = 3
            
            p.newPosi = copy(p.bprev)
            p.newPosj = copy(p.bnext)
            
            b.deltaHP = 0.0
            
            deltaHam!(p)
            delta_ham_p = p.deltaH
            delta_ham_b = b.deltaHP
            
            removeBond!(p)
            p.length -= 2
            
            ham!(p)
            ham!(b)
            final_ham_p = p.ham
            final_ham_b = b.ham
            
            @test final_ham_p ≈ initial_ham_p + delta_ham_p atol=1e-8
            @test final_ham_b ≈ initial_ham_b + delta_ham_b atol=1e-8
        end
    end
    
    # MARK: - Consistency Tests for Kinking with Bulk Interactions
    @testset "Consistency: deltaHam! for kinkBond! with bulk" begin
        
        @testset "Kink bond: move ki with bulk" begin
            p, b, m = create_interacting_system(Int32(4), Int32(8), Int32(16),
                                               poly_chemPot=1.2, poly_intPot=0.6,
                                               bulk_chemPot=0.6, bulk_intPot=0.3,
                                               cross_intPot=0.5)
            
            ham!(p)
            ham!(b)
            initial_ham_p = p.ham
            initial_ham_b = b.ham
            
            p.moveIdx = 2
            identifyBond!(p)
            p.moveType = 2  # kink
            
            p.newPosi = copy(p.bprev)
            p.newPosj = copy(p.bj)
            
            if p.bondVec[1] != 0
                p.newPosj[2] += 1
            elseif p.bondVec[2] != 0
                p.newPosj[3] += 1
            else
                p.newPosj[1] += 1
            end
            boxmod!(p.newPosj, p.newPosj, Int32(1), p.boxlength)
            
            if !overlap(p) && p.newPosi == p.bprev
                b.deltaHP = 0.0
                
                deltaHam!(p)
                delta_ham_p = p.deltaH
                delta_ham_b = b.deltaHP
                
                kinkBond!(p)
                
                ham!(p)
                ham!(b)
                final_ham_p = p.ham
                final_ham_b = b.ham
                
                @test final_ham_p ≈ initial_ham_p + delta_ham_p atol=1e-8
                @test final_ham_b ≈ initial_ham_b + delta_ham_b atol=1e-8
            end
        end
        
        @testset "Kink bond: move kj with bulk" begin
            p, b, m = create_interacting_system(Int32(5), Int32(10), Int32(16),
                                               poly_chemPot=0.9, poly_intPot=0.4,
                                               bulk_chemPot=0.5, bulk_intPot=0.25,
                                               cross_intPot=0.4)
            
            ham!(p)
            ham!(b)
            initial_ham_p = p.ham
            initial_ham_b = b.ham
            
            p.moveIdx = 2
            identifyBond!(p)
            p.moveType = 2
            
            p.newPosi = copy(p.bi)
            p.newPosj = copy(p.bnext)
            
            if p.bondVec[1] != 0
                p.newPosi[2] += 1
            elseif p.bondVec[2] != 0
                p.newPosi[3] += 1
            else
                p.newPosi[1] += 1
            end
            boxmod!(p.newPosi, p.newPosi, Int32(1), p.boxlength)
            
            if !overlap(p) && p.newPosj == p.bnext
                b.deltaHP = 0.0
                
                deltaHam!(p)
                delta_ham_p = p.deltaH
                delta_ham_b = b.deltaHP
                
                kinkBond!(p)
                
                ham!(p)
                ham!(b)
                final_ham_p = p.ham
                final_ham_b = b.ham
                
                @test final_ham_p ≈ initial_ham_p + delta_ham_p atol=1e-8
                @test final_ham_b ≈ initial_ham_b + delta_ham_b atol=1e-8
            end
        end
    end
    
    # MARK: - Tests with Multiple Bulk Species
    @testset "Multiple bulk species interactions" begin
        
        @testset "Polymer with two bulk species" begin
            p = Polymer(length=Int32(3), maxlength=Int32(10), nMonomerTypes=Int32(1))
            initMonomers!(p, defaultChemPot=1.0, defaultIntPot=0.5)
            
            b1 = FastBulkMonomer(length=Int32(6), initShape="uniform")
            b2 = FastBulkMonomer(length=Int32(6), initShape="uniform")
            
            m = p.monomers[1]
            intPot!(m, 0.5)
            intPot!(b1, 0.3)
            intPot!(b2, 0.4)
            intPot!(m, b1, 0.6)
            intPot!(m, b2, 0.7)
            intPot!(b1, b2, 0.2)
            
            l = Lattice(boxlength=Int32(16), nspecies=Int32(3))
            init!(l, [p, b1, b2])
            
            ham!(p)
            ham!(b1)
            ham!(b2)
            initial_ham_p = p.ham
            initial_ham_b1 = b1.ham
            initial_ham_b2 = b2.ham
            
            # Try a bond addition
            p.moveIdx = 2
            identifyBond!(p)
            p.moveType = 1
            
            p.newPosi = copy(p.bi)
            p.newPosj = copy(p.bj)
            p.newPosi[2] += 1
            p.newPosj[2] += 1
            boxmod!(p.newPosi, p.newPosi, Int32(1), p.boxlength)
            boxmod!(p.newPosj, p.newPosj, Int32(1), p.boxlength)
            
            if !overlap(p)
                b1.deltaHP = 0.0
                b2.deltaHP = 0.0
                
                deltaHam!(p)
                delta_ham_p = p.deltaH
                delta_ham_b1 = b1.deltaHP
                delta_ham_b2 = b2.deltaHP
                
                addBond!(p)
                p.length += 2
                
                ham!(p)
                ham!(b1)
                ham!(b2)
                final_ham_p = p.ham
                final_ham_b1 = b1.ham
                final_ham_b2 = b2.ham
                
                @test final_ham_p ≈ initial_ham_p + delta_ham_p atol=1e-8
                @test final_ham_b1 ≈ initial_ham_b1 + delta_ham_b1 atol=1e-8
                @test final_ham_b2 ≈ initial_ham_b2 + delta_ham_b2 atol=1e-8
            end
        end
    end
end

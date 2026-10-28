using Test
using Droplet

@testset "updateMonomers! re-keys monomers dictionary" begin
    
    @testset "Single polymer with 2 monomer types" begin
        l = Lattice(boxlength=32, nspecies=1)
        
        poly = Polymer(
            length=10,
            maxlength=16,
            nMonomerTypes=2,
            pattern=:block
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        
        # Initially, monomers are keyed by local indices 1, 2
        @test haskey(poly.monomers, 1)
        @test haskey(poly.monomers, 2)
        @test length(poly.monomers) == 2
        
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        initial_base_intIdx = base.intIdx
        initial_pattern_intIdx = pattern.intIdx
        
        # Initially intIdx should match keys (both are 1, 2 in the shared intMat)
        @test haskey(poly.monomers, initial_base_intIdx)
        @test haskey(poly.monomers, initial_pattern_intIdx)
        
        # Set interaction potentials with bulk to trigger potential intIdx changes
        bulk = FastBulkMonomer(length=10, initShape="uniform")
        intPot!(base, 0.5)
        intPot!(pattern, -2.0)
        intPot!(bulk, -1.0)
        intPot!(base, bulk, -0.1)
        intPot!(pattern, bulk, -0.5)
        
        # Old keys might still exist (if intIdx didn't change) or might not
        # The important thing is that updateMonomers! will fix the keys
        
        # Call updateMonomers! to re-key
        updateMonomers!(poly, base.intIdx)
        
        # After updateMonomers!, keys should match intIdx values
        @test haskey(poly.monomers, base.intIdx)
        @test haskey(poly.monomers, pattern.intIdx)
        @test !haskey(poly.monomers, 1) || base.intIdx == 1
        @test !haskey(poly.monomers, 2) || pattern.intIdx == 2
        @test length(poly.monomers) == 2
        
        # baseIntIdx should be updated
        @test poly.baseIntIdx == base.intIdx
        
        # Monomers should still be accessible by their new intIdx
        @test poly.monomers[base.intIdx] === base
        @test poly.monomers[pattern.intIdx] === pattern
    end
    
    @testset "Two polymers with different intIdx values" begin
        l = Lattice(boxlength=32, nspecies=2)
        
        poly1 = Polymer(
            length=10,
            maxlength=16,
            nMonomerTypes=2,
            pattern=:block,
            initPos=[5, 10, 10]
        )
        
        poly2 = Polymer(
            length=10,
            maxlength=16,
            nMonomerTypes=2,
            pattern=:block,
            initPos=[5, 20, 20]
        )
        
        initMonomers!(poly1, defaultChemPot=1.0)
        initMonomers!(poly2, defaultChemPot=1.0)
        
        base1 = poly1.monomers[1]
        pattern1 = poly1.monomers[2]
        base2 = poly2.monomers[1]
        pattern2 = poly2.monomers[2]
        
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        
        # Set interaction potentials - cross-polymer and with bulk to get unique intIdx
        intPot!(base1, 0.5)
        intPot!(base1, pattern1, 0.5)
        intPot!(pattern1, 0.5)
        intPot!(bulk, -1.0)
        intPot!(base1, bulk, -0.1)
        intPot!(pattern1, bulk, -0.2)
        
        intPot!(base2, 0.5)
        intPot!(base2, pattern2, 0.5)
        intPot!(pattern2, 0.5)
        intPot!(base2, bulk, -0.1)
        intPot!(pattern2, bulk, -0.2)
        
        # Cross-polymer interactions to ensure distinct intIdx
        intPot!(base1, base2, 0.0)
        intPot!(base1, pattern2, 0.0)
        intPot!(pattern1, base2, 0.0)
        intPot!(pattern1, pattern2, 0.0)
        
        # All monomers should have different intIdx values now
        intIdxs = [base1.intIdx, pattern1.intIdx, base2.intIdx, pattern2.intIdx]
        @test length(unique(intIdxs)) == 4
        
        # Update both polymers
        updateMonomers!(poly1, base1.intIdx)
        updateMonomers!(poly2, base2.intIdx)
        
        # Check poly1 monomers are keyed correctly
        @test haskey(poly1.monomers, base1.intIdx)
        @test haskey(poly1.monomers, pattern1.intIdx)
        @test poly1.monomers[base1.intIdx] === base1
        @test poly1.monomers[pattern1.intIdx] === pattern1
        @test poly1.baseIntIdx == base1.intIdx
        
        # Check poly2 monomers are keyed correctly
        @test haskey(poly2.monomers, base2.intIdx)
        @test haskey(poly2.monomers, pattern2.intIdx)
        @test poly2.monomers[base2.intIdx] === base2
        @test poly2.monomers[pattern2.intIdx] === pattern2
        @test poly2.baseIntIdx == base2.intIdx
        
        # Monomers should not overlap between polymers
        @test !haskey(poly1.monomers, base2.intIdx)
        @test !haskey(poly1.monomers, pattern2.intIdx)
        @test !haskey(poly2.monomers, base1.intIdx)
        @test !haskey(poly2.monomers, pattern1.intIdx)
    end
    
    @testset "Three monomer types" begin
        l = Lattice(boxlength=32, nspecies=1)
        
        poly = Polymer(
            length=15,
            maxlength=20,
            nMonomerTypes=3,
            pattern=:block
        )
        
        initMonomers!(poly, defaultChemPot=1.0)
        
        mono1 = poly.monomers[1]
        mono2 = poly.monomers[2]
        mono3 = poly.monomers[3]
        
        # Set interaction potentials
        intPot!(mono1, 0.1)
        intPot!(mono2, 0.2)
        intPot!(mono3, 0.3)
        
        # Update monomers
        updateMonomers!(poly, mono1.intIdx)
        
        # All three should be keyed by intIdx
        @test haskey(poly.monomers, mono1.intIdx)
        @test haskey(poly.monomers, mono2.intIdx)
        @test haskey(poly.monomers, mono3.intIdx)
        @test length(poly.monomers) == 3
        
        @test poly.baseIntIdx == mono1.intIdx
    end
    
    @testset "updateMonomers! is idempotent" begin
        l = Lattice(boxlength=32, nspecies=1)
        
        poly = Polymer(
            length=10,
            maxlength=16,
            nMonomerTypes=2,
            pattern=:block
        )
        
        initMonomers!(poly, defaultChemPot=1.0)
        
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.5)
        intPot!(pattern, -2.0)
        
        # First update
        updateMonomers!(poly, base.intIdx)
        
        keys_after_first = sort(collect(keys(poly.monomers)))
        baseIntIdx_after_first = poly.baseIntIdx
        
        # Second update (should be idempotent)
        updateMonomers!(poly, base.intIdx)
        
        keys_after_second = sort(collect(keys(poly.monomers)))
        baseIntIdx_after_second = poly.baseIntIdx
        
        @test keys_after_first == keys_after_second
        @test baseIntIdx_after_first == baseIntIdx_after_second
    end
    
    @testset "Monomers with bulk interactions" begin
        l = Lattice(boxlength=32, nspecies=2)
        
        poly = Polymer(
            length=10,
            maxlength=16,
            nMonomerTypes=2,
            pattern=:block
        )
        
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        
        initMonomers!(poly, defaultChemPot=1.0)
        
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        # Set cross-species interactions
        intPot!(base, 0.5)
        intPot!(pattern, 0.5)
        intPot!(bulk, -1.0)
        intPot!(base, bulk, -0.1)
        intPot!(pattern, bulk, -0.5)
        
        # Update polymer monomers
        updateMonomers!(poly, base.intIdx)
        
        # Check that intIdx is properly updated
        @test haskey(poly.monomers, base.intIdx)
        @test haskey(poly.monomers, pattern.intIdx)
        @test poly.baseIntIdx == base.intIdx
        
        # Verify that interaction matrix is still accessible
        @test base.intMat !== nothing
        @test pattern.intMat !== nothing
        @test base.intMat === pattern.intMat  # Should share the same matrix
    end
    
    @testset "baseIntIdx parameter is respected" begin
        l = Lattice(boxlength=32, nspecies=1)
        
        poly = Polymer(
            length=10,
            maxlength=16,
            nMonomerTypes=2,
            pattern=:block
        )
        
        initMonomers!(poly, defaultChemPot=1.0)
        
        mono1 = poly.monomers[1]
        mono2 = poly.monomers[2]
        
        intPot!(mono1, 0.5)
        intPot!(mono2, -2.0)
        
        # Set mono2 as base (non-standard choice)
        updateMonomers!(poly, mono2.intIdx)
        
        @test poly.baseIntIdx == mono2.intIdx
        
        # Set mono1 as base (standard choice)
        updateMonomers!(poly, mono1.intIdx)
        
        @test poly.baseIntIdx == mono1.intIdx
    end
end

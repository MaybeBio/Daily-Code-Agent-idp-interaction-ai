using Test
using Droplet

@testset "growPolymer! grows polymer to maxlength" begin
    
    @testset "Basic growth with single monomer type" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=10,
            maxlength=30,
            nMonomerTypes=1,
            pattern=:uniform,
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        intPot!(poly.monomers[1], 0.0)
        init!(l, [poly])
        
        @test poly.length == 10
        
        growPolymer!(poly)
        
        @test poly.length == poly.maxlength
        @test poly.length == 30
    end
    
    @testset "Growth with two monomer types" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=8,
            maxlength=24,
            nMonomerTypes=2,
            pattern=:uniform,  # Start uniform, assign pattern after growth
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        @test poly.length == 8
        
        growPolymer!(poly)
        
        @test poly.length == poly.maxlength
        @test poly.length == 24
    end
    
    @testset "Chemical potentials are restored after growth" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=10,
            maxlength=20,
            nMonomerTypes=2,
            pattern=:uniform,
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.5, defaultIntPot=0.0)
        
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(pattern, 0.0)
        
        # Set custom chemical potentials
        base.chemPot = 2.0
        pattern.chemPot = 3.0
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        growPolymer!(poly)
        
        # Check that chemical potentials are restored
        @test base.chemPot == 2.0
        @test pattern.chemPot == 3.0
    end
    
    @testset "Already at maxlength returns immediately" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=20,
            maxlength=20,
            nMonomerTypes=1,
            pattern=:uniform,
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        intPot!(poly.monomers[1], 0.0)
        init!(l, [poly])
        
        @test poly.length == poly.maxlength
        
        # Should return immediately without error
        growPolymer!(poly)
        
        @test poly.length == poly.maxlength
    end
    
    @testset "Growth with bulk species present" begin
        l = Lattice(boxlength=64, nspecies=2)
        
        # Note: length and maxlength must have same parity (both even or both odd)
        # since addBond adds 2 beads at a time
        poly = Polymer(
            length=10,
            maxlength=24,
            nMonomerTypes=2,
            pattern=:uniform,
            initShape="line"
        )
        
        bulk = FastBulkMonomer(length=100, initShape="uniform")
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(pattern, 0.0)
        intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        intPot!(pattern, bulk, -0.2)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])
        
        growPolymer!(poly)
        
        @test poly.length == poly.maxlength
    end
end

@testset "updatePattern! updates polymer pattern correctly" begin
    
    @testset "Basic pattern update" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=20,
            maxlength=20,
            nMonomerTypes=2,
            pattern=:uniform,
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Create new pattern: 8 base, 4 pattern, 8 base
        newPattern = vcat(
            fill(base.intIdx, 8),
            fill(pattern.intIdx, 4),
            fill(base.intIdx, 8)
        )
        
        updatePattern!(l, poly, newPattern)
        
        # Verify pattern was applied
        @test poly.customPattern == newPattern
        
        # Check pos structure
        key = poly.head
        position = 1
        while key != 0
            @test poly.pos[key][6] == newPattern[position]
            key = poly.pos[key][5]
            position += 1
        end
    end
    
    @testset "Lattice is updated correctly" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=15,
            maxlength=15,
            nMonomerTypes=2,
            pattern=:uniform,
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Create new pattern: 5 base, 5 pattern, 5 base
        newPattern = vcat(
            fill(base.intIdx, 5),
            fill(pattern.intIdx, 5),
            fill(base.intIdx, 5)
        )
        
        updatePattern!(l, poly, newPattern)
        
        # Check lattice at each bead position
        key = poly.head
        position = 1
        while key != 0
            x, y, z = poly.pos[key][1], poly.pos[key][2], poly.pos[key][3]
            @test l.lat[1, x, y, z] == newPattern[position]
            key = poly.pos[key][5]
            position += 1
        end
    end
    
    @testset "patternTails is recalculated" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=20,
            maxlength=20,
            nMonomerTypes=2,
            pattern=:uniform,  # Initially no patternTails
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Initially uniform, no patternTails
        @test poly.patternTails === nothing || isempty(poly.patternTails)
        
        # Create pattern with non-base region
        newPattern = vcat(
            fill(base.intIdx, 6),
            fill(pattern.intIdx, 8),
            fill(base.intIdx, 6)
        )
        
        updatePattern!(l, poly, newPattern)
        
        # patternTails should now exist
        @test poly.patternTails !== nothing
        @test haskey(poly.patternTails, pattern.intIdx)
        
        # The tail should be at the end of the pattern block
        tailKey = poly.patternTails[pattern.intIdx]
        @test poly.pos[tailKey][6] == pattern.intIdx
    end
    
    @testset "Uses customPattern when no newPattern provided" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=16,
            maxlength=16,
            nMonomerTypes=2,
            pattern=:uniform,
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Set customPattern before calling updatePattern!
        poly.customPattern = vcat(
            fill(base.intIdx, 4),
            fill(pattern.intIdx, 8),
            fill(base.intIdx, 4)
        )
        
        # Call without newPattern - should use customPattern
        updatePattern!(l, poly)
        
        # Verify pattern was applied from customPattern
        key = poly.head
        position = 1
        while key != 0
            @test poly.pos[key][6] == poly.customPattern[position]
            key = poly.pos[key][5]
            position += 1
        end
    end
    
    @testset "Error when polymer not at maxlength" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=10,
            maxlength=20,
            nMonomerTypes=2,
            pattern=:uniform,
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        newPattern = fill(base.intIdx, 10)  # Matches current length but not maxlength
        
        @test_throws ErrorException updatePattern!(l, poly, newPattern)
    end
    
    @testset "Error when pattern length mismatch" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=20,
            maxlength=20,
            nMonomerTypes=2,
            pattern=:uniform,
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        
        base = poly.monomers[1]
        
        intPot!(base, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Pattern too short
        shortPattern = fill(base.intIdx, 15)
        @test_throws ErrorException updatePattern!(l, poly, shortPattern)
        
        # Pattern too long
        longPattern = fill(base.intIdx, 25)
        @test_throws ErrorException updatePattern!(l, poly, longPattern)
    end
    
    @testset "Error when pattern contains invalid intIdx" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=20,
            maxlength=20,
            nMonomerTypes=2,
            pattern=:uniform,
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        
        base = poly.monomers[1]
        
        intPot!(base, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Pattern with non-existent intIdx
        invalidPattern = fill(Int32(999), 20)
        @test_throws ErrorException updatePattern!(l, poly, invalidPattern)
    end
    
    @testset "Error when no pattern available" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=20,
            maxlength=20,
            nMonomerTypes=1,
            pattern=:uniform,
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        intPot!(poly.monomers[1], 0.0)
        init!(l, [poly])
        
        # Clear customPattern
        poly.customPattern = nothing
        
        # Call without newPattern when customPattern is nothing
        @test_throws ErrorException updatePattern!(l, poly)
    end
end

@testset "growPolymer! and updatePattern! integration" begin
    
    @testset "Full workflow: grow then update pattern" begin
        l = Lattice(boxlength=64, nspecies=1)
        
        poly = Polymer(
            length=10,
            maxlength=30,
            nMonomerTypes=2,
            pattern=:uniform,
            initShape="line"
        )
        
        initMonomers!(poly, defaultChemPot=1.0, defaultIntPot=0.0)
        
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Step 1: Grow
        growPolymer!(poly)
        @test poly.length == poly.maxlength
        
        # Step 2: Define and apply pattern
        newPattern = vcat(
            fill(base.intIdx, 10),
            fill(pattern.intIdx, 10),
            fill(base.intIdx, 10)
        )
        
        updatePattern!(l, poly, newPattern)
        
        # Verify pattern
        key = poly.head
        position = 1
        while key != 0
            @test poly.pos[key][6] == newPattern[position]
            key = poly.pos[key][5]
            position += 1
        end
        
        # Verify patternTails
        @test haskey(poly.patternTails, pattern.intIdx)
    end
    
    @testset "Two polymers grow and pattern sequentially" begin
        l = Lattice(boxlength=64, nspecies=2)
        
        poly1 = Polymer(
            length=8,
            maxlength=24,
            nMonomerTypes=2,
            pattern=:uniform,
            initShape="line",
            initPos=[5, 16, 16]
        )
        
        poly2 = Polymer(
            length=8,
            maxlength=24,
            nMonomerTypes=2,
            pattern=:uniform,
            initShape="line",
            initPos=[5, 48, 48]
        )
        
        initMonomers!(poly1, defaultChemPot=1.0, defaultIntPot=0.0)
        initMonomers!(poly2, defaultChemPot=1.0, defaultIntPot=0.0)
        
        base1 = poly1.monomers[1]
        pattern1 = poly1.monomers[2]
        base2 = poly2.monomers[1]
        pattern2 = poly2.monomers[2]
        
        # Set interactions to get unique intIdx
        intPot!(base1, 0.0)
        intPot!(pattern1, 0.0)
        intPot!(base2, 0.0)
        intPot!(pattern2, 0.0)
        intPot!(base1, base2, 0.0)
        
        updateMonomers!(poly1, base1.intIdx)
        updateMonomers!(poly2, base2.intIdx)
        
        init!(l, [poly1, poly2])
        
        # Grow both polymers sequentially
        growPolymer!(poly1)
        growPolymer!(poly2)
        
        @test poly1.length == poly1.maxlength
        @test poly2.length == poly2.maxlength
        
        # Apply patterns
        pattern1Vec = vcat(
            fill(base1.intIdx, 8),
            fill(pattern1.intIdx, 8),
            fill(base1.intIdx, 8)
        )
        
        pattern2Vec = vcat(
            fill(base2.intIdx, 8),
            fill(pattern2.intIdx, 8),
            fill(base2.intIdx, 8)
        )
        
        updatePattern!(l, poly1, pattern1Vec)
        updatePattern!(l, poly2, pattern2Vec)
        
        # Verify both patterns applied correctly
        @test haskey(poly1.patternTails, pattern1.intIdx)
        @test haskey(poly2.patternTails, pattern2.intIdx)
    end
end

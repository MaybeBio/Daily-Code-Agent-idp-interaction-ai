using Test
using Droplet

@testset "centerMassPattern tests" begin
    
    @testset "Single bead pattern in middle" begin
        # Simple case: 5-bead polymer with 1 pattern bead in the middle
        l = Lattice(boxlength=20, nspecies=2)
        
        poly = Polymer(
            length=5,
            maxlength=5,
            initShape="line",
            nMonomerTypes=2,
            initPos=[10, 10, 10],
            pattern=:uniform
        )
        
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(base, pattern, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Pattern: base-base-pattern-base-base
        patternVec = Int32[base.intIdx, base.intIdx, pattern.intIdx, base.intIdx, base.intIdx]
        updatePattern!(l, poly, patternVec)
        
        # Get position of bead 3 (the pattern bead)
        key = poly.head
        key = poly.pos[key][5]  # bead 2
        key = poly.pos[key][5]  # bead 3
        expected_cm = Float64.(poly.pos[key][1:3])
        
        cm = centerMassPattern(poly, pattern.intIdx)
        
        @test cm ≈ expected_cm
    end
    
    @testset "Pattern at head" begin
        # Pattern beads at the beginning of polymer
        l = Lattice(boxlength=20, nspecies=2)
        
        poly = Polymer(
            length=6,
            maxlength=6,
            initShape="line",
            nMonomerTypes=2,
            initPos=[10, 10, 10],
            pattern=:uniform
        )
        
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(base, pattern, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Pattern: pattern-pattern-pattern-base-base-base
        patternVec = Int32[pattern.intIdx, pattern.intIdx, pattern.intIdx, 
                           base.intIdx, base.intIdx, base.intIdx]
        updatePattern!(l, poly, patternVec)
        
        # Collect positions of first 3 beads
        key = poly.head
        positions = []
        for i in 1:3
            push!(positions, Float64.(poly.pos[key][1:3]))
            key = poly.pos[key][5]
        end
        
        expected_cm = sum(positions) ./ 3
        cm = centerMassPattern(poly, pattern.intIdx)
        
        @test cm ≈ expected_cm
    end
    
    @testset "Pattern at tail" begin
        # Pattern beads at the end of polymer
        l = Lattice(boxlength=20, nspecies=2)
        
        poly = Polymer(
            length=6,
            maxlength=6,
            initShape="line",
            nMonomerTypes=2,
            initPos=[10, 10, 10],
            pattern=:uniform
        )
        
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(base, pattern, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Pattern: base-base-base-pattern-pattern-pattern
        patternVec = Int32[base.intIdx, base.intIdx, base.intIdx,
                           pattern.intIdx, pattern.intIdx, pattern.intIdx]
        updatePattern!(l, poly, patternVec)
        
        # Collect positions of last 3 beads (beads 4-6)
        key = poly.head
        positions = []
        for i in 1:6
            if i >= 4
                push!(positions, Float64.(poly.pos[key][1:3]))
            end
            key = poly.pos[key][5]
        end
        
        expected_cm = sum(positions) ./ 3
        cm = centerMassPattern(poly, pattern.intIdx)
        
        @test cm ≈ expected_cm
    end
    
    @testset "All beads are pattern" begin
        # Entire polymer is one pattern
        l = Lattice(boxlength=20, nspecies=2)
        
        poly = Polymer(
            length=4,
            maxlength=4,
            initShape="line",
            nMonomerTypes=2,
            initPos=[10, 10, 10],
            pattern=:uniform
        )
        
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(base, pattern, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Pattern: all pattern beads
        patternVec = Int32[pattern.intIdx, pattern.intIdx, pattern.intIdx, pattern.intIdx]
        updatePattern!(l, poly, patternVec)
        
        # Collect all positions
        key = poly.head
        positions = []
        for i in 1:4
            push!(positions, Float64.(poly.pos[key][1:3]))
            key = poly.pos[key][5]
        end
        
        expected_cm = sum(positions) ./ 4
        cm = centerMassPattern(poly, pattern.intIdx)
        
        @test cm ≈ expected_cm
    end
    
    @testset "No pattern exists" begin
        # All beads are base type, pattern type doesn't exist
        l = Lattice(boxlength=20, nspecies=2)
        
        poly = Polymer(
            length=5,
            maxlength=5,
            initShape="line",
            nMonomerTypes=2,
            initPos=[10, 10, 10],
            pattern=:uniform
        )
        
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(base, pattern, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # All base beads, no pattern
        patternVec = Int32[base.intIdx, base.intIdx, base.intIdx, base.intIdx, base.intIdx]
        updatePattern!(l, poly, patternVec)
        
        cm = centerMassPattern(poly, pattern.intIdx)
        
        # Should return zeros when pattern doesn't exist
        @test cm == [0.0, 0.0, 0.0]
    end
    
    @testset "Two separate pattern regions" begin
        # Two separate pattern regions should only count the tail region
        l = Lattice(boxlength=20, nspecies=2)
        
        poly = Polymer(
            length=9,
            maxlength=9,
            initShape="line",
            nMonomerTypes=2,
            initPos=[10, 10, 10],
            pattern=:uniform
        )
        
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(base, pattern, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Pattern: pattern-pattern-base-base-base-pattern-pattern-pattern-base
        patternVec = Int32[pattern.intIdx, pattern.intIdx, base.intIdx, 
                           base.intIdx, base.intIdx,
                           pattern.intIdx, pattern.intIdx, pattern.intIdx, base.intIdx]
        updatePattern!(l, poly, patternVec)
        
        # Function should find the LAST contiguous pattern region (beads 6-8)
        # because it walks backward from the tail
        key = poly.head
        positions = []
        for i in 1:9
            if i >= 6 && i <= 8
                push!(positions, Float64.(poly.pos[key][1:3]))
            end
            key = poly.pos[key][5]
        end
        
        expected_cm = sum(positions) ./ 3
        cm = centerMassPattern(poly, pattern.intIdx)
        
        @test cm ≈ expected_cm
    end
    
    @testset "Very short polymer (2 beads)" begin
        # Edge case: minimal polymer with pattern
        l = Lattice(boxlength=20, nspecies=2)
        
        poly = Polymer(
            length=2,
            maxlength=2,
            initShape="line",
            nMonomerTypes=2,
            initPos=[10, 10, 10],
            pattern=:uniform
        )
        
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        pattern = poly.monomers[2]
        
        intPot!(base, 0.0)
        intPot!(base, pattern, 0.0)
        intPot!(pattern, 0.0)
        
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])
        
        # Pattern: base-pattern
        patternVec = Int32[base.intIdx, pattern.intIdx]
        updatePattern!(l, poly, patternVec)
        
        # Get position of bead 2
        key = poly.head
        key = poly.pos[key][5]
        expected_cm = Float64.(poly.pos[key][1:3])
        
        cm = centerMassPattern(poly, pattern.intIdx)
        
        @test cm ≈ expected_cm
    end
end

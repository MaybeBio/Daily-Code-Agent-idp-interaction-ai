using Test
using Droplet

@testset "find_patterns function" begin
    boxlength = 128
    
    @testset "Uniform pattern (single monomer type)" begin
        # Create polymer with uniform pattern (all beads type 1)
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        poly = Polymer(;
            length = 20,
            maxlength = 32,
            nMonomerTypes = 1,
            pattern = :uniform,
            initShape = "line"
        )
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        result = find_patterns(poly)
        
        # Should have exactly one intIdx
        @test length(result) == 1
        
        # Get the intIdx (should be 1)
        intIdx = first(keys(result))
        patterns = result[intIdx]
        
        # Should have exactly one pattern spanning the whole polymer
        @test length(patterns) == 1
        @test patterns[1] == (1, 20)
    end
    
    @testset "Custom pattern: [1,1,2,2,2,1,3,3,1,1]" begin
        # Create polymer with custom pattern
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        poly = Polymer(;
            length = 10,
            maxlength = 16,
            nMonomerTypes = 3,
            pattern = :custom,
            customPattern = [1,1,2,2,2,1,3,3,1,1],
            initShape = "line"
        )
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        result = find_patterns(poly)
        
        # Should have 3 monomer types
        @test length(result) == 3
        
        # Get intIdx values (should be sorted)
        intIdxList = sort(collect(keys(result)))
        intIdx1, intIdx2, intIdx3 = intIdxList
        
        # Type 1 should appear in 3 regions: (1,2), (6,6), (9,10)
        @test length(result[intIdx1]) == 3
        @test (1, 2) in result[intIdx1]
        @test (6, 6) in result[intIdx1]
        @test (9, 10) in result[intIdx1]
        
        # Type 2 should appear in 1 region: (3,5)
        @test length(result[intIdx2]) == 1
        @test result[intIdx2][1] == (3, 5)
        
        # Type 3 should appear in 1 region: (7,8)
        @test length(result[intIdx3]) == 1
        @test result[intIdx3][1] == (7, 8)
    end
    
    @testset "Alternating pattern" begin
        # Note: The codebase requires non-base types (2+) to appear in contiguous blocks
        # So we'll test a pattern where type 1 is dispersed but type 2 is contiguous
        # Pattern: [1,1,2,2,2,2,1,1,1,1]
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        poly = Polymer(;
            length = 10,
            maxlength = 16,
            nMonomerTypes = 2,
            pattern = :custom,
            customPattern = [1,1,2,2,2,2,1,1,1,1],
            initShape = "line"
        )
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        result = find_patterns(poly)
        
        # Should have 2 monomer types
        @test length(result) == 2
        
        intIdxList = sort(collect(keys(result)))
        intIdx1, intIdx2 = intIdxList
        
        # Type 1 appears in 2 regions: (1,2) and (7,10)
        @test length(result[intIdx1]) == 2
        @test (1, 2) in result[intIdx1]
        @test (7, 10) in result[intIdx1]
        
        # Type 2 appears once: (3,6)
        @test length(result[intIdx2]) == 1
        @test result[intIdx2][1] == (3, 6)
    end
    
    @testset "Block pattern (2 types)" begin
        # Create polymer with block pattern (divides into 2 equal blocks)
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        poly = Polymer(;
            length = 20,
            maxlength = 32,
            nMonomerTypes = 2,
            pattern = :block,
            initShape = "line"
        )
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        result = find_patterns(poly)
        
        # Should have 2 monomer types
        @test length(result) == 2
        
        intIdxList = sort(collect(keys(result)))
        intIdx1, intIdx2 = intIdxList
        
        # Each type should appear once in a contiguous block
        @test length(result[intIdx1]) == 1
        @test length(result[intIdx2]) == 1
        
        # Type 1: first half (1-10)
        @test result[intIdx1][1] == (1, 10)
        
        # Type 2: second half (11-20)
        @test result[intIdx2][1] == (11, 20)
    end
    
    @testset "Block pattern (3 types)" begin
        # Create polymer with block pattern (divides into 3 equal blocks)
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        poly = Polymer(;
            length = 30,
            maxlength = 48,
            nMonomerTypes = 3,
            pattern = :block,
            initShape = "line"
        )
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        result = find_patterns(poly)
        
        # Should have 3 monomer types
        @test length(result) == 3
        
        intIdxList = sort(collect(keys(result)))
        intIdx1, intIdx2, intIdx3 = intIdxList
        
        # Each type should appear once in a contiguous block of 10 beads
        @test length(result[intIdx1]) == 1
        @test length(result[intIdx2]) == 1
        @test length(result[intIdx3]) == 1
        
        @test result[intIdx1][1] == (1, 10)
        @test result[intIdx2][1] == (11, 20)
        @test result[intIdx3][1] == (21, 30)
    end
    
    @testset "Single bead polymer" begin
        # Edge case: polymer with only 1 bead
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        poly = Polymer(;
            length = 1,
            maxlength = 4,
            nMonomerTypes = 1,
            pattern = :uniform,
            initShape = "line"
        )
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        result = find_patterns(poly)
        
        @test length(result) == 1
        intIdx = first(keys(result))
        @test length(result[intIdx]) == 1
        @test result[intIdx][1] == (1, 1)
    end
    
    @testset "Contiguous blocks with custom pattern" begin
        # Test a more complex custom pattern with multiple blocks
        # Pattern: [1,1,1,2,2,1,1,1,1,3,3,3,1,1]
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        poly = Polymer(;
            length = 14,
            maxlength = 20,
            nMonomerTypes = 3,
            pattern = :custom,
            customPattern = [1,1,1,2,2,1,1,1,1,3,3,3,1,1],
            initShape = "line"
        )
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        result = find_patterns(poly)
        
        @test length(result) == 3
        
        intIdxList = sort(collect(keys(result)))
        intIdx1, intIdx2, intIdx3 = intIdxList
        
        # Type 1 appears in 3 blocks: (1,3), (6,9), (13,14)
        @test length(result[intIdx1]) == 3
        @test (1, 3) in result[intIdx1]
        @test (6, 9) in result[intIdx1]
        @test (13, 14) in result[intIdx1]
        
        # Type 2 appears once: (4,5)
        @test length(result[intIdx2]) == 1
        @test result[intIdx2][1] == (4, 5)
        
        # Type 3 appears once: (10,12)
        @test length(result[intIdx3]) == 1
        @test result[intIdx3][1] == (10, 12)
    end
    
    @testset "Pattern coverage verification" begin
        # Verify that all beads are covered exactly once
        # Pattern respects contiguity constraint: types 2+ appear in single contiguous blocks
        # Pattern: [1,1,2,2,2,2,1,1,1,3,3,3,3,1,1,1,1,1,1,1,4,4,4,4,1]
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        poly = Polymer(;
            length = 25,
            maxlength = 32,
            nMonomerTypes = 4,
            pattern = :custom,
            customPattern = [1,1,2,2,2,2,1,1,1,3,3,3,3,1,1,1,1,1,1,1,4,4,4,4,1],
            initShape = "line"
        )
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        result = find_patterns(poly)
        
        # Collect all covered positions
        covered = Set{Int}()
        for (intIdx, patterns) in result
            for (start, finish) in patterns
                for pos in start:finish
                    # Each position should be covered exactly once
                    @test !(pos in covered)
                    push!(covered, pos)
                end
            end
        end
        
        # All positions from 1 to length should be covered
        @test length(covered) == poly.length
        @test covered == Set(1:poly.length)
    end
end

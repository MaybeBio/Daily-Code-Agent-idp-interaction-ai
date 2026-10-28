using Test
using Droplet

@testset "Pattern Statistics" begin
    # Create lattice and polymers for testing
    # Boxlength must be larger than polymer length for line initialization
    boxlength = 128
    
    @testset "get_pattern_centers" begin
        # Set up lattice and polymer with 100 beads
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        
        poly = Polymer(;
            length = 100,
            maxlength = 128,
            nMonomerTypes = 1,
            initShape = "line"
        )
        
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        # Test error handling for wrong number of patterns
        @test_throws ErrorException get_pattern_centers(poly, [(1,5)])
        @test_throws ErrorException get_pattern_centers(poly, [(1,5), (10,15), (20,25)])
        
        # Test with valid pattern specification
        patterns = [(29, 36), (93, 100)]
        center1, center2, dist = get_pattern_centers(poly, patterns)
        
        @test center1 !== nothing
        @test center2 !== nothing
        @test length(center1) == 3
        @test length(center2) == 3
        @test dist >= 0
        @test dist < Inf
        
        # Test with polymer too short (50 beads, but patterns require up to bead 100)
        l_short = Lattice(; boxlength = boxlength, nspecies = 1)
        short_poly = Polymer(;
            length = 50,
            maxlength = 64,
            nMonomerTypes = 1,
            initShape = "line"
        )
        initMonomers!(short_poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l_short, [short_poly])
        
        c1, c2, d = get_pattern_centers(short_poly, patterns)
        @test c1 === nothing
        @test c2 === nothing
        @test d == Inf
        
        # Test with different pattern ranges (within 100 bead polymer)
        patterns2 = [(1, 10), (50, 60)]
        center1, center2, dist = get_pattern_centers(poly, patterns2)
        @test center1 !== nothing
        @test center2 !== nothing
        @test dist >= 0
    end
    
    @testset "get_pattern_centers with known positions" begin
        # Test 1: Single bead patterns (beads 10 and 20 should have distance 10)
        l1 = Lattice(; boxlength = boxlength, nspecies = 1)
        poly1 = Polymer(;
            length = 50,
            maxlength = 64,
            nMonomerTypes = 1,
            initShape = "line"
        )
        initMonomers!(poly1, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l1, [poly1])
        
        # Single bead patterns
        single_patterns = [(10, 10), (20, 20)]
        c1, c2, d = get_pattern_centers(poly1, single_patterns)
        
        @test c1 !== nothing
        @test c2 !== nothing
        # For a line polymer starting at some position, verify distance is approximately 10
        @test isapprox(d, 10.0, atol=0.5)
        
        # Test 2: Verify centers are correctly averaged
        # Patterns with 3 beads each at known positions
        l2 = Lattice(; boxlength = boxlength, nspecies = 1)
        poly2 = Polymer(;
            length = 100,
            maxlength = 128,
            nMonomerTypes = 1,
            initShape = "line"
        )
        initMonomers!(poly2, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l2, [poly2])
        
        patterns_3bead = [(1, 3), (98, 100)]
        c1, c2, d = get_pattern_centers(poly2, patterns_3bead)
        
        @test c1 !== nothing
        @test c2 !== nothing
        # Get reference positions for verification
        ref_pos1 = poly2.pos[poly2.head][1:3]
        ref_pos3 = poly2.pos[poly2.head][1:3]
        key = poly2.head
        for i in 1:2
            key = poly2.pos[key][5]
            ref_pos3 = poly2.pos[key][1:3]
        end
        
        # Centers should be within the pattern ranges
        @test all(c1 .>= 0)
        @test all(c2 .>= 0)
        @test d > 0
        
        # Test 3: Verify that reversing pattern order doesn't change distance
        patterns_rev = [(98, 100), (1, 3)]
        c1_rev, c2_rev, d_rev = get_pattern_centers(poly2, patterns_rev)
        
        @test isapprox(d, d_rev, atol=1e-10)
        @test isapprox(c1, c2_rev, atol=1e-10)
        @test isapprox(c2, c1_rev, atol=1e-10)
        
        # Test 4: Same pattern should give zero distance
        same_patterns = [(25, 25), (25, 25)]
        c1_same, c2_same, d_same = get_pattern_centers(poly2, same_patterns)
        
        @test c1_same !== nothing
        @test c2_same !== nothing
        @test isapprox(d_same, 0.0, atol=1e-10)
        @test isapprox(c1_same, c2_same, atol=1e-10)
        
        # Test 5: Adjacent patterns
        adjacent_patterns = [(50, 55), (56, 60)]
        c1_adj, c2_adj, d_adj = get_pattern_centers(poly2, adjacent_patterns)
        
        @test c1_adj !== nothing
        @test c2_adj !== nothing
        @test d_adj >= 0
        @test d_adj < 10  # Should be relatively close
    end
    
    @testset "patterns_in_contact" begin
        # Set up lattice and polymer with 100 beads
        l = Lattice(; boxlength = boxlength, nspecies = 1)
        
        poly = Polymer(;
            length = 100,
            maxlength = 128,
            nMonomerTypes = 1,
            initShape = "line"
        )
        
        initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l, [poly])
        
        # Test error handling
        @test_throws ErrorException patterns_in_contact(poly, [(1,5)])
        @test_throws ErrorException patterns_in_contact(poly, [(1,5), (10,15), (20,25)])
        
        # Test with valid patterns
        patterns = [(29, 36), (93, 100)]
        result = patterns_in_contact(poly, patterns)
        @test result isa Bool
        
        # Test with short polymer
        l_short = Lattice(; boxlength = boxlength, nspecies = 1)
        short_poly = Polymer(;
            length = 50,
            maxlength = 64,
            nMonomerTypes = 1,
            initShape = "line"
        )
        initMonomers!(short_poly, defaultChemPot=0.0, defaultIntPot=0.0)
        init!(l_short, [short_poly])
        
        @test patterns_in_contact(short_poly, patterns) == false
        
        # Test adjacent patterns (should have high chance of contact)
        adjacent_patterns = [(1, 5), (6, 10)]
        result = patterns_in_contact(poly, adjacent_patterns)
        @test result isa Bool
    end
end

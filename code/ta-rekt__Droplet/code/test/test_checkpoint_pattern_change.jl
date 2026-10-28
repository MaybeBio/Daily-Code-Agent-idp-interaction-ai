using Test
using Droplet

@testset "Checkpoint pattern change functionality" begin

    @testset "Pattern change: uniform → custom (heterogeneous)" begin
        # Start with a uniform pattern, checkpoint, restore with custom pattern
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1] 
        pat = poly.monomers[2]
        intPot!(base, -0.5) 
        intPot!(pat, -0.3) 
        intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        # Run some moves
        for _ in 1:3000 
            move!(poly) 
            move!(bulk) 
        end

        # Save positions before checkpoint
        saved_positions = [(poly.pos[k][1], poly.pos[k][2], poly.pos[k][3]) for k in poly.keys_vec[1:poly.length]]

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        # Restore with a CUSTOM pattern (different from uniform)
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        
        # Create new polymer with custom pattern: first 10 base, then 10 pattern
        newPattern = vcat(fill(Int32(1), 10), fill(Int32(2), 10))
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:custom, customPattern=newPattern, initShape="line")
        bulk2 = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1] 
        pat2 = poly2.monomers[2]
        intPot!(base2, -0.5) 
        intPot!(pat2, -0.3) 
        intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.1)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        # Verify: positions preserved
        restored_positions = [(poly2.pos[k][1], poly2.pos[k][2], poly2.pos[k][3]) for k in poly2.keys_vec[1:poly2.length]]
        @test restored_positions == saved_positions

        # Verify: pattern has changed - check intIdx in pos[key][6]
        key = poly2.head
        chainPos = 1
        while key != 0 && chainPos <= poly2.length
            expectedIntIdx = newPattern[chainPos]
            @test poly2.pos[key][6] == expectedIntIdx
            key = poly2.pos[key][5]
            chainPos += 1
        end

        # Verify: lattice is consistent with new pattern
        key = poly2.head
        while key != 0
            x, y, z = poly2.pos[key][1], poly2.pos[key][2], poly2.pos[key][3]
            @test poly2.latS[x, y, z] == poly2.pos[key][6]
            key = poly2.pos[key][5]
        end

        # Verify: patternTails is correctly set for the new pattern
        @test poly2.patternTails !== nothing
        @test haskey(poly2.patternTails, pat2.intIdx)

        # Verify: simulation can continue
        nmoves = 0
        for _ in 1:2000
            move!(poly2) && (nmoves += 1)
            move!(bulk2) && (nmoves += 1)
        end
        @test nmoves > 0

        rm("$tmpfile.jls")
    end

    @testset "Pattern change: custom → different custom" begin
        # Start with one custom pattern, restore with a different custom pattern
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        
        # Initial pattern: pattern type in first half (contiguous block)
        oldPattern = vcat(fill(Int32(2), 10), fill(Int32(1), 10))
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:custom, customPattern=oldPattern, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1] 
        pat = poly.monomers[2]
        intPot!(base, -0.5) 
        intPot!(pat, -0.3) 
        intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:3000 
            move!(poly) 
            move!(bulk) 
        end

        # Save positions
        saved_positions = deepcopy(poly.pos)

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        # Restore with DIFFERENT custom pattern: pattern type in middle (contiguous block)
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        
        newPattern = vcat(fill(Int32(1), 8), fill(Int32(2), 4), fill(Int32(1), 8))
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:custom, customPattern=newPattern, initShape="line")
        bulk2 = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1] 
        pat2 = poly2.monomers[2]
        intPot!(base2, -0.5) 
        intPot!(pat2, -0.3) 
        intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.1)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        # Verify: XYZ positions preserved (first 3 coords)
        for key in keys(saved_positions)
            @test poly2.pos[key][1:3] == saved_positions[key][1:3]
        end

        # Verify: pattern is the new one
        key = poly2.head
        chainPos = 1
        while key != 0 && chainPos <= poly2.length
            @test poly2.pos[key][6] == newPattern[chainPos]
            key = poly2.pos[key][5]
            chainPos += 1
        end

        # Verify: patternTails reflects new pattern
        @test poly2.patternTails !== nothing
        
        # Verify: hamiltonian is recalculated correctly
        ham!(poly2)
        ham_before = poly2.ham
        
        # Run some moves and check deltaH consistency
        drift = 0.0
        for _ in 1:500
            ham!(poly2)
            hamBefore = poly2.ham
            if move!(poly2)
                ham!(poly2)
                drift += abs((hamBefore + poly2.deltaH) - poly2.ham)
            end
        end
        @test drift < 1e-6

        rm("$tmpfile.jls")
    end

    @testset "Pattern preserved when same customPattern on restore" begin
        # Verify that identical patterns don't trigger recalculation
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        
        pattern = vcat(fill(Int32(1), 10), fill(Int32(2), 10))
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:custom, customPattern=pattern, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1] 
        pat = poly.monomers[2]
        intPot!(base, -0.5) 
        intPot!(pat, -0.3) 
        intPot!(bulk, -0.5)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:3000 
            move!(poly) 
            move!(bulk) 
        end

        ham!(poly)
        ham_before = poly.ham
        saved_patternTails = deepcopy(poly.patternTails)

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        # Restore with SAME pattern
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:custom, customPattern=copy(pattern), initShape="line")
        bulk2 = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1] 
        pat2 = poly2.monomers[2]
        intPot!(base2, -0.5) 
        intPot!(pat2, -0.3) 
        intPot!(bulk2, -0.5)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        # Ham should be identical (no pattern change)
        @test isapprox(poly2.ham, ham_before, atol=1e-10)
        
        # patternTails should match (from saved, not recalculated)
        @test poly2.patternTails == saved_patternTails

        rm("$tmpfile.jls")
    end

    @testset "Heterogeneous polymer with pattern + interaction change" begin
        # Full test: change both pattern AND interaction strengths
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        
        oldPattern = vcat(fill(Int32(1), 15), fill(Int32(2), 5))
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:custom, customPattern=oldPattern, initShape="line")
        bulk = FastBulkMonomer(length=60, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1] 
        pat = poly.monomers[2]
        intPot!(base, -0.2)  # weak
        intPot!(pat, -0.1) 
        intPot!(bulk, -0.2)
        intPot!(base, bulk, -0.05)
        intPot!(pat, bulk, -0.05)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:5000 
            move!(poly) 
            move!(bulk) 
        end

        saved_positions = deepcopy(poly.pos)

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        # Restore with BOTH different pattern AND stronger interactions
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        
        newPattern = vcat(fill(Int32(1), 5), fill(Int32(2), 10), fill(Int32(1), 5))
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:custom, customPattern=newPattern, initShape="line")
        bulk2 = FastBulkMonomer(length=60, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1] 
        pat2 = poly2.monomers[2]
        intPot!(base2, -2.0)  # 10x stronger
        intPot!(pat2, -1.0) 
        intPot!(bulk2, -2.0)
        intPot!(base2, bulk2, -0.5)
        intPot!(pat2, bulk2, -0.5)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        # Verify: XYZ positions preserved
        for key in keys(saved_positions)
            @test poly2.pos[key][1:3] == saved_positions[key][1:3]
        end

        # Verify: new pattern applied
        key = poly2.head
        chainPos = 1
        while key != 0 && chainPos <= poly2.length
            @test poly2.pos[key][6] == newPattern[chainPos]
            key = poly2.pos[key][5]
            chainPos += 1
        end

        # Verify: lattice consistency
        key = poly2.head
        while key != 0
            x, y, z = poly2.pos[key][1], poly2.pos[key][2], poly2.pos[key][3]
            @test poly2.latS[x, y, z] == poly2.pos[key][6]
            key = poly2.pos[key][5]
        end

        # Verify: simulation continues with new parameters
        nmoves = 0
        for _ in 1:3000
            move!(poly2) && (nmoves += 1)
            move!(bulk2) && (nmoves += 1)
        end
        @test nmoves > 0

        # Verify: deltaH tracking works correctly with new pattern
        drift = 0.0
        for _ in 1:500
            ham!(poly2)
            hamBefore = poly2.ham
            if move!(poly2)
                ham!(poly2)
                drift += abs((hamBefore + poly2.deltaH) - poly2.ham)
            end
        end
        @test drift < 1e-6

        rm("$tmpfile.jls")
    end

    @testset "Error: pattern too short for saved polymer" begin
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        intPot!(poly.monomers[1], -0.5)
        updateMonomers!(poly, poly.monomers[1].intIdx)
        init!(l, [poly])

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly])

        # Try to restore with pattern that's too short
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=1)
        
        shortPattern = fill(Int32(1), 10)  # Only 10 elements, but polymer is 20
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:custom, customPattern=shortPattern, initShape="line")
        initMonomers!(poly2, defaultChemPot=1.0)
        intPot!(poly2.monomers[1], -0.5)
        updateMonomers!(poly2, poly2.monomers[1].intIdx)

        @test_throws ErrorException restore!(tmpfile, l2, [poly2])

        rm("$tmpfile.jls")
    end

    @testset "Error: pattern contains invalid intIdx" begin
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        intPot!(poly.monomers[1], -0.5)
        updateMonomers!(poly, poly.monomers[1].intIdx)
        init!(l, [poly])

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly])

        # Try to restore with pattern containing invalid intIdx (3 doesn't exist)
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=1)
        
        invalidPattern = vcat(fill(Int32(1), 10), fill(Int32(3), 10))  # intIdx=3 doesn't exist
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:custom, customPattern=invalidPattern, initShape="line")
        initMonomers!(poly2, defaultChemPot=1.0)
        intPot!(poly2.monomers[1], -0.5)
        updateMonomers!(poly2, poly2.monomers[1].intIdx)

        @test_throws ErrorException restore!(tmpfile, l2, [poly2])

        rm("$tmpfile.jls")
    end

    @testset "Multiple restore cycles with pattern changes" begin
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        
        pattern1 = fill(Int32(1), 20)  # uniform base
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:custom, customPattern=pattern1, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        pat = poly.monomers[2]
        intPot!(base, -0.5)
        intPot!(pat, -0.3)
        intPot!(bulk, -0.5)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        tmpfile = tempname()

        # Patterns for each cycle
        patterns = [
            fill(Int32(1), 20),  # uniform
            vcat(fill(Int32(1), 10), fill(Int32(2), 10)),  # half-half
            vcat(fill(Int32(2), 5), fill(Int32(1), 10), fill(Int32(2), 5)),  # sandwich
        ]

        for (cycle, newPattern) in enumerate(patterns)
            for _ in 1:2000 
                move!(poly) 
                move!(bulk) 
            end

            checkpoint!(tmpfile, [poly, bulk])

            # Restore with next pattern
            Droplet.SPECIES_ID_COUNTER[] = 0 
            Droplet.KEYCOUNT[] = 0
            l = Lattice(boxlength=32, nspecies=2)
            poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                            pattern=:custom, customPattern=copy(newPattern), initShape="line")
            bulk = FastBulkMonomer(length=50, initShape="uniform")
            initMonomers!(poly, defaultChemPot=1.0)
            base = poly.monomers[1]
            pat = poly.monomers[2]
            intPot!(base, -0.5)
            intPot!(pat, -0.3)
            intPot!(bulk, -0.5)
            updateMonomers!(poly, base.intIdx)

            restore!(tmpfile, l, [poly, bulk])

            # Verify pattern applied correctly
            key = poly.head
            chainPos = 1
            while key != 0 && chainPos <= poly.length
                @test poly.pos[key][6] == newPattern[chainPos]
                key = poly.pos[key][5]
                chainPos += 1
            end

            # Verify simulation can continue
            nmoves = 0
            for _ in 1:500
                move!(poly) && (nmoves += 1)
                move!(bulk) && (nmoves += 1)
            end
            @test nmoves > 0
        end

        rm("$tmpfile.jls")
    end

end

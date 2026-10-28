using Test
using Droplet

@testset "restore! full round-trip" begin

    @testset "checkpoint → restore produces identical positions" begin
        Droplet.SPECIES_ID_COUNTER[] = 0; Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=100, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, -0.5); intPot!(pat, -0.5); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:5000; move!(poly); move!(bulk); end

        evolved_poly_pos = deepcopy(poly.pos)
        evolved_poly_head = poly.head
        evolved_poly_length = poly.length
        evolved_bulk_pos = deepcopy(bulk.pos)

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        # Fresh system (reset counters to simulate new session)
        Droplet.SPECIES_ID_COUNTER[] = 0; Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2, pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=100, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]; pat2 = poly2.monomers[2]
        intPot!(base2, -0.5); intPot!(pat2, -0.5); intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.1)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        @test poly2.pos == evolved_poly_pos
        @test poly2.head == evolved_poly_head
        @test poly2.length == evolved_poly_length
        @test bulk2.pos == evolved_bulk_pos

        rm("$tmpfile.jls")
    end

    @testset "Restored system can continue running moves" begin
        Droplet.SPECIES_ID_COUNTER[] = 0; Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=30, nMonomerTypes=1, pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, -0.5); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:2000; move!(poly); move!(bulk); end

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        Droplet.SPECIES_ID_COUNTER[] = 0; Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=30, nMonomerTypes=1, pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]
        intPot!(base2, -0.5); intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.1)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        nmoves = 0
        for _ in 1:5000
            move!(poly2) && (nmoves += 1)
            move!(bulk2) && (nmoves += 1)
        end
        @test nmoves > 0

        rm("$tmpfile.jls")
    end

    @testset "Restore with changed intPot uses new potentials" begin
        Droplet.SPECIES_ID_COUNTER[] = 0; Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, -0.5); intPot!(pat, -0.3); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:2000; move!(poly); move!(bulk); end
        ham!(poly)
        ham_before = poly.ham
        saved_pos = deepcopy(poly.pos)

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        # Re-run with DIFFERENT potentials
        Droplet.SPECIES_ID_COUNTER[] = 0; Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2, pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]; pat2 = poly2.monomers[2]
        intPot!(base2, -2.0); intPot!(pat2, -1.5)  # changed
        intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.1)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        @test poly2.pos == saved_pos           # same positions
        @test !isapprox(poly2.ham, ham_before, atol=1e-10)  # different ham

        rm("$tmpfile.jls")
    end

    @testset "Restore with changed chemPot uses new chemPot" begin
        Droplet.SPECIES_ID_COUNTER[] = 0; Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pat, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        for _ in 1:1000; move!(poly); end
        ham!(poly)
        ham_before = poly.ham

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly])

        Droplet.SPECIES_ID_COUNTER[] = 0; Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=1)
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2, pattern=:uniform, initShape="line")
        initMonomers!(poly2, defaultChemPot=5.0)  # changed
        base2 = poly2.monomers[1]; pat2 = poly2.monomers[2]
        intPot!(base2, 0.0); intPot!(pat2, 0.0)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2])

        @test !isapprox(poly2.ham, ham_before, atol=1e-10)

        rm("$tmpfile.jls")
    end

    @testset "Heterogeneous polymer pattern survives round-trip" begin
        Droplet.SPECIES_ID_COUNTER[] = 0; Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=64, nspecies=1)
        poly = Polymer(length=10, maxlength=30, nMonomerTypes=2, pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pat, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        growPolymer!(poly)
        patternVec = vcat(fill(base.intIdx, 10), fill(pat.intIdx, 10), fill(base.intIdx, 10))
        updatePattern!(l, poly, patternVec)

        saved_tails = copy(poly.patternTails)
        saved_custom = copy(poly.customPattern)

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly])

        Droplet.SPECIES_ID_COUNTER[] = 0; Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=64, nspecies=1)
        poly2 = Polymer(length=10, maxlength=30, nMonomerTypes=2, pattern=:uniform, initShape="line")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]; pat2 = poly2.monomers[2]
        intPot!(base2, 0.0); intPot!(pat2, 0.0)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2])

        @test poly2.pattern == :custom
        @test poly2.customPattern == saved_custom
        @test poly2.patternTails == saved_tails
        @test poly2.baseIntIdx == base.intIdx

        # Verify lattice has correct intIdx values
        key = poly2.head
        while key != 0
            x, y, z = poly2.pos[key][1], poly2.pos[key][2], poly2.pos[key][3]
            @test poly2.latS[x, y, z] == poly2.pos[key][6]
            key = poly2.pos[key][5]
        end

        rm("$tmpfile.jls")
    end

    @testset "Species count mismatch raises error" begin
        Droplet.SPECIES_ID_COUNTER[] = 0; Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=10, maxlength=10, nMonomerTypes=1, pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        intPot!(poly.monomers[1], 0.0)
        updateMonomers!(poly, poly.monomers[1].intIdx)
        init!(l, [poly])

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly])

        # Try to restore with extra species
        Droplet.SPECIES_ID_COUNTER[] = 0; Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=10, maxlength=10, nMonomerTypes=1, pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=30, initShape="uniform")

        @test_throws ErrorException restore!(tmpfile, l2, [poly2, bulk2])

        rm("$tmpfile.jls")
    end
end

using Test
using Droplet

@testset "Full simulation round-trip integration" begin

    @testset "Hamiltonian continuity across checkpoint" begin
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=100, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1] 
        pat = poly.monomers[2]
        intPot!(base, -0.5) 
        intPot!(pat, -0.3) 
        intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1) 
        intPot!(pat, bulk, -0.2)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:10000 
            move!(poly) 
            move!(bulk) 
        end

        ham!(poly) 
        ham!(bulk)
        ham_poly_end = poly.ham
        ham_bulk_end = bulk.ham

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        # Restore with SAME parameters
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=2, 
                        pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=100, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1] 
        pat2 = poly2.monomers[2]
        intPot!(base2, -0.5) 
        intPot!(pat2, -0.3) 
        intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.1) 
        intPot!(pat2, bulk2, -0.2)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        @test isapprox(poly2.ham, ham_poly_end, atol=1e-10)
        @test isapprox(bulk2.ham, ham_bulk_end, atol=1e-10)

        rm("$tmpfile.jls")
    end

    @testset "deltaHam consistency after restore" begin
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1, 
                        pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, -0.5) 
        intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.2)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:3000 
            move!(poly) 
            move!(bulk) 
        end

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=1, 
                        pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1]
        intPot!(base2, -0.5) 
        intPot!(bulk2, -0.5)
        intPot!(base2, bulk2, -0.2)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        # Run moves and verify deltaH tracks actual ham changes
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

    @testset "Lattice self-consistency after restore + moves" begin
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=30, nMonomerTypes=2, 
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

        for _ in 1:3000 
            move!(poly) 
            move!(bulk) 
        end

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=30, nMonomerTypes=2, 
                        pattern=:uniform, initShape="line")
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

        # Run more moves
        for _ in 1:5000 
            move!(poly2) 
            move!(bulk2) 
        end

        # Verify polymer lattice consistency
        key = poly2.head
        while key != 0
            x, y, z = poly2.pos[key][1], poly2.pos[key][2], poly2.pos[key][3]
            @test poly2.latS[x, y, z] == poly2.pos[key][6]
            key = poly2.pos[key][5]
        end

        # Verify bulk lattice consistency
        for j in 1:bulk2.length
            x, y, z = bulk2.pos[j][1], bulk2.pos[j][2], bulk2.pos[j][3]
            @test bulk2.latS[x, y, z] == j
        end

        rm("$tmpfile.jls")
    end

    @testset "Multiple checkpoint-restore cycles" begin
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1, 
                        pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, -0.5) 
        intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        tmpfile = tempname()

        for cycle in 1:3
            for _ in 1:2000 
                move!(poly) 
                move!(bulk) 
            end

            ham!(poly) 
            ham!(bulk)
            ham_poly = poly.ham 
            ham_bulk = bulk.ham

            checkpoint!(tmpfile, [poly, bulk])

            # Fresh restore
            Droplet.SPECIES_ID_COUNTER[] = 0 
            Droplet.KEYCOUNT[] = 0
            l = Lattice(boxlength=32, nspecies=2)
            poly = Polymer(length=20, maxlength=20, nMonomerTypes=1, 
                            pattern=:uniform, initShape="line")
            bulk = FastBulkMonomer(length=50, initShape="uniform")
            initMonomers!(poly, defaultChemPot=1.0)
            base = poly.monomers[1]
            intPot!(base, -0.5) 
            intPot!(bulk, -0.5)
            intPot!(base, bulk, -0.1)
            updateMonomers!(poly, base.intIdx)

            restore!(tmpfile, l, [poly, bulk])

            # Ham should match what we had before checkpoint
            @test isapprox(poly.ham, ham_poly, atol=1e-10)
            @test isapprox(bulk.ham, ham_bulk, atol=1e-10)

            # Can continue running
            nmoves = 0
            for _ in 1:1000
                move!(poly) && (nmoves += 1)
                move!(bulk) && (nmoves += 1)
            end
            @test nmoves > 0
        end

        rm("$tmpfile.jls")
    end

    @testset "Parameter change across checkpoint preserves positions (homogeneous)" begin
        # Simulate: run with weak interactions → checkpoint → restore with 
        # strong interactions on a homogeneous polymer
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1, 
                        pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=80, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1] 
        intPot!(base, -0.2) 
        intPot!(bulk, -0.2)
        intPot!(base, bulk, -0.05)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:5000 
            move!(poly) 
            move!(bulk) 
        end

        saved_poly_pos = deepcopy(poly.pos)
        saved_bulk_pos = deepcopy(bulk.pos)

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        # Restore with STRONG interactions
        Droplet.SPECIES_ID_COUNTER[] = 0 
        Droplet.KEYCOUNT[] = 0
        l2 = Lattice(boxlength=32, nspecies=2)
        poly2 = Polymer(length=20, maxlength=20, nMonomerTypes=1, 
                        pattern=:uniform, initShape="line")
        bulk2 = FastBulkMonomer(length=80, initShape="uniform")
        initMonomers!(poly2, defaultChemPot=1.0)
        base2 = poly2.monomers[1] 
        intPot!(base2, -2.0) 
        intPot!(bulk2, -2.0)  # 10x stronger
        intPot!(base2, bulk2, -0.5)
        updateMonomers!(poly2, base2.intIdx)

        restore!(tmpfile, l2, [poly2, bulk2])

        # Positions identical
        @test poly2.pos == saved_poly_pos
        @test bulk2.pos == saved_bulk_pos

        # But system can continue evolving with new potentials
        nmoves = 0
        for _ in 1:3000
            move!(poly2) && (nmoves += 1)
            move!(bulk2) && (nmoves += 1)
        end
        @test nmoves > 0

        rm("$tmpfile.jls")
    end
end

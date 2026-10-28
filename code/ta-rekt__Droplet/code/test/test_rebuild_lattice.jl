using Test
using Droplet

@testset "_rebuildLattice! reconstructs lattice" begin

    @testset "Polymer lattice populated correctly" begin
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pat, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        for _ in 1:2000; move!(poly); end

        # Wipe lattice and rebuild
        l.lat .= 0
        Droplet._rebuildLattice!(l, [poly])

        # Verify every bead is on the lattice with correct intIdx
        key = poly.head
        while key != 0
            x, y, z = poly.pos[key][1], poly.pos[key][2], poly.pos[key][3]
            @test poly.latS[x, y, z] == poly.pos[key][6]
            key = poly.pos[key][5]
        end
    end

    @testset "FastBulkMonomer lattice populated correctly" begin
        l = Lattice(boxlength=32, nspecies=1)
        fbulk = FastBulkMonomer(length=100, initShape="uniform")
        intPot!(fbulk, -0.5)
        init!(l, [fbulk])

        for _ in 1:1000; move!(fbulk); end

        l.lat .= 0
        Droplet._rebuildLattice!(l, [fbulk])

        for j in 1:fbulk.length
            x, y, z = fbulk.pos[j][1], fbulk.pos[j][2], fbulk.pos[j][3]
            @test fbulk.latS[x, y, z] == j
        end
    end

    @testset "movableIdxs rebuilt for FastBulkMonomer" begin
        l = Lattice(boxlength=32, nspecies=1)
        fbulk = FastBulkMonomer(length=100, initShape="uniform")
        intPot!(fbulk, -0.5)
        init!(l, [fbulk])

        for _ in 1:1000; move!(fbulk); end
        expected = copy(fbulk.movableIdxs)

        fbulk.movableIdxs = Set{Int32}()
        l.lat .= 0
        Droplet._rebuildLattice!(l, [fbulk])

        @test fbulk.movableIdxs == expected
    end

    @testset "KEYCOUNT recomputed correctly" begin
        l = Lattice(boxlength=64, nspecies=1)
        poly = Polymer(length=10, maxlength=30, nMonomerTypes=1, pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        intPot!(poly.monomers[1], 0.0)
        updateMonomers!(poly, poly.monomers[1].intIdx)
        init!(l, [poly])

        growPolymer!(poly)
        maxKey = maximum(keys(poly.pos))

        Droplet.KEYCOUNT[] = 0
        l.lat .= 0
        Droplet._rebuildLattice!(l, [poly])

        @test Droplet.KEYCOUNT[] >= maxKey
    end

    @testset "Hamiltonians recomputed after rebuild" begin
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1, pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=100, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, -0.5); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:2000; move!(poly); move!(bulk); end

        ham!(poly); ham!(bulk)
        expected_poly = poly.ham
        expected_bulk = bulk.ham

        poly.ham = -999.0; bulk.ham = -999.0
        l.lat .= 0
        Droplet._rebuildLattice!(l, [poly, bulk])

        @test isapprox(poly.ham, expected_poly, atol=1e-10)
        @test isapprox(bulk.ham, expected_bulk, atol=1e-10)
    end

    @testset "Mixed species: Polymer + BulkMonomer + FastBulkMonomer" begin
        l = Lattice(boxlength=32, nspecies=3)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1, pattern=:uniform, initShape="line")
        bm = BulkMonomer(length=30, initShape="uniform")
        fbm = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, 0.0); intPot!(bm, -0.3); intPot!(fbm, -0.5)
        intPot!(base, bm, -0.1); intPot!(base, fbm, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bm, fbm])

        for _ in 1:2000; move!(poly); move!(bm); move!(fbm); end

        l.lat .= 0
        Droplet._rebuildLattice!(l, [poly, bm, fbm])

        # Spot-check polymer
        key = poly.head
        while key != 0
            x, y, z = poly.pos[key][1], poly.pos[key][2], poly.pos[key][3]
            @test poly.latS[x, y, z] != 0
            key = poly.pos[key][5]
        end

        # Spot-check bulk monomer
        for j in 1:bm.length
            x, y, z = bm.pos[j][1], bm.pos[j][2], bm.pos[j][3]
            @test bm.latS[x, y, z] == j
        end

        # Spot-check fast bulk
        for j in 1:fbm.length
            x, y, z = fbm.pos[j][1], fbm.pos[j][2], fbm.pos[j][3]
            @test fbm.latS[x, y, z] == j
        end
    end

    @testset "speciesDir rebuilt correctly" begin
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=10, maxlength=10, nMonomerTypes=1, pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=30, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        intPot!(poly.monomers[1], 0.0); intPot!(bulk, -0.5)
        intPot!(poly.monomers[1], bulk, -0.1)
        updateMonomers!(poly, poly.monomers[1].intIdx)
        init!(l, [poly, bulk])

        l.lat .= 0; l.speciesDir = nothing
        Droplet._rebuildLattice!(l, [poly, bulk])

        @test l.speciesDir[poly.id] == 1
        @test l.speciesDir[bulk.id] == 2
    end
end

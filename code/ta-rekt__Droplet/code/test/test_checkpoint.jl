using Test
using Droplet
using Serialization

@testset "checkpoint! saves mutable state" begin

    @testset "Creates file and stores correct species count" begin
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pat, 0.0); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])

        @test isfile("$tmpfile.jls")
        data = Serialization.deserialize("$tmpfile.jls")
        @test length(data[:species]) == 2
        @test data[:version] == 1
        rm("$tmpfile.jls")
    end

    @testset "Polymer fields stored correctly" begin
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pat, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        # Apply a custom pattern for non-trivial patternTails
        newpat = vcat(fill(base.intIdx, 6), fill(pat.intIdx, 8), fill(base.intIdx, 6))
        updatePattern!(l, poly, newpat)

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly])
        sp = Serialization.deserialize("$tmpfile.jls")[:species][1]

        @test sp[:type] == :Polymer
        @test sp[:id] == poly.id
        @test sp[:pos] == poly.pos
        @test sp[:head] == poly.head
        @test sp[:length] == poly.length
        @test sp[:keys_vec] == poly.keys_vec
        @test sp[:pattern] == poly.pattern
        @test sp[:customPattern] == poly.customPattern
        @test sp[:patternTails] == poly.patternTails
        @test sp[:baseIntIdx] == poly.baseIntIdx
        rm("$tmpfile.jls")
    end

    @testset "BulkMonomer fields stored correctly" begin
        l = Lattice(boxlength=32, nspecies=1)
        bulk = BulkMonomer(length=100, initShape="uniform")
        intPot!(bulk, -0.5)
        init!(l, [bulk])

        tmpfile = tempname()
        checkpoint!(tmpfile, [bulk])
        sp = Serialization.deserialize("$tmpfile.jls")[:species][1]

        @test sp[:type] == :BulkMonomer
        @test sp[:id] == bulk.id
        @test sp[:pos] == bulk.pos
        @test sp[:length] == bulk.length
        rm("$tmpfile.jls")
    end

    @testset "FastBulkMonomer fields stored correctly" begin
        l = Lattice(boxlength=32, nspecies=1)
        fbulk = FastBulkMonomer(length=100, initShape="uniform")
        intPot!(fbulk, -0.5)
        init!(l, [fbulk])

        tmpfile = tempname()
        checkpoint!(tmpfile, [fbulk])
        sp = Serialization.deserialize("$tmpfile.jls")[:species][1]

        @test sp[:type] == :FastBulkMonomer
        @test sp[:id] == fbulk.id
        @test sp[:pos] == fbulk.pos
        @test sp[:length] == fbulk.length
        rm("$tmpfile.jls")
    end

    @testset "Captures evolved state after simulation" begin
        l = Lattice(boxlength=32, nspecies=2)
        poly = Polymer(length=10, maxlength=30, nMonomerTypes=1, pattern=:uniform, initShape="line")
        bulk = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, -0.5); intPot!(bulk, -0.5)
        intPot!(base, bulk, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bulk])

        for _ in 1:5000; move!(poly); move!(bulk); end

        evolved_pos = deepcopy(poly.pos)
        evolved_length = poly.length

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bulk])
        sp = Serialization.deserialize("$tmpfile.jls")[:species][1]

        @test sp[:pos] == evolved_pos
        @test sp[:length] == evolved_length
        rm("$tmpfile.jls")
    end

    @testset "Mixed species list ordering preserved" begin
        l = Lattice(boxlength=32, nspecies=3)
        poly = Polymer(length=10, maxlength=10, nMonomerTypes=1, pattern=:uniform, initShape="line")
        bm = BulkMonomer(length=30, initShape="uniform")
        fbm = FastBulkMonomer(length=50, initShape="uniform")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, 0.0); intPot!(bm, -0.3); intPot!(fbm, -0.5)
        intPot!(base, bm, -0.1); intPot!(base, fbm, -0.1)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly, bm, fbm])

        tmpfile = tempname()
        checkpoint!(tmpfile, [poly, bm, fbm])
        data = Serialization.deserialize("$tmpfile.jls")[:species]

        @test data[1][:type] == :Polymer
        @test data[2][:type] == :BulkMonomer
        @test data[3][:type] == :FastBulkMonomer
        @test data[1][:id] == poly.id
        @test data[2][:id] == bm.id
        @test data[3][:id] == fbm.id
        rm("$tmpfile.jls")
    end
end

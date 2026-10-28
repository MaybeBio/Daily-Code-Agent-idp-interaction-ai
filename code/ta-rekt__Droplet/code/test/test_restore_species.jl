using Test
using Droplet
using Serialization
using OrderedCollections

@testset "_restoreSpecies! overwrites mutable state" begin

    @testset "Polymer pos and chain topology restored" begin
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pat, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        for _ in 1:2000; move!(poly); end

        # Snapshot evolved state
        saved = Dict{Symbol,Any}(
            :type => :Polymer, :id => poly.id,
            :pos => deepcopy(poly.pos), :head => poly.head,
            :length => poly.length, :keys_vec => copy(poly.keys_vec),
            :pattern => poly.pattern, :customPattern => poly.customPattern,
            :patternTails => poly.patternTails === nothing ? nothing : copy(poly.patternTails),
            :baseIntIdx => poly.baseIntIdx,
        )
        evolved_pos = deepcopy(poly.pos)
        evolved_head = poly.head

        # Corrupt state to verify restore overwrites it
        empty!(poly.pos)
        poly.head = Int32(-1)

        Droplet._restoreSpecies!(poly, saved)

        @test poly.pos == evolved_pos
        @test poly.head == evolved_head
        @test poly.length == saved[:length]
        @test poly.keys_vec == saved[:keys_vec]
    end

    @testset "Polymer scratch buffers re-allocated" begin
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1, pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        intPot!(poly.monomers[1], 0.0)
        updateMonomers!(poly, poly.monomers[1].intIdx)
        init!(l, [poly])

        saved = Dict{Symbol,Any}(
            :type => :Polymer, :id => poly.id,
            :pos => deepcopy(poly.pos), :head => poly.head,
            :length => poly.length, :keys_vec => copy(poly.keys_vec),
            :pattern => poly.pattern, :customPattern => nothing,
            :patternTails => nothing, :baseIntIdx => poly.baseIntIdx,
        )

        # Corrupt scratch buffers
        poly.tempVec = nothing; poly.bondVec = nothing

        Droplet._restoreSpecies!(poly, saved)

        @test length(poly.tempVec) == 3
        @test length(poly.bondVec) == 3
        @test length(poly.moveVec) == 3
        @test length(poly.newPosi) == 3
        @test length(poly.newPosj) == 3
        @test poly.moveTypes !== nothing
        @test poly.moveTypeProbs !== nothing
        @test poly.deltaHP == 0.0
        @test poly.hastingsFactor == 1.0
    end

    @testset "BulkMonomer pos restored" begin
        l = Lattice(boxlength=32, nspecies=1)
        bulk = BulkMonomer(length=50, initShape="uniform")
        intPot!(bulk, -0.5)
        init!(l, [bulk])

        for _ in 1:1000; move!(bulk); end
        evolved_pos = deepcopy(bulk.pos)

        saved = Dict{Symbol,Any}(
            :type => :BulkMonomer, :id => bulk.id,
            :pos => evolved_pos, :length => bulk.length,
        )

        # Corrupt state to verify restore overwrites it
        empty!(bulk.pos)

        Droplet._restoreSpecies!(bulk, saved)

        @test bulk.pos == evolved_pos
        @test bulk.deltaH == 0.0
        @test bulk.deltaHP == 0.0
        @test bulk.freeNNs == Set{Int32}()
    end

    @testset "FastBulkMonomer pos restored with empty movableIdxs" begin
        l = Lattice(boxlength=32, nspecies=1)
        fbulk = FastBulkMonomer(length=50, initShape="uniform")
        intPot!(fbulk, -0.5)
        init!(l, [fbulk])

        for _ in 1:1000; move!(fbulk); end
        evolved_pos = deepcopy(fbulk.pos)

        saved = Dict{Symbol,Any}(
            :type => :FastBulkMonomer, :id => fbulk.id,
            :pos => evolved_pos, :length => fbulk.length,
        )

        init!(l, [fbulk])
        Droplet._restoreSpecies!(fbulk, saved)

        @test fbulk.pos == evolved_pos
        @test fbulk.movableIdxs == Set{Int32}()  # empty, rebuilt later by _rebuildLattice!
        @test fbulk.deltaH == 0.0
        @test fbulk.deltaHP == 0.0
        @test length(fbulk.bi) == 3
        @test length(fbulk.moveVec) == 2
    end

    @testset "Type mismatch raises error" begin
        poly = Polymer(length=10, maxlength=20, nMonomerTypes=1, pattern=:uniform, initShape="line")
        saved = Dict{Symbol,Any}(:type => :BulkMonomer, :id => poly.id,
            :pos => OrderedDict{Int32,Vector{Int32}}(), :length => Int32(10))
        @test_throws ErrorException Droplet._restoreSpecies!(poly, saved)
    end

    @testset "ID mismatch raises error" begin
        poly = Polymer(length=10, maxlength=20, nMonomerTypes=1, pattern=:uniform, initShape="line")
        saved = Dict{Symbol,Any}(:type => :Polymer, :id => Int32(poly.id + 100),
            :pos => OrderedDict{Int32,Vector{Int32}}(), :head => Int32(1),
            :length => Int32(10), :keys_vec => Int32[], :pattern => :uniform,
            :customPattern => nothing, :patternTails => nothing, :baseIntIdx => nothing)
        @test_throws ErrorException Droplet._restoreSpecies!(poly, saved)
    end

    @testset "chemPot and intMat are NOT overwritten" begin
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=1, pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]
        intPot!(base, -0.5)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        saved = Dict{Symbol,Any}(
            :type => :Polymer, :id => poly.id,
            :pos => deepcopy(poly.pos), :head => poly.head,
            :length => poly.length, :keys_vec => copy(poly.keys_vec),
            :pattern => poly.pattern, :customPattern => nothing,
            :patternTails => nothing, :baseIntIdx => poly.baseIntIdx,
        )

        # Change chemPot AFTER saving (simulating parameter change)
        base.chemPot = 2.0
        intMat_before = copy(base.intMat)

        Droplet._restoreSpecies!(poly, saved)

        # chemPot and intMat should be untouched
        @test base.chemPot == 2.0
        @test base.intMat == intMat_before
    end

    @testset "Pattern fields round-trip through checkpoint + restore" begin
        l = Lattice(boxlength=32, nspecies=1)
        poly = Polymer(length=20, maxlength=20, nMonomerTypes=2, pattern=:uniform, initShape="line")
        initMonomers!(poly, defaultChemPot=1.0)
        base = poly.monomers[1]; pat = poly.monomers[2]
        intPot!(base, 0.0); intPot!(pat, 0.0)
        updateMonomers!(poly, base.intIdx)
        init!(l, [poly])

        newpat = vcat(fill(base.intIdx, 6), fill(pat.intIdx, 8), fill(base.intIdx, 6))
        updatePattern!(l, poly, newpat)

        # Checkpoint → deserialize → _restoreSpecies!
        tmpfile = tempname()
        checkpoint!(tmpfile, [poly])
        data = Serialization.deserialize("$tmpfile.jls")

        # Wipe pattern fields
        poly.pattern = :uniform; poly.customPattern = nothing; poly.patternTails = nothing

        Droplet._restoreSpecies!(poly, data[:species][1])

        @test poly.pattern == :custom
        @test poly.customPattern == newpat
        @test poly.patternTails !== nothing
        @test poly.baseIntIdx == base.intIdx

        rm("$tmpfile.jls")
    end
end

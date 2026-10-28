# MARK: -intPot self
    @testset "intPot! (self)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 1
        )

        a = FastBulkMonomer(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Test setting self interaction potential
        intPot!(a, -2.0)

        init!(l, [a])
        @test a.intPot == -2.0

        # Test with BulkMonomer type as well
        b = BulkMonomer(;
            length = 10,
            intPot = -3.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        init!(l, [b])
        
        intPot!(b, -4.0)
        @test b.intPot == -4.0

        intIdxs = map(x -> x[2], b.intSpecies)

        @test b.intIdx in intIdxs
    end

# MARK: -intPot two species
    @testset "intPot! (two species)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 2
        )

        a = FastBulkMonomer(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulkMonomer(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Set cross-interaction potential
        intPot!(a, b, -3.0)

        init!(l, [a, b])

        # Both species should have interaction structures
        @test a.intSpecies !== nothing
        @test b.intSpecies !== nothing
        @test a.intMat !== nothing
        @test b.intMat !== nothing
        @test a.intIdx !== nothing
        @test b.intIdx !== nothing

        intIdxsA = map(x -> x[2], a.intSpecies)
        intIdxsB = map(x -> x[2], b.intSpecies)

        # each species' interaction list contains its own index
        @test a.intIdx in intIdxsA
        @test b.intIdx in intIdxsB

        # They should share the same matrix
        @test a.intMat === b.intMat

        # Matrix should be 2x2
        @test size(a.intMat) == (2, 2)

        # Matrix should be symmetric
        @test a.intMat[a.intIdx, b.intIdx] == -3.0
        @test a.intMat[b.intIdx, a.intIdx] == -3.0

        # Diagonal entries should equal self intPot
        @test a.intMat[a.intIdx, a.intIdx] == a.intPot
        @test a.intMat[b.intIdx, b.intIdx] == b.intPot
    end

# MARK: -intPot three species
    @testset "intPot! (three species)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 3
        )

        a = FastBulkMonomer(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulkMonomer(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        c = FastBulkMonomer(;
            length = 10,
            intPot = -3.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Set pairwise interactions
        intPot!(a, b, -4.0)
        intPot!(a, c, -5.0)
        intPot!(b, c, -6.0)

        init!(l, [a, b, c])

        # All three should share the same matrix
        @test a.intMat === b.intMat
        @test b.intMat === c.intMat

        # Matrix should be 3x3
        @test size(a.intMat) == (3, 3)

        # Check all pairwise interactions (symmetric)
        @test a.intMat[a.intIdx, b.intIdx] == -4.0
        @test a.intMat[b.intIdx, a.intIdx] == -4.0
        @test a.intMat[a.intIdx, c.intIdx] == -5.0
        @test a.intMat[c.intIdx, a.intIdx] == -5.0
        @test a.intMat[b.intIdx, c.intIdx] == -6.0
        @test a.intMat[c.intIdx, b.intIdx] == -6.0

        # Diagonal entries should equal self intPot
        @test a.intMat[a.intIdx, a.intIdx] == a.intPot
        @test a.intMat[b.intIdx, b.intIdx] == b.intPot
        @test a.intMat[c.intIdx, c.intIdx] == c.intPot

        intIdxsA = map(x -> x[2], a.intSpecies)
        intIdxsB = map(x -> x[2], b.intSpecies)

        # each species' interaction list contains its own index
        @test a.intIdx in intIdxsA
        @test b.intIdx in intIdxsB
    end

# MARK: -intPot update
    @testset "intPot! (update existing)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 2
        )

        a = FastBulkMonomer(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulkMonomer(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Set initial interaction
        intPot!(a, b, -3.0)

        init!(l, [a, b])

        @test a.intMat[a.intIdx, b.intIdx] == -3.0

        # Update interaction - should overwrite
        intPot!(a, b, -5.0)
        @test a.intMat[a.intIdx, b.intIdx] == -5.0
        @test a.intMat[b.intIdx, a.intIdx] == -5.0

        # Update in reverse order
        intPot!(b, a, -7.0)
        @test a.intMat[a.intIdx, b.intIdx] == -7.0
        @test a.intMat[b.intIdx, a.intIdx] == -7.0
    end

# MARK: -intPot self after pair
    @testset "intPot! (self after pair interaction)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 2
        )

        a = FastBulkMonomer(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulkMonomer(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Set pairwise interaction first
        intPot!(a, b, -3.0)

        init!(l, [a, b])

        # Now update self interaction - should update diagonal
        intPot!(a, -10.0)
        @test a.intPot == -10.0
        @test a.intMat[a.intIdx, a.intIdx] == -10.0

        # Cross interaction should be unchanged
        @test a.intMat[a.intIdx, b.intIdx] == -3.0
    end

# MARK: -getIntPot
    @testset "getIntPot" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 3
        )

        a = FastBulkMonomer(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulkMonomer(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        c = FastBulkMonomer(;
            length = 10,
            intPot = -3.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Before setting interactions, getIntPot should return nothing
        @test getIntPot(a, b) === nothing
        
        # Set interactions
        intPot!(a, b, -4.0)
        intPot!(a, c, -5.0)
        intPot!(b, c, -6.0)
        
        init!(l, [a, b, c])

        # Now getIntPot should return correct values
        @test getIntPot(a, b) == -4.0
        @test getIntPot(b, a) == -4.0  # symmetric
        @test getIntPot(a, c) == -5.0
        @test getIntPot(c, a) == -5.0
        @test getIntPot(b, c) == -6.0
        @test getIntPot(c, b) == -6.0

        # Self interactions
        @test getIntPot(a, a) == a.intPot
        @test getIntPot(b, b) == b.intPot
        @test getIntPot(c, c) == c.intPot
    end

# MARK: -intPot chain merge
    @testset "intPot! (chain merge)" begin
        # Test merging two separate interaction networks
        l = Lattice(;
            boxlength = 8,
            nspecies = 4
        )

        a = FastBulkMonomer(;
            length = 5,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulkMonomer(;
            length = 5,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        c = FastBulkMonomer(;
            length = 5,
            intPot = -3.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        d = FastBulkMonomer(;
            length = 5,
            intPot = -4.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Create two separate networks: (a,b) and (c,d)
        intPot!(a, b, -10.0)
        intPot!(c, d, -20.0)
        
        init!(l, [a, b, c, d])

        # Verify separate networks
        @test a.intMat === b.intMat
        @test c.intMat === d.intMat
        @test a.intMat !== c.intMat

        # Now connect the networks
        intPot!(b, c, -30.0)

        # All should now share the same matrix
        @test a.intMat === b.intMat
        @test b.intMat === c.intMat
        @test c.intMat === d.intMat

        # Matrix should be 4x4
        @test size(a.intMat) == (4, 4)

        # All previous interactions should be preserved
        @test getIntPot(a, b) == -10.0
        @test getIntPot(c, d) == -20.0
        @test getIntPot(b, c) == -30.0

        # Diagonal entries preserved
        @test a.intMat[a.intIdx, a.intIdx] == -1.0
        @test a.intMat[b.intIdx, b.intIdx] == -2.0
        @test a.intMat[c.intIdx, c.intIdx] == -3.0
        @test a.intMat[d.intIdx, d.intIdx] == -4.0

        intIdxsA = map(x -> x[2], a.intSpecies)
        intIdxsB = map(x -> x[2], b.intSpecies)
        intIdxsC = map(x -> x[2], c.intSpecies)
        intIdxsD = map(x -> x[2], d.intSpecies)

        # each species' interaction list contains its own index
        @test a.intIdx in intIdxsA
        @test b.intIdx in intIdxsB
        @test c.intIdx in intIdxsC
        @test d.intIdx in intIdxsD
    end

# MARK: -intPot with BulkMonomer
    @testset "intPot! (mixed BulkMonomer and FastBulkMonomer)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 2
        )

        a = BulkMonomer(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        b = FastBulkMonomer(;
            length = 10,
            intPot = -2.0,
            chemPot = 1.0,
            initShape = "uniform"
        )

        # Set cross-interaction between BulkMonomer and FastBulkMonomer
        intPot!(a, b, -5.0)
        
        init!(l, [a, b])

        # Both should have interaction structures
        @test a.intMat !== nothing
        @test b.intMat !== nothing
        @test a.intMat === b.intMat

        # Check values
        @test getIntPot(a, b) == -5.0
        @test getIntPot(b, a) == -5.0
    end


# MARK: -intPot same species
    @testset "intPot! (same species twice)" begin
        l = Lattice(;
            boxlength = 8,
            nspecies = 1
        )

        a = FastBulkMonomer(;
            length = 10,
            intPot = -1.0,
            chemPot = 1.0,
            initShape = "uniform"
        )
        
        # Set interaction with itself
        intPot!(a, a, -5.0)

        init!(l, [a])

        # Should update diagonal
        @test a.intMat !== nothing
        @test size(a.intMat) == (1, 1)
        @test a.intMat[1, 1] == -5.0
    end
    

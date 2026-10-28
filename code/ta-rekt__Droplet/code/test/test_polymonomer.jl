# MARK: -PolyMonomer init
@testset "Basic PolyMonomer initialization" begin
    # Create a simple polymer with 2 monomer types
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2
    )
    
    # Initialize monomers with default potentials
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    
    # Check that monomers were created (keyed by intIdx)
    @test poly.monomers !== nothing
    @test length(poly.monomers) == 2
    @test haskey(poly.monomers, Int32(1))
    @test haskey(poly.monomers, Int32(2))
    
    # Check that each monomer has the poly reference
    @test poly.monomers[1].poly === poly
    @test poly.monomers[2].poly === poly
    
    # Check that monomers inherit chemical and interaction potentials
    @test poly.monomers[1].chemPot == 0.0
    @test poly.monomers[2].chemPot == 0.0
    @test poly.monomers[1].intPot == 0.0
    @test poly.monomers[2].intPot == 0.0
    
    # Check that interaction structures are initialized
    @test poly.monomers[1].intSpecies !== nothing
    @test poly.monomers[1].intMat !== nothing
    @test poly.monomers[1].intIdx !== nothing
    
    # Check that both monomers share the same interaction matrix
    @test poly.monomers[1].intMat === poly.monomers[2].intMat
end

# MARK: -PolyMonomer pattern
@testset "PolyMonomer pattern assignment" begin
    l = Lattice(; boxlength = 10, nspecies = 1)
    
    # Helper to extract intIdx pattern from chain
    function getPatternFromChain(p::Polymer)
        result = Int32[]
        key = p.head
        while key != 0 && length(result) < p.length
            push!(result, p.pos[key][6])
            key = p.pos[key][5]
        end
        return result
    end
    
    # Test uniform pattern (default)
    poly1 = Polymer(; length = 4, maxlength = 8, nMonomerTypes = 2, pattern = :uniform)
    initMonomers!(poly1, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly1])
    pattern1 = getPatternFromChain(poly1)
    @test all(pattern1 .== 1)  # All beads are monomer type 1 (intIdx=1)
    
    # Test block pattern (valid: contiguous blocks)
    l2 = Lattice(; boxlength = 10, nspecies = 1)
    poly2 = Polymer(; length = 4, maxlength = 8, nMonomerTypes = 2, pattern = :block)
    initMonomers!(poly2, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l2, [poly2])
    pattern2 = getPatternFromChain(poly2)
    @test pattern2 == [1, 1, 2, 2]
    
    # Test custom contiguous pattern (valid)
    l3 = Lattice(; boxlength = 10, nspecies = 1)
    poly3 = Polymer(; length = 6, maxlength = 8, nMonomerTypes = 3, pattern = :custom, customPattern = [1, 2, 2, 3, 3, 1])
    initMonomers!(poly3, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l3, [poly3])
    pattern3 = getPatternFromChain(poly3)
    @test pattern3 == [1, 2, 2, 3, 3, 1]
    
    # Test alternating pattern throws error (non-contiguous blocks)
    l4 = Lattice(; boxlength = 10, nspecies = 1)
    poly4 = Polymer(; length = 4, maxlength = 8, nMonomerTypes = 2, pattern = :alternating)
    initMonomers!(poly4, defaultChemPot=0.0, defaultIntPot=0.0)
    @test_throws ErrorException init!(l4, [poly4])
    
    # Test non-contiguous custom pattern throws error
    l5 = Lattice(; boxlength = 10, nspecies = 1)
    poly5 = Polymer(; length = 4, maxlength = 8, nMonomerTypes = 2, pattern = :custom, customPattern = [1, 2, 1, 2])
    initMonomers!(poly5, defaultChemPot=0.0, defaultIntPot=0.0)
    @test_throws ErrorException init!(l5, [poly5])
end

# MARK: -PolyMonomer intPot
@testset "Setting interaction potentials between monomers" begin
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 3
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    
    dna = poly.monomers[1]
    promoter = poly.monomers[2]
    activator = poly.monomers[3]
    
    dna.name = "dna"
    promoter.name = "promoter"
    activator.name = "activator"
    
    # Set self-interaction
    intPot!(dna, 0.5)
    @test dna.intPot == 0.5
    @test dna.intMat[dna.intIdx, dna.intIdx] == 0.5
    
    # Set cross-interaction
    intPot!(promoter, activator, 2.0)
    @test promoter.intMat[promoter.intIdx, activator.intIdx] == 2.0
    @test activator.intMat[activator.intIdx, promoter.intIdx] == 2.0
    
    # Check that all monomers still share the same matrix
    @test dna.intMat === promoter.intMat
    @test promoter.intMat === activator.intMat
end

# MARK: -PolyMono FastBulk intPot
@testset "Interaction between FastBulkMonomer and PolyMonomer" begin
    l = Lattice(;
        boxlength = 4,
        nspecies = 3
    )
    
    tf = FastBulkMonomer(;
        length = 12,
        initShape = "uniform"
    )
    
    poly = Polymer(;
        length = l.boxlength,
        maxlength = 8,
        nMonomerTypes = 3
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    
    dna = poly.monomers[1]
    promoter = poly.monomers[2]
    
    dna.name = "dna"
    promoter.name = "promoter"
    
    # Set interaction potentials between bulk monomer and poly monomers
    intPot!(tf, dna, 0.3)
    intPot!(tf, promoter, 2.0)
    
    # Check that interactions were set correctly
    @test getIntPot(tf, dna) == 0.3
    @test getIntPot(dna, tf) == 0.3
    @test getIntPot(tf, promoter) == 2.0
    @test getIntPot(promoter, tf) == 2.0
    
    # Check that all share the same interaction matrix
    @test tf.intMat === dna.intMat
    @test dna.intMat === promoter.intMat
end

# MARK: -PolyMonomer mixed
@testset "Full initialization with mixed species" begin
    l = Lattice(;
        boxlength = 8,
        nspecies = 3
    )
    
    tf = FastBulkMonomer(;
        length = 32,
        initShape = "uniform"
    )
    
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2,
        pattern = :block
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    
    dna = poly.monomers[1]
    promoter = poly.monomers[2]
    
    dna.name = "dna"
    promoter.name = "promoter"
    
    intPot!(tf, 1.0)
    intPot!(tf, dna, 0.3)
    intPot!(tf, promoter, 2.0)
    
    # Initialize the species array with mixed types
    species = [tf, dna, promoter]
    
    # This should work with the new init! implementation
    init!(l, species)
    
    # Verify lattice was initialized
    @test l.lat !== nothing
    @test l.nspecies > 0
    
    # Verify polymer has lattice view
    @test poly.latS !== nothing
    @test poly.boxlength == 8
    
    # Verify bulk monomer has lattice view
    @test tf.latS !== nothing
    @test tf.boxlength == 8
end

# MARK: -PolyMonomer helpers
@testset "Position and lattice delegation helpers" begin
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2
    )
    
    l = Lattice(;
        boxlength = 8,
        nspecies = 1
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly])
    
    monomer = poly.monomers[1]
    
    # Test delegation helpers
    @test getLatS(monomer) === poly.latS
    @test getTempVec(monomer) === poly.tempVec
    @test getBoxlength(monomer) == poly.boxlength
    
    # Test with regular species for comparison
    tf = FastBulkMonomer(;
        length = 32,
        initShape = "uniform"
    )
    
    l2 = Lattice(;
        boxlength = 8,
        nspecies = 1
    )
    
    init!(l2, [tf])
    
    @test getLatS(tf) === tf.latS
    @test getTempVec(tf) === tf.tempVec
    @test getBoxlength(tf) == tf.boxlength
end

# MARK: -PolyMonomer nn/nnp
@testset "nn and nnp with PolyMonomer" begin
    l = Lattice(;
        boxlength = 8,
        nspecies = 2
    )
    
    poly = Polymer(;
        length = 4,
        maxlength = 8,
        nMonomerTypes = 2
    )
    
    tf = FastBulkMonomer(;
        length = 32,
        initShape = "uniform"
    )
    
    initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
    init!(l, [poly, tf])
    
    monomer = poly.monomers[1]
    
    # Test nn (same-species neighbors)
    pos = poly.pos[1][1:3]
    neighbors = nn(monomer, pos)
    @test neighbors >= 0
    @test neighbors <= 6
    
    # Test nnp (cross-species neighbors)
    neighbors_cross = nnp(monomer, tf, pos)
    @test neighbors_cross >= 0
    @test neighbors_cross <= 6
end
# Tests for initialization and move identification

@testset "Init, identify" begin
    l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    # Polymer
    p = Polymer(;
        length = 16,
        maxlength = 20,
        nMonomerTypes = 1
    )
    
    initMonomers!(p, defaultChemPot=1.0, defaultIntPot=1.0)

    # dilute monomers
    b = BulkMonomer(;
        length = 256
    )
    
    # Set interaction potential after init
    intPot!(b, 1.0)

    init!(l, [p, b])

    for i in 1:16
        @test p.latS[i,8,8] == 1
        @test l.lat[2,i,8,8] == 1
    end

    @test p.pos[1][4] == 0
    @test p.pos[5][4] == 4
    @test p.pos[5][5] == 6
    @test p.pos[16][5] == 0

    p.moveIdx = 5

    identifyBond!(p)

    @test p.ki == 5
    @test p.kj == 6
    @test p.kprev == 4
    @test p.knext == 7
    @test p.bi == [5,8,8]
    @test p.bj == [6,8,8]
    @test p.bprev == [4,8,8]
    @test p.bnext == [7,8,8]

    p.moveIdx = 1

    identifyBond!(p)

    @test p.ki == 1
    @test p.kj == 2
    @test p.kprev == 0

    p.moveIdx = 15

    identifyBond!(p)

    @test p.ki == 15
    @test p.kj == 16
    @test p.knext == 0
end

@testset "Movetypes" begin
    l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    # Polymer
    p = Polymer(;
        length = 16,
        maxlength = 20,
        nMonomerTypes = 1
    )
    
    initMonomers!(p, defaultChemPot=1.0, defaultIntPot=1.0)

    # dilute monomers
    b = BulkMonomer(;
        length = 256
    )

    # Set interaction potential
    intPot!(b, 1.0)

    # initialize the system
    init!(l, [p, b])

    ### ### ### Adding a bond

    p.moveIdx = 5
    identifyBond!(p)
    p.moveType = 1  # bond addition
    p.moveVec = [0,0,1]
    p.newPosi = p.bi .+ p.moveVec
    p.newPosj = p.bj .+ p.moveVec
    
    # Original length
    origLength = p.length
    
    makeMove!(p)
    p.length = length(keys(p.pos))

    p.moveIdx = 3
    identifyBond!(p)
    moveTypes!(p)

    @test p.moveTypes == [1,1,1,1]

    p.moveIdx = 4
    identifyBond!(p)
    moveTypes!(p)

    @test p.moveTypes == [1,1,2,1]

    p.moveIdx = 5
    identifyBond!(p)
    moveTypes!(p)

    @test p.moveTypes == [2,2,1,1]

    p.moveIdx = 17
    identifyBond!(p)
    moveTypes!(p)

    @test p.moveTypes == [1,1,1,3]

    p.moveIdx = 18
    identifyBond!(p)
    moveTypes!(p)

    @test p.moveTypes == [2,2,1,1]

    
    p.moveIdx = 18
    identifyBond!(p)
    p.moveType = 1  # bond addition
    p.moveVec = [0,1,0]
    p.newPosi = p.bi .+ p.moveVec
    p.newPosj = p.bj .+ p.moveVec
    makeMove!(p)
    p.length = length(keys(p.pos))

    p.moveIdx = 18
    identifyBond!(p)
    moveTypes!(p)

    @test p.moveTypes == [1,2,1,2]
end

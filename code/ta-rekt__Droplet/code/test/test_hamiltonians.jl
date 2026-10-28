# Hamiltonian tests for various species types

# Include new comprehensive Hamiltonian tests
include("test_polymer_hamiltonian.jl")
include("test_polymer_interacting.jl")

# Legacy integration tests
@testset "Hamiltonian (Polymer)" begin
    l = Lattice(;
        boxlength = 16,
        nspecies = 1
    )

    p = Polymer(;
        length = 16,
        maxlength = 64,
        nMonomerTypes = 1,
        initPos = [1,8,8],
        pattern=:uniform
    )
    
    initMonomers!(p)
    a = p.monomers[1]

    # self and interaction potentials
    intPot!(a, -1.0)

    slist = [p]

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = [0]
    frame = [1]
    frame_joint = 1
    nframes = 10000
    step = 1
    deltaH = [0.0]
    niter = 0
    moveMade = [false]

    hamsA = Vector{Float64}(undef, 0)
    hams_ = [hamsA]
    hams = []
    accDeltaHA = Vector{Float64}(undef, 0)
    accDeltaH_ = [accDeltaHA]
    accDeltaH = []
    
    # first frame, initial energy
    push!(hams_[1], getHam(a))
    push!(accDeltaH_[1], getHam(a))

    bulkSpecies = Vector{Any}(undef, 0)

    for (sp, idx) in a.intSpecies
        push!(hams, hams_[idx])
        push!(accDeltaH, accDeltaH_[idx])

        if typeof(sp) != PolyMonomer
            push!(bulkSpecies, (sp, idx))
        end
    end

    ham!(a)
    p.deltaH = 0.0

    while frame_joint <= nframes - 2
        for (s, i) in a.intSpecies
            if move!(s)
                nmoves[i] += 1
                deltaH[i] += getDeltaHam(s)

                push!(accDeltaH[i], accDeltaH[i][end] + deltaH[i])
                deltaH[i] = 0
                
                ham!(s)
                push!(hams[i], getHam(s))

                frame[i] += 1
                frame_joint += 1
            end
        end
    end

    @test hams[1] == accDeltaH[1]
end

@testset "Hamiltonian (homogeneous interacting polymer)" begin
    l = Lattice(;
        boxlength = 16,
        nspecies = 3
    )

    p = Polymer(;
        length = 16,
        maxlength = 64,
        nMonomerTypes = 1,
        initPos = [1,8,8],
        pattern=:uniform
    )
    
    initMonomers!(p)
    a = p.monomers[1]

    b = FastBulkMonomer(;
        length = 64,
        initShape = "uniform"
    )

    c = FastBulkMonomer(;
        length = 64,
        initShape = "uniform"
    )

    # self and interaction potentials
    intPot!(a, -1.0)
    intPot!(b, -2.0)
    intPot!(c, -3.0)
    intPot!(a, b, -0.5)
    intPot!(a, c, 0.5)

    slist = [p, b, c]

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = [0,0,0]
    frame = [1,1,1]
    frame_joint = 1
    nframes = 10000
    step = 1
    deltaH = [0.0,0.0,0.0]
    niter = 0
    moveMade = [false, false, false]

    hamsA = Vector{Float64}(undef, 0)
    hamsB = Vector{Float64}(undef, 0)
    hamsC = Vector{Float64}(undef, 0)
    hams_ = [hamsA, hamsB, hamsC]
    hams = []
    accDeltaHA = Vector{Float64}(undef, 0)
    accDeltaHB = Vector{Float64}(undef, 0)
    accDeltaHC = Vector{Float64}(undef, 0)
    accDeltaH_ = [accDeltaHA, accDeltaHB, accDeltaHC]
    accDeltaH = []
    
    # first frame, initial energy
    push!(hams_[1], getHam(a))
    push!(hams_[2], getHam(b))
    push!(hams_[3], getHam(c))
    push!(accDeltaH_[1], getHam(a))
    push!(accDeltaH_[2], getHam(b))
    push!(accDeltaH_[3], getHam(c))

    bulkSpecies = Vector{Any}(undef, 0)

    for (sp, idx) in a.intSpecies
        push!(hams, hams_[idx])
        push!(accDeltaH, accDeltaH_[idx])

        if typeof(sp) != PolyMonomer
            push!(bulkSpecies, (sp, idx))
        end
    end

    ham!(a)
    ham!(b)
    ham!(c)
    p.deltaH = 0.0
    b.deltaH = 0.0
    c.deltaH = 0.0
    p.deltaHP = 0.0
    b.deltaHP = 0.0
    c.deltaHP = 0.0

    while frame_joint <= nframes - 2
        for (s, i) in a.intSpecies
            if move!(s)
                nmoves[i] += 1
                deltaH[i] += getDeltaHam(s)

                # update deltaH for the partners
                for (partner, j) in a.intSpecies
                    if j != i
                        deltaH[j] += getDeltaHP(partner)
                    end
                end

                push!(accDeltaH[i], accDeltaH[i][end] + deltaH[i])
                deltaH[i] = 0
                
                ham!(s)
                push!(hams[i], getHam(s))

                frame[i] += 1
                frame_joint += 1
            end
        end
    end

    @test hams[1] == accDeltaH[1]
    @test hams[2] == accDeltaH[2]
    @test hams[3] == accDeltaH[3]
end

@testset "Hamiltonian (heterogeneous polymer without bulk)" begin
    l = Lattice(;
        boxlength = 20,
        nspecies = 1
    )

    p = Polymer(;
        length = 10,
        maxlength = 30,
        nMonomerTypes = 2,
        initPos = [1,10,10],
        pattern = :custom,
        customPattern = [1,1,1,1,2,2,2,2,1,1]
    )
    
    initMonomers!(p, defaultChemPot=1.0, defaultIntPot=-0.3)
    
    # Set different interaction potentials for each monomer type
    m1 = p.monomers[1]
    m2 = p.monomers[2]
    intPot!(m1, -0.3)  # base type self-interaction
    intPot!(m2, -0.5)  # pattern type self-interaction

    slist = [p]

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = 0
    frame = 1
    nframes = 5000
    
    hams = Vector{Float64}(undef, 0)
    accDeltaH = Vector{Float64}(undef, 0)
    
    # first frame, initial energy
    ham!(p)
    push!(hams, p.ham)
    push!(accDeltaH, p.ham)

    while frame <= nframes - 1
        if move!(p)
            nmoves += 1
            
            push!(accDeltaH, accDeltaH[end] + p.deltaH)
            
            ham!(p)
            push!(hams, p.ham)

            frame += 1
        end
    end

    # Use isapprox for floating point comparison
    @test all(isapprox.(hams, accDeltaH, atol=1e-10))
end

@testset "Hamiltonian (heterogeneous polymer with bulk)" begin
    l = Lattice(;
        boxlength = 20,
        nspecies = 4  # m1, m2, b, c
    )

    p = Polymer(;
        length = 10,
        maxlength = 30,
        nMonomerTypes = 2,
        initPos = [1,10,10],
        pattern = :custom,
        customPattern = [1,1,1,1,2,2,2,2,1,1]
    )
    
    initMonomers!(p, defaultChemPot=1.0, defaultIntPot=-0.3)
    
    m1 = p.monomers[1]
    m2 = p.monomers[2]

    b = FastBulkMonomer(;
        length = 64,
        initShape = "uniform"
    )

    c = FastBulkMonomer(;
        length = 64,
        initShape = "uniform"
    )

    # self-interaction potentials
    intPot!(m1, -0.3)  # base type
    intPot!(m2, -0.5)  # pattern type
    intPot!(b, -0.2)
    intPot!(c, -0.1)
    
    # polymer-bulk cross interactions (different for each monomer type)
    intPot!(m1, b, -0.1)
    intPot!(m1, c, -0.05)
    intPot!(m2, b, -0.4)
    intPot!(m2, c, -0.2)

    slist = [p, b, c]

    # initialize lattice and species
    init!(l, slist)

    #### Time evolution test - track each species separately
    nmoves = 0
    frame = 1
    nframes = 10000
    
    # Track hamiltonians per species
    p_hams = Vector{Float64}(undef, 0)
    p_accDeltaH = Vector{Float64}(undef, 0)
    b_hams = Vector{Float64}(undef, 0)
    b_accDeltaH = Vector{Float64}(undef, 0)
    c_hams = Vector{Float64}(undef, 0)
    c_accDeltaH = Vector{Float64}(undef, 0)
    
    # first frame, initial energy
    ham!(p)
    ham!(b)
    ham!(c)
    push!(p_hams, p.ham)
    push!(p_accDeltaH, p.ham)
    push!(b_hams, b.ham)
    push!(b_accDeltaH, b.ham)
    push!(c_hams, c.ham)
    push!(c_accDeltaH, c.ham)
    
    # Reset deltas
    p.deltaH = 0.0
    b.deltaH = 0.0
    c.deltaH = 0.0
    p.deltaHP = 0.0
    b.deltaHP = 0.0
    c.deltaHP = 0.0

    # List of movable species: polymer, b, c
    movable = [p, b, c]
    accDeltaHs = [p_accDeltaH, b_accDeltaH, c_accDeltaH]
    hamss = [p_hams, b_hams, c_hams]
    
    while frame <= nframes - 1
        chooser = rand(1:3)
        s = movable[chooser]
        
        if move!(s)
            nmoves += 1
            
            # Accumulate deltaH for mover
            push!(accDeltaHs[chooser], accDeltaHs[chooser][end] + getDeltaHam(s))
            
            # Accumulate deltaHP for partners
            for j in 1:3
                if j != chooser
                    push!(accDeltaHs[j], accDeltaHs[j][end] + getDeltaHP(movable[j]))
                end
            end
            
            # Compute fresh ham for all species
            for j in 1:3
                ham!(movable[j])
                push!(hamss[j], getHam(movable[j]))
            end

            frame += 1
        end
    end

    # Use isapprox for floating point comparison
    @test all(isapprox.(p_hams, p_accDeltaH, atol=1e-10))
    @test all(isapprox.(b_hams, b_accDeltaH, atol=1e-10))
    @test all(isapprox.(c_hams, c_accDeltaH, atol=1e-10))
end

@testset "Hamiltonian (BulkMonomer)" begin
    l = Lattice(;
        boxlength = 16,
        nspecies = 1
    )

    # BulkMonomer species
    p = BulkMonomer(;
        length = 500
    )

    slist = [p]

    # Set interaction potential after init
    intPot!(p, -1.0)

    # initialize lattice and species
    init!(l, slist)

    # initial energy
    hamInitial = p.ham

    #### iterating multiple moves
    nmoves = 0
    frame = 1
    nframes = 1000
    step = 50
    deltaH = 0
    moveMade = false
    
    hams = zeros(nframes)
    accDeltaH = zeros(nframes)
    
    # first frame
    accDeltaH[frame] = hamInitial
    
    ham!(p)
    hams[frame] = hamInitial

    while frame <= nframes - 1
        if moveMade
            deltaH += p.deltaH

            if mod(nmoves, step) == 0
                frame += 1

                accDeltaH[frame] = accDeltaH[frame-1] + deltaH
                deltaH = 0
                
                ham!(p)
                hams[frame] = p.ham
            end
        end
        
        moveMade = move!(p)
        nmoves += moveMade
    end

    @test hams == accDeltaH
end

@testset "Hamiltonian (fast BulkMonomer)" begin
    l = Lattice(;
        boxlength = 16,
        nspecies = 1
    )

    # BulkMonomer species
    p = FastBulkMonomer(;
        length = 500
    )

    slist = [p]

    # Set interaction potential after init
    intPot!(p, -1.0)

    # initialize lattice and species
    init!(l, slist)

    # initial energy
    hamInitial = p.ham

    #### iterating multiple moves
    nmoves = 0
    frame = 1
    nframes = 1000
    step = 50
    deltaH = 0
    moveMade = false
    
    hams = zeros(nframes)
    accDeltaH = zeros(nframes)
    
    # first frame
    accDeltaH[frame] = hamInitial
    
    ham!(p)
    hams[frame] = hamInitial

    while frame <= nframes - 1
        if moveMade
            deltaH += p.deltaH

            if mod(nmoves, step) == 0
                frame += 1

                accDeltaH[frame] = accDeltaH[frame-1] + deltaH
                deltaH = 0
                
                ham!(p)
                hams[frame] = p.ham
            end
        end
        
        moveMade = move!(p)
        nmoves += moveMade
    end
    
    @test hams == accDeltaH
end

@testset "Hamiltonian (multi-component BulkMonomer)" begin
    l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    a = BulkMonomer(;
        length = 128,
        initShape = "box",
        initParams = 1,
        initPos = [8,8,8]
    )

    b = BulkMonomer(;
        length = 128,
        initShape = "uniform"
    )

    c = BulkMonomer(;
        length = 128,
        initShape = "uniform"
    )

    # self and interaction potentials
    intPot!(a, -1.0)
    intPot!(b, -2.0)
    intPot!(c, -3.0)
    intPot!(a, b, -1.5)
    intPot!(a, c, -0.5)
    intPot!(b, c, -1.0)

    slist = [a, b, c]

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = [0,0,0]
    frame = [1,1,1]
    frame_joint = 1
    nframes = 10000
    step = 1
    deltaH = [0.0,0.0,0.0]
    niter = 0
    moveMade = [false, false, false]

    hamsA = Vector{Float64}(undef, 0)
    hamsB = Vector{Float64}(undef, 0)
    hamsC = Vector{Float64}(undef, 0)
    hams_ = [hamsA, hamsB, hamsC]
    hams = []
    accDeltaHA = Vector{Float64}(undef, 0)
    accDeltaHB = Vector{Float64}(undef, 0)
    accDeltaHC = Vector{Float64}(undef, 0)
    accDeltaH_ = [accDeltaHA, accDeltaHB, accDeltaHC]
    accDeltaH = []
    
    # first frame, initial energy
    push!(hams_[1], a.ham)
    push!(hams_[2], b.ham)
    push!(hams_[3], c.ham)
    push!(accDeltaH_[1], a.ham)
    push!(accDeltaH_[2], b.ham)
    push!(accDeltaH_[3], c.ham)

    for (sp, idx) in a.intSpecies
        push!(hams, hams_[idx])
        push!(accDeltaH, accDeltaH_[idx])
    end

    ham!(a)
    ham!(b)
    ham!(c)
    a.deltaH = 0.0
    b.deltaH = 0.0
    c.deltaH = 0.0

    while frame_joint <= nframes - 2
        for (s, i) in a.intSpecies
            if move!(s)
                nmoves[i] += 1
                deltaH[i] += s.deltaH

                # update deltaH for the partners
                for (p, j) in s.intSpecies
                    if j != s.intIdx
                        deltaH[j] += p.deltaHP
                    end
                end

                push!(accDeltaH[i], accDeltaH[i][end] + deltaH[i])
                deltaH[i] = 0
                
                ham!(s)
                push!(hams[i], s.ham)

                frame[i] += 1
                frame_joint += 1
            end
        end
    end

    @test hams[1] == accDeltaH[1]
    @test hams[2] == accDeltaH[2]
    @test hams[3] == accDeltaH[3]
end

@testset "Hamiltonian (multi-component fast BulkMonomer)" begin
    l = Lattice(;
        boxlength = 16,
        nspecies = 2
    )

    a = FastBulkMonomer(;
        length = 128,
        initShape = "box",
        initParams = 1,
        initPos = [8,8,8]
    )

    b = FastBulkMonomer(;
        length = 128,
        initShape = "uniform"
    )

    c = FastBulkMonomer(;
        length = 128,
        initShape = "uniform"
    )

    # self and interaction potentials
    intPot!(a, -1.0)
    intPot!(b, -2.0)
    intPot!(c, -3.0)
    intPot!(a, b, -1.5)
    intPot!(a, c, -0.5)

    slist = [a, b, c]

    # initialize lattice and species
    init!(l, slist)

    #### iterating multiple moves
    nmoves = [0,0,0]
    frame = [1,1,1]
    frame_joint = 1
    nframes = 10000
    step = 1
    deltaH = [0.0,0.0,0.0]
    niter = 0
    moveMade = [false, false, false]

    hamsA = Vector{Float64}(undef, 0)
    hamsB = Vector{Float64}(undef, 0)
    hamsC = Vector{Float64}(undef, 0)
    hams_ = [hamsA, hamsB, hamsC]
    hams = []
    accDeltaHA = Vector{Float64}(undef, 0)
    accDeltaHB = Vector{Float64}(undef, 0)
    accDeltaHC = Vector{Float64}(undef, 0)
    accDeltaH_ = [accDeltaHA, accDeltaHB, accDeltaHC]
    accDeltaH = []
    
    # first frame, initial energy
    push!(hams_[1], a.ham)
    push!(hams_[2], b.ham)
    push!(hams_[3], c.ham)
    push!(accDeltaH_[1], a.ham)
    push!(accDeltaH_[2], b.ham)
    push!(accDeltaH_[3], c.ham)

    for (sp, idx) in a.intSpecies
        push!(hams, hams_[idx])
        push!(accDeltaH, accDeltaH_[idx])
    end

    ham!(a)
    ham!(b)
    ham!(c)
    a.deltaH = 0.0
    b.deltaH = 0.0
    c.deltaH = 0.0

    while frame_joint <= nframes - 2
        for (s, i) in a.intSpecies
            if move!(s)
                nmoves[i] += 1
                deltaH[i] += s.deltaH

                # update deltaH for the partners
                for (p, j) in s.intSpecies
                    if j != s.intIdx
                        deltaH[j] += p.deltaHP
                    end
                end

                push!(accDeltaH[i], accDeltaH[i][end] + deltaH[i])
                deltaH[i] = 0
                
                ham!(s)
                push!(hams[i], s.ham)

                frame[i] += 1
                frame_joint += 1
            end
        end
    end

    @test hams[1] == accDeltaH[1]
    @test hams[2] == accDeltaH[2]
    @test hams[3] == accDeltaH[3]
end

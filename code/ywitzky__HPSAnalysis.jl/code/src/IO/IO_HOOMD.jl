function getGSDTrajectoryFiles(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    FolderPath="/"*joinpath(split(Sim.TrajectoryFile,"/")[1:end-1])*"/"
    return ["$(FolderPath)$file" for file in readdir(FolderPath) if length(file)>4 && file[1:4]=="traj" && file[end-3:end]==".gsd"]
end

function read_HOOMD_GSD!(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    N = Sim.NAtoms
    image_old = zeros(R, (N, 3))

    images_not_saved = zeros(Bool, Sim.NSteps)# false
    step=1
    for file in getGSDTrajectoryFiles(Sim)
        gsdfileobj = GSDFormat.open(file, "r")

        for gsd_step in 2:length(gsdfileobj) ### eachindex(gsdfileobj); skip the first index since we enabled write_at_start
            xyz   = gsdfileobj[gsd_step].particles.position

            Sim.x[1:N,step] .=  xyz[1:N,1] .*10.0
            Sim.y[1:N,step] .=  xyz[1:N,2] .*10.0
            Sim.z[1:N,step] .=  xyz[1:N,3] .*10.0

            image =  gsdfileobj[gsd_step].particles.image 
            Sim.x_uw[1:N,step] .= Sim.x[1:N,step] .+ Sim.BoxLength[1].*(image[1:N,1])
            Sim.y_uw[1:N,step] .= Sim.y[1:N,step] .+ Sim.BoxLength[2].*image[1:N,2]
            Sim.z_uw[1:N,step] .= Sim.z[1:N,step] .+ Sim.BoxLength[3].*(image[1:N,3])
            image_old .= image
            step += 1
        end
    end

    return nothing
end

function read_HOOMD_Particles(Sim::SimData{R,I}, File::String) where {R<:Real, I<:Integer}
    lines = readlines(File)
    Sim.NAtoms = length(lines)-1
    Sim.NBonds = Sim.NAtoms - Sim.NChains
    Sim.NAngles = Sim.NAtoms - 2*Sim.NChains
    Sim.NDihedrals = Sim.NAtoms - 2*Sim.NChains

    Sim.IDs = zeros(I, Sim.NAtoms)
    Sim.Charges = zeros(R, Sim.NAtoms)
    Sim.Masses = zeros(R, Sim.NAtoms)

    for line in lines[2:end]
        split = Base.split(line, ",")
        index  = parse(I, split[1])
        id     = parse(I, split[2]) +1 ### +1 to be able to use id as index, and be consistent with dictionary
        charge = parse(R, split[6])
        mass   = parse(R, split[7])
        Sim.IDs[index] = id
        Sim.Charges[index] = charge
        Sim.Masses[index] = mass
    end

    Sim.ChainMasses = [sum(Sim.Masses[start:stop]) for (start,stop) in zip(Sim.ChainStart, Sim.ChainStop)  ]
    return nothing
end

function read_HOOMD_Sequences(Sim::SimData{R,I}, File::String) where {R<:Real, I<:Integer}
    lines = readlines(File)
    Sim.NChains = length(lines)
    Sim.Sequences = Vector{String}()
    for (id,line) in enumerate(lines)
        push!(Sim.Sequences, line)
    end

    Sim.ChainStart = zeros(I, Sim.NChains)
    Sim.ChainStop  = zeros(I, Sim.NChains)
    Sim.ChainStart[1] = 1
    Sim.ChainStop[1] = length(Sim.Sequences[1])
    for (i,Seq) in enumerate(Sim.Sequences)
        if i ==1 
            continue
        end
        Sim.ChainStart[i] = Sim.ChainStart[i-1]+length(Seq)
        Sim.ChainStop[i] = Sim.ChainStop[i-1]+length(Seq)
    end

    Sim.ChainLength = Sim.ChainStop .- Sim.ChainStart .+ 1
    Sim.MaxChainLength = maximum(Sim.ChainLength)
    return nothing
end

function read_HOOMD_Dictionaries(Sim::SimData{R,I}, File::String) where {R<:Real, I<:Integer}
    lines = readlines(File)
    Sim.IDToMasses  = zeros(R, length(lines)-1)
    Sim.IDToSigmas  = zeros(R, length(lines)-1)
    Sim.IDToLambdas = zeros(R, length(lines)-1)
    for line in lines[2:end]
        split = Base.split(line, ",")
        id       = parse(I, split[1])
        Sim.IDToResName[id] = strip(split[2])
        Sim.IDToMasses[id]  = parse(R, split[4])
        Sim.IDToSigmas[id]  = parse(R, split[5])
        Sim.IDToLambdas[id] = parse(R, split[6])
    end
    Sim.NAtomTypes= length(lines)-1
    Sim.NBondTypes=1
    Sim.NDihedralTypes=0
    return nothing
end

function read_HOOMD_Params(Sim::SimData{R,I}, File::String) where {R<:Real, I<:Integer}
    file = open(File, "r")
    Sim.BoxLength = zeros(R, 3) 
    Sim.BoxSize = zeros(R, 3,2 )
    while ! eof(file)
        line = readline(file)  
        split = Base.split(line)
        if split[1] =="Simname:"
            Sim.SimulationName = split[2]
        elseif split[1] == "Temp:"
            Sim.TargetTemp = parse(R, split[2])
        elseif split[1] == "NOut:"
            nout =  parse(I, split[2])
            Sim.EnergyWriteOutFreq = nout
            Sim.TrajWriteOutFreq = nout
            Sim.Reduce= Sim.EnergyWriteOutFreq//Sim.TrajWriteOutFreq
        elseif split[1] == "Lx:"
            Sim.BoxLength[1] =  parse(R, split[2])*10.0
        elseif split[1] == "Ly:"
            Sim.BoxLength[2] =  parse(R, split[2])*10.0
        elseif split[1] == "Lz:"
            Sim.BoxLength[3] =  parse(R, split[2])*10.0
        elseif split[1] == "Trajectory:"
            Sim.TrajectoryFile= Sim.BasePath*split[2]
            Sim.EnergyFile= Sim.BasePath*split[2]
        end
    end
    Sim.BoxSize[1,1] = -Sim.BoxLength[1]/2.
    Sim.BoxSize[1,2] =  Sim.BoxLength[1]/2.
    Sim.BoxSize[2,1] = -Sim.BoxLength[2]/2.
    Sim.BoxSize[2,2] =  Sim.BoxLength[2]/2.
    Sim.BoxSize[3,1] = -Sim.BoxLength[3]/2.
    Sim.BoxSize[3,2] =  Sim.BoxLength[3]/2.
    return nothing
end

function read_HOOMD_Folder(Sim::SimData{R,I}, Folder::String) where {R<:Real, I<:Integer}
    read_HOOMD_Params(Sim, Folder*"Params.txt")
    read_HOOMD_Dictionaries(Sim, Folder*"Dictionaries.txt")
    read_HOOMD_Sequences(Sim, Folder*"Sequences.txt")
    read_HOOMD_Particles(Sim, Folder*"Particles.txt")
    return nothing
end
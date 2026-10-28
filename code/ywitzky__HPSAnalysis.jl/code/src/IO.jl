#using HPSAnalysis.Polyply
using JLD2, Mmap, LoopVectorization ,  Printf, HDF5, GSDFormat

include("./IO/Xtc.jl")
include("./IO/IO_HOOMD.jl")
include("./IO/Sequence_IO.jl")
include("./IO/Formats.jl")

function parseXYZ!(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    traj=open(Sim.TrajectoryFile, "r")
    
    xnew  = zeros(R, Sim.NAtoms)
    ynew  = zeros(R, Sim.NAtoms)
    znew  = zeros(R, Sim.NAtoms)

    step=0
    atom_cnt=0
    write_step=0
    for line in eachline(traj)    
        nspaces=count(" ", line)
        if (step%Sim.StepFrequency==0)
            if atom_cnt==0 ### Atoms. Timestep: xxxxx
                if write_step!=0
                    Sim.x[:, write_step] .= xnew
                    Sim.y[:, write_step] .= ynew
                    Sim.z[:, write_step] .= znew
                end
                write_step+=1
                atom_cnt+=1
                continue
            elseif nspaces==3 
                lines = Base.split(line)
                if lines[1]=="Atoms."
                    continue
                end                
                split = parse.(R,lines)

                xnew[atom_cnt] = split[2]
                ynew[atom_cnt] = split[3]
                znew[atom_cnt] = split[4]
    
                atom_cnt+=1
            end
        end
        if nspaces==0
            atom_cnt=0
            step+=1
        end
        
    end
    close(traj)

    Sim.x = Sim.x[:,1:write_step-1 ]
    Sim.y = Sim.y[:,1:write_step-1 ]
    Sim.z = Sim.z[:,1:write_step-1 ]

    return write_step-1
end

function replaceE(String)
    if count("e", String)>0
        String="0.0"
    end
    return String
end

function parseXYZ_CSV!(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    traj=open(Sim.TrajectoryFile, "r")
    
    xnew  = zeros(eltype(Sim.x), Sim.NAtoms)
    ynew  = zeros(eltype(Sim.x), Sim.NAtoms)
    znew  = zeros(eltype(Sim.x), Sim.NAtoms)

    Sim.FrameWeights = zeros(eltype(Sim.x), Sim.NSteps)

    step=0
    atom_cnt=0
    write_step=0
    two_lines=false
    for line in eachline(traj)    
        nspaces=count(",", line)

        if two_lines
            x = tryparse(eltype(Sim.x), line) ### avoid weights that are so small that they cant be parsed
            write_step+=1
            Sim.FrameWeights[write_step] = x!== nothing ? x : 0.0
            if write_step!=1
                Sim.x[:, write_step-1] .= xnew
                Sim.y[:, write_step-1] .= ynew
                Sim.z[:, write_step-1] .= znew
            end
            atom_cnt=1
        end
        if nspaces==0
            two_lines = ~two_lines
        end
        if (atom_cnt==0 && two_lines) ### Atoms. Timestep: xxxxx
            continue
        elseif nspaces==2
            lines = Base.split(line,",")
            if lines[1]=="Atoms."
                continue
            end

            lines = replaceE.(lines)

            split = parse.(eltype(Sim.x),lines)
            xnew[atom_cnt] = split[1]
            ynew[atom_cnt] = split[2]
            znew[atom_cnt] = split[3]

            atom_cnt+=1
        end        
    end
    close(traj)

    return write_step-1
end

function readEnergyFile(Sim::SimData{R,I};EnergyFile::String, NumSteps=-1) where {R<:Real, I<:Integer}

    if Sim.HOOMD return end
    if NumSteps==-1
        step_cnt=1
        hash_cnt=0

        energy  = open(EnergyFile)
        NEntries = countlines(energy) -1 ### reads till the end
        close(energy)

        energy  = open(EnergyFile)

        line = readline(energy)
        while (line[1]=="#"[1])
            line = readline(energy)
        end
        NCols = count(" ", line)+1
        Sim.Energies = zeros(eltype(Sim.x), (NEntries, NCols))
        Sim.Energies[1,:] .= parse.(eltype(Sim.x),split(line))
        while ~eof(energy) 
            line = readline(energy)

            if(line[1]!="# 1"[1])
                step_cnt+=1
                Sim.Energies[step_cnt,:] .= parse.(eltype(Sim.x),split(line))
            else
                hash_cnt+=1
            end
        end

        println("energy file step: $(step_cnt) , number of restarts: $(hash_cnt) vs. $(Int32(floor(NEntries*Sim.Reduce))) , $(Sim.Reduce)")
        Sim.NRestarts=hash_cnt
        Sim.NSteps=step_cnt+1
        Sim.NSteps=Int32(floor(NEntries*Sim.Reduce)+hash_cnt+1000)#+hash_cnt+3000)
    else
        Sim.NSteps=NumSteps
    end
end

function parseXTC!(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    (stat, traj) = Xtc.xtc_init(Sim.TrajectoryFile)
    
    if Sim.NAtoms != traj.natoms 
        printstyled("Number of atoms in xtc does not match lammps input files.\n", color=:red)
    end

    Sim.x[:, 1] .= traj.xyz[1,:]
    Sim.y[:, 1] .= traj.xyz[2,:]
    Sim.z[:, 1] .= traj.xyz[3,:]
    prev_step=-1
    for step in 1:Sim.NSteps
        stat= Xtc.read_xtc(traj) ### xtc works in nm instead of angstroem
        Sim.x[:, step] .= traj.xyz[1,:].*10.0
        Sim.y[:, step] .= traj.xyz[2,:].*10.0
        Sim.z[:, step] .= traj.xyz[3,:].*10.0
        if prev_step==traj.step[1]
            Xtc.close_xtc(traj)
            return step
        end
        prev_step = traj.step[1]
    end
    Xtc.close_xtc(traj)
end

function readH5MD!(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    start = 1
    for rest in 0:(Sim.NRestarts)
        file = replaceLammpsVariables(Sim.TrajectoryFile, Dict("NRes"=>rest))
        c = h5open(file, "r") do file
            data = read(file, "particles")
            N= Int32(length(axes(data["all"]["position"]["value"],3)))
            stop = start + N-1
            Sim.x[:,start:stop] .= data["all"]["position"]["value"][1,:,:] 
            Sim.y[:,start:stop] .= data["all"]["position"]["value"][2,:,:] 
            Sim.z[:,start:stop] .= data["all"]["position"]["value"][3,:,:] 

            Sim.x_uw[:,start:stop] .= Sim.x[:,start:stop] .+ Sim.BoxLength[1].*data["all"]["image"]["value"][1,:,:] 
            Sim.y_uw[:,start:stop] .= Sim.y[:,start:stop] .+ Sim.BoxLength[2].*data["all"]["image"]["value"][2,:,:] 
            Sim.z_uw[:,start:stop] .= Sim.z[:,start:stop] .+ Sim.BoxLength[3].*data["all"]["image"]["value"][3,:,:] 

            start += N
        end
    end

    return size(Sim.x,2)
end

function CountStepsInGSDTrajectoryFiles(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    N=0
    for file in getGSDTrajectoryFiles(Sim)
        gsdfileobj = GSDFormat.open_gsd(file, "r";application="gsd.hoomd ", schema="hoomd", schema_version=(1, 4))
        N += convert(Int32, GSDFormat.get_nframes(gsdfileobj)) -1 ### hoomd is set to write at the start which we will throw away to not double count if we stich trajetories together
    end
    return N
end


function readXYZ!(Sim::SimData{T,I}; TrajectoryFile::String, EnergyFile::String, Minimize=false, NumSteps=-1, Delimiter=" ")  where {T<:Real,I<:Integer}
    readEnergyFile(Sim; EnergyFile=EnergyFile, NumSteps=NumSteps)
    N = Sim.NSteps

    ### Sim.StepFrequency is for data that will be reduceed, Sim.reduce sets discrepancy between energy data and .xyz created by presorting
    if Sim.TrajectoryFile[end-2:end] =="xyz"
        traj=open(TrajectoryFile, "r")
        Sim.NAtoms= parse(I,readline(traj))
    elseif Sim.TrajectoryFile[end-2:end] =="xtc"
        (_, xtc) = xtc_init(TrajectoryFile)
        Sim.NAtoms = xtc.natoms
        close(xtc)
    elseif Sim.TrajectoryFile[end-1:end] =="h5"
        N=I(0)
        for rest in 0:Sim.NRestarts
            c = h5open(replaceLammpsVariables(Sim.TrajectoryFile, Dict("NRes"=>rest)), "r") do file
                content = read(file, "particles")
                N += I(length(axes(content["all"]["position"]["value"],3)))
            end
        end
    elseif Sim.TrajectoryFile[end-2:end] =="gsd"
        Sim.NSteps= CountStepsInGSDTrajectoryFiles(Sim)
        N= Sim.NSteps
    else
        error("Can find TrajectoryFile:\"$TrajectoryFile\".")
    end

    Sim.xio= open(Sim.xFilePath,"w+")
    Sim.yio= open(Sim.yFilePath,"w+")
    Sim.zio= open(Sim.zFilePath,"w+")

    ### ask for more in case write out was incomplete for energy files
    Sim.x =  Mmap.mmap(Sim.xio, Matrix{T}, (Sim.NAtoms,N))
    Sim.y =  Mmap.mmap(Sim.yio, Matrix{T}, (Sim.NAtoms,N))
    Sim.z =  Mmap.mmap(Sim.zio, Matrix{T}, (Sim.NAtoms,N))

    ### actually read the data
    if Sim.TrajectoryFile[end-3:end] ==".xyz"
        if Delimiter==","
            Sim.NSteps= parseXYZ_CSV!(Sim) 
        else
            Sim.NSteps= parseXYZ!(Sim) 
        end
    elseif Sim.TrajectoryFile[end-3:end] ==".xtc"
        parseXTC!(Sim)

    elseif Sim.TrajectoryFile[end-2:end] ==".h5" || Sim.TrajectoryFile[end-3:end] ==".gsd"
        ### assume h5md format but could technically be any different hdf5 format
        Sim.x_uw_io= open(Sim.x_uw_FilePath,"w+")
        Sim.y_uw_io= open(Sim.y_uw_FilePath,"w+")
        Sim.z_uw_io= open(Sim.z_uw_FilePath,"w+")
    
        Sim.x_uw =  Mmap.mmap(Sim.x_uw_io, Matrix{T}, (Sim.NAtoms,N))
        Sim.y_uw =  Mmap.mmap(Sim.y_uw_io, Matrix{T}, (Sim.NAtoms,N))
        Sim.z_uw =  Mmap.mmap(Sim.z_uw_io, Matrix{T}, (Sim.NAtoms,N))

        if Sim.TrajectoryFile[end-2:end] ==".h5"
            Sim.NSteps=readH5MD!(Sim)
        elseif  Sim.TrajectoryFile[end-3:end] ==".gsd"
            read_HOOMD_GSD!(Sim)
        end

        Mmap.sync!(Sim.x)
        Mmap.sync!(Sim.y)
        Mmap.sync!(Sim.z)
        close(Sim.xio)
        close(Sim.yio)
        close(Sim.zio)
    else
        println("Can not read file format: $( Sim.TrajectoryFile[end-3:end])")
    end

    ### sync RAM to disk before closing
    Mmap.sync!(Sim.x)
    Mmap.sync!(Sim.y)
    Mmap.sync!(Sim.z)
    pos = zeros(length(Sim.x[:,1]), 3)
    pos[:, 1] = Sim.x[:, 1]
    pos[:, 2] = Sim.y[:, 1]
    pos[:, 3] = Sim.z[:, 1]
    close(Sim.xio)
    close(Sim.yio)
    close(Sim.zio)
    return pos
end

function replaceLammpsVariables(str::S, dict::Dict{String, A}) where {S<:AbstractString, A<:Any}
    for key in keys(dict)
        if (typeof(dict[key]) <:  AbstractString)
            str = replace(str, "\${$(key)}"=>String(dict[key]))
        else
            str = replace(str, "\${$(key)}"=>"$(dict[key])")
        end
    end
    return str
end

function readLMP!(Sim::SimData{R,I}, LMPFile::String) where {R<:Real, I<:Integer}
    file = open(LMPFile, "r")
    while ! eof(file)
        line = readline(file)  
        split = Base.split(line)
        if length(split)==0 continue end
        if( split[1]=="#")
            continue
        elseif( split[1]=="dump") ### works only for xyz so far
            tmp = Base.split(split[6],"./")
            Sim.TrajectoryFile = Sim.BasePath*(replaceLammpsVariables(tmp[2], Sim.LammpsVariables))
            Sim.TrajWriteOutFreq = parse(Int32, split[5])
        elseif( split[1]=="variable") ### variable assignment in lammps , mostly for REMD/Charge-HREMD 
            if split[3]=="world"
                continue ### parsed in HREMD lmp read
            else
                if split[2] !="NRes"
                Sim.LammpsVariables[split[2]] = replaceLammpsVariables(strip(split[4], ['"']), Sim.LammpsVariables)
                end
            end
        elseif (split[1]=="fix")
            if (split[4] == "langevin" || split[4] == "langevin/omp"|| split[4] == "langevin/kk")
                if split[5]!=split[6]
                    printstyled("Start and End Temperatur are not identical." ; color=:red)
                else
                    Sim.TargetTemp = parse(Float64, split[5])
                end
            elseif ( split[4]=="print")
                tmp = Base.split(replaceLammpsVariables(split[20], Sim.LammpsVariables),"./")
                Sim.EnergyWriteOutFreq = parse(eltype(Sim.EnergyWriteOutFreq), split[5])

                Sim.EnergyDict= Dict{String, eltype(Sim.NReplica)}([(strip(val,['"', '\\', '(', ')' , '$']), id) for (id, val) in enumerate(split[6:18])])
                Sim.EnergyFile = Sim.BasePath*(tmp[end])
            end
        elseif split[1]=="read_data"
            tmp = Base.split(split[2],"./")
            Sim.StartFile = Sim.BasePath*( replaceLammpsVariables(split[2], Sim.LammpsVariables))  
        end
    end
    Sim.Reduce= Sim.EnergyWriteOutFreq//Sim.TrajWriteOutFreq
end

function readMasses!(Sim::SimData{R,I}, file::IOStream) where {R<:Real, I<:Integer}
    Sim.IDToMasses=zeros(Float32, Sim.NAtomTypes)
    for i in 1:(Sim.NAtomTypes+1)
        words = split(readline(file))
        if length(words)==4 || length(words)==5
            ID = parse(Int32, words[1])
            Sim.IDToResName[ID]=strip(words[4], ['(', '\'', ' ', ])
            if length(Sim.IDToResName[ID])>1
                Sim.IDToResName[ID] = lowercase(strip(Sim.IDToResName[ID], ['\'',',']))
            end
            Sim.IDToMasses[ID]=parse(Float32, words[2])
        end
    end
end

function readAtoms!(Sim::SimData{R,I}, file::IOStream) where {R<:Real, I<:Integer}
    Sim.Sequences=[""]
    Sim.ChainStart=zeros(Sim.NChains)
    Sim.ChainStop=zeros(Sim.NChains)

    Sim.IDs=zeros(Float32, Sim.NAtoms)
    Sim.Charges=zeros(Float32, Sim.NAtoms)
    Sim.Masses=zeros(Float32, Sim.NAtoms)

    curChain=1
    Sim.ChainStart[1]=1
    for i in 1:(Sim.NAtoms+1)
        words = split(readline(file))
        if length(words)>0
            if words[1]=="#"
                continue
            end
        end
        if length(words)>=4
            ID = parse(Int32, words[1])
            atom_ID=parse(Int32, words[3]) 
            Sim.IDs[ID]=atom_ID
            Sim.Masses[ID]= Sim.IDToMasses[atom_ID]
            Sim.Charges[ID]=parse(Float32, words[4])
            ChainID = parse(Int32, words[2])
            if curChain != ChainID
                Sim.ChainStop[curChain]=i-2
                curChain= ChainID
                Sim.ChainStart[ChainID]=i-1
                push!(Sim.Sequences, "") ### appends at the end
            end
            Sim.Sequences[curChain] *= Sim.IDToResName[atom_ID]
        end
    end

    Sim.ChainStop[Sim.NChains]=Sim.NAtoms
    Sim.ChainLength = Sim.ChainStop .- Sim.ChainStart .+ 1
    Sim.MaxChainLength = maximum(Sim.ChainLength)

    ### precompute masses of a chain
    Sim.ChainMasses=zeros(Float32, Sim.NChains)
    for chain in 1:Sim.NChains
        #=
        ### first/last amonio acids dont have two amino bonds
        Sim.Masses[Sim.ChainStart[chain]] += 2.0
        Sim.Masses[Sim.ChainStop[chain]] += 16.0
        =#
        for ind in Sim.ChainStart[chain]:Sim.ChainStop[chain]
            Sim.ChainMasses[chain]+= Sim.Masses[ind]
        end
        
    end
end

function readBonds(Sim::SimData{R,I}, file::IOStream) where {R<:Real, I<:Integer}
    for i in 1:(Sim.NBonds+1)
        readline(file)
    end
end

function readDihedrals(Sim::SimData{R,I}, file::IOStream) where {R<:Real, I<:Integer}
    for i in 1:(Sim.NDihedrals+1)
        readline(file)
    end
end

function readAngles(Sim::SimData{R,I}, file::IOStream) where {R<:Real, I<:Integer}
    for i in 1:(Sim.NAngles+1)
        readline(file)
    end
end

function readStartFile!(Sim::SimData{R,I}, StartFile::String) where {R<:Real, I<:Integer}
    file = open(StartFile, "r")

    cnt = 1
    while ! eof(file) 
        words = split(readline(file))
        if length(words)==1
            if words[1]=="Masses"
                readMasses!(Sim, file)
            elseif words[1]=="Atoms"
                readAtoms!(Sim, file)
            elseif words[1]=="Bonds"
                readBonds(Sim, file) ### so far just reads the lines without saving
            elseif words[1]=="Angles"
                readAngles(Sim, file) ### so far just reads the lines without saving
            elseif words[1]=="Dihedrals"
                readDihedrals(Sim, file) ### so far just reads the lines without saving
            end
        elseif length(words)==2
            if words[2]=="atoms"
                Sim.NAtoms=parse(Int32, words[1])
            elseif words[2]=="bonds"
                Sim.NBonds=parse(Int32, words[1])
                Sim.NChains = Sim.NAtoms-Sim.NBonds
            elseif words[2]=="angles"
                Sim.NAngles=parse(Int32, words[1])
            elseif words[2]=="dihedrals"
                Sim.NDihedrals=parse(Int32, words[1])
            end
        elseif length(words)==3
            if words[2]=="atom" && words[3]=="types"
                Sim.NAtomTypes=parse(Int32, words[1])
            elseif words[2]=="bond" && words[3]=="types"
                Sim.NBondTypes=parse(Int32, words[1])
            elseif words[2]=="angle"&& words[3]=="types"
                Sim.NAngleTypes=parse(Int32, words[1])
            elseif words[2]=="dihedral"&& words[3]=="types"
                Sim.NDihedralTypes=parse(Int32, words[1])
            end
        elseif length(words)==4
            if words[3]=="xlo" && words[4]=="xhi"
                Sim.BoxSize[1,1] = parse(Float32, words[1])
                Sim.BoxSize[1,2] = parse(Float32, words[2])
                Sim.BoxLength[1] = Sim.BoxSize[1,2]-Sim.BoxSize[1,1]
            elseif words[3]=="ylo" && words[4]=="yhi"
                Sim.BoxSize[2,1] = parse(Float32, words[1])
                Sim.BoxSize[2,2] = parse(Float32, words[2])
                Sim.BoxLength[2] = Sim.BoxSize[2,2]-Sim.BoxSize[2,1]
            elseif words[3]=="zlo" && words[4]=="zhi"
                Sim.BoxSize[3,1] = parse(Float32, words[1])
                Sim.BoxSize[3,2] = parse(Float32, words[2])
                Sim.BoxLength[3] = Sim.BoxSize[3,2]-Sim.BoxSize[3,1]
            end
        end
    end
end

function initReducedData(Path::String,XYZFile::String, NumSteps::Int,NChains::Int, ChainStart::Vector{Int},ChainStop::Vector{Int};Reparse=true, LoadAll=true)
    Data = SimData()
    Data.BasePath= Path
    Data.PlotPath=Data.BasePath*"/Plots/"
    Data.DataPath=Data.BasePath*"Data/"
    Data.xFilePath = Data.DataPath*"x.bin"
    Data.yFilePath = Data.DataPath*"y.bin"
    Data.zFilePath = Data.DataPath*"z.bin"
    Data.x_uw_FilePath = Data.DataPath*"x.bin"
    Data.y_uw_FilePath = Data.DataPath*"y.bin"
    Data.z_uw_FilePath = Data.DataPath*"z.bin"
    Data.TrajectoryFile=Path*XYZFile
    mkpath(Data.PlotPath)
    mkpath(Data.DataPath)

    Data.NChains=NChains
    Data.ChainStart =ChainStart
    Data.ChainStop=ChainStop

    Data.ChainLength = Data.ChainStop .- Data.ChainStart .+ 1
    Data.NAtoms  = sum(Data.ChainLength)
    Data.IDs = ones(eltype(Data.NSteps), Data.NAtoms)


    if LoadAll
        if Reparse
            readXYZ!(Data; TrajectoryFile=Data.TrajectoryFile, EnergyFile="", NumSteps=NumSteps, Delimiter=",")
        else
            ReadAllData!(Data)
        end
    end

    ### x,y,z are stored on disk and lazyly synchronised, ONLY IN READ MODE!!!
    Data.xio= open(Data.xFilePath,"r+")
    Data.yio= open(Data.yFilePath,"r+")
    Data.zio= open(Data.zFilePath,"r+")

    Data.x =  Mmap.mmap(Data.xio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.y =  Mmap.mmap(Data.yio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.z =  Mmap.mmap(Data.zio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    
    Data.x_uw_io= open(Data.xFilePath,"r+")
    Data.y_uw_io= open(Data.yFilePath,"r+")
    Data.z_uw_io= open(Data.zFilePath,"r+")

    Data.x_uw =  Mmap.mmap(Data.x_uw_io, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.y_uw =  Mmap.mmap(Data.y_uw_io, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.z_uw =  Mmap.mmap(Data.z_uw_io, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))

    return Data

end

function initReducedData_NoConfStorage(Path::String, NumSteps::Int,NChains::Int, ChainStart::Vector{Int},ChainStop::Vector{Int}) 
    Data = SimData()
    Data.BasePath= Path
    Data.PlotPath=Data.BasePath*"/Plots/"
    Data.DataPath=Data.BasePath*"Data/"
    Data.xFilePath = Data.DataPath*""
    Data.yFilePath = Data.DataPath*""
    Data.zFilePath = Data.DataPath*""
    Data.x_uw_FilePath = Data.DataPath*""
    Data.y_uw_FilePath = Data.DataPath*""
    Data.z_uw_FilePath = Data.DataPath*""
    Data.TrajectoryFile=Path*""
    mkpath(Data.PlotPath)
    mkpath(Data.DataPath)

    Data.NChains=NChains
    Data.ChainStart =ChainStart
    Data.ChainStop=ChainStop
    Data.NSteps = NumSteps

    Data.ChainLength = Data.ChainStop .- Data.ChainStart .+ 1
    Data.NAtoms  = sum(Data.ChainLength)
    Data.Masses = ones(eltype(Data.x), Data.NAtoms)
    Data.ChainMasses = [sum(Data.Masses[start:stop]) for (start,stop) in zip(ChainStart, ChainStop)  ]

    Data.IDs = ones(eltype(Data.NSteps), Data.NAtoms)

    Data.x =  zeros( eltype(Data.x), (Data.NAtoms,Data.NSteps))
    Data.y =  zeros( eltype(Data.x), (Data.NAtoms,Data.NSteps))
    Data.z =  zeros( eltype(Data.x), (Data.NAtoms,Data.NSteps))

    Data.x_uw = Data.x
    Data.y_uw = Data.y
    Data.z_uw = Data.z

    return Data

end

#read the defines Parameters and coordinates
function initData(Path::String;LmpName="", StepFrequency=1, Reparse=true, LoadAll=true, Reduce=1, EquilibrationTime=1, LammpsVariables= Dict{String,Any}(), BasePathAdd="", HOOMD=false, ReadBig=false) #where {R<:Real,I<:Integer}
    
    #Define path where the Data should be
    Data = SimData()#::SimData{R,I}
    Data.HOOMD=HOOMD
    Data.BasePath= Path
    Data.PlotPath=Data.BasePath*"/Plots/"*BasePathAdd
    Data.DataPath=Data.BasePath*"/Data/"*BasePathAdd
    Data.xFilePath = Data.DataPath*"x.bin"
    Data.yFilePath = Data.DataPath*"y.bin"
    Data.zFilePath = Data.DataPath*"z.bin"
    Data.x_uw_FilePath = Data.DataPath*"x_uw.bin"
    Data.y_uw_FilePath = Data.DataPath*"y_uw.bin"
    Data.z_uw_FilePath = Data.DataPath*"z_uw.bin"
    Data.Reduce=Reduce ### only a factor for prereduced data.
    mkpath(Data.PlotPath)
    mkpath(Data.DataPath)
    Data.LammpsVariables = LammpsVariables


    Data.StepFrequency=StepFrequency
    #Read HoomdSetup Datas
    if Data.HOOMD
        read_HOOMD_Folder(Data,Data.BasePath*"/HOOMD_Setup/")
    else
        Data.SimulationName=split(LmpName,".")[1]
        readLMP!(Data, Path*LmpName)
        readStartFile!(Data, Data.StartFile)
    end

    if LoadAll
        if Reparse
            readXYZ!(Data; TrajectoryFile=Data.TrajectoryFile, EnergyFile=Data.EnergyFile)
            #unfoldPositions(Data)
        else
            ReadAllData!(Data;ReadBig=ReadBig)
        end
    end
    if Data.NSteps<EquilibrationTime
        if LoadAll ### is this the right choice?
            printstyled("Choosen EquilibrationTime: $(EquilibrationTime) is below the number of steps: $(Data.NSteps) in the simulation. Computations might not compute results.\n";color=:red)
        end
        Data.EquilibrationTime=Data.NSteps-1
    else
        Data.EquilibrationTime=EquilibrationTime
    end
    if ~Data.HOOMD
        readStartFile!(Data, Data.StartFile)
    end

    ### x,y,z are stored on disk and lazyly synchronised, ONLY IN READ MODE!!!
    Data.xio= open(Data.xFilePath,"r+")
    Data.yio= open(Data.yFilePath,"r+")
    Data.zio= open(Data.zFilePath,"r+")

    Data.x =  Mmap.mmap(Data.xio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.y =  Mmap.mmap(Data.yio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.z =  Mmap.mmap(Data.zio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    
    if Reparse && Data.TrajectoryFile[end-2:end] !="gsd"
        unfoldPositions(Data)
    end

    Data.x_uw_io= open(Data.x_uw_FilePath,"r+")
    Data.y_uw_io= open(Data.y_uw_FilePath,"r+")
    Data.z_uw_io= open(Data.z_uw_FilePath,"r+")

    Data.x_uw =  Mmap.mmap(Data.x_uw_io, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.y_uw =  Mmap.mmap(Data.y_uw_io, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.z_uw =  Mmap.mmap(Data.z_uw_io, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))

    return Data
end

@doc raw"""
    SaveAllData(Sim::SimData{R,I}) where {R<:Real, I<:Integer}

Save of the Field information as a jld2 file.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Creat**:
* Save File.
"""
function SaveAllData(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Dictionary= Dict()
    Sim.BigDataList = [] #Array{Symbol}()
    for fieldname in fieldnames(typeof(Sim))
        if String(fieldname)=="x" || String(fieldname)=="y" ||  String(fieldname)=="z" || String(fieldname)=="x_ci"|| String(fieldname)=="y_ci"|| String(fieldname)=="z_ci" || String(fieldname)=="x_uw"|| String(fieldname)=="y_uw"|| String(fieldname)=="z_uw" 
            continue
        end
        if Base.summarysize(getfield(Sim,fieldname))>50*1024^2 ### more than 50Mb
            push!(Sim.BigDataList,fieldname )
            continue
        end

        Dictionary[String(fieldname)]= getfield(Sim,fieldname)
    end
    SaveFields(Sim,"All",Dictionary)

    for fields in Sim.BigDataList
        JLD2.save(Sim.DataPath*String(fields)*".jld2",Dict(String(fields)=> getfield(Sim,fields)))
    end
end

function SaveFields(Sim::SimData{R,I}, DataFileName::String, FieldDict::Dict{Any,Any}) where {R<:Real, I<:Integer}
    save(Sim.DataPath*DataFileName*".jld2", FieldDict)
end

function GetBigDataRecords(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    FilesInData =readdir(Sim.DataPath; join=true)
    Handles = Vector{String}([])
    split_len=0
    for to_get in FilesInData
        tmp = split(lstrip(to_get,'.'),".")
        split_len = length(tmp)
        filetype = tmp[split_len]
        if filetype=="jld2"
            tmp2  = split(tmp[1], "/")
            split_len = length(tmp2)
            push!(Handles, tmp2[split_len])
        end
    end
    deleteat!(Handles, findall(x->x=="All", Handles))
    return Handles
end

function ReadAllData!(Sim::SimData{R,I}; ReadBig=false) where {R<:Real, I<:Integer} 
    Dictionary= Dict()
    for fieldname in fieldnames(typeof(Sim))
        if String(fieldname)!="x" && String(fieldname)!="y" &&  String(fieldname)!="z" && String(fieldname)!="x_uw" && String(fieldname)!="y_uw" &&  String(fieldname)!="z_uw" && ~(fieldname in Sim.BigDataList)
            Dictionary[fieldname]= getfield(Sim,fieldname)
        end
    end
    olddata = load(Sim.DataPath*"$("All").jld2")

    BigDataFiles = GetBigDataRecords(Sim)
    if ~ReadBig 
        if length(BigDataFiles)>0
            printstyled("Reading of big data components has to be activated by ReadBig.\n";color=:yellow)
        end
    else
        for to_get in BigDataFiles
            setfield!(Sim, Symbol(to_get) ,load(Sim.DataPath*to_get*".jld2")[to_get])
        end
    end

    for (key, field) in Dictionary
        try
            if !(String(key) in BigDataFiles)
                setfield!(Sim, key, olddata[String(key)])
            end
        catch
           printstyled("Cannot find field: $(field) with key $(key).\n"; color=:yellow)
            continue
        end
    end
end

function ReadFields!(Sim::SimData{R,I},FileName::String, FieldDict::Dict{Symbol, Any}) where {R<:Real, I<:Integer}
    for (key, field) in FieldDict
        try
            tmp =load(Sim.DataPath*"$(FileName).jld2", String(key) )
            setfield!(Sim, key, tmp)
        catch
            printstyled("Cannot find field: $(field) with key $(key).\n"; color=:yellow)
            continue
        end
    end
end

function ClearData(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    close(Sim.xio)
    close(Sim.yio)
    close(Sim.zio)

    close(Sim.x_uw_io)
    close(Sim.y_uw_io)
    close(Sim.z_uw_io)

    for fieldname in fieldnames(typeof(Sim))
        if String(fieldname)!="x" || String(fieldname)!="y" ||  String(fieldname)!="z" || String(fieldname)!="x_uw" || String(fieldname)!="y_uw" ||  String(fieldname)!="z_uw"
            setfield!(Sim, fieldname, nothing)
        end
    end
end

function SetFieldsFromData(Sim::SimData{R,I}, Strings::Vector{String}) where {R<:Real, I<:Integer}
    Dictionary=Dict()
    for fieldname in fieldnames(typeof(Sim))
        Dictionary[String(fieldname)]= fieldname
    end

    Results = []
    for to_get in Strings
        if to_get in keys(Dictionary)
            setfield!(Sim, Dictionary[to_get], load(Sim.DataPath*"$("All").jld2" , to_get))
        else
            printstyled("Cannot load query: $(to_get).\n"; color=:yellow)
        end
    end
end

function GetFieldsFromData(Sim::SimData{R,I}, Strings::Vector{String}) where {R<:Real, I<:Integer}
    FilesInData =readdir(Sim.DataPath; join=true)

    Dictionary=Dict()
    for fieldname in fieldnames(typeof(Sim))
        Dictionary[String(fieldname)]= fieldname
    end

    BigDataFiles  = GetBigDataRecords(Sim)

    Results = []
    for to_get in Strings
        if to_get in BigDataFiles
            push!(Results, load(Sim.DataPath*"$(to_get).jld2", to_get))
        elseif to_get in keys(Dictionary)
            push!(Results, load(Sim.DataPath*"$("All").jld2",to_get))
        else
            printstyled("Cannot load query: $(to_get).\n"; color=:yellow)
        end
    end
    return Tuple(Results)
end

function WriteAsPDB(Sim::SimData{T,I}; Start=convert(I,1), Stop=convert(I,0), Stepsize=convert(I,10), Wrapped=false) where {T<:Real, I<:Integer}
    if Stop==0
        Stop=Int32(Sim.NSteps)
    end
    if Wrapped
        filename=Sim.DataPath*Sim.SimulationName*"_wrapped.pdb"
        WriteAsPDB_Help(Sim, Sim.x, Sim.y, Sim.z, filename, convert(eltype(Sim.NSteps),(Start)), convert(eltype(Sim.NSteps),Stop), convert(eltype(Sim.NSteps),Stepsize))
    else
        filename=Sim.DataPath*Sim.SimulationName*".pdb"
        WriteAsPDB_Help(Sim, Sim.x_uw, Sim.y_uw, Sim.z_uw, filename, convert(eltype(Sim.NSteps),(Start)), convert(eltype(Sim.NSteps),Stop), convert(eltype(Sim.NSteps),Stepsize))
    end
end

function WriteAsPDB_Help(Sim::SimData{T,I}, x::Matrix{T}, y::Matrix{T}, z::Matrix{T}, filename::String, start::I, stop::I, stepsize::I) where {T<:Real, I<:Integer}

    file = open(filename, "w+")
    CN="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz123456789"
    for step in start:stepsize:stop
        write(file, "REMARK Frame: $(step) \n")
        for chain in 1:Sim.NChains
            chain_letter = CN[mod(chain,(length(CN)-1))+1]
            for atom in Sim.ChainStart[chain]:Sim.ChainStop[chain]
                res_name = get(BioData.OneToThree, (Sim.IDToResName[Sim.IDs[atom]])[1], "XXX")
                write(file, "ATOM  $(@sprintf("%5i", atom))  CA  $(res_name) $(chain_letter)$(@sprintf("%4i", chain))    $(@sprintf("%8.3f",x[atom, step]))$(@sprintf("%8.3f",y[atom, step]))$(@sprintf("%8.3f",z[atom, step]))                          \n")
            end
            write(file, "TER   \n")#{num:>5}      {bioh.OneToThree[str(self.Sequence[chain][cnt])]} {self.CN[chain%self.CNN]}{cnt+1:>4}\n")
        end
        write(file, "ENDMDL\n")
    end
    write(file, "END")
end


@doc raw"""
    AlphaFold_startpos(ProteintoCif::Dict{String,String}, Proteins::Vector{String}, Sequences::Vector{String})

    Reads the position of the CA for all proteins out of the the *.cif files generated by AlphaFold.

**Arguments**:
- `ProteintoCif::Dict{String,String}`: Dictionary connection Protein name to file name.
- `Proteins::Vector{String}`: List of all names of all proteins in the simulation.
- `Sequences::Vector{String}`: List of amino acid sequences corresponding to the proteins.

**Returns**:
- A dictionary linking the protein name to AlphaFold structure.
- A dictionary linking the protein name to the number of amino acids.
"""
function AlphaFold_startpos(ProteintoCif::Dict{String,String}, Proteins::Vector{String}, Sequences::Vector{String})
    LengthDict = Dict{String, Int32}()
    for (Prot, Seq) in zip(Proteins,Sequences)
        LengthDict[Prot] =length(Seq)
    end
    PositionDict =Dict{String, Matrix{Float32}}()
    for protein in Set(Proteins)
        pos = zeros(Float64, LengthDict[protein], 3)
        index = 1
        open(ProteintoCif[protein], "r") do io
            for line in eachline(io)
                field = split(strip(line))
                if field[1] == "ATOM" && field[4] == "CA"
                    field = split(strip(line))
                    x = parse(Float64, field[11])
                    y = parse(Float64, field[12])
                    z = parse(Float64, field[13])
                    pos[index, :] = [x, y, z]
                    index += 1
                end
            end
        end
        shift= sum.(collect(extrema(pos,dims=1)))./2.0
        pos[:,1].-=shift[1]
        pos[:,2].-=shift[2]
        pos[:,3].-=shift[3]

        PositionDict[protein] = pos
    end
    return PositionDict, LengthDict
end

@doc raw"""
    getLargestBoundingbox(PositionDict::Dict{String, Matrix{<:Real}})

Returns the tuple  of the smallest and highest positions of any starting configuration in the dictionary.

**Arguments**:
- `PositionDict::Dict{String,String}`: A dictionary linking the protein name to AlphaFold structure.

**Returns**:
A tuple ([x_min, y_min, z_min], [x_max, y_max, z_max]).
"""
function getLargestBoundingbox(PositionDict::Dict{String, Matrix{R}}) where {R<:Real}
    min_val = [floatmax(),floatmax(),floatmax()] .* 0.9
    max_val = [floatmax(),floatmax(),floatmax()] .* -0.9

    for pos in values(PositionDict)
        min_val_  = [minimum(pos[:,i]) for i in 1:3]
        max_val_  = [maximum(pos[:,i]) for i in 1:3]

        for i in 1:3
            min_val[i] = min(min_val[i],min_val_[i])
            max_val[i] = max(max_val[i],max_val_[i])
        end
    end
    return min_val, max_val
end

function GetLastFrameForRestart(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    ### x,y,z are stored on disk and lazyly synchronised, ONLY IN READ MODE!!!
    Sim.x_uw =  Mmap.mmap(Sim.x_uw_io, Matrix{eltype(Sim.x)}, (Sim.NAtoms,Sim.NSteps))
    Sim.y_uw =  Mmap.mmap(Sim.y_uw_io, Matrix{eltype(Sim.x)}, (Sim.NAtoms,Sim.NSteps))
    Sim.z_uw =  Mmap.mmap(Sim.z_uw_io, Matrix{eltype(Sim.x)}, (Sim.NAtoms,Sim.NSteps))

    ### periodically unwrap one Axis
    pos = zeros(eltype(Sim.x), Sim.NAtoms, 3)
    pos[:,1] .= Sim.x_uw[:,Sim.NSteps]
    pos[:,2] .= Sim.y_uw[:,Sim.NSteps]
    pos[:,3] .= Sim.z_uw[:,Sim.NSteps]

    close(Sim.x_uw_io)
    close(Sim.y_uw_io)
    close(Sim.z_uw_io)
    return pos
end

function GetRestartData(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    return (Sim.Sequences, Sim.BoxSize, GetLastFrameForRestart(Sim) )
end

function getInfoForRestart(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    pos = zeros(R, Sim.NAtoms, 3)
    pos[:,1] .= Sim.x[:,1]
    pos[:,2] .= Sim.y[:,1]
    pos[:,3] .= Sim.z[:,1]

    #=Seqs = deepcopy(Sim.Sequences)
    for (i,Seq) in enumerate(Seqs)
        Seqs[i] = replace(Seq, "'"=>"", ","=>"") ### remove the modificatiosn
    end=#


    ### reverse engineer the modified to normal resnames for HOOMD
    ResToRes = Dict([(v[1],v[1]) for (_,v) in Sim.IDToResName])
    ResNameToID =Dict([(v,u) for (u,v) in Sim.IDToResName])
    val = values(Sim.IDToResName)
    modified_res =  [x for x in val if islowercase(x[1])]
    for res in modified_res
        id = ResNameToID[res]
        for c_id in 1:length(Sim.IDToMasses)
            if id ==c_id continue end
            if Sim.IDToMasses[id]-Sim.IDToMasses[c_id]==2 || Sim.IDToMasses[id]-Sim.IDToMasses[c_id]==16
                ResToRes[Sim.IDToResName[id][1]] = Sim.IDToResName[c_id][1]
                break
            end
        end
    end

    Seqs = Vector{String}()
    for Seq in Sim.Sequences
        push!(Seqs, string(join([ResToRes[x] for x in Seq])))
    end

    ### reverse engineer for LAMMPS
    # uppercase.(Seqs)

    return pos, uppercase.(Seqs), collect(permutedims(Sim.BoxSize, (2,1)))[:], Sim.TargetTemp
end

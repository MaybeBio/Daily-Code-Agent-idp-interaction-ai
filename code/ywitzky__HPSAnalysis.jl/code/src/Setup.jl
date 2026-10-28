include("./Setup/Polyply.jl")

module Setup

export createDenseStartingPosition, writeCollectedSlurmScript

using GSDFormat, HPSAnalysis, JSON, Printf, Mmap, LinearAlgebra, Random

include("../data/BioData.jl")
include("Setup/HOOMD_Setup.jl")
include("Setup/LAMMPS_Setup.jl")
include("Setup/Calvados_Setup.jl")
include("Setup/ENM_Setup.jl")


using .BioData


BondStrength = Dict("Calvados2"=>19.19,"TestModell"=>20.0 )
EpsilonAshbaughHatch = Dict("Calvados2"=> 0.20,"TestModell"=>0.2)
Cutoff = Dict("Calvados2"=> 40.0, "HPS-Alpha"=>35.0,"TestModell"=>35.0 )

Formats=Dict([("h5md", "position image create_group yes")])

#Set start coordinates for the AA, by adding offset for x,y-position and offset in z-position if Seqindex%ProteinsPerLayer==0
function createStartingPosition(Sequences, BoxSize, LatticeSpacing=8.0, BondSpacing=3.8)
    NSeq = length(Sequences)
    MaxSeqLength =0 
    for Seq in Sequences
        MaxSeqLength = max(MaxSeqLength, length(Seq))
    end

    pos = zeros(NSeq, MaxSeqLength, 3)
    ProteinsPerLayer=ceil(sqrt(NSeq))
    NLayer = ceil(NSeq/ProteinsPerLayer)
    zoffset =(BoxSize[5]+BoxSize[6])/2. -((NLayer-1)/2.)*LatticeSpacing
    for (Seqindex, Seq) in enumerate(Sequences)
        xoffset = (BoxSize[1]+BoxSize[2])/2.
        xoffset += ( (Seqindex%ProteinsPerLayer) -(ProteinsPerLayer-1)/2.)*LatticeSpacing
        yoffset = (BoxSize[3]+BoxSize[4])/2.
        yoffset -= MaxSeqLength/2. *BondSpacing
        for (Aaid,_) in enumerate(Seq)
            pos[Seqindex, Aaid, 1] = xoffset
            pos[Seqindex, Aaid, 2] = yoffset
            pos[Seqindex, Aaid, 3] = zoffset
            yoffset += BondSpacing
        end
        if Seqindex%ProteinsPerLayer==0
            zoffset += LatticeSpacing
        end
    end
    return pos
end

function createDenseStartingPosition(Sequences, BoxSize, LatticeSpacing=8.0, BondSpacing=3.8)
    NSeq = length(Sequences)
    MaxSeqLength =0 
    for Seq in Sequences
        MaxSeqLength = max(MaxSeqLength, length(Seq))
    end

    pos = zeros(NSeq, MaxSeqLength, 3)
    xoffset= BoxSize[1]+0.5
    yoffset=BoxSize[3]+0.5
    zoffset =BoxSize[5]+0.5
    forward = 1.
    xforward = 1.
    ytmp= 0.
    for (Seqindex, Seq) in enumerate(Sequences)
        for (Aaid,_) in enumerate(Seq)
            ytmp =  yoffset+BondSpacing*forward
            if ytmp <(BoxSize[4] -0.5) && ytmp >(BoxSize[3] +0.5)
                yoffset = ytmp
            else
                xtmp =  xoffset+LatticeSpacing*xforward
                if xtmp<(BoxSize[2] -0.5) && xtmp >(BoxSize[1] +0.5)
                    xoffset = xtmp ### otherwise overlappting chains, non ideal harmonic bond will equilibrate fast
                else
                    xforward *=-1.
                    zoffset += LatticeSpacing
                    #xoffset = BoxSize[1]+0.5
                end
                forward *= -1.
            end
            pos[Seqindex, Aaid, 1] = xoffset
            pos[Seqindex, Aaid, 2] = yoffset
            pos[Seqindex, Aaid, 3] = zoffset
        end
    end
    return pos
end

function correctPositionInBounds(pos, boxSize, Sequences)
    box = boxSize/2.0
    box_neg = -1.0.*box

    for (SeqID, Seq) in enumerate(Sequences)
        for (atom,_) in enumerate(Seq)
            for i in 1:3
                if pos[SeqID, atom, i]>box[i]
                    pos[SeqID, atom, i] -= boxSize[i]
                end
                if pos[SeqID, atom, i]< box_neg[i]
                    pos[SeqID, atom, i] += boxSize[i]
                end
            end
        end
    end
    return pos
end

function rescalePositions(pos, boxSize, Sequences, image, rescale_factor=5.5/3.8)
    uw_pos = deepcopy(pos) 
    for i in 1:3
        uw_pos[:,:,i] .+= boxSize[i] .*image[:,:,i]
    end
    uw_pos ./= rescale_factor
    return correctPositionInBounds(uw_pos, boxSize, Sequences)
end

function getImageCopyNumber(pos, boxSize, Sequences)
    diff = zeros(eltype(pos), 3)
    image = zeros(Int32, size(pos))
    tmp_pos = deepcopy(pos)

    for (SeqID, Seq) in enumerate(Sequences)
        for (atom,_) in enumerate(Seq)
            if atom ==1 continue end
            for i in 1:3
                diff[i] = tmp_pos[SeqID, atom, i]-tmp_pos[SeqID, atom-1, i]

                if abs(diff[i]-boxSize[i])<abs(diff[i]) # || abs(diff[i]+boxSize[i])<abs(diff[i])
                    image[SeqID, atom,i] =  image[SeqID, atom-1,i]-1 #-1*convert(Int32,(round((diff[i]/(boxSize[i]/2.)/2.0))))  + image[SeqID, atom-1,i]
                    #tmp_pos[SeqID, atom, i] +=
                elseif abs(diff[i]+boxSize[i])<abs(diff[i])
                    image[SeqID, atom,i] =  image[SeqID, atom-1,i]+1
                    #image[SeqID, atom,i] = -1*convert(Int32,(round((diff[i]/(boxSize[i]/2.)/2.0))))  + image[SeqID, atom-1,i]
                else
                    image[SeqID, atom,i] = image[SeqID, atom-1,i]
                end
            end
        end
    end

    return image
end

@doc raw"""
    writeStartConfiguration(fileName, StartFileName, Info, Sequences, BoxSize,NSteps=100_000_000; SimulationType="Calvados2", Temperature=300,MixingRule="1-1001-1", Pos =zeros(Float32, 0),InitStyle="Slab", SaltConcentration=0.15, pH=6.0, ChargeTemperSteps=[], ChargeTemperSwapSteps=100_000, HOOMD=false, OneToChargeDef=BioData.OneToHPSCharge, OneToLambdaDef=BioData.OneToCalvados2Lambda, OneToSigmaDef=BioData.OneToHPSCalvadosSigma,WriteOutFreq=100_000, Device="GPU", yk_cut=40.0, ah_cut=20.0)

Writes the start configuration for a molecular dynamics simulation.
    
**Arguments**
- `fileName::String`: Name of the output file.
- `StartFileName::String`: Name of the initial configuration file.
- `Info::String`: Additional information about the simulation.
- `Sequences::Vector{String}`: List of amino acid sequences corresponding to the proteins.
- `BoxSize::Vector{Float}`: A vector of minmal/maximal box dimensions in each axis. ([x_min, x_max, y_min, y_max, z_min, z_max]).
- `NSteps::Int`: Number of simulation steps (default: 100,000,000).
- `SimulationType::String`: IDP model used for simulations.(default: "Calvados2").
- `Temperature::Float`: Temperature in Kelvin (default: 300).
- `MixingRule::String`: Mixing rule for optional dihedral potential for Calvados models.
- `Pos::Vector{Float}`: Initial positions. (default: empty array).
- `InitStyle::String`: Initialization style, e.g., "Slab" or "Pos".
- `SaltConcentration::Float`: Salt concentration in M (default: 0.15).
- `pH::Float`: pH level (default: 6).
- `ChargeTemperSteps`: List of charge tempering steps.
- `ChargeTemperSwapSteps::Int`: Swap steps for charge tempering.
- `HOOMD::Boolean`: Boolean to enable HOOMD compatibility.
- `OneToChargeDef::Dict`: Dictionary defining the amount of charge for each one letter atom type.
- `OneToSigmaDef::Dict`: Dictionary defining the sigma parameter in LJ potentials for each one letter atom type.
- `OneToLambdaDef::Dict`: Dictionary defining the lambda paramter in Ashbaugh-Hatch potentials for each one letter atom type.
- `WriteOutFreq::Int`: Frequency of writing output (default: 100,000).
- `Device::String`: Computational device, "CPU"/"GPU" (default: "GPU").
- `yk_cut::Float`: Cutoff distances for yukawa potential.
- `ah_cut::Float`: Cutoff distances for ashbaugh hatch potential.
- `domain::List(Int)`: Domains in which the ENM is active.
- `ENM::Tuple`: Data that are nessesary for ENM (bond name, id, group).

**Creates**:
* Writes data files with the start configuration.
"""
function writeStartConfiguration(BasePath, fileName, StartFileName, Info, Sequences, BoxSize,NSteps=100_000_000; SimulationType="Calvados2", Temperature=300,MixingRule="1-1001-1", Pos =zeros(Float32, 0),InitStyle="Slab", SaltConcentration=0.15, pH=6, ChargeTemperSteps=[], ChargeTemperSwapSteps=100_000, HOOMD=false, OneToChargeDef=BioData.OneToHPSCharge, OneToLambdaDef=BioData.OneToCalvados2Lambda, OneToSigmaDef=BioData.OneToHPSCalvadosSigma,WriteOutFreq=100_000, Device="GPU", yk_cut=40.0, ah_cut=20.0, domain=Array([]), ENM=nothing, SlabAxis=2)

    ChargeTemperSim=length(ChargeTemperSteps)>0

    NAtoms,NBonds,NAngles,NDihedrals,AlphaAddition,SimulationType,AtomTypes, LongAtomTypes, AaToId, IdToAa,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, NAtomTypes, dihedral_short_map, dihedral_long_map, dihedral_eps, dihedral_list = startConfigurationSetup(Sequences,SimulationType,pH,OneToChargeDef,OneToLambdaDef,OneToSigmaDef,MixingRule)


    #Set start coordinates for the AA, with different variations
    if InitStyle=="Slab"
        #By adding offset for x,y-position and offset in z-position if Seqindex%ProteinsPerLayer==0
        pos = createStartingPosition(Sequences, BoxSize)
    elseif (InitStyle=="Pos")
        #By getting Initial Positions
        cnt = 1
        pos_init = Pos
        pos = zeros(eltype(Pos),length(Sequences), maximum(length.(Sequences)), 3)
        for (Seqindex, Seq) in enumerate(Sequences)
            for (Aaid,_) in enumerate(Seq)
                pos[Seqindex, Aaid, 1] = pos_init[cnt,1]
                pos[Seqindex, Aaid, 2] = pos_init[cnt,2]
                pos[Seqindex, Aaid, 3] = pos_init[cnt,3]
                cnt+=1
            end
        end
    else
        #By getting Initial Positions
        pos = createDenseStartingPosition(Sequences, BoxSize)
    end

    #AltBox=[BoxLengthShort,BoxLengthLong,BoxLengthShort]
    AltBox = [BoxSize[2]-BoxSize[1], BoxSize[4]-BoxSize[3], BoxSize[6]-BoxSize[5]]
    #Periodic Boundery, if outside -> define coordinates inside and set image
    pos = Setup.correctPositionInBounds(pos, AltBox, Sequences) ### poly coordinates are not in box sometimes
    pos = Setup.correctPositionInBounds(pos, AltBox, Sequences) ### 2 times is the charm....
    pos = Setup.correctPositionInBounds(pos, AltBox, Sequences) ### 3 times is the charm....

    image = Setup.getImageCopyNumber(pos, AltBox, Sequences)

    #domain = []
    #Write all Inputs, Parameters (Yukawa Interaction with Debye-Hückle), Dictionaries and the Start-File
    if HOOMD
        writeHOOMD(BasePath, Sequences,pos,image,OneToCharge,AaToId,OneToMass,OneToSigma,OneToLambda,AlphaAddition,dihedral_long_map,dihedral_eps,SimulationType,Temperature,SaltConcentration,BoxSize,StartFileName,NSteps,WriteOutFreq,Device,yk_cut,ah_cut,pH,domain,NAtoms,NBonds,NAngles,NDihedrals,dihedral_short_map,dihedral_list, ENM, SlabAxis)
    else
        writeHPSLammps(fileName,Sequences,AtomTypes,LongAtomTypes,LongAtomTypesToRes,InitStyle,ChargeTemperSteps,ChargeTemperSwapSteps,pos,image,OneToCharge,AaToId,OneToMass,OneToSigma,OneToLambda,AlphaAddition,dihedral_long_map,dihedral_eps,SimulationType,Temperature,SaltConcentration,BoxSize,StartFileName,NSteps,WriteOutFreq,pH,NAtoms,NBonds,NAngles,NDihedrals,Info,NAtomTypes)
    end
end

function getBonds(Sequences::Vector{String}; M=2)
    N = sum(length.(Sequences).-(M-1))
    bonds = zeros(Int32, (N, M))
    bondid = 0
    atomid = 0
    for (SeqId, Sequence) in enumerate(Sequences)
        for (ResId,Res) in enumerate(Sequence)
            if ResId>(length(Sequence)-M+1)
                atomid += 1
                continue
            else
                bondid += 1
                for (id, i) in enumerate(0:M-1)
                    bonds[bondid,id] = atomid + i 
                end
                atomid += 1
            end
        end
    end 
    return bonds
end

@doc raw"""
    writeGSDStartFile(FileName::String, NAtoms::I, NBonds::I, NAngles::I, NDihedrals::I,Box::Vector{R}, Positions::Array{R}, AaToId::Dict{Char, <:Integer},Sequences,  InputImage::Array{I2}, InputMasses::Array{<:Real}, InputCharges::Array{R}, DihedralMap::Dict, DihedralList::Matrix{<:Integer}, AaToSigma::Dict{Char, <:Real}, UseAngles::Bool, SimulationType::String, Domains, ENM) where {I<:Integer, R<:Real, I2<:Integer}

A GSD data file is written, that include the parameters for the Simulation witch are given as Arguments.
    
**Arguments**
- `FileName::String`: Name of the output GSD file.
- `NAtoms::Int`: The total number of amino acids (atoms) in the system.
- `NBonds::Int`: The number of bonds between amino acids (`NAtoms - 1`).
- `NAngles::Int`: The number of angles formed between amino acids (`NAtoms - 2`).
- `NDihedrals::Int`: The number of dihedral angles between amino acids (`NAtoms - 3`).
- `Box::Vector{Float}`: 3-element vector specifying the dimensions of the simulation box
- `Positions::Array{Float}`: Array of the atomic coordinates of the Aminoacids.
- `AaToId::Dict{Char, Int}`: Dictionary mapping each amino acid type to a unique ID number.
- `Sequences::Array{String}`: A list of protein sequences, where each sequence is represented as a string of amino acids.
- `InputImage::Array{Int}`: An array used to determine periodic boundary conditions and correct atomic positions.
- `InputMasses::Array{Float}`: A 1D array specifying the mass of each amino acid.
- `InputCharges::Array{Float}`: A 1D array specifying the electric charge of each amino acid.
- `DihedralMap::Dict`: A dictionary mapping unique dihedral angle definitions (four atom indices) to dihedral types.
- `DihedralList::Matrix{Int}`: A matrix where each row defines a specific dihedral angle using atom indices.
- `AaToSigma::Dict{Char, <:Real}`: A dictionary mapping amino acid types to their Lennard-Jones σ-parameter (used in force field calculations).
- `UseAngles::Bool`: If `true`, angle and dihedral interactions are included in the GSD file.
- `SimulationType::String`: Type of simulation.
- `ENM::Tuple`: Data that are nessesary for ENM (bond name, id, group).
- `DomainRanges`::Vector{UnitRange}: Range of ENMs

**Creates**:
* Writes the GSD data files.
"""
function writeGSDStartFile(FileName::String, NAtoms::I, NBonds::I, NAngles::I, NDihedrals::I,Box::Vector{R}, Positions::Array{R}, AaToId::Dict{Char, <:Integer},Sequences,  InputImage::Array{I2}, InputMasses::Array{<:Real}, InputCharges::Array{<:Real}, DihedralMap::Dict, DihedralList::Matrix{<:Integer}, AaToSigma::Dict{Char, <:Real}, UseAngles::Bool, SimulationType::String,  ENM, groups) where {I<:Integer, R<:Real, I2<:Integer}
 
    snapshot = GSDFormat.Frame()    
    snapshot.configuration.step = 1 
    snapshot.configuration.dimensions = 3 
    snapshot.configuration.box = [Box[1],Box[2], Box[3], 0, 0, 0]./10.0 #4:6 are tilt
    snapshot.particles.N = NAtoms
    # pos = zeros(length(Sequences), maximum(length.(Sequences)), 3)
    # reshape -> [N_total, 3]
    snapshot.particles.position = reshape(permutedims(Positions, (2,1,3)), (size(Positions, 1)*size(Positions, 2), 3))./10.0 ### permute to get alignment in memory, reshape to match gsd formart
    IdToAa = Dict(value => key for (key, value) in AaToId)
    snapshot.particles.types =  [string(IdToAa[Id]) for Id in  sort(collect(values(AaToId))) ] #string.(collect(keys(AaToId)))
    snapshot.particles.typeid = [Int32(AaToId[AA])-1 for AA in join(Sequences)] ### convert to python numbering
    snapshot.particles.image = reshape(permutedims(InputImage, (2,1,3)), (size(InputImage, 1)*size(InputImage, 2), 3)) ### permute to get alignment in memory, reshape to match gsd formart
    snapshot.particles.mass = InputMasses
    snapshot.particles.charge = InputCharges
    snapshot.particles.diameter =  [Float32(AaToSigma[AA])/10.0  for AA in join(Sequences)]
    snapshot.particles.body = groups

    ### Bond_data.group = (self.N, getM(data)
    B_N, B_types, B_typeid, B_group_matrix = NBonds, ["O-O"], zeros(Int32, NBonds), getBonds(Sequences, M=2)
    # Create Bonds
    if SimulationType == "Calvados3"
        ENMB_N, ENMB_types, ENMB_typeid, ENMB_group_vector, harmonic = ENM
        B_N = ENMB_N
        B_types = ENMB_types
        B_typeid = ENMB_typeid
        ENMB_group_matrix = permutedims(hcat(collect.(ENMB_group_vector)...))
        B_group_matrix = ENMB_group_matrix
    end

    snapshot.bonds.N = B_N
    snapshot.bonds.types = B_types
    snapshot.bonds.typeid = B_typeid
    snapshot.bonds.group = B_group_matrix

    if UseAngles
        # Create Angles
        snapshot.angles.N = NAngles
        snapshot.angles.types = ["O-O-O"]
        snapshot.angles.typeid = zeros(Int32, NAngles)
        snapshot.angles.group = getBonds(Sequences, M=3)

        # Create Dihedrals
        snapshot.dihedrals.N =  NDihedrals 
        snapshot.dihedrals.types = ["$(ids[1])-$(ids[2])-$(ids[3])-$(ids[4])" for ids in collect(keys(DihedralMap))]#string.(collect(values(DihedralMap)))
        snapshot.dihedrals.typeid = [DihedralMap[DihedralList[key,:]]-1 for key in axes(DihedralList,1)] ### convert to python numbering
        snapshot.dihedrals.group = getBonds(Sequences, M=4)
    end

    file = GSDFormat.open(FileName, 'w')
    GSDFormat.append(file, snapshot)
    GSDFormat.close(file)
end

"""
    rotate_randomly(pos:, ROTATION_MATRICES_90, rng)

The function picks a random 90°‑multiple rotation from the 24 predefined rotationmatrices,
returns the according rotated positions of the (N,3) input position 

**Arguments**
- `pos::Matrix{R}`: Positions to be rotated.
- `ROTATION_MATRICES_90::Vector{Matrix{R}}`: Predefined vector of 3x3 rotation matrices.
- `rng::AbstractRNG=Random.GLOBAL_RNG`: Random number generator.

**Returns**
- Rotated positions.
"""
@inline function rotate_randomly(pos::Matrix{R}, ROTATION_MATRICES_90::Vector{Matrix{R}},
                         rng::AbstractRNG=Random.GLOBAL_RNG) where {R<:Real}
    # Pick a random rotation matrix (stored as Int8)
    RotMat = rand(rng, ROTATION_MATRICES_90)   # 3×3

    # Pick a random rotation matrix (stored as Int8)
    return copy(pos) * RotMat
end

@doc raw"""
    CreateStartConfiguration(SimulationName::String, Path::String, BoxSize::Vector{R}, Proteins::Vector{String}, Sequences::Vector{String} ; Axis=`y`, Regenerate=true,SimulationType="Calvados2",ProteinToDomain=Dict(),ProteinToCif=Dict())

Creates the file structure and initialises particle positions for the given parameters.

**Arguments**:
- `SimulationName::String`: The name of the simulation.
- `Path::String`: The base directory where simulation data will be stored.
- `BoxSize::Vector{R}`: A vector defining the box dimensions (x, y, z).
- `Proteins::Vector{String}`: List of protein names used in the simulation.
- `Sequences::Vector{String}`: List of amino acid sequences corresponding to the proteins.
**Optional Arguments**:
- `Axis::String`: The axis along which the system is unfolded.
- `Regenerate::Bool`: If true, regenerates initial positions using the Polyply package.
- `SimulationType::String`: Type of Simulation (default is Calvados2).
- `ProteinToDomain::Dict`: Dictionary of domains for the proteins.
- `ProteinToCif::Dict`: Dictionary of the AlphaFold data files for each protein.

**Returns**:
* A tuple (pos, Data) containing the initial positions and the simulation data structure.
"""
function CreateStartConfiguration(SimulationName::String, Path::String, BoxSize::Vector{R}, Proteins::Vector{String}, Sequences::Vector{String} ; Axis="y", Regenerate=true,SimulationType="Calvados2",ProteinToDomain=Dict(),ProteinToCif=Dict()) where {R<:Real}
    #Definition of Paths for the parameters
    Data = HPSAnalysis.SimData()
    Data.BasePath= Path
    Data.PlotPath=Data.BasePath*"/Plots/"
    Data.DataPath=Data.BasePath*"/Data/"
    Data.xFilePath = Data.DataPath*"x.bin"
    Data.yFilePath = Data.DataPath*"y.bin"
    Data.zFilePath = Data.DataPath*"z.bin"
    Data.x_uw_FilePath = Data.DataPath*"x_uw.bin"
    Data.y_uw_FilePath = Data.DataPath*"y_uw.bin"
    Data.z_uw_FilePath = Data.DataPath*"z_uw.bin"
    Data.Reduce=1 ### only a factor for prereduced data.
    mkpath(Data.PlotPath)
    mkpath(Data.DataPath)

    #Definition of the Box
    Data.BoxLength = [Float32(BoxSize[1]),Float32(BoxSize[2]),Float32(BoxSize[3])]
    Data.BoxSize = zeros(eltype(Data.x), 3,2 )#Matrix of FloatType
    Data.BoxSize[1,1] = -Data.BoxLength[1]/2.
    Data.BoxSize[1,2] =  Data.BoxLength[1]/2.
    Data.BoxSize[2,1] = -Data.BoxLength[2]/2.
    Data.BoxSize[2,2] =  Data.BoxLength[2]/2.
    Data.BoxSize[3,1] = -Data.BoxLength[3]/2.
    Data.BoxSize[3,2] =  Data.BoxLength[3]/2.

    #Definition: Number of Steps,Chains
    Data.NSteps = 1 ### only for the creation creation, will be changed later.
    Data.NChains = length(Sequences)
    Data.StepFrequency=1
    Data.SimulationName=SimulationName

    Data.Sequences = Sequences

    #Definition of the length of all chains
    Data.NAtoms = sum(length.(Data.Sequences))

    #All Residues get an ID, same Residue -> same ID, and definition for other way around
    Data.IDs = zeros(eltype(Data.NAtoms), Data.NAtoms)
    cnt = Int32(1)
    ind = 1
    ResNameToID = Dict{Char, Int32}()
    for chain in Data.Sequences
        for AA in chain
            if !haskey(ResNameToID,AA)
                ResNameToID[AA] = cnt
                cnt += 1
            end
            Data.IDs[ind] = ResNameToID[AA]
            ind += 1
        end
    end
    Data.IDToResName = Dict( (v => string(k)) for (k, v) in ResNameToID)

    #Matrix of the start and end values of each chain
    Data.ChainStart = zeros(eltype(Data.NSteps), Data.NChains)
    Data.ChainStop  = zeros(eltype(Data.NSteps), Data.NChains)
    Data.ChainStart[1] = 1
    Data.ChainStop[1] = length(Data.Sequences[1])
    for (i,Seq) in enumerate(Data.Sequences)
        if i ==1 
            continue
        end
        Data.ChainStart[i] = Data.ChainStart[i-1]+length(Seq)
        Data.ChainStop[i]  = Data.ChainStop[i-1] +length(Seq)
    end

    ### allocate disk space for X
    #creat path for the coordinates
    Data.xio= open(Data.xFilePath,"w+")
    Data.yio= open(Data.yFilePath,"w+")
    Data.zio= open(Data.zFilePath,"w+")
    #mmap creat 2D-Matrix for the coordinates (Lenght, Steps(=1))
    #Matrices of start coordinates (for each x,y,z-value)
    Data.x =  Mmap.mmap(Data.xio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.y =  Mmap.mmap(Data.yio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.z =  Mmap.mmap(Data.zio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))


    ######## Start Generation of Initial positions with Polyply
    InitFiles= "$(Data.BasePath)/InitFiles/"
    mkpath(InitFiles)

    ### all possible rotation matrixes by 90 degrees
    ROTATION_MATRICES_90 = Vector{Matrix{R}}([[-1 0 0; 0 -1 0; 0 0 1], [-1 0 0; 0 1 0; 0 0 -1], [1 0 0; 0 -1 0; 0 0 -1], [1 0 0; 0 1 0; 0 0 1], [-1 0 0; 0 0 -1; 0 -1 0], [-1 0 0; 0 0 1; 0 1 0], [1 0 0; 0 0 -1; 0 1 0], [1 0 0; 0 0 1; 0 -1 0], [0 -1 0; -1 0 0; 0 0 -1], [0 -1 0; 1 0 0; 0 0 1], [0 1 0; -1 0 0; 0 0 1], [0 1 0; 1 0 0; 0 0 -1], [0 -1 0; 0 0 -1; 1 0 0], [0 -1 0; 0 0 1; -1 0 0], [0 1 0; 0 0 -1; -1 0 0], [0 1 0; 0 0 1; 1 0 0], [0 0 -1; -1 0 0; 0 1 0], [0 0 -1; 1 0 0; 0 -1 0], [0 0 1; -1 0 0; 0 -1 0], [0 0 1; 1 0 0; 0 1 0], [0 0 -1; 0 -1 0; -1 0 0], [0 0 -1; 0 1 0; 1 0 0], [0 0 1; 0 -1 0; 1 0 0], [0 0 1; 0 1 0; -1 0 0]])
    rng = MersenneTwister()

    if SimulationType=="Calvados3"
        if Regenerate
            if Axis != "y"
                @error("I can only setup with slab axis being y")
            end
            ###Creat a pdb data from the AlphaFold cif data
            HPSAnalysis.RewriteCifToPDB(Data.BasePath,ProteinToCif, Proteins )

            ####  get individual bounding boxes, compute the largest one
            PositionDict, LengthDict = HPSAnalysis.AlphaFold_startpos(ProteinToCif, Proteins, Sequences) 
            min_vals, max_vals = HPSAnalysis.getLargestBoundingbox(PositionDict)

            max_length = maximum(max_vals .- min_vals .+ 20.0) ### add 20 Angstrom padding for excluded volume
            max_N = div.(Data.BoxLength, max_length)

            count = 0
            NLayers = div(Data.NChains, max_N[3]*max_N[1])+1
            yoffset = ( (NLayers-1)*max_length)/2.0 
            old = 1
            for iy in 0:max_N[2]-1 ### Y is the typicall slab axis
                for iz in 0:max_N[3]-1
                    for ix in 0:max_N[1]-1
                        count += 1
                        if count > Data.NChains
                            break
                        end
                        pos = PositionDict[Proteins[count]]
                        proteinlength = LengthDict[Proteins[count]]

                        offset = [(ix+0.5)*max_length, (iy+0.5)*max_length-yoffset, (iz+0.5)*max_length] .+10.0

                        tmp_pos = rotate_randomly(pos, ROTATION_MATRICES_90, rng)

                        Data.x[old:old+proteinlength-1, 1] = tmp_pos[:, 1] .+ offset[1]
                        Data.y[old:old+proteinlength-1, 1] = tmp_pos[:, 2] .+ offset[2] 
                        Data.z[old:old+proteinlength-1, 1] = tmp_pos[:, 3] .+ offset[3]
                        old += proteinlength
                    end
                end
            end
        end
    elseif SimulationType=="Calvados2"
        if Regenerate
            ### generate Martini ITP Files
            mkpath("$(InitFiles)ITPS_Files/")
            HPSAnalysis.Polyply.GenerateITPFilesOfSequence(Proteins, Data.Sequences, "$(InitFiles)ITPS_Files/")

            ### Generate Topology files
            TopologyFile = "$(InitFiles)$(Data.SimulationName).top"
            HPSAnalysis.Polyply.GenerateSlabTopologyFile(TopologyFile,"$(InitFiles)ITPS_Files/", Proteins, Data.SimulationName)

            ### generate coordinates
            HPSAnalysis.Polyply.GenerateCoordinates(InitFiles, Data.SimulationName, BoxSize/10.0, TopologyFile)

            ### convert to PDB
            #Polyply.ConvertGroToPDB(InitFiles, Data.SimulationName)

            ### read positons from gro
            HPSAnalysis.Polyply.readSimpleGRO("$(InitFiles)$SimulationName.gro", Data.x,Data.y,Data.z)
        end
    end


    ### sync RAM to disk before closing,
    Mmap.sync!(Data.x)
    Mmap.sync!(Data.y)
    Mmap.sync!(Data.z)
    close(Data.xio)
    close(Data.yio)
    close(Data.zio)

    ########## unfold along y-Axis so that the slab can be extended

    ### x,y,z are stored on disk and lazyly synchronised, ONLY IN READ MODE!!!
    Data.xio= open(Data.xFilePath,"r+")
    Data.yio= open(Data.yFilePath,"r+")
    Data.zio= open(Data.zFilePath,"r+")

    Data.x =  Mmap.mmap(Data.xio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.y =  Mmap.mmap(Data.yio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.z =  Mmap.mmap(Data.zio, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))

    HPSAnalysis.unfoldPositions(Data)

    Data.x_uw_io= open(Data.x_uw_FilePath,"r+")
    Data.y_uw_io= open(Data.y_uw_FilePath,"r+")
    Data.z_uw_io= open(Data.z_uw_FilePath,"r+")

    Data.x_uw =  Mmap.mmap(Data.x_uw_io, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.y_uw =  Mmap.mmap(Data.y_uw_io, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))
    Data.z_uw =  Mmap.mmap(Data.z_uw_io, Matrix{eltype(Data.x)}, (Data.NAtoms,Data.NSteps))

    ### periodically unwrap one Axis
    pos = zeros(eltype(Data.x), Data.NAtoms, 3)
    if Axis=="x"
        pos[:,1] .= Data.x_uw[:,1]
        pos[:,2] .= Data.y[:,1]
        pos[:,3] .= Data.z[:,1]
    elseif Axis=="y"
        pos[:,2] .= Data.y_uw[:,1]
        pos[:,1] .= Data.x[:,1]
        pos[:,3] .= Data.z[:,1]
    elseif Axis=="z"
        pos[:,3] .= Data.z_uw[:,1]
        pos[:,1] .= Data.x[:,1]
        pos[:,2] .= Data.y[:,1]
    end

    ### shift so that box center is at 0,0,0
    pos[:,1] .-= BoxSize[1]/2
    if SimulationType!="Calvados3"
        pos[:,2] .-= BoxSize[2]/2 # y as slab axis is centered before
    end
    pos[:,3] .-= BoxSize[3]/2

    Data.x[:,1] .=  pos[:,1]
    Data.y[:,1] .=  pos[:,2]
    Data.z[:,1] .=  pos[:,3]

    Data.x_uw[:,1] .-= BoxSize[1]/2
    Data.y_uw[:,1] .-= BoxSize[2]/2
    Data.z_uw[:,1] .-= BoxSize[3]/2

    tmp_x = Data.x_uw
    tmp_z = Data.z_uw
    Data.x_uw= deepcopy(Data.x)
    Data.z_uw = deepcopy(Data.z)
    HPSAnalysis.WriteAsPDB(Data; Wrapped=false)
    Data.x_uw = tmp_x
    Data.z_uw = tmp_z

    ### since they are not needed anymore
    close(Data.xio)
    close(Data.yio)
    close(Data.zio)
    close(Data.x_uw_io)
    close(Data.y_uw_io)
    close(Data.z_uw_io)

    return (pos, Data) 
end

end #module Setup


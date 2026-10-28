@doc raw"""
    WriteHOOMDSequences(filename::String, Sequences)

Write a data file with the sequence of each protein.
    
**Arguments**
- `filename::String`: Path where the data will be saved.
- `Sequences::Array{String}`: List of protein sequences.

**Creat**:
* Write a file with all Sequences.
"""
function WriteHOOMDSequences(filename::String, Sequences)
    io = open(filename, "w");

    for (SeqID, Seq) in enumerate(Sequences)
        write(io, "$Seq\n")
    end
    close(io);
end

@doc raw"""
    WriteDictionaries(filename::String, ToCharge, ToID, ToMass, ToDiameter, ToLambda)

Write a data file with all Dictionaries that are used in the Simulation.
    
**Arguments**
- `filename::String`: Path where the data will be saved.
- `ToCharge::Dict()`: Dictionary defining the charge for the Aminoacids.
- `ToID::Dict()`: Dictionary mapping each amino acid type to a unique ID number.
- `ToMass::Dict()`: Dictionary defining the mass for the Aminoacids.
- `ToDiameter::Dict()`: Dictionary defining the diameter for the Aminoacids.
- `ToLambda::Dict()`: Dictionary defining the Lambda for the Aminoacids.

**Creat**:
* Write a file with all Dictionaries.
"""
function WriteDictionaries(filename::String, ToCharge, ToID, ToMass, ToDiameter, ToLambda)
    io = open(filename, "w")
    write(io, "// ID,resname, Charge, Mass, λ   \n")

    for key in sort(collect(keys(ToID)), by=x->ToID[x])
        write(io, " $(ToID[key]), $key, $(ToCharge[key]), $(ToMass[key]), $(ToDiameter[key]), $(ToLambda[key])\n")
    end
    close(io);
end

@doc raw"""
    WriteHOOMDParticlesInput(filename::String, pos::Array{R}, ToCharge, ToID, Sequences, ToMass, ToDiameter, image) where {R<:Real}

Write a data file that contains position, charge, mass and diameter for the the amino acids.
    
**Arguments**
- `filename::String`: Path where the data will be saved.
- `pos::Array`: List of positions for each amino acid.
- `Sequences::Array{String}`: List of protein sequences.
- `ToMass::Dict()`: Dictionary defining the mass for the amino acids.
- `ToDiameter::Dict()`: Dictionary defining the diameter for the amino acids.
- `image::Array{<Integer}`: Number of the periodic image of each amino acid.

**Creat**:
* Write a file with the escential datas of each Aminoacid.
"""
function WriteHOOMDParticlesInput(filename::String, pos::Array{R}, ToCharge, ToID, Sequences, ToMass, ToDiameter, image) where {R<:Real}

    io = open(filename, "w");
    write(io, "### N, id, x ,y ,z , charge, m, diameter \n");
    cnt=1
    for (SeqID, Seq) in enumerate(Sequences)
        for (atom,res) in enumerate(Seq)
            write(io, "$(cnt), $(ToID[res] -1) , $(@sprintf("%.3f",pos[SeqID,atom, 1])), $(@sprintf("%.3f",pos[SeqID,atom, 2])), $(@sprintf("%.3f",pos[SeqID,atom, 3])), $(ToCharge[res]), $(ToMass[res]), $(ToDiameter[res]) , $(image[SeqID, atom,1]) , $(image[SeqID, atom,2]) , $(image[SeqID, atom,3]) \n");
            cnt += 1
        end
    end

    close(io);
end

function WriteDihedrals(filename, dihedral_map, dihedral_eps)
    io = open(filename, "w");
    write(io , "# $(dihedral_eps) \n")
    for key in keys(dihedral_map)
        write(io, "$(key[1]) $(key[2]) $(key[3]) $(key[4]) $(dihedral_map[key])\n")
    end

    close(io)
end

function WriteGroups(filename, groups)
    io = open(filename, "w");
    for val in groups
        write(io, "$(val)\n")
    end
    close(io)
end

@doc raw"""
    WriteParams(filename, SimName, Temp, NSteps, NOut, Timestep, Box, Seed; Minimise=true, TrajectoryName="traj.gsd", UseAngles=true, UseCharge=true, Alt_GSD_Start="-", Create_Start_Config=false, ϵ_r=1.73136, κ=1.0, Device="GPU", yk_cut=4.0, ah_cut=2.0, ionic=0.1, pH=7.0, SimType="Calvados2",domain=Array([[0,0]]))

Write a data file that contains all Parameters of the Simulation.
    
**Arguments**
- `filename::String`: Path where the data will be saved.
- `SimName::String`: Name of the Simulation.
- `Temp::Float`: Temperature of the Simulation.
- `NSteps::Float`: The Number of steps for the Simulation.
- `NOut::Integer`: The Number of steps to determine datas.
- `Timestep::Float`: Time step of the simulation.
- `Box::Array`: Size of the Simulation Box.
- `Seed::Float`: Seed for the pseudorandom number generator.
- `Minimise::Boolean`: Sould the energy be minimize befor the Simulation.
- `TrajectoryName::String`: Data file for the trajectorys.
- `UseAngles::Boolean`: (De-)Activate angle and dihedral potentials (default=off).
- `UseCharge::Boolean`: (De-)Activate yukawa potentials (default=on).
- `Alt_GSD_Start::String`: Name of alternative GSD start file.
- `Create_Start_Config::Boolean`: Set to true, if the python scripts should create a start configuration (numerically unstable/badly equilibrated). 
- `ϵ_r::Float`: Dielectric constant.
- `κ::Float`: Debye length of electrostatic potential.  .
- `Device::String`: Device on which simulation is running (default=GPU).
- `yk_cut::Float`: Cutoff for the Yukawa potential.
- `ah_cut::Float`: Cutoff for the Ashbaugh potential.
- `ionic::Float`: Ionic strength used for Yukawa potential.
- `pH::Float`: PH value of the simulation.
- `SimType::String`: Type of the simulation (like Calvados2, Calvados3).
- `domain::Array`: Area in with the ENM is active (only for Calvados3).
- `SlabAxis::Int`: Integer representing the long slab axis.

**Creat**:
* Write a file with all parameters of the simulation that are given from the arguments.
"""
function WriteParams(filename, SimName, Temp, NSteps, NOut, Timestep, Box, Seed; Minimise=true, TrajectoryName="traj.gsd", UseAngles=true, UseCharge=true, Alt_GSD_Start="-", Create_Start_Config=false, ϵ_r=1.73136, κ=1.0, Device="GPU", yk_cut=4.0, ah_cut=2.0, ionic=0.1, pH=7.0, SimType="Calvados2",domain=Array([[0,0]]), SlabAxis=2)
    io = open(filename, "w");
    write(io, "Simname: $SimName\n")
    write(io, "Domains: $(domain)\n")
    write(io, "Seed: $Seed\n")
    write(io, "Temp: $Temp\n")
    write(io, "ionic: $ionic\n")
    write(io, "pH: $pH\n")
    write(io, "SimulationType: $SimType\n")
    write(io, "NSteps: $NSteps\n")
    write(io, "NOut: $NOut\n")
    write(io, "dt: $Timestep\n")
    write(io, "Lx: $(Box[1])\n")
    write(io, "Ly: $(Box[2])\n")
    write(io, "Lz: $(Box[3])\n")
    write(io, "Minimise: $(Minimise)\n")
    write(io, "Trajectory: $(TrajectoryName)\n")
    write(io, "UseAngles: $(UseAngles)\n")
    write(io, "UseCharge: $(UseCharge)\n")
    write(io, "Alt_GSD_Start: $(Alt_GSD_Start)\n")
    write(io, "Create_Start_Config: $(Create_Start_Config)\n")
    write(io, "epsilon_r: $(ϵ_r)\n")
    write(io, "yk_prefactor: $(138.9315360433804/ϵ_r)\n") ### e^2/(4*pi*epsilon_r) in units of KJ/mol*nm
    #write(io, "yk_D: $(D)\n")
    write(io, "kappa: $(κ)\n")
    write(io, "Device: $(Device)\n")
    write(io, "YukawaCutoff: $(yk_cut)\n")
    write(io, "AHCutoff: $(ah_cut)\n")
    write(io, "Axis: $(SlabAxis)\n")
    close(io);
end

@doc raw"""
    WriteENM_HOOMD_Indices(filename::String, ENM)

Write a data file with all necessary datas for the ENM in HOOMD.
    
**Arguments**
- `filename::String`: Path where the data will be saved.
- `ENM::Tuple`: The necessary datas.

**Creat**:
* Write a file with all datas for the ENM.
"""
function WriteENM_HOOMD_Indices(filename::String, ENM)
    if isnothing(ENM) return nothing end
    ENMB_N, ENMB_types, ENMB_typeid, ENMB_group_vector, harmonic = ENM
    io = open(filename, "w")
    write(io, "// N, ENMB_type, ENMB_typeid, ENMB_group_vector, harmonic\n")
    for (i,id) in enumerate(ENMB_typeid)
        write(io, "$(i) , $(ENMB_types[id+1]) , $(id) , $(ENMB_group_vector[i]) , $(harmonic[ENMB_types[id+1]]) \n")
    end
    close(io);
end

@doc raw"""
    writeHOOMD(Sequences,pos,image,OneToCharge,AaToId,OneToMass,OneToSigma,OneToLambda,AlphaAddition,dihedral_long_map,dihedral_eps,SimulationType,Temperature,SaltConcentration,BoxSize,StartFileName,NSteps,WriteOutFreq,Device,yk_cut,ah_cut,pH,domain,NAtoms,NBonds,NAngles,NDihedrals,dihedral_short_map,dihedral_list, ENM)

Writes the datas for the simulation in different files with HOOMD.
    
**Arguments**
- from writeStartConfiguration()

**Creates**:
* Writes the datas for the simulation.
"""
function writeHOOMD(BasePath, Sequences,pos,image,OneToCharge,AaToId,OneToMass,OneToSigma,OneToLambda,AlphaAddition,dihedral_long_map,dihedral_eps,SimulationType,Temperature,SaltConcentration,BoxSize,StartFileName,NSteps,WriteOutFreq,Device,yk_cut,ah_cut,pH,domain,NAtoms,NBonds,NAngles,NDihedrals,dihedral_short_map,dihedral_list, ENM, SlabAxis)
        mkpath("$(BasePath)/HOOMD_Setup")
        WriteHOOMDSequences("$(BasePath)/HOOMD_Setup/Sequences.txt", Sequences)
        WriteHOOMDParticlesInput("$(BasePath)/HOOMD_Setup/Particles.txt", pos,  OneToCharge, AaToId,Sequences, OneToMass, OneToSigma, image)
        if AlphaAddition
            WriteDihedrals("$(BasePath)/HOOMD_Setup/DihedralMap.txt", dihedral_long_map, dihedral_eps)
        end
        (ϵ_r, κ) = DetermineYukawaInteractions(;SimulationType=SimulationType, Temperature=Temperature, SaltConcentration=SaltConcentration)

        BoxLength = [BoxSize[2]-BoxSize[1],BoxSize[4]-BoxSize[3],BoxSize[6]-BoxSize[5] ]
        WriteParams("$(BasePath)/HOOMD_Setup/Params.txt", StartFileName, Temperature, NSteps, WriteOutFreq, 0.01, BoxLength/10.0, rand(1:65535), UseAngles=AlphaAddition;ϵ_r=ϵ_r, κ=κ,Device=Device, yk_cut=yk_cut/10.0, ah_cut=ah_cut/10.0, ionic=SaltConcentration, pH=pH, SimType=SimulationType, domain=domain,Create_Start_Config=true, SlabAxis) ### BoxLength has to be convert to nm
        WriteDictionaries("$(BasePath)/HOOMD_Setup/Dictionaries.txt", OneToCharge, AaToId, OneToMass, OneToSigma, OneToLambda)
        WriteENM_HOOMD_Indices("$(BasePath)/HOOMD_Setup/ENM_indices.txt", ENM)
        InputMasses = [OneToMass[res] for res in join(Sequences)]
        InputCharges = [OneToCharge[res] for res in join(Sequences)]

        offset = 0
        DomainRanges = []
        ### use unofficial feature "floppy boddies" to group ENMs and turn off their interactions
        groups = fill(-1, NAtoms)
        # @TODO clean up the following part
        @warn "Need to update Domain definitions for groups incase different proteins are simulated. Code assumes just one protein within the dictionary of folded proteins."
        if !isnothing(ENM)
            @assert length(keys(domain)) == 1
        
            name = first(keys(domain)) ### wrongsly assumes that the only entry is the right proteins
            for Seq in Sequences
                append!(DomainRanges, [(start+offset):(stop+offset) for (start, stop) in deepcopy(domain[name])])
                offset += length(Seq)
            end
         ### @TODO add a proper check for Calvados2
            for (i, range) in enumerate(DomainRanges)
                groups[range] .= -i-1
            end 
        end
        WriteGroups("$(BasePath)/HOOMD_Setup/Groups.txt" ,groups)

        writeGSDStartFile("$BasePath$StartFileName.gsd", NAtoms, NBonds, NAngles, NDihedrals,Float32.(BoxLength), Float32.(pos), AaToId,Sequences,image, InputMasses, InputCharges, dihedral_short_map, dihedral_list, OneToSigma, AlphaAddition, SimulationType, ENM, groups)    
end

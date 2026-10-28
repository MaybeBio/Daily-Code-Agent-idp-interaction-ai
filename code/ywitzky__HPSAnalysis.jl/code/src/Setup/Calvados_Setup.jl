@doc raw"""
    DetermineYukawaInteractions(;SimulationType="", Temperature=300, SaltConcentration=-1)

Calculates the constants for the Yukawa potencial for given temperature and salt concentration.
    
**Arguments**
- `SimulationType::String`: Type of Simulation (e.g.: "Calvados2").

**Optional Arguments**:
- `Temperature::Float`: Temperature of the Simulation.
- `SaltConcentration::Float`: Salt concentration of the Simulation.

**Return**:
* `ϵ_r::Float`: Temperature-dependent dielectric constant.
* `κ::Float`: Invers of the Debye length.
"""
function DetermineYukawaInteractions(;SimulationType="", Temperature=300, SaltConcentration=-1)
    ### constants
	e = 1.602*10.0^-(19) ### Charge of electron
	e_0 = 8.854188*10.0^-(12) ### vacuum permitivity
	kb = 1.380*10.0^-(23) ### boltzmann constant
    NA = 6.02214086 * 10.0^23 # 1/mol Avogadro constant

    ### Saltconcentration in mol/L
    if SaltConcentration==-1 || !in(SimulationType, ["Calvados2","Calvados3"])
        κ = 10.0 # Å, Screening length
        ϵ_r = 80.0 # what ever? I think unitless, relative permitivity
    else
        ### taken from RESEARCH ARTICLE Improved predictions of phase behaviour of intrinsically disordered proteins by tuning the interaction range; Giulio Tesei , Kresten Lindorff-Larsen; 2022
        RT =  8.3145*Temperature*1e-3
	    ϵ_r= 5321.0/Temperature+233.76-0.9297*Temperature+1.417*10.0^-3*Temperature^2-8.292*10.0^-7*Temperature^3 ### relative permitivity
        #λ = e^2 /(4.0*pi*e_0*ϵ_r*kb*Temperature) ### m, Bjerrum length
        λ = 1.6021766^2/(4*pi*8.854188*ϵ_r)*6.022*1000/RT ### ?, Bjerrum length

        κ = sqrt(8.0*pi*λ*SaltConcentration*6.022/10)
        #κ = sqrt(1.0/(8.0*pi*λ*NA*SaltConcentration*1000.0))*10.0^10 ### Å Screening length, named D in the paper
    end
    return (ϵ_r, κ)
end



@doc raw"""
    CalvadosSetup(Sequences,AtomTypes,pH)

Define escential Parameters for the Simulation based on the Simulation Type.
    
**Arguments**
- `Sequences::Array{String}`: List of sequences of Proteins.
- `SimulationType::String`: Type of Simulation (e.g.: "Calvados2").
- `pH::Float`: pH-value of the system.

**Return**:
A tuple containing:
- `AtomTypes::Set{Char}`: Set of unique amino acid types in the provided sequences.
- `LongAtomTypes::Set{Char}`: Set of amino acid types where the first and last amino acid in a sequence are treated as distinct types.
- `AaToId::Dict{Char, Int32}`: Dictionary mapping each amino acid type to a unique ID number.
- `IdToAa::Dict{Int32, Char}`: Dictionary mapping each ID number to its corresponding amino acid type.
- `ResToLongAtomType::Dict{Tuple{Char, Bool}, Char}`: Dictionary mapping standard amino acids to their modified forms when at the beginning or end of a sequence.
- `LongAtomTypesToRes::Dict{Char, Tuple{Char, Bool}}`: Reverse mapping of `ResToLongAtomType`.
- `OneToCharge::Dict{Char, Float}`: Dictionary containing the charge values of amino acids, modified based on simulation type.
- `OneToMass::Dict{Char, Float}`: Dictionary containing the mass values of amino acids.
- `OneToSigma::Dict{Char, Float}`: Dictionary containing the sigma values of amino acids.
- `OneToLambda::Dict{Char, Float}`: Dictionary containing the lambda values of amino acids.
- `OneToHPSDihedral0110::Dict{Char, Any}`: Dictionary containing dihedral parameters for a specific configuration.
- `OneToHPSDihedral1001::Dict{Char, Any}`: Dictionary containing dihedral parameters for another configuration.

**Notes**:
- If `SimulationType` is "Calvados2", the first and last amino acids in each sequence are assigned different types to account for altered mass due to peptide bonding. Also the charge of histidine ('H') is adjusted based on the provided pH value using the formula: `1/(1+10^(pH-6))`.
"""
function CalvadosSetup(Sequences, AtomTypes, pH, AaToId, SimulationType)
    index_cnt = length(AtomTypes)
    LongAtomTypes=deepcopy(AtomTypes)
    ResToLongAtomType = Dict()
    NewSequences= deepcopy(Sequences)
    cnt = 96 ### use lower case ascii letters as modified residue types
    ### @TODO FIX THIS, works only on first chain
    for (id,Sequence) in enumerate(Sequences)
        if( ~((Sequence[1], true) in  keys(ResToLongAtomType)))
            ResToLongAtomType[(Sequence[1], true)] = Char(cnt+=1) 
            push!(LongAtomTypes, ResToLongAtomType[(Sequence[1], true)] )#
            AaToId[ResToLongAtomType[(Sequence[1], true)] ] = index_cnt+=1
            NewSequences[id] =  ResToLongAtomType[(Sequence[1], true)]* Sequence[2:end]
        end
        if( ~((Sequence[end], false) in  keys(ResToLongAtomType)))
            ResToLongAtomType[(Sequence[end], false)] = Char(cnt+=1)
            push!(LongAtomTypes, ResToLongAtomType[(Sequence[end], false)])
            AaToId[ResToLongAtomType[(Sequence[end], false)] ] = index_cnt+=1
            NewSequences[id] = NewSequences[id][1:end-1]* ResToLongAtomType[(Sequence[end], false)]
        end
    end

    Sequences.=NewSequences
    LongAtomTypes = union(LongAtomTypes, AtomTypes)
    LongAtomTypesToRes=Dict( (v => k) for (k, v) in ResToLongAtomType)

    OneToCharge = deepcopy(BioData.OneToHPSCharge)
    OneToMass = deepcopy(BioData.AaToWeight)
    OneToSigma  = deepcopy(BioData.OneToHPSCalvadosSigma)
    OneToHPSDihedral0110 = deepcopy(BioData.OneToHPSDihedral0110)
    OneToHPSDihedral1001 = deepcopy(BioData.OneToHPSDihedral1001)
    OneToCharge['H'] = 1.0/(1+10^(pH-6))  ### HIS is pH-Dependent for Calvados2/Calvados3

    if SimulationType=="Calvados2"
        OneToLambda = deepcopy(BioData.OneToCalvados2Lambda)
    elseif SimulationType=="Calvados3"
        OneToLambda = deepcopy(BioData.OneToCalvados3Lambda)
    end

    for e in LongAtomTypes ### added the charge modifications for the first/last amino acids
        if ~(e in keys(OneToCharge))
            (AA, front) = LongAtomTypesToRes[e]
            OneToCharge[e] = front ? OneToCharge[AA] +1 :  OneToCharge[AA] -1
            OneToMass[e] = front ? OneToMass[AA] +2.0 :  OneToMass[AA] +16
            OneToSigma[e] = OneToSigma[AA]
            OneToLambda[e] = OneToLambda[AA]
            OneToHPSDihedral0110[e] = OneToHPSDihedral0110[AA]
            OneToHPSDihedral1001[e] = OneToHPSDihedral1001[AA]
        end
    end
    return (AtomTypes, LongAtomTypes, AaToId,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, OneToHPSDihedral0110, OneToHPSDihedral1001)
end


@doc raw"""
    DetermineAtomTypes(Sequences, SimulationType, pH; OneToChargeDef=BioData.OneToHPSCharge, OneToLambdaDef=BioData.OneToCalvados2Lambda, OneToSigmaDef=BioData.OneToHPSCalvadosSigma)

Defines essential Parameters for the simulation based on the type of the simulation.
    
**Arguments**
- `Sequences::Array{String}`: List of sequences of proteins.
- `SimulationType::String`: Type of simulation (e.g.: "Calvados2").
- `pH::Float`: pH-value of the system.

**Optional Arguments**:
- `OneToChargeDef::Dict()`: Dictionary defining the charge for the Aminoacids.
- `OneToLambdaDef::Dict()`: Dictionary defining the Lambda for the Aminoacids.
- `OneToSigmaDef::Dict()`: Dictionary defining the Sigma for the Aminoacids.

**Return**:
A tuple containing:
- `AtomTypes::Set{Char}`: Set of unique amino acid types in the provided sequences.
- `LongAtomTypes::Set{Char}`: Set of amino acid types where the first and last amino acid in a sequence are treated as distinct types.
- `AaToId::Dict{Char, Int32}`: Dictionary mapping each amino acid type to a unique ID number.
- `IdToAa::Dict{Int32, Char}`: Dictionary mapping each ID number to its corresponding amino acid type.
- `ResToLongAtomType::Dict{Tuple{Char, Bool}, Char}`: Dictionary mapping standard amino acids to their modified forms when at the beginning or end of a sequence.
- `LongAtomTypesToRes::Dict{Char, Tuple{Char, Bool}}`: Reverse mapping of `ResToLongAtomType`.
- `OneToCharge::Dict{Char, Float}`: Dictionary containing the charge values of amino acids, modified based on simulation type.
- `OneToMass::Dict{Char, Float}`: Dictionary containing the mass values of amino acids.
- `OneToSigma::Dict{Char, Float}`: Dictionary containing the sigma values of amino acids.
- `OneToLambda::Dict{Char, Float}`: Dictionary containing the lambda values of amino acids.
- `OneToHPSDihedral0110::Dict{Char, Any}`: Dictionary containing dihedral parameters for a specific configuration.
- `OneToHPSDihedral1001::Dict{Char, Any}`: Dictionary containing dihedral parameters for another configuration.

**Notes**:
- If `SimulationType` is "Calvados2", the first and last amino acids in each sequence are assigned different types to account for altered mass due to peptide bonding. Also the charge of histidine ('H') is adjusted based on the provided pH value using the formula: `1/(1+10^(pH-6))`.
- If `SimulationType` is "HPS-Alpha", predefined values for charge, lambda, and sigma are used.
- If an unknown `SimulationType` is provided, the function falls back on the user-supplied dictionaries for charge, lambda, and sigma values.
"""
function DetermineAtomTypes(Sequences, SimulationType, pH; OneToChargeDef=BioData.OneToHPSCharge, OneToLambdaDef=BioData.OneToCalvados2Lambda, OneToSigmaDef=BioData.OneToHPSCalvadosSigma)
    #Define Aminoacids to ID Dict
    AtomTypes= Set(join(Sequences))
    AaToId = Dict{Char,Int32}()
    for (index, value) in enumerate(AtomTypes)
        AaToId[value]=index
    end

    LongAtomTypes=deepcopy(AtomTypes)
    ResToLongAtomType = Dict()
    LongAtomTypesToRes=Dict()
    OneToMass = deepcopy(BioData.AaToWeight)
    OneToHPSDihedral0110 = deepcopy(BioData.OneToHPSDihedral0110)
    OneToHPSDihedral1001 = deepcopy(BioData.OneToHPSDihedral1001)
    if SimulationType=="HPS-Alpha"
        OneToCharge = deepcopy(BioData.OneToHPSCharge)
        OneToLambda = deepcopy(BioData.OneToHPSUrryLambda)
        OneToSigma  = deepcopy(BioData.OneToHPSCalvadosSigma)
    elseif SimulationType=="Calvados2"||SimulationType=="Calvados3"
        (AtomTypes, LongAtomTypes, AaToId,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, OneToHPSDihedral0110, OneToHPSDihedral1001)=CalvadosSetup(Sequences, AtomTypes, pH, AaToId, SimulationType)
    else ### take the ones which are supplied
        OneToCharge = deepcopy(OneToChargeDef)
        OneToLambda = deepcopy(OneToLambdaDef)
        OneToSigma = deepcopy(OneToSigmaDef)
    end
    IdToAa=Dict( (v => k) for (k, v) in AaToId)

    return (AtomTypes, LongAtomTypes, AaToId, IdToAa,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, OneToHPSDihedral0110, OneToHPSDihedral1001)
end


function startConfigurationSetup(Sequences,SimulationType,pH,OneToChargeDef,OneToLambdaDef,OneToSigmaDef,MixingRule="1-1001-1")
    #Define Number of all Aminoacids, Bonds, Angles and Dihedrals, set AlphaAddition
    NAtoms=0
    NBonds=0
    NAngles=0
    NDihedrals=0
    for Seq in Sequences
        NAtoms += length(Seq)
        NBonds += length(Seq)-1
        NAngles += length(Seq)-2
        NDihedrals += length(Seq)-3
    end

    AlphaAddition=false
    if SimulationType=="Calvados2+Alpha"
        AlphaAddition = true
        SimulationType="Calvados2"
    elseif SimulationType=="Calvados3+Alpha"
        AlphaAddition = true
        SimulationType="Calvados3"
    else
        AlphaAddition=false
    end

    (AtomTypes, LongAtomTypes, AaToId, IdToAa,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, OneToHPSDihedral0110, OneToHPSDihedral1001) =  DetermineAtomTypes(Sequences, SimulationType, pH; OneToChargeDef=OneToChargeDef, OneToLambdaDef=OneToLambdaDef, OneToSigmaDef=OneToSigmaDef)
    #Define length of all chains
    NAtomTypes = length(LongAtomTypes)

    #if AlphaAddition then determine the Dihedral
    dihedral_short_map=Dict()
    dihedral_list = zeros(Int32, (0,0))
    dihedral_long_map=Dict()
    dihedral_eps=[]
    if AlphaAddition
        (dihedral_short_map, dihedral_long_map, dihedral_eps, dihedral_list) = determineDihedrals(Sequences, AtomTypes, AaToId, OneToHPSDihedral0110, OneToHPSDihedral1001, MixingRule)
        NDihedralsTypes = length(dihedral_eps)
    end
    return NAtoms,NBonds,NAngles,NDihedrals,AlphaAddition,SimulationType,AtomTypes, LongAtomTypes, AaToId, IdToAa,ResToLongAtomType, LongAtomTypesToRes, OneToCharge, OneToMass, OneToSigma, OneToLambda, NAtomTypes, dihedral_short_map, dihedral_long_map, dihedral_eps, dihedral_list
end


function determineDihedrals(Sequences, Types, TypeToId, OneToHPSDihedral0110, OneToHPSDihedral1001, MixingRule="1-1001-1")
    dihedral_short_map = Dict()
    dihedral_long_map = Dict()

    dihedral_eps = ones(400)*Inf
    dihedral_cnt=0
    eps=0.
    AA_ind = 1
    AA2_ind = 1

    ReverseDihedralMap = Dict()
    key=()
    M = MixingRule=="1-1001-1" ? 4 : 2
    dihedral_list = zeros(Int32, sum(length.(Sequences).-3), M)
    index = 1
    for Sequence in Sequences
        SeqLength = length(Sequence)
        for (ind, AA) in enumerate(Sequence)
            if ind> (SeqLength-3) 
                continue
            else
                if MixingRule=="1001"
                    AA2 = Sequence[ind+3]
                    AA_ind = TypeToId[AA]
                    AA2_ind = TypeToId[AA2]
                    eps = (OneToHPSDihedral1001[AA]+OneToHPSDihedral1001[AA2])/2.
                    key = (min(AA_ind,AA2_ind), max(AA_ind,AA2_ind))
                elseif MixingRule=="1-1001-1"
                    if (ind ==1 )
                        AAminus = "Z"
                        AA_min_ind = 0
                    else
                        AAminus = Sequence[ind-1]
                        AA_min_ind = TypeToId[AAminus]
                    end
                    if ( ind <=(SeqLength-4) ) 
                        AA3= Sequence[ind+4]
                        AA3_ind = TypeToId[AA3]
                    else
                        AA3="X"
                        AA3_ind =-1
                    end
                    AA2 = Sequence[ind+3]
                    AA_ind = TypeToId[AA]
                    AA2_ind = TypeToId[AA2]
                    if (ind ==1)
                        eps = (OneToHPSDihedral1001[AA]+OneToHPSDihedral1001[AA2]+OneToHPSDihedral1001[AA3])/3.
                    elseif ( ind >(SeqLength-4) ) 
                        eps = (OneToHPSDihedral1001[AA]+OneToHPSDihedral1001[AA2]+OneToHPSDihedral1001[AAminus])/3.
                    else
                        eps = (OneToHPSDihedral1001[AA]+OneToHPSDihedral1001[AA2]+OneToHPSDihedral1001[AAminus]+OneToHPSDihedral1001[AA3])/4.
                    end
                    key = (sort([AA_min_ind,AA_ind, AA2_ind, AA3_ind]))
                elseif MixingRule=="0110"
                    AA1 = Sequence[ind+1]
                    AA2 = Sequence[ind+2]
                    AA_ind = TypeToId[AA1]
                    AA2_ind = TypeToId[AA2]
                    eps = (OneToHPSDihedral0110[AA1]+OneToHPSDihedral0110[AA2])/2.
                    key = (min(AA_ind,AA2_ind), max(AA_ind,AA2_ind))
                end
                if !haskey(dihedral_short_map, key)
                    eps_short = round(eps;digits=4)
                    same_eps_ind = findfirst(x->x==eps_short, dihedral_eps)
                    if same_eps_ind===nothing
                        dihedral_cnt+=1
                        dihedral_short_map[key]=dihedral_cnt
                        dihedral_long_map[key]=dihedral_cnt

                        ReverseDihedralMap[eps_short] = key
                        dihedral_eps[dihedral_cnt] = eps_short
                    else
                        dihedral_long_map[key]=same_eps_ind
                        key = ReverseDihedralMap[eps_short]
                    end
                end
                dihedral_list[index, :] .= key 
                index += 1
            end
        end
    end
    dihedral_eps = dihedral_eps[1:dihedral_cnt]
    return (dihedral_short_map, dihedral_long_map, dihedral_eps, dihedral_list)
end

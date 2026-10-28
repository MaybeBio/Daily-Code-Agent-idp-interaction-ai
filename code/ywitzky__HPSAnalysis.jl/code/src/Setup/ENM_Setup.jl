
@doc raw"""
    BuildENMModel(Sim::HPSAnalysis.SimData{T,I}, DomainDict, Proteins, Sequences, ProteinJSON) where {T<:Real, I<:Integer}

Calculate the Indices, that are nessesary to creat a start file im HOOMD.
    
**Arguments**
- `Sim::HPSAnalysis.SimData{T,I}`: The simulation datas.
- `DomainDict`: The Domains in which the ENM is active.
- `Proteins`: List of Protein Names.
- `Sequences`: The Sequences of the Proteins.
- `ProteinJSON`: AlphaFold data of the Proteins.

**Return**:
- Number of bonds
- Vector of bond type names 
- ID connecting bond tuple to bond type
- Vector of tuples defining all bonds
- Dict{String, Dict{Symbol, Float64}} which defines the bonds as used in HOOMD.
"""
function BuildENMModel(Sim::HPSAnalysis.SimData{T,I}, DomainDict, Proteins, Sequences, ProteinJSON; plDDTcut=90.0, pae_cut=1.85, precision=1000) where {T<:Real, I<:Integer} 
    ConstraintDict, BackboneCorrectionDict = DetermineCalvados3ENMfromAlphaFold(Sim.BasePath, DomainDict, Proteins, ProteinJSON; BBProtein="CA", rcut = 9.0, plDDTcut=plDDTcut, pae_cut=pae_cut, precision=precision)
    HOOMD_Indices = ComputeHOOMD_ENMIndices(ConstraintDict, BackboneCorrectionDict, Sequences, Proteins)
    UnfoldedRegions =  GenerateUnfoldedRegions(Proteins, DomainDict, Sequences)
    All_Indices = CombineBackboneAndENM(Proteins, Sequences, HOOMD_Indices, UnfoldedRegions, BackboneCorrectionDict)

    return All_Indices
end

@doc raw"""
    CombineBackboneAndENM(Proteins, Sequences, HOOMD_Indices, UnfoldedRegions, BackboneCorrectionDict)

Combines the bond length potential of backbone and ENM.
    
**Arguments**
- `Proteins`: List of Protein Names.
- `Sequences`: The Sequences of the Proteins.
- `HOOMD_Indices`: Tuple containing the info of ENM as defined by this functions return.
- `UnfoldedRegions`: Dictionary defining the unfolded regions.
- `BackboneCorrectionDict`::Dictionary with atoms, bond length and id which will be connected via backbone.

**Return**:
- Number of bonds
- Vector of bond type names 
- ID connecting bond tuple to bond type
- Vector of tuples defining all bonds
- Dict{String, Dict{Symbol, Float64}} which defines the bonds as used in HOOMD.
"""
function CombineBackboneAndENM(Proteins, Sequences, HOOMD_Indices, UnfoldedRegions, BackboneCorrectionDict)
    offsets = vcat([0],cumsum(length.(Sequences)))

    BackboneDiffLengths = Dict([prot=> [(i,j) for (i,j,_) in BackboneCorrectionDict[prot]] for prot in Set(Proteins)]) 

    #N_Bonds = sum([sum([b-a for (a,b) in UnfoldedRegions[prot]]) for prot in Proteins])
    #BB_ID = zeros(Int32, N_Bonds)
    BB_groups = Vector{Tuple{Int32, Int32}}()

    (ENM_Bonds, ENM_types, ENM_typeid, ENM_groups, harmonic) = HOOMD_Indices
    ### ENM_types, harmonic contains [O-O] already
    for (I, prot) in enumerate(Proteins)
        for Domain in UnfoldedRegions[prot]
            for i in range(Domain[1], Domain[2]-1)
                if !( (i,i+1) in BackboneDiffLengths[prot])
                    push!(BB_groups, (offsets[I]+i-1,offsets[I]+i)) ## shift from 1 -> 0, because of julia -> python
                end
            end
        end
    end
    N_Bonds = length(BB_groups)
    BB_ID = zeros(Int32, N_Bonds)

    return (N_Bonds+ENM_Bonds, ENM_types, vcat(BB_ID, ENM_typeid), vcat(BB_groups, ENM_groups), harmonic )
end

@doc raw"""
    ComputeHOOMD_ENMIndices(ConstraintDict, BackboneCorrectionDict, Sequences, Proteins)

Calculate the Indices, that are nessesary to creat a start file im HOOMD, .
    
**Arguments**
- `ConstraintDict::Dictionary with atoms, bond length and id which will be connected via ENM.
- `BackboneCorrectionDict`::Dictionary with atoms, bond length and id which will be connected via backbone.
- `Sequences`: The Sequences of the Proteins.
- `Proteins`: List of Protein Names.

**Return**:
- Number of bonds
- Vector of bond type names 
- ID connecting bond tuple to bond type
- Vector of tuples defining all bonds
- Dict{String, Dict{Symbol, Float64}} which defines the bonds as used in HOOMD.
"""
function ComputeHOOMD_ENMIndices(ConstraintDict, BackboneCorrectionDict , Sequences, Proteins)
    offsets = cumsum(length.(Sequences))

    bondLength = 0.38 # in nm
    BB_N = 0
    ENM_N = 0
    BB_ID = Int[]
    ENM_ID = Int[]
    BB_groups = Vector{Tuple{Int32, Int32}}()
    ENM_groups = Vector{Tuple{Int32, Int32}}()
    harmonic = Dict{String, Dict{Symbol, Float64}}()
    harmonic["O-O"] = Dict(:r => bondLength, :k => 8033.0) # default backbone bond values
    off = 0

    for (i, Prot) in enumerate(Proteins)
        for (atom_1, atom_2, r0, ind) in BackboneCorrectionDict[Prot]
            if r0 != bondLength
                push!(BB_ID, ind) 
            elseif r0 == bondLength
                push!(BB_ID, 0)
            end
            push!(BB_groups, (atom_1 -1 + off, atom_2 -1 + off)) ## shift from 1 -> 0, because of julia -> python
            BB_N+= 1
        end
        off = offsets[i]
    end
    off = 0
    for (i, Prot) in enumerate(Proteins)
        for (atom_1, atom_2, r0, ind) in ConstraintDict[Prot]
            push!(ENM_ID, ind) 
            push!(ENM_groups, (atom_1 -1 + off, atom_2 -1 + off))
            ENM_N+= 1
        end
        off = offsets[i]
    end
    BB_tmp  = vcat([BackboneCorrectionDict[prot] for prot in Set(Proteins)]...)
    ENM_tmp = vcat([ConstraintDict[prot]         for prot in Set(Proteins)]...)

    BB_types  = ["BB_$(ind)"  for ind in sort(collect(Set(getindex.(BB_tmp ,4)))) if ind != 0]
    ENM_types = ["ENM_$(ind)" for ind in sort(collect(Set(getindex.(ENM_tmp,4))))]

    for (r,ind) in zip(getindex.(BB_tmp,3) , getindex.(BB_tmp,4))
        if ind!=0 ### "BB_0" is the same as "O-O"
            harmonic["BB_$(ind)"] = Dict(:r => r, :k => 8033)
        end
    end
    for (r,ind) in zip(getindex.(ENM_tmp,3), getindex.(ENM_tmp,4))
        harmonic["ENM_$(ind)"] = Dict(:r => r, :k => 700)
    end
    
    ### sort list according to types indices
    merged = vcat(BB_types, ENM_types)
    inds = sortperm(parse.(Int32,getindex.(split.(merged, "_"),2)))

    return (BB_N + ENM_N, vcat(["O-O"], merged[inds]), vcat(BB_ID, ENM_ID), vcat(BB_groups, ENM_groups), harmonic)
end

@doc raw"""
    DetermineCalvados3ENMfromAlphaFold(BasePath, DomainDict, Proteins, ProteinJSON; BBProtein="CA", rcut = 9.0, plDDTcut=90.0, pae_cut=1.85)
Return a dictionary of atoms and there distances that are nessesary for the Elastic Network Model, and a dictionary for the correction of the backbone bonds in the ENM domain.
    
**Arguments**
- `BasePath`: The base path of the simulation setup procedure.
- `DomainDict`:: The Domains in which the ENM is active.
- `Proteins`: List of Protein Names.
- `ProteinJSON`: Dictionary of AlphaFold data of the Proteins in JSON format.
- `BBProtein`: The atom from which the AlphaFold datas are set for the aminoacid.
- `rcut`: Cut of length for the ENM in Angstroem.
- `plDDTcut`: Cut of plDDT parameter of AlphaFold reference for the ENM.
- `pae_cut`: Cut of pae parameter of AlphaFold reference for the ENM.

**Return**:
* `ConstraintDict`: Dictionary that maps Protein names to a Vector containing Tuples of Indices of i,j and distance r which define additional bonds necessary for ENM.
* `BackboneCorrectionDict`: Dictionary that maps Protein names to a Vector containing Tuples of Indices of i,j and distance r which define the backbone bonds that may have different lengths in ENM regions.
"""
function DetermineCalvados3ENMfromAlphaFold(BasePath::String, DomainDict, Proteins, ProteinJSON; BBProtein="CA", rcut = 9.0, plDDTcut=90.0, pae_cut=1.85, precision=1000)
    ## distances in nm
    ciffolder = "$(BasePath)/InitFiles/CifFiles"
    ConstraintDict = Dict{String, Vector{Tuple{Int,Int, Float64, Int32}}}() ### Contains the bonds from ENM
    BackboneCorrectionDict = Dict{String, Vector{Tuple{Int,Int, Float64, Int32}}}() ### Contains the bonds from backbone with different lengths in folded regions
    bondLength = 0.38  ## coordinates are in nm

    ExistingConstraints=Dict{Float64, Int32}(bondLength=>1) ### maps existing constraint determined by distance onto integer
    ExistingBackbone   =Dict{Float64, Int32}(bondLength=>1) ### maps existing backbone   determined by distance onto integer
    cnt=1
    Masses = Dict("C"=>12.011, "N"=>14.007, "O"=>15.999, "S"=>32.065, "H"=>1.0078)
    for Prot in Set(Proteins)
        BackboneCorrectionDict[Prot] = []
        ConstraintDict[Prot] = []
        if length(DomainDict[Prot])>0
            CifPath = "$(ciffolder)/$(Prot).cif"
            lines = readlines(CifPath)
            x = zeros(length(lines))
            y = zeros(length(lines))
            z = zeros(length(lines))
            norm = zeros(length(lines))
            plDDT = zeros(length(lines))

            step = 0
            for line in readlines(CifPath)
                fields = strip.(split(line))
                if !isempty(fields) && fields[1] == "ATOM" #&& fields[4] == BBProtein
                    if fields[4] == "N"
                        step += 1
                    end
                    mass = Masses[fields[3]]
                    x[step] += parse(Float64, fields[11])*mass
                    y[step] += parse(Float64, fields[12])*mass
                    z[step] += parse(Float64, fields[13])*mass
                    plDDT[step]+= parse(Float64, fields[15])*mass
                    norm[step] += mass
                end
            end
            x = x[1:step]./norm[1:step]
            y = y[1:step]./norm[1:step]
            z = z[1:step]./norm[1:step]
            plDDT[1:step] = plDDT[1:step]./norm[1:step]
            pae =JSON.parsefile(ProteinJSON[Prot])["pae"]
            for i in 1:step-1 # up to NAtom-1
                #in_any_domain = false
                for Domain in DomainDict[Prot]
                    if Domain[1] ≤ i ≤ Domain[2] && i + 1 ≤ Domain[2]
                        dist = sqrt((x[i+1]-x[i])^2 + (y[i+1]-y[i])^2 + (z[i+1]-z[i])^2)
                        dist /= 10.0 ### convert to nm
                        dist = round(dist*precision)/precision ### reduces distance in bond length to .3f precision
                        if !(dist in keys(ExistingBackbone))
                            ExistingBackbone[dist]=cnt
                            cnt+=1
                        end
                        push!(BackboneCorrectionDict[Prot], (i, i+1,dist, ExistingBackbone[dist]))

                        #in_any_domain = true
                        if plDDT[i] ≥ plDDTcut
                            for j in i+3:Domain[2]
                                dist_sqr = (x[j]-x[i])^2 + (y[j]-y[i])^2 + (z[j]-z[i])^2
                                if dist_sqr < rcut^2 && plDDT[j] ≥ plDDTcut && pae[i][j] < pae_cut && pae[j][i] < pae_cut
                                    dist=sqrt(dist_sqr)
                                    dist /= 10.0 ### convert to nm
                                    dist = round(dist*precision)/precision ### reduces distance in bond length to .3f precision
                                    if !(dist in keys(ExistingConstraints))
                                        ExistingConstraints[dist]=cnt
                                        cnt+=1
                                    end
                                    push!(ConstraintDict[Prot], (i, j, dist, ExistingConstraints[dist]))
                                end
                            end
                        end
                        break
                    end
                end
            
                #if !in_any_domain
                #    push!(BackboneCorrectionDict[Prot], (i, i+1, bondLength))
                #end
            end
        end
    end
    return ConstraintDict, BackboneCorrectionDict
end

@doc raw"""
    GenerateUnfoldedRegions(Proteins, DomainDict, Sequences)

Generates a dictionary that defines the unfolded regions of proteins based on the definition of folded regions.
    
**Arguments**
- `Proteins`: List of Protein Names.
- `DomainDict`:: The Domains in which the ENM is active.
- `Sequences`: The Sequences of the Proteins.

**Return**:
- UnfoldedDict:: Dict{String, Vector{Tuple{Int64, Int64}}}(): Maps Protein names onto vectors of tuples that define the regions of unfolded regions.
"""
function GenerateUnfoldedRegions(Proteins, DomainDict, Sequences)
    ProtLength = Dict([prot=>length(seq) for (seq, prot) in zip(Sequences, Proteins)])
    UnfoldedDict= Dict{String, Vector{Tuple{Int64, Int64}}}()

    for Prot in Set(Proteins)
        UnfoldedDomains = []
        N = ProtLength[Prot]
        if length(DomainDict[Prot])>0
            FoldedDomains = sort(DomainDict[Prot])
            if FoldedDomains[1][1]!=1
                push!(UnfoldedDomains, (1, FoldedDomains[1][1]))
            end
            for (i, val) in enumerate(FoldedDomains[1:end-1])
                push!(UnfoldedDomains, (val[2], FoldedDomains[i+1][1]))
            end
            if FoldedDomains[end][2]!=N
                push!(UnfoldedDomains, (FoldedDomains[end][2], N))
            end
        else 
            push!(UnfoldedDomains, (1, N))
        end
        UnfoldedDict[Prot] = UnfoldedDomains
    end
    return UnfoldedDict
end
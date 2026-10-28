function allocateBondAngleVecs(Sim::SimData{T,I}) where {T<:Real, I<:Integer}
    Δx = zeros(T, (Sim.NAtoms-1))
    Δy = zeros(T,(Sim.NAtoms-1))
    Δz = zeros(T,(Sim.NAtoms-1))
    length_vec = zeros(T, (Sim.NAtoms-1))
    return Δx, Δy,Δz,length_vec
end

function computeBondAngles!(Sim::SimData{T, I}, Δx, Δy,Δz,length_vec ; step=1) where {T<:Real, I<:Integer}
    Δx .= @views Sim.x_uw[2:end, step] .-  Sim.x_uw[1:end-1, step]
    Δy .= @views Sim.y_uw[2:end, step] .-  Sim.y_uw[1:end-1, step]
    Δz .= @views Sim.z_uw[2:end, step] .-  Sim.z_uw[1:end-1, step]
    
    length_vec  .=  @. sqrt(Δx^2+Δy^2+Δz^2)

    Sim.CosBondAngles  .= @views (Δx[1:end-1] .* Δx[2:end] .+ Δy[1:end-1].*Δy[2:end]  .+ Δz[1:end-1].*Δz[2:end] )./ (length_vec[1:end-1].*length_vec[2:end])

    Sim.CosBondAngles .-= 10^-7*sign.(Sim.CosBondAngles) ### Kyras error correction, such that "DomainError with 1.0000001: acos(x) not defined for |x| > 1" doesnt occure

    Sim.BondAngles .=  @. T(π) - acos(-Sim.CosBondAngles)
    return nothing
end

@doc raw"""
    computeAvgBondAngles(Sim::SimData{T, I}; NBins = 180) where {T<:Real, I<:Integer}

Computes the average bond angles.

Computes the bond angles, the cosine thereof for each triplet of beads i, i+1, i+2 along the backbone of the proteins for each frame. The average bond angle and the average cosine of the bond angle are computed over all steps for each triplet independly. The histogram averages also over all triplets. The persistence length is computed if the average of the cosine of the bond angle is large than zero and otherwise set to 0. 

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.
- `NBins`::Int: Number of bins for the angle histogram (default: 180).

**Creates**:
* `Sim.AvgBondAngles`: Stores the averaged bond angles.
* `Sim.AvgCosBondAngles`: Stores the averaged cosine of bond angles.
* `Sim.LocalPersistenceLength`: Computes the persistence length from bond angles.
* `Sim.BondAngleHist`: Stores the histogram of bond angles.
"""
function computeAvgBondAngles(Sim::SimData{T, I}; NBins = 180) where {T<:Real, I<:Integer}
    r_0=T(-3.8)

    Sim.BondAngles=zeros(T, (Sim.NAtoms-2) )
    Sim.CosBondAngles=zeros(T, (Sim.NAtoms-2) )
    Sim.AvgBondAngles=zeros(T, (Sim.NAtoms-2))
    Sim.AvgCosBondAngles=zeros(T, (Sim.NAtoms-2) )
    Sim.LocalPersistenceLength=zeros(T, (Sim.NAtoms-2) )
    Sim.BondAngleHist=zeros(T , NBins)

    InvAngleResolution = T(NBins)/π

    Δx, Δy,Δz,length_vec = allocateBondAngleVecs(Sim)
    for step in Sim.EquilibrationTime:Sim.NSteps
        #Calculation of the Bond Angles for each Step, with error
        computeBondAngles!(Sim, Δx, Δy,Δz,length_vec; step=step)

        Sim.AvgBondAngles .+= Sim.BondAngles
        Sim.AvgCosBondAngles .+=Sim.CosBondAngles

        ### add to histogram
        for chain in 1:Sim.NChains
            for atom in Sim.ChainStart[chain]:Sim.ChainStop[chain]-2
                #Round up to next Integer (I)
                index = ceil(I, Sim.BondAngles[atom]*InvAngleResolution)
                if index != 0 
                    Sim.BondAngleHist[index]+= 1.0
                end
            end
        end
    end
    ### normalize avg
    Sim.AvgBondAngles ./= T(Sim.NSteps-Sim.EquilibrationTime)
    Sim.AvgCosBondAngles ./= T(Sim.NSteps-Sim.EquilibrationTime)
    #if val>T(0.0) true then r_0/log(val) if this is false then T(0.0)
    Sim.LocalPersistenceLength = [ val>T(0.0) ? r_0/log(val) : T(0.0) for val in Sim.AvgCosBondAngles ]

    ### normalize histogram
    NAngles = 0
    for chain in 1:Sim.NChains
        NAngles +=Sim.ChainStop[chain]-Sim.ChainStart[chain]-2
    end
    #conversion into Real -> invers
    inv = 1.0/convert(T, (Sim.NSteps-Sim.EquilibrationTime+1)*NAngles/InvAngleResolution)
    Sim.BondAngleHist .*= inv

end

function computeLocalPersistance(Sim::SimData{R,I}, r_0=1.0, MaxVariants = 5) where {R<:Real, I<:Integer}

    MaxLength =minimum(Sim.ChainLength)÷2
    NVariants  = min(MaxLength, MaxVariants)
    BondAngles=zeros(eltype(Sim.x),( Sim.NAtoms-2, Sim.NSteps))
    CosBondAngles=zeros(eltype(Sim.x),( Sim.NAtoms-2, Sim.NSteps))
    AvgBondAngles=zeros(eltype(Sim.x),( Sim.NAtoms-2))
    AvgCosBondAngles=zeros(eltype(Sim.x),( Sim.NAtoms-2))
    Sim.VarLocalPersistenceLength=zeros(eltype(Sim.x), ( Sim.NAtoms-2,  min(minimum(Sim.ChainLength)÷2, MaxVariants)))

    vec = zeros(eltype(Sim.x),(3,Sim.NAtoms-1))
    length_vec = zeros(eltype(Sim.x), (Sim.NAtoms-1))
    dot_product = zeros(eltype(Sim.x),(Sim.NAtoms-2))
    for step_size in 1:NVariants
        #AvgBondAngles .=0 
        AvgCosBondAngles .*= 0.0
        for step in 1:Sim.NSteps
            #compute connecting vectors 
            for i in 1:(Sim.NAtoms-step_size)
                vec[1,i]=Sim.x_uw[i+step_size, step]-Sim.x_uw[i,step]
                vec[2,i]=Sim.y_uw[i+step_size, step]-Sim.y_uw[i,step]
                vec[3,i]=Sim.z_uw[i+step_size, step]-Sim.z_uw[i,step]
                length_vec[i] = sqrt(vec[1,i]^2+vec[2,i]^2+vec[3,i]^2)
            end

            ### compute angles 
            for i in 1:(Sim.NAtoms-step_size-1)
                dot_product[i] = vec[1,i]*vec[1,i+1] + vec[2,i]*vec[2,i+1] + vec[3,i]*vec[3,i+1]
                CosBondAngles[i,step]= dot_product[i]/(length_vec[i]*length_vec[i+1])
                if abs(CosBondAngles[i,step])>1
                    if CosBondAngles[i,step]<1.0001 && CosBondAngles[i,step]>1.0
                        CosBondAngles[i,step]=0.99999
                    end
                    if CosBondAngles[i,step]>-1.0001 && CosBondAngles[i,step]<-1.0
                        CosBondAngles[i,step]=-0.9999
                    end
                end

                #[i, step] = pi-acos(-1.0*Sim.CosBondAngles[i,step]) ### the first vector points in the wrong direction -> angles is 180-angle
                #AvgBondAngles[i]+=Sim.BondAngles[i, step]
                AvgCosBondAngles[i]+=CosBondAngles[i, step]
            end
        end
        #AvgBondAngles ./= Sim.NSteps
        AvgCosBondAngles ./= Sim.NSteps
        for i in 1:(Sim.NAtoms-step_size-1)
            if(AvgCosBondAngles[i]>0.0)
                Sim.VarLocalPersistenceLength[i,step_size] = -r_0/log(AvgCosBondAngles[i])
            else
                Sim.VarLocalPersistenceLength[i, step_size]=0.0
            end
        end
    end

end

@doc raw"""
    computeDihedralAngles(Sim::SimData{R,I}) where {R<:Real, I<:Integer}

Computes the dihedral angles (torsion angles).

**Arguments**:
 - `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Creates**:
 * Sim.TorsionAngles with computed dihedral angles.
"""
function computeDihedralAngles(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    ### also computes the improper dihedral angles for residues between chains !!!!
    Sim.TorsionAngles=zeros(( Sim.NAtoms-3, Sim.NSteps))

    vec = zeros(eltype(Sim.x),(3,Sim.NAtoms-1))
    length_vec = zeros(eltype(Sim.x), (Sim.NAtoms-1))
    multlength=zeros(eltype(Sim.x),(Sim.NAtoms-1))
    cross = zeros(eltype(Sim.x),(3,Sim.NAtoms-2))

    vec = zeros(Float64,(3,Sim.NAtoms-1))
    length_vec = zeros(Float64, (Sim.NAtoms-1))
    multlength=zeros(Float64,(Sim.NAtoms-1))
    cross = zeros(Float64,(3,Sim.NAtoms-2))

    for step in 1:Sim.NSteps
        #compute connecting vectors 
        for i in 1:Sim.NAtoms-1
            vec[1,i]=Sim.x_uw[i+1, step]-Sim.x_uw[i,step]
            vec[2,i]=Sim.y_uw[i+1, step]-Sim.y_uw[i,step]
            vec[3,i]=Sim.z_uw[i+1, step]-Sim.z_uw[i,step]
            length_vec[i] = sqrt(vec[1,i]^2+vec[2,i]^2+vec[3,i]^2)
        end

        for i in 2:Sim.NAtoms-1 
            cross[1,i-1]=(vec[2,i-1]*vec[3,i]-vec[3,i-1]*vec[2,i])
            cross[2,i-1]=(vec[3,i-1]*vec[1,i]-vec[1,i-1]*vec[3,i])
            cross[3,i-1]=(vec[1,i-1]*vec[2,i]-vec[2,i-1]*vec[1,i])
        end

        #compute torsion angle https://en.wikipedia.org/wiki/Dihedral_angle
        for i in 1:Sim.NAtoms-3 
            arg1 = length_vec[i+1]*(vec[1,i]*cross[1,i+1]+vec[2,i]*cross[2,i+1]+vec[3,i]*cross[3,i+1])
            arg2 = cross[1,i]*cross[1,i+1]+cross[2,i]*cross[2,i+1]+cross[3,i]*cross[3,i+1]

            Sim.TorsionAngles[i, step] = atan(arg1,arg2 )
        end
    end
end

@doc raw"""
    computeDihedralHist(Sim::SimData{R,I}; N = 500) where {R<:Real, I<:Integer}

Computes a histogram of dihedral angles.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Creates**:
* Sim.TorsionHist with the invert of NDihedrals.
"""
function computeDihedralHist(Sim::SimData{R,I}; N = 500) where {R<:Real, I<:Integer}
    #Sim.TorsionHist=OffsetArray(zeros(R , 1+Sim.NAtomTypes, 360),1:1+Sim.NAtomTypes, -179:180 ) 
    Sim.TorsionHist=OffsetArray(zeros(R , 1, 360),1:1, -179:180 ) 

    AllInd = 1+Sim.NAtomTypes
    InvAngleResolution = 360.0/(2.0*π)
    for step in Sim.EquilibrationTime:Sim.NSteps #ceil(I,Sim.NSteps/)
        for chain in 1:Sim.NChains
            for atom in Sim.ChainStart[chain]:Sim.ChainStop[chain]-3
                index =  ceil(I, Sim.TorsionAngles[atom, step]*InvAngleResolution)
                if index>180
                    index=180
                end
                if index<-179
                    index=-179
                end
                #Sim.TorsionHist[Sim.IDs[atom], index]+= 1
                Sim.TorsionHist[1, index]+= 1.0
            end
        end
    end
    #=
    for type in 1:Sim.NAtomTypes
        for index in -179:180
            Sim.TorsionHist[AllInd, index] += Sim.TorsionHist[type, index] 
        end
    end=#

    NDihedrals = 0
    for chain in 1:Sim.NChains
        NDihedrals +=Sim.ChainStop[chain]-Sim.ChainStart[chain]-3
    end
    inv = 1.0/convert(R, (Sim.NSteps-Sim.EquilibrationTime+1)*NDihedrals*(2.0*π)/360.0)
    Sim.TorsionHist .*= inv
end


function decidePseudoHelicals(Sim::SimData{R,I}, Pseudohelical::Vector{Bool}, out::Vector{N}, step, lower=0.25, upper=1.3) where {R<:Real, I<:Integer, N<:Number}
    ### according to "Developing Bonded Potentials for a Coarse-Grained Model of Intrinsically Disordered Proteins" by Azamat Rizuan,..., Jeetain Mittal
    Pseudohelical .= false

    ### apply the (i+1, i+2) rule Dihedral -ID = atom-ID+1
    for chain in 1:Sim.NChains
        for atom in Sim.ChainStart[chain]:Sim.ChainStop[chain]-3
            if(Sim.TorsionAngles[atom, step]<upper && Sim.TorsionAngles[atom, step]>lower)
                Pseudohelical[atom+1]=true
                Pseudohelical[atom+2]=true
            end
        end
    end
    ### detect the helical states among the pseudo helical ones, n=1 therefore at least 3 consecutives pseudohelical states are needed and the first and last do not count as helical
    inhelix=false
    helixstart=0
    helixstop=0
    for atom in 1:Sim.NAtoms
        if inhelix && ~Pseudohelical[atom]
            helixstop=atom-1
            if helixstop-helixstart>2
                for id in helixstart+1:helixstop-1
                    out[id]+=1
                end
            end
            inhelix=false
        elseif ~inhelix && Pseudohelical[atom]
            helixstart=atom
            inhelix=true
        end
    end
end

function computeAlphaHelix(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    ### according to "Developing Bonded Potentials for a Coarse-Grained Model of Intrinsically Disordered Proteins" by Azamat Rizuan,..., Jeetain Mittal
    
    Sim.AlphaHelixProb = zeros(eltype(Sim.x), Sim.NAtoms)
    Pseudohelical = zeros(Bool,Sim.NAtoms)
    for step in 1:Sim.NSteps
        Pseudohelical .= false

        decidePseudoHelicals(Sim, Pseudohelical, Sim.AlphaHelixProb, step)
    end
    Sim.AlphaHelixProb ./=Float64(Sim.NSteps)
end

function computeBondLengthsENM(Sim::SimData{R,I}, ENM_Def_Path::String, ProtName::String) where {R<:Real,I<:Integer}
    ### uses externally defined ENM storage
    entry= load(ENM_Def_Path)[ProtName]
    N_Bonds, types, typeid, pairs,type_mapping = entry
    typed_list = [Vector{Tuple{Int32,Int32}}() for name in types]

    ### compute all pairs of interaction for one bond type
    if minimum(typeid)==0 typeid .+= 1 end
    range = Sim.ChainStart[1]:Sim.ChainStop[1]
    for (i,pair) in enumerate(pairs)
        if pair[1]∈range && pair[2]∈range
            push!(typed_list[typeid[i]], pair)
        end
    end

    step=1000
    for (c,inter) in enumerate(typed_list)
        inter_pairs = sort(inter)

        distances=[]
        for step in 1:10:Sim.NSteps
            for start in Sim.ChainStart
                for pair in inter
                    i,j= start .+ pair
                    dist = distance(Sim, [Sim.x[i,step], Sim.y[i,step], Sim.z[i,step]], [Sim.x[j,step], Sim.y[j,step], Sim.z[j,step]])
                    push!(distances, dist)
                end
            end
        end
        println("Type $c")
        println("Mean $(mean(distances)) vs $(type_mapping[types[c]][:r]*10.0)")
        println("Std  $(std(distances))")

    end
end

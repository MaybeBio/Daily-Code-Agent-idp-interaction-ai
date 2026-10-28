
#using LaTeXStrings

mutable struct HREMD_Data{FloatType<:AbstractFloat,IntType<:Integer}
    ChargeFile::String
    NReplica::IntType
    FloatType::FloatType
    ReplicaVariables::Dict{String,Any}
    ReplicaData::Vector{SimData{FloatType, IntType}}
    HREMD_Data() = new{Float32, Int32}("",1,one(Float32),  Dict{String, Any}(), [])
end

function readHREMDLmp!(HSim::HREMD_Data, LMPFile::String)
    file = open(LMPFile, "r")
    while ! eof(file)
        line = readline(file)  
        split = Base.split(line)
        if length(split)==0 continue end
        if( split[1]=="variable") 
            if split[3]=="world" ### variable assignment in lammps , mostly for REMD/Charge-HREMD 
                if occursin("." ,split[4]) ### do shady type assumption based on 
                    HSim.ReplicaVariables[split[2]] = vcat(parse.(eltype(HSim.FloatType),split[4:end]))
                else
                    HSim.ReplicaVariables[split[2]] = vcat(parse.(eltype(HSim.NReplica),split[4:end]))
                end
            end
        end
    end
    close(file)
    nothing
end

function initDataHREMD(Path::String, LmpName::String, ChargeFile::String; StepFrequency=1, Reparse=true, LoadAll=true, Reduce=1, EquilibrationTime=1)
    HSim = HREMD_Data()
    readHREMDLmp!(HSim, Path*LmpName)
    HSim.ChargeFile=ChargeFile
    ### check if HREMD is consistent
    if "id" in keys(HSim.ReplicaVariables)
        HSim.NReplica = length(HSim.ReplicaVariables["id"])
        for key in keys(HSim.ReplicaVariables)
            if length(HSim.ReplicaVariables[key]) != HSim.NReplica
                printstyled("HREMD information incomplete. Lammps world variable $(key) does not have the appropriate size.\n"; color=:red)
            end
        end
    else
        printstyled("HREMD information incomplete. No replica id with key name \"id\" has been specified.\n"; color=:red)
    end

    ### read replica information

    bla = Dict{String, Any}([(key, HSim.ReplicaVariables[key][1] ) for key in keys(HSim.ReplicaVariables)])
    HSim.ReplicaData = [initData(Path, LmpName; StepFrequency=StepFrequency, Reparse=Reparse, LoadAll=LoadAll, Reduce=Reduce, EquilibrationTime=EquilibrationTime, LammpsVariables=Dict{String, Any}([(key, HSim.ReplicaVariables[key][id] ) for key in keys(HSim.ReplicaVariables)]), BasePathAdd="Replica_$(id)/")  for id in HSim.ReplicaVariables["id"] ]
    return HSim
end

function BoltzmannWeight(E::T, Temp::T) where {T<:Real}
    return exp(-4.184*E/6.02214076/(Temp*1.38064)) ### expecting energies in kcal/mol, Temp in K

end

function ReweightingHREMD(ChargeFile::String, Sims::Vector{SimData{T, I}}, Start::I, StepWidth::I, ChargeFactor::T; Temp=280.0) where {I<:Integer, T<:Real}
    Weights =Vector{Vector{eltype(Sims[1].x)}}()
    Charges, ID= readChargeAssignments(ChargeFile)
    NSims = size(Charges)[1]

    minLength=10000000000
    for (SimID,Sim) in enumerate(Sims) 
        StepWidth= Int(1/Sim.Reduce)
        #StepWidth=1
        End= length(Charges[SimID, :])*StepWidth
        #StepWidth=1
        minLength= min(minLength, length(Sim.Energies[Start+StepWidth:StepWidth:End,1]))
    end

    Energies = zeros(T, (NSims, minLength))
    CorrectedEnergies = zeros(T, (NSims, minLength))

    for (SimID,Sim) in enumerate(Sims) 
        EtotalID = Sim.EnergyDict["etotal"]
        ECoulID = Sim.EnergyDict["ecoul"]
        StepWidth= Int(1/Sim.Reduce)
        #StepWidth=1
        End= length(Charges[SimID, :])*StepWidth
       # StepWidth=1
        Weight  =  BoltzmannWeight.(  Sim.Energies[Start+StepWidth:StepWidth:End,EtotalID] .+ Sim.Energies[Start+StepWidth:StepWidth:End,ECoulID].* (Charges[SimID,1:end-1].- ChargeFactor)./ChargeFactor, Temp)
        push!(Weights, Weight)

        Energies[SimID, :] .= Sim.Energies[Start+StepWidth:StepWidth:End,EtotalID][1:minLength]
        CorrectedEnergies[SimID, :] .= (Sim.Energies[Start+StepWidth:StepWidth:End,EtotalID] .+ Sim.Energies[Start+StepWidth:StepWidth:End,ECoulID].* (Charges[SimID,1:end-1].- ChargeFactor)./ChargeFactor)[1:minLength]
    end
    #ComputeWHAM(NSims, minLength, ID, Energies, CorrectedEnergies, convert(T, 1.0/(Temp*1.38064)))
    return Weights, ID
    LmpFiles = [[Sim.SimulationName*".lmp" for Sim in Sims]]
    Paths = [[Sim.BasePath for Sim in Sims]]
    BasePathAdd = [[Sim.DataPath[length(Sim.BasePath*"/Data/")+1:end] for Sim in Sims]]
    Fields=["RGSeries"]
    return AverageFields(Paths,LmpFiles, Fields, Start; StepWidth=StepWidth, Weights=Weights,BasePathAdd=BasePathAdd)
end

function ComputeWHAM(NSims::I, NSteps::I, IDs::Array{Int32}, Energies::Array{T}, CorrectedEnergies::Array{T},β::T) where {I<:Integer, T<:Real}
    ### expect λ*E_Coul as energies
    f     = zeros(BigFloat, NSims)
    f_old = ones(BigFloat, NSims)*2
    P = zeros(BigFloat, ( NSims, NSteps,))
    println("NSims $(NSims), NSteps, $(NSteps)")
    i =0
    while(sum(abs.(f.-f_old))>1.0 && i<500)
        for t in 1:NSteps
            weight = NSteps*sum(exp.(f[IDs[:,t]].-β.*Energies[:,t]))
            for n in 1:NSims
                P[n,t] = exp(f[IDs[n,t]]-β*CorrectedEnergies[n,t])/weight
            end
        end
        f_old .= f
        f .= -1.0.*log.(sum(P, dims=2))
        i += 1
        if i%100==1
            println("Iteration: $(i)  $(sum(abs.(f.-f_old)))")
        end
    end
end

function readChargeAssignments(fileName::String)
    file = open(fileName,"r+")
    NLines = countlines(file)-2
    close(file)

    file = open(fileName,"r+")
    readline(file)
    line = readline(file)
    Qvals = parse.(Float64, split(line)[2:end])
    NReplica = length(Qvals)

    Charges = zeros(Float64, (NReplica, NLines))
    ID  = zeros(Int32, (NReplica, NLines))
    Step = zeros(Int32, NLines)

    Restarts= zeros(Int32, 0)
    cnt = 0
    while ~eof(file) 
        line = readline(file)
        if line[1]=='#'
            push!(Restarts,parse(Int32,split(line)[3]))
            continue
        end
        line_val=parse.(Int32,split(line)[1:end])
        cnt+=1
        Step[cnt] = line_val[1]

        for i in 1:NReplica
            Charges[i, cnt] = Qvals[line_val[i+1]+1]
            ID[i, cnt ] = line_val[i+1]+1
        end
    end

    NRes=length(Restarts)
    return Step[1:end-NRes], Charges[:,1:end-NRes], ID[:,1:end-NRes], Restarts
end

function PlotEnergyRelaxation(ChargeFile::String, Sims::Vector{SimData{T, I}}, Start::I, StepWidth::I, ChargeFactor::T; Temp=280.0) where {I<:Integer, T<:Real}

    Charges, ID= readChargeAssignments(ChargeFile)

    IsSame = Charges[:,2:end].==Charges[:,1:end-1]

    AxisRanges=Vector{Vector{Tuple{I,I}}}()
    Qvals=Vector{Vector{Float64}}()

    for IRep in 1:6 #length(Sims)
        push!(Qvals, Vector{Float64}())
        push!(AxisRanges, Vector{Tuple{Int32,Int32}}())
        i = 1
        while(i<size(Charges)[2])
            old = i 
            push!(Qvals[IRep], Charges[IRep, i])
            while(i <size(Charges)[2] && IsSame[IRep,i])
                i+=1
            end
            push!(AxisRanges[IRep], (old, i))
            i+=1
        end
    end
end

function GetHistograms(Data::HREMD_Data{T,I}) where {I<:Integer, T<:Real}
    SwapSteps, Charges, ID, Restarts = readChargeAssignments(Data.ChargeFile)
    return MatchSwapAndOutputSteps(Data, SwapSteps)
end

function MatchSwapAndOutputSteps(Data::HREMD_Data{T, I}, SwapSteps::Array{I}) where {I<:Integer, T<:Real}
    Sim = Data.ReplicaData[1] # all replicas have the same data that is needed here
    Energies = Sim.Energies

    TrajID = zeros(I, Sim.NSteps)

    ### Match the Trajectoriesf
    for i in 1:Sim.NSteps
        tmp = i * Sim.TrajWriteOutFreq
        for (s,e) in enumerate(SwapSteps)
            if e > tmp
                TrajID[i] = s-1
            end
        end
    end
    return TrajID
end

function PlotSwapSteps(ChargeFile::String)
    return plotSwapSteps(ChargeFile)
end
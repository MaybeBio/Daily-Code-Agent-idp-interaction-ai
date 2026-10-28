#include("./Analysis.jl")

using Statistics

### Struct to efficiently implement equal weights
struct FakeArray <: AbstractArray{Any,1}
    val::Number
    size::Int64
end

@inline function Base.getindex(a::FakeArray, ind::Integer)
    return a.val
end

@inline function Base.size(a::FakeArray)
    return (1,a.size)
end

function AverageFields(Paths::Vector{Vector{String}}, LmpFiles::Vector{Vector{String}}, Fields::Vector{String}, Start::Integer; StepWidth::Integer=typemax(Int32), Weights=[], BasePathAdd="", HOOMD=false)

    if length(Weights)==0
        Weights = fill(FakeArray(1,1), length(Paths))
    end
    if length(BasePathAdd)==0
        BasePathAdd = [fill("", length(file)) for file in LmpFiles]
    end

    Results = Dict{String,  Tuple{Any, Vector{Any}}}()
    KeyNames = Vector{Tuple{eltype(Fields), String}}()
    Tmp = Dict{Tuple{Int, String}, Any}()#Dict{Tuple{Tuple{Int32,String}, Any}}#Vector{Array{eltype(Any)}}}}
    Axes = Dict("RGSeries" =>2 , "AlphaHelixProb"=>0, "VarLocalPersistenceLength"=>0, "Asphericity"=>2,"Acylindricity"=>2, "ShapeAsymmetry"=> 2,"REEHists"=>0 )

    ### Initialise Data Structure
    for (ProtID, PathVec) in enumerate(Paths)
        Data = initData(PathVec[1]; LmpName=LmpFiles[ProtID][1], Reparse=false, LoadAll=false,BasePathAdd=BasePathAdd[ProtID][1], HOOMD=HOOMD)
        (RGMeasureStep, End) = GetFieldsFromData(Data, ["RGMeasureStep", "NSteps"])
        LocStep = RGMeasureStep < StepWidth ? RGMeasureStep : StepWidth
        for (ID, field) in enumerate(Fields)
            Key = "Mean_"*String(field)
            push!(KeyNames, (field, Key))
            (field_data,) = GetFieldsFromData(Data, [String(field)])
            NDims = ndims(field_data)
            if Axes[field]>NDims || Axes[field]<0
                printstyled("Axes indices for field: \"$field \" is too low/high for ndims=$(NDims). Stop working!\n"; color=:red)
                return Results
            end
            
            #start_val = Start > size(Data)[Axes[field]] ?  1 : Start
            #end_val = Start > size(Data)[Axes[field]] ?  size(Data)[Axes[field]] : End

            if Axes[field]==0 ### already averages, take number of measurements as weight
                Tmp[(ProtID,Key)]=CMA_Init(field_data,weight=length(Start:LocStep:End))
            else
                Tmp[(ProtID,Key)]=CMA_Init(selectdim(field_data, Axes[field], Start), weight=Weights[ProtID][1])
               # local tmp_axis = axes()
                for Step in Start+LocStep:End
                    CMA_Update(selectdim(field_data, Axes[field], Step), Tmp[(ProtID,Key)], weight=Weights[ProtID][Step])
                end
            end
        end
    end

    for (ProtID, PathVec) in enumerate(Paths)
         ### Add other simulation runs
        for (SimRun,Path) in enumerate(PathVec[2:end])
            Data = initData(PathVec[1]; LmpName=LmpFiles[ProtID][1], Reparse=false, LoadAll=false,BasePathAdd=BasePathAdd[ProtID][1], HOOMD=HOOMD)
            (RGMeasureStep, End) = GetFieldsFromData(Data, ["RGMeasureStep", "NSteps"])
            LocStep = RGMeasureStep < StepWidth ? RGMeasureStep : StepWidth
            for (ID, field) in enumerate(Fields)
                (field,Key)=KeyNames[ID]
                ### Dont write into field! They are mutable
                (field_data,) = GetFieldsFromData(Data, [String(field)])
                if Axes[field]==0
                    CMA_Update(field_data, Tmp[(ProtID,Key)]; weight=length(Start:LocStep:End))
                else
                    for Step in Start+LocStep:End
                        CMA_Update(selectdim(field_data,Axes[field], Step), Tmp[(ProtID,Key)], weight=Weights[ProtID][Step])
                    end
                end
            end
        end

        ### compute final result
        for (ID, field) in enumerate(Fields)
            (field,Key)=KeyNames[ID]
            (mean,std, error) = CMA_Finalize(Tmp[(ProtID,Key)])
            if ProtID==1
                Results[Key] = (Vector{typeof(mean)}, Vector{eltype(mean)}(undef, length(Paths)))
                Results["$(Key)_σ"] =  (Vector{typeof(std)},  Vector{eltype(std)}(undef, length(Paths)))
                Results["$(Key)_Error"] = (Vector{typeof(error)},  Vector{eltype(error)}(undef, length(Paths)))
            end
            Results[Key][2][ProtID] = mean
            Results["$(Key)_σ"][2][ProtID] =std
            Results["$(Key)_Error"][2][ProtID] = error
        end
    end

    return Results
end

function Average(Data::Vector, Axis::Int32)   
    NDataPoints = Float64(size(Data)[Axis])
    Mean = Statistics.mean(Data, dims=Axis) *NDataPoints
    Error = Statistics.stdm(Data,mean, dims=Axis)/sqrt(NDataPoints)*NDataPoints
    return Mean, Error, NDataPoints
end


### CMA= Cumulated moving average, and error
mutable struct CMA_Obj
    N::Number
    SUM::Union{Number,Array{<:Number},Array{Array{<:Number}}, Array{Array{Array{<:Number}}}}#::Array{<:Number}
    SQR::Union{Number,Array{<:Number},Array{Array{<:Number}}, Array{Array{Array{<:Number}}}}#Array{<:Number}
    Weight::Union{Number,Array{<:Number}}
end

function CMA_Init(data; weight=BigFloat(1.0))
    if data isa Array{<:Number} ### casting doesnt work with vectors of vectors
        data = BigFloat.(data) 
    end
    if weight isa Array{<:Number}
        weight= BigFloat.(weight)
        weighted_data = data.*weight
    else
        weighted_data = data*weight
    end

    Out = CMA_Obj(BigFloat(1.0), weighted_data, data.*weighted_data, weight)
    return Out
end

function CMA_Update(data, Obj; weight=BigFloat(1.0))
    if data isa Array{<:Number}
        data = BigFloat.(data) 
    end   

    if weight isa Array{<:Number}
        weight= BigFloat.(weight)
        weighted_data = data.*weight
    else
        weighted_data = data*weight
    end
    Obj.N += 1.0
    Obj.SUM .+= weighted_data
    Obj.SQR .+= data.*weighted_data
    Obj.Weight += weight
end

function CMA_Finalize(Obj)
    CMA = Obj.SUM ./Obj.Weight
    if Obj.N==1
        CMA_σ = @. sqrt( ((Obj.SQR-2.0*CMA*Obj.SUM +CMA^2*Obj.Weight)/(Obj.Weight)))
        return Float64.(CMA), Float64.(CMA_σ), zeros(Float64, size(CMA))
    end

    if all(isapprox.(Obj.Weight,1 , atol=10^-5))
        CMA_σ = @. sqrt( ((Obj.SQR-2.0*CMA*Obj.SUM +CMA^2*Obj.Weight)/(Obj.Weight)))
        CMA_err= @. CMA_σ/sqrt((Obj.N-1))
    else
        CMA_σ = @.  sqrt( ((Obj.SQR-2.0*CMA*Obj.SUM +CMA^2*Obj.Weight)/(Obj.Weight-1)))
        CMA_err= @. CMA_σ/sqrt( (Obj.N))
    end

    return Float64.(CMA), Float64.(CMA_σ), Float64.(CMA_err)
end

using StaticArrays


function computeRGComponentSeries!(Sim::SimData{R,I}) where{R<:Real, I<:Integer}
    Sim.RGComponentSeries=zeros(R, Sim.NChains, Sim.NSteps,3)
    Sim.RGSeries=zeros(R, Sim.NChains, Sim.NSteps)

    x=zeros(R, Sim.NChains)
    y=zeros(R, Sim.NChains)
    z=zeros(R, Sim.NChains)
    for (cnt,step) in enumerate(1:Sim.NSteps)
        for (chain, start, stop) in zip(collect(1:Sim.NChains), Sim.ChainStart, Sim.ChainStop)
            x[chain] =sum((Sim.x_uw[start:stop, step] .- Sim.COM_uw[chain,1,step]).^2 .*Sim.Masses[start:stop])
            y[chain] =sum((Sim.y_uw[start:stop, step] .- Sim.COM_uw[chain,2,step]).^2 .*Sim.Masses[start:stop])
            z[chain] =sum((Sim.z_uw[start:stop, step] .- Sim.COM_uw[chain,3,step]).^2 .*Sim.Masses[start:stop])
        end
        x ./= Sim.ChainMasses
        y ./= Sim.ChainMasses
        z ./= Sim.ChainMasses

        Sim.RGComponentSeries[:, cnt,1] .= @. sqrt(x)
        Sim.RGComponentSeries[:, cnt,2] .= @. sqrt(y)
        Sim.RGComponentSeries[:, cnt,3] .= @. sqrt(z)
        Sim.RGSeries[:,cnt] .= @. sqrt((x+y+z))
    end
end

function computeRGSeries!(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    computeRGComponentSeries!(Sim)

    #Sim.RGSeries=zeros(eltype(Sim.x), Sim.NChains, Sim.NSteps-Sim.EquilibrationTime+1)

    #Sim.RGSeries[:,:] .= sqrt.(Sim.RGComponentSeries[:,:,1].^2 .+Sim.RGComponentSeries[:,:,2].^2 .+Sim.RGComponentSeries[:,:,3].^2)
 
end

function computeREESeries!(Sim::SimData{R,I}, StartStops=[[]]) where{R<:Real, I<:Integer}
    #Array{Array{Tuple{I,I}}}()
    diff=0.
    if length(StartStops[1])==0
        StartStops= [Array{Tuple{I,I}}(undef, 0,) for chain in 1:Sim.NChains]
    end

    for chain in 1:Sim.NChains
        push!( StartStops[chain] , (Sim.ChainStart[chain], Sim.ChainStop[chain]))
    end
    Sim.REESeries = zeros(R,  Sim.NSteps-Sim.EquilibrationTime+1,Sim.NChains,maximum(length.(StartStops)) )
    Sim.REEComponents = zeros(R,  Sim.NSteps-Sim.EquilibrationTime+1,Sim.NChains,maximum(length.(StartStops)),3) 

    for (cnt,step) in enumerate(Sim.EquilibrationTime:Sim.NSteps)
        for chain in 1:Sim.NChains
            for (id, (start, stop)) in enumerate(StartStops[chain])
                Sim.REEComponents[cnt,chain,id,1] =(Sim.x_uw[ start, step]- Sim.x_uw[stop, step])
                Sim.REEComponents[cnt,chain,id,2] =(Sim.y_uw[ start, step]- Sim.y_uw[stop, step])
                Sim.REEComponents[cnt,chain,id,3] =(Sim.z_uw[ start, step]- Sim.z_uw[stop, step])
                Sim.REESeries[cnt,chain,id] = sqrt.(Sim.REEComponents[cnt,chain,id, 1]^2 .+ Sim.REEComponents[cnt,chain,id,2]^2 .+Sim.REEComponents[cnt,chain,id,3]^2)
            end
        end
    end
end

function computeREEAutocorr(Sim::SimData{R,I},  NLags=5000)  where{R<:Real, I<:Integer}
    DecorrelationTimes = zeros(eltype(Sim.x),(length(axes(Sim.REESeries,2)), length(axes(Sim.REESeries,3))))
    NLags_, AutoCorr, LagTimes = InitAutoCorrelationArrays(Sim; NLags=NLags)
    _, AutoCorr2, _ = InitAutoCorrelationArrays(Sim; NLags=NLags)


    for chain in 1:Sim.NChains
        for Section in axes(Sim.REESeries,3)
            DecorrelationTimes[chain, Section] = computeAutoCorrelation(Sim, Sim.REESeries[:,chain,Section], AutoCorr,LagTimes, NLags_)
        end
    end

    NMax= min((NLags+Sim.EquilibrationTime), Sim.NSteps)

    mean_=zeros(R, length(Sim.NChains), length(axes(Sim.REESeries,3)),3)
    mean2_=zeros(R, length(Sim.NChains), length(axes(Sim.REESeries,3)))
    Δr = zeros(R, NLags)

    x = deepcopy(Sim.REEComponents[:,:,:,1])
    y = deepcopy(Sim.REEComponents[:,:,:,2])
    z = deepcopy(Sim.REEComponents[:,:,:,3])

    for chain in 1:Sim.NChains
        for Section in axes(Sim.REESeries,3)
            x[:,chain,Section] .-= StatsBase.mean(x[:,chain,Section], dims=1)
            y[:,chain,Section] .-= StatsBase.mean(y[:,chain,Section], dims=1)
            z[:,chain,Section] .-= StatsBase.mean(z[:,chain,Section], dims=1)
        end
    end

    AutoCorr2=zeros(R, NLags_, Sim.NChains,length(axes(Sim.REESeries,3) ))
    for chain in 1:Sim.NChains
        for Section in axes(Sim.REESeries,3)
            for (i,step) in enumerate(Sim.EquilibrationTime:Sim.NSteps)
                start = i 
                stop = min(Sim.NSteps-Sim.EquilibrationTime,i+NLags-1)
                N=stop-start+1

                Δr[1:N]  .= x[i,chain,Section].* x[start:stop,chain,Section]
                Δr[1:N] .+= y[i,chain,Section].* y[start:stop,chain,Section]
                Δr[1:N] .+= z[i,chain,Section].* z[start:stop,chain,Section]
                AutoCorr2[1:N, chain, Section] .+= Δr[1:N]
                #=for j in i:min(Sim.NSteps-Sim.EquilibrationTime,i+NLags-1)
                    AutoCorr2[j-i+1, chain, Section] += x[i,chain,Section]*x[j,chain,Section];
                    AutoCorr2[j-i+1, chain, Section] += y[i,chain,Section]*y[j,chain,Section]
                    AutoCorr2[j-i+1, chain, Section] += z[i,chain,Section]*z[j,chain,Section]
                end=#
                #start = i
                #stop = min(Sim.NSteps-Sim.EquilibrationTime,i+NLags-1)
                #AutoCorr2[start:stop, chain, Section] .+=  Sim.REEComponents[i,chain,Section].*Sim.REEComponents[start:stop,chain,Section]
                #mean[chain, Section,1] += Sim.REEComponents[i,chain,Section,1]
                #mean[chain, Section,2] += Sim.REEComponents[i,chain,Section,2]
                #mean[chain, Section,3] += Sim.REEComponents[i,chain,Section,3]

                #+(Sim.REEComponents[i,chain,Section],mean[chain, Section],mean[chain, Section]);
                
               
            end
           # mean2[chain, Section] = sum(  Sim.REEComponents[:,chain,Section,1]*Sim.REEComponents[i,chain,Section];
        end
    end
    
    #mean ./= length(Sim.EquilibrationTime:Sim.NSteps)
    #mean2 ./= length(Sim.EquilibrationTime:Sim.NSteps)

    #=
    autozero = zeros(R, length(Sim.NChains), length(axes(Sim.REESeries,3)))
    for chain in 1:Sim.NChains
        for Section in axes(Sim.REESeries,3)
            for (i,step) in enumerate(Sim.EquilibrationTime:Sim.NSteps)
                autozero[chain, Section] +=  (Sim.REEComponents[i,chain,Section] - mean[chain,Section])*(Sim.REEComponents[i,chain,Section] - mean[chain,Section]);
            end
        end
    end
    autozero ./= length(Sim.EquilibrationTime:Sim.NSteps)
    =#
    #AutoCorr2[:, :, :] .-=mean2
    for i in 1:NLags_
        AutoCorr2[i, :, :] ./= (Sim.NSteps-Sim.EquilibrationTime-i)
    end
   # AutoCorr2[:, :, :] .-=mean2


    for chain in 1:Sim.NChains
        for Section in axes(Sim.REESeries,3)
            #AutoCorr2[:, chain, Section] .-= mean2[chain, Section]
            AutoCorr2[:,chain,Section] ./= AutoCorr2[1,chain,Section]
        end
    end

    Sim.REEVecSeries = AutoCorr2

    return DecorrelationTimes
end

function computeREEHist(Sim::SimData{R,I}, StartStops=[[]]; ScalFac=5.0, DecorrelationTimes=[]) where {R<:Real, I<:Integer}

    if StartStops!=[[]]
        for chain in 1:Sim.NChains
            push!( StartStops[chain]  ,convert.(eltype(first(StartStops[1])),(Sim.ChainStart[chain], Sim.ChainStop[chain])))
        end
    else
        StartStops = [ [(Sim.ChainStart[chain], Sim.ChainStop[chain])] for chain in 1:Sim.NChains ]
    end

    ### check format of Decorrelation times
    if size(DecorrelationTimes)!=( Sim.NChains,maximum(length.(StartStops)))
        printstyled("Wrong shape of decorrelation times for REE-Hist. Continue with RG_MeasureStep."; color=:yellow)
        DecorrelationTimes= fill(Sim.RGMeasureStep, (Sim.NChains,maximum(length.(StartStops))))
    end


    ### compute axes of REE Hist
    maxLen=0
    Sim.REEHistsAxes = Array{UnitRange{Int32}}(undef,Sim.NChains,maximum(length.(StartStops)))
    for chain in Sim.NChains
        for (id, (start,stop)) in enumerate(StartStops[chain])
            Sim.REEHistsAxes[chain,id] =1:(convert(eltype(Sim.NSteps),ceil((stop-start)*ScalFac))+20)
            maxLen = max(maxLen,convert(eltype(Sim.NSteps),ceil((stop-start)*ScalFac))+20)
        end
    end

    ### Compute Histogram
    Sim.REEHists =  zeros(Sim.NChains,maximum(length.(StartStops)), maxLen) 
    if length(Sim.FrameWeights)>=Sim.NSteps ### If weights exist for different time steps
        for chain in 1:Sim.NChains
            for (id, (start,stop)) in enumerate(StartStops[chain]) 
                for (cnt, step) in enumerate(Sim.EquilibrationTime:DecorrelationTimes[chain, id]:Sim.NSteps)

                if Sim.FrameWeights[step]>0 ### non physical positions have a weight of 0 in Rosenbluth method
                Sim.REEHists[chain,id, Int(round(Sim.REESeries[cnt,chain,id]))+1] += Sim.FrameWeights[step]
                end
            end
        end
    end
        Sim.REEHists /= sum(Sim.FrameWeights)
    else ### unweighted frames
        for chain in 1:Sim.NChains
            for (id,_) in enumerate(StartStops[chain]) 
                for (cnt, step) in enumerate(Sim.EquilibrationTime:DecorrelationTimes[chain, id]:Sim.NSteps)

                    Sim.REEHists[chain,id, Int(round(Sim.REESeries[cnt,chain,id]))+1] += 1
                end
            end
        end
        Sim.REEHists /= length(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)
    end
end 

function computeRGCorrelationTime(Sim::SimData{R,I}, NLags=5000) where{R<:Real, I<:Integer}
    NLags = min(NLags, length(Sim.EquilibrationTime:Sim.NSteps))
    Sim.RGAutocorr = zeros(Sim.NChains, NLags)
    lags = Vector(range(0,NLags-1))


    for chain in 1:Sim.NChains
        Sim.RGAutocorr[chain,:] .= StatsBase.autocor(Sim.RGSeries[chain,Sim.EquilibrationTime:Sim.NSteps], lags;demean=true)
    end

    half_time =zeros(I, Sim.NChains)
    for chain in Sim.NChains
        for i in 1:NLags
            if(Sim.RGAutocorr[chain, i])<0.5
                half_time[chain]=i
                break
            end
        end
    end
    Sim.RGMeasureStep = max(1,ceil(I,Statistics.mean(half_time)))
end

function computeInertiaTensor(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Sim.InertiaTensorEigVals = zeros(eltype(Sim.x), (3,Sim.NChains, Sim.NSteps))

    xtmp = zeros(R, Sim.NAtoms)
    ytmp = zeros(R, Sim.NAtoms)
    ztmp = zeros(R, Sim.NAtoms)

    Matrix = zeros(eltype(Sim.x), 3,3)
    for step in 1:Sim.NSteps
        for chain in 1:Sim.NChains
            Matrix .= 0
            for atom in Sim.ChainStart[chain]:Sim.ChainStop[chain]
                @inbounds xtmp[atom] = (Sim.x_uw[atom, step]-Sim.COM_uw[chain, 1, step])
                @inbounds ytmp[atom] = (Sim.y_uw[atom, step]-Sim.COM_uw[chain, 2, step])
                @inbounds ztmp[atom] = (Sim.z_uw[atom, step]-Sim.COM_uw[chain, 3, step])
                @inbounds Matrix[1,1]+= xtmp[atom]*xtmp[atom]*Sim.Masses[atom]
                @inbounds Matrix[1,2]+= xtmp[atom]*ytmp[atom]*Sim.Masses[atom]
                @inbounds Matrix[1,3]+= xtmp[atom]*ztmp[atom]*Sim.Masses[atom]
                @inbounds Matrix[2,2]+= ytmp[atom]*ytmp[atom]*Sim.Masses[atom]
                @inbounds Matrix[2,3]+= ytmp[atom]*ztmp[atom]*Sim.Masses[atom]
                @inbounds Matrix[3,3]+= ztmp[atom]*ztmp[atom]*Sim.Masses[atom]
            end
            @inbounds Matrix[2,1]=Matrix[1,2]
            @inbounds Matrix[3,1]=Matrix[1,3]
            @inbounds Matrix[3,2]=Matrix[2,3]
            @inbounds Matrix ./= Sim.ChainMasses[chain] 
            @inbounds Sim.InertiaTensorEigVals[:, chain,step] .= eigvals(Matrix)
        end
    end
end

function computePolymerCharacteristics(Sim::SimData{R,I}, Start=100)  where {R<:Real, I<:Integer}
    computePolymerCharacteristics(Sim, Start, [Start:Sim.RGMeasureStep:Sim.NSteps for _ in 1:Sim.NChains])
end

function computePolymerCharacteristics(Sim::SimData{R,I}, Start, Range::Vector{Vector{I2}}) where {R<:Real, I<:Integer, I2<:Integer}
    λ = Sim.InertiaTensorEigVals
    Sim.ShapeAsymmetry =  1.0 .-3.0.*(λ[1,:,:].*λ[2,:,:].+λ[1,:,:].*λ[3,:,:].+λ[2,:,:].*λ[3,:,:])./(λ[1,:,:].+λ[2,:,:].+λ[3,:,:]).^2 ### arash + janka
    Sim.ParallelInertiaTensor = @. (λ[1,:,:]+λ[2,:,:])/2.0
    Sim.AspectRatio = @. Sim.ParallelInertiaTensor[:,:]/λ[3,:,:] ### Arash tanja paper
    Sim.Asphericity = @. (λ[3,:,:]-0.5*(λ[1,:,:]+λ[2,:,:]))/(λ[1,:,:]+λ[2,:,:]+λ[3,:,:]) ### wikipedia + own norm
    Sim.Acylindricity = @. (λ[2,:,:] - λ[1,:,:])/(λ[1,:,:].+λ[2,:,:].+λ[3,:,:]) ### https://journals.aps.org/pre/pdf/10.1103/PhysRevE.106.064606

    #Range = isnothing(Range) ? [Start:Sim.RGMeasureStep:Sim.NSteps for _ in 1:Sim.NChains] :  Range
    NDataPoints = convert(R, (Sim.NSteps-Start+1)÷Sim.RGMeasureStep)
    Sim.MeanShapeAsymmetry = [sum(Sim.ShapeAsymmetry[i,Range[i]])./NDataPoints for i in 1:Sim.NChains]
    Sim.MeanAspectRatio    = [sum(Sim.AspectRatio[   i,Range[i]])./NDataPoints for i in 1:Sim.NChains]
    Sim.MeanAsphericity    = [sum(Sim.Asphericity[   i,Range[i]])./NDataPoints for i in 1:Sim.NChains] 
    Sim.MeanAcylindricity  = [sum(Sim.Acylindricity[ i,Range[i]])./NDataPoints for i in 1:Sim.NChains]
end

function computeRG(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Sim.RGMeanByChain = zeros(R, Sim.NChains)
    Sim.RGErrByChain  = zeros(R, Sim.NChains)
    cnt = 0
    Sim.RGMean = zero(R)
    Sim.RGErr  = zero(R)

    for chain in 1:Sim.NChains
        cnt =0
        for step in Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps #range(1,Sim.NSteps-Sim.EquilibrationTime-1,step=Sim.RGMeasureStep)
            cnt+=1
            Sim.RGMeanByChain[chain] +=   Sim.RGSeries[chain, step]  
        end
        Sim.RGMeanByChain[chain]/=cnt
        Sim.RGMean += Sim.RGMeanByChain[chain]
        for step in Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps #range(1,Sim.NSteps-Sim.EquilibrationTime-1,step=Sim.RGMeasureStep)
            Sim.RGErrByChain[chain] += (Sim.RGSeries[chain, step] -Sim.RGMeanByChain[chain])^2
        end
        Sim.RGErrByChain[chain] = sqrt( Sim.RGErrByChain[chain]/cnt)/sqrt(cnt-1)
        Sim.RGErr += Sim.RGErrByChain[chain]
    end
    Sim.RGMean/=Sim.NChains
    Sim.RGErr /=Sim.NChains
end


@doc raw"""
    computeIntraChainScalingNaiv(Sim::SimData{R,I})

Computes the square intra chain distance |r_i-r_j| for all combinations and averages over the trajectory.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Create**:
* Sim.IntraChainScalingNaiv
"""
function computeIntraChainScalingNaiv(Sim::SimData{R,I}) where {R<:Real, I<:Integer}

    Results = [zeros(R, Sim.ChainLength[I],Sim.ChainLength[I] ) for I in 1:Sim.NChains]
    for step in Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
        for C in 1:Sim.NChains
            for (i_rel,i) in enumerate(Sim.ChainStart[C]:Sim.ChainStop[C])
                for j in i+1:Sim.ChainStop[C]
                    j_rel = j - Sim.ChainStart[C]+1
                    dist_sqr = (Sim.x_uw[i,step]-Sim.x_uw[j,step])^2+(Sim.y_uw[i,step]-Sim.y_uw[j,step])^2+(Sim.z_uw[i,step]-Sim.z_uw[j,step])^2
                    Results[C][i_rel,j_rel] += dist_sqr
                end
            end
        end
    end
    
    N =  length(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)
    for C in 1:Sim.NChains
        for (i_rel,i) in enumerate(Sim.ChainStart[C]:Sim.ChainStop[C])
            for j in i+1:Sim.ChainStop[C]
                j_rel = j - Sim.ChainStart[C]+1
                Results[C][i_rel,j_rel] /= N
            end
        end
    end
    
    Sim.IntraChainScalingNaiv=Results
end


@doc raw"""
    computeIntraChainScalingSlidingWindow(Sim::SimData{R,I})

Computes the square intra chain distance |r_i-r_j|^2 and averages in a sliding window approach.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Create**:
* Sim.IntraChainScalingSlidingWindow
"""
function computeIntraChainScalingSlidingWindow(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Results = [zeros(R, Sim.ChainLength[C]-1 ) for C in 1:Sim.NChains]
    for step in Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
        for C in 1:Sim.NChains
            for (i_rel,i) in enumerate(Sim.ChainStart[C]:Sim.ChainStop[C])
                for j in i+1:Sim.ChainStop[C]
                    τ = j-i
                    dist_sqr = (Sim.x_uw[i,step]-Sim.x_uw[j,step])^2+(Sim.y_uw[i,step]-Sim.y_uw[j,step])^2+(Sim.z_uw[i,step]-Sim.z_uw[j,step])^2
                    Results[C][τ] += dist_sqr
                end
            end
        end
    end

    N =  length(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)
    for C in 1:Sim.NChains
        for (i, _) in enumerate( Results[C])
             Results[C][i] /= (length(Results[C]) -i+1)*N
        end
    end
    Sim.IntraChainScalingSlidingWindow=Results
end


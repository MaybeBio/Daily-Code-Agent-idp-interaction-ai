#include("./IO.jl")
using LoopVectorization, StatsBase, Statistics, OffsetArrays, LinearAlgebra, LsqFit, Graphs

include("./Analysis/CellList.jl")
include("./Analysis/ChargeBonds.jl")
include("./Analysis/Cluster.jl")
include("./Analysis/Polymer.jl")
include("./Analysis/Slab.jl")
include("./Analysis/BondedAnalysis.jl")
include("./Analysis/SequenceAnalysis.jl")
include("./Analysis/ContactMatrix.jl")
include("./Analysis/HREMD/Unify.jl")
include("./Analysis/HREMD/HREMD.jl")


function symmetriesMatrix(A::Matrix{T}) where {T<:Number}
    B = zeros(T, axes(A))
    for i in axes(A,1)
        for j in axes(A,2)
            if i==j continue end
            B[i,j] += A[i,j]
            B[j,i] += A[i,j]
        end
    end
    return B
end

function isNotInSameChain(Sim::SimData, i::I,j::I) where{I<:Integer}
    if abs(j-i)>Sim.MaxChainLength
        return true
    else
        for chain in 1:Sim.NChains
            if( i<= Sim.ChainStop[chain] &&  i>= Sim.ChainStart[chain])
                if (j >= Sim.ChainStart[chain]  && j<= Sim.ChainStop[chain])
                    return false
                else
                    return true
                end
            end
        end
    end
    #return true
end

function FitSeries(Sim::SimData, Series::Vector{Float64},Func, Init::Vector{Float64}, SaveAs; xlabel="", ylabel="", label="")
    axis= collect(Float64, axes(Series)[1])
    Fit=LsqFit.curve_fit(Func , axis, Series, Init)
    if Fit.converged
        plotFitSeries(axis, Series, Func, Fit, Sim.PlotPath*Sim.SimulationName*"_$SaveAs"; xlabel=xlabel, ylabel=ylabel, label=label)
        return Fit.param
    else
        printstyled("Fit did not converge!"; color=:red)
        return nothing
    end
end

function InitAutoCorrelationArrays(Sim::SimData; NLags=5000)
    NLags = min(NLags, length(Sim.EquilibrationTime:Sim.NSteps)-1)
    AutoCorr = zeros(eltype(Sim.x), NLags)
    lagTimes = Vector(range(1,NLags))
    return NLags, AutoCorr, lagTimes
end

function computeAutoCorrelation(Sim::SimData, Series, Autocorr,LagTimes, NLags)
    Autocorr[:] .= StatsBase.autocor(Series, LagTimes; demean=true)
    for i in 1:NLags
        if(Autocorr[i])<0.5*Autocorr[1]
            return i
        end
    end
    return NLags
end

#=function computeREEHist(Sim, StartStops=[[]]; ScalFac=5.0)
    for chain in 1:Sim.NChains
        insert!( StartStops[chain]  ,1, convert.(eltype(StartStops[1][1]),(Sim.ChainStart[chain], Sim.ChainStop[chain])))
    end
    maxDiff=maximum()
    Sim.REEHistsAxes = [[1:(convert(eltype(Sim.NSteps),ceil((y[2]-y[1])*ScalFac))+20)  for y in x]  for x  in StartStops]
    Sim.REEHists =  [[  zeros(Int(ceil((y[2]-y[1])*ScalFac))+20)  for y in x]  for x in StartStops]

    for (cnt, step) in enumerate(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)
        for chain in 1:Sim.NChains
            for id in 1:length(StartStops[chain])
                Sim.REEHists[chain][id][Int(round(Sim.REESeries[cnt,chain,id]))] += 1
            end
        end
    end
    Sim.REEHists /= length(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)
end =#

function computeCOP!(Sim::HPSAnalysis.SimData{R,I}) where {R<:Real, I<:Integer}
    Sim.COP    =  zeros(Sim.NChains,3  , Sim.NSteps)
    Sim.COP_uw =  zeros(Sim.NChains,3  , Sim.NSteps)
    
    for step in 1:Sim.NSteps
        for (chain, start, stop) in zip(collect(1:Sim.NChains), Sim.ChainStart, Sim.ChainStop)
            Sim.COP[chain,1, step]+= sum(Sim.x[start:stop, step])
            Sim.COP[chain,2, step]+= sum(Sim.y[start:stop, step])
            Sim.COP[chain,3, step]+= sum(Sim.z[start:stop, step])
            Sim.COP_uw[chain,1, step]+= sum(Sim.x_uw[start:stop, step])
            Sim.COP_uw[chain,2, step]+= sum(Sim.y_uw[start:stop, step])
            Sim.COP_uw[chain,3, step]+= sum(Sim.z_uw[start:stop, step])
        end
        Sim.COP[:,1,step] ./= Sim.ChainLength
        Sim.COP[:,2,step] ./= Sim.ChainLength
        Sim.COP[:,3,step] ./= Sim.ChainLength

        Sim.COP_uw[:,1,step] ./= Sim.ChainLength
        Sim.COP_uw[:,2,step] ./= Sim.ChainLength
        Sim.COP_uw[:,3,step] ./= Sim.ChainLength
    end
end

@doc raw"""
    computeCOM!(Sim::HPSAnalysis.SimData{R,I}) where {R<:Real, I<:Integer}

Computes the center of mass (COM) for each polymer chain.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Updates**:
- `Sim.COM`: COM of each chain.
- `Sim.COM_uw`: Unwrapped COM.
"""
function computeCOM!(Sim::HPSAnalysis.SimData{R,I}) where {R<:Real, I<:Integer}
    Sim.COM    =  zeros(R, Sim.NChains,3  , Sim.NSteps)
    Sim.COM_uw =  zeros(R, Sim.NChains,3  , Sim.NSteps)
    
    tmp  = zeros(R,  Sim.NAtoms,3 )
    tmp2 = zeros(R,  Sim.NAtoms,3)
    for step in 1:Sim.NSteps
        @inbounds tmp[:,1]  .=  Sim.x_uw[:, step].*Sim.Masses
        @inbounds tmp[:,2]  .=  Sim.y_uw[:, step].*Sim.Masses
        @inbounds tmp[:,3]  .=  Sim.z_uw[:, step].*Sim.Masses
        @inbounds tmp2[:,1] .=  Sim.x[:, step].*Sim.Masses
        @inbounds tmp2[:,2] .=  Sim.y[:, step].*Sim.Masses
        @inbounds tmp2[:,3] .=  Sim.z[:, step].*Sim.Masses

        for (chain, start, stop) in zip(collect(1:Sim.NChains), Sim.ChainStart, Sim.ChainStop)
            Sim.COM[chain,1, step]+= sum(tmp2[start:stop,1])
            Sim.COM[chain,2, step]+= sum(tmp2[start:stop,2])
            Sim.COM[chain,3, step]+= sum(tmp2[start:stop,3])
            Sim.COM_uw[chain,1, step]+= sum(tmp[start:stop,1])
            Sim.COM_uw[chain,2, step]+= sum(tmp[start:stop,2])
            Sim.COM_uw[chain,3, step]+= sum(tmp[start:stop,3])
        end
        Sim.COM[:,1,step] ./= Sim.ChainMasses
        Sim.COM[:,2,step] ./= Sim.ChainMasses
        Sim.COM[:,3,step] ./= Sim.ChainMasses

        Sim.COM_uw[:,1,step] ./= Sim.ChainMasses
        Sim.COM_uw[:,2,step] ./= Sim.ChainMasses
        Sim.COM_uw[:,3,step] ./= Sim.ChainMasses
    end
end

function computeWeightDistribution(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    ### Resolution in Angstrom
    
    ### 3 maps for x-y, x-z, y-z
    Sim.OffsetInd = zeros(eltype(Sim.NSteps),3,2)
    Sim.OffsetInd[1,1]=Int32(ceil(Sim.BoxSize[1,1]/Sim.Resolution))+1
    Sim.OffsetInd[1,2]=Int32(ceil(Sim.BoxSize[1,2]/Sim.Resolution))
    Sim.OffsetInd[2,1]=Int32(ceil(Sim.BoxSize[2,1]/Sim.Resolution))+1
    Sim.OffsetInd[2,2]=Int32(ceil(Sim.BoxSize[2,2]/Sim.Resolution))
    Sim.OffsetInd[3,1]=Int32(ceil(Sim.BoxSize[3,1]/Sim.Resolution))+1
    Sim.OffsetInd[3,2]=Int32(ceil(Sim.BoxSize[3,2]/Sim.Resolution))

    Sim.WeightMapXY =OffsetArray(zeros(eltype(Sim.x), Int32(ceil((Sim.BoxSize[1,2]-Sim.BoxSize[1,1])/Sim.Resolution))  , Int32(ceil((Sim.BoxSize[2,2]-Sim.BoxSize[2,1])/Sim.Resolution))  , Int32(Sim.NSteps)), Int32(ceil(Sim.BoxSize[1,1]/Sim.Resolution))+1:Int32(ceil(Sim.BoxSize[1,2]/Sim.Resolution)),   Int32(ceil(Sim.BoxSize[2,1]/Sim.Resolution))+1:Int32(ceil(Sim.BoxSize[2,2]/Sim.Resolution))  ,1:Sim.NSteps)
    Sim.WeightMapXZ =OffsetArray(zeros(eltype(Sim.x), Int32(ceil((Sim.BoxSize[1,2]-Sim.BoxSize[1,1])/Sim.Resolution))  , Int32(ceil((Sim.BoxSize[3,2]-Sim.BoxSize[3,1])/Sim.Resolution)) ,  Int32(Sim.NSteps)), Int32(ceil(Sim.BoxSize[1,1]/Sim.Resolution))+1:Int32(ceil(Sim.BoxSize[1,2]/Sim.Resolution)),   Int32(ceil(Sim.BoxSize[3,1]/Sim.Resolution))+1:Int32(ceil(Sim.BoxSize[3,2]/Sim.Resolution)),  1:Sim.NSteps);
    Sim.WeightMapYZ =OffsetArray(zeros(eltype(Sim.x), Int32(ceil((Sim.BoxSize[2,2]-Sim.BoxSize[2,1])/Sim.Resolution))  , Int32(ceil((Sim.BoxSize[3,2]-Sim.BoxSize[3,1])/Sim.Resolution)) , Int32(Sim.NSteps) ),  Int32(ceil(Sim.BoxSize[2,1]/Sim.Resolution))+1:Int32(ceil(Sim.BoxSize[2,2]/Sim.Resolution)),   Int32(ceil(Sim.BoxSize[3,1]/Sim.Resolution))+1:Int32(ceil(Sim.BoxSize[3,2]/Sim.Resolution))  , 1:Sim.NSteps);

    for step in 1:Sim.NSteps
        for atom in 1:Sim.NAtoms
            xind = Sim.x[atom,step] > Sim.BoxSize[1,1] && Sim.x[atom,step] < Sim.BoxSize[1,2] ? Int32(ceil(Sim.x[atom,step]/Sim.Resolution)) : continue
            yind = Sim.y[atom,step] > Sim.BoxSize[2,1] && Sim.y[atom,step] < Sim.BoxSize[2,2] ?  Int32(ceil(Sim.y[atom,step]/Sim.Resolution)) : continue
            zind = Sim.z[atom,step] > Sim.BoxSize[3,1] && Sim.z[atom,step] < Sim.BoxSize[3,2] ? Int32(ceil(Sim.z[atom,step]/Sim.Resolution)) :  continue

            Sim.WeightMapXY[xind , yind, step] += Sim.Masses[atom]
            Sim.WeightMapXZ[xind , zind, step] += Sim.Masses[atom]
            Sim.WeightMapYZ[yind , zind, step] += Sim.Masses[atom]
        end
    end
end

function computeRGComponentSlabHist(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    printstyled("Function computeRGComponentSlabHist returns the result and doesnt save it!\n"; color=:red)
    if Sim.SlabAxis==1
        SlabCoord = Sim.COM
    elseif Sim.SlabAxis==2
        SlabCoord = Sim.COM
    elseif Sim.SlabAxis == 3
        SlabCoord = Sim.COM
    else 
        ArgumentError("SlabAxis is not properly specified.")
    end
     AllAxis = [1,2,3]
    deleteat!(AllAxis, Sim.SlabAxis)
    
    RGCompHist = OffsetArray(zeros(eltype(Sim.x), Int32(ceil((Sim.BoxSize[Sim.SlabAxis,2]-Sim.BoxSize[Sim.SlabAxis,1])/Sim.Resolution)) , 3), Int32(ceil(Sim.BoxSize[Sim.SlabAxis,1]/Sim.Resolution))+1:Int32(ceil(Sim.BoxSize[Sim.SlabAxis,2]/Sim.Resolution)) ,  1:3)

    lowestind = Int32(ceil(Sim.BoxSize[Sim.SlabAxis,1]/Sim.Resolution))+1
    highestind = Int32(ceil(Sim.BoxSize[Sim.SlabAxis,2]/Sim.Resolution)) 
    
    volume =((Sim.BoxSize[AllAxis[1],2]-Sim.BoxSize[AllAxis[1],1])*(Sim.BoxSize[AllAxis[2],2]-Sim.BoxSize[AllAxis[2],1]))
    conversion = 1.66053906660/volume/Sim.Resolution
    for id in 1:3
        for step in Sim.EquilibrationTime:Sim.NSteps
            for chain in 1:Sim.NChains
                ind = Int32(ceil(SlabCoord[chain,2,step]/Sim.Resolution))
                ind = ind >= lowestind && ind<= highestind ?  ind : continue
                RGCompHist[ind, id]+=Sim.RGComponentSeries[chain, step, id]
            end
        end
    end

    RGCompHist ./= (Sim.NChains*(Sim.NSteps-Sim.EquilibrationTime))
    return RGCompHist
end 

function computeMSDofChains!(Sim::HPSAnalysis.SimData{R,I}; MaxDelta=100_000) where {R<:Real, I<:Integer}
    N = min(Sim.NSteps-Sim.EquilibrationTime, MaxDelta)
    Sim.MSD = zeros(R, Sim.NChains, N,3 )
    if length(Sim.COM)==0
        computeCOM!(Sim)
    end

    xdiff =zeros(R, Sim.NChains)
    ydiff =zeros(R, Sim.NChains)
    zdiff =zeros(R, Sim.NChains)


    step = Sim.EquilibrationTime

    for step in Sim.EquilibrationTime:Sim.NSteps
        for step2 in step+1:min(Sim.NSteps, step+MaxDelta)
            Δt = step2 - step
            xdiff .= (Sim.COM_uw[:,1,step2].- Sim.COM_uw[:,1,step]).^2
            ydiff .= (Sim.COM_uw[:,2,step2].- Sim.COM_uw[:,2,step]).^2
            zdiff .= (Sim.COM_uw[:,3,step2].- Sim.COM_uw[:,3,step]).^2
            
            Sim.MSD[:,Δt,1] .+= xdiff
            Sim.MSD[:,Δt,2] .+= ydiff
            Sim.MSD[:,Δt,3] .+= zdiff
        end
    end
    ### old technique
    #=for step in Sim.EquilibrationTime:Sim.NSteps
        for step2 in step+1:min(Sim.NSteps, step+MaxDelta)
            Δt = step2 - step
            for chain in 1:Sim.NChains
                xdiff[chain] = (Sim.COM[chain,1,step2]- Sim.COM[chain,1,step])%Sim.BoxLength[1]
                ydiff[chain] = (Sim.COM[chain,2,step2]- Sim.COM[chain,2,step])%Sim.BoxLength[2]
                zdiff[chain] = (Sim.COM[chain,3,step2]- Sim.COM[chain,3,step])%Sim.BoxLength[3]

                xdiff[chain] = xdiff[chain] - Sim.BoxLength[1]*(xdiff[chain]÷(Sim.BoxLength[1]/2.))
                ydiff[chain] = ydiff[chain] - Sim.BoxLength[2]*(ydiff[chain]÷(Sim.BoxLength[2]/2.))
                zdiff[chain] = zdiff[chain] - Sim.BoxLength[3]*(zdiff[chain]÷(Sim.BoxLength[3]/2.))

                xdiff[chain] = xdiff[chain]^2
                ydiff[chain] = ydiff[chain]^2
                zdiff[chain] = zdiff[chain]^2
            end
            
            Sim.MSD[:,Δt,1] .+= xdiff
            Sim.MSD[:,Δt,2] .+= ydiff
            Sim.MSD[:,Δt,3] .+= zdiff
        end
    end =#
    for step in axes(Sim.MSD,2)
        #for chain in axes(Sim.MSD, 1)
            Sim.MSD[:,step,1] ./=  (N-step+1) #Sim.MSD[chain,step,1]/(N-step+1)
            Sim.MSD[:,step,2] ./=  (N-step+1) #Sim.MSD[chain,step,2]/(N-step+1)
            Sim.MSD[:,step,3] ./=  (N-step+1) # Sim.MSD[chain,step,3]/(N-step+1)
        #end
    end
end

function unfoldPositions(Sim::SimData{R,I}; CheckFiles=true) where {R<:Real, I<:Integer}

    if CheckFiles && length(Sim.TrajectoryFile)>=3 
        if Sim.TrajectoryFile[end-2:end] ==".h5" || Sim.TrajectoryFile[end-2:end] =="gsd"
            return nothing 
        end 
    end

    ### allocate data on hard drive
    Sim.x_uw_io= open(Sim.x_uw_FilePath,"w+")
    Sim.y_uw_io= open(Sim.y_uw_FilePath,"w+")
    Sim.z_uw_io= open(Sim.z_uw_FilePath,"w+")

    Sim.x_uw =  Mmap.mmap(Sim.x_uw_io, Matrix{R}, (Sim.NAtoms,Sim.NSteps))
    Sim.y_uw =  Mmap.mmap(Sim.y_uw_io, Matrix{R}, (Sim.NAtoms,Sim.NSteps))
    Sim.z_uw =  Mmap.mmap(Sim.z_uw_io, Matrix{R}, (Sim.NAtoms,Sim.NSteps))

    ### initialise with wrapped data
    Sim.x_uw[:,1] .= Sim.x[:,1] # = deepcopy(Sim.x)
    Sim.y_uw[:,1] .= Sim.y[:,1] # .= deepcopy(Sim.y)
    Sim.z_uw[:,1] .= Sim.z[:,1] # .= deepcopy(Sim.z)

    ### do the actual work
    Sim.CellStep = zeros(1)
    for step in 1:Sim.NSteps ### changed that one to 1 .... was 2 before
        Sim.CellStep[1]=step
        checkCIPositions(Sim)
    end
    #unfoldPositions_SD(Sim)

    ### sync RAM to disk before closing
    Mmap.sync!(Sim.x_uw)
    Mmap.sync!(Sim.y_uw)
    Mmap.sync!(Sim.z_uw)
    close(Sim.x_uw_io)
    close(Sim.y_uw_io)
    close(Sim.z_uw_io)
end

function checkCIPositions(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    step = Sim.CellStep[1]
    too_large= 15.
 
    diff = zeros(eltype(Sim.x), 3)
    diff2 = zeros(eltype(Sim.x), 3)

    COM_new = zeros(eltype(Sim.x), 3)
    COM_old = zeros(eltype(Sim.x), 3)
    dist = 0

    for chain in 1:Sim.NChains
        start = Sim.ChainStart[chain]
        stop = Sim.ChainStop[chain]

        Sim.x_uw[start, step]=Sim.x[start, step]
        Sim.y_uw[start, step]=Sim.y[start, step]
        Sim.z_uw[start, step]=Sim.z[start, step]

        for atom in start:(stop-1)
            diff[1] = Sim.x[atom+1, step]- Sim.x_uw[atom, step]
            diff[2] = Sim.y[atom+1, step]- Sim.y_uw[atom, step]
            diff[3] = Sim.z[atom+1, step]- Sim.z_uw[atom, step]

            diff2.= diff

            if abs(diff[1]-Sim.BoxLength[1])<abs(diff[1]) || abs(diff[1]+Sim.BoxLength[1])<abs(diff[1])
                Sim.x_uw[atom+1, step] = Sim.x[atom+1, step] - Sim.BoxLength[1]*(convert(eltype(Sim.NChains),(round((diff[1]/(Sim.BoxLength[1]/2.)))/2.)))
                diff2[1] = Sim.x_uw[atom+1, step]- Sim.x_uw[atom, step]
            else
                Sim.x_uw[atom+1, step]=Sim.x[atom+1, step]
            end            
            
            if abs(diff[2]-Sim.BoxLength[2])<abs(diff[2]) || abs(diff[2]+Sim.BoxLength[2])<abs(diff[2])
                Sim.y_uw[atom+1, step] = Sim.y[atom+1, step]- Sim.BoxLength[2]*(convert(eltype(Sim.NChains),(round((diff[2]/(Sim.BoxLength[2]/2.)))/2.)))
                diff2[2] = Sim.y_uw[atom+1, step]- Sim.y_uw[atom, step]
            else
                Sim.y_uw[atom+1, step]=Sim.y[atom+1, step]
            end      

            if abs(diff[3]-Sim.BoxLength[3])<abs(diff[3]) || abs(diff[3]+Sim.BoxLength[3])<abs(diff[3])
                Sim.z_uw[atom+1, step] = Sim.z[atom+1, step] - Sim.BoxLength[3]*(convert(eltype(Sim.NChains),(round((diff[3]/(Sim.BoxLength[3]/2.)))/2.)))
                diff2[3] = Sim.z_uw[atom+1, step]- Sim.z_uw[atom, step]
            else
                Sim.z_uw[atom+1, step]=Sim.z[atom+1, step]
            end         

            dist = sqrt(diff2[1]^2+diff2[2]^2+diff2[3]^2)
            if dist >too_large
                printstyled("CI Positions didnt work!!!!\n"; color=:red)
                println(Sim.BoxLength)
                println("Atom: $atom, Step: $step,  Diff: $diff,Diff2: $diff2, Dist: $dist ")
                println("x $(Sim.x_uw[atom, step]), $(Sim.x_uw[atom+1, step]) ")
                println("y $(Sim.y_uw[atom, step]), $(Sim.y_uw[atom+1, step]) ")
                println("z $(Sim.z_uw[atom, step]), $(Sim.z_uw[atom+1, step]) \n")

                println("x $(Sim.x[atom, step]), $(Sim.x[atom+1, step]) ")
                println("y $(Sim.y[atom, step]), $(Sim.y[atom+1, step]) ")
                println("z $(Sim.z[atom, step]), $(Sim.z[atom+1, step]) ")
                #return
            end
        end

        if step>1
        COM_old[1] = sum(Sim.x_uw[start:stop, step-1])
        COM_old[2] = sum(Sim.y_uw[start:stop, step-1])
        COM_old[3] = sum(Sim.z_uw[start:stop, step-1])
        COM_old ./=(Sim.ChainStop[chain]-Sim.ChainStart[chain])

        COM_new[1] = sum(Sim.x_uw[start:stop, step])
        COM_new[2] = sum(Sim.y_uw[start:stop, step])
        COM_new[3] = sum(Sim.z_uw[start:stop, step])
        COM_new ./=(Sim.ChainStop[chain]-Sim.ChainStart[chain])

        diff .= COM_new.-COM_old
        if abs(diff[1]-Sim.BoxLength[1])<abs(diff[1]) || abs(diff[1]+Sim.BoxLength[1])<abs(diff[1]) 
            #Sim.x_uw[start:stop, step] .-=  Sim.BoxLength[1]*(convert(eltype(Sim.NChains),(round( (diff[1]/(Sim.BoxLength[1]/2.)) /2.))))
            Sim.x_uw[start:stop, step] .-=  Sim.BoxLength[1]*(round(diff[1]/(Sim.BoxLength[1]/2.0))÷2)
            COM_new[1] = sum(Sim.x_uw[start:stop, step])
        end            
        
        if abs(diff[2]-Sim.BoxLength[2])<abs(diff[2]) || abs(diff[2]+Sim.BoxLength[2])<abs(diff[2])
            #Sim.y_uw[start:stop, step] .-=  Sim.BoxLength[2]*(convert(eltype(Sim.NChains),(round( (diff[2]/(Sim.BoxLength[2]/2.)) /2.))))
            Sim.y_uw[start:stop, step] .-=  Sim.BoxLength[2]*(round(diff[2]/(Sim.BoxLength[2]/2.)) ÷2)
            COM_new[2] = sum(Sim.y_uw[start:stop, step])
        end         
        
        if abs(diff[3]-Sim.BoxLength[3])<abs(diff[3]) || abs(diff[3]+Sim.BoxLength[3])<abs(diff[3])
            #Sim.z_uw[start:stop, step] .-=  Sim.BoxLength[3]*(convert(eltype(Sim.NChains),(round((diff[3]/(Sim.BoxLength[3]/2.))/2.))))
            Sim.z_uw[start:stop, step] .-=  Sim.BoxLength[3]*(round(diff[3]/(Sim.BoxLength[3]/2.)) ÷2)
            COM_new[3] = sum(Sim.z_uw[start:stop, step])
        end 
      
        COM_new ./=(Sim.ChainStop[chain]-Sim.ChainStart[chain])
        
        diff .= COM_new.-COM_old

        dist = sqrt(diff[1]^2+diff[2]^2+diff[3]^2)
        end
    end
end

function computeBondsHists(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Sim.x_uw_io= open(Sim.x_uw_FilePath,"r+")
    Sim.y_uw_io= open(Sim.y_uw_FilePath,"r+")
    Sim.z_uw_io= open(Sim.z_uw_FilePath,"r+")

    Sim.x_uw =  Mmap.mmap(Sim.x_uw_io, Matrix{eltype(Sim.x)}, (Sim.NAtoms,Sim.NSteps))
    Sim.y_uw =  Mmap.mmap(Sim.y_uw_io, Matrix{eltype(Sim.x)}, (Sim.NAtoms,Sim.NSteps))
    Sim.z_uw =  Mmap.mmap(Sim.z_uw_io, Matrix{eltype(Sim.x)}, (Sim.NAtoms,Sim.NSteps))

    Resolution=100.
    Sim.BondHist = zeros(eltype(Sim.x), 1000)
    diff= 0.
    ind = Int(0)
    for step in 1:Sim.NSteps
        for chain in 1:Sim.NChains
            for atom in Sim.ChainStart[chain]:(Sim.ChainStop[chain]-1)
                diff =(Sim.x_uw[atom+1, step]- Sim.x_uw[atom, step])^2
                diff+=(Sim.y_uw[atom+1, step]- Sim.y_uw[atom, step])^2
                diff+=(Sim.z_uw[atom+1, step]- Sim.z_uw[atom, step])^2

                ind = Int32(ceil(sqrt(diff)*Resolution))

                if ind==0
                    continue
                end
                if ind<1000
                    Sim.BondHist[ind] += 1
                end
            end
        end
    end
    for ind in axes(Sim.BondHist,1)
        Sim.BondHist[ind] /= ((Sim.NSteps-1)*(Sim.NAtoms-Sim.NChains))/Resolution
    end
end

@inline function sqr_distance(Sim::SimData{R,I}, a::Vector{R}, b::Vector{R}) where {R<:Real, I<:Integer}
    inv = 1.0 ./Sim.BoxLength
    Δ = a.-b
    Δ .-= floor.(I, Δ.*inv .+0.5).*Sim.BoxLength
    return sum(Δ.^2)
end

distance(Sim::SimData{R,I}, a::Vector{R}, b::Vector{R}) where {R<:Real, I<:Integer} =  sqrt(sqr_distance(Sim, a, b))

#### probably deprecated check sequence analysis and the sequence monte carlo
function computeSequenceHydropathyDecoration(Sim::SimData{R,I}; SimType="Calvados2", β=-1.0) where {R<:Real, I<:Integer}
    Sim.HydropathyDecoration = zeros(R, Sim.NChains)
    Sim.MeanHydropathy = zeros(R, Sim.NChains)
    IDToLambda = Dict()
    if SimType=="Calvados2"
        for id in keys(Sim.IDToResName)
            IDToLambda[id]= BioData.OneToCalvados2Lambda[Sim.IDToResName[id][1]]
        end
    end
    λ = zeros(R, Sim.NAtoms)
    for res in 1:Sim.NAtoms
        λ[res] = IDToLambda[Sim.IDs[res]]
    end
    for chain in 1:Sim.NChains
        for res in (Sim.ChainStart[chain]+1):Sim.ChainStop[chain]
            λ_i =  λ[res]
            for res2 in 1:res-1
                Sim.HydropathyDecoration[chain] += (λ_i+λ[res2])*abs(res2-res)^β
            end
        end
        Sim.HydropathyDecoration[chain] /= Sim.ChainLength[chain] 
        Sim.MeanHydropathy[chain] = sum(λ[Sim.ChainStart[chain]:Sim.ChainStop[chain]])/Sim.ChainLength[chain]
    end
end



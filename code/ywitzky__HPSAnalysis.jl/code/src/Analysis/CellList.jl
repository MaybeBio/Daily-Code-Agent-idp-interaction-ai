const SCALE_FACTOR=1+10^-6

function initCellLists(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Sim.CellStep = zeros(I, 1)
    Sim.CellDimensions= zeros(3,2)
    for i in 1:3
            if !(-Sim.BoxSize[i,2] ≈ Sim.BoxSize[i,1])
                @error "Simulation box has to be centered at (0,0,0) to use CellLists!"
            end
            ### increase CellResolution slightly such that boxes always end at multiple of box size
            n = floor(I, Sim.BoxSize[i,2]/Sim.CellResolution[i])
            Sim.CellResolution[i] =  Sim.BoxSize[i,2]*SCALE_FACTOR/n ### increase ever so slightly such that hash(BoxSize) lands in the correct cell.
            Sim.CellDimensions[i,1] = -n+1 ### hash function is assymetric
            Sim.CellDimensions[i,2] = n
    end

    CellVolume=prod(Sim.CellResolution)
    Sim.MaxParticlesPerCell= I(ceil(CellVolume/(4/3*pi*(4.5/2.)^3) )+1)# volume of the smallest residues

    xrange = (Sim.CellDimensions[1,1]-1):(Sim.CellDimensions[1,2]+1) 
    yrange = (Sim.CellDimensions[2,1]-1):(Sim.CellDimensions[2,2]+1)
    zrange = (Sim.CellDimensions[3,1]-1):(Sim.CellDimensions[3,2]+1)

    Sim.CellList         = OffsetArray([Vector{I}() for x in xrange, y in yrange, z in zrange], xrange, yrange, zrange)
    Sim.PositiveCellList = OffsetArray([Vector{I}() for x in xrange, y in yrange, z in zrange], xrange, yrange, zrange)
    Sim.NegativeCellList = OffsetArray([Vector{I}() for x in xrange, y in yrange, z in zrange], xrange, yrange, zrange)
    Sim.ChainNumCellList = OffsetArray([Vector{I}() for x in xrange, y in yrange, z in zrange], xrange, yrange, zrange)

    map(x->sizehint!(x, Sim.MaxParticlesPerCell),  Sim.CellList)
    map(x->sizehint!(x, Sim.MaxParticlesPerCell),  Sim.PositiveCellList)
    map(x->sizehint!(x, Sim.MaxParticlesPerCell),  Sim.NegativeCellList)
    map(x->sizehint!(x, Sim.MaxParticlesPerCell),  Sim.ChainNumCellList)

    ### allocate work arrays for later function calls
    Sim.CenterBox = zeros(3*Sim.MaxParticlesPerCell  ,3)
    Sim.NeighBox  = zeros(3*Sim.MaxParticlesPerCell  ,3)
    Sim.CL_Dist   = zeros(3*Sim.MaxParticlesPerCell^2,3)

    resetCellLists(Sim)
end

function resetCellLists(Sim::SimData{T,I}) where {T<:Real, I<:Integer} 
    map(x->resize!(x,0), Sim.CellList)
    map(x->resize!(x,0), Sim.PositiveCellList)
    map(x->resize!(x,0), Sim.NegativeCellList)
    map(x->resize!(x,0), Sim.ChainNumCellList)
end

function destroyCellLists(Sim::SimData{T,I}) where {T<:Real, I<:Integer}
    Sim.CellList         = zeros(0)
    Sim.PositiveCellList = zeros(0)
    Sim.NegativeCellList = zeros(0)
    Sim.ChainNumCellList = zeros(0)
end

function computeCellLists(Sim::SimData{T,I}; ComputeCharges=true::Bool, ComputeChains=false::Bool) where {T<:Real, I<:Integer}
    step = Sim.CellStep[1]
    resetCellLists(Sim)

    x_low, x_high, y_low, y_high, z_low, z_high = computeCellMaxima(Sim)

    xdiff =convert(I,1+x_high-x_low)
    ydiff =convert(I,1+y_high-y_low)
    zdiff =convert(I,1+z_high-z_low)

    ind  = 1
    xind = 1
    yind = 1
    zind = 1

    xind_arr = LoopVectorization.@avx ceil.(I,Sim.x[:,step]/Sim.CellResolution[1])
    yind_arr = LoopVectorization.@avx ceil.(I,Sim.y[:,step]/Sim.CellResolution[2])
    zind_arr = LoopVectorization.@avx ceil.(I,Sim.z[:,step]/Sim.CellResolution[3])
    
    atom_cnt = 0
    for (C, (start,stop)) in enumerate(zip(Sim.ChainStart, Sim.ChainStop))
        @inbounds for atom in start:stop #1:Sim.NAtoms
                
            xind = xind_arr[atom]

            yind = yind_arr[atom]
            
            zind = zind_arr[atom]

            ### all atoms in cell 
            push!(Sim.CellList[xind, yind, zind], atom)

            if ComputeCharges
                ### all negative residues
                if Sim.Charges[atom]<0
                    push!(Sim.NegativeCellList[xind, yind, zind], atom)
                end

                ### all positives residues
                if Sim.Charges[atom]>0
                    push!(Sim.PositiveCellList[xind, yind, zind], atom)
                end
            end 

            if ComputeChains
                push!(Sim.ChainNumCellList[xind, yind, zind],C)
            end
        end
    end


    ### copy borders on the other side...
    CellLists = Vector{OffsetArray{Vector{I}}}([Sim.CellList])

    if ComputeCharges
        append!(CellLists,[Sim.NegativeCellList,Sim.PositiveCellList])
    end

    if ComputeChains
        append!(CellLists,[Sim.ChainNumCellList])
    end

    if ComputeCharges
        @error "ComputeCellLists is currently not correct for ComputeCharges"
    end

    CellList = CellLists[1]
    @inbounds begin 
        for yind in y_low-1:y_high+1
            for zind in z_low-1:z_high+1
                for (i,e) in enumerate(CellList[x_low, yind, zind])
                    ### copy -x to x and back
                    if ~(e in CellList[x_high+1, yind, zind])
                        for CL in CellLists ### iterating over all at once allows to detect duplicates based on atom number while still being able to have duplicates of chain numbers
                            push!(CL[x_high+1, yind, zind], CL[x_low, yind, zind][i] )
                        end
                    end
                end

                for (i,e) in enumerate(CellList[x_high, yind, zind])
                    ### copy x to -x and back
                    if ~(e in CellList[x_low-1, yind, zind])
                        for CL in CellLists
                            push!( CL[x_low-1, yind, zind],  CL[x_high, yind, zind][i])
                        end
                    end
                end
            end
        end

        for xind in x_low-1:x_high+1
            for zind in z_low-1:z_high+1
                for (i,e) in enumerate(CellList[xind, y_low, zind])
                    ### copy -y to y and back
                    if ~(e in CellList[xind, y_high+1, zind])
                        for CL in CellLists
                            push!( CL[xind, y_high+1, zind], CL[xind, y_low, zind][i])
                        end
                    end
                end

                for (i,e) in enumerate(CellList[xind, y_high, zind])
                    ### copy y to -y and back
                    if ~(e in CellList[xind, y_low-1, zind])
                        for CL in CellLists
                            push!(CL[xind, y_low-1, zind], CL[xind, y_high, zind][i])
                        end
                    end
                end
            end
        end

        for xind in x_low-1:x_high+1
            for yind in y_low-1:y_high+1
                for (i,e) in enumerate(CellList[xind, yind, z_low])
                    ### copy -z to z and back
                    if ~(e in CellList[xind, yind, z_high+1])
                        for CL in CellLists
                            push!(CL[xind, yind, z_high+1], CL[xind, yind, z_low][i])
                        end
                    end
                end

                for (i,e) in enumerate(CellList[xind, yind, z_high])
                    ### copy z to -z and back
                    if ~(e in CellList[xind, yind, z_low-1])
                        for CL in CellLists
                            push!(CL[xind, yind, z_low-1], CL[xind, yind, z_high][i])
                        end
                    end
                end
            end
        end
    end
end

function computeCellMaxima(Sim::SimData{T,I}) where {T<:Real, I<:Integer}
    x_low  = ceil(I,Sim.BoxSize[1,1]/Sim.CellResolution[1])#+1
    x_high = ceil(I,Sim.BoxSize[1,2]/Sim.CellResolution[1])
    y_low  = ceil(I,Sim.BoxSize[2,1]/Sim.CellResolution[2])#+1
    y_high = ceil(I,Sim.BoxSize[2,2]/Sim.CellResolution[2])
    z_low  = ceil(I,Sim.BoxSize[3,1]/Sim.CellResolution[3])#+1
    z_high = ceil(I,Sim.BoxSize[3,2]/Sim.CellResolution[3])

    return x_low, x_high, y_low, y_high, z_low, z_high
end

function iterateThroughCellList(Sim::SimData{T,I}, func::Function; CellMaxima=nothing) where {T<:Real, I<:Integer}
    step = Sim.CellStep[1]
    if isnothing(CellMaxima)
        x_low, x_high, y_low, y_high, z_low, z_high = computeCellMaxima(Sim)
    else
        x_low, x_high, y_low, y_high, z_low, z_high = CellMaxima
    end

    if x_high-x_low < 3 || y_high-y_low < 3 || z_high-z_low < 3
        printstyled("Function iterateThroughCellList wont properly iterate through CellList because of too little cells in at least one direction. \n $(Sim.CellDimensions)"; color=:red)
        return
    end

    fill!(Sim.CenterBox,0)
    fill!(Sim.NeighBox ,0)
    fill!(Sim.CL_Dist  ,0)
    border=false

    for xind in x_low-1:x_high ### -1 allows do avoid double counting in subfunctions
        for yind in y_low-1:y_high
            for zind in z_low-1:z_high
                ### check whether follow up functions has to consider wrap behaviour
                border = xind==x_low || xind == x_high || yind==y_low || yind == y_high || zind==z_low || zind == z_high 
                func(Sim, xind, yind, zind,border)
                #=if border
                    if xind==x_low 
                        func(Sim,x_high, yind, zind,border)
                    elseif xind==x_high
                        func(Sim, x_low, yind, zind,border)
                    elseif yind==y_high
                        func(Sim, xind, y_low, zind,border)                
                    elseif yind==y_low
                        func(Sim,xind, y_high, zind,border)
                    elseif zind==z_high
                        func(Sim, xind, yind, z_low,border)
                    elseif zind==z_low
                        func(Sim,xind, yind, z_high,border)
                    end
                end=#
            end
        end
    end
end

function computeDistancesForCellNeighbour(BoxLength::Vector{T}, CenterBox::Matrix{T}, NeighBox::Matrix{T}, CL_Dist::Matrix{T}, IMax::I, JMax::I) where {I<:Integer, T<:Real}
    ### IMax = CenterBox, JMax = NeighBox
    inv = 1.0./BoxLength

    cnt = 1 
    #=if  border
        for i in 1:IMax
            for j in 1:JMax
                for dim in 1:3
                    CL_Dist[cnt,dim] = ((NeighBox[j,dim]-CenterBox[i,dim])%BoxLength[dim])
                    #CL_Dist[cnt,dim] = CL_Dist[cnt,dim] - BoxLength[dim]*CL_Dist[cnt,dim]÷(BoxLength[dim]/2.)
                end
                CL_Dist[cnt,1] = CL_Dist[cnt,1]^2 + CL_Dist[cnt,2]^2 + CL_Dist[cnt,3]^2
                CL_Dist[cnt,1] = sqrt(CL_Dist[cnt,1])
                cnt+=1
            end
        end
        else
        for i in 1:IMax
            for j in 1:JMax
                CL_Dist[cnt,1] = (NeighBox[j,1]-CenterBox[i,1])^2
                CL_Dist[cnt,1] = (NeighBox[j,2]-CenterBox[i,2])^2 + CL_Dist[cnt,1]
                CL_Dist[cnt,1] = (NeighBox[j,3]-CenterBox[i,3])^2 + CL_Dist[cnt,1]
                CL_Dist[cnt,1] = sqrt(CL_Dist[cnt,1])
                cnt+=1
            end
        end
    end =#
    #https://en.wikipedia.org/wiki/Periodic_boundary_conditions
    for i in 1:IMax
        @inbounds for j in 1:JMax
            LoopVectorization.@avx CL_Dist[cnt,1] = NeighBox[j,1]-CenterBox[i,1]
            LoopVectorization.@avx CL_Dist[cnt,2] = NeighBox[j,2]-CenterBox[i,2]
            LoopVectorization.@avx CL_Dist[cnt,3] = NeighBox[j,3]-CenterBox[i,3]

            LoopVectorization.@avx CL_Dist[cnt,1] -= floor(I, CL_Dist[cnt,1] * inv[1]+0.5)*BoxLength[1]
            LoopVectorization.@avx CL_Dist[cnt,2] -= floor(I, CL_Dist[cnt,2] * inv[2]+0.5)*BoxLength[2]
            LoopVectorization.@avx CL_Dist[cnt,3] -= floor(I, CL_Dist[cnt,3] * inv[3]+0.5)*BoxLength[3]

            LoopVectorization.@avx CL_Dist[cnt,1] = sqrt(CL_Dist[cnt,1]^2 +CL_Dist[cnt,2]^2 + CL_Dist[cnt,3]^2)
            cnt+=1
        end
    end
end


function computeSqrDistancesForCellNeighbour(BoxLength::Vector{T}, CenterBox::AbstractArray{T}, NeighBox::AbstractArray{T}, CL_Dist::Matrix{T}, IMax::I, JMax::I) where {I<:Integer, T<:Real}
    ### IMax = CenterBox, JMax = NeighBox
    inv = one(T)./BoxLength

    cnt = 1 

    #https://en.wikipedia.org/wiki/Periodic_boundary_conditions
    for i in 1:IMax
         @inbounds for j in 1:JMax
            CL_Dist[cnt,1] = NeighBox[j,1]-CenterBox[i,1]
            CL_Dist[cnt,2] = NeighBox[j,2]-CenterBox[i,2]
            CL_Dist[cnt,3] = NeighBox[j,3]-CenterBox[i,3]

            CL_Dist[cnt,1] -= floor(I, CL_Dist[cnt,1] * inv[1]+0.5)*BoxLength[1]
            CL_Dist[cnt,2] -= floor(I, CL_Dist[cnt,2] * inv[2]+0.5)*BoxLength[2]
            CL_Dist[cnt,3] -= floor(I, CL_Dist[cnt,3] * inv[3]+0.5)*BoxLength[3]

            CL_Dist[cnt,1] = CL_Dist[cnt,1]^2 +CL_Dist[cnt,2]^2 + CL_Dist[cnt,3]^2

            cnt+=1
        end
    end
end

function getUniqueCellsOfChain(Sim::SimData{T,I}, Chain::I2)where {T<:Real, I<:Integer, I2<:Integer}
    start = Sim.ChainStart[Chain]
    stop  = Sim.ChainStop[Chain]
    step = Sim.CellStep[1]
    ### results in out of bounds if  Sim.x[i,s]==Sim.BoxLength[..]
    @views xind_arr = ceil.(I,Sim.x[start:stop,step]/Sim.CellResolution[1])
    @views yind_arr = ceil.(I,Sim.y[start:stop,step]/Sim.CellResolution[2])
    @views zind_arr = ceil.(I,Sim.z[start:stop,step]/Sim.CellResolution[3])

    return unique!(collect(zip(xind_arr, yind_arr, zind_arr)))
end
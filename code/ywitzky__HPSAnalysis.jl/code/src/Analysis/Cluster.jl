@doc raw"""
    computeClustersByCOM(Sim::HPSAnalysis.SimData{T,Int}; Cutoff=50.0) where{T<:Real, Int<:Integer}

Computes clusters based on distances between protein COMs and cutoff.

Computes graph network for each frame in EquilibrationTime:RGMeasureStep:NSteps by naivly computing minimal distances between all proteins center of masses based on the given cutoff. Returns the concatenation of groups of IDs belonging to weekly connected graph components of all frames.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.
- `Cutoff::Float`: Distance criterion for a Cluster.

**Returns**:
- `Cluster::Vector{Vector{Vector{Int}}}`: List of Clusters of the chains.
"""
function computeClustersByChainCOM(Sim::HPSAnalysis.SimData{T,Int}; Cutoff=50.0) where{T<:Real, Int<:Integer}
    ### cutoff from dignon et al. Sequence determinants of protein phase behavior from a coarse-grained model
    Clusters = Vector{Vector{Vector{Int}}}()

    for (i,step) in enumerate(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)
        G = Graphs.SimpleGraph{Int}()
        add_vertices!(G, Sim.NChains)

        for C in 1:Sim.NChains
            dx = Sim.COM_uw[C+1:Sim.NChains,1, step] .- Sim.COM_uw[C,1, step]
            dy = Sim.COM_uw[C+1:Sim.NChains,2, step] .- Sim.COM_uw[C,2, step]
            dz = Sim.COM_uw[C+1:Sim.NChains,3, step] .- Sim.COM_uw[C,3, step]

            dx -= Sim.BoxLength[1].*round.(Int, dx./Sim.BoxLength[1])
            dy -= Sim.BoxLength[2].*round.(Int, dy./Sim.BoxLength[2])
            dz -= Sim.BoxLength[3].*round.(Int, dz./Sim.BoxLength[3])

            dist_sqr = dx.^2 .+dy.^2 .+dz.^2
            for (k,K) in enumerate(C+1:Sim.NChains)
                if dist_sqr[k]<Cutoff^2
                    add_edge!(G, C,K)
                end
            end
        end

        fill!(Sim.ChargeChainContactMatrix, 0)
        push!(Clusters,  weakly_connected_components(G))

    end

    Sim.ClusterRange = Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
    Sim.Clusters = Clusters
    return Clusters
end

function computeChargeClusters(Sim::HPSAnalysis.SimData{T,Int}; Cutoff=5.0) where{T<:Real, Int<:Integer}
    printstyled("Function computeClusters does only detect charge contacts for clusters.\n"; color=:yellow)
    ### abuse cellList mechanism developed for Charge contacts
    Sim.ChargeResidueContactMatrix = zeros(eltype(Sim.x), Sim.MaxChainLength,Sim.MaxChainLength)
    Sim.ChargeChainContactMatrix = zeros(eltype(Sim.x), Sim.NChains,Sim.NChains)
    Clusters = Vector{Vector{Vector{Int}}}()

    pre = Sim.ChargeAnalysisCutoff 
    Sim.ChargeAnalysisCutoff = Cutoff
    initCellLists(Sim)
    for (i,step) in enumerate(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)
        Sim.ChargeContactValid = Dict()
        
        G = Graphs.SimpleGraph{Int}()
        add_vertices!(G, Sim.NChains)
        
        Sim.CellStep[1] = step
        HPSAnalysis.resetCellLists(Sim) ### empty all cells
        if i%100==0
            println("Step: $(Sim.CellStep) of $(Sim.NSteps)")
        end
        HPSAnalysis.computeCellLists(Sim) ### compute CellLists  for step

        iterateThroughCellList(Sim, computeChargeBondTimesForCell )

        for key in keys(Sim.ChargeContactValid)
            atom_i = key[1] ### always positive charge
            atom_j = key[2] ### always negative charge

            ### compute charge bond Matrix

            ### compute chain numbers
            chain_num_i = Int32(ceil(atom_i/Sim.MaxChainLength))
            chain_num_j = Int32(ceil(atom_j/Sim.MaxChainLength))

            ### compute residue number in chain
            id_in_chain_i = atom_i - (chain_num_i-1)*Sim.MaxChainLength
            id_in_chain_j = atom_j - (chain_num_j-1)*Sim.MaxChainLength
            
            ### increase Contactmatrices
            Sim.ChargeResidueContactMatrix[id_in_chain_i, id_in_chain_j] += 1
            Sim.ChargeChainContactMatrix[chain_num_i, chain_num_j] += 1
        end

        for C in 1:Sim.NChains
            for K in C+1:Sim.NChains
                if Sim.ChargeChainContactMatrix[C,K]>0
                    add_edge!(G, C,K)
                end
            end
        end

        fill!(Sim.ChargeChainContactMatrix, 0)
        push!(Clusters,  weakly_connected_components(G))

    end
    Sim.ChargeAnalysisCutoff = pre


    Sim.ClusterRange = Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
    Sim.Clusters = Clusters
    return Clusters
end


@doc raw"""
    getLargestClusterIDs(Sim::SimData{T,I}) where {I<:Integer, T<:Real}

Computes IDs for the largest clusters.

Warning is emitted if less than 80% of the proteins are within the largest cluster. The cluster definition might not be reasonable when being used to center slab simulations or other tasks.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Returns**:
- `LClustID::Vector{Int}`: List of IDs for largest clusters.
"""
function getLargestClusterIDs(Sim::SimData{T,I} ) where {I<:Integer, T<:Real}
    LClustID = argmax.([length.(C) for C in Sim.Clusters]) ### returns index for clusters with most proteins inside

    LargeClusterLengths = [length(Sim.Clusters[step][LId]) for (step, LId) in enumerate(LClustID)]

    if sum(LargeClusterLengths.< Sim.NChains*0.8)>0
        @warn "⚠️ The number of proteins within the Largest Cluster is less than 80% of the proteins in the simulation. Be cautious when using this cluster definition for slab centering etc. ."
    end

    return LClustID
end

@doc raw"""
    computeCOMOfLargestCluster(Sim::SimData{T,I} ) where {I<:Integer, T<:Real}

Computes COM of the largest cluster along the axis Sim.SlabAxis in the unwrapped positions.

For the algorithm to work, it is important that the dense phase of the first frame is not crossing the periodic boundaries. The algorithm take dense phase center of the previous step to pre center the data to about issues  at the boundaries. Therefore it works only when the displacement along Sim.SlabAxis between consecutive frames is smaller than half of the box length.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Returns**:
- `Coms::Vector{Float}`: COM of largest cluster along the Sim.SlabAxis for each step.
"""
function computeCOMOfLargestCluster(Sim::SimData{T,I} ) where {I<:Integer, T<:Real}
    Coms= zeros(T, length(Sim.ClusterRange))#Vector{Int}
    LClustID = getLargestClusterIDs(Sim)

    Len = Sim.BoxLength[Sim.SlabAxis]
    Len_inv = one(T)/Sim.BoxLength[Sim.SlabAxis]

    mass= zero(T)
    LId = LClustID[1]
    for CID in Sim.Clusters[1][LId] ### compute first where everything is (hopefully) nicely wrapped
        Coms[1] += Sim.COM[CID,Sim.SlabAxis,1]*Sim.ChainMasses[CID]
        mass+= Sim.ChainMasses[CID]
    end
    Coms/= mass

    for (j,step) in enumerate(Sim.ClusterRange)
        LId=LClustID[j]
        if j==1 continue end
        mass= zero(T)
        for CID in Sim.Clusters[j][LId] 
            ### center according to history; avoids centering issues when dense phase cross pbc
            ### wrap into central image, once issues are avoided
            tmp = (Sim.COM_uw[CID,Sim.SlabAxis,step]-Coms[j-1])
            tmp -= Len*round(I, tmp*Len_inv)

            Coms[j] += tmp  *Sim.ChainMasses[CID]

            mass+= Sim.ChainMasses[CID]
        end
        Coms[j] /= mass 
        Coms[j] += Coms[j-1] ### add previous position get unwrapped position
    end

    return Coms
end

@doc raw"""
    computeCOMsOfCluster(Sim::SimData{T,I} ) where {I<:Integer, T<:Real}

Computes COM for each cluster in each frame.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Returns**:
- `Coms::Matrix{Float}`: Matrix of COM for each cluster and time step.
"""
function computeCOMsOfCluster(Sim::SimData{T,I} ) where {I<:Integer, T<:Real}
    Coms= Vector{Matrix{T}}()
    for (Cstep,C) in enumerate(Sim.Clusters)
        step = Sim.ClusterRange[Cstep]
        tmp = zeros(length(C),3 )
        for (CID,Cluster) in enumerate(C)
            mass=0
            for id in Cluster
                tmp[CID,:].+= Sim.COM[id,:,step]*Sim.ChainMasses[CID]
                mass+= Sim.ChainMasses[CID]
            end
            tmp[CID,:]./=mass
        end
        push!(Coms, tmp)
    end
    return Coms
end

@inline function isOverlapping1D(xmin1,xmax1,xmin2,xmax2)
    return xmax1 >= xmin2 && xmax2 >= xmin1
end

@inline function isOverlapping3D(minima::Array{R}, maxima::Array{R}, I::Int,J::Int, BoxLength::Vector{R}) where {Int<:Integer, R<:Real}
    res  = true
    for i in 1:3
        res &= (isOverlapping1D(minima[i,I],maxima[i,I],minima[i,J],maxima[i,J])|| isOverlapping1D(minima[i,I],maxima[i,I],minima[i,J]+BoxLength[i],maxima[i,J]+BoxLength[i]) || isOverlapping1D(minima[i,I],maxima[i,I],minima[i,J]-BoxLength[i],maxima[i,J]-BoxLength[i]))
    end 
    return res 
end

@doc raw"""
    computeClustersByBeadDistance(Sim::HPSAnalysis.SimData{T,Int}; Cutoff=5.0) where{T<:Real, Int<:Integer}

Computes clusters based on minimal distances between amino acids of protein and cutoff.

Computes graph network for each frame in EquilibrationTime:RGMeasureStep:NSteps by naivly comparing all pairs of proteins, computing the minimal distance between any two amino acids of those proteins and adding an edge in the network if the distance is below the given cutoff. Returns the concatenation of groups of IDs belonging to weekly connected graph components of each frame.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.
- `Cutoff::Float`: Distance criterion for a Cluster.

**Returns**:
- `Cluster::Vector{Vector{Vector{Int}}}`: List of Clusters of the chains for each step.
"""
function computeClustersByBeadDistance(Sim::SimData{T,Int}; Cutoff=5.0) where {Int<:Integer, T<:Real}

    Clusters = Vector{Vector{Vector{Int}}}()

    minima = zeros(T,3, Sim.NChains)
    maxima = zeros(T,3, Sim.NChains)

    Cut_half = Cutoff/2.0f0
    for (i,step) in enumerate(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)
        println("Step $step/$(Sim.NSteps)")
        G = Graphs.SimpleGraph{Int}()
        add_vertices!(G, Sim.NChains)
        #=for I in 1:Sim.NChains
            Range_I=  Sim.ChainStart[I]:Sim.ChainStop[I]

            minima[1,I] = minimum(Sim.x_uw[Range_I, step])-Cut_half
            maxima[1,I] = maximum(Sim.x_uw[Range_I, step])+Cut_half
            minima[2,I] = minimum(Sim.y_uw[Range_I, step])-Cut_half
            maxima[2,I] = maximum(Sim.y_uw[Range_I, step])+Cut_half
            minima[3,I] = minimum(Sim.z_uw[Range_I, step])-Cut_half
            maxima[3,I] = maximum(Sim.z_uw[Range_I, step])+Cut_half

            minima[1,I] -= Sim.BoxLength[1].*round.(Int,  minima[1,I]./Sim.BoxLength[1])
            minima[2,I] -= Sim.BoxLength[2].*round.(Int,  minima[2,I]./Sim.BoxLength[2])
            minima[3,I] -= Sim.BoxLength[3].*round.(Int,  minima[3,I]./Sim.BoxLength[3])

            maxima[1,I] -= Sim.BoxLength[1].*round.(Int,  maxima[1,I]./Sim.BoxLength[1])
            maxima[2,I] -= Sim.BoxLength[2].*round.(Int,  maxima[2,I]./Sim.BoxLength[2])
            maxima[3,I] -= Sim.BoxLength[3].*round.(Int,  maxima[3,I]./Sim.BoxLength[3])
        end
        =#

        for I in 1:Sim.NChains
            Range_I = Sim.ChainStart[I]:Sim.ChainStop[I]

            for J in I+1:Sim.NChains
                for j in  Sim.ChainStart[J]:Sim.ChainStop[J]
                ### check if there is already overlap between the bounding boxes of the proteins
                #if isOverlapping3D(minima, maxima, I,J, Sim.BoxLength)

                    dx = Sim.x_uw[Range_I, step] .- Sim.x_uw[j, step]
                    dy = Sim.y_uw[Range_I, step] .- Sim.y_uw[j, step]
                    dz = Sim.z_uw[Range_I, step] .- Sim.z_uw[j, step]

                    dx -= Sim.BoxLength[1].*round.(Int, dx./Sim.BoxLength[1])
                    dy -= Sim.BoxLength[2].*round.(Int, dy./Sim.BoxLength[2])
                    dz -= Sim.BoxLength[3].*round.(Int, dz./Sim.BoxLength[3])

                    dist_sqr = dx.^2 .+dy.^2 .+dz.^2

                    if any(dist_sqr.<Cutoff^2)
                        add_edge!(G, I,J)
                        continue
                    end
                #end
                end
            end
        end

        fill!(Sim.ChargeChainContactMatrix, 0)
        push!(Clusters,  weakly_connected_components(G))

    end

    Sim.ClusterRange = Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
    Sim.Clusters = Clusters
    return Clusters
end
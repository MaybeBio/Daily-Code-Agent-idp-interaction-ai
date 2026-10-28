@doc raw"""
    Only created to debugging purposes use *computeInterChainContactMatrix* in production code. May be faster for very small systems.
"""
function computeInterChainContactMatrix_Naiv(Sim::SimData{R,I}; bounds::Union{Nothing, Vector{Tuple{R, R}}}=nothing, CutoffDict=Sim.IDToSigmas, rel_cutoff=1.0) where {R<:Real, I<:Integer}
    #Sim.ContactMatrices = [zeros(R, Sim.ChainLength[C], Sim.ChainLength[C]) for C in 1:Sim.NChains]
    Sim.ContactMatrices = Matrix(zeros(R, Sim.ChainLength[1], Sim.ChainLength[1]))

    ### allocate square cutoffs, indexed by id 
    cutoffs = [CutoffDict[id]*rel_cutoff for id in sort(collect(keys(CutoffDict)))]

    UsedChains = 0

    max_cutoff_sqr = maximum(cutoffs)^2

    for step in Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
        println("Step: $step")
        ### detect chains which are fully within bounds
        if isnothing(bounds)
            chains = collect(1:Sim.NChains)
        else
            chains = Vector{I}()
            for (C, (start, stop)) in enumerate(zip(Sim.ChainStart, Sim.ChainStop))
                (x_min, x_max) = extrema(Sim.x[start:stop,step])
                (y_min, y_max) = extrema(Sim.y[start:stop,step])
                (z_min, z_max) = extrema(Sim.z[start:stop,step])

                if  x_min>bounds[1][1] && x_max < bounds[1][2] && 
                    y_min>bounds[2][1] && y_max < bounds[2][2] && 
                    z_min>bounds[3][1] && z_max < bounds[3][2]
                    push!(chains,C)
                end
            end
        end
        UsedChains += length(chains)

        ### iterate naively over all pairs of chains
        for (Cid, C) in enumerate(chains)
            start = Sim.ChainStart[C]
            stop  = Sim.ChainStop[C]
            @inbounds for (i_rel,i) in enumerate(start:stop)
                i_cutoff = cutoffs[Sim.IDs[i]]

                ### periodic distance
                dist_x = abs.(Sim.x[:,step].-Sim.x[i,step]) 
                dist_y = abs.(Sim.y[:,step].-Sim.y[i,step]) 
                dist_z = abs.(Sim.z[:,step].-Sim.z[i,step]) 

                dist_x -= @. floor(I, dist_x * 1.0/Sim.BoxLength[1]+0.5)*Sim.BoxLength[1]
                dist_y -= @. floor(I, dist_y * 1.0/Sim.BoxLength[2]+0.5)*Sim.BoxLength[2]
                dist_z -= @. floor(I, dist_z * 1.0/Sim.BoxLength[3]+0.5)*Sim.BoxLength[3]

                dist_sqr = @. dist_x^2+dist_y^2+dist_z^2

                for j in findall(dist_sqr.<max_cutoff_sqr)
  
                    if  (j< start || j>stop) # C != C2
                        ChainOfj = searchsortedfirst(Sim.ChainStart, j) -1 
                        ChainOfj = ChainOfj != size(Sim.ChainStart,1) && Sim.ChainStart[ChainOfj+1]==j ? ChainOfj+1 : ChainOfj

                        j_rel = j - Sim.ChainStart[ChainOfj] +1

                        if dist_sqr[j]< ((cutoffs[Sim.IDs[j]]+i_cutoff)/2.0)^2
                            Sim.ContactMatrices[i_rel,j_rel] += 1
                        end
                    end
                end
            end
        end
    end
    println("UsedChains $UsedChains")
    Sim.ContactMatrices = (Sim.ContactMatrices.+transpose(Sim.ContactMatrices))./2.0
    UsedChains = isnan(UsedChains) ? R(1) : UsedChains ### zeros cause nans in the following
    Sim.ContactMatrices ./= UsedChains

    return Sim.ContactMatrices
end

@doc raw"""
    computeInterChainContactMatrix(Sim::SimData{R,I})

Computes the inter chain contact histrogram at each step of the *ClusterRange* by computing the distances between the residues of all proteins that are within the *bounds* in the *SlabAxis* and the neighbours the amino acids. The distance cutoffs for the histograms are defined by *rel_cutoff* ⋅ *CutoffDict*[type].

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.
- `bounds::Union{Nothing, Vector{Tuple{R2, R2}}}`: Bound alongs the SlabAxis at which particles are considered to be in the bulk of the dense phase.
- `CutoffDict::Dict{Char,Float}=Sim.IDToSigmas`: Dictionary that defines the type specific factor for the cutoff.
- `rel_cutoff=1.0`: Global scaling factor for the cutoff.


**Create**:
* Sim.ContactMatrices
"""
function computeInterChainContactMatrix(Sim::SimData{R,I}; bounds::Union{Nothing, Vector{Tuple{R, R}}}=nothing, CutoffDict=Sim.IDToSigmas, rel_cutoff=1.0, NSubMatrices=5,fac=5)::Matrix{R}  where {R<:Real, I<:Integer}
    ContactMatrices = [Matrix(zeros(R, Sim.ChainLength[1], Sim.ChainLength[1])) for _ in 1:NSubMatrices]
    UsedChains = zeros(NSubMatrices)

    chain_of_atom = zeros(I, Sim.NAtoms)
    for C in 1:Sim.NChains
        for atom in Sim.ChainStart[C]:Sim.ChainStop[C]
            chain_of_atom[atom] = C
        end
    end


    ### allocate square cutoffs, indexed by id 
    cutoffs = [CutoffDict[id]*rel_cutoff for id in sort(collect(keys(CutoffDict)))]

    Sim.CellResolution = ones(3).*maximum(cutoffs)
    initCellLists(Sim)
    ### Elastic Network models screw up the estimation for the MaxParticlesPerCell
    Sim.CenterBox = zeros(R, (fac*3 *Sim.MaxParticlesPerCell,3))
    Sim.NeighBox  = zeros(R, (fac*27*Sim.MaxParticlesPerCell,3))
    Sim.CL_Dist   = zeros(R, (fac*27*Sim.MaxParticlesPerCell^2,3)) ### 3 rows intermediate steps, first row contains the results

    ChainAtomsNeighbour = zeros(I, fac*27*Sim.MaxParticlesPerCell)
    ChainsNeighbour     = zeros(I, fac*27*Sim.MaxParticlesPerCell)

    ChainAtomsInCenter  = zeros(I, fac*5*27*Sim.MaxParticlesPerCell)

    COM=zeros(R, 3)
    AxisCOM = computeCOMOfLargestCluster(Sim)

    max_id = maximum(values(Sim.IDs))  # Assuming IDs are integers
    cutoff_combinations = zeros(R, max_id, max_id)
    for i in 1:max_id
        for j in 1:max_id
            cutoff_combinations[i,j] = ((cutoffs[i] + cutoffs[j])/2)^2
        end
    end

    invBoxLength= R.(1.0./Sim.BoxLength)
    for (k,step) in enumerate(Sim.ClusterRange)### ≈ startstep:stepwidth:NSteps
        if step < Sim.EquilibrationTime || step %Sim.RGMeasureStep != 0
            continue 
        end ### ensures that there exist a clustering and while downsampling the evaluated frames
        #println("Step $k of $(length(Sim.ClusterRange))")

        COM[Sim.SlabAxis] = AxisCOM[k]
        chains = isnothing(bounds) ? collect(1:Sim.NChains) : getChainsInBounds(Sim, bounds,COM, step, Sim.BoxLength, invBoxLength)

        #UsedChains += length(chains)

        ### compute CellLists for step
        Sim.CellStep[1] = step
        resetCellLists(Sim)
        computeCellLists(Sim; ComputeCharges=false, ComputeChains=true) 

        ### iterate naively over all pairs of chains
        for (Cid, C) in enumerate(chains)
            C_subMatId = mod1(Cid, NSubMatrices)
            UsedChains[C_subMatId] += 1
            start = Sim.ChainStart[C]
            stop  = Sim.ChainStop[C]

            ### iterate over all cell lists where chain C is present
            for (xind, yind, zind) in getUniqueCellsOfChain(Sim, C)
                IMax = 0
                @inbounds for id in Sim.CellList[xind, yind, zind]
                    if start <= id <= stop
                        IMax += 1
                        ChainAtomsInCenter[IMax] = id
                        Sim.CenterBox[IMax,1] = Sim.x[id, step]
                        Sim.CenterBox[IMax,2] = Sim.y[id, step]
                        Sim.CenterBox[IMax,3] = Sim.z[id, step]
                    end
                end
                #ChainAtomsInCenter = [id for id in Sim.CellList[xind, yind, zind] if id>=start && id <= stop]
                #IMax = length(ChainAtomsInCenter)
                if IMax ==0 continue end
                ### get positions for center box
                #Sim.CenterBox[1:IMax,1] .= Sim.x[ChainAtomsInCenter[1:IMax], step]
                #Sim.CenterBox[1:IMax,2] .= Sim.y[ChainAtomsInCenter[1:IMax], step]
                #Sim.CenterBox[1:IMax,3] .= Sim.z[ChainAtomsInCenter[1:IMax], step]


                ### have to compute all octants since we do not count all chains; aggregating all of them at once is quicker than individual computation -> SIMD 
                JMax = 0
                for zi in -1:1
                    zi_res = zind+zi
                    for yi in -1:1
                        yi_res = yind+yi
                        for xi in -1:1
                            xi_res = xind+xi
                            ### filter all particles that dont belong to chain C
                            #=
                            for (i,id) in enumerate(Sim.CellList[xi_res, yi_res, zi_res])
                                if id<start || id > stop
                                    ChainAtomsNeighbour[JMax] = id 
                                    ChainsNeighbour[JMax] = chain_of_atom[id]#Sim.ChainNumCellList[xi_res, yi_res, zi_res][i]
                                    JMax += 1
                                end
                            end
                            =#
                            cell   = Sim.CellList[xi_res, yi_res, zi_res]          # Vector{I}
                            chainv = Sim.ChainNumCellList[xi_res, yi_res, zi_res] # Vector{I} (same length)

                            # ---- fast loop over the cell (no enumerate, no bounds checks) -----
                            @simd for i in eachindex(cell)
                                id = cell[i]               # atom index
                                # “id not belonging to the current chain”
                                if id < start || id > stop
                                    JMax += 1
                                    ChainAtomsNeighbour[JMax] = id
                                    #ChainsNeighbour[JMax]    = chainv[id]   # or Sim.chain_of_atom[id] if you stored it globally
                                    
                                end
                            end
                        end
                    end
                end
                #println("B")

                if JMax ==0 continue end

                ### get positions for comparison box
                Sim.NeighBox[1:JMax,1] .= Sim.x[ChainAtomsNeighbour[1:JMax], step]
                Sim.NeighBox[1:JMax,2] .= Sim.y[ChainAtomsNeighbour[1:JMax], step]
                Sim.NeighBox[1:JMax,3] .= Sim.z[ChainAtomsNeighbour[1:JMax], step]

                ### compute distances in an optimized manor-> saved in Sim.CL_Dist[:,1]
                #CB_tmp = @views Sim.CenterBox[1:IMax,:]
                #NB_tmp = @views Sim.NeighBox[1:JMax,:]
                computeSqrDistancesForCellNeighbour(Sim.BoxLength, Sim.CenterBox, Sim.NeighBox, Sim.CL_Dist,IMax, JMax)

                cnt = 1
                for i in 1:IMax
                    atom_i = ChainAtomsInCenter[i]
                    rel_i = atom_i - start+1
                    cutoff_i = cutoffs[Sim.IDs[atom_i]]
                    id_i = Sim.IDs[ChainAtomsInCenter[i]]

                    #res_i = Sim.IDs[atom_i]
                    for j in 1:JMax
                        #println("j $j")
                        atom_j = ChainAtomsNeighbour[j] ### here we take neighbour cell
                        id_j = Sim.IDs[ChainAtomsNeighbour[j]]
                        #cutoff_combinations[id_i,id_j]  #
                        #println("$j, $id_i, $id_j, $cnt")
                        if Sim.CL_Dist[cnt,1]< cutoff_combinations[id_i,id_j]  #((cutoff_i+ cutoffs[Sim.IDs[atom_j]])/2.0)^2
                            #rel_j = atom_j - Sim.ChainStart[ChainsNeighbour[j]] +1
                            rel_j = atom_j - Sim.ChainStart[chain_of_atom[atom_j]] +1
                            ContactMatrices[C_subMatId][rel_i, rel_j] += 1
                        end
                        cnt += 1
                    end
                end
            end
        end
    end

    mean = Matrix(zeros(R, Sim.ChainLength[1], Sim.ChainLength[1]))
    error = Matrix(zeros(R, Sim.ChainLength[1], Sim.ChainLength[1]))

    for i in 1:NSubMatrices
        println("UsedChains $(UsedChains[i])")
        ContactMatrices[i] = (ContactMatrices[i].+transpose(ContactMatrices[i]))./2.0
        ContactMatrices[i] ./= UsedChains[i]
        mean .+= ContactMatrices[i] 
    end
    mean ./= NSubMatrices

    for i in 1:NSubMatrices
        error .+= (mean .- ContactMatrices[i]).^2
    end
    error .= sqrt.(error)./sqrt(NSubMatrices)

    Sim.ContactMatrices = mean
    Sim.ContactMatricesError = error
end

function computeInterChainMatrixHelper(Sim::SimData{R,I}, xind::I2, yind::I2, zind::I2, border::Bool, cutoffs::Vector{R2}, SqrCutoffMatrix::Matrix{R}, ValidChain::Vector{Bool}, maxcutff::R2,chain_of_atom::Vector{I}, ChainAtomsInCenter::Vector{I},ChainsInCenter::Vector{I}, ChainAtomsNeighbour::Vector{I}) where {R<:Real,R2<:Real, I<:Integer,I2<:Integer}

    #ChainAtomsInCenter =  [id for (C,id) in zip(Sim.ChainNumCellList[xind, yind, zind], Sim.CellList[xind,yind,zind]) if ValidChain[C]]
    #ChainsInCenter     =  [C  for (C,id) in zip(Sim.ChainNumCellList[xind, yind, zind], Sim.CellList[xind,yind,zind]) if ValidChain[C]]
    step = Sim.CellStep[1]


    cell_chains = Sim.ChainNumCellList[xind, yind, zind]
    cell_atoms  = Sim.CellList[xind,yind,zind]
    IMax=0
    @inbounds @simd for i in eachindex(cell_atoms)   # avoids bounds checks & enables SIMD
        c = cell_chains[i]               # chain number of the i‑th atom in the cell
        if ValidChain[c]                # branch is SIMD‑friendly (mask)
            IMax += 1
            id = cell_atoms[i]
            ChainAtomsInCenter[IMax] = id    # atom id
            ChainsInCenter[IMax] = c         # chain id (you could also use `Sim.chain_of_atom[atom]`)
            Sim.CenterBox[IMax,1] = Sim.x[id, step]
            Sim.CenterBox[IMax,2] = Sim.y[id, step]
            Sim.CenterBox[IMax,3] = Sim.z[id, step]
        end
    end

    #IMax = length(ChainAtomsInCenter)
    if IMax ==0
        return
    end

    ### get positions for center box
    #Sim.CenterBox[1:IMax,1] .= Sim.x[ChainAtomsInCenter[1:IMax], Sim.CellStep[1]]
    #Sim.CenterBox[1:IMax,2] .= Sim.y[ChainAtomsInCenter[1:IMax], Sim.CellStep[1]]
    #Sim.CenterBox[1:IMax,3] .= Sim.z[ChainAtomsInCenter[1:IMax], Sim.CellStep[1]]

    JMax =0
    atom_i = 0
    atom_j = 0
    res_i=0
    #cnt = 1

    ### compute all octants since not every chain is taken
    for xi in -1:1
        xi_res = xind+xi
        for yi in -1:1
            yi_res = yind+yi
            for zi in -1:1
                zi_res = zind+zi
                cnt = 1

                #ChainAtomsNeighbour = Sim.CellList[xi_res,yi_res,zi_res]  #[id for (i,(C,id)) in enumerate(zip(Sim.ChainNumCellList[xi_res, yi_res, zi_res], Sim.CellList[xi_res,yi_res,zi_res])) ]
                #JMax = length(ChainAtomsNeighbour)

                #if JMax == 0 continue end
                ### get positions for center box
                #Sim.NeighBox[1:JMax,1] .= Sim.x[ChainAtomsNeighbour, step]
                #Sim.NeighBox[1:JMax,2] .= Sim.y[ChainAtomsNeighbour, step]
                #Sim.NeighBox[1:JMax,3] .= Sim.z[ChainAtomsNeighbour, step]

                cell   = Sim.CellList[xi_res, yi_res, zi_res]          # Vector{I}
                chainv = Sim.ChainNumCellList[xi_res, yi_res, zi_res] # Vector{I} (same length)

                # ---- fast loop over the cell (no enumerate, no bounds checks) -----
                @simd for i in eachindex(cell)
                    id = cell[i]               # atom index

                        JMax += 1
                        ChainAtomsNeighbour[JMax] = id
                end
            end
        end
    end
    Sim.NeighBox[1:JMax,1] .= Sim.x[ChainAtomsNeighbour[1:JMax], step]
    Sim.NeighBox[1:JMax,2] .= Sim.y[ChainAtomsNeighbour[1:JMax], step]
    Sim.NeighBox[1:JMax,3] .= Sim.z[ChainAtomsNeighbour[1:JMax], step]

    ### compute distances in an optimized manor
    computeSqrDistancesForCellNeighbour(Sim.BoxLength, Sim.CenterBox, Sim.NeighBox, Sim.CL_Dist,IMax, JMax)

    cnt=0
    for i in 1:IMax
        @inbounds @simd  for j in 1:JMax
            cnt += 1
            if Sim.CL_Dist[cnt,1] <  maxcutff
                atom_j = ChainAtomsNeighbour[j] ### here we take neighbour cell
                C_j= chain_of_atom[atom_j]
                #C_j = Sim.ChainNumCellList[xi_res, yi_res, zi_res][j]
                if ChainsInCenter[i] != C_j
                    atom_i = ChainAtomsInCenter[i]
                    id_i= Sim.IDs[atom_i]
                    rel_i = atom_i - Sim.ChainStart[ChainsInCenter[i]]+1
                    cutoff_i = cutoffs[Sim.IDs[atom_i]]
                    id_j= Sim.IDs[atom_j]
                    if Sim.CL_Dist[cnt,1]< SqrCutoffMatrix[id_j,id_i] #(cutoff_i+ cutoffs[Sim.IDs[atom_j]])/2.0
                        rel_j = atom_j - Sim.ChainStart[C_j]+1

                        Sim.ContactMatrices[rel_i, rel_j] += 1
                    end
                end
            end
        end
    end

    return nothing
end

@doc raw"""
    Expected to be faster than *computeInterChainContactMatrix* but so far still slower for the tested system sizes.
"""
function computeInterChainContactMatrix_fast(Sim::SimData{R,I}; bounds::Union{Nothing, Vector{Tuple{R, R}}}=nothing, CutoffDict=Sim.IDToSigmas, rel_cutoff=1.0::R2)::Matrix{R} where {R<:Real, I<:Integer} 
    #Sim.ContactMatrices = [zeros(R, Sim.ChainLength[C], Sim.ChainLength[C]) for C in 1:Sim.NChains]
    Sim.ContactMatrices = Matrix(zeros(R, Sim.ChainLength[1], Sim.ChainLength[1]))

    ### allocate square cutoffs, indexed by id 
    cutoffs = [CutoffDict[id]*rel_cutoff for id in sort(collect(keys(CutoffDict)))]

    chain_of_atom = zeros(I, Sim.NAtoms)
    for C in 1:Sim.NChains
        for atom in Sim.ChainStart[C]:Sim.ChainStop[C]
            chain_of_atom[atom] = C
        end
    end

    UsedChains = 0
    Sim.CellResolution = ones(R,3).*maximum(cutoffs)
    initCellLists(Sim)
    Sim.CenterBox = zeros(R, (Sim.MaxParticlesPerCell,3))
    Sim.NeighBox  = zeros(R, (27*Sim.MaxParticlesPerCell,3))
    Sim.CL_Dist   = zeros(R, (27*Sim.MaxParticlesPerCell^2,3)) ### 3 rows intermediate steps, first row contains the results

    x_low  = ceil(I,bounds[1][1]/Sim.CellResolution[1])
    x_high = ceil(I,bounds[1][2]/Sim.CellResolution[1]) 
    y_low  = ceil(I,bounds[2][1]/Sim.CellResolution[2])
    y_high = ceil(I,bounds[2][2]/Sim.CellResolution[2]) 
    z_low  = ceil(I,bounds[3][1]/Sim.CellResolution[3])
    z_high = ceil(I,bounds[3][2]/Sim.CellResolution[3]) 

    ChainAtomsInCenter  = Vector{I}(zeros(I,Sim.MaxParticlesPerCell))
    ChainsInCenter      = Vector{I}(zeros(I,Sim.MaxParticlesPerCell))
    ChainAtomsNeighbour = Vector{I}(zeros(I,27*Sim.MaxParticlesPerCell))

    chains = Vector{Bool}([true for _ in 1:Sim.NChains])

    COM=zeros(R, 3)
    AxisCOM = computeCOMOfLargestCluster(Sim)

    #ChainAtomsInCenter  = zeros(I, 5*27*Sim.MaxParticlesPerCell)
    #ChainAtomsInCenter  = zeros(I, 5*27*Sim.MaxParticlesPerCell)

    max_id = maximum(values(Sim.IDs))  # Assuming IDs are integers
    cutoff_combinations = zeros(R, max_id, max_id)
    for i in 1:max_id
        for j in 1:max_id
            cutoff_combinations[i,j] = ((cutoffs[i] + cutoffs[j])/2)^2
        end
    end


    for (k,step) in enumerate(Sim.ClusterRange)### ≈ startstep:stepwidth:NSteps
        println("Step $k of $(length(Sim.ClusterRange))")
        COM[Sim.SlabAxis] = AxisCOM[k]
        chains_tmp = getChainsInBounds(Sim, bounds,COM, step, Sim.BoxLength, R.(1.0./Sim.BoxLength))
        chains.=true
        for id in chains_tmp
            chains[id]=true
        end
        UsedChains += sum(chains)

        ### compute CellLists for step
        Sim.CellStep[1] = step
        resetCellLists(Sim)
        computeCellLists(Sim; ComputeCharges=false, ComputeChains=true) 

        ### iterate over all pairs of chains
        computeInterChainMatrixHelper_(Sim::SimData{R,I}, xind::I2, yind::I2, zind::I2, border::Bool) where {R<:Real, I<:Integer,I2<:Integer}  = computeInterChainMatrixHelper(Sim, xind, yind, zind, border, cutoffs,cutoff_combinations, chains, maximum(cutoffs),chain_of_atom, ChainAtomsInCenter, ChainsInCenter,ChainAtomsNeighbour )#, ChainAtomsInCenter, ChainsInCenter)

        iterateThroughCellList(Sim, computeInterChainMatrixHelper_; CellMaxima=[x_low+1, x_high, y_low+1, y_high, z_low+1, z_high])
    end
    Sim.ContactMatrices = (Sim.ContactMatrices+transpose(Sim.ContactMatrices))./2.0
    Sim.ContactMatrices ./= UsedChains

    Sim.ContactMatrices
end

@doc raw"""
    computeIntraChainContactMatrix(Sim::SimData{R,I})

Computes the intra chain contact map by averaging over all pair distances first and taking the logarithm thereof. 

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Create**:
* Sim.IntraChainContactMatrix
"""
function computeIntraChainContactMatrix(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Sim.IntraChainContactMatrix = [zeros(R, Sim.ChainLength[C], Sim.ChainLength[C]) for C in 1:Sim.NChains]

    for step in Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
        for (C, (start, stop)) in enumerate(zip(Sim.ChainStart, Sim.ChainStop))
            @inbounds for (i_rel,i) in enumerate(start:stop)
                for j in i+1:stop
                    j_rel = j-Sim.ChainStart[C]+1
                    dist_sqr = (Sim.x_uw[i,step]-Sim.x_uw[j,step])^2+(Sim.y_uw[i,step]-Sim.y_uw[j,step])^2+(Sim.z_uw[i,step]-Sim.z_uw[j,step])^2
                    dist = sqrt(dist_sqr)
                    Sim.IntraChainContactMatrix[C][i_rel,j_rel] += dist
                end
            end
        end
    end
    
    for (C, (start, stop)) in enumerate(zip(Sim.ChainStart, Sim.ChainStop))
        Sim.IntraChainContactMatrix[C] ./= length( Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)
        @inbounds for (i_rel,i) in enumerate(start:stop)
            for j in i+1:stop 
                j_rel = j-Sim.ChainStart[C]+1
                Sim.IntraChainContactMatrix[C][i_rel,j_rel] = log(Sim.IntraChainContactMatrix[C][i_rel,j_rel])
                Sim.IntraChainContactMatrix[C][j_rel,i_rel] = Sim.IntraChainContactMatrix[C][i_rel,j_rel]
            end
        end
    end

    Sim.IntraChainContactMatrix
end

@inline function yukawa_potential(x::R, q1::R,q2::R, κ::R, prefac::R)::R where {R<:Real}
    return  prefac* q1 * q2 / x * exp(-x * κ)
end

@inline function ashbaugh_hatch_potential(x::R, σ::R,λ::R, ϵ::R=1.0)::R  where {R<:Real}
    lj = 4 * ϵ * ((σ / x)^12 - (σ / x)^6)
    ah = x <  2^(1/6)*σ ? lj + ϵ * (1 - λ) : λ* lj
    return ah
end

@inline function ashbaugh_hatch_potential_fast(r::R, σ::R,λ::R, ϵ::R=1.0)::R  where {R<:Real}
    inv_r = 1.0 / r
    σ_r = σ * inv_r
    σ_r2 = σ_r * σ_r
    σ_r6 = σ_r2 * σ_r2 * σ_r2
    σ_r12 = σ_r6 * σ_r6
    
    lj = 4 * ϵ * (σ_r12 - σ_r6)

    #lj = 4 * ϵ * ((σ / x)^12 - (σ / x)^6)
    ah = r <  2^(1/6)*σ ? lj + ϵ * (1 - λ) : λ* lj
    return ah
end

function relative_distance(dx::Vector{R}, dy::Vector{R}, dz::Vector{R}, Sim::SimData{R,I}, C2::I2, step::I, xi::R, yi::R, zi::R,dist_sqr::Vector{R})::Nothing where {R<:Real, I<:Integer, I2<:Integer} 
    @inbounds @fastmath for (i, idx) in enumerate(Sim.ChainStart[C2] : Sim.ChainStop[C2])
        dx[i] = abs(Sim.x[idx, step] - xi)
        dy[i] = abs(Sim.y[idx, step] - yi)
        dz[i] = abs(Sim.z[idx, step] - zi)
        dist_sqr[i] = 0.0
    end
    nothing
end

function apply_periodic!(dist::Vector{R}, dist_sqr::Vector{R},inv_bl::R, bl::R, half=R(0.5))::Vector{R} where {R<:Real} 
    @inbounds @fastmath  @simd for i in eachindex(dist)
        n = floor( dist[i] * inv_bl + half)
        dist[i] -= n * bl
        dist_sqr[i] += dist[i]^2
    end
    return dist_sqr
end

function read_FF_Interaction_Parameters(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    file = open("$(Sim.BasePath)/HOOMD_Setup/Params.txt", "r")
    yu_pref, κ,rcut_yu, rcut_ah = missing, missing, missing, missing
    while !eof(file)
        line = readline(file)  
        split = Base.split(line)
        if split[1] == "yk_prefactor:"
            yu_pref = parse(R,split[2])*R(10.0) #convert from nm to Å
        end
        if split[1] == "kappa:"
            κ = parse(R,split[2])*R(1/10) #convert from nm^-1 to Å^-1 
        end
        if split[1] == "YukawaCutoff:"
            rcut_yu = parse(R,split[2])*R(10.0) #convert from nm to Å
        end
        if split[1] == "AHCutoff:" 
            rcut_ah = parse(R,split[2])*R(10.0) #convert from nm to Å
        end
    end
    close(file)

    return yu_pref, κ, rcut_yu, rcut_ah
end

function allocateEnergyMaps(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    # energy maps (initially zero)
    numTypes = maximum(collect(keys(Sim.IDToResName)))
    emap_inter_YU  = zeros(R, Sim.ChainLength[1], Sim.ChainLength[1])
    emap_inter_AH  = zeros(R, Sim.ChainLength[1], Sim.ChainLength[1])
    emap_inter_AA_YU = zeros(R,numTypes, numTypes)
    emap_inter_AA_AH = zeros(R,numTypes, numTypes)

    sigmas  = [Sim.IDToSigmas[id]  for id in Sim.IDs]
    lambdas = [Sim.IDToLambdas[id] for id in Sim.IDs]

    dist_x   = zeros(R, Sim.ChainLength[1])
    dist_y   = zeros(R, Sim.ChainLength[1])
    dist_z   = zeros(R, Sim.ChainLength[1])
    dist_sqr           = zeros(R, Sim.ChainLength[1])

    inv_box_length = R.(1.0./Sim.BoxLength)
    inv_x = inv_box_length[1]
    inv_y = inv_box_length[2]
    inv_z = inv_box_length[3]

    bl_x = R(Sim.BoxLength[1])
    bl_y = R(Sim.BoxLength[2])
    bl_z = R(Sim.BoxLength[3])

    return emap_inter_YU,emap_inter_AH,emap_inter_AA_YU,emap_inter_AA_AH,sigmas, lambdas, dist_x,dist_y,dist_z, dist_sqr, inv_x,inv_y,inv_z, bl_x, bl_y,bl_z
end

function symmetrizeEMaps(Sim::SimData{R,I},emap_YU::Matrix{R},emap_AH::Matrix{R},emap_AA_YU::Matrix{R},emap_AA_AH::Matrix{R}, chain_count::I2;diagonal=true) where {R<:Real, I<:Integer, I2<:Integer}
    chain_count = chain_count == 0 ? 1 : chain_count
    #  Average over frames
    emap_YU    ./= R(chain_count)
    emap_AH    ./= R(chain_count)
    emap_AA_YU ./= R(chain_count)
    emap_AA_AH ./= R(chain_count)

    #  Symmetrize (add transpose)
    emap_YU    .+= transpose(emap_YU)
    emap_AH    .+= transpose(emap_AH)
    emap_AA_YU .+= transpose(emap_AA_YU)
    emap_AA_AH .+= transpose(emap_AA_AH)

    if diagonal
        #  Halve the diagonal entries (counted twice)
        for i in 1:Sim.ChainLength[1]
            emap_YU[i, i] /= 2.0
            emap_AH[i, i] /= 2.0
        end
        for i in axes(emap_AA_YU,1)
            emap_AA_YU[i, i] /= 2.0
            emap_AA_AH[i, i] /= 2.0
        end
    end

    #normalize by the number of chain pairs
    emap_AA_YU = normalize_energy_CMs(Sim, emap_AA_YU)
    emap_AA_AH = normalize_energy_CMs(Sim, emap_AA_AH)
    return emap_YU,emap_AH,emap_AA_YU,emap_AA_AH
end

@doc raw"""
    computeMeanPairEnergyMatrix_naiv(Sim::SimData{R,I}, bounds; ϵ_ah)

Computes the mean yukawa and ashbaugh hatch interaction energy seperatly for all chains that are within the bounds at each step in the Sim.ClusterRange. Also returns the mean interaction energy per amino acids for both interactions.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.
- `bounds`::Vector{Tuple{Float32, Float32}}: Vector with the upper and lower bounds accept as bulk for each dimension.
- `ϵ_ah` = 0.8368: interaction prefactor for the ashbaugh hatch potential 

**Create**:
- Sim.MeanYukawaEnergy
- Sim.MeanYukawaEnergy_perAA
- Sim.MeanAshbaughHatchEnergy
- Sim.MeanAshbaughHatchEnergy_perAA

**Returns**:
- Sim.MeanYukawaEnergy
- Sim.MeanYukawaEnergy_perAA
- Sim.MeanAshbaughHatchEnergy
- Sim.MeanAshbaughHatchEnergy_perAA
"""
function computeMeanPairEnergyMatrix_naiv(Sim::SimData{R,I}, bounds::Vector{Tuple{Float32, Float32}}; ϵ_ah = R(0.8368))  where {R<:Real, I<:Integer}
    @assert Sim.SlabAxis == 2 "computeMeanPairEnergyMatrix_naiv could not work if SlabAxis != y"
    ### energies in KJ
    yu_pref, yu_κ, rcut_yu, rcut_ah = read_FF_Interaction_Parameters(Sim)

    rcut_yu_sqr = rcut_yu^2 ### precompute sqr cutoff
    rcut_ah_sqr = rcut_ah^2 ### precompute sqr cutoff

    emap_inter_YU,emap_inter_AH,emap_inter_AA_YU,emap_inter_AA_AH,sigmas, lambdas, dist_x,dist_y,dist_z, dist_sqr, inv_x,inv_y,inv_z, bl_x, bl_y,bl_z = allocateEnergyMaps(Sim::SimData{R,I})
    inv_box_length = R.(1.0./Sim.BoxLength)

    COM=zeros(R, 3)
    AxisCOM = computeCOMOfLargestCluster(Sim)
    
    chain_count = 0
    @inbounds for (k,step) in enumerate(Sim.ClusterRange) ### ≈ startstep:stepwidth:NSteps
        ### check which chains are within bounds at the current step
        COM[Sim.SlabAxis] = AxisCOM[k]
        chains = isnothing(bounds) ? collect(1:Sim.NChains) : getChainsInBounds(Sim, bounds,COM, step, Sim.BoxLength, inv_box_length)
        chain_count+= length(chains)

        for (C_count,C1) in enumerate(chains)
            for (i1_rel,i1) in enumerate(Sim.ChainStart[C1]:Sim.ChainStop[C1])
                ### do the particle look ups once
                i1_id  = Sim.IDs[i1]
                i1_σ   = sigmas[ i1]
                i1_λ   = lambdas[i1]
                q1 = Sim.Charges[i1]
                xi1= Sim.x[i1,step]
                yi1= Sim.y[i1,step]
                zi1= Sim.z[i1,step]

                ### we have to check against all particles since we count all interactions for chains within the bounds
                for C2 in 1:Sim.NChains 
                    if C1==C2 continue end

                    ### computes relative distance dist_x,... and zeros dist_sqr
                    relative_distance(dist_x, dist_y, dist_z, Sim, C2, step, xi1, yi1, zi1,dist_sqr)

                    ### wraps back to minimal distance and sums up in dist_sq;,no allocations
                    dist_sqr = apply_periodic!(dist_x,dist_sqr,inv_x, bl_x)
                    dist_sqr = apply_periodic!(dist_y,dist_sqr,inv_y, bl_y)
                    dist_sqr = apply_periodic!(dist_z,dist_sqr,inv_z, bl_z)
        
                    range = Sim.ChainStart[C2]:Sim.ChainStop[C2]
                    @fastmath @simd for i2_rel in eachindex(dist_sqr)
                        if dist_sqr[i2_rel] < rcut_yu_sqr
                            r_sqr = dist_sqr[i2_rel]
                            r = sqrt(r_sqr)
                            
                            q2 = Sim.Charges[range[i2_rel]]
                            en_yu = yukawa_potential(r, q1, q2, yu_κ, yu_pref)
                            i2_id = Sim.IDs[range[i2_rel]]
                            
                            emap_inter_YU[i2_rel, i1_rel]  += en_yu
                            emap_inter_AA_YU[i2_id, i1_id] += en_yu
                            
                            # Branch-free Ashbaugh-Hatch condition (using squared distance)
                            if r_sqr < rcut_ah_sqr
                                σ = (i1_σ + sigmas[ range[i2_rel]]) * R(0.5)
                                λ = (i1_λ + lambdas[range[i2_rel]]) * R(0.5)

                                en_ah = ashbaugh_hatch_potential_fast(r, σ,λ, ϵ_ah)

                                emap_inter_AH[i2_rel, i1_rel]  += en_ah
                                emap_inter_AA_AH[i2_id, i1_id] += en_ah
                            end
                        end
                    end
                end
            end
        end
    end

    emap_inter_YU,emap_inter_AH,emap_inter_AA_YU,emap_inter_AA_AH = symmetrizeEMaps(Sim,emap_inter_YU,emap_inter_AH,emap_inter_AA_YU,emap_inter_AA_AH, chain_count)

    Sim.MeanYukawaEnergy = emap_inter_YU
    Sim.MeanYukawaEnergy_perAA = emap_inter_AA_YU
    Sim.MeanAshbaughHatchEnergy = emap_inter_AH
    Sim.MeanAshbaughHatchEnergy_perAA = emap_inter_AA_AH

    return emap_inter_YU, emap_inter_AA_YU, emap_inter_AH, emap_inter_AA_AH
end

function count_possible_AApairs_interOnly(Sim::SimData{R,I})::Matrix{R}  where {R<:Real, I<:Integer}
    id = Sim.IDs[Sim.ChainStart[1]:Sim.ChainStop[1]]
    NTypes = maximum(keys(Sim.IDToResName))
    AAcount_per_chain = R.([sum(id.==i) for i in 1:NTypes])

    chainCount = Sim.NChains

    ### dont add chainCount part, since we normalize by chains anyway ; chainCount * (chainCount - 1)
    pair_matrix = [i == j  ?  
        AAcount_per_chain[i] * AAcount_per_chain[i] *  0.5 : 
        AAcount_per_chain[i] * AAcount_per_chain[j] 
        for i in 1:NTypes, j in 1:NTypes
    ]

    return pair_matrix
end

function normalize_energy_CMs(Sim::SimData{R,I}, cont_map_AA::Matrix{R};contact_cutoff::R=R(1.5), normVolume::Bool = false)::Matrix{R} where {R<:Real, I<:Integer}
    pair_matrix = count_possible_AApairs_interOnly(Sim)

    #  Normalise by the number of possible pairs
    cont_map_AA_norm = [ p!=0.0 ? c/p : 0.0 for (p,c) in zip(pair_matrix,cont_map_AA)]

    if normVolume
        volume_matrix = volumeNorm_matrix(Sim, contact_cutoff=contact_cutoff)
        cont_map_AA_norm .= cont_map_AA_norm ./ volume_matrix
    end

    return cont_map_AA_norm
end

function volumeNorm_matrix(Sim::SimData{R,I}; contact_cutoff=1.5) where {R<:Real, I<:Integer}
    NTypes= maximum(key(Sim.IDToSigmas))
    sigmas = [Sim.IDToSigmas[i] for i in 1:NTypes]
    sigma_matrix  = [0.5*(σ1+σ2) for σ1 in sigmas, σ2 in sigmas]

    volume_matrix = R.((4.0/3.0 * π) .* ( (contact_cutoff * 2^(1/6)) .* sigma_matrix).^3)
    return volume_matrix
end

function getHOOMDBody(Sim::SimData{R,I})::Vector{I} where {R<:Real, I<:Integer}
    traj = GSDFormat.open(Sim.TrajectoryFile,"r")[1]
    return Vector{I}(traj.particles.body)
end


@doc raw"""
    getIntraChainEnergyMap(Sim::SimData{R,I}, bounds; ϵ_ah)

Computes the mean yukawa and ashbaugh hatch interaction energy within each chain that is within the bounds at each step in the Sim.ClusterRange. Also returns the mean interaction energy per amino acids for both interactions.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.
- `bounds`::Vector{Tuple{Float32, Float32}}: Vector with the upper and lower bounds accept as bulk for each dimension.
- `ϵ_ah` = 0.8368: interaction prefactor for the ashbaugh hatch potential 

**Create**:
- Sim.MeanIntraYukawaEnergy
- Sim.MeanIntraYukawaEnergy_perAA
- Sim.MeanIntraAshbaughHatchEnergy
- Sim.MeanIntraAshbaughHatchEnergy_perAA

**Returns**:
- Sim.MeanIntraYukawaEnergy
- Sim.MeanIntraYukawaEnergy_perAA
- Sim.MeanIntraAshbaughHatchEnergy
- Sim.MeanIntraAshbaughHatchEnergy_perAA
"""
function getIntraChainEnergyMap(Sim::SimData{R,I}, bounds::Vector{Tuple{Float32, Float32}}; ϵ_ah = R(0.8368)) where {R<:Real, I<:Integer}
    body = getHOOMDBody(Sim)
    mask = body .== -1 ### particles in body -1 interact via LJ and Yukawa potentials
    mask_rel = mask[Sim.ChainStart[1]:Sim.ChainStop[1]]

    ### energies in KJ
    yu_pref, yu_κ, rcut_yu, rcut_ah = read_FF_Interaction_Parameters(Sim)

    rcut_yu_sqr = rcut_yu^2 ### precompute sqr cutoff
    rcut_ah_sqr = rcut_ah^2 ### precompute sqr cutoff

    emap_YU,emap_AH,emap_AA_YU,emap_AA_AH,sigmas, lambdas, dist_x,dist_y,dist_z, dist_sqr, inv_x,inv_y,inv_z, bl_x, bl_y,bl_z = allocateEnergyMaps(Sim)
    inv_box_length = R.(1.0./Sim.BoxLength)

    COM=zeros(R, 3)
    AxisCOM = computeCOMOfLargestCluster(Sim)
    
    chain_count = 0
    @inbounds for (k,step) in enumerate(Sim.ClusterRange)### ≈ startstep:stepwidth:NSteps
        ### check which chains are within bounds at the current step
        COM[Sim.SlabAxis] = AxisCOM[k]
        chains = isnothing(bounds) ? collect(1:Sim.NChains) : getChainsInBounds(Sim, bounds,COM, step, Sim.BoxLength, inv_box_length)
        chain_count+= length(chains)
        for (C_count,C1) in enumerate(chains) ### only take chains which are fully in the bounds
            for (i1_rel,i1) in enumerate(Sim.ChainStart[C1]:Sim.ChainStop[C1])
                if !mask[i1] continue end ### skip amino acids that are part of RRM
                ### do the particle look ups once
                i1_id  = Sim.IDs[i1]
                i1_σ   = sigmas[i1]
                i1_λ   = lambdas[i1]
                q1 = Sim.Charges[i1]
                xi1= Sim.x[i1,step]
                yi1= Sim.y[i1,step]
                zi1= Sim.z[i1,step]

                ### computes relative distance dist_x,... and zeros dist_sqr
                relative_distance(dist_x, dist_y, dist_z, Sim, C1, step, xi1, yi1, zi1,dist_sqr)

                ### wraps back to minimal distance and sums up in dist_sq;,no allocations
                dist_sqr = apply_periodic!(dist_x,dist_sqr,inv_x, bl_x)
                dist_sqr = apply_periodic!(dist_y,dist_sqr,inv_y, bl_y)
                dist_sqr = apply_periodic!(dist_z,dist_sqr,inv_z, bl_z)
    
                range = Sim.ChainStart[C1]:Sim.ChainStop[C1]
                @fastmath @simd for i2_rel in eachindex(dist_sqr)
                    if (dist_sqr[i2_rel] < rcut_yu_sqr && i2_rel>i1_rel+1) && mask_rel[i2_rel]
                        r_sqr = dist_sqr[i2_rel]
                        r = sqrt(r_sqr)
                        
                        q2 = Sim.Charges[range[i2_rel]]
                        en_yu = yukawa_potential(r, q1, q2, yu_κ, yu_pref)
                        i2_id = Sim.IDs[range[i2_rel]]

                        emap_YU[i2_rel, i1_rel]  += en_yu
                        emap_AA_YU[i2_id, i1_id] += en_yu
                        
                        # Branch-free Ashbaugh-Hatch condition (using squared distance)
                        if r_sqr < rcut_ah_sqr
                            σ = (i1_σ + sigmas[ range[i2_rel]]) * R(0.5)
                            λ = (i1_λ + lambdas[range[i2_rel]]) * R(0.5)

                            en_ah = ashbaugh_hatch_potential_fast(r, σ,λ, ϵ_ah)

                            emap_AH[i2_rel, i1_rel]  += en_ah
                            emap_AA_AH[i2_id, i1_id] += en_ah
                        end
                    end
                end
            end
        end
    end

    emap_YU,emap_AH,emap_AA_YU,emap_AA_AH = symmetrizeEMaps(Sim,emap_YU,emap_AH,emap_AA_YU,emap_AA_AH, chain_count;diagonal=false)

    Sim.MeanIntraYukawaEnergy = emap_YU
    Sim.MeanIntraYukawaEnergy_perAA = emap_AA_YU
    Sim.MeanIntraAshbaughHatchEnergy = emap_AH
    Sim.MeanIntraAshbaughHatchEnergy_perAA = emap_AA_AH

    return emap_YU, emap_AA_YU, emap_AH, emap_AA_AH
end

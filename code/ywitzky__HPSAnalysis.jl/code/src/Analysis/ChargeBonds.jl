function computeChargeBondTimesForCell(Sim::HPSAnalysis.SimData{T,I2}, xind::I, yind::I, zind::I, border::Bool) where {I<:Integer, T<:Real, I2<:Integer}
    IMax = length(Sim.PositiveCellList[xind, yind, zind]) 
    if IMax ==0
        return
    end

    ### get positions for center box
    Sim.CenterBox[1:IMax,1] .= Sim.x[Sim.PositiveCellList[xind, yind, zind], Sim.CellStep[1]]
    Sim.CenterBox[1:IMax,2] .= Sim.y[Sim.PositiveCellList[xind, yind, zind], Sim.CellStep[1]]
    Sim.CenterBox[1:IMax,3] .= Sim.z[Sim.PositiveCellList[xind, yind, zind], Sim.CellStep[1]]

    JMax =0
    atom_i = 0
    atom_j = 0
    res_i=0
    #cnt = 1


    ### only compute "upper" octance to avoid double computation
    for xi in xind:xind+1
        for yi in yind:yind+1
            for zi in zind:zind+1
                cnt = 1
                #JMax= Sim.CellListCounters[2, xi, yi, zi]
                JMax= length(Sim.NegativeCellList[xi, yi, zi]) 

                Sim.NeighBox[1:JMax,1] .= Sim.x[Sim.NegativeCellList[xi, yi, zi], Sim.CellStep[1]]
                Sim.NeighBox[1:JMax,2] .= Sim.y[Sim.NegativeCellList[xi, yi, zi], Sim.CellStep[1]]
                Sim.NeighBox[1:JMax,3] .= Sim.z[Sim.NegativeCellList[xi, yi, zi], Sim.CellStep[1]]

                ### computes distances into CL_Dist[:,1]
                computeDistancesForCellNeighbour(Sim.BoxLength, Sim.CenterBox, Sim.NeighBox, Sim.CL_Dist,IMax, JMax)
                
                for i in 1:IMax
                    atom_i = (Sim.PositiveCellList[xind, yind, zind])[i]
                    #res_i = Sim.IDs[atom_i]
                    for j in 1:JMax
                        atom_j = (Sim.NegativeCellList[xi, yi, zi])[j] ### here we take neighbour cell
                        #res_j = Sim.IDs[atom_j]
                        if Sim.CL_Dist[cnt,1]< Sim.ChargeAnalysisCutoff
                            if isNotInSameChain(Sim, atom_i, atom_j)
                                ### compute Charge times
                                key = (atom_i,atom_j)
                                if (haskey(Sim.ChargeContactValid, key))
                                    Sim.ChargeContactTime[key] += 1
                                    Sim.ChargeContactValid[key] = true
                                else
                                    Sim.ChargeContactValid[key] = true
                                    Sim.ChargeContactTime[key] = 1
                                end
                            end
                        end 
                        cnt += 1
                    end
                end
            end
        end
    end
    return nothing
end

function computeChargeBonds(Sim::HPSAnalysis.SimData{T,I2}, xind::I, yind::I, zind::I, border::Bool) where {I<:Integer, T<:Real, I2<:Integer}
    IMax = length(Sim.PositiveCellList[xind, yind, zind]) 
    if IMax ==0
        return
    end

    ### get positions for center box
    Sim.CenterBox[1:IMax,1] .= Sim.x[Sim.PositiveCellList[xind, yind, zind], Sim.CellStep[1]]
    Sim.CenterBox[1:IMax,2] .= Sim.y[Sim.PositiveCellList[xind, yind, zind], Sim.CellStep[1]]
    Sim.CenterBox[1:IMax,3] .= Sim.z[Sim.PositiveCellList[xind, yind, zind], Sim.CellStep[1]]

    JMax =0
    atom_i = 0
    atom_j = 0
    res_i=0
    #cnt = 1


    ### only compute "upper" octance to avoid double computation
    for xi in xind:xind+1
        for yi in yind:yind+1
            for zi in zind:zind+1
                cnt = 1
                #JMax= Sim.CellListCounters[2, xi, yi, zi]
                JMax= length(Sim.NegativeCellList[xi, yi, zi]) 

                Sim.NeighBox[1:JMax,1] .= Sim.x[Sim.NegativeCellList[xi, yi, zi], Sim.CellStep[1]]
                Sim.NeighBox[1:JMax,2] .= Sim.y[Sim.NegativeCellList[xi, yi, zi], Sim.CellStep[1]]
                Sim.NeighBox[1:JMax,3] .= Sim.z[Sim.NegativeCellList[xi, yi, zi], Sim.CellStep[1]]

                ### computes distances into CL_Dist[:,1]
                computeDistancesForCellNeighbour(Sim.BoxLength, Sim.CenterBox, Sim.NeighBox, Sim.CL_Dist,IMax, JMax)
                
                for i in 1:IMax
                    atom_i = (Sim.PositiveCellList[xind, yind, zind])[i]
                    #res_i = Sim.IDs[atom_i]
                    for j in 1:JMax
                        atom_j = (Sim.NegativeCellList[xi, yi, zi])[j] ### here we take neighbour cell
                        #res_j = Sim.IDs[atom_j]
                        if Sim.CL_Dist[cnt,1]<Sim.GofR_MaxRange
                            if isNotInSameChain(Sim, atom_i, atom_j)
                                ind = ceil(I, Sim.CL_Dist[cnt,1]/Sim.GofR_Resolution)
                                Sim.GofR[ind]+=one(T)
                                if Sim.CL_Dist[cnt,1]< Sim.ChargeAnalysisCutoff
                                    key = (atom_i,atom_j)
                                    Sim.ChargeContactValid[key] = true
                                end
                            end 
                        end
                        cnt += 1
                    end
                end
            end
        end
    end
    return nothing
end

function computeChargeContacts(Sim::SimData{T,I}, MaxTimeStep::I2) where {T<:Real, I<:Integer, I2<:Integer}
    id_in_chain_i = 0 
    id_in_chain_j = 0
    chain_num_i = 0
    chain_num_j = 0
    for key in keys(Sim.ChargeContactValid)
        atom_i = key[1] ### always positive charge
        atom_j = key[2] ### always negative charge

        ### compute charge bond Matrix

        ### compute chain numbers
        chain_num_i = ceil(I,atom_i/Sim.MaxChainLength)
        chain_num_j = ceil(I,atom_j/Sim.MaxChainLength)

        ### compute residue number in chain
        id_in_chain_i = atom_i - (chain_num_i-1)*Sim.MaxChainLength
        id_in_chain_j = atom_j - (chain_num_j-1)*Sim.MaxChainLength
        
        ### increase Contactmatrices
        Sim.ChargeResidueContactMatrix[id_in_chain_i, id_in_chain_j] += 1
        Sim.ChargeChainContactMatrix[chain_num_i, chain_num_j] += 1

        ### compute bonds times
        if( ~Sim.ChargeContactValid[key]) ### Bond is broken now, time was updated in that step
            i_id = Sim.PosChargeToHistID[Sim.Charges[atom_i]]
            j_id = Sim.NegChargeToHistID[Sim.Charges[atom_j]]
            
            if (Sim.ChargeContactTime[key]<MaxTimeStep)
                Sim.ChargeContactTimeHist[i_id, j_id, Sim.ChargeContactTime[key]] += 1
            else
                Sim.ChargeContactTimeHist[i_id, j_id, MaxTimeStep] += 1
            end
            ### remove broken bonds to speed up searches
            delete!(Sim.ChargeContactTime, key) 
            delete!(Sim.ChargeContactValid, key) 
        else 
            Sim.ChargeContactValid[key]=false ### contacts that arent set to true are broken in the next step
        end
    end
    return nothing
end

function computeHREMDChargeContacts(Sim::SimData{T,I}, ResidueCM::Array{T,3}, ChainCM::Array{T,3}, ChargeID::I2) where {T<:Real, I<:Integer, I2<:Integer}
    id_in_chain_i = 0 
    id_in_chain_j = 0
    chain_num_i = 0
    chain_num_j = 0
    for key in keys(Sim.ChargeContactValid)
       ### compute bonds times
        atom_i = key[1] ### always positive charge
        atom_j = key[2] ### always negative charge

        ### compute chain numbers
        chain_num_i = ceil(I,atom_i/Sim.MaxChainLength)
        chain_num_j = ceil(I,atom_j/Sim.MaxChainLength)

        ### compute residue number in chain
        id_in_chain_i = atom_i - (chain_num_i-1)*Sim.MaxChainLength
        id_in_chain_j = atom_j - (chain_num_j-1)*Sim.MaxChainLength
        
        ### increase Contactmatrices
        ResidueCM[id_in_chain_i, id_in_chain_j, ChargeID] += one(T)
        ChainCM[chain_num_i, chain_num_j, ChargeID] += one(T)

        Sim.ChargeContactValid[key]=false ### contacts that arent set to true are broken in the next step
    end
    Sim.ChargeContactValid = Dict{Tuple{I, I}, I}()
end

function countRelevantCharges(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    pos_cnt = 1
    neg_cnt = 1
    for charge in Sim.Charges
        if (charge>0)
            if(~haskey(Sim.PosChargeToHistID, charge))
                Sim.PosChargeToHistID[charge] = pos_cnt
                pos_cnt+=1
            end
        elseif (charge<0)
            if(~haskey(Sim.NegChargeToHistID, charge))
                Sim.NegChargeToHistID[charge] = neg_cnt
                neg_cnt += 1
            end
        end
    end
    pos_cnt -= 1
    neg_cnt -= 1
    return pos_cnt, neg_cnt
end

function initChargeDicts(Sim::SimData{T,I}) where {T<:Real,I<:Integer}
    Sim.ChargeContactValid = Dict()
    Sim.ChargeContactTime  = Dict()
    Sim.PosChargeToHistID = Dict()
    Sim.NegChargeToHistID = Dict()

    printstyled("ChargeContactMatrix is only valid if all proteins have the same length.\n", color=:yellow)
end

function computeChargeCorrelations(Sim::SimData{T,I};MaxTimeStep=250) where {T<:Real, I<:Integer}

    initChargeDicts(Sim)

    ### determine number of positive and negative charges and create indices
    (pos_cnt, neg_cnt) = countRelevantCharges(Sim::SimData)

    ### to avoid allocation in wrapper through StructDispatch
    Sim.ChargeResidueContactMatrix = zeros(T, Sim.MaxChainLength,Sim.MaxChainLength)
    Sim.ChargeChainContactMatrix = zeros(T, Sim.NChains,Sim.NChains)
    Sim.ChargeContactTimeHist= zeros(T,(pos_cnt, neg_cnt, MaxTimeStep)) 
    
    initCellLists(Sim)
    for step in 1:Sim.NSteps
        Sim.CellStep[1] = step
        resetCellLists(Sim) ### empty all cells
        if Sim.CellStep[1]%10==0
            println("Step: $(Sim.CellStep) of $(Sim.NSteps)")
        end
        computeCellLists(Sim) ### compute CellLists  for step

        iterateThroughCellList(Sim, computeChargeBondTimesForCell ) ### compute ChargeContacts for all cells with func computeChargeBondTimesForCell

        #iterateThroughCellList(Sim, computeChargeBonds ) ### compute ChargeContacts for all cells with func computeChargeBondTimesForCell

        computeChargeContacts(Sim, MaxTimeStep)
    end


    atom_i = 0
    atom_j = 0
    ### count all the bonds that are still in tact
    for key in keys(Sim.ChargeContactValid)
        atom_i = key[1] ### always positive charge
        atom_j = key[2] ### always negative charge
        i_id = Sim.PosChargeToHistID[Sim.Charges[atom_i]]
        j_id = Sim.NegChargeToHistID[Sim.Charges[atom_j]]
        if (Sim.ChargeContactTime[key]<MaxTimeStep)
            Sim.ChargeContactTimeHist[i_id, j_id, Sim.ChargeContactTime[key]] += 1
        else
            Sim.ChargeContactTimeHist[i_id, j_id, MaxTimeStep] += 1
        end
    end

    ### incorporate symmetry in Contactmatrices
    Sim.ChargeResidueContactMatrix .= symmetriesMatrix(Sim.ChargeResidueContactMatrix)
    Sim.ChargeChainContactMatrix .= symmetriesMatrix(Sim.ChargeChainContactMatrix)

    ### Norm
    Sim.ChargeResidueContactMatrix /= Sim.NSteps
    Sim.ChargeChainContactMatrix /= Sim.NSteps
end

function computeHREMDChargeAnalysis(Sims::Vector{SimData{T,I}}, ID::Matrix{I2};MaxDistance=25, Resolution=0.1) where {T<:Real, I<:Integer, I2<:Integer}
    printstyled("computeHREMDChargeAnalysis is still heavily manipulated"; color=:yellow)
    NSims = length(Sims)
    ### to avoid allocation in wrapper through StructDispatch
    ChargeResidueContactMatrix = zeros(T,  (Sims[1].MaxChainLength, Sims[1].MaxChainLength, NSims))
    ChargeChainContactMatrix = zeros(T,   (Sims[1].NChains, Sims[1].NChains, NSims))

    MaxDistance *=2
    GofRRange = 0:Resolution:MaxDistance
    GofR = zeros(T,(ceil(I,MaxDistance/Resolution), NSims))
    NStep =0

    pos=0
    neg=0
    Range = 10:10:20 
    for (SID,Sim) in enumerate(Sims[3:end]) ### faulty first 2 sims, should change in general case
        println("Sim $SID")
        Sim.CellResolution =MaxDistance
        Sim.GofR_MaxRange=MaxDistance
        Sim.GofR_Resolution=Resolution
        Sim.GofR= zeros(T,ceil(I,MaxDistance/Resolution))

        initChargeDicts(Sim)

        ### determine number of positive and negative charges and create indices
        (pos, neg) = countRelevantCharges(Sim::SimData)
        
        initCellLists(Sim)
        for step in Range
            id= ID[SID]

            Sim.CellStep[1] = step
            resetCellLists(Sim) ### empty all cells
            if Sim.CellStep[1]%10==0
                println("Step: $(Sim.CellStep) of $(Sim.NSteps)")
            end
            computeCellLists(Sim) ### compute CellLists  for step

            #if maximum(length.(Sim.PositiveCellList))> 329 continue  end

            iterateThroughCellList(Sim, computeChargeBonds  ) ### compute ChargeContacts for all cells with func computeChargeBonds

            ### add to g(r) of correct charge
            GofR[:,id] += Sim.GofR[:]
            Sim.GofR .= 0.0

            computeHREMDChargeContacts(Sim, ChargeResidueContactMatrix, ChargeChainContactMatrix, id) 

        end
       # NStep +=length(Range)
    end

    ### normalisation of GofR
    n_pos = sum(Sims[1].Charges.>0.0)
    n_neg = sum(Sims[1].Charges.<0.0)
    ρ = (n_pos*n_neg)/prod(Sims[1].BoxLength)^2

    
    Norm =  ρ*4.0/3.0*π*(  GofRRange[2:end].^3-GofRRange[1:end-1].^3)*length(Range)

    for (i, Sim) in enumerate(Sims)
        ### incorporate symmetry in Contactmatrices, norm
        ChargeResidueContactMatrix[:,:,i] .= symmetriesMatrix(ChargeResidueContactMatrix[:,:,i] ) ./= NStep
        ChargeChainContactMatrix[:,:,i] .= symmetriesMatrix(ChargeChainContactMatrix[:,:,i] ) ./= NStep

        ### assign to simulations in order for saving and later use
        Sim.ChargeResidueContactMatrix = ChargeResidueContactMatrix[:,:,i]
        Sim.ChargeResidueContactMatrix = ChargeResidueContactMatrix[:,:,i]

        GofR[:,i] ./= Norm
        Sim.GofR = GofR[:,i]
        Sim.GofRRange =  (GofRRange[1:end-1].+GofRRange[2:end])./2
    end



    end


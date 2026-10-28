using Base.Iterators, HDF5, LsqFit, LaTeXStrings,Integrals, Statistics 

function getSlabCoordinate(Sim::SimData{R,I};Unwrapped=true) where {R<:Real, I<:Integer}
    if Sim.SlabAxis==1
        SlabCoord = Unwrapped ? Sim.x_uw : Sim.x
    elseif Sim.SlabAxis==2
        SlabCoord = Unwrapped ? Sim.y_uw : Sim.y
    elseif Sim.SlabAxis == 3
        SlabCoord = Unwrapped ? Sim.z_uw : Sim.z
    else 
        ArgumentError("SlabAxis is not properly specified.")
    end
    return SlabCoord
end

@inline function getRecenteredPositions(coord, com , Len, Len_inv)
    pos = (coord-com)                       ### center slab at 0
    pos -= Len*round(Int32, pos*Len_inv)    ### wrap back to central images
end

@inline function getRecenteredPositions(SlabCoord::Array{R}, atom ,j, step, AxisCOM, Len, Len_inv) where {R<:Real}
    getRecenteredPositions(SlabCoord[atom,step],AxisCOM[j], Len, Len_inv )
end

@inline function getChainsInBounds(Sim::SimData{R,I}, bounds::Union{Nothing, Vector{Tuple{R, R}}}, COM::Vector{R}, step, Len::Vector{R}, Len_inv::Vector{R}) where {R<:Real, I<:Integer}
    chains = Vector{I}()
    if !isnothing(bounds)
        for (C, (start, stop)) in enumerate(zip(Sim.ChainStart, Sim.ChainStop))
            ### computing the extrema first works only most of the times...
            (x_min, x_max) = extrema(getRecenteredPositions.(Sim.x_uw[start:stop, step], COM[1], Len[1], Len_inv[1]))
            (y_min, y_max) = extrema(getRecenteredPositions.(Sim.y_uw[start:stop, step], COM[2], Len[2], Len_inv[2]))
            (z_min, z_max) = extrema(getRecenteredPositions.(Sim.z_uw[start:stop, step], COM[3], Len[3], Len_inv[3]))

            if  x_min>bounds[1][1] && x_max < bounds[1][2] && 
                y_min>bounds[2][1] && y_max < bounds[2][2] && 
                z_min>bounds[3][1] && z_max < bounds[3][2] 

                push!(chains,C)
            end
        end
    end
    return chains
end

function getTimeStepsWhereChainsInBounds(Sim::SimData{R,I}, bounds::Vector{Tuple{R, R}}) where {R<:Real, I<:Integer} 

    N = length(Sim.ClusterRange)
    ranges = sizehint!([Vector{I}() for _ in 1:Sim.NChains],N)
    COM=zeros(R, 3)
    AxisCOM = computeCOMOfLargestCluster(Sim)
    invBoxLength= R.(1.0./Sim.BoxLength)

    for (k,step) in enumerate(Sim.ClusterRange)### ≈ startstep:stepwidth:NSteps
        COM[Sim.SlabAxis] = AxisCOM[k]

        chains = getChainsInBounds(Sim, bounds,COM, step, Sim.BoxLength, invBoxLength)

        for c in chains
            push!(ranges[c],step)
        end
    end
    return ranges
end

@doc raw"""
    computeSlabHistogram(Sim::SimData{R,I}; Use_Alpha=false, Use_Types=false) where {R<:Real, I<:Integer}
Computes centered slab histograms along Sim.SlabAxis.

Centered slab histograms are computed for all amino acids, only the positive and only the negatives at all times. **Use\_Alpha** enables the computation of alpha helices and their own slab histogram. Amino acids specific histograms are enabled through **Use\_Types**. 

Results are not return but stored in Sim.SlabHistogramSeries as an Offset array where the first dimension ranges from -boxwidth/Sim.Resolution:Sim.Resolution:boxwidth/Sim.Resolution. Default Sim.Resolution is set to ``1~\AA``. The second index are the steps at which clusters and slab histogram were computed according to Sim.ClusterRange. The third index are the amino acids which have been used the order in which they are mentioned above. If used, amino acid specific histograms have the index 4+Sim.IDs. 

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Creat**:
* `Sim.SlabHistogramSeries`: Stores mass densities across slabs for different atom types.
"""
function computeSlabHistogram(Sim::SimData{R,I}; Ranges::Union{Nothing, Vector{UnitRange{I}}}=nothing, Use_Alpha=false, Use_Types=false) where {R<:Real, I<:Integer}
    if Use_Alpha && sum(Sim.TorsionAngles[:,1])==0 
        computeDihedralAngles(Sim) ### ensure that TorsionAngles have been computed
    end

    AxisCOM = computeCOMOfLargestCluster(Sim)

    SlabCoord =  getSlabCoordinate(Sim)

    Len = deepcopy(Sim.BoxLength[Sim.SlabAxis])
    Len_inv = 1.0/Len

    # all residue types, one for positive charged, one for negative charged, one for alpha helices,  one for each type, and according to range definition
    NHists = 1 + 2 +1 + Sim.NAtomTypes*Use_Types 
    rangeoffset = 1 + 2 +1 + Sim.NAtomTypes*Use_Types 
    if !isnothing(Ranges)
        NHists += length(Ranges)
    end

    ### array with 1:N steps for -boxwidth:boxwitdh in the direction of the slab
    Sim.SlabHistogramSeries = OffsetArray(zeros(R, Int32(ceil((Sim.BoxSize[Sim.SlabAxis,2]-Sim.BoxSize[Sim.SlabAxis,1])/Sim.Resolution)) , Int32(Sim.NSteps), NHists), Int32(ceil(Sim.BoxSize[Sim.SlabAxis,1]/Sim.Resolution))+1:Int32(ceil(Sim.BoxSize[Sim.SlabAxis,2]/Sim.Resolution)) , 1:Sim.NSteps, 1:NHists)

    AllAxis = [1,2,3]
    deleteat!(AllAxis, Sim.SlabAxis)

    ### convert from dalton/AA^3 to kg/L=g/mL  divided by Volume element per bin
    volume =((Sim.BoxSize[AllAxis[1],2]-Sim.BoxSize[AllAxis[1],1])*(Sim.BoxSize[AllAxis[2],2]-Sim.BoxSize[AllAxis[2],1]))
    conversion = 1.66053906660/volume/Sim.Resolution

    if Use_Alpha
        Pseudohelical = zeros(Bool,Sim.NAtoms)
        AlphaHelixProb = zeros(eltype(Sim.x), Sim.NAtoms)
    end

    lowestind = Int32(ceil(Sim.BoxSize[Sim.SlabAxis,1]/Sim.Resolution))+1
    highestind =Int32(ceil(Sim.BoxSize[Sim.SlabAxis,2]/Sim.Resolution))

    for (j,step) in enumerate(Sim.ClusterRange)### ≈ startstep:stepwidth:NSteps
        if j %200 == 0 println("step $j/$(length(Sim.ClusterRange))") end

        if Use_Alpha ### compute helical states
            Pseudohelical .= false
            AlphaHelixProb .= 0
            decidePseudoHelicals(Sim, Pseudohelical, AlphaHelixProb, step)
        end

        for C in 1:Sim.NChains
            for (atom_in_chain, atom) in enumerate(Sim.ChainStart[C]:Sim.ChainStop[C])
                pos = getRecenteredPositions(SlabCoord, atom,j , step, AxisCOM, Len, Len_inv)
                ind = ceil(Int32,((pos-0.5)/Sim.Resolution))    ### get index according to resolution
                if ind == lowestind-1
                    ind = highestind
                end

                Sim.SlabHistogramSeries[ind , step,1]+= Sim.Masses[atom]
                if Sim.Charges[atom] > 0
                    Sim.SlabHistogramSeries[ind , step,2]+= Sim.Masses[atom]
                elseif Sim.Charges[atom] < 0
                    Sim.SlabHistogramSeries[ind , step,3]+= Sim.Masses[atom]
                end

                if Use_Alpha && AlphaHelixProb[atom]>0
                    Sim.SlabHistogramSeries[ind, step, 4] += Sim.Masses[atom]
                end

                if Use_Types
                    Sim.SlabHistogramSeries[ind, step, Sim.IDs[atom]+4]+= Sim.Masses[atom] ### lowest ID is 1
                end

                if !isnothing(Ranges)
                    atom_in_chain .∈ Ranges 
                    for (k, in_range) in enumerate(atom_in_chain .∈ Ranges)
                        if in_range
                            Sim.SlabHistogramSeries[ind, step, rangeoffset+k]+= Sim.Masses[atom]
                        end
                    end
                end
            end
        end
    end
    Sim.SlabHistogramSeries *= conversion 

    nothing
end

@doc raw"""
    computeDensityHistogram(Sim::SimData{R,I}, DivLength=I(10)) where {R<:Real, I<:Integer}
Computes a logarithmic histogram of densities of the subcubes of the simulation box.

The simulation box is divided into *DivLength*^4 subboxes where each dimension is divided into *DivLength* many subsections and the axes according to Sim.SlabAxis is divided into *DivLength*^2 many subsections.

Subboxes that are empty or those that contain the boundary of the simulation box are not considered for the histogram.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.
- `DivLength=I(10))`: Division factor determining the length of each subbox. 

**Creat**:
* `Sim.DensityHist`: Stores subcube densities with logarithmic indices.
"""
function computeDensityHistogram(Sim::HPSAnalysis.SimData{R,I}, DivLength=I(10))  where {R<:Real, I<:Integer}
    if !(Sim.SlabAxis in [1,2,3])
        ArgumentError("SlabAxis is not properly specified.")
    end

    dims = DivLength*ones(I, 3)
    dims[Sim.SlabAxis] *= DivLength### subdivide long axis for finer grid
    offset  = Sim.BoxSize[:,1]
    divider = Sim.BoxLength./R.(dims)
    BoxHist = zeros(R,  (I.(dims)...) ) ### atom , step in x,y,z makes this 
    ind = zeros(I, 3)

    bools = zeros(Bool,3)
    cnt = 0 

    ### convert from dalton/AA^3 to kg/L=g/mL  divided by Volume element per bin
    volume = prod(divider)
    res = 80 ### take 80 minorticks per power of 10
    conversion = R(1.66053906660/volume)
    Sim.DensityHist = zeros(R, 8*res) ### overexpect indices from 10^-6 to 10^1
    xoff = I(dims[3]*dims[2])
    yoff = I(dims[3])
    for step in  Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps
        for atom in 1:Sim.NAtoms
            
            ind[1] = ceil(I,(Sim.x[atom,step]-offset[1])/ divider[1])
            ind[2] = ceil(I,(Sim.y[atom,step]-offset[2])/ divider[2])
            ind[3] = ceil(I,(Sim.z[atom,step]-offset[3])/ divider[3])
            
            bools[1] = ind[1]<=0 || ind[1]>dims[1]
            bools[2] = ind[2]<=0 || ind[2]>dims[2]
            bools[3] = ind[3]<=0 || ind[3]>dims[3]
            
            if bools[1] || bools[2] || bools[3]
                cnt += 1
                continue
            end

            @inbounds BoxHist[ind[1], ind[2], ind[3]] += Sim.Masses[atom]
        end
        BoxHist .*=  conversion
        indices = ceil.(I,log10.(filter(x->x!=0.0 ,BoxHist[:])).*res).+I(6*res)

        for i in indices
            Sim.DensityHist[i] += R(1)
        end
        fill!(BoxHist, 0.0)
    end
    Sim.DensityHist ./= length(Sim.EquilibrationTime:Sim.RGMeasureStep:Sim.NSteps)*prod(dims)
end

function computeSlabDensities_Help!(SlabHistogramSeries::OffsetArrays.OffsetArray{R, 1, Array{R, 1}}, density::Array{R}) where {R<:Real}
    fill!(density,zero(R))
    N = extrema(axes(SlabHistogramSeries,1))[2]-1

    for i in 0:N 
        density[i+1] += SlabHistogramSeries[i]
        density[i+1] += SlabHistogramSeries[-i]
    end
    density ./= 2.0
    return density
end

function tanh_profile(z, params)
    ρ_l, ρ_v, W, z0 = params
    Δρ = ρ_l - ρ_v
    return (ρ_l + ρ_v)/2 .- (Δρ/2) .* tanh.((z .- z0) ./ W)
end

function determineSurfaceWidth(Sim::HPSAnalysis.SimData{R,I}; Windowlength=300, size=(400,300)) where{R<:Real, I<:Integer} 
    N = extrema(axes(Sim.SlabHistogramSeries,1))[2]-1
    density = zeros(R,N+1)
    NMeasurements= sum(Sim.SlabHistogramSeries[1,Sim.NSteps-Windowlength+1:Sim.NSteps,1].!=0.0)

    xaxis = axes(Sim.SlabHistogramSeries)[1]
    AvgHist = (sum(Sim.SlabHistogramSeries[xaxis,Sim.NSteps-Windowlength:Sim.NSteps,:], dims=2)./(NMeasurements))[:,1,1]

    HPSAnalysis.computeSlabDensities_Help!(AvgHist,density);

    xvals = collect(xaxis[1:1000])
    yvals = collect(density[1:1000])

    W0 = xvals[findfirst(x->x<0.2,yvals)] ### rough guess where the center of the tanh is
    model = curve_fit(tanh_profile, xvals, yvals, [0.5, 0.0,30.0, W0])

    plotSurfaceWidthFit(xvals, yvals, model, Sim.PlotPath*"SurfaceWidthFit"; size=size)

    return model.param[3]
end



@doc raw"""
    computeSlabHistogram(Sim::SimData; Windowlength=200, DenseWidth=25, DiluteWidth=500, MaxVal=0.9, Surface_fac=0.8, safety=75)

Computes average density within dense phase and dilute phase as well as the indices below/above which Sim.SlabHistogramSeries is in dense/dilute phase.

A mirror symmetric density around zero is computed from which dense phase approximation ρ\_app is defined as the mean value of the first **Width** steps. The dense phase boundary is **Surface_fac** times the distance **r_dense** at which the density drops below **MaxVal** times ρ\_app. Similarly the dilute phase boundary is the distance at which the density drops below **1-MaxVal** times ρ\_app plus **r_dense** times (1-**Surface_fac**).
"""
function computeSlabDensities(Sim::HPSAnalysis.SimData{R,I}; Windowlength=200, DenseWidth=25, DiluteWidth=500, MaxVal=0.9, Surface_fac=0.8, safety=75)  where {R<:Real, I<:Integer}

    safety_ = ceil(Int32, safety/Sim.Resolution)

    ### compute mirror symmetric density of centered histogram
    N = extrema(axes(Sim.SlabHistogramSeries,1))[2]-1
    density = zeros(R,N+1)
    NSteps = size(Sim.SlabHistogramSeries,2)

    ρ_dense  = zeros(R,NSteps)
    ρ_dilute = zeros(R, NSteps)

    ### newer iterations dont measure at every frame. Detect which frames actually contain data
    NMeasurements= sum(Sim.SlabHistogramSeries[1,Sim.NSteps-Windowlength+1:Sim.NSteps,1].!=0.0)
    xaxis = axes(Sim.SlabHistogramSeries)[1]
    AvgHist = (sum(Sim.SlabHistogramSeries[xaxis,Sim.NSteps-Windowlength:Sim.NSteps,:], dims=2)./(NMeasurements))[:,1,1]
    
    computeSlabDensities_Help!(AvgHist, density)

    ### take average at interval that is  certainly the center
    avg_dense = mean(density[1:DenseWidth])

    tmp = findlast(density.>avg_dense*MaxVal)
    ind_dense = round(I, tmp*Surface_fac)
    ρ_dense_tmp = mean(density[1:ind_dense]) 

    DiluteWidth = DiluteWidth >= NSteps ? NSteps-1 : DiluteWidth ### avoid negative start index <= 0
    avg_dilute = mean(density[NSteps-DiluteWidth:end])
    ind_dilute = ceil(I,findfirst(density.-avg_dilute .<(1.0-MaxVal)*(ρ_dense_tmp-avg_dilute)))+safety

    for step in axes(Sim.SlabHistogramSeries,2)
        d = @view Sim.SlabHistogramSeries[:, step, 1]

        computeSlabDensities_Help!(Sim.SlabHistogramSeries[:, step, 1], density)
        ρ_dense[step] = mean(density[1:ind_dense]) 
        ρ_dilute[step] = mean(density[ind_dilute:end])
    end

    return ρ_dense, ρ_dilute, ind_dense, ind_dilute
end

@inline function getVoxelIndex(pos, res, off)
    tmp = ceil(Int32,((pos)/res))+off
    #tmp > res ? tmp-2*off :  (tmp< 1 ? tmp+2*off : tmp)
end


function computeBinderCumulantsOfSlabDensities(Sim::HPSAnalysis.SimData{R,I}, indices_dilute, indices_dense)  where {R<:Real, I<:Integer}

    #DenseCumulantBoxes, DiluteCumulantBoxes = computeBinderCumulantsSubBoxes(Sim, indices_dilute, indices_dense)

    order_param = Sim.SlabDenseCumulantBoxes .- Sim.SlabDiluteCumulantBoxes

    Sim.BinderSlabCumulants = zeros(5)
    norm = zeros(5)

    for (wid,width) in enumerate([1,2,3,4,6])
        for ax1 in partition(1:12, width)
            for ax2 in partition(1:12, width)
                tmp = order_param[ax1,ax2,:].^2
                Sim.BinderSlabCumulants[wid] += mean(tmp)^2/mean(tmp.^2)
                norm[wid] += 1
            end
        end
    end
    Sim.BinderSlabCumulants ./= norm
end

function computeBinderCumulantsSubBoxes(Sim::HPSAnalysis.SimData{R,I}, indices_dilute, indices_dense)  where {R<:Real, I<:Integer}
    AxisCOM = computeCOMOfLargestCluster(Sim)

    AllAxis = [1,2,3]
    deleteat!(AllAxis, Sim.SlabAxis)

    ### convert from dalton/AA^3 to kg/L=g/mL  divided by Volume element per bin
    d1 = Sim.BoxLength[AllAxis[1]]/12.0
    d2 = Sim.BoxLength[AllAxis[2]]/12.0
    volume =d1*d2
    conversion = 1.66053906660/volume/Sim.Resolution

    ### get Long Axis and short axis depending on selection
    SlabCoord =  getSlabCoordinate(Sim;Unwrapped=false )
    tmp = deepcopy(Sim.SlabAxis)
    Sim.SlabAxis= AllAxis[1]
    Axis1 =  getSlabCoordinate(Sim;Unwrapped=false)
    Sim.SlabAxis= AllAxis[2]
    Axis2 =  getSlabCoordinate(Sim;Unwrapped=false)
    Sim.SlabAxis= tmp

    Len = deepcopy(Sim.BoxLength[Sim.SlabAxis])
    Len_inv = 1.0/Len

    ### divide short axes into 12 sub boxes
    DenseCumulantBoxes  = zeros(R, 12,12, Sim.NSteps)
    DiluteCumulantBoxes = zeros(R, 12,12, Sim.NSteps)

    dilute_cutoff = indices_dilute .* Sim.Resolution
    dense_cutoff  = indices_dense  .* Sim.Resolution
    for (j,step) in enumerate(Sim.ClusterRange)### ≈ startstep:stepwidth:NSteps
        if j %200 == 0 println("step $j/$(length(Sim.ClusterRange))") end

        for atom in 1:Sim.NAtoms
            ### recenter histogram according to Cluster COMs
            pos = getRecenteredPositions(SlabCoord, atom ,j, step, AxisCOM, Len, Len_inv)

            if abs(pos)< dense_cutoff
                ind1 = getVoxelIndex(Axis1[atom,step], d1, 6)
                ind2 = getVoxelIndex(Axis2[atom,step], d2, 6)
                if ind1 ∈ 1:12 && ind2 ∈ 1:12 ### catch exceptions where position is exactly a boundary
                    DenseCumulantBoxes[ind1, ind2, j] += Sim.Masses[atom]
                end
            end

            if abs(pos)> dilute_cutoff
                ind1 = getVoxelIndex(Axis1[atom,step], d1, 6)
                ind2 = getVoxelIndex(Axis2[atom,step], d2, 6)
                if ind1 ∈ 1:12 && ind2 ∈ 1:12 ### catch exceptions where position is exactly a boundary
                    DiluteCumulantBoxes[ind1, ind2, j] += Sim.Masses[atom]
                end 
            end
        end
        DenseCumulantBoxes[:,:,j] /= dense_cutoff*2 
        DiluteCumulantBoxes[:,:,j] /= (Sim.BoxLength[Sim.SlabAxis]-dilute_cutoff*2 )
    end
    Sim.SlabDenseCumulantBoxes  = DenseCumulantBoxes .*conversion
    Sim.SlabDiluteCumulantBoxes = DiluteCumulantBoxes.*conversion

    return Sim.SlabDenseCumulantBoxes, Sim.SlabDiluteCumulantBoxes
end

@doc raw"""
    computeSurfaceTension(Sim::HPSAnalysis.SimData{R,I},start::I;filename::String = "pressure.h5",tau::I = I(1),NSub::I= I(10))

Computes the **mechanical (pressure‑tensor) surface tension** of a slab
simulation that has been written out by HOOMD‑Blue.  The routine reads the
pressure‑tensor time series from an HDF5 file, extracts the normal and
tangential components with respect to the slab normal (`Sim.SlabAxis`),
applies the standard formula  

\[
\gamma = \frac{L_{\mathrm{slab}}}{2}\,
\bigl\langle P_{norm.} - P_{tang.} \bigr\rangle ,
\]

and returns the mean value together with an estimate of its statistical
uncertainty obtained by block‑averaging.

**Arguments**

- `Sim::HPSAnalysis.SimData{R,I}`:  A simulation data structure containing the Simulation information.
- `start::I`: index (in multiples of `Sim.TrajWriteOutFreq`) of the  first frame to be included in the analysis.

**Keyword arguments**

- `filename::String = "pressure.h5"`: name of the HDF5 file that contains the pressure tensor (`hoomd-data/md/compute/ThermodynamicQuantities/pressure_tensor`).
- `tau::I = I(1)`: stride for sub‑sampling the pressure‑tensor series
  (e.g. `tau=10` uses every 10‑th stored frame).
- `NSub::I = I(10)`: number of sub‑samples used for the block‑averaging
  error estimate.

**Returns**

A tuple `(γ, Δγ)` where  

- `γ` (`Sim.SurfaceTension`) is the mean surface tension,
- `Δγ` (`Sim.SurfaceTensionError`) is the estimated statistical error."""
function computeSurfaceTension(Sim::HPSAnalysis.SimData{R,I},Ranges::Vector{<:AbstractRange{I2}};filename::String="pressure.h5", tau::I=I(1), NSub::I=I(10)) where {R<:Real, I<:Integer,I2<:Integer}
    longfilename = Sim.BasePath*filename 

    if isfile(longfilename)
        file = h5open(longfilename) 
        pressure_tensor = collect(file["hoomd-data/md/compute/ThermodynamicQuantities/pressure_tensor"])
        timestep = collect(file["hoomd-data/Simulation/timestep"])

        pressure_tensor = reshape(pressure_tensor, (6,I(length(pressure_tensor)/6)))

        ### start is assumed to be in multiples of the frequency of the frame outputs
        ### pressure data is written out way more

        #start_time  = start * Sim.TrajWriteOutFreq
        conv_steps(time) = findfirst(x->x>=time* Sim.TrajWriteOutFreq, timestep)
        indices = vcat(collect.([conv_steps(first(range)):tau:conv_steps(last(range)) for range in Ranges])...)

        #start_index = findfirst(x->(x)>=start_time, timestep[:] )
        #https://hoomd-blue.readthedocs.io/en/v6.0.0/hoomd/md/compute/thermodynamicquantities.html#hoomd.md.compute.ThermodynamicQuantities.pressure_tensor
        p_xx = [pressure_tensor[1,i] for i in  indices]
        p_yy = [pressure_tensor[4,i] for i in  indices]
        p_zz = [pressure_tensor[6,i] for i in  indices]

        all_axis = [p_xx, p_yy, p_zz]
        normal_p = all_axis[Sim.SlabAxis]
        tangential_p = sum(deleteat!(all_axis,Sim.SlabAxis))./2.0
        
        n_measure = floor(I, length(normal_p)/NSub)
        if n_measure==0
            Sim.SurfaceTension =R(0)
            Sim.SurfaceTensionError=R(0)
        else
            γ = Sim.BoxLength[Sim.SlabAxis]/2.0 .*(normal_p.-tangential_p)
            Sim.SurfaceTension = mean(R.(γ)) 
            y_mean_sub = [mean(sub) for sub in Iterators.partition(γ, n_measure)][1:end-1] ### last partition is not full...
            Sim.SurfaceTensionError = Statistics.std(y_mean_sub, mean=Sim.SurfaceTension)
            #sqrt(sum((y_mean_sub.-Sim.SurfaceTension).^2)/NSub)
        end
    else
        Sim.SurfaceTension =0.0
        Sim.SurfaceTensionError=0.0
    end

    return Sim.SurfaceTension, Sim.SurfaceTensionError
end


"""
    computeSurfaceTensionCorrection(Sim::HPS.SimData{R,I}, Δρ::R2, W::R2, λ_eff::R2, charge_frac::R2) where {R<:Real, I<:Integer, R2<:Real}

Compute surface tension correction terms (AH and YK) from simulation data.

# Arguments
- `Sim`: Simulation data object (contains force field parameters).
- `Δρ`: Density difference in **kg/L**.
- `W`: Interfacial width in **Å**.
- `λ_eff`: Effective length scale, **unitless**.
- `charge_frac`: Fraction of charge contributing to YK term, **unitless**.

# Returns
- `Δγ_AH`, `Δγ_YK`: Surface tension corrections in units of **10⁻³ kJ/mol/Å²**.

# Notes
Internally calls `computeSurfaceTensionCorrection(...)` with derived parameters from `Sim`.
"""
function computeSurfaceTensionCorrection(Sim::SimData{R,I},Δρ::R2, W::R2, λ_eff::R2,charge_frac::R2,Δf::R2) where {R<:Real,I<:Integer, R2<:Real}
    yu_pref, κ, rcut_yu, rcut_ah = read_FF_Interaction_Parameters(Sim)

    m_monomer = sum(Sim.ChainMasses)./Sim.NAtoms

    Δγ_AH, Δγ_YK, Δγ_ExQ  = computeSurfaceTensionCorrection(R(Δρ), R(W), κ,m_monomer,R(λ_eff), yu_pref, R(charge_frac), R(Δf); r_c_AH=rcut_ah,r_c_DH=rcut_yu)

    return Δγ_AH,Δγ_YK, Δγ_ExQ
end

"""
    computeSurfaceTensionCorrection(Δρ::R, W::R, κ::R, m_monomer::R, λ_eff::R, yk_prefactor::R, charge_frac::R; r_c_AH=R(20.0), r_c_DH=R(40), ϵ_AH=R(0.8368)) where {R<:Real}

Compute surface tension corrections (AH and YK) using integral formulations.

# Arguments
- `Δρ`: Density difference in **kg/L**.
- `W`: Interfacial width in **Å**.
- `κ`: Inverse Debye length in **1/Å**.
- `m_monomer`: Monomer mass in **atomic mass units (Da)**.
- `λ_eff`: Effective length scale, **unitless**.
- `yk_prefactor`: YK prefactor in **kJ/mol·Å**.
- `charge_frac`: Fraction of charge contributing to YK term, **unitless**.
- `r_c_AH`: Cutoff for AH term in **Å** (default: 20.0).
- `r_c_DH`: Cutoff for DH term in **Å** (default: 40.0).
- `ϵ_AH`: AH energy parameter in **kJ/mol** (default: 0.8368).

# Returns
- `Δγ_AH`, `Δγ_YK`: Surface tension corrections in units of **10⁻³ kJ/mol/Å²**.

# Notes
- Integrals are evaluated in reduced coordinates.
- AH term uses Lennard-Jones-like kernel; YK term uses Yukawa-like kernel.
"""
function computeSurfaceTensionCorrection(Δρ::R, W::R, κ::R, m_monomer::R, λ_eff::R, yk_prefactor::R, charge_frac::R, Δf::R; r_c_AH=R(20.0), r_c_DH=R(40), ϵ_AH=R(0.8368)) where {R<:Real}
    ### integrals in reduced coordinates
    Δγ_AH((s,r),(W   )) = (r^-3 -2*r^-9)        *(3*s^3-s)*coth(s*r/W)
    Δγ_DH((s,r),(W,D))  = r^2*(D+r)/D*exp(-r/D) *(3*s^3-s)*coth(s*r/W)

    σ=5.0
    domain = ([0,r_c_AH/σ], [1,Inf]) # (lb, ub)
    prob = IntegralProblem(Δγ_AH, domain,(W/σ))
    integral_AH =solve(prob, HCubatureJL(), reltol = 1e-8, abstol = 1e-8)

    l = 1.0 # Angstroem
    prob = IntegralProblem(Δγ_DH, ([0,r_c_DH/l], [1,Inf]) ,(W/l,1.0/κ/l))
    integral_DH =solve(prob, HCubatureJL(), reltol = 1e-8, abstol = 1e-8)

    l = 1.0 # Angstroem
    prob = IntegralProblem(Δγ_DH, ([0,σ/l], [1,Inf]) ,(W/l,1.0/κ/l))
    integral_ExQ =solve(prob, HCubatureJL(), reltol = 1e-8, abstol = 1e-8)

    conv_to_number_density = 1.0/(m_monomer*1.6605) ### converts to number density per Å^3
    Δn = Δρ * conv_to_number_density

    Δγ_AH  = Δn^2 * 12.0*π* λ_eff*ϵ_AH*σ^4*integral_AH[1]*1000.0
    Δγ_YK  = Δn^2 *charge_frac^2*yk_prefactor*π/2.0*l^3 *integral_DH[1]*1000.0
    Δγ_ExQ = Δn^2 *Δf           *yk_prefactor*π/2.0*l^3 *integral_ExQ[1]*1000.0
    
    return Δγ_AH,Δγ_YK,Δγ_ExQ
end

export plotRGAutocorrComp

function HPSAnalysis.plotBondAngleHistComp(Sims::Vector{HPSAnalysis.SimData{R,I}}, Path::String) where {R <: Real , I<: Integer}
    BondAnglePlot = Plots.plot(xlabel="Bond Angle Ψ [Degree]", ylabel="P(Ψ) [%]", title="Bond angle histogram")
    for Sim in Sims
        Plots.plot!(axes(Sim.BondAngleHist,1)*180.0/size(Sim.BondAngleHist,1), Sim.BondAngleHist[:]*100, label=Sim.SimulationName)
    end
    Plots.savefig(BondAnglePlot, Path*"Comp_BondAngleHist.png")
    Plots.savefig(BondAnglePlot, Path*"Comp_BondAngleHist.pdf")
    return BondAnglePlot
end

function HPSAnalysis.plotRGAutocorrComp(Sims::Vector{HPSAnalysis.SimData{R,I}},  Path::String) where {R<:Real, I<:Integer}
    RG_Auto= Plots.plot(xlabel="lag time [τ]", ylabel="autocorr []", xlims=(0,500),ylims=(-0.5, 1.0))
    for Sim in Sims
        axis = axes(Sim.RGAutocorr[1,:])
        avg = sum(Sim.RGAutocorr, dims=1)[1,:]./Sim.NChains
        err = zeros_like(avg)
        for chain in 1:Sim.NChains
            err .+= (Sim.RGAutocorr[chain,:] .- avg).^2
        end
        err .= sqrt.(err./Sim.NChains)./(Sim.NChains-1)

        Plots.plot!(axis, avg,yerr=err,  label=Sim.SimulationName)
    end

    Plots.savefig(RG_Auto, Path*"_RGAutocorrComp.png")
    Plots.savefig(RG_Auto, Path*"_CompRGAutocorrComp.pdf")
    return RG_Auto
end

function HPSAnalysis.plotDihedralAngleHistComp(Sims::Vector{HPSAnalysis.SimData{R,I}}, Path::String) where {R<:Real, I<:Integer}
    DihedralAnglePlot= Plots.plot(xlabel="Diehedral Angle Ψ [Degree]", ylabel="P(Ψ) [%]", title="Dihedral angle histogram", label="",legend_position=:topleft)
    for Sim in Sims
        DihedralAnglePlot= Plots.plot!(axes(Sim.TorsionHist,2), Sim.TorsionHist[end,:]*100,  label=Sim.SimulationName)
    end

    Plots.savefig(DihedralAnglePlot, Path*"Comp_DihedralAngleHist.png")
    Plots.savefig(DihedralAnglePlot, Path*"Comp_DihedralAngleHist.pdf")
    return DihedralAnglePlot
end

function HPSAnalysis.computeContactMatrixColorMatrix(Sims::Vector{HPSAnalysis.SimData{R,I}},ranges::Vector{Vector{UnitRange{I2}}}=Vector{Vector{UnitRange{Int32}}}(),color_inter=:thermal,color_intra=:hawaii, xspace::Int64=50, amino_acids_volume::Vector{<:Vector{<:Real}}=Vector{Float64}[]) where {R<:Real, I<:Integer,I2<:Integer}
   
   # Prepare the contact matrices (unchanged from the original code)
    inter_CMs = [Sim.ContactMatrices for Sim in Sims]
    addon = ""

    @assert all(map(x->all(0 .<= x .<=1.0 ),  amino_acids_volume)) "'amino_acids_volume' needs to be in range (0,1]."
    if length(amino_acids_volume) == length(Sims) && all(map(x->!(isempty(x)), amino_acids_volume)) 
        addon *= length(ranges) ==  length(Sims) ? "_VolumeWeight" : ""
    else
        amino_acids_volume = [Vector{Float32}([]) for _ in enumerate(inter_CMs)]
    end

    if length(ranges) > 0
        for (i, (CM, Sim)) in enumerate(zip(inter_CMs, Sims))
            _, inter_CMs[i],_,_ = HPSAnalysis.compute_submatrix(CM, Matrix{eltype(CM)}(Sim.ContactMatricesError), ranges[i], amino_acids_volume[i])
            addon = "_with_coarse"
        end
    end

    clims_inter = (0, maximum(maximum.(inter_CMs)))
    intra_CMs   = [sum(Sim.IntraChainContactMatrix) ./ length(Sim.IntraChainContactMatrix) for Sim in Sims]
    clims_intra = (minimum(minimum.(intra_CMs)), maximum(maximum.(intra_CMs)))

    data_inter = collect(range(clims_inter..., 100))
    data_intra = collect(range(clims_intra..., 100))

    # Build colour matrices, limits, ticks for *each* simulation
    Color_Mats = Vector{Matrix{RGB{Float64}}}()
    Lims       = Vector{Tuple{Float64,Float64}}()
    Ticks      = Vector{Vector{Int}}()
    NMats      = Vector{Int}()

    for (intra_CM, inter_CM) in zip(intra_CMs, inter_CMs)
        matrix = triu(inter_CM) + tril(intra_CM, -1)
        is_triu = triu(trues(size(inter_CM)))
        colour_mat = [istriu_ ? get_color(val, clims_inter, color_inter) : get_color(val, clims_intra, color_intra) for (istriu_, val) in zip(is_triu, matrix) ]

        push!(Color_Mats, colour_mat)
        NMat = size(matrix, 1)
        push!(NMats, NMat)
        push!(Lims, (0.5, NMat + 0.5))
        push!(Ticks, vcat([1], collect(0:xspace:NMat)[2:end]))
    end

    return inter_CMs , intra_CMs, clims_inter, clims_intra, data_inter, data_intra, Color_Mats, NMats, Lims, Ticks, addon

end

# blanks out any grid subplots left over when nSims < ncols*nrows
function fillEmptyGridCells!(fig,subplot_range, ncols::Integer, nrows::Integer)
    for i in subplot_range
        row = ((i - 1) ÷ ncols) + 1
        col = ((i - 1) % ncols) + 1
        xlbl = (row == nrows) ? "Amino Acids i " : ""
        ylbl = (col == 1)     ? "Amino Acids j " : ""
        plot!(fig,subplot=i, xlabel = xlbl, ylabel = ylbl, framestyle = :none)
    end
end

# optional red boxes marking Points[i] on top of heat‑map subplot j, for j in subplot_range
function overlayPoints!(fig,Points::Vector{Vector{R2}}, NMats, subplot_range, ::Type{IntType}) where {R2<:Real, IntType<:Integer}
    length(Points) == 0 && return
    for (i, j) in enumerate(subplot_range)
        len = NMats[i]
        if length(Points[i]) > 0
            for point in Points[i]
                plot!(fig,[point,point],[floor(IntType,len*0.975), len], label="", c=:red, subplot=j, linewidth=2.0)
                plot!(fig,[floor(IntType,len*0.975), len],[point,point], label="", c=:red, subplot=j, linewidth=2.0)
            end
        end
    end
end

# one heat‑map subplot per simulation, starting at subplot index offset+1
function plotHeatmapGrid!(fig, Color_Mats, Lims, Ticks, NMats, Labels::Vector{String}, ncols::Integer, nrows::Integer, LabelSize, bspace, lspace; offset::Integer=1, xlbl::AbstractString="Amino Acids i ", ylbl::AbstractString="Amino Acids j ")
    for i in 1:length(Color_Mats)
        subplot_idx = offset + i

        row = ((i - 1) ÷ ncols) + 1
        col = ((i - 1) % ncols) + 1

        xlbl_ = (row == nrows) ? xlbl : ""
        ylbl_ = (col == 1)     ? ylbl : ""
        bs = (row == nrows) ? bspace  : 0mm
        ls = (col ==1) ? lspace : 0mm

        heatmap!(fig, Color_Mats[i],yflip = false,lims = Lims[i],ticks = Ticks[i],subplot = subplot_idx, xlabel = xlbl_, ylabel = ylbl_, AspectRatio = true, yrotation = 90.0,tick_direction = :out,colorbar=false, guidefontsize=LabelSize,ytickfontvalign=:top,bottom_margin=bs,left_margin=ls)

        # label in the upper‑right corner of each panel
        annotate!(fig, [NMats[i] / 2.0],[NMats[i] * 0.95],text(get(Labels, i, ""), :black, :center, LabelSize), subplot = subplot_idx)
    end
end

# box-filter smoothing: each pixel becomes the mean of its (2*radius+1)^2 neighborhood,
# clipped (shrunk, not wrapped/padded) at the matrix edges
function smooth_neighbors(matrix::AbstractMatrix, radius::Integer)
    radius <= 0 && return matrix
    nrows, ncols = size(matrix)
    smoothed = similar(matrix, Float64)
    for j in 1:ncols, i in 1:nrows
        i_lo, i_hi = max(1, i-radius), min(nrows, i+radius)
        j_lo, j_hi = max(1, j-radius), min(ncols, j+radius)
        smoothed[i,j] = mean(@view matrix[i_lo:i_hi, j_lo:j_hi])
    end
    return smoothed
end

function determineGrid(Grid, nSims)
    ncols = isnothing(Grid) ? ceil(Int, sqrt(nSims)) : Grid[1]
    nrows = isnothing(Grid) ? ceil(Int, nSims/ncols) : Grid[2] 
    return Plots.grid(nrows,ncols), ncols,nrows
end

function usePercentage!(inPercent, clims, data; inPermille=false)
    fac,sign = inPermille ? (1000.0,"‰") : (100.0, "\\%")
    label = inPercent ?  latexstring("P(d_{ij}<1.1\\cdotσ\\cdot2^{(1/6)})\\quad ["*sign*"]") :  L"P(d_{ij}<1.1\cdotσ\cdot2^{(1/6)})" 
    clims = inPercent ? clims.*fac : clims
    data .*= inPercent ? fac : 1.0
    return label,clims, data
end

function HPSAnalysis.plotContactMatrixComp(Sims::Vector{HPSAnalysis.SimData{R,I}},Path::String,Labels::Vector{String};Size=(600,600),color_inter=:thermal,color_intra=:hawaii,LabelSize= 10,Points::Vector{Vector{R2}}=Vector{Vector{Int32}}(),xspace::Int64=50,ranges::Vector{Vector{UnitRange{I2}}}=Vector{Vector{UnitRange{Int32}}}(),Grid=nothing,topsepspace=-7mm, bspace=0mm,rightsepspace=-7mm, tspace=0mm,lspace=0mm,inPercent=false, amino_acids_volume::Vector{<:Vector{<:Real}}=Vector{Float64}[]) where {R<:Real, I<:Integer,I2<:Integer, R2<:Real}

    nSims = length(Sims)
    grid,  ncols,nrows = determineGrid(Grid, nSims)  # Determine the grid that will hold the heat‑maps
    comb_layout = @layout [
        a{0.975w}  b{0.025w};
        [grid]    c{0.975h,0.025w}
    ]
    nGrid=ncols*nrows

    inter_CMs , intra_CMs, clims_inter, clims_intra, data_inter, data_intra, Color_Mats, NMats, Lims, Ticks,addon = HPSAnalysis.computeContactMatrixColorMatrix(Sims,ranges,color_inter,color_intra, xspace,amino_acids_volume)

    # Figure and colour‑bars
    fig = plot(layout = comb_layout,margins = -3mm, left_margin = 0mm, right_margin = 0mm,size = Size, top_margin=0mm, bottom_margin=0mm, minorticks=10)

    # intra‑chain colour bar (left)
    heatmap!(data_intra, [1],transpose([data_intra;;]),color = color_intra,colorbar = false,subplot = 1,yaxis = nothing,xmirror = true,xlabel = L"\ln(\langle d_{ij}\rangle/\AA)",ticks = true,tick_direction = :out, labelfontsize = LabelSize, xminorticks=10, bottom_margin=topsepspace, clims=clims_intra,top_margin=tspace)

    # empty placeholder (keeps the left bar thin)
    plot!(subplot = 2, framestyle = :none, bottom_margin=topsepspace, left_margin=rightsepspace,right_margin=2mm)

    ylabel_,clims_inter, data_inter = usePercentage!(inPercent, clims_inter, data_inter)

    heatmap!([1], data_inter, [data_inter;;],colorbar = false, clims = clims_inter,color = color_inter,xticks = nothing,ymirror = true,yguide = nothing,subplot = nGrid + 3,yrotation = 90.0,xmirror = true,ylabel = ylabel_,xlabel = nothing,yticks = true,tick_direction = :out,right_margin = 2.5mm,labelfontsize = LabelSize, left_margin=rightsepspace,yminorticks=10)

    # Plot the heat‑maps - one subplot per simulation
    plotHeatmapGrid!(fig, Color_Mats, Lims, Ticks, NMats, Labels, ncols, nrows, LabelSize, bspace, lspace; offset=2)

    fillEmptyGridCells!(fig, 3+nSims:2+nGrid, ncols, nrows)

    # (Optional) plot point‑sets on top of each heat‑map
    overlayPoints!(fig, Points, NMats, 3:(nSims+2), I)

    Plots.savefig(fig, Path * "Comp_ContactMatrix$(addon).png")
    Plots.savefig(fig, Path * "Comp_ContactMatrix$(addon).pdf")
    return fig
end

function HPSAnalysis.plotContactMatrixCompTopBottomBars(Sims::Vector{HPSAnalysis.SimData{R,I}},Path::String,Labels::Vector{String};Size=(600,600),color_inter=:thermal,color_intra=:hawaii,LabelSize= 10,Points::Vector{Vector{R2}}=Vector{Vector{Int32}}(),xspace::Int64=50,ranges::Vector{Vector{UnitRange{I2}}}=Vector{Vector{UnitRange{Int32}}}(),Grid=nothing,topsepspace=-7mm, bspace=0mm,rightsepspace=0mm, tspace=0mm,lspace=0mm,inPercent=false, amino_acids_volume::Vector{<:Vector{<:Real}}=Vector{Float64}[]) where {R<:Real, I<:Integer,I2<:Integer, R2<:Real}

    nSims = length(Sims)
    grid, ncols,nrows  = determineGrid(Grid, nSims)  # Determine the grid that will hold the heat‑maps
    # colour‑bars stacked above and below the grid instead of beside it
    comb_layout = @layout [
        a{0.05h};
        [grid];
        c{0.05h}
    ]
    nGrid=ncols*nrows

    inter_CMs , intra_CMs, clims_inter, clims_intra, data_inter, data_intra, Color_Mats, NMats, Lims, Ticks,addon = HPSAnalysis.computeContactMatrixColorMatrix(Sims,ranges,color_inter,color_intra, xspace,amino_acids_volume)

    # Figure and colour‑bars
    fig = plot(layout = comb_layout, left_margin = 0mm, right_margin=0.5mm, size = Size, top_margin=0mm, bottom_margin=0mm, minorticks=10)

    # intra‑chain colour bar (top)
    heatmap!(data_intra, [1],transpose([data_intra;;]),color = color_intra,colorbar = false,subplot = 1,yaxis = nothing,xmirror = true,xlabel = L"\ln(\langle d_{ij}\rangle/\AA)",ticks = true,tick_direction = :out, labelfontsize = LabelSize, xminorticks=10, bottom_margin=topsepspace, clims=clims_intra,top_margin=tspace,xlims=clims_intra)

    xlabel_,clims_inter, data_inter = usePercentage!(inPercent, clims_inter, data_inter)

    heatmap!(data_inter,[1], transpose([data_inter;;]),colorbar = false, clims = clims_inter,color = color_inter, xlabel = xlabel_,tick_direction = :out,labelfontsize = LabelSize, left_margin=rightsepspace,xminorticks=10,subplot = nGrid + 2, bottom_margin=3mm,yaxis = nothing, top_margin=topsepspace,xlims=clims_inter) 

    # Plot the heat‑maps - one subplot per simulation
    plotHeatmapGrid!(fig, Color_Mats, Lims, Ticks, NMats, Labels, ncols, nrows, LabelSize, bspace, lspace; offset=1, xlbl="", ylbl="")

    fillEmptyGridCells!(fig,2+nSims:1+nGrid, ncols, nrows)

    # (Optional) plot point‑sets on top of each heat‑map
    overlayPoints!(fig,Points, NMats, 2:(nSims+1), I)

    Plots.savefig(fig, Path * "Comp_ContactMatrixWide$(addon).png")
    Plots.savefig(fig, Path * "Comp_ContactMatrixWide$(addon).pdf")
    return fig
end


function HPSAnalysis.plotContactMatrixDifference(Sims::Vector{HPSAnalysis.SimData{R,I}}, filename::String;color_inter=:thermal,color_intra=:hawaii,ranges::Vector{UnitRange{I2}}=Vector{UnitRange{Int32}}(),xspace::Int64=50,Size=(1050,400), LabelSize=10, bspace=-3mm,lspace=-6mm, rspace=-3mm, tspace=7mm, Labels::Vector{String}=["",""],CheckSequence=true, inPercent=false,CoarseRanges::Vector{Vector{UnitRange{I3}}}=Vector{Vector{UnitRange{Int32}}}()) where {R<:Real, I<:Integer,I2<:Integer,I3<:Integer}
    @assert length(Sims)>=2 "Number of Simulations needs to be atleast 2"

    inter_CMs , intra_CMs, clims_inter, clims_intra, data_inter, data_intra, Color_Mats, NMats, Lims, Ticks,addon = HPSAnalysis.computeContactMatrixColorMatrix(Sims,CoarseRanges,color_inter,color_intra, xspace)
    grid = Plots.grid(1,3)

    comb_layout = @layout [
       [a{0.5w} b{0.5w,0.06h};  
       e{0.5w} f{0.5w,0.9h}];
    ]

    ### SUBPLOT 1

    fig1 = Plots.plot(size=Size,bottom_margin=5mm, layout=comb_layout)

    heatmap!(data_intra, [1],transpose([data_intra;;]),color = color_intra,colorbar = false,subplot = 1,yaxis = nothing,xmirror = true,xlabel = L"\ln(<d_{ij}>)",ticks = true,tick_direction = :out, labelfontsize = LabelSize, xminorticks=10, bottom_margin=bspace, clims=clims_intra, top_margin=tspace)

    label,clims_inter, data_inter = usePercentage!(inPercent, clims_inter, data_inter)

    heatmap!(data_inter, [1],transpose([data_inter;;]),color = color_inter,colorbar = false,subplot = 2,yaxis = nothing,xmirror = true,xlabel =label,ticks = true,tick_direction = :out, labelfontsize = LabelSize, xminorticks=10, bottom_margin=bspace, clims = clims_inter, top_margin=tspace)

    xlbl = "Amino Acids i "
    ylbl = "Amino Acids j "

    Plots.heatmap!(Color_Mats[1], subplot=3,yflip = false,lims = Lims[1],ticks = Ticks[1], xlabel = xlbl, ylabel = ylbl, AspectRatio = true, yrotation = 90.0,tick_direction = :out,colorbar=false, guidefontsize=LabelSize,ytickfontvalign=:top,left_margin=4mm, minorticks=10)

    #Ticks[2]
    Plots.heatmap!(Color_Mats[2], subplot=4, yflip = false,lims = Lims[2],ticks = Ticks[2], xlabel = xlbl, ylabel = "", AspectRatio = true, yrotation = 90.0,tick_direction = :out,colorbar=false, guidefontsize=LabelSize,ytickfontvalign=:top, minorticks=10)

    # label in the upper‑right corner of each panel
    annotate!([NMats[1] / 2.0],[NMats[1] * 0.95],text(get(Labels, 1, ""), :black, :center, LabelSize), subplot = 3)
    annotate!([NMats[2] / 2.0],[NMats[2] * 0.95],text(get(Labels, 2, ""), :black, :center, LabelSize), subplot = 4)

    ### SUBPLOT 2

    comb_layout = @layout [[ c{0.95w} d{0.06h,0.025w};
                             g{0.95w} h{0.9h ,0.025w}];]
    fig2=Plots.plot(size=Size,bottom_margin=5mm, layout=comb_layout)#

    ranges = length(ranges) > 0 ? ranges :  [1:size(inter_CMs[1])[1], 1:size(inter_CMs[2])[1]]
    @assert length(ranges[1]) == length(ranges[2]) "Matrix areas specified with arguemtns \"ranges\" need to be of same size"
    if CheckSequence
        @assert Sims[1].Sequences[1][ranges[1]]==Sims[2].Sequences[1][ranges[2]] "Ranges do not contain the same sequence:\n $(Sims[1].Sequences[1][ranges[1]]) \n $(Sims[2].Sequences[1][ranges[2]])"
    end

    inter_A = inter_CMs[1][ranges[1],ranges[1]]
    inter_B = inter_CMs[2][ranges[2],ranges[2]]

    intra_A = intra_CMs[1][ranges[1],ranges[1]]
    intra_B = intra_CMs[2][ranges[2],ranges[2]]

    intra_diff = intra_A.- intra_B 
    inter_diff = inter_A.- inter_B 

    matrix = triu(inter_diff) + tril(intra_diff, -1)
    is_triu = triu(trues(size(inter_diff)))

    cmax =maximum(maximum(abs.(inter_diff)))
    clims_inter = (-cmax, cmax)
    cmax =maximum(maximum(abs.(intra_diff)))
    clims_intra = (-cmax, cmax)

    color_inter = :redgreensplit #:romaO
    color_intra = :tarn

    colour_mat = [istriu_ ? get_color(-val, clims_inter, color_inter) : get_color(val, clims_intra, color_intra) for (istriu_, val) in zip(is_triu, matrix) ]
    
    ticks = vcat([1], collect(0:xspace:size(inter_diff)[1])[2:end])
    lims= (1,size(inter_diff)[1])

    data_inter = collect(range(clims_inter..., 100))
    data_intra = collect(range(clims_intra..., 100))

    rr_space=-5mm
    sep_space=-3mm
    heatmap!(data_intra, [1],transpose([data_intra;;]),color = color_intra,colorbar = false, subplot = 1,yaxis = nothing,xmirror = true,xlabel = L"\Delta\ln(<d_{ij}>)",ticks = true,tick_direction = :out, labelfontsize = LabelSize, xminorticks=10, bottom_margin=bspace, clims=clims_intra, top_margin=tspace,left_margin=sep_space,right_margin=rr_space)

    Plots.heatmap!(colour_mat, subplot=3, yflip = false, ticks=ticks, xlabel = xlbl, ylabel = "", AspectRatio = true, yrotation = 90.0,tick_direction = :out,colorbar=false, guidefontsize=LabelSize,ytickfontvalign=:top,lims=lims, left_margin=sep_space,right_margin=rr_space, minorticks=10)

    # empty placeholder (keeps the right bar thin)
    plot!(subplot = 2, framestyle = :none, bottom_margin=bspace, left_margin=lspace,right_margin=rspace, top_margin=tspace)


    label,clims_inter, data_inter = usePercentage!(inPercent, clims_inter, data_inter; inPermille=true)

    heatmap!([1], data_inter, [data_inter;;],colorbar = false, clims = clims_inter,color = cgrad(color_inter, rev=true),xticks = nothing,ymirror = true,yguide = nothing,subplot = 4,yrotation = 90.0,xmirror = true,ylabel = label ,xlabel = nothing,tick_direction = :out,labelfontsize = LabelSize, left_margin=lspace,yminorticks=10,right_margin=rspace,)

    #if any(NMats .!= last.(ranges)) ### dont plot squares if equally long proteins are compared.
        for (i,range) in zip([3,4],ranges)
            low  = first(range)
            high =  last(range)
            plot!(fig1,[low, low,high, high,low], [low, high, high, low,low], linewidth=3, linecolor=:red, label="", linestyle=:solid, subplot=i)
        end
    #end

    ### COMBINE PLOTS
    fig3 =vline([1], linewidth=10,label="", c=:black, framestyle=:none,left_margin=-9mm)
    comb_layout= @layout [ j{0.63w}  b{0.025w} c{0.3w,1.0h}]

    fig = Plots.plot(fig1,fig3,fig2,layout=comb_layout)

    Plots.savefig(fig, filename*".png")
    Plots.savefig(fig, filename*".pdf")

    return fig
end

"""
    plotContactMatrixDifferenceGroups(Groups, filename; color_inter=:redgreensplit, color_intra=:tarn, ranges=..., amino_acids_volume=..., asym=false, smooth_inter=0, Labels=String[], xspace=50, Size=(600,600), LabelSize=10, CheckSequence=true, inPercent=false, Grid=nothing, topsepspace=-7mm, bspace=0mm, lspace=0mm)

Plots one differential contact matrix per pair (upper triangle = inter-chain
contact-probability diff, lower triangle = intra-chain distance diff), each
half on its own scale/colorbar, laid out in a grid.

`Groups` holds `2·nPairs` groups, consumed as consecutive pairs
`(Groups[2k-1], Groups[2k])`. Each group is a `Vector{SimData}` of one or more
replicates; replicates go through `computeContactMatrixColorMatrix` (reusing
its coarsening/volume-weighting) and get averaged before the `a - b` diff.

- `ranges[g]`, if given, is a coarsening partition applied to every replicate
  in `Groups[g]` (one per group, not per pair). Empty by default = no coarsening.
  Since matrices aren't size-reduced or cropped, both sides of a pair must
  share the same full sequence (checked via `CheckSequence`).
- `amino_acids_volume[g]`, if given, likewise holds one volume vector per group.
- `asym=true` uses each half's true (min,max) instead of a symmetric ±cmax,
  with the neutral colour still pinned at 0 (`zero_centered_cgrad`).
- `smooth_inter` > 0 box-filters `inter_diffs` with that neighbor radius before
  scaling/coloring.
"""
function HPSAnalysis.plotContactMatrixDifferenceGroups(Groups::Vector{Vector{HPSAnalysis.SimData{R,I}}}, filename::String, Labels::Vector{String}=String[];color_inter=:redgreensplit,color_intra = :tarn,ranges::Vector{Vector{UnitRange{I2}}}=Vector{Vector{UnitRange{Int32}}}(),
amino_acids_volume::Vector{<:Vector{<:Real}}=Vector{Float64}[],
xspace::Int64=50, Size=(600,600), LabelSize=10, CheckSequence=true, inPercent=false, Grid=nothing, asym=false, smooth_inter::Integer=0,tspace=0mm, sepspace=-7mm, bspace=0mm, lspace=0mm) where {R<:Real, I<:Integer, I2<:Integer}
    addon=""
    @assert iseven(length(Groups)) "Groups must come in pairs (a,b), one pair per differential comparison"
    nPairs  = length(Groups) ÷ 2
    nGroups = length(Groups)

    @assert isempty(ranges) || length(ranges) == nGroups "Need either no coarsening ranges, or one coarsening partition per group"
    @assert isempty(amino_acids_volume) || length(amino_acids_volume) == nGroups "Need either no amino_acids_volume, or one volume vector per group"

    if CheckSequence
        for k in 1:nPairs
            @assert Groups[2k-1][1].Sequences[1] == Groups[2k][1].Sequences[1] "Group pair $k does not share the same sequence"
        end
    end

    # each group's replicates go through computeContactMatrixColorMatrix (coarsening + volume-weighting reused from there), then get averaged
    inter_group = Vector{Matrix{Float64}}(undef, nGroups)
    intra_group = Vector{Matrix{Float64}}(undef, nGroups)
    for (g, sims) in enumerate(Groups)
        group_ranges = isempty(ranges) ? Vector{Vector{UnitRange{Int32}}}() : fill(ranges[g], length(sims))
        group_volume = isempty(amino_acids_volume) ? Vector{Float64}[] : fill(amino_acids_volume[g], length(sims))
        inter_CMs, intra_CMs, _, _, _, _, _, _, _, _, _ = HPSAnalysis.computeContactMatrixColorMatrix(sims, group_ranges, color_inter, color_intra, xspace, group_volume)
        inter_group[g] = mean(inter_CMs)
        intra_group[g] = mean(intra_CMs)
    end

    inter_diffs = [inter_group[2k-1] .- inter_group[2k] for k in 1:nPairs]
    intra_diffs = [intra_group[2k-1] .- intra_group[2k] for k in 1:nPairs]
    smooth_inter > 0 && (inter_diffs = [smooth_neighbors(m, smooth_inter) for m in inter_diffs])
    smooth_inter > 0 && (intra_diffs = [smooth_neighbors(m, smooth_inter) for m in intra_diffs])

    # separate scales: inter-chain (contact-probability) and intra-chain (distance) diffs live on different magnitudes
    if asym
        clims_inter = extrema(collect(Iterators.flatten(extrema.(inter_diffs)))*-1.0)
        clims_intra = extrema(collect(Iterators.flatten(extrema.(intra_diffs))))
    else
        cmax_inter = maximum(k -> maximum(abs.(inter_diffs[k])), 1:nPairs)
        cmax_intra = maximum(k -> maximum(abs.(intra_diffs[k])), 1:nPairs)
        clims_inter = (-cmax_inter, cmax_inter)
        clims_intra = (-cmax_intra, cmax_intra)
    end


    Color_Mats = Vector{Matrix{RGBA{Float64}}}(undef, nPairs)
    Lims  = Vector{Tuple{Float64,Float64}}(undef, nPairs)
    Ticks = Vector{Vector{Int}}(undef, nPairs)
    NMats = Vector{Int}(undef, nPairs)
    for k in 1:nPairs
        ### negate inter to match plotContactMatrixDifference's sign convention between contact-probability and distance differences
        combined = triu(-inter_diffs[k]) + tril(intra_diffs[k], -1)
        is_triu = triu(trues(size(combined)))
        Color_Mats[k] = [istriu_ ? get_color(val, clims_inter, color_inter; zero_centered=asym) : get_color(val, clims_intra, color_intra; zero_centered=asym) for (istriu_, val) in zip(is_triu, combined)]
        NMat = size(combined,1)
        NMats[k] = NMat
        Lims[k]  = (0.5, NMat+0.5)
        Ticks[k] = vcat([1], collect(0:xspace:NMat)[2:end])
    end

    grid, ncols, nrows = determineGrid(Grid, nPairs)
    nGrid = ncols*nrows
    comb_layout = @layout [ a{0.05h}; [grid]; c{0.05h} ]

    fig = Plots.plot(layout=comb_layout, size=Size, minorticks=10)

    plotHeatmapGrid!(fig, Color_Mats, Lims, Ticks, NMats, Labels, ncols, nrows, LabelSize, bspace, lspace; offset=1, xlbl="", ylbl="")
    fillEmptyGridCells!(fig, nPairs+2:nGrid+1, ncols, nrows)

    # intra‑chain colour bar (top)
    data_intra = collect(range(clims_intra..., 1000))
    heatmap!(fig, data_intra, [1], transpose([data_intra;;]), color=(asym ? zero_centered_cgrad(clims_intra, color_intra) : color_intra), colorbar=false, subplot=1, yaxis=nothing, xmirror=true, xlabel=L"\Delta\ln(\langle d_{ij}\rangle/\AA)", tick_direction=:out, labelfontsize=LabelSize, xminorticks=10, clims=clims_intra, top_margin=tspace, xlims=clims_intra,bottom_margin=sepspace)

    # inter‑chain colour bar (bottom)
    data_inter = collect(range(clims_inter..., 1000))
    label,clims_inter_scaled, data_inter = usePercentage!(inPercent, clims_inter, data_inter; inPermille=true)
    heatmap!(fig, data_inter, [1], transpose([data_inter;;]), color=(asym ? zero_centered_cgrad(clims_inter_scaled, color_inter) : color_inter), colorbar=false, subplot=nGrid+2, yaxis=nothing, xlabel=label, tick_direction=:out, labelfontsize=LabelSize, xminorticks=10, clims=clims_inter_scaled, bottom_margin=3mm, top_margin=sepspace, xlims=clims_inter_scaled)

    addon *= smooth_inter > 0 ? "_smooth$(smooth_inter)" : ""
    Plots.savefig(fig, filename*"$(addon).png")
    Plots.savefig(fig, filename*"$(addon).pdf")

    return fig
end


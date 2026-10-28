
#

#=
function plotRGAutocorr(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    x = axes(Sim.RGAutocorr,2)
    avg= reduce(+,Sim.RGAutocorr, dims=1)[1,:]./Sim.NChains

    Plots.plot(xlabel="Lag time [steps]", ylabel="AutoCorr [steps]")
    for chain in 1:Sim.NChains
        Plots.plot!(x, Sim.RGAutocorr[chain,:], label="")#, label="$(chain)")
    end
    Plots.vline!([Sim.RGMeasureStep], label="")
    Plots.hline!([0.0], label="")


    Plots.plot!(x, avg, color=:black, label="avg")
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_RG_Autocorr.png")
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_RG_Autocorr.pdf")
end
=#

function HPSAnalysis.plotREEVecAutocorr(Sim::HPSAnalysis.SimData{R,I}, conf=0.1) where {R<:Real, I<:Integer}
    x = (collect(R, axes(Sim.REEVecSeries,1)).-convert(R,1.0)) * convert(R,conf)

    Plots.plot(xlabel="Lag time [ns]", ylabel="AutoCorr [steps]")
    avg = sum(Sim.REEVecSeries[:, :, 1], dims=2)./Sim.NChains
    for chain in 1:Sim.NChains
        Plots.plot!(x, Sim.REEVecSeries[:, chain, 1], label="")
    end
    Plots.plot!(x,avg[:, 1, 1], label="avg.", c=:black)

    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_REE_Vec_Autocorr.png")
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_REE_Vec_Autocorr.pdf")
    Plots.plot!(xlim=(0,3.0))
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_REE_Vec_Autocorr_zoom.png")
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_REE_Vec_Autocorr_zoom.pdf")
end

#=function plotREEAutocorr(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    x = axes(Sim.REEAutocorr,2)
    avg= reduce(+,Sim.REEAutocorr, dims=1)[1,:]./Sim.NChains

    Plots.plot(xlabel="Lag time [steps]", ylabel="AutoCorr [steps]")
    for chain in 1:Sim.NChains
        Plots.plot!(x, Sim.REEAutocorr[chain,:], label="")#, label="$(chain)")
    end
    Plots.vline!([Sim.RGMeasureStep], label="")
    Plots.hline!([0.0], label="")

    Plots.plot!(x, avg, color=:black, label="avg")
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_RG_Autocorr.png")
    Plots.savefig(Sim.PlotPath*Sim.SimulationName*"_RG_Autocorr.pdf")
end=#

function HPSAnalysis.plotChargeBondTime(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Sim.CellResolution=20.
    Sim.ChargeAnalysisCutoff=20.
    if length(Sim.ChargeContactTimeHist)==0
        SetFieldsFromData(Sim, ["ChargeContactTimeHist", "NegChargeToHistID", "PosChargeToHistID"])
       # computeChargeCorrelations(Sim)
    end

    #
   # Plots.plotlyjs()

    ChargeBondTime =  Plots.plot(size = (1150, 500),  yticks=([10^-8,10^-7,10^-6,10^-5,10^-4,0.001,0.01,0.1,1,10,100], ["  10^-8","  10^-7","  10^-6","  10^-5","  10^-4","  10^-3","  10^-2","  10^-1","  10^0","  10^1","  10^2"]),yaxis=:log, xaxis=:log, xlabel="Time ns", ylabel="P(r<2nm) [%]", grid=true, gridalpha=0.25, legend=:topright)

	step_to_time_in_ns = 10*100_000*10^-15/10^-9
	range = (axes(Sim.ChargeContactTimeHist,3)*step_to_time_in_ns)
	time = collect(0:1.0:length(axes(Sim.ChargeContactTimeHist,3))-1)

    for pos in keys(Sim.PosChargeToHistID)
		pos_id = Sim.PosChargeToHistID[pos]
		for neg in keys(Sim.NegChargeToHistID)
			neg_id = Sim.NegChargeToHistID[neg]
			hist = Sim.ChargeContactTimeHist[pos_id, neg_id,:]/sum(Sim.ChargeContactTimeHist[pos_id, neg_id,:])*100
			time_avg = sum(time.*hist)*step_to_time_in_ns/100
            Plots.plot!(range[1:150],hist[1:150], label="Charges: $(pos) , $(neg)  \n avg time: $(@sprintf("%.3f", time_avg)) ns", markershape=:cross )
		end
	end
    Plots.savefig(ChargeBondTime, Sim.PlotPath*Sim.SimulationName*"_ChargeBonds.png")

    Plots.savefig(ChargeBondTime, Sim.PlotPath*Sim.SimulationName*"_ChargeBonds.pdf")

	return ChargeBondTime
end

function HPSAnalysis.plotChargeBonds(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    if length(Sim.ChargeChainContactMatrix)==0
        SetFieldsFromData(Sim, ["ChargeChainContactMatrix", "MaxChainLength", "ChargeResidueContactMatrix"])
    #    computeChargeCorrelations(Sim)
    end
    

    Plots.gr() ### change Plots backend to get indepenent colorbars
	ChainMatrix = Plots.heatmap(1:Sim.NChains,1:Sim.NChains,Sim.ChargeChainContactMatrix, xlabel="Chain i", ylabel="Chain j" , title="Chain Contacts", colorbar=true)
	ResidueMatrix = Plots.heatmap(1:Sim.MaxChainLength,1:Sim.MaxChainLength,Sim.ChargeResidueContactMatrix, xlabel="Residue i", ylabel="Residue j" , title="Residue Contacts", colorbar=true)
	
	ChargeMatrix  = Plots.plot(ChainMatrix, ResidueMatrix, layout=2, size=(1200,500), left_margin=8mm, bottom_margin=5mm)

    Plots.savefig(ChargeMatrix, Sim.PlotPath*Sim.SimulationName*"_ChargeBondMatrix.png")
    Plots.savefig(ChargeMatrix, Sim.PlotPath*Sim.SimulationName*"_ChargeBondMatrix.pdf")

    return ChargeMatrix
end

function HPSAnalysis.plotMSDofChains(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Plots.gr()
    MSD_Plot=Plots.plot()
    MSD_Diff=Plots.plot()


    #MinLength = minimum(Sim.BoxLength)/2. ### visualisation of a improper behaviour of the Analysis after certain diffusion values
    #MaxDeviation = MinLength^2
    MSD_all = Plots.plot()
    for chain in 1:Sim.NChains			
        data = @views Sim.MSD[chain,:,:]
        Plots.plot!(axes(data,1), data[:,1]+data[:,2]+data[:,3], label="") 
    end


    N = length(axes(Sim.MSD,1))
    avg = sum(Sim.MSD, dims=1)[1,:,:]/Sim.NChains
    steps = length(axes(avg,1))
    diff_avg  = avg[2:steps,:]-avg[1:steps-1,:]

    avg_err =zero(avg) ### works like zeros_like in python
    diff_avg_err =zero(diff_avg) ### works like zeros_like in python
    diff = Sim.MSD[:,2:steps,:]-Sim.MSD[:,1:steps-1,:]
    for j in axes(Sim.MSD,3)
        for i in axes(Sim.MSD,2)
            avg_err[i,j] =  sum((avg[i,j] .- Sim.MSD[:,i,j]).^2)
        end
        for i in axes(diff,2)
            diff_avg_err[i,j] =  sum((diff_avg[i,j] .- diff[:,i,j]).^2)
        end
    end
    avg_err .= sqrt.(avg_err./N)/sqrt(N-1)
    diff_avg_err .= sqrt.(diff_avg_err./N)/sqrt(N-1)

    MaxVal = maximum(avg[1: ceil(I,Sim.NSteps*0.95)])
    ylims = [0, MaxVal]

    width=avg_err[:,1].+avg_err[:,2].+avg_err[:,3]
    Plots.plot!(MSD_all,axes(avg,1), avg[:,1]+avg[:,2]+avg[:,3]+width,fillrange=avg[:,1]+avg[:,2]+avg[:,3]-width, color=:gray, label="uncertainty")
    Plots.plot!(MSD_all,axes(avg,1), avg[:,1]+avg[:,2]+avg[:,3], color=:black, label="avg.", title="X&Y&Z-component", ylim=ylims.*3, legend=:bottomright)
    

    Plots.plot()
    for chain in 1:Sim.NChains
        data = @views Sim.MSD[chain,:,:]
        Plots.plot!(axes(data,1), data[:,1], label="") #+data[:,2]
    end
    Plots.plot!(axes(avg,1), avg[:,1]+avg_err[:,1],fillrange=avg[:,1]-avg_err[:,1], color=:gray, label="")
    MSD_XY = Plots.plot!( axes(avg,1), avg[:,1], color=:black, title="X-component" , label="", ylim=ylims) #+avg[:,2], title="X-Y"
    

    Plots.plot()
    for chain in 1:Sim.NChains
        data = @views Sim.MSD[chain,:,:]
        Plots.plot!(axes(data,1), data[:,2], label="") #+data[:,3]
    end
    Plots.plot!(axes(avg,1), avg[:,2]+avg_err[:,2],fillrange=avg[:,2]-avg_err[:,2], color=:gray, label="")
    MSD_XZ = Plots.plot!( axes(avg,1), avg[:,2], color=:black, title="Y-component" , label="", ylim=ylims)

    
    Plots.plot()
    for chain in 1:Sim.NChains
        data = @views Sim.MSD[chain,:,:]
        Plots.plot!(axes(data,1), data[:,3], label="") #data[:,2]+
    end
    Plots.plot!(axes(avg,1), avg[:,3]+avg_err[:,3],fillrange=avg[:,3]-avg_err[:,3], color=:gray, label="")
    MSD_YZ = Plots.plot!(axes(avg,1), avg[:,3], color=:black, title="Z-component" , label="" , ylim=ylims)
    MSD_Plot = Plots.plot(MSD_all, MSD_XY,MSD_XZ,MSD_YZ, layout=4, plot_title="MSD") # xscale=:log, yscale=:log,


 
    #max_step=Int32(min(500., Float32(length(axes(avg,1)))))
    #lims = [minimum(diff_avg[1:max_step,:])*0.95, maximum(diff_avg[1:max_step,:])*1.05]
    diff_ylims = [0, maximum(diff_avg[1: ceil(I,Sim.NSteps*0.75)])]
    Diff_all = Plots.plot()
    
    val=diff_avg[:,1]+diff_avg[:,2]+diff_avg[:,3]
    width =diff_avg_err[:,1]+diff_avg_err[:,2]+diff_avg_err[:,3]
    Diff_all = Plots.plot!(1:steps-1, val+width, fillrange=val-width, color=:gray, label="")
    Plots.plot!(1:steps-1, val, color=:black, label="" , title="X&Y&Z-component", ylims= diff_ylims.*3)


    Diff_XY  = Plots.plot(1:steps-1, diff_avg[:,1]+diff_avg_err[:,1], fillrange=diff_avg[:,1]-diff_avg_err[:,1], color=:gray, label="" , title="X-component", ylims= diff_ylims)
    Plots.plot!(1:steps-1, diff_avg[:,1], color=:black, label="" , title="X", ylims= diff_ylims)

    Diff_XZ  = Plots.plot(1:steps-1, diff_avg[:,2]+diff_avg_err[:,2], fillrange=diff_avg[:,2]-diff_avg_err[:,2], color=:gray, label="" , title="Y-component", ylims= diff_ylims)
    Plots.plot!(1:steps-1, diff_avg[:,2], color=:black, label="" , title="Y", ylims= diff_ylims)

    Diff_YZ  = Plots.plot(1:steps-1, diff_avg[:,3]+diff_avg_err[:,3], fillrange=diff_avg[:,3]-diff_avg_err[:,3], color=:gray, label="" , title="Z-component", ylims= diff_ylims)
    Plots.plot!(1:steps-1, diff_avg[:,3], color=:black, label="" , title="Z", ylims= diff_ylims)
    annotate!(Diff_all,(-0.2, 1.3), text("$(Sim.SimulationName)")) ### hacky way to get plot title that isnt set for subplots

    
    MSD_Diff = Plots.plot(Diff_all, Diff_XY, Diff_XZ, Diff_YZ, layout=4, plot_title="d/dt MSD")
    MSD_Master = Plots.plot(MSD_Plot , MSD_Diff, layout=2, size=(1600, 600))#,yscale=:log,xscale=:log)

    Plots.savefig(MSD_Master, Sim.PlotPath*Sim.SimulationName*"_MSD_Chains.png")
    Plots.savefig(MSD_Master, Sim.PlotPath*Sim.SimulationName*"_MSD_Chains.pdf")
    return MSD_Master
end

function HPSAnalysis.plotBondHist(Sim::SimData{R,I}) where {R<:Real, I<:Integer}

    #max_val = maximum(axes(Sim.BondHist),1)
    BondPlot= Plots.plot(axes(Sim.BondHist,1)/100, Sim.BondHist*100, xlabel="Bond distance r [AA]", ylabel="P(r) [%]", label="", xlim=(3,4.5))

    Plots.savefig(BondPlot, Sim.PlotPath*Sim.SimulationName*"_BondHist.png")
    Plots.savefig(BondPlot, Sim.PlotPath*Sim.SimulationName*"_BondHist.pdf")
    return BondPlot
end

@doc raw"""
    plotDihedralAngleHist(Sim::SimData{R,I})
Plot of the befor computed histogram of dihedral angles.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Create**:
* A histogram of dihedral angles.
"""
function HPSAnalysis.plotDihedralAngleHist(Sim)
    kα1 = 11.4
    kα2 = 0.15
    θα1 = 0.90
    θα2 = 1.02
    kβ1 = 1.80
    kβ2 = 0.65
    θβ1 = -1.55
    θβ2 = -2.50
    e0  =0.27
    e1  = 0.14
    e2  = 0.40
    Dihed_a(x, ed) = @. exp(-kα1*(x-θα1)^2-ed) + exp(-kα2*(x-θα2)^4+e0) + exp(-kα2*(x-θα2+2*π)^4+e0)
    Dihed_b(x, ed) = @. exp(-kβ1*(x-θβ1)^2+e1+ed) + exp(-kβ1*(x-θβ1-2*π)^2+e1+ed) + exp(-kβ2*(x-θβ2)^4+e2) + exp(-kβ2*(x-θβ2-2*π)^4+e2)
    U(x,ed) = @. -log(Dihed_a(x, ed)+Dihed_b(x, ed))

    NA = 6.022
    kb = 1.38
    kcal = 4184
    kt = kb *Sim.TargetTemp *NA/kcal

    x = LinRange(-π,π,500)
    E = U(x, -1.0)
    W = exp.(-1.0*E/kt)#.*sin.( x)
    P = W./sum(W)*(2*π/500.0)^-1

    E = U(x, 1.0)
    W = exp.(-1.0*E/kt)#.*sin.( x)
    P2 = W./sum(W)*(2*π/500.0)^-1


    DihedralAnglePlot= Plots.plot(axes(Sim.TorsionHist,2), Sim.TorsionHist[end,:]*100, xlabel="Diehedral Angle Ψ [Degree]", ylabel="P(Ψ) [%]", title="Dihedral angle histogram", label="")
    Plots.plot!(( x)*180.0/π, (P+P2)*100/2.0, label="Probability")

    #DihedralAnglePlot= Plots.plot( xlabel="Diehedral Angle Ψ [Degree]", ylabel="P(Ψ) [%]", title="Dihedral angle histogram", label="")
    #=for i in 2:length(axes(Sim.TorsionHist,1))
        Plots.plot!(axes(Sim.TorsionHist,2), Sim.TorsionHist[i,:]*100, label="")
    end=#
    Plots.savefig(DihedralAnglePlot, Sim.PlotPath*Sim.SimulationName*"_DihedralAngleHist.png")
    Plots.savefig(DihedralAnglePlot, Sim.PlotPath*Sim.SimulationName*"_DihedralAngleHist.pdf")
    return DihedralAnglePlot
end

@doc raw"""
    plotBondAngleHist(Sim::SimData{R,I})

Plot of the histogram of bond angles.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Create**:
* A histogram of bond angles.
"""
function HPSAnalysis.plotBondAngleHist(Sim)
    BondAnglePlot= Plots.plot(axes(Sim.BondAngleHist,1)*180.0/size(Sim.BondAngleHist,1), Sim.BondAngleHist[:]*100, xlabel="Bond Angle Ψ [Degree]", ylabel="P(Ψ) [%]", title="Bond angle histogram", label="")

    Plots.savefig(BondAnglePlot, Sim.PlotPath*Sim.SimulationName*"_BondAngleHist.png")
    Plots.savefig(BondAnglePlot, Sim.PlotPath*Sim.SimulationName*"_BondAngleHist.pdf")
    return BondAnglePlot
end

function HPSAnalysis.plotRGHist(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Plots.plot()
    RG_Hist = Plots.histogram([(Sim.RGSeries...)...], normalize=:pdf, xlabel="RG [AA]", ylabel="P(RG)", label="") #,, bins=b_range color=:gray)
    Plots.savefig(RG_Hist, Sim.PlotPath*Sim.SimulationName*"_RGHist.png")
    Plots.savefig(RG_Hist, Sim.PlotPath*Sim.SimulationName*"_RGHist.pdf")
    return RG_Hist
end

function HPSAnalysis.plotRGSeries(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    RG_Hist= Plots.plot()
    axis = axes(Sim.RGSeries[1,:])
    for chain in 1:Sim.NChains
        RG_Hist = Plots.plot!(axis, Sim.RGSeries[chain,:], xlabel="t [?]", ylabel="RG [AA]", label="")
    end
    Plots.vline!([Sim.EquilibrationTime], label="Choosen EquilibrationTime")
    Plots.savefig(RG_Hist, Sim.PlotPath*Sim.SimulationName*"_RGSeries.png")
    Plots.savefig(RG_Hist, Sim.PlotPath*Sim.SimulationName*"_RGSeries.pdf")
    return RG_Hist
end

function HPSAnalysis.plotRGAutocorr(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    RG_Auto= Plots.plot(xlabel="lag time [τ]", ylabel="autocorr []", xlims=(0,500),ylims=(-0.5, 1.0))
    axis = 0:length(Sim.RGAutocorr[1,:])-1
    avg = sum(Sim.RGAutocorr, dims=1)[1,:]./Sim.NChains

    for chain in 1:Sim.NChains
        RG_Auto = Plots.plot!(axis, Sim.RGAutocorr[chain,:], label="")
    end
    Plots.plot!(axis, avg,  label="avg. over chains", c=:black, xlim=(0,length(Sim.RGAutocorr[1,:])))

    #Plots.vline!([Sim.EquilibrationTime], label="Choosen EquilibrationTime")
    Plots.savefig(RG_Auto, Sim.PlotPath*Sim.SimulationName*"_RGAutocorr.png")
    Plots.savefig(RG_Auto, Sim.PlotPath*Sim.SimulationName*"_RGAutocorr.pdf")
    return RG_Auto
end

function HPSAnalysis.plotREEHist(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    Plots.plot()
    REE_Hist = Plots.histogram([(Sim.REESeries...)...], normalize=:pdf, xlabel="REE [AA]", ylabel="P(REE)", label="") #,, bins=b_range color=:gray)
    Plots.savefig(REE_Hist, Sim.PlotPath*Sim.SimulationName*"_REEHist.png")
    Plots.savefig(REE_Hist, Sim.PlotPath*Sim.SimulationName*"_REEHist.pdf")
    return REE_Hist
end

function HPSAnalysis.animateDensityGif(Sim::SimData{R,I}, NBins=100) where {R<:Real, I<:Integer}

    conversion = 1.66053906660/(Sim.BoxLength[1]*Sim.BoxLength[2]/Float32(NBins)^2*Sim.BoxLength[3])
    weights= Sim.Masses*conversion
    anim = @animate for step ∈ 1:10:Sim.NSteps
        fig=Plots.plot(layout=@layout[a{0.35w} b{0.2w} c{0.45w}] , size=(1750/1.6,500/2), left_margin=5mm, right_margin=2.5mm, margin=0mm, bottom_margin=6mm)
        Plots.histogram2d!(Sim.y[:,step], Sim.x[:, step], weights=weights, nbins = NBins, show_empty_bins=true, clims=(0,1.5),  xlabel="y", ylabel="x", subplot=1, colorbar=:false,xlims=Tuple(Sim.BoxSize[2,:]), ylims=Tuple(Sim.BoxSize[1,:]))
        Plots.histogram2d!(Sim.x[:, step], Sim.z[:,step], weights=weights, nbins = NBins, show_empty_bins=true, clims=(0,1.5), xlabel="x", ylabel="z", subplot=2, colorbar=:false, xlims=Tuple(Sim.BoxSize[1,:]), ylims=Tuple(Sim.BoxSize[3,:]))
        Plots.histogram2d!(Sim.y[:, step], Sim.z[:,step], weights=weights, nbins = NBins, show_empty_bins=true, clims=(0,1.5), xlabel="y", ylabel="z", subplot=3, colorbar=:right,colorbar_title="kg/L", xlims=Tuple(Sim.BoxSize[2,:]), ylims=Tuple(Sim.BoxSize[3,:]))
        Plots.annotate!(1.55*Sim.BoxSize[2,2], 1.15*Sim.BoxSize[3,1], "$(step/1000.0) μs", subplot=3)
        fig
    end
    gif(anim, "/localscratch/HPS_DATA/HPS-Alpha/SLAB/DensityGifs/$(Sim.SimulationName)_$(Sim.TargetTemp).gif", fps = 30)
end


@doc raw"""
    plotAvgSlabDensity(Sim::SimData; Windowlength=100)

Plots the average slab histogram of the frames in range Sim.NSteps-Windowlength:Sim.NSteps. Result is empty if the Windowlength is larger then the step size at which histgrams here computed.
"""
function HPSAnalysis.plotAvgSlabDensity(Sim::HPSAnalysis.SimData{R,I};size=(400, 400), IndexDict=Dict(1=>"all"), Windowlength=100, title=nothing, normalize=false,xlims=nothing, filename=nothing) where {R<:Real, I<:Integer}
    if Windowlength > Sim.NSteps
        Windowlength=Sim.NSteps
    end

    filename = isnothing(filename) ?  Sim.SimulationName : filename

    xaxis = axes(Sim.SlabHistogramSeries)[1]
    ### newer iterations dont measure at every frame. Detect which frames actually contain data
    NMeasurements= sum(Sim.SlabHistogramSeries[1,Sim.NSteps-Windowlength+1:Sim.NSteps,1].!=0.0)
    AvgHist = sum(Sim.SlabHistogramSeries[xaxis,Sim.NSteps-Windowlength:Sim.NSteps,:], dims=2)./(NMeasurements)

    if isnothing(title)
        split_ =  split(Sim.SimulationName[3:end],"_")
        title= "$(split_[1]) @ $(split_[2])"
    end

    ylabel_ = normalize ?  L"\langle \rho(z)\rangle/\langle \rho(z=0)\rangle" : L"\langle \rho(z)\rangle\quad[\textrm{kg/L}]"
    fig = Plots.plot(dpi=300, size=size, ylabel= ylabel_ , xlabel= L"z\quad[\textrm{Å}]", minorticks=10, gridalpha=1.0) 
    for (index, label) in IndexDict
        mean_ = normalize ? mean(AvgHist[-10:10,1,index]) : 1
        Plots.plot!(xaxis,AvgHist[:,1,index]./mean_,  label=label, title= title, legendtitlefontsize=9, c= index==1 ? :black : index, legend= normalize ? :bottom : :best, legend_column= normalize ?  2 : 1) # color=:black,
    end

    addon= normalize ?  "_normed" : ""
    if !isnothing(xlims)
        Plots.plot!(xlims=xlims,xticks=xlims[1]:100:xlims[2])
        addon *= "_lims"
    end

    Plots.savefig(fig, Sim.PlotPath*filename*"_AvgSlabHist_$(Windowlength)$(addon)_LastFrames.png")
    Plots.savefig(fig, Sim.PlotPath*filename*"_AvgSlabHist_$(Windowlength)$(addon)_LastFrames.pdf")
    return fig
end

@doc raw"""
    plotAvgSlabDensityEvolution(Sim::SimData{R,I}; Windowlength=100) where {R<:Real, I<:Integer}

Plot of the average density in the simulation box.

Plots the average slab histogram in range n\*Windowlength:(n+1)\*Windowlength starting at Sim.EquilibrationTime. Some slices will be empty if the Windowlength is larger then the step size at which histgrams were computed.

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.
- `Windowlength::Int`: How much steps are processed.

**Create**:
* A histogram of the average density.
"""
function HPSAnalysis.plotAvgSlabDensityEvolution(Sim::SimData{R,I}; Windowlength=100, dilute_pos=0, dense_pos=0) where {R<:Real, I<:Integer}
    if Windowlength > Sim.NSteps
        Windowlength=Sim.NSteps
    end
    
    xval= LinRange(Sim.BoxSize[Sim.SlabAxis,1], Sim.BoxSize[Sim.SlabAxis,2], 9)
    xticks = (xval.-Sim.BoxSize[Sim.SlabAxis,1], string.(xval))
    xaxis = axes(Sim.SlabHistogramSeries)[1]

    tmp = ceil(I,Sim.NSteps/Windowlength)
    AvgHist = zeros(R,(size(Sim.SlabHistogramSeries)[1], size(Sim.SlabHistogramSeries)[3], tmp) )
    for i in 1:tmp-1
        start = 1 + (i-1)*Windowlength
        stop = start + Windowlength-1
        ### newer iterations dont measure at every frame. Detect which frames actually contain data
        NMeasurements= sum(Sim.SlabHistogramSeries[0,start:stop,1].!=0.0)
        AvgHist[:,:,i] .= inter_vals.no_offset_view((reduce(+, Sim.SlabHistogramSeries[xaxis,start:stop,:], dims=2)./(NMeasurements)))[:,1,:] ### convert offsetarray back so that we can assign it to AvgHist
    end
    start = (tmp-1)*Windowlength+1
    stop=Sim.NSteps
    NMeasurements= sum(Sim.SlabHistogramSeries[0,start:stop,1].!=0.0)
    AvgHist[:,:,end] .= OffsetArrays.no_offset_view((reduce(+, Sim.SlabHistogramSeries[xaxis,start:stop,:], dims=2)./(NMeasurements))[:,1,:])

    fig = Plots.plot(dpi=300,  ylabel= "avg. density\n"* "  [kg/L]" , xlabel= "z-Axis [Å]"  ,  title="$(Sim.SimulationName)",palette=:darktest, seriescolor=:darktest, color=:darktest, markercolor=:darktest, linecolor=:darktest, zcolor=:darktest, xticks=xticks, colorbar_title="window No.", levels=tmp)

    for i in 1:tmp
        if (i)*Windowlength < Sim.EquilibrationTime continue end
        Plots.plot!(axes(AvgHist)[1][2:end-1], AvgHist[:,1,i][2:end-1],  label="",  line_z=i, marker_z=i, levels=tmp, c=:darktest)
    end

    center = maximum(axes(Sim.SlabHistogramSeries,1))
    if dilute_pos!=0
        Plots.vline!([center-dilute_pos,center+dilute_pos], c=:black, label="")
    end

    if dense_pos!=0
        Plots.vline!([center-dense_pos, center+dense_pos], c=:black, label="")
    end

    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_$(Sim.TargetTemp)_AvgSlabHist_Evolution.png")
    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_$(Sim.TargetTemp)_AvgSlabHist_Evolution.pdf")

    if dilute_pos!=0
        pos = (center-dilute_pos, center+dilute_pos)
        start = ceil(-dilute_pos/100)*100
        xval= start:100:dilute_pos
        xticks = (xval.-Sim.BoxSize[Sim.SlabAxis,1], [@sprintf("%.0f", val) for val in xval])
        Plots.plot!(xlim=pos, xticks=xticks, xminorticks=10, grid=true, minorgrid=true, gridalpha=0.8, minorgridalpha=0.4)#, xrotation=90)
        Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_$(Sim.TargetTemp)_AvgSlabHist_Evolution_zoom.png")
        Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_$(Sim.TargetTemp)_AvgSlabHist_Evolution_zoom.pdf")
    end


    return fig
end

function HPSAnalysis.plotSurfaceWidthFit(xvals, yvals, model, SaveAs; size=(400,300))
    fig = Plots.plot(xvals,yvals, size=size,label="sim. data", xlabel=L"|z-z_\textrm{drop.}|\quad[\AA]", ylabel=L"\rho(z)\quad[\textrm{kg/L}]", minorticks=10,gridalpha=1.0)
    Plots.plot!(xvals,HPSAnalysis.tanh_profile(xvals, model.param),label="tanh. fit", xlim=(0,2*model.param[4]))

    Plots.savefig(fig, SaveAs*".png")
    Plots.savefig(fig, SaveAs*".pdf")

    return fig
end

function HPSAnalysis.plotFitSeries(axis, Series, Func, Fit, SaveAs; xlabel="", ylabel="", label="")
    FitPlot = Plots.plot()
    Plots.plot!(axis, Func(axis, Fit.param), xlabel=xlabel, ylabel=ylabel, label=label)
    Plots.plot!(axis, Series, label="data")
    Plots.savefig(FitPlot, SaveAs*".png")
    Plots.savefig(FitPlot, SaveAs*".pdf")
    return FitPlot
end

function HPSAnalysis.plotSwapSteps(ChargeFile::String)
    Steps, Charges, ID, Restarts = HPSAnalysis.readChargeAssignments(ChargeFile)

    fig = Plots.plot(xlabel="Steps", ylabel="q_i"*" at Sim i")
    for i in axes(ID,1) Plots.plot!(Steps[1:end-1],  Charges[i,1:end-1], label="") end
    Plots.vline!(Restarts, c=:black, label="")
    return fig
end

function HPSAnalysis.plotHREMDMSDofChains(Data::HPSAnalysis.HREMD_Data{R,I}) where {R<:Real, I<:Integer}
    Plots.gr()
    MSD_Plot= Plots.plot()
    MSD_Diff= Plots.plot(title="X-Y-Z")
    MSD_all = Plots.plot()

    Diff_all = Plots.plot( title="X-Y-Z")


    #MinLength = minimum(Sim.BoxLength)/2. ### visualisation of a improper behaviour of the Analysis after certain diffusion values
    #MaxDeviation = MinLength^2
    MSD_XY = Plots.plot()
    MSD_XZ = Plots.plot()
    MSD_YZ = Plots.plot()

    Diff_XY = Plots.plot(title="X")
    Diff_XZ = Plots.plot(title="Y")
    Diff_YZ = Plots.plot(title="Z")
    for Sim in Data.ReplicaData

        avg = sum(Sim.MSD, dims=1)[1,:,:]/Sim.NChains
        MaxVal = maximum(avg)
        ylims = [1, MaxVal]

        Plots.plot!(MSD_all, axes(avg,1), avg[:,1]+avg[:,2]+avg[:,3], color=:black, label="", ylim=ylims.*3)

    #  Plots.hline!([MaxDeviation], label="", linewidth=1; color=:red)
        Plots.plot!(MSD_XY, axes(avg,1), avg[:,1], color=:black, title="X" , label="", ylim=ylims) #+avg[:,2], title="X-Y"

        #Plots.hline!([MaxDeviation], label="", linewidth=1; color=:red)
        Plots.plot!(MSD_XZ, axes(avg,1), avg[:,2], color=:black, title="Y" , label="", ylim=ylims)


        #Plots.hline!([MaxDeviation], label="", linewidth=1; color=:red)
        Plots.plot!(MSD_YZ, axes(avg,1), avg[:,3], color=:black, title="Z" , label="" , ylim=ylims)


        steps = length(axes(avg,1))
        diff_avg  = avg[2:steps,:]-avg[1:steps-1,:]
        #max_step=Int32(min(500., Float32(length(axes(avg,1)))))
        #lims = [minimum(diff_avg[1:max_step,:])*0.95, maximum(diff_avg[1:max_step,:])*1.05]
        lims = [0,5.0]

        Plots.plot!(Diff_all, 1:steps-1, diff_avg[:,1]+diff_avg[:,2]+diff_avg[:,3], color=:black, label="" , ylims= lims)
        Plots.plot!(Diff_XY, 1:steps-1, diff_avg[:,1], color=:black, label="" , ylims= lims)
        Plots.plot!(Diff_XZ, 1:steps-1, diff_avg[:,2], color=:black, label="" , ylims= lims)
        Plots.plot!(Diff_YZ, 1:steps-1, diff_avg[:,3], color=:black, label="" ,  ylims= lims)


    end
    Sim = Data.ReplicaData[1]

    Plots.annotate!(Diff_all,(-0.2, 1.3), Plots.text("$(Sim.SimulationName)")) ### hacky way to get plot title that isnt set for subplots

    MSD_Plot = Plots.plot(MSD_all, MSD_XY,MSD_XZ,MSD_YZ, layout=4, plot_title="MSD") # xscale=:log, yscale=:log,
    MSD_Diff = Plots.plot(Diff_all, Diff_XY, Diff_XZ, Diff_YZ, layout=4, plot_title="d/dt MSD")
    Plots.vline!([50, 100,150,200], label="restarts")
    MSD_Master = Plots.plot(MSD_Plot , MSD_Diff, layout=2, size=(1600, 600),)


    Plots.savefig(MSD_Master, Sim.PlotPath*"HREMD_MSD_Chains.png")
    Plots.savefig(MSD_Master, Sim.PlotPath*"HREMD_MSD_Chains.pdf")

    return MSD_Master
end


@doc raw"""
    plotDensityHistogram(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
Plots the logarithmic density histogram computed by [`computeDensityHistogram`](@ref).

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.

**Create**:
* A time averaged histogram of the average density per subcube of the simulation box.
"""
function HPSAnalysis.plotDensityHistogram(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    res = 80
    x=10.0.^((collect(axes(Sim.DensityHist,1)).-(6*res))./res)#./NRes

    fig = Plots.plot(x[2:end], Sim.DensityHist[2:end], label="", xlim=(10^-3,2), ylim=(0,maximum(Sim.DensityHist[10:end])), xscale=:log10, xlabel="m [kg/L]", ylabel="P(m) []", minorticks=true)

    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_$(Sim.TargetTemp)_DensityHist.png")
    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_$(Sim.TargetTemp)_DensityHist.pdf")
end


@doc raw"""
    plotIntraChainScaling(Sim::SimData{R,I}) where {R<:Real, I<:Integer}

Plots the square intra chain distance |r_i-r_j| on logarithmic scales as computed by [`computeIntraChainScalingNaiv`](@ref) and [`computeIntraChainScalingSlidingWindow`](@ref).

**Arguments**:
- `Sim::SimData{R,I}`: A simulation data structure containing the Simulation information.
"""
function HPSAnalysis.plotIntraChainScaling(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    fig = Plots.plot( xlim=(1, maximum(length.(Sim.IntraChainScalingSlidingWindow))), ylim=(1, maximum(maximum.(Sim.IntraChainScalingSlidingWindow))*1.2),  xlabel="i []", ylabel="|r_i-r_j|^2 [AA]", xscale=:log10, yscale=:log10, minorticks=true)
    
    for (S, Sequence) in enumerate(Sim.IntraChainScalingNaiv)
        for i in 1:size(Sequence, 1)
            Plots.plot!(Sequence[i,i+1:end], label= S==1 && i==1 ? "Naiv" : "", c=:gray)
        end
    end

    for (S, Sequence) in enumerate(Sim.IntraChainScalingSlidingWindow)
        Plots.plot!(Sequence, label= S==1 ? "Window" : "", c=:black)
    end

    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_IntraChainScaling.png")
    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_IntraChainScaling.pdf")
end

function HPSAnalysis.plotHOOMDTemperature(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    try
        file = h5open("$(Sim.BasePath)pressure.h5", "r")
        temperature = collect(file["hoomd-data"]["md"]["compute"]["ThermodynamicQuantities"]["kinetic_temperature"])
        timestep = collect(file["hoomd-data"]["Simulation"]["timestep"])
        close(file)
        kb = 0.00831446262
        fig = Plots.plot(timestep, temperature/kb, xlabel="step", ylabel="temperature", ylim=(Sim.TargetTemp*0.95, Sim.TargetTemp*1.05), label="Sim. Temp.")
        Plots.hline!([Sim.TargetTemp], label="Target Temp.")
        Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_KineticTemperature.png")
        Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_KineticTemperature.pdf")
    catch e 
        println(e)
    end
end

function HPSAnalysis.plotIntraChainContactMatrix(Sim::SimData{R,I}) where {R<:Real, I<:Integer}
    N = ceil(Int32, sqrt(Sim.NChains))
    fig=Plots.plot(layout=(N,N), xlabel="Amino Acids i [-]", ylabel="Amino Acids j [-]")

    for C in 1:Sim.NChains
        Plots.heatmap!(Sim.IntraChainContactMatrix[C], subplot=C, zscale=:ln)
    end
    
    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_IntraChainContactMatrix.png")
    Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_IntraChainContactMatrix.pdf")
end


function get_color(value, clims, scheme; zero_centered::Bool=false)
    if isnan(value)
        return RGBA(0.0, 0.0, 0.0, 0.0)
    end
    if zero_centered
        lo, hi = clims
        norm_val = value >= 0 ? 0.5 + 0.5*clamp(hi > 0 ? value/hi : 0.0, 0.0, 1.0) :
                                 0.5 - 0.5*clamp(lo < 0 ? value/lo : 0.0, 0.0, 1.0)
    else
        norm_val = clamp((value - clims[1]) / (clims[2] - clims[1]), 0.0, 1.0)
    end
    return colorschemes[scheme][norm_val]
end

# a colour gradient whose neutral/middle colour is pinned at value 0, even for asymmetric
# clims, so Plots' own linear clims normalization (heatmap `color=`/`clims=`) puts the
# scheme's neutral colour exactly at the 0 tick instead of at the midpoint of clims.
# Finely resamples the whole scheme through the same piecewise remap get_color(zero_centered=true)
# uses, rather than interpolating between just the 3 anchor colours (which would skip any
# intermediate hues the scheme has between its ends and its middle).
function zero_centered_cgrad(clims, scheme; nstops::Integer=101)
    lo, hi = clims
    pos0 = clamp((0.0 - lo) / (hi - lo), 0.0, 1.0)
    orig_pos(p) = p <= pos0 ? (pos0 > 0 ? 0.5*p/pos0 : 0.5) :
                              0.5 + 0.5*(p - pos0)/(pos0 < 1 ? 1 - pos0 : 1)
    ps = range(0.0, 1.0, length=nstops)
    return cgrad([colorschemes[scheme][orig_pos(p)] for p in ps], ps)
end

function HPSAnalysis.compute_submatrix(matrix::Matrix{R},matrix_err::Matrix{R}, ranges, volume::Vector{<:Real}=Float64[]) where {R<:Real}
    N = size(matrix, 1)
    N_sub = length(ranges)
    small_matrix = zeros(R, (N_sub,N_sub))
    large_matrix = zeros(R, (N,N))
    small_matrix_err = zeros(R, (N_sub,N_sub))
    large_matrix_err = zeros(R, (N,N))

    weighted = !isempty(volume)
    weighted && @assert length(volume) == N "'volume' needs one entry per row/column of 'matrix'."

    for (i,range_i) in enumerate(ranges)
        for (j, range_j) in enumerate(ranges)
            if weighted
                w = volume[range_i] * transpose(volume[range_j])
                mean_val = sum(matrix[    range_i, range_j] .* w) / sum(w)
                err_val  = 0.0 # sum(matrix_err[range_i, range_j] .* w) / sum(w)
            else
                mean_val = mean(matrix[   range_i, range_j])
                err_val = sqrt(sum(matrix_err[range_i, range_j].^2))/length(matrix_err[range_i, range_j])
            end
            large_matrix[    range_i,range_j] .= mean_val
            large_matrix_err[range_i,range_j] .= err_val
            small_matrix[    i,j] = mean_val
            small_matrix_err[i,j] = err_val
        end
    end
    return small_matrix, large_matrix, small_matrix_err, large_matrix_err
end

function HPSAnalysis.plotInterChainContactMatrix(Sim::SimData{R,I}; Size=500, ranges=Vector{UnitRange{<:Integer}}(), labels=["",""],inPercent=false) where {R<:Real, I<:Integer}

    IntraChainMatrix = sum(Sim.IntraChainContactMatrix)./length(Sim.IntraChainContactMatrix)
    
    ### replace IntraChainContactMatrix in case that 
    
    color_inter=:thermal

    ly = @layout([a{0.9w, 0.03h} b{0.1w, 0.03h}; c{0.95w, 0.95h} d{0.03w, 0.95h}])

    InterChainMatrix = Matrix{R}(Sim.ContactMatrices)
    InterChainMatrixError = Matrix{R}(Sim.ContactMatricesError)

    ### plot once with errors and once with normal values
    for (i, (avg, iserror, label_)) in enumerate(zip([IntraChainMatrix , Sim.ContactMatricesError], ["", "_Error"],labels))
        fig=Plots.plot(xlabel="Amino Acids i ", ylabel="Amino Acids j ", colorbar=:legend, legend=:outertop,)

        addon=""
        if length(ranges)>0
            _, InterChainMatrix,_,_ = HPSAnalysis.compute_submatrix(InterChainMatrix,InterChainMatrixError, ranges)
            addon ="_with_coarse"
        end
        if any(isnan.(InterChainMatrix)) || any(isnan.(avg))
            println("Skip plot.")
            continue
        end

        clims_inter = (0,maximum(InterChainMatrix)) # a
        clims_intra = iserror=="_Error" ? clims_inter : extrema(avg) #b
        toplabel    = iserror=="_Error" ? (inPercent ? L"ΔP(d_{ij}<1.1\cdotσ\cdot2^{(1/6)})\quad[\%]" :  L"ΔP(d_{ij}<1.1\cdotσ\cdot2^{(1/6)})") : L"\ln(\langle d_{ij}\rangle)"
        color_intra = iserror=="_Error" ? color_inter : :viridis

        rightlabel = inPercent ? L"P(d_{ij}<1.1\cdotσ\cdot2^{(1/6)})\quad[\%]" : L"P((d_{ij}<1.1\cdotσ\cdot2^{(1/6)})"

        matrix  = triu(InterChainMatrix)+tril(avg,-1);
        is_triu = triu(trues(size(InterChainMatrix)), 1)
        color_mat = [istriu ? get_color(val, clims_inter, color_inter) : get_color(val, clims_intra, color_intra) for (istriu, val) in zip(is_triu, matrix)]

        clims_inter = inPercent ? clims_inter.*100.0 : clims_inter
        clims_intra = (inPercent && iserror=="_Error") ? clims_intra.*100.0 : clims_intra

        data_inter = collect(range(clims_intra..., 100))
        data_intra = collect(range(clims_inter..., 100))

        clims = extrema(matrix)
        NMat = size(matrix,1)
        ticks = vcat([1], collect(0:50:NMat)[2:end])#, [NMat])

        fig = plot(layout=ly, size=(Size,Size), margin=-2.5mm, AspectRatio=true, minorticks=10, tick_direction=:out)#, link=:y )

        heatmap!(color_mat,yflip=false,lims=(0.5,NMat+0.5), ticks=ticks, colorbar=true, color = color_inter, clims=clims_intra, subplot=3, xlabel="Amino Acids i ", ylabel="Amino Acids j ", AspectRatio=false,yrotation=90.0, fontsize=10)#, yaxis = :flip)
        heatmap!(data_inter,[1],  transpose([data_inter;;]),color=color_intra,   colorbar=false, subplot=1, yaxis=nothing, xmirror=true, xlabel=toplabel ,labelfontsize=8)
        heatmap!(subplot=2, framestyle=:none)
        heatmap!([1], data_intra, [data_intra;;], colorbar=false, clims=clims_inter, color=color_inter, xlabel="",xticks=nothing,   ymirror=true,  yguide=nothing, subplot=4, labelfontsize=8, yrotation=90.0, ylabel=rightlabel)

        Plots.annotate!([-0.075*NMat], [1.1*NMat], [label_])

        Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_InterChainContactMatrix$(addon)$(iserror).png")
        Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_InterChainContactMatrix$(addon)$(iserror).pdf")
        
        if i==2
            return fig
        end
    end
end

zerotonan(x) = isapprox(x,0,atol=10^-7,rtol=10^-7) ? NaN : x

function getTypeNameDictWithEndingNames(Sim::SimData{R,I}; RenamePhosAA=Dict('#'=>"pS",'&'=> "pT", '@'=>"pY"))::Dict{I,String} where {R<:Real,I<:Integer} 
    NIds = maximum(keys(Sim.IDToResName))
    local_IDToResName = copy(Sim.IDToResName)
    non_real_aa = [id for id in 1:NIds if islowercase(Sim.IDToResName[id][1])]
    NIds_total = BioData.AaToWeight
    types =collect(keys(BioData.AaToWeight))
    for aa in non_real_aa      
        key = findfirst(x->(BioData.AaToWeight[x]+2≈Sim.IDToMasses[aa] || BioData.AaToWeight[x]+16≈Sim.IDToMasses[aa]) ,types)
        if !(key ∈ ['#','&','@'])
            local_IDToResName[aa] = lowercase("$(types[key])")
        end
        if key ∈ ['#','&','@']
            local_IDToResName[aa] = RenamePhosAA[aa]
        end
    end
    return local_IDToResName
end

function HPSAnalysis.plotMeanEnergyPerID(Sim::SimData{R,I}; inter=true)::Plots.Plot{Plots.GRBackend} where {R<:Real,I<:Integer} 
    fig =  Plots.plot(xlabel="Amino Acid i", ylabel="Amino Acid j",title=L"\quad\langle E \rangle~\textrm{per~Sequence~Index}~[\textrm{J/mol}]")

    if inter
        matrix  = (triu(Sim.MeanYukawaEnergy)+tril(Sim.MeanAshbaughHatchEnergy,-1)).*1000.0
    else
        matrix  = (triu(Sim.MeanIntraYukawaEnergy)+tril(Sim.MeanIntraAshbaughHatchEnergy,-1)).*1000.0
    end
    max = maximum(abs.(matrix); init=0.0)
    if max==0 return fig end
    NMat = size(matrix,1)
    ticks = vcat([1], collect(0:50:NMat)[2:end])
    Plots.heatmap!(zerotonan.(matrix),c=:vik,clims=(-max,max),size=(400,400), ticks=ticks, AspectRatio=true, tick_direction=:out, minorticks=10)
    return fig
end

function HPSAnalysis.plotMeanEnergyPerType(Sim::SimData{R,I};inter=true)::Plots.Plot{Plots.GRBackend} where {R<:Real,I<:Integer} 
    IDToResNames = getTypeNameDictWithEndingNames(Sim)

    NIds = maximum(keys(IDToResNames))
    ticks= [IDToResNames[id] for id in 1:NIds ]

    fig =  Plots.plot(xlabel="Amino Acid Type",ylabel="Amino Acid Type", title=L"\quad\langle E \rangle~\textrm{per~Amino~Acid~Type}~[\textrm{J/mol}]", ticks=(1:NIds, ticks)) 
    if inter
        matrix  = (triu(Sim.MeanYukawaEnergy_perAA)+tril(Sim.MeanAshbaughHatchEnergy_perAA,-1)).*1000
    else
        matrix  = (triu(Sim.MeanIntraYukawaEnergy_perAA)+tril(Sim.MeanIntraAshbaughHatchEnergy_perAA,-1)).*1000
    end
    max = maximum(abs.(matrix); init=0.0)
    if max==0 return fig end

    Plots.heatmap!(zerotonan.(matrix),c=:vik,clims=(-max,max),AspectRatio=true, size=(400,400),right_margin=1mm,colorbar_title_location=:top,colorbar_tickfontrotation=0, tick_direction=:out)
    return fig
end

function HPSAnalysis.plotMeanEnergies(Sim::SimData{R,I}; annotatePlot=true) where {R<:Real,I<:Integer} 
    figs =[]
    for (inter,addon) in [(true,""), (false,"_intra")]
        figa= plotMeanEnergyPerID(Sim; inter=inter)
        xlims=Plots.xlims(figa)
        if annotatePlot annotate!(figa,[-xlims[2]/5.0],[xlims[2]],["a)"], subplot=1) end

        figb= plotMeanEnergyPerType(Sim; inter=inter)
        xlims=Plots.xlims(figb)
        if annotatePlot annotate!(figb,[-xlims[2]/6.5],[xlims[2]],["b)"], subplot=1) end

        fig  = Plots.plot(figa,figb , size=(800,400), bottom_margin=2mm, left_margin=2mm)

        Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_MeanEnergies$(addon).png" )
        Plots.savefig(fig, Sim.PlotPath*Sim.SimulationName*"_MeanEnergies$(addon).pdf" )
        push!(figs, fig)
    end
    return Plots.plot(figs..., layout=(2,1), size=(800,800))
end
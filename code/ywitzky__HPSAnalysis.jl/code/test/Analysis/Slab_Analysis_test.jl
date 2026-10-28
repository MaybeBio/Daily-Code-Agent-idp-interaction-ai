using OffsetArrays


### construct SimData for 2 chains of different length inside slab box
Sim = HPSAnalysis.SimData()

Sim.NAtoms=61+6
Sim.NChains=2
Sim.IDs = ones(Sim.NAtoms)
Sim.Masses = ones(Float32, Sim.NAtoms)
Sim.ChainStart= [1, 62]
Sim.ChainStop = [61, Sim.NAtoms]
Sim.ChainMasses = [61, 6]
TotalMass= sum(Sim.ChainMasses)
Sim.Charges = ones(Sim.NAtoms)
Sim.BoxSize=[-10.0 10.0; -100.0 100.0; -10.0 10.0]
Sim.BoxLength=[20.0,200.0,20.0]
Sim.NSteps = 100
Sim.SlabAxis=2
Sim.x = zeros(Sim.NAtoms, Sim.NSteps)
Sim.y = zeros(Sim.NAtoms, Sim.NSteps)
Sim.z = zeros(Sim.NAtoms, Sim.NSteps)

Sim.x_uw = zeros(Sim.NAtoms, Sim.NSteps)
Sim.y_uw = zeros(Sim.NAtoms, Sim.NSteps)
Sim.z_uw = zeros(Sim.NAtoms, Sim.NSteps)


pos_C1 = Float32.(collect(-30:30))
pos_C2 = Float32.(collect(-45:-40))

offset = collect(2.0*(0:Sim.NSteps-1))
for i in 1:Sim.NSteps
    Sim.y_uw[1:61,i] = offset[i] .+ pos_C1
    Sim.y_uw[62:end,i] = offset[i] .+ pos_C2
end

Sim.EquilibrationTime=1
Sim.RGMeasureStep = 1

### definition how true slab histogram looks like
hist_def = OffsetArray(zeros(200,3), -99:100, 1:3)
for i in 1:Sim.NAtoms
    hist_def[ceil(Int32, Sim.y_uw[i, 1]),1] += 1.0/TotalMass
end
### defitnion of slab histogram for range 1 & 2
for y_pos in vcat(Sim.y_uw[1:3, 1], Sim.y_uw[62:64, 1])
    hist_def[ceil(Int32,y_pos),2] += 1.0/6.0    ### normalized just this time for comparison
end
for y_pos in vcat(Sim.y_uw[3:6, 1], Sim.y_uw[64:end, 1])
    hist_def[ceil(Int32,y_pos),3] += 1.0/8.0    ### normalized just this time for comparison
end

sigmoid(x) =  1/(1+exp(-x)) # abs(x)<10.0 ? 1/(1+exp(-x)) : 0.0
sigmoid_profile(x, w) = sigmoid(x+w)-sigmoid(x-w)

@testset "slab histograms" begin 
    ### This will emit a warning, which is irrelevant in constructed test case.
    HPSAnalysis.computeCOM!(Sim)
    HPSAnalysis.computeClustersByChainCOM(Sim;Cutoff=2.0)

    data = HPSAnalysis.computeCOMOfLargestCluster(Sim)
    @test (all(data .≈ offset)) ### check whether wether COM of largest Cluster matches intended value

    HPSAnalysis.computeSlabHistogram(Sim, Ranges=Vector{UnitRange{Int32}}([1:3,3:6 ]))
    axis = axes(Sim.SlabHistogramSeries, 1)
    ### test for the true histogram and the two subranges
    for (i,ind) in enumerate([1,5,6])
        hist = sum(Sim.SlabHistogramSeries[:,:,ind], dims=2)/Sim.NSteps
        hist /= sum(hist) ### normalise histogram u/A
        @test (all(hist.≈hist_def[:,i]))
    end


    ### decaying box slab density over time
    values = collect(1.06:-0.01:0.07)
    fill!(Sim.SlabHistogramSeries, 0.05) ### have some dilute phase
    for step in axes(Sim.SlabHistogramSeries,2)
        Sim.SlabHistogramSeries[-10:10, step,1] .= values[step]
    end

    ρ_dense, ρ_dilute, ind_dense, ind_dilute = HPSAnalysis.computeSlabDensities(Sim; Windowlength=50, DiluteWidth=50)

    @test all(ρ_dense .≈ values)
    @test all(ρ_dilute .≈ 0.05)

    ### decaying sigmoidal slab density over time
    values = collect(1.06:-0.01:0.07)
    xaxis = collect(axes(Sim.SlabHistogramSeries,1))

    ### vary the definitions of the surfaces
    for (MaxVal, Surface_fac) in zip([0.9, 0.8], [0.75, 0.5]) ### surface_fac relatively conservative here, such that analytical approximations in test case are valid;
        
        ### varying width feature was removed
        width = 35.0 #.*(1.0.+collect(axes(Sim.SlabHistogramSeries,2))./100)*Surface_fac
        for (i,step) in enumerate(axes(Sim.SlabHistogramSeries,2))
            bla = sigmoid_profile.(xaxis, width).*values[step].+0.05
            Sim.SlabHistogramSeries[:, step,1] = bla
        end
        safety_ =25
        ρ_dense, ρ_dilute, ind_dense, ind_dilute = HPSAnalysis.computeSlabDensities(Sim;DenseWidth=5,MaxVal=MaxVal, Surface_fac=Surface_fac, safety=safety_, Windowlength=50, DiluteWidth=50)

        @test all(ρ_dense .≈ (values.+0.05))
        @test all(ρ_dilute.≈ 0.05)

        ### analytical solution to sigmoid = MaxVal; index_dense is negavtive
        index_dense = -width .- log.(values ./((ρ_dense).*MaxVal.-0.05).-1) 
        index_dense = ceil.(Int32,collect(index_dense))

        ### internal index starts at 1 compared to offset array used here.
        diff = (1-Surface_fac) .*(index_dense.-1)
        index_dense = round.(Int32, (index_dense.-1).*Surface_fac)

        index_dilute = -width .- log.(values./((1.0.-MaxVal).*(ρ_dense.-0.05))).-1
        index_dilute = floor.(Int32, floor.(Int32, index_dilute) ) #-abs.(diff))

        index_dense .= minimum(abs.(index_dense))
        safety=ceil(Int32, safety_/Sim.Resolution)
        index_dilute .= maximum(abs.(index_dilute)) + safety

        @test all(abs.(index_dense) .== ind_dense)
        @test all(abs.(index_dilute) .== ind_dilute)
    end
end






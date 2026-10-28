
### construct SimData for 2 chains of different length inside slab box
Sim = HPSAnalysis.SimData()

Sim.NAtoms=40
Sim.NChains=8
Sim.IDs = ones(Sim.NAtoms)
Sim.Masses = ones(Float32, Sim.NAtoms)
Sim.ChainStart= [1, 6,11,16,21,26,31,36]
Sim.ChainStop = collect(5:5:40)
Sim.ChainMasses = 5*ones(Sim.NChains)
TotalMass= sum(Sim.ChainMasses)
#Sim.Charges = ones(Sim.NAtoms)
Sim.BoxSize=[-10.0 10.0; -10.0 10.0; -10.0 10.0]
Sim.BoxLength=[20.0,20.0,20.0]
Sim.NSteps = 1


offset = collect(1:5)
Sim.x = zeros(Sim.NAtoms, Sim.NSteps)
Sim.y = zeros(Sim.NAtoms, Sim.NSteps)
Sim.z = zeros(Sim.NAtoms, Sim.NSteps)

# chain 1
Sim.x[1:5,1] += offset  

# chain 2 incontact with 1 if cutdist>=1
Sim.x[6:10,1] .+= 1
Sim.y[6:10,1] += offset

# chain 3; incontact with 1-2 if cutdist>=1
Sim.x[11:15,1] += offset
Sim.y[11:15,1] .+= 6

# chain 4 incontact with 1-3 if cutdist>=2
Sim.x[16:20,1] += offset
Sim.y[16:20,1] .+= 8

#chain 5
Sim.x[21:25,1] += -11 .+ offset
Sim.y[21:25,1] .+= -10
Sim.z[21:25,1] .+= -10

#chain 6; incontact with 5  if cutdist>=1
Sim.x[26:30,1] .+=  10
Sim.y[26:30,1] +=  10 .- offset
Sim.z[26:30,1] .+=  10

#chain 7; incontact with 5 & 6 if cutdist>=2
Sim.x[31:35,1] .+=  10
Sim.y[31:35,1] +=  10 .- offset
Sim.z[31:35,1] .+=  -8

#chain 8; always alone
Sim.x[36:40,1] .+=  -5.0
Sim.y[36:40,1] .+=  -5.0
Sim.z[36:40,1] .+=  -5.0


Sim.x_uw = deepcopy(Sim.x)
Sim.y_uw = deepcopy(Sim.y)
Sim.z_uw = deepcopy(Sim.z)

Sim.EquilibrationTime=1
Sim.RGMeasureStep=1


@testset "computeClustersByBeadDistance" begin
    @test Vector{Vector{Vector{Int32}}}([[[1],[2],[3],[4],[5],[6],[7],[8]]])== HPSAnalysis.computeClustersByBeadDistance(Sim; Cutoff=0.8)

    @test Vector{Vector{Vector{Int32}}}([[[1,2,3],[4],[5,6],[7],[8]]]) == HPSAnalysis.computeClustersByBeadDistance(Sim; Cutoff=1.1)

    @test Vector{Vector{Vector{Int32}}}([[[1,2,3,4],[5,6,7],[8]]]) == HPSAnalysis.computeClustersByBeadDistance(Sim; Cutoff=2.1)
end

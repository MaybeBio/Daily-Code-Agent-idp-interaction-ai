using HDF5, StatsBase

@testset "Surface Tension" begin

Sim = HPSAnalysis.SimData()
Sim.SlabAxis=3
Sim.BoxLength = [10,10,10]
Sim.BasePath = TestPath * "SurfaceTensionTest/"

testname ="testfile.h5"
mkpath(TestPath * "SurfaceTensionTest/")
filename=TestPath * "SurfaceTensionTest/$(testname)"
testfile=fid = h5open(filename, "w")

key = "hoomd-data/md/compute/ThermodynamicQuantities/pressure_tensor"
N = 10_202
pressure_tensor = rand(6,N)
testfile[key] = pressure_tensor

timestep = (1:N) .* 100
key = "hoomd-data/Simulation/timestep"
testfile[key] = collect(timestep)


start = 20
NFrames = 1000
Sim.TrajWriteOutFreq =NFrames
τ=2
NSub=5
γ, Δγ = HPSAnalysis.computeSurfaceTension(Sim,[start:NFrames];filename=testname, tau=Int32(τ), NSub=Int32(NSub))

range_ = start*10:τ:NFrames*10

gamma_vec = Sim.BoxLength[Sim.SlabAxis]/2.0 *(pressure_tensor[6,range_].-(pressure_tensor[1,range_].+pressure_tensor[4,range_])./2.0)
gamma_mean = mean(gamma_vec)

gamma_sub = [mean(gamma_vec[a:e]) for (a,e) in [(1, 980),(981, 1960),(1961, 2940),(2941, 3920),(3921, 4900)]]
gamma_error = sqrt(sum((gamma_sub .- gamma_mean).^2)/4)

    @test γ  ≈ gamma_mean
    @test Δγ ≈ gamma_error
end
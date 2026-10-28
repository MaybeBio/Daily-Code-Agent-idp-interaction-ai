using GSDFormat, HDF5, Plots

SetupTestPath="$(TestPath)/Restart_test/"
rm(SetupTestPath; force=true, recursive=true)
mkpath(SetupTestPath)
pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")
sim = pyimport("Submit_HOOMD")
###############


pH=8.0
NChains = 2
Names=["A","B"]
Sequences=["EEEEE", "DDDDD"] 
BoxLengthShort=200.0
BoxLengthLong=2000.0
BoxSize = [-BoxLengthShort/2., BoxLengthShort/2.,-BoxLengthLong/2., BoxLengthLong/2.,-BoxLengthShort/2., BoxLengthShort /2.]
Info="Test Sims"

pos = zeros(Float32, 10, 3)

pos[1:5, 2] .= collect(1:5)*3.8
pos[6:10, 2] .= collect(1:5)*3.8.+50

ENM = (0,[],[],[], Dict())
cd(SetupTestPath)
HPSAnalysis.Setup.writeStartConfiguration(SetupTestPath, "./Test_slab","./Test_Start_slab.txt",Info, Sequences, BoxSize , 10_000, HOOMD=true; SimulationType="Calvados2" , Temperature=300,  InitStyle="Pos", Pos=pos , pH=pH, SaltConcentration=0.2, Device="CPU", WriteOutFreq=1_000, ENM=ENM)

@testset "Restarts" begin
    TrajectoryNumber , NStepsOld = sim.CountNumberOfTrajectoryFiles(SetupTestPath)
    @test TrajectoryNumber == 0 
    @test NStepsOld == 0 

    sim.run(SetupTestPath) ### first simulation runs additionally 8000 steps for to relax energies 

    TrajectoryNumber , NStepsOld = sim.CountNumberOfTrajectoryFiles(SetupTestPath)
    @test TrajectoryNumber == 1
    @test NStepsOld == 18

    sim.restart(SetupTestPath; ExtendedSteps=20_000)
    TrajectoryNumber , NStepsOld = sim.CountNumberOfTrajectoryFiles(SetupTestPath)
    @test TrajectoryNumber == 2
    @test NStepsOld == 20

    sim.restart(SetupTestPath; ExtendedSteps=25_000)
    TrajectoryNumber , NStepsOld = sim.CountNumberOfTrajectoryFiles(SetupTestPath)
    @test TrajectoryNumber == 3
    @test NStepsOld == 25

    sim.restart(SetupTestPath; ExtendedSteps=30_000)
    sim.restart(SetupTestPath; ExtendedSteps=35_000)

    traj_1 = GSDFormat.open("$(SetupTestPath)/traj.gsd")
    traj_2 = GSDFormat.open("$(SetupTestPath)/traj_1.gsd")
    traj_3 = GSDFormat.open("$(SetupTestPath)/traj_2.gsd")
    traj_4 = GSDFormat.open("$(SetupTestPath)/traj_3.gsd")
    traj_5 = GSDFormat.open("$(SetupTestPath)/traj_4.gsd")


    @test traj_1[1].particles.typeid==traj_2[1].particles.typeid
    @test traj_1[1].particles.types ==traj_2[1].particles.types
    @test traj_1[1].particles.mass  ==traj_2[1].particles.mass
    @test traj_1[1].particles.charge==traj_2[1].particles.charge

    @test traj_1[1].bonds.types ==traj_2[1].bonds.types
    @test traj_1[1].bonds.typeid==traj_2[1].bonds.typeid

    @test traj_1[1].angles.types ==traj_2[1].angles.types
    @test traj_1[1].angles.typeid==traj_2[1].angles.typeid

    @test traj_1[1].dihedrals.types ==traj_2[1].dihedrals.types
    @test traj_1[1].dihedrals.typeid==traj_2[1].dihedrals.typeid

    ### test simulation extension by restarts
    data = HPSAnalysis.initData(SetupTestPath; Reparse=true, LoadAll=true, HOOMD=true, ReadBig=true)
    @test data.NSteps == 35

    file = h5open("$(SetupTestPath)/pressure.h5")
    steps = file["hoomd-data"]["Simulation"]["timestep"][:]
    @test length(steps) == 355

    expected_steps = vcat(collect(100:100:8_000), collect(8_000:100:18_000), collect(18_000:100:20_000), collect(20_000:100:25_000), collect(25_000:100:30_000), collect(30_000:100:35_000))

    @test all(steps .== expected_steps)
    close(file)

    ### restart "crashed" simulation
    HPSAnalysis.Setup.WriteParams("$(SetupTestPath)/HOOMD_Setup/Params.txt", "./Test_Start_slab.txt", 300, 40_000, 1_000, 0.01, [20.0,200.0,20.0], 11111; UseAngles=false, Device="CPU")

    sim.restart(SetupTestPath; ExtendedSteps=0) ### runs a new traj for additional 10_000 steps up to 40_000
    file = h5open("$(SetupTestPath)/pressure.h5")
    steps = file["hoomd-data"]["Simulation"]["timestep"][:]
    traj_6 = GSDFormat.open("$(SetupTestPath)/traj_5.gsd")

    @test length(steps) == 406 
    @test all(steps .== vcat(expected_steps, collect(35_000:100:40_000)))
    data = HPSAnalysis.initData(SetupTestPath; Reparse=true, LoadAll=true, HOOMD=true, ReadBig=true)

    ### test if data is correctly read by initData
    alltrue = true ### reduce the number of tests since its basically one test
    for (traj, range) in zip([traj_1, traj_2, traj_3, traj_4, traj_5,traj_6], [1:18, 19:20, 21:25, 26:30, 31:35, 36:40]) 
        for (step, frame) in zip(range,[frame for (i,frame) in enumerate(traj) if i >1])
            alltrue &= all(frame.particles.position[:, 1].*10.0 .≈ data.x[:, step])
            alltrue &= all(frame.particles.position[:, 2].*10.0 .≈ data.y[:, step])
            alltrue &= all(frame.particles.position[:, 3].*10.0 .≈ data.z[:, step])
        end
    end
    @test alltrue

    ### check if last frame of old trajectory is first frame of next trajectroy
    @test all(traj_1[19].particles.position .≈ traj_2[1].particles.position)
    @test all( traj_2[3].particles.position .≈ traj_3[1].particles.position)
    @test all( traj_3[6].particles.position .≈ traj_4[1].particles.position)
    @test all( traj_4[6].particles.position .≈ traj_5[1].particles.position)
    @test all( traj_5[6].particles.position .≈ traj_6[1].particles.position)

    close(traj_1)
    close(traj_2)
    close(traj_3)
    close(traj_4)
    close(traj_5)
    close(traj_6)

    cd(TestPath)
end
SimName="test"

import HPSAnalysis.BioData as BioData
using GSDFormat 

if isdir("$SetupTestPath/HOOMD_Setup/")
    rm("$SetupTestPath/HOOMD_Setup/"; force=true, recursive=true)
end
mkpath("$SetupTestPath/HOOMD_Setup/")


Seq=["DEGHKDEGHK"]
HPSAnalysis.Setup.WriteHOOMDSequences("$SetupTestPath/HOOMD_Setup/Sequences.txt", Seq)

N=10
NChains=1
InputBonds=[[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9]]
InputAngles=[]
InputDihedrals=[]
AtomTypes= Set(join(Seq))
AaToId = Dict{Char,Int32}()
for (index, value) in enumerate(AtomTypes)
    AaToId[value]=index
end


position=[0.0 0.0 1.0 0.5 0.0 0.0 -0.5 -1.0 0.0 0.0;;;-4.0 -3.0 -2.0 -1.0 0.0 0.0 1.0 2.0 3.0 4.0;;;0.0 0.0 -1.0 -0.5 -1.0 1.0 0.5 1.0 0.0 0.0]


coor=fill(0.0,1,maximum(length(seq) for seq in Seq),3)

ChargeDict= BioData.OneToHPSCharge
WeightDict = BioData.AaToWeight
SigmaDict = BioData.OneToHPSCalvadosSigma
LambdaDict = BioData.OneToCalvados3Lambda


HPSAnalysis.Setup.WriteHOOMDParticlesInput("$SetupTestPath/HOOMD_Setup/Particles.txt",position,ChargeDict,AaToId,Seq,WeightDict,SigmaDict,coor)
HPSAnalysis.Setup.WriteDictionaries("$SetupTestPath/HOOMD_Setup/Dictionaries.txt", ChargeDict, AaToId,WeightDict,SigmaDict, LambdaDict)
HPSAnalysis.Setup.WriteDihedrals("$SetupTestPath/HOOMD_Setup/DihedralMap.txt",[],0)
HPSAnalysis.Setup.WriteParams("$SetupTestPath/HOOMD_Setup/Params.txt",SimName,300, 10, 1, 0.01, Array([10,101,10]), rand(1:65535), UseAngles=false,domain=Array([[3,8]]),Device="CPU", UseCharge=false, Create_Start_Config=true,SimType="Calvados3")

harmonic = Dict{String, Dict{Symbol, Float64}}()
r0 = 5
harmonic["O-O"]   = Dict(:r => 0.1, :k => 700)
harmonic["BB_1"]  = Dict(:r => 0.2, :k => 700)
harmonic["ENM_2"] = Dict(:r => 0.3, :k => 1800)
harmonic["ENM_3"] = Dict(:r => 0.4, :k => 1800)
harmonic["ENM_4"] = Dict(:r => 0.5, :k => 1800)

bondid=UInt32[0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 3, 4]

HPSAnalysis.Setup.WriteENM_HOOMD_Indices("$SetupTestPath/HOOMD_Setup/ENM_indices.txt", (12, ["O-O", "BB_1", "ENM_2", "ENM_3", "ENM_4"], bondid, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9), (3, 5), (3, 6), (4, 6)], harmonic))

sim.run(SetupTestPath) ### generates the following file that is the input of the simulation.

data=GSDFormat.open("$(SetupTestPath)$(SimName)_300.0_Start_slab.gsd","r")
frame = data[1]

typesid=UInt32[1, 3, 4, 2, 0, 1, 3, 4, 2, 0]
coor=Float32[0.0 -0.4 0.0; 0.0 -0.3 0.0; 0.1 -0.2 -0.1; 0.05 -0.1 -0.05; 0.0 0.0 -0.1; 0.0 0.0 0.1; -0.05 0.1 0.05; -0.1 0.2 0.1; 0.0 0.3 0.0; 0.0 0.4 0.0]
bond_group=Int32[0 1; 1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 3 5; 3 6; 4 6]
image=Int32[0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0]
types=["K", "D", "H", "E", "G"]
bondtypes=["O-O\0\0", "BB_1\0", "ENM_2", "ENM_3", "ENM_4"]

@testset "Calvados3" begin
    @test frame.particles.N==10
    @test frame.particles.position==coor
    @test frame.particles.types==types
    @test frame.particles.typeid==typesid
    @test frame.particles.image==image
    
    @test frame.bonds.N==UInt32(12)
    @test frame.bonds.types==bondtypes
    @test frame.bonds.typeid==bondid
    @test frame.bonds.group==bond_group
end

close(data)


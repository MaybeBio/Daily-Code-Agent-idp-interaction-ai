### Run different proteins for which there are reverense values for the R_g in the Calvados3 paper, to verify our implementation. As multi domain proteins we use Tia1, Ubq3, Ubq2, Ubq4 and Gal3.
using Distributed
#addprocs(5)

@everywhere using HPSAnalysis
@everywhere using PyCall
#=
BasePath = "$(SetupTestPath)/implementation_test/"
if isdir(BasePath)
    rm(BasePath; force=true, recursive=true)
end
mkpath(BasePath)
=#
@everywhere begin
    BasePath = "/localscratch/test/Calvados3/NinethTry/"

    Proteins = ["tia1", "ubq2", "ubq3", "ubq4", "gal3"]
    DomainDict = Dict("tia1" => [(6, 82),(95, 172),(190, 275)], "ubq2" => [(11, 82),(87, 158)], "ubq3" => [(1, 72),(77, 148),(153, 224)], "ubq4" => [(1, 72),(77, 148),(153, 224),(229, 300)], "gal3" => [(117, 250)])

    PP = "/localscratch/test/"
    PP ="/localscratch/test/Calvados3/AlphaFold3/"

    ProteinJSON = Dict("tia1" => "$(PP)fold_tia1/fold_tia1_full_data_0.json", "ubq2" => "$(PP)fold_ubq2/fold_ubq2_full_data_0.json", "ubq3" => "$(PP)fold_ubq3/fold_ubq3_full_data_0.json", "ubq4" => "$(PP)fold_ubq4/fold_ubq4_full_data_0.json", "gal3" => "$(PP)fold_gal3/fold_gal3_full_data_0.json")

    #ProteinCif  = Dict("tia1" => "$(PP)fold_tia1/fold_tia1_model_0.cif", "ubq2" => "$(PP)fold_ubq2/fold_ubq2_model_0.cif", "ubq3" => "$(PP)fold_ubq3/fold_ubq3_model_0.cif", "ubq4" => "$(PP)fold_ubq4/fold_ubq4_model_0.cif", "gal3" => "$(PP)fold_gal3/fold_gal3_model_0.cif")

    PP = "/localscratch/test/Calvados3/Calcados3_Alphafold2_data/CIF/"
    ProteinCif  = Dict("tia1" => "$(PP)TIA1_rank0_relax.pdb.cif", "ubq2" => "$(PP)Ubq2_rank0_relax.pdb.cif", "ubq3" => "$(PP)Ubq3_rank0_relax.pdb.cif", "ubq4" => "$(PP)Ubq4_rank0_relax.pdb.cif", "gal3" => "$(PP)Gal3_rank0_relax.pdb.cif")

    RgDict = Dict("tia1" => 2.645, "ubq2" => 2.076, "ubq3" => 2.649, "ubq4" => 3.230, "gal3" => 3.026)

    pH = Dict("ubq4" => 8.0, "ubq3" => 8.0, "ubq2" => 8.0, "tia1"=>6.0, "gal3"=> 7.0)
    width_multiplier = 1.0
    Temperatures = Dict("ubq4"=>293.0, "ubq3" => 293.0, "ubq2" => 293.0, "tia1"=>293.15, "gal3"=> 303.0)
    ionic = Dict("ubq4"=>0.33,"ubq3"=>0.33,"ubq2"=>0.33, "tia1"=>0.1, "gal3"=> 0.04)
end


@everywhere function Rg(protein, Path, Sim)
    HPSAnalysis.computeCOM!(Sim)
    HPSAnalysis.computeRGComponentSeries!(Sim)
    HPSAnalysis.computeRGCorrelationTime(Sim)
    #Sim.RGMeasureStep
    #Sim.RGSeries
    #println(Sim.RGSeries)
    #println(Sim.RGMeasureStep)
    mkpath("$Path/Rg/")
    io = open("$Path/Rg/$(protein)_Rg.txt", "w")
    write(io, "$(Sim.RGSeries)")
    close(io)
end

@everywhere function run_sim_prot(protein, BasePath, FoldedDomains, ProteinToJSON, ProteinToCif, pH, width_multiplier, Temperatures, precision=1000)
    PkgSourcePath = "/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-1])
    EnvironmentPath = HPSAnalysis.getPythonEnvironment(PkgSourcePath)
    ENV["PYCALL_JL_RUNTIME_PYTHON"] = "$(EnvironmentPath)/bin/python"
    pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")
    sim  = pyimport("Submit_HOOMD")
    func = pyimport("PythonFuncs")

    pad = "001"
    temp = Temperatures[protein]   

    mkpath(BasePath*"$(protein)")
    mkpath(BasePath*"$(protein)/$(temp)K/")
    mkpath(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/")
    mkpath(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/Restart/")
    Path = BasePath*"$(protein)/$(temp)K/RUN_$(pad)/"
    
    println(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/")
    Seq = HPSAnalysis.ProteinSequences.NameToSeq[protein]
    NChains = 1
    Sequences = [deepcopy(Seq) for _ in 1:NChains]
    Proteins  = [deepcopy(protein) for _ in 1:NChains]

    Info = "SLAB Simulation script for $protein.\n\n"
    BoxLength = Float32(250.0)      
    BoxSize = [-BoxLength/2., BoxLength/2.,-BoxLength/2., BoxLength/2.,-BoxLength/2., BoxLength /2.]

    SimulName = "$(protein)_$temp"

    (pos, Data) = HPSAnalysis.Setup.CreateStartConfiguration(SimulName,Path , Float32.([BoxLength,BoxLength*width_multiplier , BoxLength]), Proteins, Sequences, Regenerate=true; Axis="y", SimulationType="Calvados3",ProteinToDomain=FoldedDomains,ProteinToCif=ProteinToCif)

    ENM = HPSAnalysis.Setup.BuildENMModel(Data, FoldedDomains, Proteins, Sequences, ProteinToJSON; plDDTcut=0, pae_cut=10000.0, precision=precision) ### turnoff cut offs

    HPSAnalysis.Setup.writeStartConfiguration(Path, "/$(protein)_slab","/$(SimulName)_Start_slab", Info, Sequences, BoxSize , 100_000_000, HOOMD=true ; SimulationType="Calvados3" , Temperature=temp,  InitStyle="Pos", Pos=pos, SaltConcentration=ionic[protein], pH=pH[protein],domain=FoldedDomains,Device="CPU",WriteOutFreq=100_000, ENM, ah_cut=22)


    #=
    sim.run("$(Path)/")
    _, _, _, IDToMass, _, _ = func.readDictionaries("$(Path)/HOOMD_Setup/Dictionaries.txt")
    Masses = func.compute_Mass_List(Data.IDs, IDToMass)
    Data.Masses = Masses
    Data.ChainMasses = [sum(Masses)]
    Rg(protein, Path, Data)
    =#
end

@everywhere function runStuff(protein)
    run_sim_prot(protein, BasePath, DomainDict, ProteinJSON, ProteinCif, pH, width_multiplier, Temperatures)
end

#pmap(runStuff, Proteins)

runStuff(Proteins[1])
runStuff(Proteins[2])
runStuff(Proteins[3])
runStuff(Proteins[4])
runStuff(Proteins[5])

#run_sim_prot("tia1", BasePath, DomainDict, ProteinJSON, ProteinCif, pH, width_multiplier, Temperatures)

using HPSAnalysis
using PyCall

BasePath = "$(SetupTestPath)/Barostat_test/"
#if isdir(BasePath)
#    rm(BasePath; force=true, recursive=true)
#end
mkpath(BasePath)

testcif = "_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1    N N   . MET A 1 1   ? -23.095 28.656  -30.598 1.00 35.20 1   A 1 
ATOM 2    C CA  . MET A 1 1   ? -23.073 27.290  -30.037 1.00 45.41 1   A 1 
ATOM 3    C C   . MET A 1 1   ? -22.970 27.444  -28.533 1.00 47.08 1   A 1 
ATOM 4    O O   . MET A 1 1   ? -21.966 27.964  -28.073 1.00 44.37 1   A 1 
ATOM 5    C CB  . MET A 1 1   ? -21.887 26.489  -30.586 1.00 43.95 1   A 1 
ATOM 6    C CG  . MET A 1 1   ? -22.103 26.073  -32.040 1.00 41.34 1   A 1 
ATOM 7    S SD  . MET A 1 1   ? -20.636 25.331  -32.799 1.00 36.61 1   A 1 
ATOM 8    C CE  . MET A 1 1   ? -21.326 24.783  -34.373 1.00 30.08 1   A 1 
ATOM 9    N N   . ALA A 1 2   ? -24.039 27.132  -27.796 1.00 41.37 2   A 1 
ATOM 10   C CA  . ALA A 1 2   ? -24.013 27.248  -26.348 1.00 48.39 2   A 1 
ATOM 11   C C   . ALA A 1 2   ? -23.078 26.152  -25.827 1.00 49.29 2   A 1 
ATOM 12   O O   . ALA A 1 2   ? -23.384 24.969  -25.945 1.00 46.27 2   A 1 
ATOM 13   C CB  . ALA A 1 2   ? -25.442 27.123  -25.804 1.00 46.75 2   A 1 
ATOM 14   N N   . SER A 1 3   ? -21.920 26.553  -25.328 1.00 42.07 3   A 1 
ATOM 15   C CA  . SER A 1 3   ? -21.064 25.649  -24.585 1.00 47.53 3   A 1 "

io = open("$(BasePath)test.cif", "w")
write(io, testcif)
close(io)
TestDict = Dict("test"=> "$(BasePath)test.cif")
PosDict, LenDict = HPSAnalysis.AlphaFold_startpos(TestDict, ["test"],["abc"])

pos = Float32.([-23.073 27.290  -30.037; -24.013 27.248  -26.348;-21.064 25.649  -24.585 ])
shift =  ((maximum(pos, dims=1)+minimum(pos, dims=1))/2.0)[1:3]
pos[:,1] .-= shift[1]
pos[:,2] .-= shift[2]
pos[:,3] .-= shift[3]

@test LenDict == Dict{String, Int32}("test"=>Int32(3))
@test collect(keys(PosDict)) == ["test"]
@test all(isapprox.(PosDict["test"], pos, rtol=10^-3))

(min_vals, max_vals) = HPSAnalysis.getLargestBoundingbox(PosDict)
@test isapprox(min_vals , [-24.013, 25.649, -30.037].-shift, rtol=10^-3)
@test isapprox(max_vals , [-21.064, 27.290, -24.585].-shift, rtol=10^-3)

### potentially test CreateStartConfiguration 


if PythonTests
    Proteins = ["tia1", "ubq2", "ubq3", "ubq4", "gal3"]
    DomainDict = Dict("tia1" => [(6, 82),(95, 172),(190, 275)], "ubq2" => [(11, 82),(87, 158)], "ubq3" => [(1, 72),(77, 148),(153, 224)], "ubq4" => [(1, 72),(77, 148),(153, 224),(229, 300)], "gal3" => [(117, 250)])

    PP ="/localscratch/test/Calvados3/AlphaFold3/"

    ProteinJSON = Dict("tia1" => "$(PP)fold_tia1/fold_tia1_full_data_0.json", "ubq2" => "$(PP)fold_ubq2/fold_ubq2_full_data_0.json", "ubq3" => "$(PP)fold_ubq3/fold_ubq3_full_data_0.json", "ubq4" => "$(PP)fold_ubq4/fold_ubq4_full_data_0.json", "gal3" => "$(PP)fold_gal3/fold_gal3_full_data_0.json")

    ProteinCif  = Dict("tia1" => "$(PP)fold_tia1/fold_tia1_model_0.cif", "ubq2" => "$(PP)fold_ubq2/fold_ubq2_model_0.cif", "ubq3" => "$(PP)fold_ubq3/fold_ubq3_model_0.cif", "ubq4" => "$(PP)fold_ubq4/fold_ubq4_model_0.cif", "gal3" => "$(PP)fold_gal3/fold_gal3_model_0.cif")


    pH = 7.0
    width_multiplier = 1.0
    Temperatures = 300.0
    function run_sim_prot(protein, BasePath, FoldedDomains, ProteinToJSON, ProteinToCif, pH, width_multiplier, Temperatures)
        PkgSourcePath = "/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-1])
        EnvironmentPath = HPSAnalysis.getPythonEnvironment(PkgSourcePath)
        ENV["PYCALL_JL_RUNTIME_PYTHON"] = "$(EnvironmentPath)/bin/python"
        pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")
        sim = pyimport("Submit_HOOMD")

        RunsPerProtein = 1
        temp = Temperatures

        pad = "001"
        mkpath(BasePath*"$(protein)")
        mkpath(BasePath*"$(protein)/$(temp)K/")
        mkpath(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/")
        mkpath(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/Restart/")
        Path = BasePath*"$(protein)/$(temp)K/RUN_$(pad)/"
        
        Seq = HPSAnalysis.ProteinSequences.NameToSeq[protein]
        NChains = 150
        Sequences = [deepcopy(Seq) for _ in 1:NChains]
        Proteins  = [deepcopy(protein) for _ in 1:NChains]

        Info = "SLAB Simulation script for $protein.\n\n"
        BoxLengthShort = Float32(350.0)      
        BoxLengthLong  = BoxLengthShort*width_multiplier
        BoxSize = [-BoxLengthShort/2., BoxLengthShort/2.,-BoxLengthLong/2., BoxLengthLong/2.,-BoxLengthShort/2., BoxLengthShort /2.]

        SimulName = "$(protein)_$temp"

        (pos, Data) = HPSAnalysis.CreateStartConfiguration(SimulName,Path , Float32.([BoxLengthShort,BoxLengthLong , BoxLengthShort]), Proteins, Sequences, Regenerate=true; Axis="y", SimulationType="Calvados3",ProteinToDomain=FoldedDomains,ProteinToCif=ProteinCif)
        
        
        HPSAnalysis.RewriteCifToPDB(Path, ProteinToCif, Proteins)

        ENM = HPSAnalysis.Setup.BuildENMModel(Data, FoldedDomains, Proteins, Sequences, ProteinToJSON)

        HPSAnalysis.Setup.writeStartConfiguration(Path, "/$(protein)_slab","/$(SimulName)_prerun_Start_slab", Info, Sequences, BoxSize , 100_000, HOOMD=true ; SimulationType="Calvados3" , Temperature=temp,  InitStyle="Pos", Pos=pos , pH=pH,domain=FoldedDomains,Device="GPU",WriteOutFreq=1_000, ENM, SlabAxis=Data.SlabAxis)

        #=
        sim.preRun("$(Path)/", prerunSteps=100_000)
        @info "End Prerun"

        ### side step to do unwrap because old box -> do preRun as new start conficuration with unwrap in old box
        barostat_pos = HPSAnalysis.readXYZ!(Data, TrajectoryFile="$Path/prerun_traj.gsd")

        
        Sequences = [deepcopy(Seq) for _ in 1:NChains] #-> somthing has chainged the first and last letter of the first Sequenze to a and b
        
        HPSAnalysis.Setup.writeStartConfiguration(Path, "/$(protein)_slab","/$(SimulName)_Start_slab", Info, Sequences, BoxSize , 10_000, HOOMD=true ; SimulationType="Calvados3" , Temperature=temp,  InitStyle="Pos", Pos=barostat_pos , pH=pH,domain=FoldedDomains,Device="GPU",WriteOutFreq=1_000, ENM, SlabAxis=Data.SlabAxis)

        sim.run("$(Path)/")
        =#
    end

    run_sim_prot(Proteins[1], BasePath, DomainDict, ProteinJSON, ProteinCif, 7.0,0.8,300)

end
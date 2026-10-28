using HPSAnalysis

PkgSourcePath="/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-1])
EnvironmentPath= HPSAnalysis.getPythonEnvironment(PkgSourcePath)
ENV["PYCALL_JL_RUNTIME_PYTHON"]="$(EnvironmentPath)/bin/python"

using PyCall
BasePath="$SetupTestPath/RS_Prot/"


rm(BasePath; force=true, recursive=true)
mkpath(BasePath)
mkpath("$BasePath/HOOMD_Setup/")


ToCreate =  ["RS31"]
FoldedDomains = Dict("RS31" => [(1,70),(90,155)])
ProteinToJSON= Dict("RS31" =>"$(PkgPath)/data/TestData/fold_rs31_full_data_0.json","RS31a" =>"$(PkgPath)/data/TestData/fold_rs31a_full_data_0.json" )
ProteinToCif= Dict("RS31" =>"$(PkgPath)/data/TestData/fold_rs31_model_0.cif","RS31a" =>"$(PkgPath)/data/TestData/fold_rs31a_model_0.cif" )

Temperatures=300
pH=7.0

if PythonTests
    pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")
    sim = pyimport("Submit_HOOMD")
end

### implement calvados 3 parameters b hand
### taken from https://github.com/KULL-Centre/CALVADOS/blob/main/examples/slab_IDR_MDP/input/residues_CALVADOS3.csv
Lambda= Dict("R"=>0.7407902764839954, "D"=>0.092587557536158, "N"=>0.3706962163690402, "E"=>0.000249590539426, "K"=>0.1380602542039267, "H"=> 0.4087176216525476, "Q"=> 0.3143449791669133, "S"=>0.4473142572693176	, "C"=>0.5922529084601322, "G"=> 0.7538308115197386	, "T"=> 0.2672387936544146, "A"=>0.3377244362031627, "M"=>0.5170874160398543, "Y"=>0.950628687301107, "V"=>0.2936174211771383,"W"=>1.033450123574512, "L"=>0.5548615312993875, "I"=>0.5130398874425708, "P"=>0.3469777523519372, "F"=>0.8906449355499866)
Sigma= Dict("R"=>0.656, "D"=>0.558, "N"=>0.568, "E"=>0.592, "K"=>0.636, "H"=>0.608, "Q"=>0.602, "S"=>0.518, "C"=>0.548, "G"=>0.45, "T"=> 0.562, "A"=>0.504, "M"=>0.618, "Y"=>0.646	, "V"=>0.586,"W"=>0.678, "L"=>0.618, "I"=>0.618, "P"=>0.556, "F"=>0.636)
Charge= Dict{String, Float64}("R"=>1, "D"=>-1, "N"=>0, "E"=>-1, "K"=>1, "H"=>0, "Q"=>0, "S"=>0, "C"=>0, "G"=>0, "T"=> 0, "A"=>0, "M"=>0, "Y"=>0, "V"=>0,"W"=>0,"L"=>0, "I"=>0, "P"=>0, "F"=>0)
Mass= Dict("R"=>156.19,"D"=>115.09,"N"=>114.1,"E"=>129.11,"K"=>128.17,"H"=>137.14,"Q"=>128.13,"S"=>87.08,"C"=>103.14,"G"=>57.05,"T"=>101.11,"A"=>71.07,"M"=>131.2,"Y"=>163.18,"V"=>99.13,"W"=>186.22,"L"=>113.16,"I"=>113.16,"P"=>97.12,"F"=>147.18)

Charge["H"] = 1. / ( 1 + 10^(pH-6) )

function readPositionFromCif(filename)
    x_arr = []
    y_arr = []
    z_arr = []

    for line in readlines(filename)
        fields = strip.(split(line))
        if !isempty(fields) && fields[1] == "ATOM"
            ID = parse(Int, fields[2])
            symbole = fields[3]
            label_atom = fields[4]
            if label_atom =="CA"
                push!(x_arr, parse(Float64, fields[11]))
                push!(y_arr, parse(Float64, fields[12]))
                push!(z_arr, parse(Float64, fields[13]))
            end
        end
    end
    xyz = zeros(length(x_arr),3)
    xyz[:,1] .= x_arr
    xyz[:,2] .= y_arr
    xyz[:,3] .= z_arr
    return xyz
end

function parseDictionary(filename)
    Charge = Dict()
    Mass   = Dict()
    Lambda = Dict()
    Sigma  = Dict()

    lines = readlines(filename)
    for line in lines[2:end]
        fields = split(line,",")
        residue = strip(fields[2])
        Charge[residue] = parse(Float64, fields[3])
        Mass[residue]   = parse(Float64, fields[4])
        Sigma[residue]  = parse(Float64, fields[5])/10.0
        Lambda[residue] = parse(Float64, fields[6])
    end
    return Charge, Mass, Sigma, Lambda
end

### Test the full setup routine for an actual protein
for (protID, protein) in enumerate(ToCreate)
    mkpath(BasePath*"$(protein)/")
    for temp in Temperatures
        pad = "001"
        mkpath(BasePath*"$(protein)/$(temp)K/")
        mkpath(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/")
        mkpath(BasePath*"$(protein)/$(temp)K/RUN_$(pad)/Restart/")
        Path = BasePath*"$(protein)/$(temp)K/RUN_$(pad)/"
        cd(Path)

        Seque = HPSAnalysis.ProteinSequences.NameToSeq[protein]
        NChain=1

        Seque = HPSAnalysis.ProteinSequences.NameToSeq[protein]
        local Sequences= [deepcopy(Seque) for _ in 1:NChain]
        Proteins = [deepcopy(protein) for _ in 1:NChain]

        ###FoldedDomain -> NChain * FoldedDomain

        local Info ="SLAB Simulation script for $protein.\n\n"
        BoxLS=Float32(350.0)
        BoxLL=Float32(1500.)
        local BoxSize = [-BoxLS/2., BoxLS/2.,-BoxLL/2., BoxLL/2.,-BoxLS/2., BoxLS /2.]

        SimulName = "$(protein)_$temp"

        (_, Data) = HPSAnalysis.Setup.CreateStartConfiguration(SimulName,Path , Float32.([BoxLS,BoxLS , BoxLS]), Proteins, Sequences, Regenerate=true; Axis="y", SimulationType="Calvados3",ProteinToDomain=FoldedDomains[protein],ProteinToCif=ProteinToCif)

        pos = readPositionFromCif(ProteinToCif["RS31"])
        local ENM = HPSAnalysis.Setup.BuildENMModel(Data, FoldedDomains, Proteins, Sequences, ProteinToJSON)

        HPSAnalysis.Setup.writeStartConfiguration(Path, "/$(protein)_slab","/$(SimulName)_Start_slab", Info, Sequences, BoxSize , 1, HOOMD=true ; SimulationType="Calvados3" , Temperature=temp,  InitStyle="Pos", Pos=pos , pH=pH,domain=Dict(protein=>FoldedDomains[protein]),Device="CPU",WriteOutFreq=100, ENM)

        ### test if it crashes
        if PythonTests
            sim.run("$(Path)/")
        end

        Charge_test, Mass_test, Sigma_test, Lambda_test = parseDictionary("$(Path)/HOOMD_Setup/Dictionaries.txt")

        ### Add N and C termini
        first = "$(Seque[1])"
        Charge["a"] = Charge[first] +1.0
        Mass["a"]   = Mass[first] +2
        Sigma["a"]  = Sigma[first]
        Lambda["a"] = Lambda[first]

        last = "$(Seque[end])"
        Charge["b"] = Charge[last] - 1.0
        Mass["b"]   = Mass[last] +16.0
        Sigma["b"]  = Sigma[last]
        Lambda["b"] = Lambda[last]

        @testset "Calvados3 Param" begin
            @test all(map(x-> Charge[x]≈Charge_test[x], collect(keys(Charge_test))))
            @test all(map(x-> Sigma[x] ≈Sigma_test[x] , collect(keys(Sigma_test))))
            @test all(map(x-> Lambda[x]≈Lambda_test[x], collect(keys(Lambda_test))))
            @test all(map(x-> Charge[x]≈Charge_test[x], collect(keys(Charge_test))))
            @test all(map(x->isapprox(Mass[x],Mass_test[x],atol=0.01) , collect(keys(Mass_test)))) ## we have one more digit for the masses

        end
    end
end

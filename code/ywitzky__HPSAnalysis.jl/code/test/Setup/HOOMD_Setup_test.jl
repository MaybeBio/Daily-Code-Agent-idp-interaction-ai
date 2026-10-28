using HPSAnalysis,Test,Printf, PyCall

filename_test="$SetupTestPath/HOOMD_write_test.csv"
filename="$SetupTestPath/HOOMD_write.csv"

if PythonTests
    ### define python based conversion from numpy to pure python
    py"""import numpy
    def castToInt(x):
        return int(x)"""
end

@testset "WriteHOOMDSequences" begin
    filename_test="$SetupTestPath/HOOMD_sequence_test.txt"
    filename="$SetupTestPath/sequence.txt"
    sequences=["MRPVFV","MRPVF","MRPV","MRP"]
    io=open(filename_test,"w")
    for (_,seq) in enumerate(sequences)
        write(io,"$seq\n")
    end
    close(io)
    
    HPSAnalysis.Setup.WriteHOOMDSequences(filename,sequences)
    @test files_are_equal(filename_test,filename)

    if PythonTests
        Seqs, NBeads, NChains, InputBonds, InputAngles, InputDihedrals = sim.readSequences(filename)

        Bonds = []; Angles=[]; Dihedrals=[];
        cb =0;    ca =0;    cd =0;    cnt = 0
        for Seq in sequences
            append!(Bonds,    [(cnt+i-1,cnt+i) for (i, _) in enumerate(Seq[1:end-1])])
            append!(Angles,   [(cnt+i-1,cnt+i, cnt+i+1) for (i, _) in enumerate(Seq[1:end-2])])
            append!(Dihedrals,[(cnt+i-1,cnt+i, cnt+i+1, cnt+i+2) for (i, _) in enumerate(Seq[1:end-3])])
            cnt += length(Seq)
        end

        @test Seqs==sequences
        @test NBeads ==sum(length.(sequences))
        @test NChains ==length(sequences)
        @test permutedims(hcat(collect.(Bonds)...)) == InputBonds
        @test permutedims(hcat(collect.(Angles)...)) == InputAngles
        @test permutedims(hcat(collect.(Dihedrals)...)) == InputDihedrals
    end
end

@testset "WriteDictionaries" begin
    io=open(filename_test,"w")
    ToCharge_test=Dict('A'=>1.0,'B'=>-1.0,'C'=>2.0,'D'=>-2.0)
    ToID_test=Dict('A'=>1,'B'=>2,'C'=>3,'D'=>4)
    ToMass_test=Dict('A'=>1.0,'B'=>2.0,'C'=>3.0,'D'=>4.0)
    ToDiameter_test=Dict('A'=>1.0,'B'=>2.0,'C'=>3.0,'D'=>4.0)
    ToLambda_test=Dict('A'=>1.0,'B'=>2.0,'C'=>3.0,'D'=>4.0)

    write(io, "// ID,resname, Charge, Mass, λ   \n")
    keys=['A','B','C','D']
    for key in keys
        write(io, " $(ToID_test[key]), $key, $(ToCharge_test[key]), $(ToMass_test[key]), $(ToDiameter_test[key]), $(ToLambda_test[key])\n")
    end

    close(io)
    HPSAnalysis.Setup.WriteDictionaries(filename,ToCharge_test,ToID_test,ToMass_test,ToDiameter_test,ToLambda_test)
    @test files_are_equal(filename_test,filename)

    if PythonTests
        ### now test the python read scripts
        ID, IDToResName, IDToCharge, IDToMass, IDToSigma, IDToLambda = sim.readDictionaries(filename)

        ### pyCall doesnt autoconvert numpy ints if part of a dictionary
        ### additionally switch from 0 index to 1 indexing
        IDToResName = Dict((py"castToInt"(key)+1, Char(value[1])) for (key,value) in IDToResName)
        IDToCharge  = Dict((py"castToInt"(key)+1, value         ) for (key,value) in IDToCharge)
        IDToMass    = Dict((py"castToInt"(key)+1, value         ) for (key,value) in IDToMass)
        IDToSigma   = Dict((py"castToInt"(key)+1, value         ) for (key,value) in IDToSigma)
        IDToLambda  = Dict((py"castToInt"(key)+1, value         ) for (key,value) in IDToLambda)

        @test [ToID_test[key] for key in keys] == ID.+1 ### python routine uses C indices
        @test Dict(ToID_test[key]=>key                  for key in keys) == Dict{Int64, Char}(IDToResName)
        @test Dict(ToID_test[key]=>ToCharge_test[key]   for key in keys) == Dict{Int64, Float64}(IDToCharge)  
        @test Dict(ToID_test[key]=>ToMass_test[key]     for key in keys) == IDToMass 
        @test Dict(ToID_test[key]=>ToDiameter_test[key] for key in keys) == IDToSigma 
        @test Dict(ToID_test[key]=>ToLambda_test[key]   for key in keys) == IDToLambda 
    end
end

if PythonTests
    @testset "WriteHOOMDParticlesInput" begin
        filename_test="$SetupTestPath/particles_test.csv"
        filename="$SetupTestPath/particles.csv"

        sequences=["MRPVFV","MRPVF","MRPV","MRP"]
        io=open(filename_test,"w")
        set=Set(join(sequences))
        coor=fill(0.5,4,maximum(length(seq) for seq in sequences),3)
        ToCharge=Dict(atom=>num*0.1 for (num, atom) in enumerate(set))
        ToID=Dict(atom=>num for (num, atom) in enumerate(set))
        ToMass=Dict(atom=>50+num for (num, atom) in enumerate(set))
        ToDiameter=Dict(atom=>10.0+num for (num, atom) in enumerate(set))
        pos=coor
        image=coor

        write(io, "### N, id, x ,y ,z , charge, m, diameter \n");
        cnt=1
        for (SeqID, Seq) in enumerate(sequences)
            for (atom,res) in enumerate(Seq)
                write(io, "$(cnt), $(ToID[res] -1) , $(@sprintf("%.3f",pos[SeqID,atom, 1])), $(@sprintf("%.3f",pos[SeqID,atom, 2])), $(@sprintf("%.3f",pos[SeqID,atom, 3])), $(ToCharge[res]), $(ToMass[res]), $(ToDiameter[res]) , $(image[SeqID, atom,1]) , $(image[SeqID, atom,2]) , $(image[SeqID, atom,3]) \n");
                cnt += 1
            end
        end
        close(io)
        
        HPSAnalysis.Setup.WriteHOOMDParticlesInput(filename, pos, ToCharge, ToID, sequences, ToMass, ToDiameter, image)
        @test files_are_equal(filename_test,filename)

        InputPositions, InputTypes, InputCharges, InputMasses, Types, Diameter, InputImage = sim.readParticleData(filename, 18, sequences)

        @test [ToID[res]-1 for res in prod(sequences)] == InputTypes
        @test [ToCharge[res] for res in prod(sequences)] == InputCharges
        @test [ToMass[res] for res in prod(sequences)] == InputMasses
        @test [ToDiameter[res] for res in prod(sequences)] == Diameter .*10.0 ### convert from nm to AA
    end

    @testset "Write/read-Params" begin
        filename_test="$SetupTestPath/params_test.csv"
        filename="$SetupTestPath/params.csv"
        sequences=["MRPVFV","MRPVF","MRPV","MRP"]
        io=open(filename_test,"w")
        epsilon_r = 1.73136
        Params=Dict("Simname"=>"C2", "Temp"=>300.0, "NSteps"=>10000, "NOut"=>10, "dt"=>1000.0,
        "Lx"=>10, "Ly"=>100, "Lz"=>11, "Seed"=>10, "Minimise"=>true, "Trajectory"=>"testtraj.gsd", "UseAngles"=>true, "UseCharge"=>true, "Alt_GSD_Start"=>"-", "Create_Start_Config"=>false, "epsilon_r"=>epsilon_r, "kappa"=>1.0, "Device"=>"GPU", "YukawaCutoff"=>4.0, "AHCutoff"=>2.0, "ionic"=>0.1, "pH"=>7.0, "SimulationType"=>"Calvados2","Domains"=>Array([[0,0]]), "yk_prefactor"=>138.9315360433804/epsilon_r , "Axis"=>"2")
        
        for (key, value) in Params
            write(io, "$(key): $(value)\n")
        end
        close(io)
        
        HPSAnalysis.Setup.WriteParams(filename, Params["Simname"], Params["Temp"], Params["NSteps"], Params["NOut"], Params["dt"], [Params["Lx"], Params["Ly"], Params["Lz"]], Params["Seed"]; Minimise=Params["Minimise"], TrajectoryName=Params["Trajectory"], UseAngles=Params["UseAngles"], UseCharge=Params["UseCharge"], Alt_GSD_Start=Params["Alt_GSD_Start"], Create_Start_Config=Params["Create_Start_Config"], ϵ_r=Params["epsilon_r"], κ=Params["kappa"], Device=Params["Device"], yk_cut=Params["YukawaCutoff"], ah_cut=Params["AHCutoff"], ionic=Params["ionic"], pH=Params["pH"], SimType=Params["SimulationType"], domain=Params["Domains"])

        #@test files_are_equal(filename_test,filename) ### new method doesnt write everything in the correct order

        ParamDict = sim.readParam(filename)
        Params["Use_Minimised_GSD"] = ParamDict["Alt_GSD_Start"]=="-" && ParamDict["Minimise"] ? true : false
        Params["Domains"] = repr(Params["Domains"]) ### convert to String

        @test ParamDict == Params
    end

    @testset "Write/read-ENM_HOOMD_Indices" begin
        filename_test = "$SetupTestPath/enmhoomd_test.csv"
        filename = "$SetupTestPath/enmhoomd.csv"
        N = 4
        types = ["B1","B2","B3","B4"]
        id = [0,1,2,3] ### use c indexing
        group = [(0, 1), (1, 2), (2, 3), (3, 4)]
        harmonic = Dict{String, Dict{Symbol, Float64}}()
        harmonic["B1"] = Dict(:r => 0.5, :k => 700)
        harmonic["B2"] = Dict(:r => 0.5, :k => 700)
        harmonic["B3"] = Dict(:r => 0.5, :k => 700)
        harmonic["B4"] = Dict(:r => 0.5, :k => 700)
        ENM = (N, types, id, group, harmonic)

        open(filename_test, "w") do io
            write(io, "1 , B1 , 0 , (0, 1) , $(harmonic["B1"]) \n")
            write(io, "2 , B2 , 1 , (1, 2) , $(harmonic["B2"]) \n")
            write(io, "3 , B3 , 2 , (2, 3) , $(harmonic["B3"]) \n")
            write(io, "4 , B4 , 3 , (3, 4) , $(harmonic["B4"]) \n")
        end

        HPSAnalysis.Setup.WriteENM_HOOMD_Indices(filename, ENM)

        read_ENM = sim.read_ENM_HOOD_indices(filename)
        read_ENM_test = sim.read_ENM_HOOD_indices(filename_test)

        @test N == read_ENM[1]
        @test all(types == read_ENM[2])
        @test id == read_ENM[3]
        @test group == read_ENM[4]
        @test all(map(key-> (harmonic[key][:k]== read_ENM[5][key]["k"])&& (harmonic[key][:r]≈read_ENM[5][key]["r"]), collect(keys(harmonic))))

        @test all(read_ENM[1] == read_ENM_test[1])
        @test all(read_ENM[2] == read_ENM_test[2])
        @test all(read_ENM[3] == read_ENM_test[3])
        @test all(read_ENM[4] == read_ENM_test[4])
        @test all(read_ENM[5] == read_ENM_test[5])
    end
end

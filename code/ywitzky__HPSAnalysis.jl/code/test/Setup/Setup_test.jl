SetupTestPath="$(TestPath)/Setup/"
mkpath(SetupTestPath)

if PythonTests
    pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")
    sim = pyimport("Submit_HOOMD")
end
mkpath("$SetupTestPath/HOOMD_Setup/")

function files_are_equal(file1::String, file2::String)::Bool
    return read(file1, String) == read(file2, String)
end

@testset "Setup Tests" begin
    include("./HOOMD_Setup_test.jl")
    include("./GSD_test.jl")
    if PythonTests
        include("./Restart_test.jl")
    end
end

@testset "Calvados Tests" begin
    include("./Calvados/Calvados_Helper_test.jl")
    if PythonTests
        include("./Calvados/Calvados2_test.jl")
        include("./Calvados/Calvados3_test.jl")
    end
    include("./Calvados/Calvados3_ENM_test.jl")
    include("./Calvados/RS_Prot.jl") 
  
    include("./Calvados/Prerun_resize_test.jl")

    #include("./Calvados/Implementation/Calvados3_implementation_test.jl")
    #include("./Calvados/Implementation/Calvados3_implementation_analysis.jl")
end

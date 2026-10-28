using HPSAnalysis
PkgSourcePath="/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-1])
PkgPath="/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-2])

using Test, Scratch, Aqua

PythonTests = false ### turn off for Github Actions to work; python with environment not available online
if PythonTests
    EnvironmentPath= HPSAnalysis.getPythonEnvironment(PkgSourcePath)
    ENV["PYCALL_JL_RUNTIME_PYTHON"]="$(EnvironmentPath)/bin/python3"
    # might need to run the following lines once
    #ENV["PYTHON"]="$(EnvironmentPath)/bin/python3"
    # using Pkg
    # Pkg.build("PyCall")

    using PyCall
    pushfirst!(pyimport("sys")."path", "$(PkgSourcePath)/Setup/")
end

TestPath = Scratch.get_scratch!(HPSAnalysis, "test") 

@testset "Aqua" begin
    Aqua.test_all(HPSAnalysis; deps_compat=(ignore=[:Printf, :Mmap, :Libdl, :LinearAlgebra, :Statistics,:Test, :Distributed, :Random], ), project_extras=false, )
    ### ignore standard libraries, not sure how to deal with them im PackageCompatUI/add compats manually
    ### Test fails in project extras since they dont get excluded their normaly
end

include("./Analysis/Analysis_test.jl")
include("./Setup/Setup_test.jl")

module HPSAnalysis

include("./IO/PythonEnv.jl")
PkgSourcePath="/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-1])
PkgPath="/"*joinpath(split(pathof(HPSAnalysis),"/")[1:end-2])
EnvironmentPath= HPSAnalysis.getPythonEnvironment(PkgSourcePath)


### helping stuff which probably should be seperate
include("../data/BioData.jl")
include("../data/ProteinSequences.jl")
include("./Data.jl")

include("./Setup.jl")

 
include("./IO.jl")
include("./Analysis.jl")

### load empty functions to be extended in HPSAnalysisPlots
include("../ext/Plot/EmptyFuncs.jl")

end # module HPSAnalysis

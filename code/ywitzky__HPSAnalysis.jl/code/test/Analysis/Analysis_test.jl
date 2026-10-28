
@testset "Analysis Tests" begin
    include("./Slab_Analysis_test.jl")
    include("./Clustering_test.jl")
    include("./Slab_Cumulant_test.jl")
    include("./IntraChainScaling_test.jl")
    include("./IntraChainContactMatrix_test.jl")
    include("./InterChainContactMatrix_test.jl")
    include("./SurfaceTension_test.jl")
    include("./MeanPairEnergyMatrix_test.jl")
    include("./MeanPairIntraEnergyMatrix_test.jl")### needs to be after MeanPairEnergyMatrix
end
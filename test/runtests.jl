using PhylogeneticFactorAnalysis
using Test, SafeTestsets

@time @safetestset "PFA test" begin include("test.jl") end
@time @safetestset "Post processing test" begin include("test_postprocessing.jl") end

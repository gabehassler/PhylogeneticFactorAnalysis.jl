using PhylogeneticFactorAnalysis
using Test, SafeTestsets

@time @safetestset "First test" begin include("test.jl") end

using PhylogeneticFactorAnalysis

using BEASTXMLConstructor, Random, BEASTTreeUtils, CSV, DataFrames

Random.seed!(666)





function throw_error_if_exists(path::String)
    if isfile(path)
        error("File $path already exists. Delete and re-run.")
    end
end

data_path = joinpath(@__DIR__, "test_data.csv")
newick_path = joinpath(@__DIR__, "test_newick.txt")

dummy_files = [data_path, newick_path]

throw_error_if_exists.(dummy_files)



N = 100
K = 5
P = 20

Y = randn(N, P)
taxa = ["taxon$i" for i = 1:N]

data = DataFrame(taxon = taxa)
for i = 1:P
    data[!, "trait$i"] = Y[:, i]
end

CSV.write(data_path, data)

write(newick_path, writeTopology(rtree(taxa, ultrametric=true)))

msp = PhylogeneticFactorAnalysis.ModelSelectionProvider([1, 2, 3], Float64[], 5)
prior = PhylogeneticFactorAnalysis.IIDPrior()
tasks = PhylogeneticFactorAnalysis.PipelineTasks()



selection_mcmc = MCMCOptions()
final_mcmc = MCMCOptions()

pipeline_input = PhylogeneticFactorAnalysis.PipelineInput(
            data_path, newick_path, "", msp, prior, tasks,
            selection_mcmc, final_mcmc)

rm.(dummy_files);

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

# throw_error_if_exists.(dummy_files)



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

iid_msp = PhylogeneticFactorAnalysis.ModelSelectionProvider([1, 2, 3], Float64[], 2, ["CLPD"])
iid_prior = PhylogeneticFactorAnalysis.IIDPrior("none")

shrinkage_msp = PhylogeneticFactorAnalysis.ModelSelectionProvider([3, 3, 3], [10.0, 100.0, 1000.0], 2, ["CLPD"])
shrinkage_prior = PhylogeneticFactorAnalysis.ShrinkagePrior(NaN, "shape", true, false, true)
tasks = PhylogeneticFactorAnalysis.PipelineTasks()



selection_mcmc = MCMCOptions()
final_mcmc = MCMCOptions()

iid_input = PhylogeneticFactorAnalysis.PipelineInput(
            "iid",
            data_path, newick_path, iid_msp, iid_prior)

shrink_input = PhylogeneticFactorAnalysis.PipelineInput(
            "shrink",
            data_path, newick_path, shrinkage_msp, shrinkage_prior)


PhylogeneticFactorAnalysis.run_pipeline(iid_input)
PhylogeneticFactorAnalysis.run_pipeline(shrink_input)


rm.(dummy_files);

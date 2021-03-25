using PhylogeneticFactorAnalysis

using BEASTXMLConstructor, Random, BEASTTreeUtils, CSV, DataFrames

Random.seed!(666)

const CLEAN = true




function throw_error_if_exists(path::String)
    if isfile(path)
        error("File $path already exists. Delete and re-run.")
    end
end

data_path = joinpath(@__DIR__, "test_data.csv")
newick_path = joinpath(@__DIR__, "test_newick.txt")
labels_path = joinpath(@__DIR__, "test_labels.csv")

dummy_files = [data_path, newick_path, labels_path]

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

@assert iseven(P) # needed for pretty label test
ordered_labels = Vector{String}(undef, P)
erng = (div(P, 2) + 1):P
orng = 1:div(P, 2)
ordered_labels[orng] .= ["trait$i" for i = 1:2:P]
ordered_labels[erng] .= ["trait$i" for i = P:-2:2]
cats = fill("", P)
cats[erng] .= "even"
cats[orng] .= "odd"
labels_df = DataFrame(trait = ordered_labels,
                    pretty = [x * " pretty" for x in ordered_labels],
                    cat = cats)
CSV.write(labels_path, labels_df)


write(newick_path, writeTopology(rtree(taxa, ultrametric=true)))

iid_msp = PhylogeneticFactorAnalysis.ModelSelectionProvider([1, 2, 3], Float64[], 2)
iid_prior = PhylogeneticFactorAnalysis.IIDPrior("none")

shrinkage_msp = PhylogeneticFactorAnalysis.ModelSelectionProvider([3, 3, 3], [10.0, 100.0, 1000.0], 2)
shrinkage_prior = PhylogeneticFactorAnalysis.ShrinkagePrior(NaN, "shape", true, false, true)
tasks = PhylogeneticFactorAnalysis.PipelineTasks()



selection_mcmc = MCMCOptions(chain_length = 100)
final_mcmc = MCMCOptions()

iid_input = PhylogeneticFactorAnalysis.PipelineInput(
            "iid",
            data_path, newick_path, iid_msp, iid_prior,
            selection_mcmc = selection_mcmc, overwrite = true,
            labels_path = labels_path)

shrink_input = PhylogeneticFactorAnalysis.PipelineInput(
            "shrink",
            data_path, newick_path, shrinkage_msp, shrinkage_prior,
            selection_mcmc = selection_mcmc, overwrite=true)

try
    # PhylogeneticFactorAnalysis.run_pipeline(iid_input)
    # PhylogeneticFactorAnalysis.run_pipeline(shrink_input)
    # iid_input.tasks.make_selection_xml = false
    # iid_input.tasks.run_selection_xml = false
    # iid_input.tasks.run_final_xml = false
    # iid_input.tasks.record_selection_stats = false
    # iid_input.tasks.process_final_log = false
    # iid_input.tasks.plot_loadings = false
    PhylogeneticFactorAnalysis.run_pipeline(iid_input)
catch e
    @error "Something went wrong" exception=(e, catch_backtrace())
    cd(@__DIR__)
end

cd(@__DIR__)

if CLEAN
    rm.(dummy_files);
    for file in readdir()
        if file[end-2:end] in ["xml", "txt", "log"]
            rm(file)
        end
    end
end


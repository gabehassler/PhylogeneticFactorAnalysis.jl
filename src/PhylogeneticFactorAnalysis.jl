module PhylogeneticFactorAnalysis

using BEASTXMLConstructor, BeastUtils.DataStorage, BeastUtils.MatrixUtils, UnPack

struct ModelSelectionProvider
    n_factors::Vector{Int}
    shrinkage_mults::Vector{Float64} # only relevant if using shrinkage prior
    reps::Int
    statistics::Vector{String}

    function ModelSelectionProvider(n_factors::Vector{Int},
                                    shrinkage_mults::Vector{Float64},
                                    reps::Int,
                                    statistics::Vector{String})
        if length(n_factors) != length(shrinkage_mults) && !isempty(shrinkage_mults)
            throw(ArgumentError("The length of the shrinkage multipliers must" *
                  " be 0 or the length of the factors ($(length(n_factors)))"))
        end

        return new(n_factors, shrinkage_mults, reps, statistics)
    end

end

import Base.length
function length(m::ModelSelectionProvider)
    return length(m.n_factors)
end

struct GeneralModelingDecisions
    standardize_data::Bool
end

abstract type PriorParameters end

struct IIDPrior <: PriorParameters
    constraint::String
end

struct ShrinkagePrior <: PriorParameters
    multiplier::Float64 #TODO: remove?
    shrink_by::String
    fix_first::Bool
    shrink_first::Bool
    force_ordered::Bool
end

struct PipelineTasks
    make_selection_xml::Bool
    run_selection_xml::Bool
    record_selection_stats::Bool
    make_final_xml::Bool
    run_final_xml::Bool
    process_final_log::Bool
    plot_loadings::Bool
    plot_factors::Bool

    function PipelineTasks(;
                            make_selection_xml::Bool = true,
                            run_selection_xml::Bool = true,
                            record_selection_stats::Bool = true,
                            make_final_xml::Bool = true,
                            run_runal_xml::Bool = true,
                            process_final_log::Bool = true,
                            plot_loadings::Bool = true,
                            plot_factors::Bool = true
                          )
        return new(make_selection_xml,
                    run_selection_xml,
                    record_selection_stats,
                    make_final_xml,
                    run_runal_xml,
                    process_final_log,
                    plot_loadings,
                    plot_factors)
    end
end


struct PipelineInput
    name::String
    data_path::String
    tree_path::String
    labels_path::String
    trait_data::TraitData
    newick::String


    model_selection::ModelSelectionProvider
    prior::PriorParameters
    tasks::PipelineTasks

    selection_mcmc::MCMCOptions
    final_mcmc::MCMCOptions
    standardize_data::Bool


    function PipelineInput(name::String,
                           data_path::String,
                           tree_path::String,
                           model_selection::ModelSelectionProvider,
                           prior::PriorParameters;
                           labels_path::String = "",
                           tasks::PipelineTasks = PipelineTasks(),
                           selection_mcmc::MCMCOptions = MCMCOptions(),
                           final_mcmc::MCMCOptions = MCMCOptions(),
                           standardize_data::Bool = true
                           )
        td = csv_to_traitdata(data_path)
        newick = read(tree_path, String)
        return new(name, data_path, tree_path, labels_path, td, newick,
                   model_selection,
                   prior,
                   tasks,
                   selection_mcmc,
                   final_mcmc,
                   standardize_data)
    end
end

include("make_xml.jl")

function run_pipeline(input::PipelineInput)
    @unpack tasks = input

    if tasks.make_selection_xml
        make_selection_xml(input)
    end
    if tasks.run_selection_xml
        #TODO
    end
    if tasks.record_selection_stats
        #TODO
    end
    if tasks.make_final_xml
        #TODO
    end
    if tasks.run_final_xml
        #TODO
    end
    if tasks.process_final_log
        #TODO
    end
    if tasks.plot_loadings
        #TODO
    end
    if tasks.plot_factors
        #TODO
    end
end

function make_selection_xml(input::PipelineInput)
    @unpack model_selection, trait_data, name, prior, selection_mcmc = input
    @unpack reps = model_selection
    @unpack data = trait_data

    validation_data = copy(trait_data.data)
    if input.standardize_data
        standardize_data!(validation_data)
    end

    training_data = copy(validation_data)

    rng = 1:reps

    assignments = rand(rng, length(validation_data))

    newick = read(input.tree_path, String)

    for r = 1:reps
        training_data .= validation_data
        training_data[findall(x -> x == r, assignments)] .= NaN
        for m = 1:length(model_selection)
            make_training_xml(input, training_data, validation_data, m, r, standardize = false)
        end
    end
end

end
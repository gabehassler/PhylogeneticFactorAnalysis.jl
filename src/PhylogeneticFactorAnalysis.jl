module PhylogeneticFactorAnalysis

using BEASTXMLConstructor

struct ModelSelectionProvider
    n_factors::Vector{Int}
    shrinkage_mults::Vector{Float64} # only relevant if using shrinkage prior
    reps::Int
end


abstract type PriorParameters end

struct IIDPrior <: PriorParameters
end

struct ShrinkagePrior <: PriorParameters
    multiplier::Float64
    shrink_by::String
    fix_first::Bool
    shrink_first::Bool
end

struct PipelineTasks
    make_selection_xml::Bool
    run_selection_xml::Bool
    record_selection_stats::Bool
    make_final_xml::Bool
    run_runal_xml::Bool
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

struct



struct PipelineInput
    data_path::String
    labels_path::String


    model_selection::ModelSelectionProvider
    prior::PriorParameters
    tasks::PipelineTasks

    selection_mcmc::MCMCOptions
    final_mcmc::MCMCOptions
end


end
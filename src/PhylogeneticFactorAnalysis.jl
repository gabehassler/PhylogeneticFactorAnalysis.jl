module PhylogeneticFactorAnalysis

using BEASTXMLConstructor, BeastUtils.DataStorage, BeastUtils.MatrixUtils, BeastUtils.RunBeast, BeastUtils.Logs,
      UnPack, Random, DataFrames, CSV, Statistics


include("PostProcessing.jl")
using PhylogeneticFactorAnalysis.PostProcessing


ModelStat = NamedTuple{(:model, :statistics),Tuple{Int64,Array{String,1}}}

struct ModelSelectionProvider
    n_factors::Vector{Int}
    shrinkage_mults::Vector{Float64} # only relevant if using shrinkage prior
    reps::Int
    statistics::Vector{String}
    burnin::Float64
    final_names::Dict{String, ModelStat}

    function ModelSelectionProvider(n_factors::Vector{Int},
                                    shrinkage_mults::Vector{Float64},
                                    reps::Int;
                                    statistics::Vector{String} = [LPD_COND],
                                    burnin::Float64 = 0.25)
        if length(n_factors) != length(shrinkage_mults) && !isempty(shrinkage_mults)
            throw(ArgumentError("The length of the shrinkage multipliers must" *
                  " be 0 or the length of the factors ($(length(n_factors)))"))
        end

        return new(n_factors, shrinkage_mults, reps, statistics, burnin, Dict{String, ModelStat}())
    end

end

import Base.length
function length(m::ModelSelectionProvider)
    return length(m.n_factors)
end

import Base.size
function size(m::ModelSelectionProvider)
    return (length(m.n_factors), m.reps)
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

mutable struct PipelineTasks
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
    directory::String
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
    julia_seed::Int
    beast_seed::Int

    overwrite::Bool


    function PipelineInput(name::String,
                           data_path::String,
                           tree_path::String,
                           model_selection::ModelSelectionProvider,
                           prior::PriorParameters;
                           labels_path::String = "",
                           tasks::PipelineTasks = PipelineTasks(),
                           selection_mcmc::MCMCOptions = MCMCOptions(),
                           final_mcmc::MCMCOptions = MCMCOptions(),
                           standardize_data::Bool = true,
                           julia_seed::Int = Int(rand(UInt32)),
                           beast_seed::Int = -1,
                           directory = pwd(),
                           overwrite::Bool = false
                           )
        td = csv_to_traitdata(data_path)
        newick = read(tree_path, String)
        return new(name, directory,
                   data_path, tree_path, labels_path, td, newick,
                   model_selection,
                   prior,
                   tasks,
                   selection_mcmc,
                   final_mcmc,
                   standardize_data,
                   julia_seed,
                   beast_seed,
                   overwrite)
    end
end

include("paths.jl")
include("make_xml.jl")

function run_pipeline(input::PipelineInput)

    @unpack tasks, julia_seed = input
    Random.seed!(julia_seed)

    dir = basedir(input)
    if isdir(dir) && !isempty(readdir(dir)) && !input.overwrite
        error("Director \"$dir\" is not empty. Set TODO to overwrite.")
    else
        mkpath(dir)
        mkpath(statistics_dir(input))
        mkpath(selection_log_dir(input))
        mkpath(selection_xml_dir(input))
        mkpath(timer_dir(input))
    end

    cd(dir)

    if tasks.make_selection_xml
        make_selection_xml(input)
    end

    if tasks.run_selection_xml
        run_selection_xml(input)
    end

    if tasks.record_selection_stats && !tasks.run_selection_xml
        m, r = size(input.model_selection)
        for i = 1:r
            for j = 1:m
                process_selection_statistics(input, j, i)
            end
        end
    end

    if tasks.make_final_xml
        make_final_xml(input)
    end
    if tasks.run_final_xml
        run_final_xml(input)
    end
    if tasks.process_final_log
        process_final_logs(input)
    end
    if tasks.plot_loadings
        #TODO
    end
    if tasks.plot_factors
        #TODO
    end
end


function process_final_logs(input::PipelineInput)
    log_paths = final_log_paths(input)
    processed_paths = processed_log_paths(input)

    for i = 1:length(log_paths)
        svd_logs(log_paths[i], processed_paths[i])
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
        #TODO: initialize with "L" for shrinkage?
        training_data .= validation_data
        training_data[findall(x -> x == r, assignments)] .= NaN
        for m = 1:length(model_selection)
            xml = make_training_xml(input, training_data, validation_data, m, r, standardize = false)
            xml_path = selection_xml_path(input, model = m, rep = r)
            mv(xml, xml_path, force=input.overwrite)
        end
    end
end

function run_selection_xml(input::PipelineInput)
    @unpack model_selection = input

    check_stats_exist(input)

    for r = 1:model_selection.reps
        for m = 1:length(model_selection)
            xml_path = selection_xml_path(input, model = m, rep = r)
            RunBeast.run_beast(xml_path, seed = input.beast_seed, overwrite = input.overwrite)

            log_filename = log_name(input, model = m, rep = r)
            log_path = selection_log_path(input, model = m, rep = r)
            mv(log_filename, log_path, force = input.overwrite)

            timer_filename = timer_name(input, model = m, rep = r)
            tp = timer_path(input, model = m, rep = r)
            mv(timer_filename, tp, force = input.overwrite)

            if input.tasks.record_selection_stats
                process_selection_statistics(input, m, r)
            end
        end
    end
end



function check_stats_exist(input::PipelineInput)
    for stat in input.model_selection.statistics
        path = get_stat_path(input, stat)
        if !isfile(path)
            m, r = size(input.model_selection)
            stats = fill(NaN, m, r)
            df = DataFrame(stats, ["rep$i" for i = 1:r])
            CSV.write(path, df)
        end
    end
end

function process_selection_statistics(input::PipelineInput, model::Int, rep::Int)
    @unpack model_selection = input

    log_path = selection_log_path(input, model = model, rep = rep)
    stats = compute_selection_statistics(log_path, model_selection)

    for i = 1:length(stats)
        stat = model_selection.statistics[i]
        stat_path = get_stat_path(input, stat)
        df = CSV.read(stat_path, DataFrame)
        df[model, rep] = stats[i]
        CSV.write(stat_path, df)
    end
end



function make_final_xml(input::PipelineInput)
    @unpack model_selection = input
    @unpack statistics, final_names = model_selection
    best_models = find_best_models(input)

    n = length(best_models)

    if n == 1 # just 1 selection statistic
        xml = make_final_xml(input, best_models[1])
        final_names[sans_extension(xml)] =
                (model = best_models[1], statistics = statistics)
    else
        condensed = condense_models(best_models, statistics)
        for model in keys(condensed)
            model_stats = condensed[model]
            xml = make_final_xml(input, model, statistic = join(model_stats, "and"))
            final_names[sans_extension(xml)] = (model = model, statistics = model_stats)
        end
    end
end


function run_final_xml(input::PipelineInput)
    xml_paths = final_xml_paths(input)
    for path in xml_paths
        RunBeast.run_beast(path, seed = input.beast_seed, overwrite = input.overwrite)
    end
end

function condense_models(models::Vector{Int}, statistics::Vector{String})
    @assert length(models) == length(statistics)
    u = unique(models)

    d = Dict{Int, Vector{String}}()

    for model in u
        inds = findall(isequal(model), models)
        d[model] = statistics[inds]
    end

    return d
end


function find_best_models(input::PipelineInput)
    @unpack model_selection = input
    @unpack statistics = model_selection

    n = length(statistics)
    best_models = fill(0, n)

    for i = 1:n
        stat = statistics[i]
        X = Matrix(CSV.read(get_stat_path(input, stat), DataFrame))
        u = vec(mean(X, dims=2)) * mult_dict[stat]
        best_models[i] = findmax(u)[2]
    end
    return best_models
end

function compute_selection_statistics(log_path::String, model_selection::ModelSelectionProvider)
    cols, data = Logs.get_log(log_path, burnin = model_selection.burnin)
    n = length(model_selection.statistics)

    means = zeros(n)

    for i = 1:n
        stat = model_selection.statistics[i]
        ind = findall(startswith(label_dict[stat]), cols)
        @assert length(ind) == 1
        ind = ind[1]

        stat_data = data[:, ind]

        if stat == LPD_COND #TODO: implement directly in BEAST?
            like_ind = findall(isequal("likelihood"), cols)
            @assert length(like_ind) == 1
            like_ind = like_ind[1]
            stat_data .-= @view data[:, like_ind]
        end

        means[i] = mean(stat_data)
    end

    return means
end

end
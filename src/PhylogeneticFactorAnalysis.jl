module PhylogeneticFactorAnalysis

using Pkg.Artifacts


using BEASTXMLConstructor
using BeastUtils, BeastUtils.DataStorage,
      BeastUtils.MatrixUtils, BeastUtils.RunBeast, BeastUtils.Logs,
      BeastUtils.PosteriorSummary
using UnPack, Random, DataFrames, CSV
using Statistics
using EzXML
using JLD
using LinearAlgebra
using RCall

import BeastUtils.DataStorage.TraitData
import BEASTXMLConstructor.MCMCOptions

export load_jld,
       pfa,
       parse_xml,
       run_pipeline,
       start_from,
       run_only,
       check_beast,
       check_r,
       import_example,
       run_example,
       get_data,
       plot_correlation,
       plot_loadings

post_dir = "PostProcessing"
pfa_dir = "PhylogeneticFactorAnalysis"
plotting_dir = "plotting"

# import standalone modules
include(joinpath(post_dir, "PostProcessing.jl"))
using PhylogeneticFactorAnalysis.PostProcessing

include(joinpath(post_dir, "BEASTPostProcessing.jl"))
using PhylogeneticFactorAnalysis.BEASTPostProcessing


const BEAST_HOME = "BEAST_HOME"


include(joinpath(pfa_dir, "dependencies.jl"))


function find_beast()
    if haskey(ENV, BEAST_HOME)
        return RunBeast.find_beast()
    else
        return RunBeast.find_beast(
                        beast_home = joinpath(artifact"beast_jar", "beast.jar"))
    end
end

function check_beast()
    return RunBeast.check_beast(beast_jar = find_beast())
end

ModelStat = NamedTuple{(:model, :statistics),Tuple{Int64,Array{String,1}}}

mutable struct Bookkeeper
    started::Int
    completed::Int
    total::Int
    function Bookkeeper(n::Int)
        return new(0, 0, n)
    end
end

mutable struct ModelSelectionProvider
    n_factors::Vector{Int}
    shrinkage_mults::Vector{Float64} # only relevant if using shrinkage prior
    reps::Int
    statistics::Vector{String}
    burnin::Float64
    mcmc_options::MCMCOptions
    final_names::Dict{String, ModelStat}

    function ModelSelectionProvider(n_factors::Vector{Int},
                                    shrinkage_mults::Vector{Float64},
                                    reps::Int;
                                    mcmc_options::MCMCOptions = MCMCOptions(),
                                    statistics::Vector{String} = [LPD_COND],
                                    burnin::Float64 = 0.25)
        if length(n_factors) != length(shrinkage_mults) && !isempty(shrinkage_mults)
            throw(ArgumentError("The length of the shrinkage multipliers must" *
                  " be 0 or the length of the factors ($(length(n_factors)))"))
        end

        return new(n_factors, shrinkage_mults, reps, statistics, burnin,
                   mcmc_options, Dict{String, ModelStat}())
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

abstract type PriorParameters end

struct IIDPrior <: PriorParameters
    constraint::String
end

struct ShrinkagePrior <: PriorParameters
    shrink_by::String
    fix_first::Bool
    shrink_first::Bool
    scale_first::Bool
    force_ordered::Bool
    spacing::Float64

    function ShrinkagePrior(;
                            shrink_by::String = SHAPE,
                            fix_first::Bool = true,
                            shrink_first::Bool = false,
                            scale_first::Bool = true,
                            force_ordered::Bool = true,
                            spacing::Float64 = 0.9)
        return new(shrink_by, fix_first, shrink_first, scale_first, force_ordered, spacing)
    end
    #TODO: add spacing
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

mutable struct PlotAttributes
    labels_path::String
    class_path::String
    burnin::Float64
    hpd_alpha::Float64
    scale_loadings_by_factors::Bool

    function PlotAttributes(;
                            labels_path::String = "",
                            class_path::String = "",
                            burnin::Float64 = 0.1,
                            hpd_alpha::Float64 = 0.05,
                            scale_loadings_by_factors::Bool = true
                            )
        return new(labels_path, class_path, burnin, hpd_alpha,
                   scale_loadings_by_factors)
    end
end

struct TraitsAndTree
    data_path::String
    tree_path::String
    trait_data::TraitData
    newick::String
    discrete_inds::Vector{Int}

    function TraitsAndTree(data_path::String, tree_path::String; discrete_inds::Vector{Int} = Int[])
        traits = parse_traitdata(data_path)
        newick = read(tree_path, String)

        return new(data_path, tree_path, traits, newick, discrete_inds)
    end
end

mutable struct ModelOptions # miscellaneous modeling options
    root_sample_size::Float64
    standardize_data::Bool
end

mutable struct PostProcessingOptions # instructions for post-processing
    rotation_plan::RotationPlan
end

function PostProcessingOptions()
    return PostProcessingOptions(RotationPlan(SVDRotation, SignRotation))
end

function ModelOptions(;
                      root_sample_size::Float64 = 1e-3,
                      standardize_data::Bool = false
                      )
    return ModelOptions(root_sample_size,
                        standardize_data)
end


mutable struct PipelineInput
    name::String
    directory::String

    data::TraitsAndTree

    model_selection::ModelSelectionProvider
    prior::PriorParameters
    model_options::ModelOptions
    processing_options::PostProcessingOptions
    tasks::PipelineTasks

    final_mcmc::MCMCOptions
    julia_seed::Int
    beast_seed::Int

    overwrite::Bool
    plot_attrs::PlotAttributes
    initialize_parameters::Bool
    merged_xml::String
    warnings::Vector{String}
    jar_path::String


    function PipelineInput(name::String,
                           data::TraitsAndTree,
                           model_selection::ModelSelectionProvider,
                           prior::PriorParameters;
                           model_options = ModelOptions(),
                           processing_options = PostProcessingOptions(),
                           tasks::PipelineTasks = PipelineTasks(),
                           final_mcmc::MCMCOptions = MCMCOptions(chain_length = 100_000),
                           julia_seed::Int = Int(rand(UInt32)),
                           beast_seed::Int = -1,
                           directory = pwd(),
                           overwrite::Bool = false,
                           plot_attrs::PlotAttributes = PlotAttributes(),
                           initialize_parameters::Bool = false,
                           merged_xml::String = "",
                           jar_path::String = find_beast()
                           )
        return new(name, directory,
                   data,
                   model_selection,
                   prior,
                   model_options,
                   processing_options,
                   tasks,
                   final_mcmc,
                   julia_seed,
                   beast_seed,
                   overwrite,
                   plot_attrs,
                   initialize_parameters,
                   merged_xml,
                   [""],
                   jar_path
                   )
    end
end

function dimensions(input::PipelineInput, model::Int)
    n, p = size(input.data.trait_data.data)
    k = input.model_selection.n_factors[model]
    return (n = n, p = p, k = k)
end

function get_data(input::PipelineInput)
    return input.data.trait_data
end

include(joinpath(pfa_dir, "paths.jl"))
include(joinpath(pfa_dir, "make_xml.jl"))
include(joinpath(pfa_dir, "writer.jl"))
include(joinpath(pfa_dir, "parsers.jl"))
include(joinpath(pfa_dir, "ui.jl"))
include(joinpath(pfa_dir, "display.jl"))
include(joinpath(pfa_dir, "tasks.jl"))
include("examples.jl")
include(joinpath(pfa_dir, "logs.jl"))
include(joinpath(plotting_dir, "plotting.jl"))
include(joinpath(plotting_dir, "outside_pipeline.jl"))
include(joinpath(pfa_dir, "post_processing.jl"))


function run_pipeline(input::PipelineInput)

    empty!(input.warnings)

    @unpack tasks, julia_seed = input
    Random.seed!(julia_seed)

    dir = basedir(input)
    if isdir(dir) && !isempty(readdir(dir)) && !input.overwrite
        error("Director \"$dir\" is not empty. Either manually delete the " *
              "directory or set $OVERWRITE=\"true\" in the $PFA xml block to " *
              "automatically overwrite.")
    else
        mkpath(dir)
        mkpath(statistics_dir(input))
        mkpath(selection_log_dir(input))
        mkpath(selection_xml_dir(input))
        mkpath(timer_dir(input))
    end
    original_dir = pwd()
    cd(dir)

    write_jld(input.name * "_backup.jld", input)

    check_tasks!(input)

    if tasks.make_selection_xml
        make_selection_xml(input)
    end

    if tasks.run_selection_xml
        run_selection_xml(input)
    end

    if tasks.record_selection_stats && !tasks.run_selection_xml
        check_stats_exist(input)
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
        plot_loadings(input)
    end
    if tasks.plot_factors
        plot_factors(input)
    end

    for warning in input.warnings
        @warn warning
    end

    cd(original_dir)
end

function check_tasks!(input::PipelineInput)
    @unpack tasks, model_selection = input
    if length(model_selection) == 1
        tasks.make_selection_xml = false
        tasks.run_selection_xml = false
        tasks.record_selection_stats = false

        if tasks.make_final_xml
            best_models = fill(1, length(model_selection.statistics))
            make_final_xml(input, best_models)
            tasks.make_final_xml = false
        end
    end
end

const L_HEADER = "L"
const PREC_HEADER = "factorPrecision"

function check_spacing!(::Vector{Float64}, ::IIDPrior)
    # do nothing
end

function check_spacing!(x::Vector{Float64}, prior::ShrinkagePrior)
    original_sum_squares = dot(x, x)

    spacing = prior.spacing
    for i = 2:length(x)
        max_val = 0.8 * spacing * abs(x[i - 1])
        if abs(x[i]) > max_val
            x[i] = max_val
        end
    end
    new_sum_squares = dot(x, x)

    return x .* sqrt(original_sum_squares / new_sum_squares)
end


function default_parameters(k::Int, p::Int)
    L = zeros(k, p)
    for i = 1:k
        L[i, i] = 1.0
    end

    return (L = L, precs = ones(p))
end

function initialize_parameters(input::PipelineInput,
                               data::Matrix{Float64},
                               model::Int)
    filename = make_init_xml(input, data, model, standardize = false)
    xml_path = init_xml_path(input)
    mv(filename, xml_path, force = true)

    run_beast(xml_path, seed = input.beast_seed, overwrite = true,
              beast_jar = input.jar_path)
    log_filename = log_name(input, stat=INIT)
    log_path = init_log_path(input)
    mv(log_filename, log_path, force = true)

    cols, data = get_log(log_path)
    L_cols = findall(startswith(L_HEADER), cols)
    F_cols = findall(startswith(FAC_HEADER), cols)


    @unpack n, k, p = dimensions(input, model)
    L = reshape(data[end, L_cols], p, k)'
    F = reshape(data[end, F_cols], k, n)
    d = (1 / n) * diag(F * F')
    Lsvd = svd(Diagonal(d) * L)

    @unpack S, Vt = Lsvd
    check_spacing!(S, input.prior)

    L = Diagonal(S) * Vt

    prec_cols = findall(startswith(PREC_HEADER), cols)
    precs = data[end, prec_cols]

    return (L = L, precs = precs)
end

function plot_factors(input::PipelineInput)
    log_paths = processed_log_paths(input)
    stat_paths = factors_statistic_paths(input)
    plot_paths = factors_plot_paths(input)

    for i = 1:length(log_paths)
        prep_factors(log_paths[i], stat_paths[i])
        factor_plot(plot_paths[i], stat_paths[i], input.data.tree_path,
                    input.plot_attrs.class_path)
    end
end

function standardize_continuous!(X::Matrix{Float64}, discrete_inds::Vector{Int})
    continuous_inds = setdiff(1:size(X, 2), discrete_inds)
    continuous_data = X[:, continuous_inds]
    standardize_data!(continuous_data)
    X[:, continuous_inds] .= continuous_data
    return X
end

function make_selection_xml(input::PipelineInput)
    @unpack model_selection, data, name, prior, model_options = input
    @unpack trait_data, newick, discrete_inds = data
    @unpack reps = model_selection
    @unpack data = trait_data

    validation_data = copy(trait_data.data)
    if model_options.standardize_data
        standardize_continuous!(validation_data, discrete_inds)
    end

    training_data = copy(validation_data)

    rng = 1:reps

    assignments = rand(rng, length(validation_data))

    for r = 1:reps
        L_inits = Dict{Int, ModelParams}()

        training_data .= validation_data
        training_data[findall(x -> x == r, assignments)] .= NaN
        for m = 1:length(model_selection)

            @unpack k, p = dimensions(input, m)
            if !haskey(L_inits, k)
                if input.initialize_parameters
                    params = initialize_parameters(input, training_data, m)
                    L_inits[k] = params
                else
                    L_inits[k] = default_parameters(k, p)
                end
            end

            xml = make_training_xml(input, training_data, validation_data, m, r,
                                    L_inits[k], standardize = false)
            xml_path = selection_xml_path(input, model = m, rep = r)
            mv(xml, xml_path, force=input.overwrite)
        end
    end
end

function run_selection_xml(input::PipelineInput, rep::Int, model::Int,
                           lk::ReentrantLock, bookkeeper::Bookkeeper)
    xml_path = selection_xml_path(input, model = model, rep = rep)

    begin
        lock(lk)
        try
            bookkeeper.started += 1
            println("starting $(basename(xml_path)) ($(bookkeeper.started)" *
                " of $(bookkeeper.total) started)")
        finally
            unlock(lk)
        end
    end

    capture_output = Threads.nthreads() > 1

    RunBeast.run_beast(xml_path, seed = input.beast_seed,
                        overwrite = input.overwrite,
                        beast_jar = input.jar_path,
                        capture_output = capture_output)

    begin
        lock(lk)
        try
            bookkeeper.completed += 1
            println("completed $(basename(xml_path)) " *
                "($(bookkeeper.completed) of $(bookkeeper.total) completed)")
        finally
            unlock(lk)
        end
    end

    log_filename = log_name(input, model = model, rep = rep)
    log_path = selection_log_path(input, model = model, rep = rep)
    mv(log_filename, log_path, force = input.overwrite)

    timer_filename = timer_name(input, model = model, rep = rep)
    tp = timer_path(input, model = model, rep = rep)
    mv(timer_filename, tp, force = input.overwrite)

    if input.tasks.record_selection_stats
        process_selection_statistics(input, model, rep, lk = lk)
    end
end


function run_selection_xml(input::PipelineInput)
    @unpack model_selection = input

    check_stats_exist(input)
    stats_lock = ReentrantLock()
    total_xml = model_selection.reps * length(model_selection)
    bookkeeper = Bookkeeper(total_xml)

    if Threads.nthreads() > 1
        println("Model selection using multiple threads. " *
        "BEAST screen output will be suppressed for model selection only.\n")
    end


    Threads.@threads for r = 1:model_selection.reps
        @sync for m = 1:length(model_selection)
            run_selection_xml(input, r, m, stats_lock, bookkeeper)
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

function process_selection_statistics(input::PipelineInput, model::Int,
                                      rep::Int;
                                      lk::ReentrantLock = ReentrantLock())
    @unpack model_selection = input

    log_path = selection_log_path(input, model = model, rep = rep)
    stats = compute_selection_statistics(input, log_path, model_selection)

    for i = 1:length(stats)
        stat = model_selection.statistics[i]
        stat_path = get_stat_path(input, stat)
        begin
            lock(lk)
            try
                df = CSV.read(stat_path, DataFrame)
                df[model, rep] = stats[i]
                CSV.write(stat_path, df)
            finally
                unlock(lk)
            end
        end

    end
end



function make_final_xml(input::PipelineInput)
    @unpack model_selection = input
    @unpack statistics, final_names = model_selection
    best_models = find_best_models(input)

    return make_final_xml(input, best_models)
end

function make_final_xml(input::PipelineInput, best_models::Vector{Int})
    @unpack model_selection = input
    @unpack statistics, final_names = model_selection

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
        RunBeast.run_beast(path, seed = input.beast_seed,
                           overwrite = input.overwrite,
                           beast_jar = input.jar_path)
    end
end


function process_final_logs(input::PipelineInput)
    log_paths = final_log_paths(input)
    processed_paths = processed_log_paths(input)
    models = input.model_selection.final_names

    for i = 1:length(log_paths)
        best_model = models[sans_extension(log_paths[i])].model
        @unpack k, p, n = dimensions(input, best_model)
        @unpack rows, cols = rotation_rows_and_cols(input, best_model)

        postprocess_log(log_paths[i], processed_paths[i], k, p, n,
                 rows, cols,
                 input.processing_options.rotation_plan)
    end
end

function rotation_rows_and_cols(input::PipelineInput, model::Int)
    @unpack k, p = dimensions(input, model)
    if typeof(input.prior) == IIDPrior
        constraint = input.prior.constraint
        if constraint == HYBRID
            return (rows = collect(2:k), cols = collect(2:p))
        elseif constraint == UPPER_TRIANGULAR
            return (rows = Int[], cols = Int[])
        elseif constraint == ORTHOGONAL
            return (rows = collect(1:k), cols = collect(1:p))
        else
            error("unrecognized constraint '$constraint'")
        end
    end

    return (rows = collect(1:k), cols = collect(1:p))
end

function plot_loadings(input::PipelineInput)
    log_paths = processed_log_paths(input)
    stat_paths = loadings_statistic_paths(input)
    plot_paths = loadings_plot_paths(input)
    for i = 1:length(log_paths)
        prep_loadings(input, log_paths[i], stat_paths[i])
        load_plot(plot_paths[i], stat_paths[i], labels_path = input.plot_attrs.labels_path)
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

function compute_selection_statistics(input::PipelineInput, log_path::String,
                                model_selection::ModelSelectionProvider)
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

        ess = effective_sample_size(stat_data)
        if ess < 100
            warning = "Selection statistic $stat in file $(basename(log_path)) " *
                      "has low effective sample size ($ess)."
            @warn warning

            push!(input.warnings, warning)
        end

        means[i] = mean(stat_data)
    end

    return means
end

end
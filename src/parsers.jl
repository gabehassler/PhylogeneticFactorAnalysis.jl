function parse_xml(xml::EzXML.Document)
    pfa = xml.root
    parse_pfa(pfa)
end

function parse_pfa(node::EzXML.Node)
    @assert node.name == PFA
    julia_seed = attr(node, JULIA_SEED, Int, default=rand(UInt32))
    directory = attr(node, DIRECTORY, String, default= pwd())

    data = parse_data(get_child_by_name(node, DATA),
                      alternative_directory = directory)

    default_name = sans_extension(basename(data.data_path))
    nm = attr(node, NAME, String, default = default_name)

    standardize = attr(node, STANDARDIZE, Bool, default=false)
    overwrite = attr(node, OVERWRITE, Bool, default=false)



    model_selection = parse_child(node, MODEL_SELECTION, parse_model_selection)
    prior_node = find_prior(node)
    if prior_node.name == IID_PRIOR
        prior = parse_iid_prior(find_prior(node))
    elseif prior_node.text == SHRINKAGE_PRIOR
        error("not yet implemented")
    end

    tasks = parse_child(node, TASKS, parse_tasks, default = PipelineTasks())
    plots = parse_child(node, PLOTS, parse_plots, default = PlotAttributes())
    final_mcmc = parse_child(node, MCMC, parse_mcmc, default = MCMCOptions(chain_length = 100_000))

    return PipelineInput(nm, data, model_selection, prior,
                         tasks = tasks,
                         julia_seed = julia_seed,
                         directory = directory,
                         plot_attrs = plots,
                         standardize_data = standardize,
                         overwrite = overwrite,
                         final_mcmc = final_mcmc)
    # plots = parse_plots(child_nodes)

end

function parse_child(node::EzXML.Node, child_name::String, parser::Function; default = nothing)
    optional = !isnothing(default)
    child_node = get_child_by_name(node, child_name, optional = optional)
    if isnothing(child_node) && optional
        return default
    end

    return parser(child_node)
end



function check_path(path::String, dir::String)
    if isfile(path)
        return abspath(path)
    elseif isfile(joinpath(dir, path))
        return abspath(joinpath(dir, path))
    else
        if dir == pwd()
            error("Cannot locate file $(abspath(path))")
        else
            error("Cannot locate file $(abspath(path)) or $(abspath(joinpath(dir, path)))")
        end
    end
end


const PFA = "phylogeneticFactorAnalysis"

const NAME = "name"
const DIRECTORY = "directory"
const DATA = "data"
const DATA_PATH = "traits"
const TREE_PATH = "tree"
const JULIA_SEED = "partitionSeed"
const BEAST_SEED = "mcmcSeed"
const STANDARDIZE = "standardizeData"
const OVERWRITE = "overwrite"
const CHAIN_LENGTH = "chainLength"
const MCMC = "mcmcOptions"

const TASKS = "tasks"
const PLOTS = "plotting"

const KEYS = Dict{String, Vector{Pair{String, Symbol}}}(
    PFA =>
        [NAME => :name,
        DIRECTORY => :directory,
        JULIA_SEED => :julia_seed],
)

import Base.parse
function parse(::Type{T}, s::AbstractString) where T <: AbstractString
    return convert(T, s)
end

function attr(node::EzXML.Node, key::String, type::Type; default = nothing)
    if haskey(node, key)
        return parse(type, node[key])
    elseif !isnothing(default)
        return convert(type, default)
    else
        error("Node $(node.name) does not have required attribute '$key'")
    end
end


function get_child_by_name(node::EzXML.Node, nm::String; optional::Bool = false)
    children = nodes(node)
    ind = findfirst(x -> x.name == nm, children)
    if isnothing(ind)
        if optional
            return nothing
        else
            error("$(node.name) must have child with name $nm")
        end
    end
    return children[ind]
end

function parse_text_array(node::EzXML.Node, type::Type)
    @assert node.type == EzXML.TEXT_NODE

    return parse.(type, split(string(node)))
end

function parse_text_array_child(node::EzXML.Node, type::Type)
    children = nodes(node)

    if length(children) != 1 || children[1].type != EzXML.TEXT_NODE
        error("Element $(node.name) must have exactly one text child")
    end

    return parse_text_array(children[1], type)
end

function parse_data(node::EzXML.Node; alternative_directory::String)
    data_path = attr(node, DATA_PATH, String)
    tree_path = attr(node, TREE_PATH, String)

    data_path = check_path(data_path, alternative_directory)
    tree_path = check_path(tree_path, alternative_directory)
    return TraitsAndTree(data_path, tree_path)
end


function parse_model_selection(node::EzXML.Node)
    reps = attr(node, REPEATS, Int, default=10)
    stats = attr(node, SELECTION_STAT, String, default = "CLPD")
    stats = string.(split(stats))

    for stat in stats
        if !(stat in STATS)
            error("Statistic $stat not recognized. Must be one of the followin: " *
                   join(STATS, ", "))
        end
    end

    burnin = attr(node, BURNIN, Float64, default=0.25)
    factors = get_child_by_name(node, N_FACTORS)
    factors = parse_text_array_child(factors, Int)

    shrinkage = get_child_by_name(node, SHRINKAGE, optional = true)
    if isnothing(shrinkage)
        shrinkage = Float64[]
    else
        shrinkage = parse_text_array_child(shrinkage, Float64)
    end

    mcmc = get_child_by_name(node, MCMC, optional = true)
    if isnothing(shrinkage)
        mcmc = MCMCOptions()
    else
        mcmc = parse_mcmc(mcmc)
    end

    return ModelSelectionProvider(factors, shrinkage, reps,
                                  statistics = stats, burnin = burnin,
                                  mcmc_options = mcmc)
end

const CONSTRAINTS = ["orthogonal", "upperTriangular", "hybrid"]

function find_prior(node::EzXML.Node)
    children = nodes(node)
    prior_inds = findall(x -> x.name in [SHRINKAGE_PRIOR, IID_PRIOR], children)
    if length(prior_inds) > 1
        error("More that one prior supplied (" *
            join([x.name for x in children[prior_inds]], ", ") * "). " *
            "There must be exactly one prior.")
    elseif length(prior_inds) == 0
        error("No prior supplied.")
    end

    return children[prior_inds[1]]
end

function parse_iid_prior(iid::EzXML.Node)
    constraint = attr(iid, CONSTRAINT, String)
    if !(lowercase(constraint) in lowercase.(CONSTRAINTS))
        error("Constraint '$constraint' not recognized. Must be one of: " *
            join(CONSTRAINTS, ", "))
    end
    return IIDPrior(constraint)
end

function parse_tasks(node::EzXML.Node)
    error("not implemented")
end

function parse_mcmc(node::EzXML.Node)
    chain_length = attr(node, CHAIN_LENGTH, Int, default = 10_000)
    return MCMCOptions(chain_length = chain_length)
end

function parse_plots(node::EzXML.Node)
    error("not implemented")
end
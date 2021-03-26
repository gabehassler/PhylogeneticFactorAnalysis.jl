function parse_xml(xml::EzXML.Document)
    pfa = xml.root
    parse_pfa(pfa)
end

function parse_pfa(node::EzXML.Node)
    @assert node.name == PFA
    julia_seed = attr(node, JULIA_SEED, Int, default=rand(UInt32))
    directory = attr(node, DIRECTORY, String, default= pwd())
    data_path = attr(node, DATA_PATH, String)
    tree_path = attr(node, TREE_PATH, String)

    data_path = check_path(data_path, directory)
    tree_path = check_path(data_path, directory)

    default_name = sans_extension(basename(data_path))
    nm = attr(node, NAME, String, default = default_name)


    #TODO: other attributes

    model_selection = parse_model_selection(get_child_by_name(node, MODEL_SELECTION))
    prior_node = find_prior(node)
    if prior_node.name == IID_PRIOR
        prior = parse_iid_prior(find_prior(node))
    elseif prior_node.text == SHRINKAGE_PRIOR
        error("not yet implemented")
    end

    task_node = get_child_by_name(node, TASKS, optional = true)
    if isnothing(task_node)
        tasks = PipelineTasks()
    else
        tasks = parse_tasks(task_node)
    end

    plot_node = get_child_by_name(node, PLOTS, optional = true)
    if isnothing(plot_node)
        plots = PlotAttributes()
    else
        plots = parse_plots(plot_node)
    end

    return PipelineInput(nm, data_path, tree_path, model_selection, prior,
                         tasks = tasks,
                         julia_seed = julia_seed,
                         directory = directory,
                         plot_attrs = plots)
    # plots = parse_plots(child_nodes)

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
const DATA_PATH = "data"
const TREE_PATH = "tree"
const JULIA_SEED = "juliaSeed"
const BEAST_SEED = "beastSeed"

const TASKS = "tasks"
const PLOTS = "plotting"

const KEYS = Dict{String, Vector{Pair{String, Symbol}}}(
    PFA =>
        [NAME => :name,
        DIRECTORY => :directory,
        JULIA_SEED => :julia_seed],
)

import Base.parse
function parse(type::Type{T}, s::AbstractString) where T <: AbstractString
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

    return ModelSelectionProvider(factors, shrinkage, reps,
                                  statistics = stats, burnin = burnin)
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
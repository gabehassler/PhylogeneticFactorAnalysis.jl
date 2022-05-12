function load_jld(x::String)
    # jldopen("iid.jld") do file
    #     addrequire(file, PhylogeneticFactorAnalysis)
    #     addrequire(file, BeastUtils.DataStorage)
    #     write(file, "pfa_input", input)
    # end
    input = jldopen(x, "r") do file
        read(file, "pfa_input")
    end
    # return load(x)["pfa_input"]
    return input
    # return true
end

function write_jld(x::String, input::PipelineInput)
    jldopen(x, "w") do file
        addrequire(file, PhylogeneticFactorAnalysis)
        addrequire(file, BeastUtils)
        addrequire(file, BEASTXMLConstructor)
        write(file, "pfa_input", input)
    end
end

function xml_doc(input::PipelineInput)
    doc = XMLDocument()
    root = make_element(input)
    setroot!(doc, root)
    return doc
end


function make_element(obj)
    node = ElementNode(name(obj))

    sub_components = xml_components(obj)
    for component in sub_components
        el = make_element(component)
        link!(node, el)
    end
    for attr in xml_attributes(obj)
        node[attr[1]] = attr[2]
    end

    return node
end

function make_element(x::EzXML.Node)
    return x
end

# function make_element(text::TextNode)
#     return text
# end

# PipelineInput
function xml_attributes(input::PipelineInput)
    return [NAME => input.name,
            DATA_PATH => input.data_path,
            TREE_PATH => input.tree_path,
            DIRECTORY => input.directory,
            STANDARDIZE => string(input.standardize_data),
            JULIA_SEED => string(input.julia_seed),
            BEAST_SEED => string(input.beast_seed),
            OVERWRITE => string(input.overwrite)]
end

function xml_components(input::PipelineInput)
    return [input.model_selection, input.prior]
end

function name(::PipelineInput)
    return "phylogeneticFactorAnalysis"
end

# iid prior

const CONSTRAINT = "constraint"

function xml_attributes(prior::IIDPrior)
    return [CONSTRAINT => prior.constraint]
end

function xml_components(::IIDPrior)
    return Nothing[]
end

function name(::IIDPrior)
    return IID_PRIOR
end

const IID_PRIOR = "iidPrior"
const SHRINKAGE_PRIOR = "orthogonalShrinkagePrior"

# model selection
const REPEATS = "repeats"
const SELECTION_STAT = "selectionStatistic"
const BURNIN = "burnin"
const MODEL_SELECTION = "modelSelection"
const N_FACTORS = "nFactors"
const SHRINKAGE = "shrinkageStrength"

function nv(nm::String, val::Array{T}) where T
    return (name = nm, value = val)
end

function xml_attributes(msp::ModelSelectionProvider)
    return [REPEATS => string(msp.reps),
            SELECTION_STAT => join(msp.statistics, ' '),
            BURNIN => string(msp.burnin)
            ]
end

function xml_components(msp::ModelSelectionProvider)
    comps = [nv(N_FACTORS, msp.n_factors)]
    if length(msp.shrinkage_mults) > 0
        push!(comps, nv(SHRINKAGE, msp.shrinkage_mults))
    end
    return comps
end

function name(::ModelSelectionProvider)
    return MODEL_SELECTION
end

# other
NT = NamedTuple{(:name, :value),Tuple{String,Array{T,1}}} where T


function name(x::NT)
    return x.name
end

function xml_components(x::NT)
    return [TextNode(join(x.value, ' '))]
end

function xml_attributes(::NT)
    return Nothing[]
end





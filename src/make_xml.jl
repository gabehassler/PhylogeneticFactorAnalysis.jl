const MSE_STAT = "MSE"
const LPD_COND = "CLPD"
# const LPD_MARG = "MLPD"
const NO_STAT = ""
const JOINT = "joint"
# const REMOVED = "removed"
# const partition_dict = Dict(LPD_COND => false, MSE_STAT => false)
const trait_dict = Dict(LPD_COND => JOINT, MSE_STAT => JOINT)
const label_dict = Dict(LPD_COND => "$(trait_dict[LPD_COND]).", MSE_STAT => "traitValidation.TotalSum")
const mult_dict = Dict(LPD_COND => 1, MSE_STAT => -1)

const SHAPE = "shape"
const SCALE = "scale"
const BOTH = "both"

const LOADINGS_WEIGHT = 3.0

function set_common_options(bx::BEASTXMLElement, options::MCMCOptions; standardize::Bool = true)

    set_options!(bx, options)
    lgo = BEASTXMLConstructor.get_loadings_op(bx, component="matrix")
    lgo.weight = LOADINGS_WEIGHT

    facs = BEASTXMLConstructor.get_integratedFactorModel(bx)
    facs.standardize_traits = standardize
end


function make_final_xml(input::PipelineInput, model::Int; statistic::String = "")

    @unpack trait_data, newick, model_selection, prior, final_mcmc = input
    @unpack data, taxa = trait_data

    bx = make_initial_xml(data, taxa, newick, model_selection, prior, model, log_factors = true)
    set_common_options(bx, final_mcmc, standardize = input.standardize_data)

    filename = xml_name(input, stat=statistic)
    path = BEASTXMLConstructor.save_xml(filename, bx)
    return filename
end


function make_training_xml(input::PipelineInput, training_data::Matrix{Float64},
                           validation_data::Matrix{Float64},
                           model::Int, rep::Int; standardize::Bool = true)

    @unpack name, trait_data, newick, model_selection, prior, selection_mcmc = input
    bx = make_initial_xml(training_data, trait_data.taxa, newick, model_selection, prior, model, log_factors = false)

    set_common_options(bx, selection_mcmc, standardize = standardize)

    add_validation(bx, validation_data, model_selection.statistics)

    filename = xml_name(input, model = model, rep = rep)
    path = BEASTXMLConstructor.save_xml(filename, bx)
    return filename
end

function add_validation(bx::BEASTXMLElement, validation_data::Matrix{Float64}, statistics::Vector{String})
    BEASTXMLConstructor.add_trait!(bx, validation_data, JOINT)

    like = BEASTXMLConstructor.get_traitLikelihood(bx)
    facs = BEASTXMLConstructor.get_integratedFactorModel(bx)

    for stat in statistics
        if stat == LPD_COND
            flpd = BEASTXMLConstructor.FactorLogPredictiveDensity(facs, like, trait_name = JOINT)
            BEASTXMLConstructor.add_loggable(bx, flpd, already_made = false)
        elseif stat == MSE_STAT
            tree_model = BEASTXMLConstructor.get_treeModel(bx)

            trait_validation =
                BEASTXMLConstructor.TraitValidationXMLElement(tree_model, like)

            cross_validation =
                BEASTXMLConstructor.CrossValidationXMLElement(trait_validation)

            BEASTXMLConstructor.set_validation_type!(cross_validation,
                                                BeastNames.SQUARED_ERROR)

            BEASTXMLConstructor.set_log_sum!(cross_validation, true)

            BEASTXMLConstructor.add_loggable(bx, cross_validation,
                                        already_made=false)
        else
            error("unknown statistics: $stat")
        end
    end
end


function shrinkage_shapes_and_scales(k::Int, shrink::Float64, prior::ShrinkagePrior)
    means = [shrink for i = 1:k]

    if !prior.shrink_first
        means[1] = 1.0
    end

    scales = zeros(k)
    shapes = zeros(k)
    if prior.shrink_by == SCALE
        shapes .= 1.0
        scales .= means
    elseif prior.shrink_by == SHAPE
        shapes .= means
        scales .= 1.0
    elseif prior.shrink_by == BOTH
        shapes .= sqrt.(means)
        scales .= shapes
    else
        error("unrecognized `shrink_by` field: \"$(prior.shrink_by)\"")
    end

    return (shapes = shapes, scales = scales)
end

function make_initial_xml(data::Matrix{Float64},
                          taxa::Vector{String},
                          newick::String,
                          selection_vars::ModelSelectionProvider,
                          prior::ShrinkagePrior,
                          model::Int;
                          log_factors::Bool = false) # TODO: get from input

    k = selection_vars.n_factors[model]

    bx = BEASTXMLConstructor.make_orthogonal_pfa_xml(data, taxa, newick,
                                k,
                                fix_first = prior.fix_first,
                                shrink_first = prior.shrink_first,
                                timing=true,
                                log_factors = log_factors,
                                force_ordered  = prior.force_ordered)

    facs = BEASTXMLConstructor.get_integratedFactorModel(bx)

    set_scale = true # TODO: set based on whether we initialize with loadings

    shrink = selection_vars.shrinkage_mults[model]
    @unpack shapes, scales = shrinkage_shapes_and_scales(k, shrink, prior)


    BEASTXMLConstructor.set_shrinkage_mults!(facs,
                                             shapes = shapes,
                                             scales = scales,
                                             set_scale = set_scale,
                                             set_mults = true)

    ops = BEASTXMLConstructor.get_operators(bx)

    if prior.fix_first
        indices = collect(2:k)
        BEASTXMLConstructor.set_muliplicative_gamma_indices(bx, indices)
    end

    return bx
end

function make_initial_xml(data::Matrix{Float64},
                          taxa::Vector{String},
                          newick::String,
                          selection_vars::ModelSelectionProvider,
                          prior::IIDPrior,
                          model::Int;
                          log_factors::Bool = false)
    bx = BEASTXMLConstructor.make_pfa_xml(data, taxa, newick,
                                        selection_vars.n_factors[model],
                                        log_factors = log_factors,
                                        useHMC = false,
                                        timing = true)

    lgo = BEASTXMLConstructor.get_loadings_op(bx)
    lgo.sparsity = prior.constraint

    return bx
end


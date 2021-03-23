# struct TreeData
#     taxa::Vector{String}
#     data::Matrix{Float64}
#     newick::String
# end

const MSE_STAT = "MSE"
const LPD_COND = "CLPD"
# const LPD_MARG = "MLPD"
const NO_STAT = ""
const JOINT = "joint"
# const REMOVED = "removed"
const partition_dict = Dict(LPD_COND => false, MSE_STAT => false)
const trait_dict = Dict(LPD_COND => JOINT, MSE_STAT => JOINT)
const label_dict = Dict(LPD_COND => "$(trait_dict[LPD_COND]).", MSE_STAT => "traitValidation.TotalSum")
const mult_dict = Dict(LPD_COND => 1, MSE_STAT => -1)

const LOADINGS_WEIGHT = 3.0

function make_training_xml(input::PipelineInput, training_data::Matrix{Float64},
                           validation_data::Matrix{Float64},
                           model::Int, rep::Int; standardize::Bool = true)

    @unpack name, trait_data, newick, model_selection, prior, selection_mcmc = input
    bx = make_initial_xml(training_data, trait_data.taxa, newick, model_selection.n_factors[model], prior, log_factors = false)

    set_options!(bx, selection_mcmc)

    lgo = BEASTXMLConstructor.get_loadings_op(bx, component="matrix")
    lgo.weight = LOADINGS_WEIGHT

    facs = BEASTXMLConstructor.get_integratedFactorModel(bx)
    facs.standardize_traits = standardize

    add_validation(bx, validation_data, model_selection.statistics)

    filename = join([name, "model$model", "r$rep"], '_')
    BEASTXMLConstructor.save_xml(filename, bx)
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


function make_initial_xml(data::Matrix{Float64},
                          taxa::Vector{String},
                          newick::String,
                          k::Int,
                          prior::ShrinkagePrior;
                          log_factors::Bool = false) # TODO: get from input

    bx = BEASTXMLConstructor.make_orthogonal_pfa_xml(data, taxa, newick,
                                k,
                                fix_first = prior.fix_first,
                                shrink_first = prior.shrink_first,
                                timing=true,
                                log_factors = log_factors,
                                force_ordered  = prior.force_ordered)



end

function make_initial_xml(data::Matrix{Float64},
                          taxa::Vector{String},
                          newick::String,
                          k::Int,
                          prior::IIDPrior;
                          log_factors::Bool = false)
    bx = BEASTXMLConstructor.make_pfa_xml(data, taxa, newick,
                                        k,
                                        log_factors = log_factors,
                                        useHMC = false,
                                        timing = true)

    lgo = BEASTXMLConstructor.get_loadings_op(bx, component="matrix")
    lgo.sparsity = prior.constraint

    return bx
end


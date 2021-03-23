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

function make_xml(base_name::String,
                trait_data::TraitData,
                validation_data::Matrix{Float64},
                newick::String,
                selection::ModelSelectionProvider,
                prior::ShrinkagePrior,
                mcmc_options::MCMCOptions,
                model::Int,
                rep::Int;
                log_factors::Bool = false,
                standardize::Bool = true) # TODO: get from input

    @unpack data, taxa = trait_data
    bx = BEASTXMLConstructor.make_orthogonal_pfa_xml(data, taxa, newick,
                                selection.n_factors[model],
                                fix_first = prior.fix_first,
                                shrink_first = prior.shrink_first,
                                timing=true,
                                log_factors = log_factors,
                                force_ordered  = prior.force_ordered)
    # TODO: shrinkage specific stuff

end

function make_xml(base_name::String,
                  trait_data::TraitData,
                  validation_data::Matrix{Float64},
                  newick::String,
                  selection::ModelSelectionProvider,
                  prior::IIDPrior,
                  mcmc_options::MCMCOptions,
                  model::Int,
                  rep::Int;
                  log_factors::Bool = false,
                  standardize::Bool = true) # TODO: get from input
    @unpack data, taxa = trait_data
    bx = BEASTXMLConstructor.make_pfa_xml(data, taxa, newick,
                                        selection.n_factors[model],
                                        log_factors = log_factors,
                                        useHMC = false,
                                        timing = true)

    set_options!(bx, mcmc_options)
    lgo = BEASTXMLConstructor.get_loadings_op(bx, component="matrix")
    lgo.weight = LOADINGS_WEIGHT

    lgo.sparsity = prior.constraint

    facs = BEASTXMLConstructor.get_integratedFactorModel(bx)
    facs.standardize_traits = standardize



    if size(validation_data) == size(data) #TODO: better way to check whether to add extra data
        BEASTXMLConstructor.add_trait!(bx, validation_data, JOINT)

        like = BEASTXMLConstructor.get_traitLikelihood(bx)

        for stat in selection.statistics
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

    filename = join([base_name, "model$model", "r$rep"], '_')
    BEASTXMLConstructor.save_xml(filename, bx)
end


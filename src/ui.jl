function pfa(input_file::String)
    if endswith(input_file, ".xml")
        input = parse_xml(input_file)
    elseif endswith(input_file, ".jld")
        input = load_jld(input_file)
    elseif startswith(input_file, "<?xml")
        input = parse_xml(input_file, is_string = true)
    else
        error("Unknown file extension. Must be either an 'xml' file or 'jld' file.")
    end

    old_wd = pwd()
    try
        run_pipeline(input)
    catch e
        cd(old_wd)
        @error "Something went wrong" exception=(e, catch_backtrace())
    end
end

function pfa(;name::String,
            data_path::String,
            newick_path::String,
            partition_seed::Int = -1,
            mcmc_seed::Int = -1,
            overwrite::Bool = true,
            standardize_traits::Bool = false,
            loadings_prior::String = "iid",
            loadings_constraint::String = "orthogonal",
            factors::Vector{Int} = [2],
            discrete_indices = Int[])
    data = TraitsAndTree(data_path, newick_path, discrete_inds = discrete_indices)
    model_options = ModelOptions(standardize_data = standardize_traits)
    if loadings_prior == "iid"
        prior = IIDPrior(loadings_constraint)
    else
        error("Unknown or unimplemented prior type")
    end


    model_selection = ModelSelectionProvider(factors, nothing, 1)

    input = PipelineInput(name, data, model_selection, prior,
                         model_options = model_options,
                         julia_seed = partition_seed,
                         beast_seed = mcmc_seed,
                         overwrite = overwrite
                         )

    old_wd = pwd()
    try
        run_pipeline(input)
    catch e
        cd(old_wd)
        @error "Something went wrong" exception=(e, catch_backtrace())
    end
end


const R_PLOT_SCRIPT = joinpath(@__DIR__, "..", "R", "plots.R")
const LOAD_HEADER = "L"
const SV_HEADER = "sv"
const FAC_HEADER = "factors."


function prep_loadings(input::PipelineInput, log_path::String,
                            csv_path::String)

    @unpack plot_attrs, data = input
    @unpack trait_data = data
    @unpack burnin, hpd_alpha, scale_loadings_by_factors = plot_attrs #TODO: just pass plot_attrs?

    original_labels = trait_data.trait_names

    # L_header =  LOAD_HEADER
    # sv_header = SV_HEADER
    # fac_header = FAC_HEADER
    prep_loadings(log_path, csv_path,
                  hpd_alpha = hpd_alpha,
                  burnin=burnin,
                  scale_loadings_by_factors = scale_loadings_by_factors,
                  original_labels = original_labels)
end

function prep_loadings(log_path::BeastLog, csv_path::String;
                       burnin::Float64 = 0.1,
                       hpd_alpha::Float64 = 0.05,
                       L_header::String = LOAD_HEADER,
                       sv_header::String = SV_HEADER,
                       fac_header::String = FAC_HEADER,
                       scale_loadings_by_factors::Bool = true,
                       original_labels::Vector{<:AbstractString} = String[],
                       k::Int = -1,
                       n_traits::Int = -1,
                       offset::Int = 0)


    cols, data = get_log(log_path, burnin=burnin)

    L_inds = findall(x -> startswith(x, L_header), cols)
    sv_inds = findall(x -> startswith(x, sv_header), cols)
    L_cols = cols[L_inds]
    L_data = @view data[:, L_inds]

    k = k == -1 ? length(sv_inds) : k
    n_traits = n_traits == -1 ? k : n_traits

    p, r = divrem(length(L_cols), k)
    @assert r == 0

    original_labels = isempty(original_labels) ?
                                    ["trait$i" for i = 1:p] : original_labels


    n = size(data, 1)

    df = DataFrame()
    types = [String, Int, Float64, Float64, Float64, Float64]
    nms = ["trait", "factor", "L", "perc", "hpdu", "hpdl"]
    dim = k * p

    for i = 1:length(types)
        df[!, nms[i]] = Vector{types[i]}(undef, dim)
    end


    row_counts = zeros(Int, k)

    stdev_adjustments = ones(n, k)

    if scale_loadings_by_factors
        fac_inds = findall(startswith(fac_header), cols)
        if length(fac_inds) == 0
            error("Cannot find any logged factor values. " *
                "Consider adding a 'fac_header' keywork argument.")
        end

        n_taxa, r = divrem(length(fac_inds), n_traits)

        @assert r == 0
        F_data = @view data[:, fac_inds]

        if k != n_traits
            @warn "The number of factors does not equal the number of traits. " *
                  "Assuming indices $(offset + 1) to $(offset + k) " *
                  "correspond to the factors. " *
                  "Add an 'offset' keywork argument to change this behavior."
        end

        for i = 1:n
            F = reshape(F_data[i, :], n_traits, n_taxa)[1:k, :]
            stds = std(F, dims=2, corrected=false)
            stdev_adjustments[i, :] .= vec(stds)
        end
    end


    for i = 1:k
        for j in 1:p
            ind = (i - 1) * p + j
            @assert L_cols[ind] == "$L_header$i$(j)" || L_cols[ind] == "$L_header.$i.$j"

            df.trait[ind] = original_labels[j]
            df.factor[ind] = i

            vals = @view(L_data[:, ind]) .* @view(stdev_adjustments[:, i])

            μ = mean(vals)
            df.L[ind] = μ

            hpd = hpd_interval(vals, alpha = hpd_alpha)
            df.hpdu[ind] = hpd.upper
            df.hpdl[ind] = hpd.lower

            perc_same = count(x -> sign(μ) * x > 0.0, vals) / n
            df.perc[ind] = perc_same
            if hpd.upper < 0.0 || hpd.lower > 0.0
                row_counts[i] += 1
            end
        end
    end

    keep_rows = findall(x -> x > 1, row_counts)
    k_effective = length(keep_rows)

    CSV.write(csv_path, df)

    return k_effective
end

function load_prep_and_plot(plot_name::String, log_path::BeastLog, stats_path::String;
        labels_path::String = "",
        width_scale::Float64 = 1.0,
        kw_args...)
    prep_loadings(log_path, stats_path; kw_args...)
    load_plot(plot_name, stats_path, labels_path = labels_path, width_scale = width_scale)
end

function load_plot(plot_name::String, statistics::DataFrame; kwargs...)
    tmp_path = "tmp.csv"
    @assert !isfile(tmp_path)
    CSV.write(tmp_path, statistics)
    load_plot(plot_name, tmp_path; kwargs...)
    rm(tmp_path)
end


function load_plot(plot_name::String, statistics_path::String;
        labels_path::String = "",
        width_scale::Float64 = 1.0)
    @rput statistics_path
    @rput plot_name

    plots_path = R_PLOT_SCRIPT

    if isempty(labels_path)
        labels_path = missing
    end

    @rput plots_path
    labels_array = [labels_path] # can't @rput missing directly
    @rput labels_array
    @rput width_scale
    R"""
    source(plots_path)
    plot_loadings(statistics_path, plot_name, labels_path = labels_array[[1]],
        width_scale = width_scale)
    """
end

function prep_factors(svd_path::BeastLog, out_path::String;
                      fac_header::String = FAC_HEADER,
                      sv_header::String = SV_HEADER,
                      k::Int = -1,
                      burnin::Float64 = 0.1)
    cols, data = get_log(svd_path, burnin = burnin)
    if k == -1
        k = length(findall(x -> startswith(x, sv_header), cols))
    end
    fac_inds = findall(x -> startswith(x, fac_header), cols)
    fac_cols = cols[fac_inds]
    fac_means = vec(mean(data[:, fac_inds], dims = 1))
    nk = length(fac_cols)
    n, r = divrem(nk, k)
    @assert r == 0

    F = zeros(n, k)
    taxa = Vector{String}(undef, n)

    ind = 1
    f = length(fac_header) + 1
    for i = 1:n
        s = split(fac_cols[ind][f:end], '.')
        taxon = join(s[1:(end - 1)], '.')

        taxa[i] = taxon
        for j = 1:k
            s = split(fac_cols[ind][f:end], '.')
            @assert join(s[1:(end - 1)], '.') == taxon
            @assert s[end] == "$j"
            F[i, j] = fac_means[ind]
            ind += 1
        end
    end

    df = DataFrame(taxon = taxa)

    for i = 1:k
        df[!, Symbol("f$i")] = @view F[:, i]
    end

    CSV.write(out_path, df)
    return taxa, F
end

function prep_optional_arguments(x::Array)
    return [isempty(y) ? missing : y for y in x]
end

function factor_plot(args...; kwargs...)
    prep_r_factors(args...; kwargs...)
    run_r_factors()
end

function prep_r_factors(plot_path::String, stats_path::String, tree_path::String,
                     class_path::String;
                     fac_names::Vector{String} = String[],
                     layout::String = "rectangular",
                     tip_labels::Bool = true,
                     line_width::Real = 1.0,
                     include_only::AbstractArray{<:AbstractString} = String[],
                     relabel::DataFrame = DataFrame()
                     )
    @rput plot_path
    @rput stats_path
    @rput tree_path
    @rput R_PLOT_SCRIPT

    if isempty(class_path)
        class_path = missing
    end
    class_array = [class_path]
    @rput class_array

    optional_arguments = prep_optional_arguments([fac_names, include_only, relabel])
    @rput optional_arguments


    @rput layout
    @rput tip_labels
    @rput line_width
end

function run_r_factors()

    R"""
    source(R_PLOT_SCRIPT)

    fac_names <- optional_arguments[[1]]
    include_only <- optional_arguments[[2]]
    relabel <- optional_arguments[[3]]


    plot_factor_tree(plot_path, tree_path, stats_path, class_path=class_array[[1]],
                     fac_names = fac_names, layout = layout,
                     tip_labels = tip_labels, line_width = line_width,
                     include_only = include_only, relabel = relabel)
    """
end

function factor_prep_and_plot(plot_path::String, log_path::BeastLog,
                              stats_path::String, tree_path::String;
                              class_path::String = "",
                              check_stats::Bool = false,
                              fac_names::AbstractArray{<:AbstractString} = String[],
                              layout::String = "rectangular",
                              tip_labels::Bool = true,
                              line_width::Real = 1.0,
                              include_only::AbstractArray{<:AbstractString} = String[],
                              relabel::DataFrame = DataFrame(),
                              kwargs...)
    if !(check_stats && isfile(stats_path))
        prep_factors(log_path, stats_path; kwargs...)
    end

    factor_plot(plot_path, stats_path, tree_path, class_path,
                layout = layout, fac_names = fac_names,
                tip_labels = tip_labels,
                line_width = line_width,
                include_only = include_only,
                relabel = relabel)
end


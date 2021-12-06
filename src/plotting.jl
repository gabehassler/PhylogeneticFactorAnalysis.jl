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

function prep_loadings(log_path::String, csv_path::String;
                       burnin::Float64 = 0.1,
                       hpd_alpha::Float64 = 0.05,
                       L_header::String = LOAD_HEADER,
                       sv_header::String = SV_HEADER,
                       fac_header::String = FAC_HEADER,
                       scale_loadings_by_factors::Bool = true,
                       original_labels::Vector{<:AbstractString} = String[],
                       k::Int = -1,
                       n_traits::Int = -1)


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

        n_taxa, r = divrem(length(fac_inds), n_traits)
        @show length(fac_inds)
        @show p
        @assert r == 0
        F_data = @view data[:, fac_inds]

        if k != n_traits
            @warn "The number of factors does not equal the number of traits. " *
                  "Assuming indices 1 to $k correspond to the factors."
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

function load_prep_and_plot(plot_name::String, log_path::String, stats_path::String;
        labels_path::String = "",
        width_scale::Float64 = 1.0,
        kw_args...)
    prep_loadings(log_path, stats_path; kw_args...)
    load_plot(plot_name, stats_path, labels_path = labels_path, width_scale = width_scale)
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

function prep_factors(svd_path::String, out_path::String)
    cols, data = get_log(svd_path)
    k = length(findall(x -> startswith(x, SV_HEADER), cols))
    fac_inds = findall(x -> startswith(x, FAC_HEADER), cols)
    fac_cols = cols[fac_inds]
    fac_means = vec(mean(data[:, fac_inds], dims = 1))
    nk = length(fac_cols)
    n, r = divrem(nk, k)
    @assert r == 0

    F = zeros(n, k)
    taxa = Vector{String}(undef, n)

    ind = 1
    for i = 1:n
        s = split(fac_cols[ind], '.')
        taxon = join(s[2:(end - 1)], '.')
        taxa[i] = taxon
        for j = 1:k
            s = split(fac_cols[ind], '.')
            @assert join(s[2:(end - 1)], '.') == taxon
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

function factor_plot(plot_path::String, stats_path::String, tree_path::String,
                     class_path::String)
    @rput plot_path
    @rput stats_path
    @rput tree_path
    @rput R_PLOT_SCRIPT

    if isempty(class_path)
        class_path = missing
    end
    class_array = [class_path]
    @rput class_array
    R"""
    source(R_PLOT_SCRIPT)
    plot_factor_tree(plot_path, tree_path, stats_path, class_path=class_array[[1]])
    """
end

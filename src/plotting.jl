const R_PLOT_SCRIPT = joinpath(@__DIR__, "plots.R")
const LOAD_HEADER = "L"
const SV_HEADER = "sv"
const FAC_HEADER = "factors."


function prep_loadings(input::PipelineInput, log_path::String,
                            csv_path::String)

    @unpack plot_attrs, data = input
    @unpack trait_data = data
    @unpack burnin, hpd_alpha, scale_loadings_by_factors = plot_attrs #TODO: just pass plot_attrs?

    cat_dict = Dict{String, String}()
    original_labels = trait_data.trait_names

    cols, data = get_log(log_path, burnin=burnin)
    L_header =  LOAD_HEADER
    sv_header = SV_HEADER
    fac_header = FAC_HEADER

    L_inds = findall(x -> startswith(x, L_header), cols)
    sv_inds = findall(x -> startswith(x, sv_header), cols)
    L_cols = cols[L_inds]
    L_data = @view data[:, L_inds]

    k = length(sv_inds)
    p, r = divrem(length(L_cols), k)
    @assert r == 0

    n = size(data, 1)

    df = DataFrame([String, Int, Float64, Float64, Float64, Float64],
                    ["trait", "factor", "L", "perc", "hpdu", "hpdl"],
                    k * p)

    row_counts = zeros(Int, k)

    stdev_adjustments = ones(n, k)

    if scale_loadings_by_factors
        fac_inds = findall(startswith(fac_header), cols)
        n_taxa, r = divrem(length(fac_inds), k)
        @assert r == 0
        F_cols = @view cols[fac_inds]
        F_data = @view data[:, fac_inds]

        for i = 1:n
            F = reshape(F_data[i, :], k, n_taxa)
            stds = std(F, dims=2, corrected=false)
            stdev_adjustments[i, :] .= vec(stds)
        end
    end


    for i = 1:k
        for j in 1:p
            ind = (i - 1) * p + j
            @assert L_cols[ind] == "$L_header$i$(j)"

            df.trait[ind] = original_labels[j]
            df.factor[ind] = i

            vals = @view(L_data[:, ind]) .* @view(stdev_adjustments[:, i])
            μ = mean(vals)
            df.L[ind] = μ

            hpd = hpd_interval(vals)
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

function load_plot(plot_name::String, statistics_path::String; labels_path::String = "")
    @rput statistics_path
    @rput plot_name

    plots_path = R_PLOT_SCRIPT

    if isempty(labels_path)
        labels_path = missing
    end

    @rput plots_path
    labels_array = [labels_path] # can't @rput missing directly
    @rput labels_array
    R"""
    source(plots_path)
    plot_loadings(statistics_path, plot_name, labels_path = labels_array[[1]])
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
        taxon = split(fac_cols[ind], '.')[2]
        taxa[i] = taxon
        for j = 1:k
            s = split(fac_cols[ind], '.')
            @assert s[2] == taxon
            @assert s[3] == "$j"
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

function factor_plot(plot_path::String, stats_path::String, tree_path::String)
    @rput plot_path
    @rput stats_path
    @rput tree_path
    @rput R_PLOT_SCRIPT
    R"""
    source(R_PLOT_SCRIPT)
    plot_factor_tree(plot_path, tree_path, stats_path)
    """
end

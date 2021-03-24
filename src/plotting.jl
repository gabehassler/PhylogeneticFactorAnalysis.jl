function prep_for_plotting(input::PipelineInput, log_path::String,
                            csv_path::String)

    @unpack plot_attrs, labels_path, trait_data = input
    @unpack burnin, hpd_alpha, scale_loadings_by_factors = plot_attrs #TODO: just pass plot_attrs?

    cat_dict = Dict{String, String}()
    original_labels = trait_data.trait_names

    if isempty(labels_path)
        for label in original_labels
            cat_dict[label] = "NA"
        end
    else
        labels_df = CSV.read(labels_path, DataFrame)
        for label in original_labels
            ind = findall(isequal(label), labels_df.label)
            @assert length(ind) == 1
            ind = ind[1]

            cat_dict[label] = df.category[ind]
        end
    end



    cols, data = get_log(log_path, burnin=burnin)
    L_header = "L"
    sv_header = "sv"
    fac_header = "factors."

    L_inds = findall(x -> startswith(x, L_header), cols)
    sv_inds = findall(x -> startswith(x, sv_header), cols)
    L_cols = cols[L_inds]
    L_data = @view data[:, L_inds]

    k = length(sv_inds)
    p, r = divrem(length(L_cols), k)
    @assert r == 0

    n = size(data, 1)

    df = DataFrame([String, String, Int, Float64, Float64, Float64, Float64],
                    ["trait", "cat", "factor", "L", "perc", "hpdu", "hpdl"],
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
            df.cat[ind] = cat_dict[original_labels[j]]
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
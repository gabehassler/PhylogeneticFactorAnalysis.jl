function plot_loadings(plot_paths::Vector{String},
                       log_path::BeastLog,
                       trait_names::Vector{String},
                       trait_dims::Vector{Int};
                       original_labels::Vector{<:Vector{<:AbstractString}} = Vector{String}[],
                       factor_partitions::AbstractVector{Int} = 1:(length(plot_paths)),
                       kwargs...)

    offset = 0
    n = length(trait_names)
    joint_name = join(trait_names, '.') * ".joint"
    n_traits = sum(trait_dims)
    for i in factor_partitions
        trait_name = trait_names[i]
        labels = isempty(original_labels) ? String[] : original_labels[i]
        plot_loadings(plot_paths[i], log_path,
              L_header = "$trait_name.L.", fac_header = joint_name,
              prec_header = "$trait_name.factorPrecision",
              k = trait_dims[i], n_traits = n_traits,
              offset = offset,
              original_labels = labels;
              kwargs...)
        offset += trait_dims[i]
    end
end

function plot_loadings(log_path::BeastLog,
                       trait_names::Vector{String},
                       trait_dims::Vector{Int};
                       dir::String = pwd(),
                       file_base::String = "",
                       factor_partitions::AbstractVector{Int} = 1:(length(trait_names)),
                       kwargs...)
    nm = isempty(file_base) ? "" : file_base * "_"
    plot_paths = [joinpath(dir, "$(nm)$(x)_loadings.pdf") for x in trait_names]
    plot_loadings(plot_paths, log_path, trait_names, trait_dims, factor_partitions = factor_partitions; kwargs...)
end

function plot_loadings(plot_path::String, log_path::BeastLog; kwargs...)
    max_i = 1000
    stats_path = "tmp.csv"
    i = 0
    while isfile(stats_path) && i < max_i
        i += 1
        stats_path = "tmp$i.csv"
    end

    if i == max_i
        error("Cannot find unused temporary file name. " *
            "Consider deleting 'tmp<i>.csv' files.")
    end



    load_prep_and_plot(plot_path, log_path, stats_path; kwargs...)
    rm(stats_path)
end


function process_variance(log_path::BeastLog; var_start::String = "mbd.variance")
    df = import_log(log_path, burnin = 0.1)
    var_data = Matrix(df[!, startswith.(names(df), Ref(var_start))])

    n, p2 = size(var_data)
    p = Int(sqrt(p2))

    cor_data = zeros(n, p2)

    for i = 1:n
        V = reshape(var_data[i, :], p, p)
        C = cov2corr(V)
        cor_data[i, :] .= vec(C)
    end

    return (V = process_mat(var_data), C = process_mat(cor_data))
end


function process_mat(X::Matrix{Float64}, d1::Int, d2::Int)
    μ = mean(X, dims = 1)
    n, p = size(X)
    @assert p == d1 * d2

    probs = [count(sign.(X[:, i]) .== sign(μ[i])) / n for i = 1:p]
    hpds = hpd_interval(X)
    hpdls = [hpds[i].lower for i = 1:p]
    hpdus = [hpds[i].upper for i = 1:p]

    return (μ = reshape(μ, d1, d2),
            prob = reshape(probs, d1, d2),
            hpdl = reshape(hpdls, d1, d2),
            hpdu = reshape(hpdus, d1, d2))
end

function process_mat(X::Matrix{Float64})
    n, p = size(X)
    d = isqrt(p)
    return process_mat(X, d, d)
end


# struct CorrPlotParameters
#     factor_inds::Vector{Vector{Int}}
#     optimization_inds::Vector{Vector{Int}}
# end

function plot_correlation(path::String,
                          corrs::Matrix{Float64},
                          hpdl::Matrix{Float64},
                          hpdu::Matrix{Float64},
                          labels::Vector{String};
                          width::Real = 1 + 0.75 * size(corrs, 1),
                          annotations::Union{Nothing, DataFrame} = nothing
                          )

    @rput corrs
    @rput hpdl
    @rput hpdu
    @rput labels
    @rput path
    @rput width

    @rput annotations

    R"""
    library(ggplot2)
    library(dplyr)
    library(ggforce)
    library(wesanderson)



    prep_cor <- function(C) {
        # Convert the correlation matrix to a long format
        levs = rownames(C)
        C_long <- as.data.frame(as.table(C))
        colnames(C_long) <- c("Var1", "Var2", "value")
        C_long$Var1 = factor(C_long$Var1, levels = levs, ordered = TRUE)
        C_long$Var2 = factor(C_long$Var2, levels = levs, ordered = TRUE)
        C_long
    }



    # annotations = data.frame(
    #     ind1 = c(1, 3, 4),
    #     ind2 = c(2, 10, 5),
    #     rank = c(1, 2, 3)
    # )

    annotate <- function(var1, var2, annotations) {
        var1 = as.numeric(var1)
        var2 = as.numeric(var2)
        inds = which(annotations$ind1 == var2 & annotations$ind2 == var1)
        if (length(inds) > 0) {
            return(annotations$rank[inds])
        } else {
            return(NA)
        }
    }

    prep_cor2 <- function(cdf, inds = "all", annotations = NULL) {
        to_offset = c("value", "lower", "upper")
        cdf = cdf %>%
            mutate(Var2r = factor(Var2, levels = rev(levels(Var2)))) %>%
            mutate(v1 = as.numeric(Var1),
                v2 = as.numeric(Var2r),
                v1a = v1,
                v2a = v2,
                v1b = v1,
                v2b = v2,
                v1_center = ifelse(Var1 == Var2, v1, NA),
                v2_center = ifelse(Var1 == Var2, v2, NA),
                v1_top = ifelse(Var2 == levels(Var2)[1], v1, NA),
                v2_top = ifelse(Var2 == levels(Var2)[1] & Var1 != levels(Var1)[1], v2, NA) + 1) %>%
            mutate(across(all_of(to_offset), ~ .x * 0.5 + v2, .names = "{col}_shift")) %>%
            mutate(y0 = v2,
                x0start = v1 - 0.5,
                x0end = x0start + 1)

        to_null_strict = c(paste0(to_offset, "_shift"), "y0", "x0start", "x0end", "v1a", "v2a")
        to_null = c("v1b", "v2b")
        if (inds == "upper") {
            cdf = cdf %>%
                mutate(across(all_of(to_null_strict), ~ ifelse(Var1 <= Var2, NA, .x))) %>%
                mutate(across(all_of(to_null), ~ ifelse(Var1 < Var2, NA, .x)))
        } else if (inds == "lower") {
            cdf = cdf %>%
                mutate(across(all_of(to_null_strict), ~ ifelse(Var1 >= Var2, NA, .x))) %>%
                mutate(across(all_of(to_null), ~ ifelse(Var1 > Var2, NA, .x)))
        } else if (inds == "all") {
            cdf = cdf %>%
                mutate(across(all_of(to_null_strict), ~ ifelse(Var1 == Var2, NA, .x)))
        }

        if (is.null(annotations)) {
            cdf = cdf %>%
                mutate(marker = NA)
        } else {
            cdf = cdf %>%
                rowwise() %>%
                mutate(marker = annotate(Var1, Var2, annotations))
        }

        return(cdf)
    }


    wpal <- wes_palette("Royal1")

    pal <- c(wpal[1], "white", wpal[2])
    rownames(corrs) <- labels
    colnames(corrs) <- labels
    rownames(hpdl) <- labels
    colnames(hpdl) <- labels
    rownames(hpdu) <- labels
    colnames(hpdu) <- labels

    C = prep_cor(corrs)
    L = prep_cor(hpdl)
    U = prep_cor(hpdu)
    C$lower = L$value
    C$upper = U$value


    # cdf = prep_cor(C) %>%
    #     # dplyr::mutate(Var1 = ifelse(Var1 <= Var2, NA, Var1)) %>%
    #     rowwise() %>%
    #     mutate(marker = annotate(Var1, Var2, annotations),
    #         lower = ifelse(value == 1, 1, value - 0.3 * abs(runif(1))),
    #         upper = ifelse(value == 1, 1, value + 0.3 * abs(runif(1)))
    #     )


    print(annotations)

    C = prep_cor2(C, inds = "upper", annotations = annotations)

    plt = C %>%
        ggplot(aes(x = v1, y = v2)) +
        # geom_point(aes(size = 200 * value)) +
        # ggforce::geom_circle(aes(r = abs(value) * 0.5, fill = value, x0 = as.numeric(Var1), y0 = as.numeric(Var2))) +
        geom_crossbar(aes(y = value_shift, ymin = lower_shift, ymax = upper_shift, fill = value), width = 0.5, linewidth=0.2) +
        coord_fixed() +
        geom_tile(color = "black", fill = NA, aes(x = v1a, y = v2a)) +
        geom_text(aes(x = v1_center + 0.4, y = v2_center, label = Var1), hjust = 1) +
        geom_text(aes(x = v1_top, y = v2_top - 0.4, label = Var1), angle = 45, hjust = 0, vjust = 0) +
        geom_text(aes(x = v1 + 0.35, y = v2 + 0.35, label = marker)) +
        geom_segment(aes(x = x0start, y = y0, xend = x0end, yend = y0), color = "black", linetype="dashed") +
        theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(angle = 0, hjust = 1)) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
        theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank()) +
        labs(fill = "Correlation") +
        ylim(1, length(labels) + 2)

    ggsave(path, plot = plt, width = width, height = width, units = "in", dpi = 300)
    # library(corrplot)
    # library(wesanderson)
    # wpal <- wes_palette("Royal1")

    # pdf(path, width=width, height=width)
    # corrplot(corrs, low=hpdl, upp=hpdu, plot="rect", type="upper", diag=FALSE,
    #          tl.col="black",
    #          col=colorRampPalette(pal)(200))
    # dev.off()
    """
end

function plot_correlation(bn::String, labels::Vector{String}; kwargs...)
    return plot_correlation("$bn.log", "$bn.pdf", labels)
end

function plot_correlation(log_path::BeastLog, plot_path::String, labels::Vector{String}; kwargs...)
    @unpack C = process_variance(log_path)
    @unpack μ, hpdl, hpdu = C
    plot_correlation(plot_path, μ, hpdu, hpdl, labels; kwargs...)
end


function cor_mat(M1::AbstractMatrix{Float64}, M2::AbstractMatrix{Float64})
    k = size(M1, 2)
    C = zeros(k, k)
    for i = 1:k
        for j = 1:k
            C[i, j] = cor(M1[:, i], M2[:, j])
        end
    end
    return C
end
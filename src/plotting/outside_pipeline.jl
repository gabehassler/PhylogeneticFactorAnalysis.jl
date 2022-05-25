function plot_loadings(plot_paths::Vector{String},
                       log_path::BeastLog,
                       trait_names::Vector{String},
                       trait_dims::Vector{Int};
                       factor_partitions::AbstractVector{Int} = 1:(length(plot_paths)),
                       kwargs...)

    offset = 0
    n = length(trait_names)
    joint_name = join(trait_names, '.') * ".joint"
    n_traits = sum(trait_dims)
    for i in factor_partitions
        trait_name = trait_names[i]
        plot_loadings(plot_paths[i], log_path,
              L_header = "$trait_name.L.", fac_header = joint_name,
              k = trait_dims[i], n_traits = n_traits,
              offset = offset;
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

function plot_correlation(path::String,
                          corrs::Matrix{Float64},
                          hpdl::Matrix{Float64},
                          hpdu::Matrix{Float64},
                          labels::Vector{String})

    @rput corrs
    @rput hpdl
    @rput hpdu
    @rput labels
    @rput path

    R"""
    library(corrplot)
    rownames(corrs) <- labels
    colnames(corrs) <- labels
    rownames(hpdl) <- labels
    colnames(hpdl) <- labels
    rownames(hpdu) <- labels
    colnames(hpdu) <- labels

    pdf(path)
    corrplot(corrs, low=hpdl, upp=hpdu, plot="rect")
    dev.off()
    """
end

function plot_correlation(bn::String, labels::Vector{String})
    return plot_correlation("$bn.log", "$bn.pdf", labels)
end

function plot_correlation(log_path::BeastLog, plot_path::String, labels::Vector{String})
    @unpack C = process_variance(log_path)
    @unpack μ, hpdl, hpdu = C
    plot_correlation(plot_path, μ, hpdu, hpdl, labels)
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
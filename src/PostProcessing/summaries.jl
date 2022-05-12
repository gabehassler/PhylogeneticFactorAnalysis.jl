
# required functions:
# distance(value::AbstractVector{Float64}, cluster::Int, summary::ClusterStatistics)
# constructor(data::Array{Float64, 3}, reference::Int)

struct MultivariateNormal <: ClusterStatistics
    means::Vector{Vector{Float64}}
    precisions::Vector{Matrix{Float64}}
    logdeterminants::Vector{Float64}

    function MultivariateNormal(means::Vector{Vector{Float64}},
                                          precisions::Vector{Matrix{Float64}},
                                          logdeterminants::Vector{Float64})
        @assert length(means) == length(precisions) == length(logdeterminants)
        return new(means, precisions, logdeterminants)
    end
end

function MultivariateNormal(means::Vector{Vector{Float64}},
                                      precisions::Vector{Matrix{Float64}})
    logdets = [-logdet(p) for p in precisions]
    return MultivariateNormal(means, precisions, logdets)
end

function MultivariateNormal(data::Array{Float64, 3},
                                      reference::Int)
    p, k = size(data, (1, 2))
    means = [data[:, i, reference] for i = 1:k]
    precisions = [Matrix(Diagonal(ones(p))) for _ = 1:k]
    logdets = zeros(k)

    return MultivariateNormal(means, precisions, logdets)
end


function distance(value::AbstractVector{Float64},
                  cluster::Int,
                  summary::MultivariateNormal)

    @unpack means, precisions, logdeterminants = summary

    errs = value - means[cluster]

    return errs' * precisions[cluster] * errs - logdeterminants[cluster]
end

function clear_stats!(stats::ClusterStatistics)
    k = length(stats.means)
    for i = 1:k
        fill!(stats.means[i], 0.0)
        fill!(stats.precisions[i], 0.0)
    end
    fill!(stats.logdeterminants, 0.0)
end

################################################################################
## IndependentNormal
################################################################################

struct IndependentNormal <: ClusterStatistics
    means::Vector{Vector{Float64}}
    precisions::Vector{Vector{Float64}}
    logdeterminants::Vector{Float64}
end

function IndependentNormal(data::Array{Float64, 3};
                           reference::Int = size(data, 3))
    p, k, _ = size(data)
    means = [data[:, i, reference] for i = 1:k]
    precisions = [ones(p) for _ = 1:k]
    logdets = zeros(k)

    return IndependentNormal(means, precisions, logdets)
end


function distance(value::AbstractVector{Float64},
                  cluster::Int,
                  summary::IndependentNormal)

    @unpack means, precisions, logdeterminants = summary
    errs = value - means[cluster]
    return errs' * Diagonal(precisions[cluster]) * errs - logdeterminants[cluster]
end

function recompute_statistics!(assignments::ClusterAssignments,
                               data::AbstractArray{Float64, 3},
                               stats::IndependentNormal)

    p, k, n = size(data)
    clear_stats!(stats)
    @unpack means, precisions, logdeterminants = stats

    for i = 1:n
        for j = 1:k
            y = value_at_cluster(data, assignments, i, j)
            means[j] += y
            precisions[j] += y.* y
        end
    end

    for j = 1:k
        precisions[j] .-= means[j]
        means[j] ./= n

        precisions[j] .= n ./ precisions[j]
        logdeterminants[j] = log(prod(precisions[j]))
    end
end

################################################################################
## Geodesics
################################################################################

struct Geodesic <: ClusterStatistics
    means::Vector{Vector{Float64}}
end

function Geodesic(data::AbstractArray{Float64, 3}, reference::Int)
    return Geodesic(@view data[:, :, reference])
end

function Geodesic(data::AbstractMatrix{Float64})
    p, k = size(data)
    means = [zeros(p) for _ = 1:k]
    for i = 1:k
        y = data[:, i] + 1e-6 * randn(p) #need to add a little jitter
        means[i] .= y / norm(y)
    end
    return Geodesic(means)
end

function distance(value::AbstractVector{Float64}, cluster::Int, stats::Geodesic)
    normalized_value = value ./ norm(value)
    # @show norm(normalized_value)
    # @show norm(stats.means[cluster])
    return acos(dot(normalized_value, stats.means[cluster]))
end

function recompute_statistics!(assignments::ClusterAssignments,
                               data::AbstractArray{Float64, 3},
                               stats::Geodesic)
    p, k,n = size(data)
    @unpack means = stats
    for i = 1:k
        fill!(means[i], 0.0)
    end

    for i = 1:n
        for j = 1:k
            y = value_at_cluster(data, assignments, i, j)
            means[j] .+= y ./ norm(y)
        end
    end

    for j = 1:k
        means[j] ./= norm(means[j])
    end
end


################################################################################
## Procrustes
################################################################################

struct Procrustes <: ClusterStatistics
    means::Matrix{Float64}
    # TODO
end

function Procrustes(data::Array{Float64, 3},
                    reference::Int)

    p, k, n = size(data)
    X = reshape(mean(data, dims=3), p, k)

    return Procrustes(X)
end



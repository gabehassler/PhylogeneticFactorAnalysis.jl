abstract type ClusterStatistics <: AbstractTransformer end
function compute_transform!(t::Type{<:ClusterStatistics}, X::Array{Float64, 3},
                            reference_ind::Int)
    return cluster!(X, t(X[:, :, reference_ind]))
end

struct ClusterAssignments <: Transform
    clusters::Matrix{Int}
    signs::Matrix{Int}
    # data::Array{Float64, 3}
    # means::Matrix{Float64}
    # precisions::Vector{<:AbstractArray{Float64, 2}}
    # statistics::ClusterStatistics


    function ClusterAssignments(::Int, k::Int, n::Int)
        clusters = collect(1:k) * ones(n)'
        signs = ones(k, n)

        return new(clusters, signs)
    end
end

function reassign_clusters!(assignments::ClusterAssignments,
                            data::Array{Float64, 3},
                            statistics::ClusterStatistics)
    @unpack clusters, signs = assignments

    p, k, n = size(data)
    perms = permutations(1:k)
    min_perm = zeros(Int, k)
    min_signs = zeros(Int, k)
    current_signs = zeros(Int, k)

    n_updated = 0


    M = zeros(k, k)
    S = zeros(Int, k, k)

    for i = 1:n

        for clust = 1:k
            for point = 1:k
                y = @view data[:, point, i]

                dist_pos = distance(y, clust, statistics)
                dist_neg = distance(-1 * y, clust, statistics)

                dist = 0.0
                sgn = 0

                if dist_pos > dist_neg
                    dist = dist_neg
                    sgn = -1
                else
                    dist = dist_pos
                    sgn = 1
                end

                M[point, clust] = dist
                S[point, clust] = sgn
            end
        end

        min_ss = Inf

        for perm in perms
            ss = 0.0
            for point = 1:k
                ss += M[perm[point], point]
                current_signs[point] = S[perm[point], point]
            end


            if ss < min_ss
                min_ss = ss
                min_perm .= perm
                min_signs .= current_signs
            end
        end


        if min_perm != clusters[:, i] || min_signs != signs[:, i]
            n_updated += 1
        end
        clusters[:, i] .= min_perm
        signs[:, i] .= min_signs
    end
    return n_updated
end

function cluster(data::Array{Float64, 3}, statistics::ClusterStatistics;
                 max_iterations::Int = 100
                 )


    clusters = ClusterAssignments(size(data)...)
    ind = 0
    n_updated = 1

    while ind < max_iterations && n_updated != 0
        ind += 1
        @show ind
        n_updated = reassign_clusters!(clusters, data, statistics)
        @show n_updated
        recompute_statistics!(clusters, data, statistics)

    end
    return clusters
end

function cluster!(data::Array{Float64, 3}, statistics::ClusterStatistics; kwargs...)
    clusters = cluster(data, statistics; kwargs...)
    data .= apply_transform(data, clusters)
    return clusters
end


function apply_transform(data::AbstractArray{Float64, 3},
                        assignments::ClusterAssignments)
    @unpack clusters, signs = assignments

    arranged_data = copy(data)
    p, k, n = size(data)
    for i = 1:n
        arranged_data[:, :, i] .= data[:, clusters[:, i], i] * Diagonal(signs[:, i])
    end
    return arranged_data
end
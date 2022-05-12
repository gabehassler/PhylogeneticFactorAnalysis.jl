function value_at_cluster(data::AbstractArray{Float64, 3},
                          assignments::ClusterAssignments,
                          state::Int,
                          clust::Int;
                          buffer::Vector{Float64} = zeros(size(data, 1)))

    @unpack clusters, signs = assignments
    buffer .= data[:, clusters[clust, state], state]
    buffer .*= signs[clust, state] # TODO: clusters[clust]?

    return buffer
end

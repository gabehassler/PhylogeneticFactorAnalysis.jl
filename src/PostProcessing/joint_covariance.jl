function induce_block_independence(V::AbstractArray{Float64},
                                   dims::AbstractVector{Int}...)
    n = length(dims)
    T = diagm(0 => ones(size(V, 1)))

    for i = 1:n
        dimsi = dims[i]
        V_sub = V[dimsi, dimsi]
        chol = cholesky(V_sub)
        T[dimsi, dimsi] .= inv(chol.U)
    end

    return T
end

function induce_block_independence!(V::AbstractArray{Float64},
                                    dims::AbstractVector{Int}...)
    T = induce_block_independence(V, dims...)
    V .= T' * V * T #TODO: do in block form
    return T
end






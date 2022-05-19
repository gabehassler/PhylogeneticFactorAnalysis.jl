function procrustes(data::Array{Float64, 3};
                    max_iterations::Int = 20,
                    tol::Float64 = 1e-8)
    Y = copy(data)
    p, k, n = size(data)
    max_diff = Inf
    counter = 0
    # H = 10 # TODO: set to N?
    J = Diagonal(ones(p))
    # J = I - 1/H * ones(p, p)

    X_bar = reshape(mean(Y, dims=3), p, k)

    while max_diff > tol && counter < max_iterations
        rotate_to(Y, X_bar, J = J)

        X_bar_old = X_bar
        X_bar = reshape(mean(Y, dims=3), p, k)
        max_diff = maximum(abs.(X_bar - X_bar_old))

        counter += 1
        @show max_diff
        @show counter
    end

    return Y
end

function rotate_to(X::AbstractMatrix{Float64},
                   reference::AbstractMatrix{Float64};
                   J::AbstractMatrix{Float64} = I)
    C = X' * J * reference
    s = svd(C)
    H = s.U * s.Vt
    @assert isortho(H)
    return X * H
end

function rotate_to(X::AbstractArray{Float64, 3},
                   reference::AbstractMatrix{Float64}; kwargs...)
    for i = 1:size(X, 3)
        X_sub = @view X[:, :, i]
        X_sub .= rotate_to(X_sub, reference; kwargs...)
    end
end



function absmax(x::AbstractArray{T}) where T <: Real
    return maximum(abs.(x))
end

function isortho(x::AbstractMatrix{Float64})
    return absmax(I - x * x') < 1e-10
end


################################################################################
## Procrustes
################################################################################

struct ProcrustesRotation <: AbstractRotation
    reference::AbstractArray{Float64, 2}
end

function ProcrustesRotation(data::AbstractArray{Float64, 3})
    p, k, _ = size(data)
    return ProcrustesRotation(reshape(mean(data, dims=3), p, k))
end

function update_rotation!(rotations::Rotations, data::Array{Float64, 3},
                          procrustes::ProcrustesRotation)
    p, k, n = size(data)

    for i = 1:n
        R = view_rotation(rotations, i)
        Y = @view data[:, :, i]
        Yr = Y * R # TODO: cache

        C = Yr' * procrustes.reference
        s = svd(C)
        H = s.U * s.Vt
        update_rotation!(rotations, H, i)
    end
end

function update_statistics!(rotations::Rotations, data::Array{Float64, 3},
                            procrustes::ProcrustesRotation)
    p, k, n = size(data)

    Y_bar = zeros(p, k)

    for i = 1:n
        Y_bar .+=  rotate_sample(data, rotations, i) # TODO: use cache
    end

    Y_bar ./= n

    diff = absmax(Y_bar - procrustes.reference)
    procrustes.reference .= Y_bar

    return diff
end


################################################################################
## SVD
################################################################################

struct SVDRotation <: AbstractRotation
    # doesn't need to store any statistics
end

function SVDRotation(::Matrix{Float64})
    return SVDRotation()
end

function update_rotation!(rotations::Rotations, data::Array{Float64, 3},
                          ::SVDRotation)

    p, k, n = size(data)

    for i = 1:n
        Y = rotate_sample(data, rotations, i)
        s = svd(Y)
        update_rotation!(rotations, s.Vt', i)
    end
end


################################################################################
## Procrustes
################################################################################

struct ProcrustesRotation <: AbstractRotation
    reference::AbstractArray{Float64, 2}
    dims::Vector{Int}
end

function ProcrustesRotation(data::AbstractArray{Float64, 2})
    k = size(data, 2)
    return ProcrustesRotation(data, collect(1:k))
end

function ProcrustesRotation(data::AbstractArray{Float64, 3})
    p, k, _ = size(data)
    dims = collect(1:k)
    return ProcrustesRotation(reshape(mean(data, dims=3), p, k), dims)
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
    dims::Vector{Int}
end

function SVDRotation(X::Matrix{Float64})
    k = size(X, 2)
    return SVDRotation(collect(1:k))
end

function update_rotation!(rotations::Rotations, data::Array{Float64, 3},
                          svd_rot::SVDRotation)

    p, k, n = size(data)
    @unpack dims = svd_rot
    Y_dims = zeros(p, length(dims))

    for i = 1:n
        Y = rotate_sample(data, rotations, i)
        Y_dims .= @view Y[:, dims]
        s = svd(Y_dims)
        R_dims = s.Vt'
        if length(dims) == k
            R = R_dims
        else
            R = Matrix(Diagonal(ones(k)))
            R[dims, dims] .= R_dims
        end
        update_rotation!(rotations, R, i)
    end
end

################################################################################
## Signs
################################################################################

struct SignRotation <: AbstractRotation
    dims::Vector{Int}
end

function SignRotation(X::Matrix{Float64})
    k = size(X, 2)
    return SignRotation(collect(1:k))
end

function update_rotation!(rotations::Rotations, data::Array{Float64, 3},
                          sign_rot::SignRotation)

    p, k, n = size(data)
    @unpack dims = sign_rot


    means = zeros(p, k)
    sum_squares = zeros(p, k)


    for i = 1:n #not using built-in mean/var for memory efficiency w/ abs value
        Y = rotate_sample(data, rotations, i) # TODO: buffer this?
        for j in dims
            for l = 1:p
                x = abs(Y[l, j])
                means[l, j] += x
                sum_squares[l, j] += x * x
            end
        end
    end

    means ./= n
    sum_squares ./= n
    zs = zeros(p, k)
    for j in dims
        for l = 1:p
            var = sum_squares[l, j] - means[l, j] * means[l, j]
            z = 0.0
            if means[l, j] > 0
                z = means[l, j] / sqrt(var)
            end
            zs[l, j] = z
        end
    end

    ref_inds = zeros(Int, k)
    for i = 1:k
        ref_inds[i] = findmax(@view zs[:, i])[2]
    end

    for i = 1:n
        r = ones(k)
        Y = rotate_sample(data, rotations, i)
        for j in dims
            if Y[ref_inds[j], j] < 0
                r[j] = -1
            end
        end

        update_rotation!(rotations, Diagonal(r), i)
    end
end

function minimal_variance(R::Matrix{Float64}, X::Array{Float64, 3})
    k, p, n = size(X)
    X_rotated = naive_rotate(R, X) #TODO: use buffer + sum_squares instead
    return -sum(var(X_rotated, dims=3))
end

function four_norm(R::Matrix{Float64}, X::Array{Float64, 3})
    X_rotated = naive_rotate(R, X) #TODO: use buffer + sum_squares instead
    μ = mean(X_rotated, dims=3)
    return norm(μ, 100)
end

function max_significance(R::Matrix{Float64}, X::Array{Float64, 3})
    X_rotated = naive_rotate(R, X) #TODO: use buffer + sum_squares instead
    k, p, n = size(X)

    percs = [count(x -> x > 0, X_rotated[i, j, :]) / n for i in 1:k, j in 1:p]

    for i = 1:length(percs)
        if percs[i] < 0.5
            percs[i] = 1.0 - percs[i]
        end
    end
    return maximum(percs)
end

function total_significance(R::Matrix{Float64}, X::Array{Float64, 3})
    X_rotated = naive_rotate(R, X) #TODO: use buffer + sum_squares instead
    k, p, n = size(X)

    percs = [count(x -> x > 0, X_rotated[i, j, :]) / n for i in 1:k, j in 1:p]

    for i = 1:length(percs)
        if percs[i] < 0.5
            percs[i] = 1.0 - percs[i]
        end
    end
    return sum(percs)
end

function maximum_correlation(R::Matrix{Float64}, X::Array{Float64, 3})
    X_rotated = naive_rotate(R, X) #TODO: use buffer + sum_squares instead
    k, p, n = size(X)
    # @show k, p
    # @assert k == p

    μ = abs.(mean(X_rotated, dims=3))[:, :, 1]
    # for i = 1:k
    #     μ[i, i] = 0.0
    # end

    return maximum(μ)
end






function naive_rotate(R::Matrix{Float64}, X::Array{Float64, 3})
    k, p, n = size(X)
    X_rotated = zeros(k, p, n) #TODO: use buffer + sum_squares instead
    for i = 1:n
        X_rotated[:, :, i] .= R * @view X[:, :, i]
    end
    return X_rotated
end


function optimize_2d(C::Array{Float64, 3}, f::Function; n_steps::Int = 1000) #TODO: this is bad, but easy
    k, q, n = size(C)
    @assert k == 2

    R = zeros(k, k)
    θs = collect(0:(π / (2 * n_steps)):(π / 2)) # should be symmetric across π / 4.

    n = length(θs)
    x = zeros(n)
    for i = 1:n
        fill_rotation!(R, θs[i])
        x[i] = f(R, C)
    end

    max_ind = findmax(x)[2]
    fill_rotation!(R, θs[max_ind])

    return R
end

function optimize_3d(C::Array{Float64, 3}, f::Function; n_steps::Int = 20) #TODO: this is bad, but easy
    k, q, n = size(C)
    @assert k == 3

    R = zeros(k, k)
    θs = collect(0:(π / (2 * n_steps)):(π / 2)) # should be symmetric across π / 4.

    n = length(θs)
    x = zeros(n^3)

    max_value = -Inf
    θ_max = (NaN, NaN, NaN)
    ind = 1
    for i = 1:n
        for j in 1:n
            for k in 1:n
                θ = (θs[i], θs[j], θs[k])
                fill_rotation!(R, θ)
                value = f(R, C)
                if value > max_value
                    max_value = value
                    θ_max = θ
                end
                ind += 1
            end
        end
    end

    fill_rotation!(R, θ_max)

    return R
end

function optimize(C::Array{Float64, 3}, f::Function; kwargs...)
    k = size(C, 1)
    if k == 1
        return ones(1, 1)
    elseif k == 2
        return optimize_2d(C, f; kwargs...)
    elseif k == 3
        return optimize_3d(C, f; kwargs...)
    end
    error("Not implemented for k = $k")
end

function fill_rotation!(R::AbstractMatrix{Float64}, θ::Float64)
    sinθ = sin(θ)
    cosθ = cos(θ)

    R[1, 1] = cosθ
    R[1, 2] = -sinθ
    R[2, 1] = sinθ
    R[2, 2] = cosθ
end

function fill_rotation!(R::Matrix{Float64}, θ::Tuple{Float64, Float64, Float64})
    Rs = [basic_3d(θ[i], i) for i = 1:3] # TODO: very memory inefficient
    R .= Rs[1] * Rs[2] * Rs[3]
end

function basic_3d(θ::Float64, fixed::Int)
    R = zeros(3, 3)
    R[fixed, fixed] = 1.0
    inds = setdiff(1:3, fixed)
    R_sub = @view R[inds, inds]
    fill_rotation!(R_sub, θ)
    return R
end



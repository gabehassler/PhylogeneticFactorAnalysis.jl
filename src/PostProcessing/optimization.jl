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

# function rotated_corr(Ri::AbstractVector{Float64}, Σrr::AbstractMatrix{Float64}, Σrj::AbstractVector{Float64})

#     Σ̂ii = Ri' * Σrr * Ri
#     Σ̂ij = Ri' * Σrj
#     return Σ̂ij / sqrt(Σ̂ii)
# end

# function rotated_corr(Ri::AbstractVector{Float64}, Σrr::AbstractArray{Float64, 3}, Σrj::AbstractMatrix{Float64})
#     k, _, n = size(Σrr)
#     sum = zeros(k)
#     for i = 1:n
#         sum .+= rotated_corr(Ri, @view(Σrr[:, :, i]), @view(Σrj[:, i]))
#     end
#     return sum / n
# end
# function rotated_corr_single(Ri, Σrr, Σrj, Σjj)
#     # @show size(Ri)
#     # @show Σrr
#     # @show Σrj
#     Σ̂ii = Ri' * Σrr * Ri
#     Σ̂ij = Ri' * Σrj
#     return Σ̂ij / sqrt(Σ̂ii * Σjj)
# end

# function rotated_corr(Ri, Σrr, Σrj, Σjj)
#     _, _, n = size(Σrr)
#     sum = 0.0
#     for i = 1:n
#         Σrri = @view Σrr[:, :, i]
#         Σrji = @view Σrj[:, i]
#         Σjji = Σjj[i]
#         sum += rotated_corr_single(Ri, Σrri, Σrji, Σjji)
#     end
#     return -sum / n
# end

function optimize_correlation_from_variance(Σ::Array{Float64, 3},
        rot_inds::AbstractVector{Int},
        target_inds::AbstractVector{Int})

    if !isempty(intersect(rot_inds, target_inds))
        error("not yet implemented") # is probably fine but need to check
    end

    k, _, n = size(Σ)

    if length(target_inds) > k
        error("cannot target more indices than there are dimensions")
    end

    kr = length(rot_inds)

    R = zeros(kr, 0)
    for ind in target_inds
        if size(R, 2) == kr - 1
            r = nullspace(R')
            @show mean(r' * Σ[rot_inds, ind, i] for i = 1:n)
            if mean(r' * Σ[rot_inds, ind, i] for i = 1:n)[1] < 0
                r = -r
            end
        else
            r = optimize_correlation_from_variance(Σ, rot_inds, ind, R)
        end
        R = [R r]
    end

    if size(R, 2) < kr
        R = [R nullspace(R')]
    end

    @assert R' * R ≈ I

    R
end

function optimize_correlation_from_variance(Σ::Array{Float64, 3},
        rot_inds::AbstractVector{Int},
        target_ind::Int,
        R::Matrix{Float64})

    k, _, n = size(Σ)

    @show rot_inds
    @show target_ind

    model = Model(Ipopt.Optimizer)
    set_silent(model)
    kr = length(rot_inds)

    @variable(model, r[1:kr])
    @variable(model, -1 <= c[1:n] <= 1)
    @variable(model, Σii[1:n] >= 0 + 1e-5)
    @variable(model, Σij[1:n])
    @constraint(model, Σii == [r' * Σ[rot_inds, rot_inds, i] * r for i in 1:n])
    @constraint(model, Σij == [r' * Σ[rot_inds, target_ind, i] for i in 1:n])
    @constraint(model, c == [Σij[i] / sqrt(Σii[i] * Σ[target_ind, target_ind, i]) for i in 1:n])
    @constraint(model, sum(r.^2) == 1)

    m = size(R, 2)
    for i = 1:m
        @constraint(model, r' * R[:, i] == 0)
    end

    @objective(model, Max, sum(c))
    optimize!(model)
    assert_is_solved_and_feasible(model)
    return(value.(r))
end


# function optimize_correlation_from_variance(Σrr::Array{Float64, 3},
#         Σrj::Array{Float64, 2}, Σjj::Vector{Float64})
#     k, n = size(Σrj)
#     @show size(Σrj)

#     optf = OnceDifferentiable(r -> rotated_corr(r, Σrr, Σrj, Σjj), ones(k), autodiff=:forward)
#     manif = Optim.Sphere()
#     u0 = randn(k)
#     # @show u0
#     # error()
#     u0 = u0 / norm(u0)
#     opt = Optim.optimize(optf, u0, Optim.ConjugateGradient(manifold=manif))
#     display(opt)
#     r = Optim.minimizer(opt)
#     R = zeros(k, k)
#     R[1, :] .= r
#     if k > 1
#         R[2:k, :] .= nullspace(R[1, :]')'
#     end
#     return R
# end



function optimize_one_trait(C::Array{Float64, 2}; kwargs...)
    @warn "Ignoring provided optimization function"
    c_mean = mean(C, dims=2)
    k = length(c_mean)

    R = zeros(k, k)
    R[1, :] .= c_mean / norm(c_mean)
    if k > 1
        R[2:k, :] .= nullspace(R[1, :]')'
    end

    return R
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



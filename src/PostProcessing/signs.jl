using Base: Float64
struct Signs <: Transform
    signs::Matrix{Int}
end

function apply_transform(data::Array{Float64, 3}, signs::Signs)
    p, k, n = size(data)
    new_data = copy(data)
    for i = 1:n
        for j = 1:k
            new_data[:, j, i] .*= signs.signs[j, i]
        end
    end

    return new_data
end

function apply_inverse(data::Array{Float64, 3}, signs::Signs)
    apply_transform(data, signs)
end


function naive_constrain_signs(data::Array{Float64, 3})
    p, k, n = size(data)
    abs_data = abs.(data) # todo: probably don't need to re-allocate
    μ = mean(abs_data, dims=3)
    σ = std(abs_data, dims=3, mean = μ)
    Z = μ ./ σ

    pivots = zeros(Int, k)
    for i = 1:k
        z = @view Z[:, i]
        pivots[i] = findmax(z)[2]
    end

    signs = zeros(Int, k, n)

    for i = 1:n
        for j = 1:k
            signs[j, i] = sign(data[pivots[j], j, i])
        end
    end

    return Signs(signs)
end

function naive_constrain_signs!(data::Array{Float64, 3})
    signs = naive_constrain_signs(data)
    data .= apply_transform(data, signs)
    return signs
end
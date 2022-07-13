struct Rotations <: Transform
    rotations::Array{Float64, 3}
end

function Rotations(::Int, k::Int, n::Int)
    return Rotations(k, n)
end

function Rotations(k::Int, n::Int)
    R = zeros(k, k, n)
    for i = 1:n
        for j = 1:k
            R[j, j, i] = 1.0
        end
    end
    return Rotations(R)
end

function view_rotation(r::Rotations, n::Int)
    return @view r.rotations[:, :, n]
end

function update_rotation!(r::Rotations, R::AbstractMatrix{Float64}, n::Int)
    r.rotations[:, :, n] .= r.rotations[:, :, n] * R
end

function rotate_sample(data::AbstractArray{Float64, 3}, r::Rotations, n::Int)
    R = view_rotation(r, n)
    Y = @view data[:, :, n]
    return Y * R
end

function apply_transform(data::Array{Float64, 3}, rotations::Rotations)
    new_data = zeros(size(data)...)
    for i = 1:size(data, 3)
        new_data[:, :, i] = rotate_sample(data, rotations, i)
    end
    return new_data
end

function apply_inverse(data::Array{Float64, 3}, rotations::Rotations;
                           transpose::Bool = false)

    if transpose
        return apply_transform(data, rotations) # R' = inv(R) for rotation matrices
    else
        p, k, n = size(data)
        new_data = zeros(p, k, n)
        for i = 1:n
            Y = @view data[:, :, i]
            new_data[:, :, i] .= Y * view_rotation(rotations, i)'
        end

        return new_data
    end
end

abstract type AbstractRotation <: AbstractTransformer end
# interface:
#   update_rotation!
#   update_statistics!

function compute_transform!(t::Type{<:AbstractRotation}, X::Array{Float64, 3},
                            reference_ind::Int)
    return rotate!(X, t(X[:, :, reference_ind]))
end

function compute_transform!(t::Type{<:AbstractRotation}, X::Array{Float64, 3})
    return rotate!(X, t(mean(X, dims=3)[:, :, 1]))
end


function update_statistics!(::Rotations, ::Array{Float64, 3},
                            ::AbstractRotation)
    return 0.0 # default is non-iterative (returns 0)
end

include("procrustes.jl")


function rotate(data::Array{Float64, 3}, rotator::AbstractRotation;
                rotations = Rotations(size(data)...),
                max_iterations::Int = 100,
                tolerance::Float64 = 1e-8)

    @assert tolerance > 0

    iteration = 0
    diff = Inf
    while diff > tolerance && iteration < max_iterations
        update_rotation!(rotations, data, rotator)
        diff = update_statistics!(rotations, data, rotator) # returns 0 if non-iterative
        iteration += 1
    end

    return rotations
end

import LinearAlgebra.rotate!

function rotate!(data::Array{Float64, 3}, rotator::AbstractRotation; kwargs...)
    rotation = rotate(data, rotator; kwargs...)
    data .= apply_transform(data, rotation) #TODO: make more memory efficient
    return rotation
end
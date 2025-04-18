module PostProcessing

export cluster,
       cluster!,
       rearrange_data,
       MultivariateNormal,
       IndependentNormal,
       Geodesic,
       rotate,
       rotate!,
       Rotations,
       ProcrustesRotation,
       Permutation,
       SVDRotation,
       SignRotation,
       apply_transform,
       apply_inverse,
       rotate_sample,
       Signs,
       naive_constrain_signs,
       naive_constrain_signs!,
       compose,
       check_valid_rotation,
       check_transform_results,
       expand_rotation,
       RotationPlan,
       do_rotations!,
       optimize,
       induce_block_independence,
       induce_block_independence!

using LinearAlgebra, UnPack, Combinatorics, Statistics
# using Optim
# using ForwardDiff
using JuMP
using Ipopt


abstract type Transform end
# interface:
#   apply_transform
#   apply_inverse
#   compose

abstract type AbstractTransformer end
# interface:
#   compute_transform!


struct RotationPlan
    transformers::Vector{DataType}

    function RotationPlan(transformers::Vector{DataType})
        for t in transformers
            if !(t <: AbstractTransformer)
                throw(ArgumentError("$t is not an valid transformer"))
            end
        end
        return new(transformers)
    end
end

function RotationPlan(transformers::Type{<:AbstractTransformer}...)
    return RotationPlan([t for t in transformers])
end


function do_rotations!(r::RotationPlan, X::Array{Float64, 3}, reference_ind::Int)
    n = length(r.transformers)
    rs = Vector{Transform}(undef, n)
    for i = 1:n
        rs[i] = reference_ind == 0 ? compute_transform!(r.transformers[i], X) :
                compute_transform!(r.transformers[i], X, reference_ind)
    end
    return compose(rs...)
end


include("clustering.jl")
include("summaries.jl")
include("rotation.jl")
include("utils.jl")
include("signs.jl")
include("conversion.jl")
include("testing.jl")
include("optimization.jl")
include("joint_covariance.jl")


end

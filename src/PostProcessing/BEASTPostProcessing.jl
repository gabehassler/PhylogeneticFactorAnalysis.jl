module BEASTPostProcessing
#want to keep BEAST-specific code seperate from more general code

export post_process

using Statistics: cov2cor!
using DataFrames
using CSV
using Statistics
using LinearAlgebra
using UnPack
using BeastUtils.Logs
using BeastUtils.MatrixUtils
using PhylogeneticFactorAnalysis.PostProcessing

const L_HEADER = "L"
const FAC_HEADER = "factors."
const PREC_HEADER = "factorPrecision"

struct JointParameters
    tree_dims::Vector{Int}
    data_dims::Vector{Int}
    trait_names::Vector{String}
    joint_name::String
    n_taxa::Int
end

function subset_startswith(df::DataFrame, s::AbstractString)
    inds = findall(x -> startswith(x, s), names(df))
    return df[!, inds]
end

function count_startswith(S::Array{String}, s::String)
    return count(x -> startswith(x, s), S)
end

import LinearAlgebra.adjoint
function adjoint(X::AbstractMatrix{<:AbstractString})
    k, p = size(X)
    Y = [X[i, j] for j = 1:p, i = 1:k]
end

function absmax(x::Array{<:Real})
    return maximum(abs.(x))
end


function collect_matrix(df::DataFrame, header::String, n_rows::Int, n_cols::Int;
                        transpose::Bool = false)

    df_sub = subset_startswith(df, header)
    if size(df_sub, 2) == 0
        @warn "could not find rows starting with $header in DataFrame. skipping"
        return zeros(0, 0, 1), Matrix{String}(undef, 0, 0)
    end
    labels = names(df_sub)

    labels = reshape(labels, n_rows, n_cols)
    labels = transpose ? labels' : labels

    @assert size(df_sub, 2) == n_rows * n_cols

    n_states = size(df_sub, 1)

    final_rows, final_cols = transpose ? (n_cols, n_rows) : (n_rows, n_cols)

    data = zeros(final_rows, final_cols, n_states)
    M = Matrix(df_sub)
    for i = 1:n_states
        X = reshape(M[i, :], n_rows, n_cols)
        data[:, :, i] .= transpose ? X' : X
    end

    return data, labels
end


function make_2d(data::AbstractMatrix{Float64})
    return data
end

function make_2d(data::AbstractArray{Float64, 3})
    p, k, n = size(data)

    return reshape(data, p * k, n)
end

function update_df!(df::DataFrame,
        data_sets::Vector{<:Pair{<:AbstractArray{<:AbstractString},
                                 <:AbstractArray{Float64}}
                         }
                   )
    for set in data_sets
        labels = set[1]
        data = set[2]

        reshaped_data = make_2d(data)

        for i = 1:length(labels)
            df[!, labels[i]] = reshaped_data[i, :]
        end

    end
    return nothing
end

function transpose_samples(X::Array{Float64, 3})
    p, k, n = size(X)
    Y = zeros(k, p, n)
    for i = 1:n
        Y[:, :, i] .= (@view X[:, :, i])'
    end
    return Y
end





function rotate_other_parameters!(F::Array{Float64, 3}, C::Array{Float64, 3},
                                  V::Array{Float64, 3}, rotations::Rotations,
                                  model_dims::AbstractVector{Int}, joint_dim::Int)

    expanded_rotation = expand_rotation(rotations, model_dims, joint_dim)

    if size(F, 1) > 0
        F .= apply_inverse(F, expanded_rotation, transpose = true)
    end

    other_dims = setdiff(1:joint_dim, model_dims)
    V_sub = V[other_dims, model_dims, :]
    V_sub = apply_inverse(V_sub, rotations, transpose = true)

    for i = 1:size(V, 3)
        Vr = @view V_sub[:, :, i]
        V[other_dims, model_dims, i] .= Vr
        V[model_dims, other_dims, i] .= Vr'
        if size(C, 1) > 0
            C[:, :, i] .= cov2corr(V[:, :, i])
        end
    end
end


function rotate_submodel!(df::DataFrame, parameters::JointParameters,
        plan::RotationPlan,
        model::Int;
        map_ind::Int = 1,
        double_check::Bool = false,
        optimization::Union{Nothing, Function} = nothing)
    dim_factor = parameters.tree_dims[model]
    dim_trait = parameters.data_dims[model]
    dim_joint = sum(parameters.tree_dims)
    trait_name = parameters.trait_names[model]
    joint_name = parameters.joint_name
    n_taxa = parameters.n_taxa

    L, L_labels = collect_matrix(df, "$trait_name.$L_HEADER", dim_trait, dim_factor)
    F, F_labels = collect_matrix(df, joint_name, dim_joint, n_taxa,
                            transpose = true)
    C, C_labels = collect_matrix(df, "correlation.", dim_joint, dim_joint)
    V, V_labels = collect_matrix(df, "mbd.variance", dim_joint, dim_joint)
    @assert size(C)[1:2] == size(V)[1:2] == (dim_joint, dim_joint)


    prec_header = "$trait_name.factorPrecision"
    prec_df = subset_startswith(df, prec_header)
    prec = Matrix(prec_df)'
    prec_labels = names(prec_df)
    @assert length(prec_labels) == dim_trait


    prop_header = "$trait_name.proportion"
    prop_df = subset_startswith(df, prop_header)
    @assert size(prop_df, 2) == 2 + dim_factor * 2


    final_rotation  = do_rotations!(plan, L, map_ind)

    offset = sum(parameters.tree_dims[1:(model - 1)])
    model_inds = (offset + 1):(offset + dim_factor)
    other_inds = setdiff(1:dim_joint, model_inds)

    rotate_other_parameters!(F, C, V, final_rotation, model_inds, dim_joint)

    F_model = @view F[:, model_inds, :]
    update_proportions!(prop_df, F_model, L, prec, prop_header)

    n = size(C, 3)

    if !isnothing(optimization)
        if count(parameters.tree_dims .!= parameters.data_dims) > 1
            @warn "optimization with multiple factor models has not been tested"
        end
        C_sub = C[model_inds, other_inds, :]
        R = PostProcessing.optimize_one_trait(C_sub, optimization)
        R = R'
        optimized_rotation = zeros(dim_factor, dim_factor, n) # TODO: don't need to repeat same matrix over and over
        for i = 1:n
            optimized_rotation[:, :, i] .= R
        end

        optimized_rotation = Rotations(optimized_rotation)

        L = apply_transform(L, optimized_rotation)
        rotate_other_parameters!(F, C, V, optimized_rotation, model_inds, dim_joint)
    end


    prop_labels = names(prop_df)
    props = Matrix(prop_df)

    update_df!(df,
             [L_labels => L, F_labels => F, C_labels => C,
             V_labels => V, prop_labels => props']
            )
end




import Base:length
function length(j::JointParameters)
    return length(j.tree_dims)
end

function is_factor(j::JointParameters, ind::Int)
    return j.tree_dims[ind] != j.data_dims[ind]
end

function rotate_multi_sem(log_path::String, rotate_path::String,
        plan::RotationPlan, parameters::JointParameters;
        burnin::Float64 = 0.0,
        minimum_map::Float64 = max(0.1 - burnin, 0.0),
        kw_args...)

    @unpack tree_dims, data_dims = parameters
    q = sum(tree_dims)
    m = length(parameters)

    df = import_log(log_path, burnin=burnin)
    n = size(df, 1)

    map_offset = Int(round(n * minimum_map))
    map_ind = findmax(df.joint[(map_offset + 1):end])[2] + map_offset


    for i = 1:m
        if is_factor(parameters, i)
            rotate_submodel!(df, parameters, plan, i; map_ind = map_ind,
                    kw_args...)
        end
    end
    CSV.write(rotate_path, df, delim = '\t')
end

function post_process(log_path::String,
                      rotated_path::String,
                      traits::Vector{Pair{String, Tuple{Int, Int}}},
                      n_taxa::Int;
                      optimize::Bool = false,
                      rotation_plan = RotationPlan(SVDRotation,
                                                   ProcrustesRotation)
                      )

    n = length(traits)
    tree_dims = [x[2][1] for x in traits]
    data_dims = [x[2][2] for x in traits]
    trait_names = [x[1] for x in traits]
    joint_name = join(trait_names, '.') * ".joint"

    params = JointParameters(tree_dims, data_dims, trait_names, joint_name, n_taxa)

    optimization = optimize ? PostProcessing.maximum_correlation : nothing

    rotate_multi_sem(log_path, rotated_path,
           rotation_plan,
           params,
           double_check = true,
           optimization = optimization)
    return nothing
end


function factor_proportion_statistics(F::AbstractMatrix{Float64},
                                      L::AbstractMatrix{Float64},
                                      λ::AbstractVector{Float64})
    FtF = F' * F
    LLt = L' * L

    f = mean(F, dims=1)
    n = size(F, 1)
    FtF -= n * f' * f

    H = LLt .* FtF #hadamard product
    h = diag(H)
    fac_var = sum(H)
    ind_var = sum(h)

    res_var = (n - 1) * sum(1/x for x in λ)

    total_var = sum(H) + res_var

    fac_prop = fac_var / total_var
    ind_prop = ind_var / fac_var

    fac_props = h ./ fac_var
    total_props = fac_props .* fac_prop

    return (factor_proportion = fac_prop,
            absolute_proportions = total_props,
            relative_proportions = fac_props,
            marginal_proportion = ind_prop
            )

end

function update_proportions!(df::DataFrame,
                             F::AbstractArray{Float64, 3},
                             L::AbstractArray{Float64, 3},
                             λ::AbstractArray{Float64, 2},
                             header::AbstractString)
    n = size(F, 3)
    k = size(F, 2)
    for i = 1:n
        Fi = @view F[:, :, i]
        Li = @view L[:, :, i]
        λi = @view λ[:, i]
        fp, aps, rps, mp = factor_proportion_statistics(Fi, Li, λi)
        df[i, header * ".factorProportion"] = fp
        df[i, header * ".relativeMarginalProportion"] = mp
        for j = 1:k
            df[i, header * ".absoluteProportion.$j"] =  aps[j]
            df[i, header * ".relativeProportion.$j"] =  rps[j]
        end
    end
end


end
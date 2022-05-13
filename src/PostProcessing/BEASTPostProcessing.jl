module BEASTPostProcessing
#want to keep BEAST-specific code seperate from more general code

export post_process

using Statistics: cov2cor!
# using PhylogeneticFactorAnalysis.PostProcessing: ClusterStatistics, optimize_2d, minimal_variance
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

function collect_loadings(log_path::String; kwargs...)
    df = CSV.read(log_path, DataFrame, header=4, delim='\t')
    return collect_loadings(df; kwargs...)
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

function collect_dimensions(df::DataFrame;
                            L_header::String = L_HEADER,
                            fac_header::String = FAC_HEADER,
                            prec_header::String = PREC_HEADER,
                            extra_traits::Int = 0)

    nms = names(df)
    n_traits = count_startswith(nms, prec_header)

    n_factors, r = divrem(count_startswith(nms, L_header), n_traits)
    @assert r == 0

    n_taxa, r = divrem(count_startswith(nms, fac_header), n_factors + extra_traits)
    @assert r == 0

    return (n_traits = n_traits, n_factors = n_factors, n_taxa = n_taxa)
end


function collect_loadings(df::DataFrame;
                          L_header = "L",
                          fac_header = "factors.",
                          prec_header = "factorPrecision",
                          normalization = "none")

    @unpack n_traits, n_factors, n_taxa = collect_dimensions(df)
    k = n_factors

    n_states = size(df, 1)

    L_data = subset_startswith(df, L_header)
    fac_data = subset_startswith(df, fac_header)

    L = zeros(n_traits, k, n_states)
    L_mat = Matrix(L_data)
    F_mat = Matrix(fac_data)
    for state = 1:n_states
        L_sub = reshape(L_mat[state, :], n_traits, k)

        if normalization == "none"
            L[:, :, state] .= L_sub
        elseif normalization == "factor_sd"

            F_sub = reshape(F_mat[state, :], k, n_taxa)

            F_std = std(F_sub, dims=2, corrected=false)
            L[:, :, state] .= L_sub * Diagonal(F_std)
        elseif normalization == "loadings_norm"
            L_norms = [norm(L_sub[:, i]) for i = 1:k]
            L[:, :, state] .= L_sub * inv(Diagonal(L_norms))
        else
            error("Unknown normalization")
        end
    end

    return L
end

function make_labels(header::String, x::Vector{Float64})
    return header
end

function make_labels(header::String, x::Matrix{Float64})
    return ["header$i" for i = 1:size(x, 1)]
end

function make_labels(header::String, x::Array{Float64, 3})
    p, k, _ = size(x)
    return vec(["$header$j$i" for i = 1:p, j = 1:k])
end

function make_labels(labels::Array{String}, x::Array{Float64, 3})
    s = size(@view x[:, :, 1])
    if size(labels) == s
        return labels
    elseif size(labels') == s
        return labels'
    else
        @show s
        @show size(labels)
        error("incompatible dimensions")
    end
end

function make_2d(data::Matrix{Float64})
    return data
end

function make_2d(data::Array{Float64, 3})
    p, k, n = size(data)

    return reshape(data, p * k, n)
end

T = Union{AbstractString, AbstractArray{<:AbstractString}}

function save_log(log_path::String, states::Vector{Int},
                  data_sets::Vector{<:Pair{<:T, <:Array{Float64}}})


    df = DataFrame(state = states)
    save_log(log_path, df, data_sets)
end


function save_log(log_path::String, df::DataFrame,
                  data_sets::Vector{<:Pair{<:T, <:Array{Float64}}})

    update_df!(df, data_sets)
    CSV.write(log_path, df, delim='\t')
end

function update_df!(df::DataFrame,
        data_sets::Vector{<:Pair{<:T, <:Array{Float64}}})
    for set in data_sets
        header = set[1]
        data = set[2]

        labels = make_labels(header, data)

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

function loadings_and_factors(df::DataFrame; extra_traits::Int = 0)
    @unpack n_traits, n_factors, n_taxa = collect_dimensions(df,
                                                extra_traits = extra_traits)
    n_factors_plus = n_factors + extra_traits
    L, L_cols = collect_matrix(df, L_HEADER, n_traits, n_factors)
    F, F_cols = collect_matrix(df, FAC_HEADER, n_factors_plus,
                               n_taxa, transpose=true)

    C, C_labels = collect_matrix(df, "traits.cor", n_factors_plus, n_factors_plus)
    V, V_labels = collect_matrix(df, "traits.variance.matrix", n_factors_plus, n_factors_plus)
    return (L = L, F = F, C = C, V = V,
            L_labels = L_cols, F_labels = F_cols,
            C_labels = C_labels, V_labels = V_labels)
end

function loadings_and_factors(path::String; burnin::Float64 = 0.0, kwargs...)
    df = import_log(path, burnin=burnin)
    return loadings_and_factors(df; kwargs...)
end

import PhylogeneticFactorAnalysis.PostProcessing.check_valid_rotation
function check_valid_rotation(p1::String, p2::String;
                              b1::Float64 = 0.0, b2::Float64 = 0.0,
                              tol::Float64 = 1e-10,
                              extra_traits::Int = 0)
    lf1 = loadings_and_factors(p1, burnin = b1, extra_traits = extra_traits)
    lf2 = loadings_and_factors(p2, burnin = b2, extra_traits = extra_traits)
    L1, F1, C1, V1 = lf1.L, lf1.F, lf1.C, lf1.V
    L2, F2, C2, V2 = lf2.L, lf2.F, lf2.C, lf2.V

    p, k, n = size(L1)
    _, q, m = size(F1)
    @assert m == n
    @assert q == k + extra_traits

    for i = 1:n
        Y1 = F1[:, 1:k, i] * L1[:, :, i]'
        Y2 = F2[:, 1:k, i] * L2[:, :, i]'
        diff = absmax(Y1 - Y2)
        if diff > tol
            error("The maximum difference between the two Y_hat values is " *
                  "$diff. Maximum allowable difference is $tol.")
        end
    end

    diff = absmax(F1[:, (k + 1):q, :] - F2[:, (k + 1):q, :])
    if diff > tol
        error("The maximum difference between the unrotated values is " *
              "$diff. Maximum allowable difference is $tol.")
    end

    for i = 1:n
        P1 = inv(V1[:, :, i])
        P2 = inv(V2[:, :, i])
        F1i = F1[:, :, i]
        S1 = F1i * P1 * F1i'

        F2i = F2[:, :, i]
        S2 = F2i * P2 * F2i'

        diff = absmax(S1 - S2)
        if diff > tol
            error("The maximum difference between the two precision-adjusted " *
                  "inner products in $diff. Maximum allowable difference is $tol.")
        end
    end

    for i = 1:n
        C1_sub = cov2corr(V1[:, :, i])
        C2_sub = cov2corr(V2[:, :, i])
        diff1 = absmax(C1_sub - C1[:, :, i])
        diff2 = absmax(C2_sub - C2[:, :, i])
        diff = max(diff1, diff2)
        if diff > tol
            error("The maximum difference between the correlations is $diff. " *
                  "Maximum allowable difference is $tol.")
        end
    end
end



# function rotate_other_parameters!(F::Array{Float64, 3}, C::Array{Float64, 3},
#                                   V::Array{Float64, 3}, rotations::Rotations,
#                                   n_factors::Int, extra_traits::Int)

#     expanded_rotation = expand_rotation(rotations, extra_traits)
#     if (size(F, 1)) > 0
#         F .= apply_inverse(F, expanded_rotation, transpose = true)
#     end

#     n_factors_plus = n_factors + extra_traits

#     r1 = (n_factors + 1):n_factors_plus
#     r2 =  1:n_factors
#     V_sub = V[r1, r2, :]
#     V_sub = apply_inverse(V_sub, rotations, transpose = true)

#     for i = 1:size(V, 3)
#         Vr = @view V_sub[:, :, i]
#         V[r1, r2, i] .= Vr
#         V[r2, r1, i] .= Vr'

#         C[:, :, i] .= cov2corr(V[:, :, i])
#     end
# end

function rotate_other_parameters!(F::Array{Float64, 3}, C::Array{Float64, 3},
                                  V::Array{Float64, 3}, rotations::Rotations,
                                  n_factors::Int, extra_traits::Int)

    expanded_rotation = expand_rotation(rotations, extra_traits)

    F .= apply_inverse(F, expanded_rotation, transpose = true)

    n_factors_plus = n_factors + extra_traits

    r1 = (n_factors + 1):n_factors_plus
    r2 =  1:n_factors
    V_sub = V[r1, r2, :]
    V_sub = apply_inverse(V_sub, rotations, transpose = true)

    for i = 1:size(V, 3)
        Vr = @view V_sub[:, :, i]
        V[r1, r2, i] .= Vr
        V[r2, r1, i] .= Vr'

        C[:, :, i] .= cov2corr(V[:, :, i])
    end
end



function rotate_sem(in_path::String, out_path::String, plan::RotationPlan;
                    extra_traits::Int = 0,
                    burnin::Float64 = 0.0,
                    double_check::Bool = false,
                    optimization::Union{Nothing, Function} = nothing)
    df = import_log(log_path, burnin=burnin)
    n = size(df, 1)
    @unpack n_traits, n_factors, n_taxa = collect_dimensions(df,
                                                extra_traits = extra_traits)
    map_ind = findmax(df.joint)[2]

    n_factors_plus = n_factors + extra_traits
    L, L_labels = collect_matrix(df, L_HEADER, n_traits, n_factors)
    F, F_labels = collect_matrix(df, FAC_HEADER, n_factors_plus, n_taxa,
                               transpose = true)
    C, C_labels = collect_matrix(df, "correlation.", n_factors_plus, n_factors_plus)
    V, V_labels = collect_matrix(df, "variance", n_factors_plus, n_factors_plus)
    @assert size(C) == size(V) == (n_factors_plus, n_factors_plus)


    final_rotation  = do_rotations!(plan, L, map_ind)



    if double_check
        L_original = loadings_and_factors(df, extra_traits = extra_traits).L
        check_valid_rotation(final_rotation)
        check_transform_results(L_original, L, final_rotation)
    end

    rotate_other_parameters!(F, C, V, final_rotation, n_factors, extra_traits)

    if !isnothing(optimization)
        C_sub = C[1:n_factors, (n_factors + 1):end, :]
        R = optimize(C_sub, optimization)
        R = R'
        optimized_rotation = zeros(n_factors, n_factors, n) # TODO: don't need to repeat same matrix over and over
        for i = 1:n
            optimized_rotation[:, :, i] .= R
        end

        optimized_rotation = Rotations(optimized_rotation)

        L = apply_transform(L, optimized_rotation)
        rotate_other_parameters!(F, C, V, optimized_rotation, n_factors, extra_traits)
    end

    save_log(out_path,
             df,
             [L_labels => L, F_labels => transpose_samples(F), C_labels => C, V_labels => V]
            )

    if double_check
        check_valid_rotation(in_path, out_path, b1 = burnin, tol=1e-3,
                             extra_traits = extra_traits)
    end
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


    final_rotation  = do_rotations!(plan, L, map_ind)

    offset = sum(parameters.tree_dims[1:(model - 1)])
    model_inds = (offset + 1):(offset + dim_factor)
    other_inds = setdiff(1:dim_joint, model_inds)

    rotate_other_parameters!(F, C, V, final_rotation, model_inds, dim_joint)

    n = size(C, 3)

    if !isnothing(optimization)
        if count(parameters.tree_dims .!= parameters.data_dims) > 1
            @warn "optimization with multiple factor models has not been tested"
        end
        C_sub = C[model_inds, other_inds, :]
        R = optimize(C_sub, optimization)
        R = R'
        optimized_rotation = zeros(dim_factor, dim_factor, n) # TODO: don't need to repeat same matrix over and over
        for i = 1:n
            optimized_rotation[:, :, i] .= R
        end

        optimized_rotation = Rotations(optimized_rotation)

        L = apply_transform(L, optimized_rotation)
        rotate_other_parameters!(F, C, V, optimized_rotation, model_inds, dim_joint)
    end

    update_df!(df,
             [L_labels => L, F_labels => transpose_samples(F), C_labels => C, V_labels => V]
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
        kw_args...)

    @unpack tree_dims, data_dims = parameters
    q = sum(tree_dims)
    m = length(parameters)

    df = import_log(log_path, burnin=burnin)
    n = size(df, 1)

    map_ind = findmax(df.joint)[2]


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
                      optimize::Bool = false
                      )

    n = length(traits)
    tree_dims = [x[2][1] for x in traits]
    data_dims = [x[2][2] for x in traits]
    trait_names = [x[1] for x in traits]
    joint_name = join(trait_names, '.') * ".joint"

    params = JointParameters(tree_dims, data_dims, trait_names, joint_name, n_taxa)

    optimization = optimize ? PostProcessing.maximum_correlation : nothing

    rotate_multi_sem(log_path, rotated_path,
           RotationPlan(SVDRotation, ProcrustesRotation),
           params,
        #    extra_traits = extra_traits,
           burnin = 0.1,
           double_check = true,
           optimization = optimization)
    return nothing
end

end
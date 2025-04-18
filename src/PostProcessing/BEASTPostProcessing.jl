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
    # Construct the regular expression pattern to match `s` followed by numbers
    pattern = r"^" * s * r".*\d+"
    # Find the indices of column names that match the pattern
    inds = findall(x -> occursin(pattern, x), names(df))
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
    M = Matrix{Float64}(df_sub)
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


function rotate_other_parameters!(matrices::Vector{Array{Float64, 3}},
                                  rotations::Rotations,
                                  inverses::Vector{Bool},
                                  sides::Vector{Symbol},
                                  expand::Vector{Bool},
                                  model_dims::AbstractVector{Int},
                                  joint_dim::Int)

    expanded_rotations = expand_rotation(rotations, model_dims, joint_dim)
    n = length(matrices)
    for i = 1:n
        if minimum(size(matrices[i])) > 0
            rotation = expand[i] ? expanded_rotations : rotations
            matrices[i] .= apply_transform(matrices[i], rotation, side = sides[i], inverse = inverses[i])
        end
    end
end





# function rotate_other_parameters!(F::Array{Float64, 3}, C::Array{Float64, 3},
#                                   V::Array{Float64, 3}, rotations::Rotations,
#                                   model_dims::AbstractVector{Int}, joint_dim::Int)

#     expanded_rotation = expand_rotation(rotations, model_dims, joint_dim)

#     if size(F, 1) > 0
#         F .= apply_inverse(F, expanded_rotation, transpose = true)
#     end

#     other_dims = setdiff(1:joint_dim, model_dims)

#     if size(V, 1) > 0
#         V_sub = V[other_dims, model_dims, :]
#         V_sub = apply_inverse(V_sub, rotations, transpose = true)

#         for i = 1:size(V, 3)
#             Vr = @view V_sub[:, :, i]
#             V[other_dims, model_dims, i] .= Vr
#             V[model_dims, other_dims, i] .= Vr'
#             if size(C, 1) > 0
#                 C[:, :, i] .= cov2corr(V[:, :, i]) #TOCO: just update correlation
#             end
#         end
#     else
#         if size(C, 1) > 0
#             error("not implemented when correlation is present but no variance")
#         end
#     end

# end


function rotate_submodel!(df::DataFrame, parameters::JointParameters,
        plan::RotationPlan,
        model::Int;
        map_ind::Int = 0,
        double_check::Bool = false,
        optimization_inds::Vector{Vector{Int}},
        variance_header::AbstractString = "mbd.variance",
        correlation_header::AbstractString = "correlation.",
        L_header::AbstractString = "$(parameters.trait_names[model]).$L_HEADER",
        F_header::AbstractString = parameters.joint_name,
        prec_header = "$(parameters.trait_names[model]).factorPrecision",
        prop_header = "$(parameters.trait_names[model]).proportion",
        resvar_header = "$(parameters.trait_names[model]).variance"
        )
    dim_factor = parameters.tree_dims[model]
    dim_trait = parameters.data_dims[model]
    dim_joint = sum(parameters.tree_dims)
    trait_name = parameters.trait_names[model]
    joint_name = parameters.joint_name
    n_taxa = parameters.n_taxa

    L, L_labels = collect_matrix(df, L_header, dim_trait, dim_factor)
    F, F_labels = collect_matrix(df, F_header, dim_joint, n_taxa,
                            transpose = true)
    C, C_labels = collect_matrix(df, correlation_header, dim_joint, dim_joint)
    V, V_labels = collect_matrix(df, variance_header, dim_joint, dim_joint)
    C_dims = size(C)[1:2]
    if !(C_dims == (dim_joint,dim_joint) || C_dims == (0, 0))
        error("Correlation matrix with header $correlation_header has size " *
              "$C_dims. Should be ($dim_joint, $dim_joint)")
    end
    V_dims = size(V)[1:2]
    if !(V_dims == (dim_joint,dim_joint) || V_dims == (0, 0))
        error("Variance matrix with header $variance_header has size " *
              "$V_dims. Should be ($dim_joint, $dim_joint)")
    end
    if C_dims == (0, 0) || V_dims == (0, 0)
        @warn "Could not find either correlation or variance matrix for " *
              "post-processing. If this is a standard latent factor model " *
              "then this is expected. If the variance is not fixed to the " *
              "identity, then the output variance/correlation will not be " *
              "post-processed."
    end

    Vr, Vr_labels = collect_matrix(df, resvar_header, dim_factor, dim_factor)


    prec_df = subset_startswith(df, prec_header)
    prec = Matrix(prec_df)'
    prec_labels = names(prec_df)
    if length(prec_labels) != dim_trait
        error("Unknown trait dimension. $(length(prec_labels)) based on parsing the log file, but $dim_trait supplied")
    end




    final_rotation  = do_rotations!(plan, L, map_ind)

    offset = sum(parameters.tree_dims[1:(model - 1)])
    model_inds = (offset + 1):(offset + dim_factor)
    other_inds = setdiff(1:dim_joint, model_inds)

    rotate_other_parameters!([F, C, V, Vr], final_rotation,
                             [false, false, false, false],
                             [:right, :both, :both, :both],
                             [true, true, true, false],
                             model_inds, dim_joint)


    n = size(C, 3)

    if !isempty(optimization_inds[model])
        if !isempty(intersect(abs.(optimization_inds[model]),  model_inds))
            error("cannot optimize correlation within a model")
        end

        # @assert length(optimization_inds[model]) == 1 # TODO: generalize

        if count(parameters.tree_dims .!= parameters.data_dims) > 1
            @warn "optimization with multiple factor models has not been tested"
        end
        # C_sub = C[model_inds, abs.(optimization_inds[model]), :]
        # c = zeros(dim_factor, n)
        # signs = sign.(optimization_inds[model])
        # @show size(C_sub)
        # @show signs
        # for i = 1:n
        #     c[:, i] .= C_sub[:, :, i] * signs
        # end

        # Vrr = V[model_inds, model_inds, :]
        # Vrj = V[model_inds, optimization_inds[model][1], :]
        # Vjj = V[optimization_inds[model][1], optimization_inds[model][1], :]
        # R = PostProcessing.optimize_correlation_from_variance(Vrr, Vrj, Vjj)
        R = PostProcessing.optimize_correlation_from_variance(V, model_inds, optimization_inds[model])

        # R = PostProcessing.optimize_one_trait(c); @warn "NEED TO FIX THIS TO DEAL WITH NON-IDENTITY SUBMATRIX"
        # R = R'
        optimized_rotation = zeros(dim_factor, dim_factor, n) # TODO: don't need to repeat same matrix over and over
        for i = 1:n
            optimized_rotation[:, :, i] .= R
        end

        optimized_rotation = Rotations(optimized_rotation)

        L = apply_transform(L, optimized_rotation)
        rotate_other_parameters!([F, C, V, Vr], optimized_rotation,
                                 [false, false, false, false],
                                 [:right, :both, :both, :both],
                                 [true, true, true, false],
                                 model_inds, dim_joint)
    end
    adjust_correlation!(C, V)


    F_model = @view F[:, model_inds, :]


    labs = [L_labels => L, F_labels => F, C_labels => C,
            V_labels => V, Vr_labels => Vr]

    prop_df = subset_startswith(df, prop_header)

    if ncol(prop_df) > 0

        @assert size(prop_df, 2) == 2 + dim_factor * 2
        update_proportions!(prop_df, F_model, L, prec, prop_header)

        prop_labels = names(prop_df)
        props = Matrix(prop_df)
        push!(labs, prop_labels => props')
    end


    update_df!(df, labs)
end

function adjust_correlation!(C::AbstractArray{Float64, 3},
                             V::AbstractArray{Float64, 3})
    n = size(C, 3)
    for i = 1:n
        C[:, :, i] .= cov2cor(V[:, :, i])
    end
end

function cov2cor(Σ::AbstractMatrix{Float64})
    D = sqrt.(inv(Diagonal(Σ)))
    return D * Σ * D
end




import Base:length
function length(j::JointParameters)
    return length(j.tree_dims)
end

function is_factor(j::JointParameters, ind::Int)
    return j.tree_dims[ind] != j.data_dims[ind]
end

function rotate_multi_sem(log_path::String, args...; burnin::Float64 = 0.0, kwargs...)
    df = import_log(log_path, burnin=burnin)
    rotate_multi_sem!(df, args...; kwargs...)
end

function rotate_multi_sem!(df::DataFrame, rotate_path::String,
        plan::RotationPlan, parameters::JointParameters;
        burnin::Float64 = 0.0,
        minimum_map::Float64 = max(0.1 - burnin, 0.0),
        use_map::Bool = true, # maximum a posteriori estimate
        kw_args...)

    @unpack tree_dims, data_dims = parameters
    q = sum(tree_dims)
    m = length(parameters)

    # df = import_log(log_path, burnin=burnin)
    n = size(df, 1)

    map_ind = 0
    if use_map
        map_offset = Int(round(n * minimum_map))
        map_ind = findmax(df.joint[(map_offset + 1):end])[2] + map_offset
    end


    for i = 1:m
        if is_factor(parameters, i)
            rotate_submodel!(df, parameters, plan, i; map_ind = map_ind,
                    kw_args...)
        end
    end
    CSV.write(rotate_path, df, delim = '\t')
    # df = nothing
end

function factor_summaries(df::DataFrame, parameters::JointParameters;
        F_header::AbstractString = parameters.joint_name,
        kwargs...)

    m = length(parameters)
    s = size(df, 1)


    dim_joint = sum(parameters.tree_dims)
    n_taxa = parameters.n_taxa

    fac_inds = findall(parameters.tree_dims .!= parameters.data_dims)

    n = sum(parameters.tree_dims[fac_inds] .* parameters.data_dims[fac_inds])

    F, F_labels = collect_matrix(df, F_header, dim_joint, n_taxa,
    transpose = true)

    taxon_regex = Regex("$(F_header)\\.(.+)\\.(\\d+)\$")
    taxa = [match(taxon_regex, x)[1] for x in F_labels[:, 1]]
    for i in 1:dim_joint
        for j in 1:n_taxa
            @assert F_labels[j, i] == "$F_header.$(taxa[j]).$i"
        end
    end

    ydf = DataFrame(model = Vector{String}(undef, n),
                    factor = fill(0, n),
                    trait = fill(0, n),
                    mean = fill(NaN, n),
                    max = fill(NaN, n),
                    min = fill(NaN, n)
                    )

    fac_offset = 0
    y_offset = 0



    for i in fac_inds
        trait_name = parameters.trait_names[i]
        L_header = "$(trait_name).$L_HEADER" # TODO: allow for different headers

        dim_factor = parameters.tree_dims[i]

        factor_inds = fac_offset .+ (1:dim_factor)

        dim_trait = parameters.data_dims[i]

        fac_offset += dim_factor

        L, L_labels = collect_matrix(df, L_header, dim_trait, dim_factor)
        @show size(L)
        @show size(F)

        Fi = @view F[:, factor_inds, :]
        F_mean = reshape(mean(Fi, dims = 1), dim_factor, s)
        F_max = reshape(maximum(Fi, dims = 1), dim_factor, s)
        F_min = reshape(minimum(Fi, dims = 1), dim_factor, s)


        for j = 1:dim_factor

            for k = 1:dim_trait
                y_offset += 1

                μ = mean(L[k, j, :] .* F_mean[j, :])
                mx = mean(L[k, j, :] .* F_max[j, :])
                mn = mean(L[k, j, :] .* F_min[j, :])

                ydf[y_offset, :model] = trait_name
                ydf[y_offset, :factor] = j
                ydf[y_offset, :trait] = k
                ydf[y_offset, :mean] = μ
                ydf[y_offset, :max] = mx
                ydf[y_offset, :min] = mn

            end
        end





    end

    return ydf


end

function post_process(log_path::String,
                      rotated_path::String,
                      traits::Vector{Pair{String, Tuple{Int, Int}}},
                      n_taxa::Int;
                      optimization_inds::Vector{Vector{Int}} = fill(Int[], length(traits)),
                      rotation_plan = RotationPlan(SVDRotation,
                                                   ProcrustesRotation),
                      joint_name::String = "",
                      use_map::Bool = true, # MAP = maximum a posteriori estimate
                      burnin::Float64 = 0.0,
                      kwargs...
                      )

    n = length(traits)
    tree_dims = [x[2][1] for x in traits]
    data_dims = [x[2][2] for x in traits]
    trait_names = [x[1] for x in traits]
    joint_name = isempty(joint_name) ? join(trait_names, '.') * ".joint" : joint_name

    params = JointParameters(tree_dims, data_dims, trait_names, joint_name, n_taxa)

    df = import_log(log_path, burnin = burnin)

    @show kwargs


    rotate_multi_sem!(df, rotated_path,
           rotation_plan,
           params,
           double_check = true,
           use_map = use_map;
           optimization_inds = optimization_inds,
           burnin = burnin,
           kwargs...)

    return factor_summaries(df, params)

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
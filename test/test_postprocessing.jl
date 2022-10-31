using Test
using PhylogeneticFactorAnalysis
using PhylogeneticFactorAnalysis.BEASTPostProcessing
using PhylogeneticFactorAnalysis.PostProcessing
using UnPack
using BeastUtils.MatrixUtils
using BeastUtils.Logs
using LinearAlgebra
using DataFrames

BPP = PhylogeneticFactorAnalysis.BEASTPostProcessing
absmax = BPP.absmax

const T = Vector{Pair{String, Tuple{Int, Int}}}


function loadings_and_factors(df::DataFrame, dims::T, N::Int)

    ks = [x[2][1] for x in dims]
    ps = [x[2][2] for x in dims]
    K = sum(ks)
    P = sum(ps)
    S = size(df, 1)
    L = zeros(P, K, S)
    F = zeros(N, K, S)

    m = length(dims)
    k_offset = 0
    p_offset = 0
    for i = 1:m
        trait = dims[i][1]
        k, p = dims[i][2]
        if k == p
            Li = Diagonal(ones(k))
            for j = 1:S
                L[p_offset .+ (1:p), k_offset .+ (1:k), j] .= Li
            end
        else
            Li, _ = BPP.collect_matrix(df, "$trait.L", p, k)
            L[p_offset .+ (1:p), k_offset .+ (1:k), :] .= Li
        end
        k_offset += k
        p_offset += p
    end

    trait_names = [x[1] for x in dims]
    F, _ = BPP.collect_matrix(df, join(trait_names, '.') * ".joint", K,
                               N, transpose=true)

    C, _ = BPP.collect_matrix(df, "correlation.mbd.variance", K, K)
    V, _ = BPP.collect_matrix(df, "mbd.variance", K, K)
    return (L = L, F = F, C = C, V = V)
end

function loadings_and_factors(path::String, dims::T, n::Int;
                              burnin::Float64 = 0.0)
    df = import_log(path, burnin=burnin)
    return loadings_and_factors(df, dims, n)
end

import PhylogeneticFactorAnalysis.PostProcessing.check_valid_rotation
function check_valid_rotation(p1::String, p2::String, dims::T, n::Int;
                              b1::Float64 = 0.0, b2::Float64 = 0.0,
                              tol::Float64 = 1e-9)
    lf1 = loadings_and_factors(p1, dims, n, burnin = b1)
    lf2 = loadings_and_factors(p2, dims, n, burnin = b2)
    L1, F1, C1, V1 = lf1.L, lf1.F, lf1.C, lf1.V
    L2, F2, C2, V2 = lf2.L, lf2.F, lf2.C, lf2.V

    p, k, n = size(L1)
    _, q, m = size(F1)
    @assert m == n
    @assert q == k

    for i = 1:n
        Y1 = F1[:, 1:k, i] * L1[:, :, i]'
        Y2 = F2[:, 1:k, i] * L2[:, :, i]'
        diff = BPP.absmax(Y1 - Y2)
        if diff > tol
            error("The maximum difference between the two Y_hat values is " *
                  "$diff. Maximum allowable difference is $tol.")
        end
    end

    offset = 0
    for dim in dims
        kd, pd = dim[2]
        if kd == pd
            diff = absmax(F1[:, offset .+ (1:kd), :] - F2[:, offset .+ (1:kd), :])
            if diff > tol
                error("The maximum difference between the unrotated values is " *
                      "$diff. Maximum allowable difference is $tol.")
            end
        end
        offset += kd
    end


    for i = 1:n
        P1 = inv(V1[:, :, i])
        P2 = inv(V2[:, :, i])
        F1i = F1[:, :, i]
        S1 = F1i * P1 * F1i' / n

        F2i = F2[:, :, i]
        S2 = F2i * P2 * F2i' / n

        diff = absmax(S1 - S2)
        if diff > tol
            display(S1)
            display(S2)
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

    return true
end


log_path = joinpath(@__DIR__, "..", "examples", "post_processing", "test.log")
rotated_path = joinpath(@__DIR__, "rotated.log")


dims = ["fac" => (3, 10), "trait" => (1, 1)]
n_taxa = 20


post_process(log_path, rotated_path, dims, n_taxa)
@test check_valid_rotation(log_path, rotated_path, dims, n_taxa)

post_process(log_path, rotated_path, dims, n_taxa, optimization_inds=[[4]])
@test check_valid_rotation(log_path, rotated_path, dims, n_taxa)

# post_process(log_path, rotated_path, dims, n_taxa,
#              rotation_plan = RotationPlan(SVDRotation, Signs))
# @test check_valid_rotation(log_path, rotated_path, dims, n_taxa)

GC.gc() #TODO: this should get unlinked earlier
rm(rotated_path)
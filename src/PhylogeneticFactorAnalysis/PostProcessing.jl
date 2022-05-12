module PostProcessing

export svd_logs

using LinearAlgebra, Statistics, LinearAlgebra.BLAS
using BeastUtils.Logs, BeastUtils.MatrixUtils, BeastUtils.PosteriorSummary


const DEBUG = false

function my_show(x, s::String)
    println("$s:")
    display(x)
    println()
end

function process_log(log_path::String, header::Union{Regex, String}, k::Int,
                    p::Int;
                    fac_header::String = "factors.",
                    proportion_header::String = "factorProportion.",
                    rotate_factors::Bool = false,
                    find_best::Bool = true,
                    do_svd::Bool = true,
                    rotation_cols::Vector{Int} = zeros(Int, k),
                    transposed::Bool = false,
                    relevant_cols::Vector{Int} = collect(1:p),
                    relevant_rows::Vector{Int} = collect(1:k))

    cols, data = get_log(log_path, burnin = 0.0)
    L_inds = findall(x -> startswith(x, header), cols)

    @assert length(L_inds) == k * p

    f_inds = rotate_factors ? findall(x -> startswith(x, fac_header), cols) : Int[]
    prop_inds = findall(x -> startswith(x, proportion_header), cols)
    prop_inds = isnothing(prop_inds) ? Int[] : prop_inds
    compute_prop = !isnothing(prop_inds) && rotate_factors

    n_taxa = div(length(f_inds), k)
    @assert length(f_inds) == n_taxa * k

    L_cols = cols[L_inds]
    L_data = data[:, L_inds]

    f_cols = cols[f_inds]
    f_data = data[:, f_inds]

    prop_cols = cols[prop_inds]
    prop_data = data[:, prop_inds]
    all_prop_ind = 1
    abs_prop_inds = 2:(k + 1)
    rel_prop_inds = (k + 2):(2 * k + 1)
    marg_prop_ind = 2 * k + 2

    if compute_prop
        @assert length(prop_inds) == marg_prop_ind
        @assert prop_cols[all_prop_ind] == proportion_header * "factorProportion"
        for i = 1:k
            @assert prop_cols[abs_prop_inds[i]] == proportion_header * "absoluteProportion.$i"
            @assert prop_cols[rel_prop_inds[i]] == proportion_header * "relativeProportion.$i"
        end
        @assert prop_cols[marg_prop_ind] == proportion_header * "relativeMarginalProportion"
    end

    n = size(data, 1)


    d_storage = zeros(n, k)
    v_storage = zeros(n, k * p)
    l_storage = zeros(n, k * p)
    f_storage = zeros(n, n_taxa * k)
    prop_storage = zeros(n, 2 * k + 2)

    ke = length(relevant_rows)
    pe = length(relevant_cols)

    do_svd = do_svd && ke != 0

    other_rows = setdiff(1:k, relevant_rows)
    other_cols = setdiff(1:k, relevant_cols)

    L_full = zeros(k, p)
    Ft_full = zeros(k, n_taxa)
    UtFt_full = zeros(k, n_taxa)

    L = zeros(ke, pe)
    Ft = zeros(ke, n_taxa)
    UtFt = zeros(ke, n_taxa)

    fac_inds = vec([(i - 1) * k + row for row in relevant_rows, i in 1:n_taxa])

    S = zeros(ke)
    Vt = zeros(ke, pe)

    for i = 1:n
        fill_L!(L_full, L_data, i, transposed)
        for j = 1:pe
            for l = 1:ke
                L[l, j] = L_full[relevant_rows[l], relevant_cols[j]] # TODO: add inbounds
            end
        end

        if DEBUG
            my_show(L_full, "L_full")
            my_show(L, "L")
        end


        if do_svd
            s = svd(L)
            S .= s.S
            Vt .= s.Vt
        else
            for i = 1:ke
                S[i] = norm(@view L[i, :])
                for j = 1:pe
                    Vt[i, j] = L[i, j] / S[i]
                end
            end
        end

        if do_svd
            for j = 2:ke
                @assert S[j - 1] >= S[j]
            end
        end
        # rotate_svd!(s, cols = rotation_cols)
        if DEBUG
            L2 = U * Diagonal(S) * Vt
            @assert maximum(abs.(L2 - L)) < 1e-12
        end

        for j = 1:ke
            d_storage[i, relevant_rows[j]] = S[j]
            offset = (relevant_rows[j] - 1) * p
            for l = 1:pe
                col = relevant_cols[l]
                ind = offset + col
                v_storage[i, ind] = Vt[j, l]
                l_storage[i, ind] = Vt[j, l] * S[j]
            end
        end
        for j = 1:(k - ke)
            row = other_rows[j]
            l = @view L_full[row, :]
            nrm = norm(l)
            d_storage[i, row] = nrm

            rng = ((row - 1) * p + 1):(row * p)
            l_storage[i, rng] .= l
            if nrm != 0
                v_storage[i, rng] .= l / nrm # TODO: memory efficiency
            end
        end

        #
        # d_storage[i, :] .= s.S
        # v_storage[i, :] .= vec(s.Vt')
        #
        # l_storage[i, :] .= vec(s.Vt' * Diagonal(s.S))

        if rotate_factors
            U = Matrix(Diagonal(ones(ke)))
            if do_svd
                U = s.U
            end
            copyto!(Ft, @view f_data[i, fac_inds])
            gemm!('T', 'N', 1.0, U, Ft, 0.0, UtFt)
            for j = 1:n_taxa
                offset = (j - 1) * k
                for l = 1:ke
                    f_storage[i, offset + relevant_rows[l]] = UtFt[l, j]
                end
            end
            for j = 1:n_taxa
                offset = (j - 1) * k
                for l = 1:(k - ke)
                    ind = offset + other_rows[l]
                    f_storage[i, ind] = f_data[i, ind]
                end
            end
        end

        if compute_prop
            fill_L!(L_full, l_storage, i, transposed)
            copyto!(Ft_full, @view f_storage[i, fac_inds])

            LLt = L_full * L_full' #TODO cache these and below
            FtF = Ft_full * Ft_full'
            V = LLt .* FtF
            v = diag(V)
            sv = sum(v)
            rel_props = v ./ sv
            abs_props = rel_props * prop_data[i, all_prop_ind]


            prop_storage[i, rel_prop_inds] .= rel_props
            prop_storage[i, abs_prop_inds] .= abs_props
            prop_storage[i, all_prop_ind] = prop_data[i, all_prop_ind]

            sV = sum(V)
            prop_storage[i, marg_prop_ind] = sv / sV
        end





        # if i > 2
        #     error("break")
        # end
    end

    force_rotation = false
    if any(rotation_cols .!= 0)
        force_rotation = true
    end
    if !force_rotation
        if !find_best
            fill!(rotation_cols, 1)
        else
            rotation_cols .= find_rotation_cols(v_storage, k, p)
        end
    end
    reflect!(v_storage, l_storage, f_storage, rotation_cols, k, p)
    return d_storage, v_storage, l_storage, f_storage, prop_storage
end

function fill_L!(L::Matrix{Float64}, data::Matrix{Float64}, row::Int,
                 transposed::Bool)
    k, p = size(L)
    if transposed
        for i = 1:p
            offset = k * (i - 1)
            for j = 1:k
                L[j, i] = data[row, offset + j]
            end
        end
    else
        for i = 1:k
            offset = p * (i - 1)
            for j = 1:p
                L[i, j] = data[row, offset + j]
            end
        end
    end
end

function rotate_svd!(s::SVD{Float64,Float64,Array{Float64,2}}; cols::Vector{Int} = ones(Int, length(s.S)))
    k = length(s.S)
    rev = ones(k)
    for i = 1:k
        if s.Vt[i, cols[i]] < 0.0
            rev[i] = -1.0
        end
    end
    Q = Diagonal(rev)
    s.U .= s.U * Q #TODO: make memory efficient
    s.Vt .= Q * s.Vt #TODO: see above
end

function find_rotation_cols(data::Matrix{Float64}, k::Int, p::Int)
    abs_mat = abs.(data)
    means = vec(mean(abs_mat, dims = 1))
    vars = vec(var(abs_mat, dims = 1))

    zs = means ./ sqrt.(vars)
    nan_inds = findnans(zs)
    zs[nan_inds] .= 0.0
    cols = zeros(Int, k)
    for i = 1:k
        l = ((i - 1) * p) + 1
        u = i * p
        zmax = findmax(@view zs[l:u])
        cols[i] = zmax[2]
    end
    return cols
end

function process_loadings(log_path::String, header::Union{Regex, String}, k::Int, p::Int; find_best::Bool = true)
    cols, data = get_log(log_path, header, burnin = 0.0)

    @assert length(cols) == k * p

    rotation_cols = ones(Int, k)
    if find_best
        rotation_cols .= find_rotation_cols(data, k, p)
    end
    reflect!(data, rotation_cols, k, p)
    return data, cols
end

function reflect!(data::Matrix{Float64}, loadings::Matrix{Float64},
                    factors::Matrix{Float64}, cols::Vector{Int}, k::Int, p::Int)
    @assert length(cols) == k
    n, kp = size(data)
    @assert kp == k * p

    @assert size(factors, 1) == n
    n_taxa = div(size(factors, 2), k)
    @assert n_taxa * k == size(factors, 2)

    for i = 1:k
        offset = (i - 1) * p
        check_ind = offset + cols[i]
        fac_inds = i:k:(n_taxa * k)
        for j = 1:n
            if data[j, check_ind] < 0.0
                for k = (offset + 1):(i * p)
                    data[j, k] = -data[j, k]
                    loadings[j, k] = -loadings[j, k]
                end
                for ind in fac_inds
                    factors[j, ind] = -factors[j, ind]
                end
            end
        end
    end
    return data
end


function svd_logs(path::String, new_path::String, k::Int, p::Int;
        cols::Vector{Int} = zeros(Int, k),
        Lid::String = "L",
        fid::String = "factors.",
        propid::String = "factorProportion.",
        rotate_factors::Bool = false,
        transposed::Bool = false,
        do_svd::Bool = true,
        relevant_rows::Vector{Int} = collect(1:k),
        relevant_cols::Vector{Int} = collect(1:p))

    d, v, l, f, props = process_log(path, Lid, k, p, rotation_cols = cols,
                                transposed = transposed,
                                fac_header = fid,
                                proportion_header = propid,
                                rotate_factors = rotate_factors,
                                relevant_rows = relevant_rows,
                                relevant_cols = relevant_cols,
                                do_svd = do_svd)
    d_labels = ["sv$i" for i = 1:k]
    v_labels = vec(["V$i$j" for j = 1:p, i = 1:k])
    l_labels = vec(["L$i$j" for j = 1:p, i = 1:k])
    cols = get_cols(path)
    f_inds = rotate_factors ? findall(x -> startswith(x, fid), cols) : Int[]
    f_labels = string.(cols[f_inds])

    prop_inds = findall(x -> startswith(x, propid), cols)
    prop_inds = rotate_factors && !isnothing(prop_inds) ? prop_inds : Int[]
    prop_labels = string.(cols[prop_inds])


    states = get_log(path, "state", burnin = 0.0)[2]


    labels = ["state"; d_labels; v_labels; l_labels; f_labels; prop_labels]
    data = [states d v l f props]


    make_log(new_path, data, labels, includes_states = true)
end

# function svd_logs(path::String, new_path::String;
#     prec_start = "factorPrecision",
#     Lid::String = "L",
#     fid::String = "factors.",
#     rotate_factors::Bool = false,
#     transposed::Bool = false)

#     cols = Logs.get_cols(path)
#     p = length(findall(x -> startswith(x, prec_start), cols))
#     kp = length(findall(x -> startswith(x, Lid), cols))
#     k, r = divrem(kp, p)
#     if r != 0
#         error("Unable to determine the number of factors and traits.")
#     end

#     return svd_logs(path, new_path, k, p, Lid = Lid, fid = fid,
#             rotate_factors = rotate_factors, transposed = transposed)
# end



# k = 2
# p = 100
#
# path = joinpath(Directories.desktop, "simFactor_N100_P100_K2_gibbs.log")
# new_path = joinpath(Directories.desktop, "gibbs.log")
#
# svd_logs(path, new_path, k, p)

end

function find_starts(x::AbstractArray{<:AbstractString}, s::AbstractString)
    return findall(a -> startswith(a, s), x)
end

function process_log(log_path::String; kw_args...)
    df = import_log(log_path; kw_args...)
    n = size(df, 1)
    nms = names(df)

    sv_inds = find_starts(nms, SV_HEADER)
    k = length(sv_inds)

    L_inds = find_starts(nms, L_HEADER)
    p = div(length(L_inds), k)
    L = zeros(p, k, n)

    F_inds = find_starts(nms, FAC_HEADER)
    n_taxa = div(length(F_inds), k)
    F = zeros(n_taxa, k, n)

    for i = 1:n
        L[:, :, i] .= reshape(Vector(df[i, L_inds]), p, k)
        F[:, :, i] .= reshape(Vector(df[i, F_inds]), k, n_taxa)'
    end

    return (L = L, F = F)
end





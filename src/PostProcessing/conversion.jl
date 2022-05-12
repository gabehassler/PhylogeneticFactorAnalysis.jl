function compose(args::Transform...)
    n = length(args)
    t = args[1]
    for i = 2:n
        t = pairwise_compose(t, args[i])
    end

    return t
end

function pairwise_compose(t1::Transform, t2::Transform)
    pairwise_compose(Rotations(t1), Rotations(t2))
end

function pairwise_compose(r1::Rotations, r2::Rotations)
    R = apply_transform(r1.rotations, r2)
    return Rotations(R)
end

function Rotations(r::Rotations)
    return r
end

function Rotations(signs::Signs)
    k, n = size(signs.signs)
    R = zeros(k, k, n)
    for i = 1:n
        for j = 1:k
            R[j, j, i] = signs.signs[j, i]
        end
    end
    return Rotations(R)
end

function Rotations(ca::ClusterAssignments)
    k, n = size(ca.clusters)
    R = zeros(k, k, n)
    I = Matrix(Diagonal(ones(k)))
    for i = 1:n
        p = @view ca.clusters[:, i]
        signs = Diagonal(@view ca.signs[:, i])
        R[:, :, i] = I[:, p] * signs
    end

    return Rotations(R)
end

function expand_rotation(r::Rotations, dim::Int)
    # _, k, n = size(r.rotations)
    # k2 = dim + k
    # R = zeros(k2, k2, n)
    # for i = 1:n
    #     R[1:k, 1:k, i] .= r.rotations[:, :, i]
    #     for j = (k + 1):k2
    #         R[j, j, i] = 1.0
    #     end
    # end

    # return Rotations(R)
    _, k, _ = size(r.rotations)

    expand_rotation(r, 1:k, k + dim)
end

function expand_rotation(r::Rotations, model_inds::AbstractVector{Int},
        joint_dim::Int)

    _, _, n = size(r.rotations)
    R = zeros(joint_dim, joint_dim, n)
    other_inds = setdiff(1:joint_dim, model_inds)
    for i = 1:n
        R[model_inds, model_inds, i] .= r.rotations[:, :, i]
        for j in other_inds
            R[j, j, i] = 1.0
        end
    end

    return Rotations(R)
end


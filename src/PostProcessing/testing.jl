const TOL = 1e-12

function absmax(X::Array{<:Real})
    return maximum(abs.(X))
end

function is_ortho(x::AbstractMatrix{Float64}; tol::Float64 = TOL)
    return absmax(I - x' * x) < tol
end

function check_valid_rotation(t::Transform)
    return check_valid_rotation(Rotations(t))
end

function check_valid_rotation(r::Rotations)
    for i = 1:size(r.rotations, 3)
        x = @view r.rotations[:, :, i]
        if !is_ortho(x)
            return false
        end
    end

    return true
end

function check_transform_results(orignal_data::Array{Float64, 3},
                                 rotated_data::Array{Float64, 3},
                                 t::Transform;
                                 tol::Float64 = TOL)
    return absmax(rotated_data - apply_transform(orignal_data, t)) < tol
end
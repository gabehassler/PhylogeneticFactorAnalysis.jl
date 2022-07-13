function postprocess_log(input::AbstractString, output::AbstractString,
                     k::Int, p::Int, n::Int,
                     rows::AbstractVector{Int}, cols::AbstractVector{Int},
                     rotation_plan::RotationPlan)
    if rows != 1:k || cols != 1:p
        error("New rotation shceme not implemented for hybrid constraint")
    end

    traits = ["traits" => (k, p)]

    BEASTPostProcessing.post_process(input, output, traits, n,
                                     rotation_plan = rotation_plan,
                                     prec_header = "factorPrecision",
                                     F_header = "factors.",
                                     L_header = "L",
                                     prop_header = "factorProportion",
                                     use_map=false)
end

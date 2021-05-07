const example_dir = joinpath(@__DIR__, "..", "examples")

const aqui_data_continuous = joinpath(example_dir, "aquilegia_continuous.csv")
const aqui_data_discrete = joinpath(example_dir, "aquilegia_with_discrete.csv")
const aqui_newick = joinpath(example_dir, "aquilegia_newick.txt")

const xml0 = joinpath(example_dir, "example_0.xml")
const xml1 = joinpath(example_dir, "example_1.xml")
const xml2 = joinpath(example_dir, "example_2.xml")
const xml3 = joinpath(example_dir, "example_3.xml")
const xml4 = joinpath(example_dir, "example_4.xml")


const EXAMPLES = Dict(
        0 => [xml0, aqui_data_continuous, aqui_newick],
        1 => [xml1, aqui_data_continuous, aqui_newick],
        2 => [xml2, aqui_data_discrete, aqui_newick],
        3 => [xml3, aqui_data_continuous, aqui_newick],
        4 => [xml4, aqui_data_continuous, aqui_newick]
        )


function import_example(n::Int; dir::String = pwd(), force::Bool = false)
    new_dir = joinpath(dir, "pfa_example_$n")
    mkpath(new_dir)
    for file in EXAMPLES[n]
        cp(file, joinpath(new_dir, basename(file)), force = force)
    end
    return new_dir
end

function run_example(n::Int; args...)
    example_dir = import_example(n; args...)
    pfa(joinpath(example_dir, "example_$n.xml"))
end

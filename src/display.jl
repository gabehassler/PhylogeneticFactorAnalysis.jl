import Base.show

function Base.show(io::IO, input::PipelineInput)
    output = "$(typeof(input))"
    output *= "Phylogenetic factor analysis with the following instructions:"
    components = [
        "name" => input.name,
        "directory" => input.directory,
        "prior" => typeof(input.prior)
    ]

    for comp in components
        output *= "\n\t$(comp[1]): $(comp[2])"
    end

    println(io, output)
end
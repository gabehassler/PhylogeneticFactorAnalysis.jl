import Base.show

function Base.show(io::IO, input::PipelineInput)
    output = "$(typeof(input)) - "
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

function Base.show(io::IO, tasks::PipelineTasks)
    fns = fieldnames(PipelineTasks)
    n = length(fns)
    components = Vector{String}(undef, n)
    for i = 1:n
        fn = fns[i]
        val = getfield(tasks, fn)
        components[i] = "\t" * string(fn) * ": $val"
    end
    output = "$(typeof(tasks)) - Tasks to be run in this analysis:\n" *
            join(components, '\n')
    println(io, output)
end


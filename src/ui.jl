function pfa(input_file::String)
    if endswith(input_file, ".xml")
        input = parse_xml(input_file)
    elseif endswith(input_file, ".jld")
        input = load_jld(input_file)
    else
        error("Unknown file extension. Must be either an 'xml' file or 'jld' file.")
    end
    run_pipeline(input)
end
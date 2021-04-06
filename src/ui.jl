function pfa(input_file::String)
    if endswith(input_file, ".xml")
        input = parse_xml(input_file)
    elseif endswith(input_file, ".jld")
        input = load_jld(input_file)
    else
        error("Unknown file extension. Must be either an 'xml' file or 'jld' file.")
    end

    old_wd = pwd()
    try
        run_pipeline(input)
    catch e
        cd(old_wd)
        @error "Something went wrong" exception=(e, catch_backtrace())
    end
end
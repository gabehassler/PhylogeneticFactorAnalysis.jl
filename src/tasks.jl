function check_completion(input::PipelineInput)
    selection_xml_made = check_selection_xml(input)
end


function check_selection_xml(input::PipelineInput)
    error("Not implemented")
end


const BEGINNING = "beginning"

const TASK_FIELDNAMES = [:make_selection_xml,
                         :run_selection_xml,
                         :record_selection_stats,
                         :make_final_xml,
                         :run_final_xml,
                         :process_final_log,
                         :plot_loadings,
                         :plot_factors]

const FROM_FIELDS = Dict{String, Symbol}(
                                         "plots" => :plot_loadings
                                        )



function start_from(s::String, input::PipelineInput)
    @unpack tasks = input

    @assert all(fieldnames(typeof(tasks)) .== TASK_FIELDNAMES)

    start_field = FROM_FIELDS[s]
    start_ind = findfirst(isequal(start_field), TASK_FIELDNAMES)
    for i = 1:(start_ind - 1)
        setfield!(tasks, TASK_FIELDNAMES[i], false)
    end
    for i = start_ind:length(TASK_FIELDNAMES)
        setfield!(tasks, TASK_FIELDNAMES[i], true)
    end

    tasks.make_final_xml = true; @warn "TODO: de-couple names from make_final_xml"
end
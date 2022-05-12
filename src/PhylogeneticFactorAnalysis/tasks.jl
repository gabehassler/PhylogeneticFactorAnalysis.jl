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
                            "plots" => :plot_loadings,
                            "selection_stats" => :record_selection_stats,
                            "process_final" => :process_final_log,
                            "run_final_xml" => :run_final_xml
                                        )

const ONLY_FIELDS = Dict{String, Vector{Symbol}}(
                            "plots" => [:plot_loadings, :plot_factors],
                            "initial_xml" => [:make_selection_xml],
                            "model_selection" => [:make_selection_xml,
                                                  :run_selection_xml,
                                                  :record_selection_stats],
                            "prepare_final_xml" => [:record_selection_stats,
                                                    :make_final_xml]
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

function run_only(s::String, input::PipelineInput)
    @unpack tasks = input

    @assert all(fieldnames(typeof(tasks)) .== TASK_FIELDNAMES)
    do_fields = ONLY_FIELDS[s]
    dont_fields = setdiff(TASK_FIELDNAMES, do_fields)
    for field in do_fields
        setfield!(tasks, field, true)
    end
    for field in dont_fields
        setfield!(tasks, field, false)
    end
end

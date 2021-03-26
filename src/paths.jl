function basedir(input::PipelineInput)
    return joinpath(input.directory, input.name)
end

function statistics_dir(input::PipelineInput)
    return basedir(input)
end

function get_stat_path(input::PipelineInput, stat::String)
    return joinpath(statistics_dir(input), "$stat.csv")
end

function selection_xml_dir(input::PipelineInput)
    return joinpath(basedir(input), "selection_xml")
end

function selection_xml_path(input::PipelineInput; kwargs...)
    return joinpath(selection_xml_dir(input), xml_name(input; kwargs...))
end

function selection_log_dir(input::PipelineInput)
    return joinpath(basedir(input), "selection_logs")
end

function selection_log_path(input::PipelineInput; kwargs...)
    return joinpath(selection_log_dir(input), log_name(input; kwargs...))
end

function timer_dir(input::PipelineInput)
    return joinpath(basedir(input), "timing_files")
end

function timer_path(input::PipelineInput; kwargs...)
    return joinpath(timer_dir(input), timer_name(input; kwargs...))
end


function final_paths_helper(input::PipelineInput, extension::String)
    basenames = keys(input.model_selection.final_names)
    return [x * extension for x in basenames]
end

function final_xml_paths(input::PipelineInput)
    final_paths_helper(input, ".xml")
end

function final_log_paths(input::PipelineInput)
    final_paths_helper(input, ".log")
end

function processed_log_paths(input::PipelineInput)
    final_paths_helper(input, "_processed.log")
end

function loadings_statistic_paths(input::PipelineInput)
    final_paths_helper(input, "_loadingsStatistics.csv")
end

function loadings_plot_paths(input::PipelineInput)
    final_paths_helper(input, "_loadings.pdf")
end

function factors_statistic_paths(input::PipelineInput)
    final_paths_helper(input, "_factorMeans.csv")
end

function factors_plot_paths(input::PipelineInput)
    final_paths_helper(input, "_factors.pdf")
end

import Base.basename
function basename(input::PipelineInput; model::Int = 0, rep::Int = 0, stat="")
    name = input.name

    if model != 0
        name = name * "_model$model"
    end
    if !isempty(stat)
        name = name * "_$stat"
    end
    if rep != 0
        name = name * "_rep$rep"
    end

    return name
end

function xml_name(input::PipelineInput; kwargs...)
    return basename(input; kwargs...) * ".xml"
end

function log_name(input::PipelineInput; kwargs...)
    return basename(input; kwargs...) * ".log"
end

function timer_name(input::PipelineInput; kwargs...)
    return basename(input; kwargs...) * "_timer.txt"
end


function sans_extension(x::String)
    l = findlast(isequal('.'), x)
    return x[1:(l - 1)]
end
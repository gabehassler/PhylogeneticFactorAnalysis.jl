using PhylogeneticFactorAnalysis, Test



check_beast()
@test check_r()


cl = 100
examples = 1:4
for example in examples
    input_dir = "pfa_example_$(example)"
    output_dir = input_dir * "_output"
    @assert !isdir(output_dir)
    @assert !isdir(input_dir)
    run_example(example, override_chainlength = cl)
    @test isdir(input_dir)
    @test isdir(output_dir)
    rm(output_dir, recursive = true)
    rm(input_dir, recursive = true)
end

run_example(0)
input_dir = "pfa_example_0"
output_dir = input_dir * "_output"

input = parse_xml(joinpath(input_dir, "example_0.xml"))
input.overwrite = true

start_from("plots", input)
run_pipeline(input)
tasks = input.tasks
@test tasks.make_selection_xml == tasks.run_selection_xml ==
        tasks.record_selection_stats == tasks.run_final_xml == false

run_only("initial_xml", input)
run_pipeline(input)

@test tasks.make_selection_xml = true
@test tasks.run_selection_xml == tasks.record_selection_stats ==
        tasks.make_final_xml == tasks.run_final_xml ==
        tasks.process_final_log == tasks.plot_loadings == tasks.plot_factors ==
        false

display(input)
display(tasks)

log_path = joinpath(output_dir, output_dir * "_processed.log")
PhylogeneticFactorAnalysis.process_log(log_path)

@test isdir(input_dir)
@test isdir(output_dir)
rm(output_dir, recursive = true)
rm(input_dir, recursive = true)
using PhylogeneticFactorAnalysis, Test


example_dir = "pfa_example_0"
test_dir = "pfa_test"

check_beast()
@test check_r()

run_example(0)

@test isdir(example_dir)
@test isdir(test_dir)

input = parse_xml(joinpath(example_dir, "example_0.xml"))
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
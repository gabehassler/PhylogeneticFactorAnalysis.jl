# Partial Analyses

You don't need to run the entire pipeline every time you run an analysis.
The instructions below have instructions for only running part of the pipeline.

__NOTE: THIS CODE HAS NOT BEEN EXHAUSTIVELY TESTED AND MAY HAVE UNINTENDED CONSEQUENCES IF TRYING TO RESUME WORK FROM AN EXISTING ANALYSIS (LIKE DELETING FILES IN THE SPECIFIED FOLDER).
YOU SHOULD BACK UP YOUR WORK BEFORE TRYING TO START AN ANALYSIS FROM MIDWAY THROUGH__

## Example: Model Selection Only
The following example shows how to set the pipeline to only do model selection and make the final BEAST xml, but not the running the final xml or plotting:

```
julia> import_example(1); # this is just so the example will work, you don't need to do this

julia> input = parse_xml("pfa_example_1/example_1.xml") #path to your own xml file
PhylogeneticFactorAnalysis.PipelineInputPhylogenetic factor analysis with the following instructions:
        name: aquilegia1
        directory: C:\Users\gabeh\OneDrive\Desktop
        prior: PhylogeneticFactorAnalysis.IIDPrior

julia> tasks = input.tasks
PhylogeneticFactorAnalysis.PipelineTasks - Tasks to be run in this analysis:
        make_selection_xml: true
        run_selection_xml: true
        record_selection_stats: true
        make_final_xml: true
        run_final_xml: true
        process_final_log: true
        plot_loadings: true
        plot_factors: true


julia> tasks.run_final_xml = false; tasks.process_final_log = false; tasks.plot_loadings = false; tasks.plot_factors = false
false

julia> tasks
PhylogeneticFactorAnalysis.PipelineTasks - Tasks to be run in this analysis:
        make_selection_xml: true
        run_selection_xml: true
        record_selection_stats: true
        make_final_xml: true
        run_final_xml: false
        process_final_log: false
        plot_loadings: false
        plot_factors: false

julia> run_pipeline(input)
...
```

## Example: Plotting only

You can also use this to re-run or start from an existing analysis.
If you want to re-run the final analysis with a greater chain length, you can adjust the `chanLength` attribute in your xml then do the following.
Note that you must have the atrribute `overwrite="true"` in the `phylogeneticFactorAnalysis` element (i.e. `<phylogeneticFactorAnalysis overwtie="true" ...>`).
Also, (and I'm planning to fix this at some point) you have to have `make_final_xml = true` in the `tasks` object if you want to do anything after that point like plotting.

For example, to just re-run the plots from an existing analysis, you could do:
```
julia> input = parse_xml("pfa_example_1/example_1.xml") #path to your own xml file
PhylogeneticFactorAnalysis.PipelineInputPhylogenetic factor analysis with the following instructions:
        name: aquilegia1
        directory: C:\Users\gabeh\OneDrive\Desktop
        prior: PhylogeneticFactorAnalysis.IIDPrior

julia> tasks = input.tasks
PhylogeneticFactorAnalysis.PipelineTasks - Tasks to be run in this analysis:
        make_selection_xml: true
        run_selection_xml: true
        record_selection_stats: true
        make_final_xml: true
        run_final_xml: true
        process_final_log: true
        plot_loadings: true
        plot_factors: true


julia> tasks.make_selection_xml = false; tasks.run_selection_xml = false; tasks.record_selection_stats = false; tasks.run_final_xml = false;

julia> tasks
PhylogeneticFactorAnalysis.PipelineTasks - Tasks to be run in this analysis:
        make_selection_xml: false
        run_selection_xml: false
        record_selection_stats: false
        make_final_xml: true
        run_final_xml: false
        process_final_log: true
        plot_loadings: true
        plot_factors: true


julia> run_pipeline(input)
...
```
# Examples

Example data and XML files can be examined and run through the `import_example` function.
The examples are numbered 1-4.
To access the first set of example files, run:
```
import_example(1)
```
This will create a directory `pfa_example_1` in your working directory and populate it with the relevant files.
The files associated with Example 1 are `aquilegia_continuous.csv` containing the trait data, `aquilegia_newick.txt` containing the phylogenetic tree, and `example_1.xml` containing the instructions for running PFA.
To run this example, enter
```
pfa("pfa_example_1/example_1.xml")
```

The table below provides a quick summary of the instructions provided in each example XML file.

| Example | Data Type | Prior | Model Selection |
|---|---|---|---|
1 | continuous | i.i.d. Gaussian | yes |
2 | continuous & discrete | i.i.d. Gaussian | yes |
3 | continuous | shrinkage | yes |
4 | continuous | i.i.d. Gaussian | no |

# PhylogeneticFactorAnalysis.jl Documentation

View this repository on GitHub: <https://github.com/gabehassler/PhylogeneticFactorAnalysis.jl>

If you use this software in published research, please cite:

Hassler, Gabriel W., Brigida Gallone, Leandro Aristide, William L. Allen, Max R. Tolkoff, Andrew J. Holbrook, Guy Baele, Philippe Lemey, and Marc A. Suchard. "Principled, practical, flexible, fast: a new approach to phylogenetic factor analysis." _Methods in Ecology and Evolution_ (2022).

You can find the paper here: <https://doi.org/10.1111/2041-210X.13920>


## Installation Instructions

### Install Java

PhylogeneticFactorAnalysis will not automatically install Java, which is required to run BEAST.
To check if you have Java installed, run `java -version` from your computer's command line (not in Julia).
If Java is installed, you will see something like:
```
$ java -version
openjdk version "1.8.0_282"
OpenJDK Runtime Environment (build 1.8.0_282-8u282-b08-0ubuntu1~18.04-b08)
OpenJDK 64-Bit Server VM (build 25.282-b08, mixed mode)
```

If Java is not installed, follow the instructions on the [BEAST website](https://beast.community/installing) for installin Java.
__Note that you ONLY need to install Java. You do NOT need to install BEAST, which comes with PhylogeneticFactorAnalysis.__



### Optional: Set up R Environment
PhylogeneticFactorAnalysis uses R to generate automated plots via the Julia package [RCall](https://juliainterop.github.io/RCall.jl/stable/).
If you are interested in producing these plots, we recommend setting up your R environment with the appropriate packages and ensuring RCall can locate your current R installation.

If you do not have R installed already, you can download and install it from the [R website](https://www.r-project.org/).

Once you have installed R, ensure the appropriate packages are installed by typing the following into the R console:
```
install.packages(c("ggplot2", "wesanderson", "colorspace", "tidyr", "ggtree", "phytools", "tidytree", "aplot", "RColorBrewer", "ggnewscale", "phyclust", "treeio"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("treeio")
BiocManager::install("ggtree")
```

### Add `PhylogeneticFactorAnalysis` to your Julia environment
From the Julia REPL, run
```
import Pkg; Pkg.add("PhylogeneticFactorAnalysis")
```




### Check dependencies

Once you have added `PhylogeneticFactorAnalysis`, check that the dependencies are properly installed.

To check Java and BEAST, enter into the Julia REPL
```
using PhylogeneticFactorAnalysis

check_beast()
```

To check R, add the following line to the code above:
```
check_r()
```
It's OK (and likely inevitable) if `check_r()` results in several R warnings.
As long as it does not throw an error and terminates with `all necessary R packages installed` your R environment should be set up appropriately.

## Instructions

PhylogeneticFactorAnalysis takes a relatively simple xml file and uses it to build and run BEAST xml files.
Users must specify the following:
1. a csv file containing the trait data
2. newick-formatted tree file
3. prior and constraint on the loadings matrix

Here is the simplest xml you can supply
```
<?xml version="1.0" encoding="UTF-8"?>
<phylogeneticFactorAnalysis
    name="test" <!-- this will be the header of all folders and files created -->
    standardizeTraits="true" <!-- we recommend standardizing traits, but you may set to false if desired-->
    >
  <data traits="test_traits.csv" <!-- path to your traits file -->
        tree="test_newick.txt" <!-- path to your newick file -->
        />
  <iidPrior constraint="upperTriangular" <!-- this is the specific prior you use on the loadings. the options are "upperTriangular", "orthogonal", or "hybrid". -->
        />
  <modelSelection>
    <nFactors>5</nFactors> <!-- this is the number of factors you want to use. you may supply as many as you want, (e.g. <nFactors>1 2 3 4 5</nFactors>). If you only supply one it will skip model selection and simply choose that number of factors for the final analysis-->
  </modelSelection>
</phylogeneticFactorAnalysis>
```
The above xml performs a single BEAST run with a 5-factor model and upper-triangular constraint on the loadings.

Here is a more complex xml that can handle discrete traits:
```
<?xml version="1.0" encoding="UTF-8"?>
<phylogeneticFactorAnalysis name="aquilegiaBinary" partitionSeed="666" mcmcSeed="666" overwrite="true" standardizeTraits="true">
  <data traits="aquilegia_binary.csv" tree="aquilegia_newick.txt" discreteIndices="11 12 13 14"/>
  <iidPrior constraint="orthogonal"/>
  <modelSelection repeats="5" selectionStatistic="CLPD" burnin="0.25">
    <nFactors>1 2 3 4 5</nFactors>
    <mcmcOptions chainLength="100000"/>
  </modelSelection>
  <mcmcOptions chainLength="100000"/>
</phylogeneticFactorAnalysis>
```


### Explanation of elements and attributes
Element: `phylogeneticFactorAnalysis`
 - This is the root element that encompasses all other elements.
 - Attributes
   - `name`: the name of the folder and files that will be created. This can be set to anything
   - `directory`: the directory you want to store the results in. You can also put your data files here and the pipeline will find them.
   - `partitionSeed` (optional): the random number seed used for creating the training and test data sets. Can be any positive integer.
   - `mcmcSeed` (optional): the random number seed for the BEAST MCMC simulation
   - `overwrite` (optional): allows you to overwrite any existing folders. Possible values are "true" or "false". __If this is set to "true", it will automatically delete everything in the folder specified by `name`.__
   - `standardizeTraits` (optional): specifies whether data should be standardized (i.e. centered and scaled) prior to analysis. Possible values are "true" or "false".

Element: `data`
 - Locations and instructions related to the input data.
 - Attributes
   - `traits`: path to csv file where your trait data is stored. The csv file must have the first column corresponding to the taxa with column header "taxon".
   - `tree`: path to the newick-formatted tree file.
   - `discreteIndices` (optional): if any traits are discrete, list the indices that correspond to the discrete traits.

Element: `iidPrior`
 - This is the prior. You must specify one prior, but you may specify `orthogonalShrinkagePrior`  (below) instead of `iidPrior`.
 - Attributes
   - `constraint`: the constraint on the loadings matrix. Options are "upperTriangular", "orthogonal", or "hybrid".

Element: `orthogonalShrinkagePrior`
 - This is the prior. You must specify one prior, but you may specify `iidPrior` (above) instead of `orthogonalShrinkagePrior`.
 - Attributes: __TODO__

Element: `modelSelection`
 - Supplies the models you want to compare to each other.
 - Attributes
   - `repeats`: Specifies the "k" value for k-fold cross validation. Should be a positive integer > 1 (although we wouldn't recommend going below 5).
   - `selectionStatistic` (optional): the specific selection statistic you would like to use. Values are "CLPD" or "MSE".
   - `burnin` (optional): the MCMC burnin before evaluating the selection statistic. Should be any number between 0 and 1. Defaults to 0.25.
 - Sub-element: `nFactors`
   - specify the number of factors for each model
 - Sub-element: `mcmcOptions` (optional)
   - Meta-parameters related to the MCMC simulation
   - Attributes
     - `chainLength`: the number of states to run the MCMC chain. Defaults to `10000`.

Element: `mcmcOptions`
 - Meta-parameters related to the final MCMC simulation


## Running analyses
To load the package, type `using PhylogeneticFactorAnalysis` into the Julia REPL.

After loading the package, the only function required to run analyses is `pfa`.
If `test.xml` is the path to your xml file, simply enter `pfa("test.xml")` into the Julia REPL to run the pipeline.


## Multithreading
This package uses multi-threading to speed up analyses.
Essentially, it runs many BEAST runs in parallel rather than sequentially.
To take advantage of the multi-threading capabilities, you must start Julia with multiple threads.

By default, Julia starts with 1 thread.
To check the number of threads, call `Threads.nthreads()` from within Julia.
```
julia> Threads.nthreads()
1
```

To start with more threads, either set the `JULIA_NUM_THREADS` environment variable on your local machine to the number of threads you want, or start Julia from the command line as follows:
```
$ julia --threads 12
```
Note that you can replace the `12` with the number of threads you want.
Using `--threads auto` lets Julia decide how many threads to use based on your local computer.
You can find more information on setting the number of available threads here: <https://docs.julialang.org/en/v1/manual/multi-threading/>.

To avoid memory issues (and otherwise overburdening your machine), the number of active threads is limited to the number of cross-validation training/test sets.
If you're doing 5-fold cross validation, for example, then there will be a maximum of 5 active threads.
You can set this via the `repeats` attribute in the `modelSelection` xml element.



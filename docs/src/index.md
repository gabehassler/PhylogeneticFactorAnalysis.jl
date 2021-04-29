# PhylogeneticFactorAnalysis.jl Documentation

## Installation instructions

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
```

### Add `PhylogeneticFactorAnalysis` to your Julia environment
From the Julia REPL, run
```
import Pkg; Pkg.add(url="https://github.com/gabehassler/PhylogeneticFactorAnalysis.jl.git")
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


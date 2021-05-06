var documenterSearchIndex = {"docs":
[{"location":"#PhylogeneticFactorAnalysis.jl-Documentation","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"","category":"section"},{"location":"#Installation","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"Installation","text":"","category":"section"},{"location":"#Install-Java","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"Install Java","text":"","category":"section"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"PhylogeneticFactorAnalysis will not automatically install Java, which is required to run BEAST. To check if you have Java installed, run java -version from your computer's command line (not in Julia). If Java is installed, you will see something like:","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"$ java -version\nopenjdk version \"1.8.0_282\"\nOpenJDK Runtime Environment (build 1.8.0_282-8u282-b08-0ubuntu1~18.04-b08)\nOpenJDK 64-Bit Server VM (build 25.282-b08, mixed mode)","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"If Java is not installed, follow the instructions on the BEAST website for installin Java. Note that you ONLY need to install Java. You do NOT need to install BEAST, which comes with PhylogeneticFactorAnalysis.","category":"page"},{"location":"#Optional:-Set-up-R-Environment","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"Optional: Set up R Environment","text":"","category":"section"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"PhylogeneticFactorAnalysis uses R to generate automated plots via the Julia package RCall. If you are interested in producing these plots, we recommend setting up your R environment with the appropriate packages and ensuring RCall can locate your current R installation.","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"If you do not have R installed already, you can download and install it from the R website.","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Once you have installed R, ensure the appropriate packages are installed by typing the following into the R console:","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"install.packages(c(\"ggplot2\", \"wesanderson\", \"colorspace\", \"tidyr\", \"ggtree\", \"phytools\", \"tidytree\", \"aplot\", \"RColorBrewer\", \"ggnewscale\", \"phyclust\", \"treeio\"))","category":"page"},{"location":"#Add-PhylogeneticFactorAnalysis-to-your-Julia-environment","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"Add PhylogeneticFactorAnalysis to your Julia environment","text":"","category":"section"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"From the Julia REPL, run","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"import Pkg; Pkg.add(\"PhylogeneticFactorAnalysis\")","category":"page"},{"location":"#Check-dependencies","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"Check dependencies","text":"","category":"section"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Once you have added PhylogeneticFactorAnalysis, check that the dependencies are properly installed.","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"To check Java and BEAST, enter into the Julia REPL","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"using PhylogeneticFactorAnalysis\n\ncheck_beast()","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"To check R, add the following line to the code above:","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"check_r()","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"It's OK (and likely inevitable) if check_r() results in several R warnings. As long as it does not throw an error and terminates with all necessary R packages installed your R environment should be set up appropriately.","category":"page"},{"location":"#Instructions","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"Instructions","text":"","category":"section"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"PhylogeneticFactorAnalysis takes a relatively simple xml file and uses it to build and run BEAST xml files. Users must specify the following:","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"a csv file containing the trait data\nnewick-formatted tree file\nprior and constraint on the loadings matrix","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Here is the simplest xml you can supply","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<phylogeneticFactorAnalysis\n    name=\"test\" <!-- this will be the header of all folders and files created -->\n    standardizeTraits=\"true\" <!-- we recommend standardizing traits, but you may set to false if desired-->\n    >\n  <data traits=\"test_traits.csv\" <!-- path to your traits file -->\n        tree=\"test_newick.txt\" <!-- path to your newick file -->\n        />\n  <iidPrior constraint=\"upperTriangular\" <!-- this is the specific prior you use on the loadings. the options are \"upperTriangular\", \"orthogonal\", or \"hybrid\". -->\n        />\n  <modelSelection>\n    <nFactors>5</nFactors> <!-- this is the number of factors you want to use. you may supply as many as you want, (e.g. <nFactors>1 2 3 4 5</nFactors>). If you only supply one it will skip model selection and simply choose that number of factors for the final analysis-->\n  </modelSelection>\n</phylogeneticFactorAnalysis>","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"The above xml performs a single BEAST run with a 5-factor model and upper-triangular constraint on the loadings.","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Here is a more complex xml that can handle discrete traits:","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<phylogeneticFactorAnalysis name=\"aquilegiaBinary\" partitionSeed=\"666\" mcmcSeed=\"666\" overwrite=\"true\" standardizeTraits=\"true\">\n  <data traits=\"aquilegia_binary.csv\" tree=\"aquilegia_newick.txt\" discreteIndices=\"11 12 13 14\"/>\n  <iidPrior constraint=\"orthogonal\"/>\n  <modelSelection repeats=\"5\" selectionStatistic=\"CLPD\" burnin=\"0.25\">\n    <nFactors>1 2 3 4 5</nFactors>\n    <mcmcOptions chainLength=\"100000\"/>\n  </modelSelection>\n  <mcmcOptions chainLength=\"100000\"/>\n</phylogeneticFactorAnalysis>","category":"page"},{"location":"#Explanation-of-elements-and-attributes","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"Explanation of elements and attributes","text":"","category":"section"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Element: phylogeneticFactorAnalysis","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"This is the root element that encompasses all other elements.\nAttributes\nname: the name of the folder and files that will be created. This can be set to anything\npartitionSeed (optional): the random number seed used for creating the training and test data sets. Can be any positive integer.\nmcmcSeed (optional): the random number seed for the BEAST MCMC simulation\noverwrite (optional): allows you to overwrite any existing folders. Possible values are \"true\" or \"false\". If this is set to \"true\", it will automatically delete everything in the folder specified by name.\nstandardizeTraits (optional): specifies whether data should be standardized (i.e. centered and scaled) prior to analysis. Possible values are \"true\" or \"false\".","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Element: data","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Locations and instructions related to the input data.\nAttributes\ntraits: path to csv file where your trait data is stored. The csv file must have the first column corresponding to the taxa with column header \"taxon\".\ntree: path to the newick-formatted tree file.\ndiscreteIndices (optional): if any traits are discrete, list the indices that correspond to the discrete traits.","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Element: iidPrior","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"This is the prior. You must specify one prior, but you may specify orthogonalShrinkagePrior  (below) instead of iidPrior.\nAttributes\nconstraint: the constraint on the loadings matrix. Options are \"upperTriangular\", \"orthogonal\", or \"hybrid\".","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Element: orthogonalShrinkagePrior","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"This is the prior. You must specify one prior, but you may specify iidPrior (above) instead of orthogonalShrinkagePrior.\nAttributes: TODO","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Element: modelSelection","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Supplies the models you want to compare to each other.\nAttributes\nrepeats: Specifies the \"k\" value for k-fold cross validation. Should be a positive integer > 1 (although we wouldn't recommend going below 5).\nselectionStatistic (optional): the specific selection statistic you would like to use. Values are \"CLPD\" or \"MSE\".\nburnin (optional): the MCMC burnin before evaluating the selection statistic. Should be any number between 0 and 1. Defaults to 0.25.\nSub-element: nFactors\nspecify the number of factors for each model\nSub-element: mcmcOptions (optional)\nMeta-parameters related to the MCMC simulation\nAttributes\nchainLength: the number of states to run the MCMC chain. Defaults to 10000.","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Element: mcmcOptions","category":"page"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"Meta-parameters related to the final MCMC simulation","category":"page"},{"location":"#Running-analyses","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"Running analyses","text":"","category":"section"},{"location":"","page":"PhylogeneticFactorAnalysis.jl Documentation","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"The only function required to run analyses is pfa. If test.xml is the path to your xml file, simply enter pfa(\"test.xml\") into the Julia REPL to run the pipeline.","category":"page"}]
}

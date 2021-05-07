var documenterSearchIndex = {"docs":
[{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Example data and XML files can be examined and run through the import_example function. The examples are numbered 1-4. To access the first set of example files, run:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"import_example(1)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"This will create a directory pfa_example_1 in your working directory and populate it with the relevant files. The files associated with Example 1 are aquilegia_continuous.csv containing the trait data, aquilegia_newick.txt containing the phylogenetic tree, and example_1.xml containing the instructions for running PFA. To run this example, enter","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"pfa(\"pfa_example_1/example_1.xml\")","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The table below provides a quick summary of the instructions provided in each example XML file.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Example Data Type Prior Model Selection\n1 continuous i.i.d. Gaussian yes\n2 continuous & discrete i.i.d. Gaussian yes\n3 continuous shrinkage yes\n4 continuous i.i.d. Gaussian no","category":"page"},{"location":"#PhylogeneticFactorAnalysis.jl-Documentation","page":"Home","title":"PhylogeneticFactorAnalysis.jl Documentation","text":"","category":"section"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#Install-Java","page":"Home","title":"Install Java","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PhylogeneticFactorAnalysis will not automatically install Java, which is required to run BEAST. To check if you have Java installed, run java -version from your computer's command line (not in Julia). If Java is installed, you will see something like:","category":"page"},{"location":"","page":"Home","title":"Home","text":"$ java -version\nopenjdk version \"1.8.0_282\"\nOpenJDK Runtime Environment (build 1.8.0_282-8u282-b08-0ubuntu1~18.04-b08)\nOpenJDK 64-Bit Server VM (build 25.282-b08, mixed mode)","category":"page"},{"location":"","page":"Home","title":"Home","text":"If Java is not installed, follow the instructions on the BEAST website for installin Java. Note that you ONLY need to install Java. You do NOT need to install BEAST, which comes with PhylogeneticFactorAnalysis.","category":"page"},{"location":"#Optional:-Set-up-R-Environment","page":"Home","title":"Optional: Set up R Environment","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PhylogeneticFactorAnalysis uses R to generate automated plots via the Julia package RCall. If you are interested in producing these plots, we recommend setting up your R environment with the appropriate packages and ensuring RCall can locate your current R installation.","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you do not have R installed already, you can download and install it from the R website.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Once you have installed R, ensure the appropriate packages are installed by typing the following into the R console:","category":"page"},{"location":"","page":"Home","title":"Home","text":"install.packages(c(\"ggplot2\", \"wesanderson\", \"colorspace\", \"tidyr\", \"ggtree\", \"phytools\", \"tidytree\", \"aplot\", \"RColorBrewer\", \"ggnewscale\", \"phyclust\", \"treeio\"))","category":"page"},{"location":"#Add-PhylogeneticFactorAnalysis-to-your-Julia-environment","page":"Home","title":"Add PhylogeneticFactorAnalysis to your Julia environment","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"From the Julia REPL, run","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg; Pkg.add(\"PhylogeneticFactorAnalysis\")","category":"page"},{"location":"#Check-dependencies","page":"Home","title":"Check dependencies","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Once you have added PhylogeneticFactorAnalysis, check that the dependencies are properly installed.","category":"page"},{"location":"","page":"Home","title":"Home","text":"To check Java and BEAST, enter into the Julia REPL","category":"page"},{"location":"","page":"Home","title":"Home","text":"using PhylogeneticFactorAnalysis\n\ncheck_beast()","category":"page"},{"location":"","page":"Home","title":"Home","text":"To check R, add the following line to the code above:","category":"page"},{"location":"","page":"Home","title":"Home","text":"check_r()","category":"page"},{"location":"","page":"Home","title":"Home","text":"It's OK (and likely inevitable) if check_r() results in several R warnings. As long as it does not throw an error and terminates with all necessary R packages installed your R environment should be set up appropriately.","category":"page"},{"location":"#Instructions","page":"Home","title":"Instructions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PhylogeneticFactorAnalysis takes a relatively simple xml file and uses it to build and run BEAST xml files. Users must specify the following:","category":"page"},{"location":"","page":"Home","title":"Home","text":"a csv file containing the trait data\nnewick-formatted tree file\nprior and constraint on the loadings matrix","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here is the simplest xml you can supply","category":"page"},{"location":"","page":"Home","title":"Home","text":"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<phylogeneticFactorAnalysis\n    name=\"test\" <!-- this will be the header of all folders and files created -->\n    standardizeTraits=\"true\" <!-- we recommend standardizing traits, but you may set to false if desired-->\n    >\n  <data traits=\"test_traits.csv\" <!-- path to your traits file -->\n        tree=\"test_newick.txt\" <!-- path to your newick file -->\n        />\n  <iidPrior constraint=\"upperTriangular\" <!-- this is the specific prior you use on the loadings. the options are \"upperTriangular\", \"orthogonal\", or \"hybrid\". -->\n        />\n  <modelSelection>\n    <nFactors>5</nFactors> <!-- this is the number of factors you want to use. you may supply as many as you want, (e.g. <nFactors>1 2 3 4 5</nFactors>). If you only supply one it will skip model selection and simply choose that number of factors for the final analysis-->\n  </modelSelection>\n</phylogeneticFactorAnalysis>","category":"page"},{"location":"","page":"Home","title":"Home","text":"The above xml performs a single BEAST run with a 5-factor model and upper-triangular constraint on the loadings.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here is a more complex xml that can handle discrete traits:","category":"page"},{"location":"","page":"Home","title":"Home","text":"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<phylogeneticFactorAnalysis name=\"aquilegiaBinary\" partitionSeed=\"666\" mcmcSeed=\"666\" overwrite=\"true\" standardizeTraits=\"true\">\n  <data traits=\"aquilegia_binary.csv\" tree=\"aquilegia_newick.txt\" discreteIndices=\"11 12 13 14\"/>\n  <iidPrior constraint=\"orthogonal\"/>\n  <modelSelection repeats=\"5\" selectionStatistic=\"CLPD\" burnin=\"0.25\">\n    <nFactors>1 2 3 4 5</nFactors>\n    <mcmcOptions chainLength=\"100000\"/>\n  </modelSelection>\n  <mcmcOptions chainLength=\"100000\"/>\n</phylogeneticFactorAnalysis>","category":"page"},{"location":"#Explanation-of-elements-and-attributes","page":"Home","title":"Explanation of elements and attributes","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Element: phylogeneticFactorAnalysis","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is the root element that encompasses all other elements.\nAttributes\nname: the name of the folder and files that will be created. This can be set to anything\npartitionSeed (optional): the random number seed used for creating the training and test data sets. Can be any positive integer.\nmcmcSeed (optional): the random number seed for the BEAST MCMC simulation\noverwrite (optional): allows you to overwrite any existing folders. Possible values are \"true\" or \"false\". If this is set to \"true\", it will automatically delete everything in the folder specified by name.\nstandardizeTraits (optional): specifies whether data should be standardized (i.e. centered and scaled) prior to analysis. Possible values are \"true\" or \"false\".","category":"page"},{"location":"","page":"Home","title":"Home","text":"Element: data","category":"page"},{"location":"","page":"Home","title":"Home","text":"Locations and instructions related to the input data.\nAttributes\ntraits: path to csv file where your trait data is stored. The csv file must have the first column corresponding to the taxa with column header \"taxon\".\ntree: path to the newick-formatted tree file.\ndiscreteIndices (optional): if any traits are discrete, list the indices that correspond to the discrete traits.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Element: iidPrior","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is the prior. You must specify one prior, but you may specify orthogonalShrinkagePrior  (below) instead of iidPrior.\nAttributes\nconstraint: the constraint on the loadings matrix. Options are \"upperTriangular\", \"orthogonal\", or \"hybrid\".","category":"page"},{"location":"","page":"Home","title":"Home","text":"Element: orthogonalShrinkagePrior","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is the prior. You must specify one prior, but you may specify iidPrior (above) instead of orthogonalShrinkagePrior.\nAttributes: TODO","category":"page"},{"location":"","page":"Home","title":"Home","text":"Element: modelSelection","category":"page"},{"location":"","page":"Home","title":"Home","text":"Supplies the models you want to compare to each other.\nAttributes\nrepeats: Specifies the \"k\" value for k-fold cross validation. Should be a positive integer > 1 (although we wouldn't recommend going below 5).\nselectionStatistic (optional): the specific selection statistic you would like to use. Values are \"CLPD\" or \"MSE\".\nburnin (optional): the MCMC burnin before evaluating the selection statistic. Should be any number between 0 and 1. Defaults to 0.25.\nSub-element: nFactors\nspecify the number of factors for each model\nSub-element: mcmcOptions (optional)\nMeta-parameters related to the MCMC simulation\nAttributes\nchainLength: the number of states to run the MCMC chain. Defaults to 10000.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Element: mcmcOptions","category":"page"},{"location":"","page":"Home","title":"Home","text":"Meta-parameters related to the final MCMC simulation","category":"page"},{"location":"#Running-analyses","page":"Home","title":"Running analyses","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To load the package, type using PhylogeneticFactorAnalysis into the Julia REPL.","category":"page"},{"location":"","page":"Home","title":"Home","text":"After loading the package, the only function required to run analyses is pfa. If test.xml is the path to your xml file, simply enter pfa(\"test.xml\") into the Julia REPL to run the pipeline.","category":"page"}]
}
using Documenter
using PhylogeneticFactorAnalysis

makedocs(
    sitename = "PhylogeneticFactorAnalysis",
    format = Documenter.HTML(),
    modules = [PhylogeneticFactorAnalysis],
    pages = ["Home" => "index.md",
             "Examples" => "examples.md"]
)

deploydocs(
    repo = "github.com/gabehassler/PhylogeneticFactorAnalysis.jl.git",
    devbranch = "main"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

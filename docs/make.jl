using Documenter
using PhylogeneticFactorAnalysis

makedocs(
    sitename = "PhylogeneticFactorAnalysis",
    format = Documenter.HTML(),
    modules = [PhylogeneticFactorAnalysis]
)

deploydocs(
    repo = "https://github.com/gabehassler/PhylogeneticFactorAnalysis.jl.git",
    devbranch = "main"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

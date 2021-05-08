repos <- c("remotes",
                   "ggplot2",
                   "wesanderson",
                   "colorspace",
                   "tidyr",
                   "ggtree",
                   "phytools",
                   "tidytree",
                   "aplot",
                   "RColorBrewer",
                   "ggnewscale",
                   "phyclust",
                   "treeio",
                   "BiocManager")
bio_repos <- c("ggtree", "treeio")

for (repo in repos) {
   if (!require(repo)) install.packages(repo, repos="https://cloud.r-project.org/")
}

for (repo in bio_repos) {
   BiocManager::install(repo)
}

# saveRDS(remotes::dev_package_deps(dependencies = TRUE), ".github/depends.Rds", version = 2)

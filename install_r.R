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
                   "treeio")

for (repo in repos) {
   print("======================================================================")
   print(paste("installing", repo, "..."))
   install.packages(repo, repos="https://cloud.r-project.org/")
   print(paste(repo, "installed"))
}

# saveRDS(remotes::dev_package_deps(dependencies = TRUE), ".github/depends.Rds", version = 2)

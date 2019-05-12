pkgs = c("BiocManager",
         "tidyverse", "knitr", "rmarkdown",
         "kableExtra", "cowplot", "RColorBrewer",
         "pheatmap", "viridis" )
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)

BiocManager::install("edgeR")
BiocManager::install("Glimma")
BiocManager::install("limma")
BiocManager::install("DESeq2")
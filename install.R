pkgs = c("tidyverse", "knitr", "rmarkdown",
         "kableExtra", "cowplot", "RColorBrewer",
         "pheatmap", "viridis" )
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("edgeR")
BiocManager::install("Glimma")
BiocManager::install("limma")
BiocManager::install("DESeq2")
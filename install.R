pkgs = c("tidyverse", "knitr", "rmarkdown",
         "kableExtra", "cowplot", "RColorBrewer",
         "pheatmap", "viridis" )
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("Glimma")
biocLite("limma")
biocLite("DESeq2")
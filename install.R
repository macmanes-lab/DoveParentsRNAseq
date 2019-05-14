pkgs = c("tidyverse", "knitr", "rmarkdown",
         "kableExtra", "cowplot", "RColorBrewer",
         "pheatmap", "viridis" )
install.packages(pkgs)

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("Glimma")
biocLite("limma")
biocLite("DESeq2")
pkgs = c("tidyverse", "knitr", "rmarkdown",
         "kableExtra", "cowplot", "RColorBrewer",
         "pheatmap", "viridis",
         "ggfortify", "factoextra",
         "corrplot",
         "readxl", "modelr", "lubridate",
         "WGCNA", "magrittr",  "forcats",
         "ggimage")
         
install.packages(pkgs)

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
biocLite("Glimma")
biocLite("limma")
biocLite("DESeq2")
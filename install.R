pkgs = c("tidyverse", "knitr", "rmarkdown",
          "cowplot", "viridis", "readr", "ggsignif")
         
install.packages(pkgs)

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("BiocParallel")
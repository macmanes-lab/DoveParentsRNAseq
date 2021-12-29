pkgs = c("tidyverse", "knitr", "rmarkdown",
          "cowplot", "viridis", "readr", "ggsignif",
        "Rtsne")
         
install.packages(pkgs)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("BiocParallel")
BiocManager::install("recount3")
BiocManager::install("biomaRt")

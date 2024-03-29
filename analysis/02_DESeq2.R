# 02_DESeq2.R"

library(tidyverse)
library(DESeq2)
library(BiocParallel)
register(MulticoreParam(6))

source("R/functions.R")  # load custom functions 
source("R/themes.R")  # load custom themes and color palletes
source("R/genelists.R") # load candidate genes

# DESeq2 was not designed to run on 100 + samples. But, I really like it, so I do it anyways. 
# Some of these commands take like 15 min to run using 6 cores. 

# count data
countData <- read_csv("results/00_countData.csv") %>%
  column_to_rownames(var = "X1")

#### uncomment this to subset the data for quick analysis
#print("subset for quick run")
#countData  <- head(countData, 1550)

# col data or variable informaiton
colData <- read.csv("metadata/00_colData.csv", header = T, row.names = 1) %>%
  mutate(tissue = factor(tissue, levels = tissuelevels),
         treatment <- factor(treatment, levels = alllevels)) %>%
  mutate(sextissue = as.factor(paste(sex, tissue, sep = "_"))) 
colData <- colData %>% drop_na()
head(colData)
dim(colData)
# 984 samples, 11 variables

print(ncol(countData) == nrow(colData))

head(names(countData))
head(row.names(colData))

savealltheDEGs <- function(whichdds, whichgroup){
  
  createDEGdfs(whichdds, whichgroup, "bldg", "control")
  createDEGdfs(whichdds, whichgroup, references, charlevelsnocontrol)
  
  # sequential
  createDEGdfs(whichdds, whichgroup, "lay", "inc.d3")
  createDEGdfs(whichdds, whichgroup, "inc.d3", "inc.d9")
  createDEGdfs(whichdds, whichgroup, "inc.d9", "inc.d17")
  createDEGdfs(whichdds, whichgroup, "inc.d17", "hatch")
  createDEGdfs(whichdds, whichgroup, "hatch", "n5")
  createDEGdfs(whichdds, whichgroup, "n5", "n9")
  
  # removal
  createDEGdfs(whichdds, whichgroup, "inc.d3", "m.inc.d3")
  createDEGdfs(whichdds, whichgroup, "inc.d9", "m.inc.d9")
  createDEGdfs(whichdds, whichgroup, "inc.d17", "m.inc.d17")
  createDEGdfs(whichdds, whichgroup, "hatch", "m.n2")
  
  # replacement
  createDEGdfs(whichdds, whichgroup, manipcontrols, replacements)
}

######## differential gene expression


ddsHf <- returndds2(c("female_hypothalamus"))
savealltheDEGs(ddsHf, "female_hypothalamus")

ddsHm <- returndds2(c("male_hypothalamus"))
savealltheDEGs(ddsHm, "male_hypothalamus")


ddsGm <- returndds2(c("male_gonads"))
savealltheDEGs(ddsGm, "male_gonads")

ddsGf <- returndds2(c("female_gonads"))
savealltheDEGs(ddsGf, "female_gonads")

ddsPf <- returndds2(c("female_pituitary"))
savealltheDEGs(ddsPf, "female_pituitary")

ddsPm <- returndds2(c("male_pituitary"))
savealltheDEGs(ddsPm, "male_pituitary")




# 02_DESeq2"

library(tidyverse)
library(DESeq2)
library(BiocParallel)
register(MulticoreParam(6))
library(cowplot)

source("R/functions.R")  # load custom functions 
source("R/themes.R")  # load custom themes and color palletes

# DESeq2 was not designed to run on 100 + samples. But, I really like it, so I do it anyways. Some of these commands take like 15 min to run using 6 cores. 
# Also, I don't run every chuck every time. When I want new analysese, I add them to the bottom and set `eval = F` for chunks I don't want to rerurn 

# count data
countData <- read.csv("results/00_counts.csv", header = T, row.names = 1)

#### uncomment this to subset the data for quick analysis
#print("subset for quick run")
#countData  <- head(countData, 1500)

# col data or variable informaiton
colData <- read.csv("metadata/00_colData.csv", header = T, row.names = 1) %>%
  mutate(tissue = factor(tissue, levels = tissuelevels),
         treatment <- factor(treatment, levels = alllevels)) %>%
  mutate(sextissue = as.factor(paste(sex, tissue, sep = "_"))) 
colData <- colData %>% drop_na()
# 987 samples, 11 variables

print(ncol(countData) == nrow(colData))

######## differential gene expression

# contrastst part 1. versus control and bldg
references <- c("control", "bldg")
allotherstages   <- c(charlevelsnocontrol, maniplevels)

# contrastst part 2. manips versus internal or other manips
manipcontrols <- c("inc.d9", "inc.d17", "hatch", "n5",  "prolong", "early")
replacements   <- c("early", "prolong", "extend")

ddsMG <- returndds("male_hypothalamus")
vsdMG <- returnvsd(ddsMG, "male_hypothalamus")
createDEGdftreatmentvcontrols(ddsMG, "male_hypothalamus", references, allotherstages)
createDEGdftreatmentvcontrols(ddsMG, "male_hypothalamus", manipcontrols, replacements)

ddsFH <- returndds("female_hypothalamus")
vsdFH <- returnvsd(ddsFG, "female_hypothalamus")
createDEGdftreatmentvcontrols(ddsFH, "female_hypothalamus", references, allotherstages)
createDEGdftreatmentvcontrols(ddsFH, "female_hypothalamus", manipcontrols, replacements)

ddsMP <- returndds("male_pituitary")
vsdMP <- returnvsd(ddsMP, "male_pituitary")
createDEGdftreatmentvcontrols(ddsMG, "male_pituitary", references, allotherstages)
createDEGdftreatmentvcontrols(ddsMG, "male_pituitary", manipcontrols, replacements)

ddsFP <- returndds("female_pituitary")
vsdFP <- returnvsd(ddsFP, "female_pituitary")
createDEGdftreatmentvcontrols(ddsMP, "male_hypothalamus", references, allotherstages)
createDEGdftreatmentvcontrols(ddsMP, "male_hypothalamus", manipcontrols, replacements)

ddsMG <- returndds("male_gonads")
vsdMG <- returnvsd(ddsMG, "male_gonads")
createDEGdftreatmentvcontrols(ddsMG, "male_gonads", references, allotherstages)
createDEGdftreatmentvcontrols(ddsMG, "male_gonads", manipcontrols, replacements)

ddsFG <- returndds("female_gonads")
vsdFG <- returndds(ddsFG, "female_gonads")
createDEGdftreatmentvcontrols(ddsFG, "female_gonads", references, allotherstages)
createDEGdftreatmentvcontrols(ddsFG, "female_gonads", manipcontrols, replacements)

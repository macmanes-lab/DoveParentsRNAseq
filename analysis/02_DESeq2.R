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

#ddsMP <- returndds("male_pituitary")
#vsdMP <- returnvsd(ddsMP, "male_pituitary")
#createDEGdftreatmentvcontrols(ddsMP, "male_pituitary", references, allotherstages)
#createDEGdftreatmentvcontrols(ddsMP, "male_pituitary", manipcontrols, replacements)
#createDEGdftreatmentvcontrols(ddsMP, "male_pituitary", "lay", "inc.d3")
#createDEGdftreatmentvcontrols(ddsMP, "male_pituitary", "inc.d3", "inc.d9")
#createDEGdftreatmentvcontrols(ddsMP, "male_pituitary", "inc.d9", "inc.d17")
#createDEGdftreatmentvcontrols(ddsMP, "male_pituitary", "inc.d17", "hatch")
#createDEGdftreatmentvcontrols(ddsMP, "male_pituitary", "hatch", "n5")
#createDEGdftreatmentvcontrols(ddsMP, "male_pituitary", "n5", "n9")

#ddsFP <- returndds("female_pituitary")
#vsdFP <- returnvsd(ddsFP, "female_pituitary")
#createDEGdftreatmentvcontrols(ddsFP, "female_pituitary", references, allotherstages)
#createDEGdftreatmentvcontrols(ddsFP, "female_pituitary", manipcontrols, replacements)
#createDEGdftreatmentvcontrols(ddsFP, "female_pituitary", "lay", "inc.d3")
#createDEGdftreatmentvcontrols(ddsFP, "female_pituitary", "inc.d3", "inc.d9")
#createDEGdftreatmentvcontrols(ddsFP, "female_pituitary", "inc.d9", "inc.d17")
#createDEGdftreatmentvcontrols(ddsFP, "female_pituitary", "inc.d17", "hatch")
#createDEGdftreatmentvcontrols(ddsFP, "female_pituitary", "hatch", "n5")
#createDEGdftreatmentvcontrols(ddsFP, "female_pituitary", "n5", "n9")

#ddsMG <- returndds("male_gonads")
#vsdMG <- returnvsd(ddsMG, "male_gonads")
#createDEGdftreatmentvcontrols(ddsMG, "male_gonads", references, allotherstages)
#createDEGdftreatmentvcontrols(ddsMG, "male_gonads", manipcontrols, replacements)
#createDEGdftreatmentvcontrols(ddsMG, "male_gonads", "lay", "inc.d3")
#createDEGdftreatmentvcontrols(ddsMG, "male_gonads", "inc.d3", "inc.d9")
#createDEGdftreatmentvcontrols(ddsMG, "male_gonads", "inc.d9", "inc.d17")
#createDEGdftreatmentvcontrols(ddsMG, "male_gonads", "inc.d17", "hatch")
#createDEGdftreatmentvcontrols(ddsMG, "male_gonads", "hatch", "n5")
#createDEGdftreatmentvcontrols(ddsMG, "male_gonads", "n5", "n9")

ddsFG <- returndds("female_gonads")
vsdFG <- returnvsd(ddsFG, "female_gonads")
#createDEGdftreatmentvcontrols(ddsFG, "female_gonads", references, allotherstages)
#createDEGdftreatmentvcontrols(ddsFG, "female_gonads", manipcontrols, replacements)
#createDEGdftreatmentvcontrols(ddsFG, "female_gonads", "lay", "inc.d3")
#createDEGdftreatmentvcontrols(ddsFG, "female_gonads", "inc.d3", "inc.d9")
#createDEGdftreatmentvcontrols(ddsFG, "female_gonads", "inc.d9", "inc.d17")
#createDEGdftreatmentvcontrols(ddsFG, "female_gonads", "inc.d17", "hatch")
#createDEGdftreatmentvcontrols(ddsFG, "female_gonads", "hatch", "n5")
#createDEGdftreatmentvcontrols(ddsFG, "female_gonads", "n5", "n9")

#ddsMH <- returndds("male_hypothalamus")
#vsdMH <- returnvsd(ddsMH, "male_hypothalamus")
#createDEGdftreatmentvcontrols(ddsMH, "male_hypothalamus", references, allotherstages)
#createDEGdftreatmentvcontrols(ddsMH, "male_hypothalamus", manipcontrols, replacements)
#createDEGdftreatmentvcontrols(ddsMH, "male_hypothalamus", "lay", "inc.d3")
#createDEGdftreatmentvcontrols(ddsMH, "male_hypothalamus", "inc.d3", "inc.d9")
#createDEGdftreatmentvcontrols(ddsMH, "male_hypothalamus", "inc.d9", "inc.d17")
#createDEGdftreatmentvcontrols(ddsMH, "male_hypothalamus", "inc.d17", "hatch")
#createDEGdftreatmentvcontrols(ddsMH, "male_hypothalamus", "hatch", "n5")
#createDEGdftreatmentvcontrols(ddsMH, "male_hypothalamus", "n5", "n9")

ddsFH <- returndds("female_hypothalamus")
vsdFH <- returnvsd(ddsFH, "female_hypothalamus")
#createDEGdftreatmentvcontrols(ddsFH, "female_hypothalamus", references, allotherstages)
#createDEGdftreatmentvcontrols(ddsFH, "female_hypothalamus", manipcontrols, replacements)
#createDEGdftreatmentvcontrols(ddsFH, "female_hypothalamus", "lay", "inc.d3")
#createDEGdftreatmentvcontrols(ddsFH, "female_hypothalamus", "inc.d3", "inc.d9")
#createDEGdftreatmentvcontrols(ddsFH, "female_hypothalamus", "inc.d9", "inc.d17")
#createDEGdftreatmentvcontrols(ddsFH, "female_hypothalamus", "inc.d17", "hatch")
#createDEGdftreatmentvcontrols(ddsFH, "female_hypothalamus", "hatch", "n5")
#createDEGdftreatmentvcontrols(ddsFH, "female_hypothalamus", "n5", "n9")


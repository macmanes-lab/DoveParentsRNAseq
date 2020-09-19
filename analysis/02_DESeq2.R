# 02_DESeq2"

library(tidyverse)
library(DESeq2)
library(BiocParallel)
register(MulticoreParam(6))

source("R/functions.R")  # load custom functions 
source("R/themes.R")  # load custom themes and color palletes
source("R/genelists.R") # load candidate genes

# DESeq2 was not designed to run on 100 + samples. But, I really like it, so I do it anyways. 
# Some of these commands take like 15 min to run using 6 cores. 
# Also, I don't run every chuck every time. When I want new analyses, 
# I add them to the bottom and set `eval = F` for chunks I don't want to rerurn 

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


### same code, for female and male HPG

#ddsMH <- returndds("male_hypothalamus")
#vsdMH <- returnvsd(ddsMH, "male_hypothalamus")
#createDEGdfs(ddsMH, "male_hypothalamus", references, allotherstages)
#createDEGdfs(ddsMH, "male_hypothalamus", manipcontrols, replacements)

# sequential
#createDEGdfs(ddsMH, "male_hypothalamus", "lay", "inc.d3")
#createDEGdfs(ddsMH, "male_hypothalamus", "inc.d3", "inc.d9")
#createDEGdfs(ddsMH, "male_hypothalamus", "inc.d9", "inc.d17")
#createDEGdfs(ddsMH, "male_hypothalamus", "inc.d17", "hatch")
#createDEGdfs(ddsMH, "male_hypothalamus", "hatch", "n5")
#createDEGdfs(ddsMH, "male_hypothalamus", "n5", "n9")
#createDEGdfs(ddsMH, "male_hypothalamus", "control", "bldg")

# removal
#createDEGdfs(ddsMH, "male_hypothalamus",  "inc.d3", "m.inc.d3")
#createDEGdfs(ddsMH, "male_hypothalamus",  "inc.d9", "m.inc.d9")
#createDEGdfs(ddsMH, "male_hypothalamus", "inc.d17", "m.inc.d17")
#createDEGdfs(ddsMH, "male_hypothalamus", "hatch", "m.n2")


#ddsFH <- returndds("female_hypothalamus")
#vsdFH <- returnvsd(ddsFH, "female_hypothalamus")
#createDEGdfs(ddsFH, "female_hypothalamus", references, allotherstages)
#createDEGdfs(ddsFH, "female_hypothalamus", manipcontrols, replacements)
#createDEGdfs(ddsFH, "female_hypothalamus", "lay", "inc.d3")
#createDEGdfs(ddsFH, "female_hypothalamus", "inc.d3", "inc.d9")
#createDEGdfs(ddsFH, "female_hypothalamus", "inc.d9", "inc.d17")
#createDEGdfs(ddsFH, "female_hypothalamus", "inc.d17", "hatch")
#createDEGdfs(ddsFH, "female_hypothalamus", "hatch", "n5")
#createDEGdfs(ddsFH, "female_hypothalamus", "n5", "n9")
#createDEGdfs(ddsFH, "female_hypothalamus", "control", "bldg")
#createDEGdfs(ddsFH, "female_hypothalamus",  "inc.d3", "m.inc.d3")
#createDEGdfs(ddsFH, "female_hypothalamus", "inc.d9", "m.inc.d9")
#createDEGdfs(ddsFH, "female_hypothalamus", "inc.d17", "m.inc.d17")
#createDEGdfs(ddsFH, "female_hypothalamus", "hatch", "m.n2")


ddsMP <- returndds("male_pituitary")
vsdMP <- returnvsd(ddsMP, "male_pituitary")
#createDEGdfs(ddsMP, "male_pituitary", references, allotherstages)
createDEGdfs(ddsMP, "male_pituitary", manipcontrols, replacements)
createDEGdfs(ddsMP, "male_pituitary", "lay", "inc.d3")
createDEGdfs(ddsMP, "male_pituitary", "inc.d3", "inc.d9")
createDEGdfs(ddsMP, "male_pituitary", "inc.d9", "inc.d17")
createDEGdfs(ddsMP, "male_pituitary", "inc.d17", "hatch")
createDEGdfs(ddsMP, "male_pituitary", "hatch", "n5")
createDEGdfs(ddsMP, "male_pituitary", "n5", "n9")
createDEGdfs(ddsMP, "male_pituitary", "control", "bldg")
createDEGdfs(ddsMP, "male_pituitary", "inc.d3", "m.inc.d3")
createDEGdfs(ddsMP, "male_pituitary", "inc.d9", "m.inc.d9")
createDEGdfs(ddsMP, "male_pituitary", "inc.d17", "m.inc.d17")
createDEGdfs(ddsMP, "male_pituitary", "hatch", "m.n2")


ddsFP <- returndds("female_pituitary")
vsdFP <- returnvsd(ddsFP, "female_pituitary")
createDEGdfs(ddsFP, "female_pituitary", references, allotherstages)
createDEGdfs(ddsFP, "female_pituitary", manipcontrols, replacements)
createDEGdfs(ddsFP, "female_pituitary", "lay", "inc.d3")
createDEGdfs(ddsFP, "female_pituitary", "inc.d3", "inc.d9")
createDEGdfs(ddsFP, "female_pituitary", "inc.d9", "inc.d17")
createDEGdfs(ddsFP, "female_pituitary", "inc.d17", "hatch")
createDEGdfs(ddsFP, "female_pituitary", "hatch", "n5")
createDEGdfs(ddsFP, "female_pituitary", "n5", "n9")
createDEGdfs(ddsFP, "female_pituitary", "control", "bldg")
createDEGdfs(ddsFP, "female_pituitary", "inc.d3", "m.inc.d3")
createDEGdfs(ddsFP, "female_pituitary", "inc.d9", "m.inc.d9")
createDEGdfs(ddsFP, "female_pituitary", "inc.d17", "m.inc.d17")
createDEGdfs(ddsFP, "female_pituitary", "hatch", "m.n2")


ddsMG <- returndds("male_gonads")
vsdMG <- returnvsd(ddsMG, "male_gonads")
createDEGdfs(ddsMG, "male_gonads", references, allotherstages)
createDEGdfs(ddsMG, "male_gonads", manipcontrols, replacements)
createDEGdfs(ddsMG, "male_gonads", "lay", "inc.d3")
createDEGdfs(ddsMG, "male_gonads", "inc.d3", "inc.d9")
createDEGdfs(ddsMG, "male_gonads", "inc.d9", "inc.d17")
createDEGdfs(ddsMG, "male_gonads", "inc.d17", "hatch")
createDEGdfs(ddsMG, "male_gonads", "hatch", "n5")
createDEGdfs(ddsMG, "male_gonads", "n5", "n9")
createDEGdfs(ddsMG, "male_gonads", "control", "bldg")
createDEGdfs(ddsMG, "male_gonads", "inc.d3", "m.inc.d3")
createDEGdfs(ddsMG, "male_gonads", "inc.d9", "m.inc.d9")
createDEGdfs(ddsMG, "male_gonads", "inc.d17", "m.inc.d17")
createDEGdfs(ddsMG, "male_gonads", "hatch", "m.n2")


ddsFG <- returndds("female_gonads")
vsdFG <- returnvsd(ddsFG, "female_gonads")
createDEGdfs(ddsFG, "female_gonads", references, allotherstages)
createDEGdfs(ddsFG, "female_gonads", manipcontrols, replacements)
createDEGdfs(ddsFG, "female_gonads", "lay", "inc.d3")
createDEGdfs(ddsFG, "female_gonads", "inc.d3", "inc.d9")
createDEGdfs(ddsFG, "female_gonads", "inc.d9", "inc.d17")
createDEGdfs(ddsFG, "female_gonads", "inc.d17", "hatch")
createDEGdfs(ddsFG, "female_gonads", "hatch", "n5")
createDEGdfs(ddsFG, "female_gonads", "n5", "n9")
createDEGdfs(ddsFG, "female_gonads", "control", "bldg")
createDEGdfs(ddsFG, "female_gonads",  "inc.d3", "m.inc.d3")
createDEGdfs(ddsFG, "female_gonads",  "inc.d9", "m.inc.d9")
createDEGdfs(ddsFG, "female_gonads",  "inc.d17", "m.inc.d17")
createDEGdfs(ddsFG, "female_gonads", "hatch", "m.n2")
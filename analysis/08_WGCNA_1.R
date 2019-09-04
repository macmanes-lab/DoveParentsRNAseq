

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Load the packages
library(tidyverse)
library(WGCNA)
library(dplyr)
library(magrittr)
library(forcats)


# read the count data
countData <- read.csv("../results/00_countData_characterization.csv", row.names = 1)

# read the sample meta data or column data
colData <- read.csv("../metadata/00_colData_characterization.csv", row.names = 1, stringsAsFactors = T)


## rename rownames and colnames for better vizualiation
colData <- colData %>%
  mutate(sex = fct_recode(sex,
                          "F" = "female",
                          "M" = "male"),
         tissue = fct_recode(tissue,
                             "pit" = "pituitary",
                             "hyp" = "hypothalamus",
                             "gon" = "gonad"),
         hypothesis = fct_recode(treatment,
                                 "anticipation" = "control",
                                 "anticipation" = "bldg",
                                 "incubation" = "lay",
                                 "incubation" = "inc.d3",
                                 "incubation" = "inc.d9",
                                 "incubation" = "inc.d17",
                                 "hatchling.care" = "hatch",
                                 "hatchling.care" = "n5",
                                 "hatchling.care" = "n9"))

# wgcna needs numeric identifiers
colData$ID <- as.numeric(colData$bird)

# create a short sample desctiption, make it the row and columns names of colData and countData
colData$sample <- paste(colData$treatment, colData$sex, colData$tissue,  colData$ID, sep = ".")
  
# set row and count names to be the same
row.names(colData) <- colData$sample
colnames(countData) <- colData$sample

head(colData)
head(countData)

# create new grouping for subsets
colData$sextissue <-  as.factor(paste(colData$sex, colData$tissue, sep = "."))
levels(colData$sextissue)

# subset to samples
colData <- colData %>% dplyr::filter(sextissue %in% c("F.pit", "M.pit")) 

# which counts to save
savecols <- as.character(colData$sample) 
savecols <- as.vector(savecols) 
  
# read counts, save counts that match colData
countData <- countData %>% dplyr::select(one_of(savecols)) 

  
#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


datExpr0 <- as.data.frame(t(countData))

head(names(datExpr0))  # columns are genes
head(rownames(datExpr0)) # rows are samples



#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


gsg = goodSamplesGenes(datExpr0, verbose = 0);
gsg$allOK



#=====================================================================================
#
#  Code chunk 4 removing genes
#
#=====================================================================================


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  #if (sum(!gsg$goodGenes)>0) 
  #  printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  #if (sum(!gsg$goodSamples)>0) 
  #  printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


sampleTree = hclust(dist(datExpr0), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.

plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="", cex.main = 1, cex = 0.4)
abline(h = 2000000, col = "red");



#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# Plot a line to show the cut

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 2000000, minSize = 3)
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================

traitData <- colData %>% select(tissue,sex,treatment,hypothesis)

# numeric
head(traitData)
traitData %<>% mutate_if(is.factor,as.numeric)
head(traitData)


# remove columns or add that hold information we do not need, if necessary

allTraits <- traitData
allTraits$sample <- colData$sample
row.names(allTraits) <- allTraits$sample
# subset to keep only the "good samples"
Samples <- rownames(datExpr)
traitRows <-  match(Samples, allTraits$sample)
datTraits <-  allTraits[traitRows, ]
datTraits$sample <- NULL


collectGarbage()


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Pituitary dendrogram with trait heatmap",
                    cex.dendroLabels = 0.4
                    )




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
countData <- read.csv("../results/00_countData_manipluation.csv", row.names = 1)

# read the sample meta data or column data
colData <- read.csv("../metadata/00_colData_manipluation.csv", row.names = 1, stringsAsFactors = T)

geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)

## rename rownames and colnames for better vizualiation
colData <- colData %>%
  mutate(sex = fct_recode(sex,
                          "F" = "female",
                          "M" = "male"),
         tissue = fct_recode(tissue,
                             "pit" = "pituitary",
                             "hyp" = "hypothalamus",
                             "gon" = "gonad"))

# wgcna needs numeric identifiers
colData$ID <- as.numeric(colData$bird)

# create a short sample desctiption, make it the row and columns names of colData and countData
colData$sample <- paste(colData$treatment, colData$sex, colData$tissue,  colData$ID, sep = ".")
  
# set row and count names to be the same
row.names(colData) <- colData$sample
colnames(countData) <- colData$sample

head(countData)
head(colData)

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

traitData <- colData %>% select(tissue,sex,treatment)

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


## part 4 combined 

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================



# Load the WGCNA package
library(WGCNA)
library(magrittr) # for %<>%

# The following setting is important, do not omit.

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.60,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


enableWGCNAThreads()


# this is my attempt to solve an error
datExpr %<>% mutate_if(is.integer,as.numeric)
head(str(datExpr))

net = blockwiseModules(datExpr, power = 10,
                       verbose = 5)
names(net)
net$colors

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
mergedColors
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE)

# save data frame with genes and their module colors
genes_modules <- as.data.frame(net$unmergedColors)

# find candidate genes, PRL and PRLR
# PRL = NP_990797.2
# PRLR = XP_015132722.1

genes_modules$entrezid <- row.names(genes_modules)

prolatin_modules <- genes_modules %>%
  filter(entrezid %in% c("NP_990797.2", "XP_015132722.1"))
prolatin_modules

PRL_associated <- genes_modules %>% filter(`net$unmergedColors` %in% c("black"))
PRL_associated <- left_join(PRL_associated, geneinfo, by = "entrezid")
PRL_associated <- PRL_associated %>% arrange(Name)
str(PRL_associated)
write.csv(PRL_associated, "../results/08_PRL_manipulated.csv", row.names = F)

PRLR_associated <- genes_modules %>% filter(`net$unmergedColors` %in% c("red"))
PRLR_associated <- left_join(PRLR_associated, geneinfo, by = "entrezid")
PRLR_associated <- PRLR_associated %>% arrange(Name)
str(PRLR_associated)
write.csv(PRLR_associated, "../results/08_PRLR_manipulated.csv", row.names = F)


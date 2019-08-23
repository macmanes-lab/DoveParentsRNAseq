

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Load the packages
library(WGCNA)
library(dplyr)
library(magrittr)
library(forcats)

subsetWGCNA <- function(whichgroups){

#Read in the data set
colData <- read.csv("../metadata/00_colData_characterization.csv", row.names = 1, stringsAsFactors = T)
countData <- read.csv("../results/00_countData_characterization.csv", row.names = 1)

# create new grouping for subsets
colData$sextissue <-  paste(colData$sex, colData$tissue, sep = ".")

## rename rownames and colnames for better vizualiation
colData <- colData %>%
  mutate(sex = fct_recode(sex,
                              "F" = "female",
                              "M" = "male"),
         tissue = fct_recode(tissue,
                             "pit" = "pituitary",
                             "hyp" = "hypothalamus",
                             "gon" = "gonad"))

colData$ID <- as.numeric(colData$bird)
colData$sample <- paste(colData$sex, colData$tissue, colData$treatment, colData$ID, sep = ".")
head(colData$sample)

row.names(colData) <- colData$sample
colnames(countData) <- colData$sample



# subset to look within one tissue in one sex
colData <- colData %>%
  dplyr::filter(sextissue %in% whichgroups) %>%
  droplevels()
row.names(colData) <- colData$sample
head(colData)
  
# which counts to save
savecols <- as.character(colData$sample) 
savecols <- as.vector(savecols) 
  
# save counts that match colData
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

#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
  #   cex.axis = 1.5, cex.main = 2)
#abline(h = 1000000, col = "red");



#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# Plot a line to show the cut

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 1000000, minSize = 3)
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


traitData <- colData[c(5,3,4)]

#traitData$treatment <- factor(traitData$treatment, 
#                              levels =c("control",  "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", "hatch", "n5", "n9"))
#traitData$tissue <- factor(traitData$tissue, 
 #                             levels =c("hypothalamus", "pituitary", "gonad"))


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
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap",
                    cex.dendroLabels = 0.4
                    )

}


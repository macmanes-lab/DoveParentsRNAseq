#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Load the packages
library(WGCNA)
library(dplyr)
library(magrittr)

#Read in the data set
colData <- read.csv("../metadata/00_colData_characterization.csv", row.names = 1, stringsAsFactors = T)
countData <- read.csv("../results/00_countData_characterization.csv", row.names = 1)

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


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Flagging genes and samples with too many missing values...
# ..step 1
#..Excluding 15 genes from the calculation due to too many missing samples or zero variance.
#..step 2

# Removing genes: NP_001005571.1, NP_001292076.1, NP_989761.2, NP_990385.1, XP_015129157.1, XP_015129381.1, XP_015130369.1, XP_015130427.1, XP_015130613.1, XP_015136658.1, XP_015145985.1, XP_015149350.1, XP_015151860.1, XP_015152548.1, XP_425714.3

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# Plot a line to show the cut
abline(h = 20000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 2000000, minSize = 10)
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


traitData <- colData[c(2:6)]

traitData$treatment <- factor(traitData$treatment, 
                              levels =c("control",  "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", "hatch", "n5", "n9"))
traitData$tissue <- factor(traitData$tissue, 
                              levels =c("hypothalamus", "pituitary", "gonad"))



# solution
str(traitData)
traitData %<>% mutate_if(is.factor,as.numeric)
str(traitData)


# remove columns that hold information we do not need, if necessary
allTraits <- traitData
row.names(allTraits) <- colData$V1
allTraits$samples <- colData$V1

head(allTraits)

# subset to keep only the "good samples"
Samples <- rownames(datExpr)
traitRows <-  match(Samples, allTraits$samples)
datTraits <-  allTraits[traitRows, ]
rownames(datTraits) <-  allTraits[traitRows, 6];  #using column 5 as row name
datTraits$samples <- NULL

collectGarbage()


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


save(datExpr, datTraits, file = "08_WGCNA_1.RData")




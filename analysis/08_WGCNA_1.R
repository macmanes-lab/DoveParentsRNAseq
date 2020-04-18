# WGCNA pipeline

# Load the packages
library(tidyverse)
library(WGCNA)
library(dplyr)
library(magrittr)
library(forcats)

source("../R/functions.R")

# read the sample meta data or column data
colData <- read_csv("../metadata/00_colData_characterization.csv") %>%
  mutate(sex = factor(sex),
         tissue = factor(tissue),
         treatment = factor(treatment)) %>%
  mutate(sex = fct_recode(sex, "F" = "female", "M" = "male"),
         tissue = fct_recode(tissue,
                             "P" = "pituitary",
                             "H" = "hypothalamus",
                             "G" = "gonad"),
         H1 = fct_collapse(treatment,
                                 "NP" = c("control", "bldg"),
                                 "eggs" = c("lay","inc.d3","inc.d9","inc.d17"),
                                 "chicks" = c("hatch","n5", "n9")),
         H2 = fct_collapse(treatment,
                         "NP" = c("control", "bldg"),
                         "early" = c("lay","inc.d3","inc.d9"),
                         "late" = c("inc.d17", "hatch","n5", "n9"))) %>%
  mutate(id = as.numeric(factor(bird))) %>% # wgcna needs numeric identifiers
  select(id, bird, group, sex, tissue, treatment, H1, H2) %>%
  mutate(sample = paste(treatment, sex, tissue, id, sep = "."),
         sextissue = paste(sex, tissue, sep = "."))
colData <- as.data.frame(colData)
row.names(colData) <- colData$sample
head(colData)


# read the count data
countData <- read.csv("../results/00_countData_characterization.csv", 
                      row.names = 1) 
colnames(countData) <- colData$sample

# subset col and count data for each tissue
colDataHyp <- subsetcolData3(colData, c("M.H", "F.H"))
colDataPit <- subsetcolData3(colData, c("M.P", "F.P"))
colDataGon <- subsetcolData3(colData, c("M.G", "F.G"))

countDataHyp <- subsetcountData3(colDataHyp)
countDataPit <- subsetcountData3(colDataPit)
countDataGon <- subsetcountData3(colDataGon)

##### first pass, pituitary

# WGCNA

datExpr0 <- as.data.frame(t(countDataGon))
head(names(datExpr0))  # columns are genes
head(rownames(datExpr0)) # rows are samples

gsg = goodSamplesGenes(datExpr0, verbose = 0); # check good genes
gsg$allOK 

# remove bad genes
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


# detect outlier samples
sampleTree = hclust(dist(datExpr0), method = "average");
plot(sampleTree, main = "Sample clustering to detect outliers", 
     sub="", xlab="", cex.main = 1, cex = 0.4)
abline(h = 2000000, col = "red");

clust = cutreeStatic(sampleTree, cutHeight = 2000000, minSize = 3)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# keep only the "good samples" and important variables

traitData <- colData %>% select(tissue,sex,treatment,H1, H2)
traitData %<>% mutate_if(is.factor,as.numeric)
allTraits <- traitData
allTraits$sample <- colData$sample
row.names(allTraits) <- allTraits$sample
Samples <- rownames(datExpr)
traitRows <-  match(Samples, allTraits$sample)
datTraits <-  allTraits[traitRows, ]
datTraits$sample <- NULL
head(datTraits)

collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Gonad dendrogram with trait heatmap",
                    cex.dendroLabels = 0.4
                    )


# Choose a set of soft-thresholding powers	
powers = c(c(1:10), seq(from = 12, to=20, by=2))	
# Call the network topology analysis function	
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)	
# Plot the results:	

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


enableWGCNAThreads()	

# this is my attempt to solve an error	
datExpr %<>% mutate_if(is.integer,as.numeric)	
head(str(datExpr))	
net = blockwiseModules(datExpr, power = 8,	
                       verbose = 5)	
names(net)	
head(net$colors)
head(net$unmergedColors)

# Convert labels to colors for plotting	
mergedColors = labels2colors(net$colors)	
unmergedColors	= labels2colors(net$unmergedColors)	

# Plot the dendrogram and the module colors underneath	
plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]],	
                    "Module colors",	
                    dendroLabels = FALSE)	

# save data frame with genes and their module colors	
genes_modules <- as.data.frame(net$colors)	

genes_modules %>%
  group_by(`net$colors`) %>%
  summarize(n = n()) %>% arrange(n) 

genes_modules$gene <- row.names(genes_modules)	

candidategenes <- read_csv("../results/candidategenes.csv") %>% pull(x)

candidategenemodules <- genes_modules %>%	
  filter(gene %in% candidategenes)	 %>%
  arrange(`net$colors`, gene) %>%
  distinct(`net$colors`) %>% pull(`net$colors`) %>%
  droplevels()
candidategenemodules

candidategenesdf <- genes_modules %>%	
  filter(gene %in% candidategenes)	 %>%
  arrange(`net$colors`, gene) 
candidategenesdf

candidategeneassociated <- genes_modules %>% 
  filter(`net$colors` %in% candidategenemodules)  %>% 
  filter(!`net$colors` %in% c("turquoise", "grey")) %>%
  filter(!grepl('LOC', gene)) %>%
  arrange(`net$colors`, gene) %>%
  group_by(`net$colors`) %>%
  summarize(genes = str_c(gene, collapse = ", "))
head(candidategeneassociated)	


### save files
write.csv(candidategeneassociated, "../results/candidategeneassociated.csv", row.names = F)
write.csv(genes_modules, "../results/08_genes_modules.csv", row.names = F)
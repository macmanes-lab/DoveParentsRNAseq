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

PRL_associated <- genes_modules %>% filter(`net$unmergedColors` %in% c("greenyellow"))
PRL_associated <- left_join(PRL_associated, geneinfo, by = "entrezid")
PRL_associated <- PRL_associated %>% arrange(Name)
PRL_associated
write.csv(PRL_associated, "../results/08_PRL_associated.csv", row.names = F)

PRLR_associated <- genes_modules %>% filter(`net$unmergedColors` %in% c("magenta"))
PRLR_associated <- left_join(PRLR_associated, geneinfo, by = "entrezid")
PRLR_associated <- PRLR_associated %>% arrange(Name)
PRLR_associated
write.csv(PRLR_associated, "../results/08_PRLR_associated.csv", row.names = F)


VIPR1_associated <- genes_modules %>% filter(`net$unmergedColors` %in% c("yellow"))
VIPR1_associated <- left_join(VIPR1_associated, geneinfo, by = "entrezid")
VIPR1_associated <- VIPR1_associated %>% arrange(Name)
VIPR1_associated
write.csv(VIPR1_associated, "../results/08_VIPR1_associated.csv", row.names = F)

AR_associated <- genes_modules %>% filter(`net$unmergedColors` %in% c("green"))
AR_associated <- left_join(AR_associated, geneinfo, by = "entrezid")
AR_associated <- AR_associated %>% arrange(Name)
AR_associated
write.csv(AR_associated, "../results/08_AR_associated.csv", row.names = F)

POMC_associated <- genes_modules %>% filter(`net$unmergedColors` %in% c("darkred"))
POMC_associated <- left_join(POMC_associated, geneinfo, by = "entrezid")
POMC_associated <- POMC_associated %>% arrange(Name)
POMC_associated
write.csv(POMC_associated, "../results/08_POMC_associated.csv", row.names = F)



#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]

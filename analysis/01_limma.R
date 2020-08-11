# Limma

# DESeq2 is _not_ recommended for experiments with more than 100 samples ([see Mike Love's post](https://mikelove.wordpress.com/2016/09/28/deseq2-or-edger/)), so I decided to try the limma package. I followed [this tutorial](https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html).

library(limma)
library(Glimma)
library(edgeR)

# This script follows 00_datawrangling.R

# import "colData" which contains sample information and "countData" which contains read counts
colData <- read.csv("metadata/00_colData.csv", header = T, row.names = 1)
countData <- read.csv("results/00_counts.csv", header = T, row.names = 1)
head(colData)
head(countData[1:3])

print("subset for quick run")
countData  <- head(countData, 1000) # suset for quick analysis

geneinfo <- row.names(countData)

# Then, I followed the steps from <https://github.com/macmanes-lab/RockDove/blob/master/parental_care/parental_analysis.Rmd>.

head(colData$group)
print("create a large DGEList with 3 elements")
parentalobject <- DGEList(counts=countData, genes=geneinfo, group=colData$group)

print("transform raw counts to countspermillion")
cpms <- cpm(parentalobject)

print("calculate number of lowly lowly expressed genes and remove them")
table(rowSums(parentalobject$counts==0)==10)
keep_genes <- rowSums(cpms >= 1) >= 10
dge <- parentalobject[keep_genes, ]

print("specific the design")
parentaldesign <- model.matrix(~ colData$group )
colnames(parentaldesign) <- levels(colData$group)

print("The TMM normalization")
parentalobject <- calcNormFactors(parentalobject)
parentalobject <- estimateCommonDisp(parentalobject)
parentalobject <- estimateTagwiseDisp(parentalobject)
parentalobject <- estimateDisp(parentalobject, parentaldesign)
parentalobject <- estimateGLMCommonDisp(parentalobject, parentaldesign, verbose=TRUE)
parentalobject <- estimateGLMTrendedDisp(parentalobject, parentaldesign)
parentalobject <- estimateGLMTagwiseDisp(parentalobject, parentaldesign)

print("find and print data")
names(parentalobject)
head(parentalobject$counts[1:3])
head(parentalobject$pseudo.counts[1:3])

print("save output")
write.csv(parentalobject$pseudo.counts, "results/01_limma.csv")
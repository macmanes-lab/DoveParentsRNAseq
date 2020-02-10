# data wrangle for multiple figures

# this is for tsne and pca
# Note, pseudocounts is too big for storage on github :(
pseudocounts <- read_csv("../results/01_pseudo.counts.csv")
head(pseudocounts[1:3])
pseudocounts <- as.data.frame(pseudocounts)
row.names(pseudocounts) <- pseudocounts$X1
pseudocounts$X1 <- NULL
# prep count data for all samples
countData <- as.data.frame(t(pseudocounts))
head(countData[1:3])

# this is for tsne and pca
# prep col data for all samples
colData <- read.csv("../metadata/00_samples.csv", header = T, row.names = 1)
colData$treatment <- factor(colData$treatment, levels = alllevels)
colData <- colData %>% mutate(tissue = fct_recode(tissue, "gonads" = "gonad"))
colData$tissue <- factor(colData$tissue, levels = tissuelevels)
row.names(colData) <- colData$V1

# check ready for analysis
#row.names(countData) == row.names(colData)

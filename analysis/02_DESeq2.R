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

print("subset for quick run")
countData  <- head(countData, 1000) # suset for quick analysis

# gene information
geneinfo <- read.csv("metadata/00_geneinfo.csv", row.names = 1)

# col data or variable informaiton
colData <- read.csv("metadata/00_colData.csv", header = T, row.names = 1) %>%
  mutate(tissue = factor(tissue, levels = tissuelevels),
         treatment <- factor(treatment, levels = alllevels)) %>%
  mutate(sextissue = as.factor(paste(sex, tissue, sep = "_"))) 
colData <- colData %>% drop_na()
# 987 samples, 8 variables

print(ncol(countData) == nrow(colData))


# differential gene expression


createDEGdftreatment <- function(up, down, mytissue){
  
  res <- results(dds, contrast = c("treatment", up, down), 
                 independentFiltering = T, alpha = 0.1)
  
  DEGs <- data.frame(gene = row.names(res),
                        padj = res$padj, 
                        logpadj = -log10(res$padj),
                        lfc = res$log2FoldChange,
                        sextissue = mytissue)
  DEGs <- na.omit(DEGs)
  DEGs <- DEGs %>%
    dplyr::mutate(direction = ifelse(DEGs$lfc > 0 & DEGs$padj < 0.1, 
                                     yes = up, no = ifelse(DEGs$lfc < 0 & DEGs$padj < 0.1, 
                                                           yes = down, no = "NS"))) %>% 
    dplyr::arrange(desc(lfc)) %>%
    dplyr::mutate(direction = factor(DEGs$direction, levels = c(down, "NS", up)))
  
  # write DEGsframe of only significant genes
  DEGs <- DEGs %>% dplyr::filter(direction != "NS")
  print(str(DEGs))
  
  partialfilename = paste("_", down, "_", up, sep = "")
  myfilename = paste0("results/DESeq2/", mytissue, partialfilename, "_DEGs.csv")
  
  write.csv(DEGs, myfilename, row.names = F)
  # return DEGs frome with all data, included NS genes
  # print(head(DEGs))
}  


for(i in levels(colData$sextissue)){
  
  newcolData <- subsetcolData2(colData, i)
  
  # save counts that match colData
  savecols <- as.character(newcolData$V1) 
  savecols <- as.vector(savecols) 
  
  newcountData <- countData %>% dplyr::select(one_of(savecols)) 
  
  dds <- DESeqDataSetFromMatrix(countData = newcountData,
                                colData = newcolData,
                                design = ~ treatment )
  dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
  print(dds)
  print(dim(dds))
  dds <- DESeq(dds, parallel = TRUE) # Differential expression analysis
  
  vsd <- as.data.frame(assay(vst(dds, blind=FALSE)))
  
  myfilename = paste0("results/DEseq2/", i, "_vsd.csv")
  write.csv(vsd, myfilename)
  
  #return(dds)
  #return(vsd)
  #print(head(vsd))

  # save differential gene expression results
  
  control.bldg <- createDEGdftreatment("bldg", "control", i)
  
}


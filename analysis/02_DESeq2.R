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


#### comment this out to run the whole thing
#print("subset for quick run")
#countData  <- head(countData, 1500) # suset for quick analysis

# col data or variable informaiton
colData <- read.csv("metadata/00_colData.csv", header = T, row.names = 1) %>%
  mutate(tissue = factor(tissue, levels = tissuelevels),
         treatment <- factor(treatment, levels = alllevels)) %>%
  mutate(sextissue = as.factor(paste(sex, tissue, sep = "_"))) 
colData <- colData %>% drop_na()
# 987 samples, 8 variables

print(ncol(countData) == nrow(colData))

# differential gene expression


returndds <- function(whichgroup){
  
  newcolData <- colData %>%
    dplyr::filter(sextissue %in% whichgroup) %>%      
    droplevels()
  row.names(newcolData) <- newcolData$V1
  
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
  return(dds)
}


returnvsd <- function(whichdds, whichgroup){

  vsd <- as.data.frame(assay(vst(whichdds, blind=FALSE)))
  myfilename = paste0("results/DEseq2/treatment/", whichgroup, "_vsd.csv")
  write.csv(vsd, myfilename)
  return(vsd)
  print(head(vsd))
  
}

ddsMG <- returndds("male_hypothalamus")
#ddsFG <- returndds("female_hypothalamus")
#ddsMG <- returndds("male_gonads")
#ddsFH <- returndds("female_gonads")
#ddsMP <- returndds("male_pituitary")
#ddsFP <- returndds("female_pituitary")

vsdMG <- returnvsd(ddsMG, "male_hypothalamus")


# DEGs

downs <- c("control", "bldg")
ups   <- c(charlevelsnocontrol, maniplevels)

createDEGdftreatmentvcontrols <- function(whichdds, whichgroup){
  
  for(down in downs){
    for(up in ups){
      if (up != down) {
        
        k <- paste(down, up, sep = " vs. ") #assigns usique rownames
        print(k)
  
        res <- results(whichdds, contrast = c("treatment", up, down), 
                       independentFiltering = T, alpha = 0.1)
        
        DEGs <- data.frame(gene = row.names(res),
                           padj = res$padj, 
                           logpadj = -log10(res$padj),
                           lfc = res$log2FoldChange,
                           sextissue = whichgroup)
        DEGs <- na.omit(DEGs)
        
        DEGs <- DEGs %>%
          dplyr::mutate(direction = ifelse(lfc > 0 & padj < 0.1, yes = up, 
                                           ifelse(lfc < 0 & padj < 0.1, yes = down, 
                                                  no = "NS"))) %>% 
          dplyr::arrange(desc(lfc)) %>%
          dplyr::mutate(direction = factor(direction, levels = c(down, "NS", up)))
        
        # write DEGsframe of only significant genes
        DEGs <- DEGs %>% dplyr::filter(direction != "NS")
        print(str(DEGs))
        
        partialfilename = paste("_", down, "_", up, sep = "")
        myfilename = paste0("results/DESeq2/treatment/", 
                            whichgroup, partialfilename, "_DEGs.csv")
        
        write.csv(DEGs, myfilename, row.names = F)
        # return DEGs frome with all data, included NS genes
        # print(head(DEGs))
        
        print(head(DEGs))
      
      }
    }  
  }
  
  up <- up[-1] 
  down <- down[-1] 
}

createDEGdftreatmentvcontrols(ddsMG, "male_hypothalamus")



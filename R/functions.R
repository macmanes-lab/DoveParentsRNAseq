# run DESeq on subset of data
# e.g. dds <- subsetDESeq("male_hypothalamus")
subsetDESeq <- function(eachgroup){
  
  # subset to look within one tissue in one sex
  colData <- m.colData %>%
    dplyr::filter(sextissue == eachgroup) %>%
    droplevels()
  row.names(colData) <- colData$V1
  
  # which counts to save
  savecols <- as.character(colData$V1) 
  savecols <- as.vector(savecols) 
  
  # save counts that match colData
  countData <- m.countData %>% dplyr::select(one_of(savecols)) 
  
  # check that row and col lenghts are equal
  print(ncol(countData) == nrow(colData))
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ treatment )
  
  print(dds)
  dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
  print(dim(dds))
  
  dds <- DESeq(dds) # Differential expression analysis
  return(dds)
}


# print total number of differntially expressed genes
# numDEGs('m.inc.d3', 'm.inc.d9')


numDEGs <- function(dds, group1, group2){
  res <- results(dds, contrast = c("treatment", group1, group2), independentFiltering = T)
  sumpadj <- sum(res$padj < 0.01, na.rm = TRUE)
  return(sumpadj)
}


## plot DEGs 

plottotalDEGs <- function(dds, mysubtitle){
  
  a <- group1
  b <- group2
    
    # comapre all contrasts, save to datafrmes
    totalDEGS=data.frame()
    for (i in a){
      for (j in b){
        if (i != j) {
          k <- paste(i,j, sep = ".") #assigns usique rownames
          #print(k)
          totalDEGS[k,1]<-i               
          totalDEGS[k,2]<-j
          totalDEGS[k,3]<- numDEGs(dds, i,j) #caluculates number of DEGs
        }
      }
      b <- b[-1]  # drop 1st element of second string to not recalculate DEGs
    }
    
  
  totalDEGS$V1 <- factor(totalDEGS$V1, levels = 
                           c("m.inc.d3",  "m.inc.d8",
                             "m.inc.d9", "m.inc.d17",
                             "prolong", "extend", "m.n2"))
  totalDEGS$V2 <- factor(totalDEGS$V2, levels = 
                           c("m.inc.d3",  "m.inc.d8",
                             "m.inc.d9", "m.inc.d17",
                             "prolong", "extend", "m.n2"))
  
  allcontrasts <- totalDEGS %>%
    ggplot( aes(V1, V2)) +
    geom_tile(aes(fill = V3)) +
    scale_fill_viridis(na.value="#440154", 
                       limits = c(0, 3000),
                       breaks = c(0, 1000, 2000, 3000)) + 
    xlab(" ") + ylab("Timepoint") +
    labs(fill = "# of DEGs",
         subtitle = mysubtitle) +
    theme_minimal(base_size = 8) + 
    theme(axis.text.x = element_text(angle = 90))
  
  return(allcontrasts)
}

# plot DEGs just for characterization study
plottotalDEGschar <- function(dds, mysubtitle){
  
  a <- group1
  b <- group2
  
  # comapre all contrasts, save to datafrmes
  totalDEGS=data.frame()
  for (i in a){
    for (j in b){
      if (i != j) {
        k <- paste(i,j, sep = ".") #assigns usique rownames
        #print(k)
        totalDEGS[k,1]<-i               
        totalDEGS[k,2]<-j
        totalDEGS[k,3]<- numDEGs(dds, i,j) #caluculates number of DEGs
      }
    }
    b <- b[-1]  # drop 1st element of second string to not recalculate DEGs
  }
  
  totalDEGS$V1 <- factor(totalDEGS$V1, levels = 
                           c("control",  "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", "hatch", "n5", "n9"))
  totalDEGS$V2 <- factor(totalDEGS$V2, levels = 
                           c("control",  "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", "hatch", "n5", "n9"))
  
  allcontrasts <- totalDEGS %>%
    ggplot( aes(V1, V2)) +
    geom_tile(aes(fill = V3)) +
    scale_fill_viridis(na.value="#440154", 
                       limits = c(0, 3000),
                       breaks = c(0, 1000, 2000, 3000)) + 
    xlab(" ") + ylab("Timepoint") +
    labs(fill = "# of DEGs",
         subtitle = mysubtitle) +
    theme_minimal(base_size = 8) + 
    theme(axis.text.x = element_text(angle = 90))
  
  return(allcontrasts)
}


# resturn pvalues for all genes
returnpadj <- function(group1, group2){
  res <- results(dds, contrast = c("treatment", group1, group2), independentFiltering = T)
  pvals <- as.data.frame(res$padj)
  padjcolname <- as.character(paste("padj", group1, group2, sep=""))
  colnames(pvals) <- c(padjcolname)
  return(pvals)
}

## calculate principal components
# I've taken the pca function from DESeq2 and elaborated it so that I could extract up to 6 PCs
pcadataframe <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4],PC5 = pca$x[, 5],PC6 = pca$x[, 6],group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:6]
    return(d)
  }
}

# plot pcs 
# eg. plotPCAs(dds.female_hypothalamus, "female hypothalamus")
plotPCAs <- function(dds, mysubtitle){
  
  vsd <- vst(dds, blind=FALSE) # variance stabilized 
  
  # create the dataframe using my function pcadataframe
  pcadata <- pcadataframe(vsd, intgroup=c("treatment"), returnData=TRUE)
  percentVar <- round(100 * attr(pcadata, "percentVar"))
  print(percentVar)
  
  print(summary(aov(PC1 ~ treatment, data=pcadata)))
  print(TukeyHSD(aov(PC1 ~ treatment, data=pcadata), which = "treatment"))
  
  print(summary(aov(PC2 ~ treatment, data=pcadata))) 
  print(summary(aov(PC3 ~ treatment, data=pcadata))) 
  print(summary(aov(PC4 ~ treatment, data=pcadata))) 
  
  pca1 <- ggplot(pcadata, aes(treatment, PC1,color = treatment)) + 
    geom_boxplot() +
    geom_point() +
    ylab(paste0("PC1: ", percentVar[1],"% variance")) +
    xlab(NULL) +
    theme_cowplot(font_size = 8, line_size = 0.25) +
    labs(subtitle = " ") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90))
  
  
  pca2 <- ggplot(pcadata, aes(treatment, PC2,color = treatment)) + 
    geom_boxplot() +
    geom_point() +
    ylab(paste0("PC2: ", percentVar[2],"% variance")) +
    xlab(NULL) +
    theme_cowplot(font_size = 8, line_size = 0.25) +
    labs(subtitle = " ") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90))
  
  pca12 <- ggplot(pcadata, aes(PC1, PC2,color = treatment)) + 
    geom_point() +
    stat_ellipse() +
    ylab(paste0("PC2: ", percentVar[2],"% variance")) +
    xlab(paste0("PC1: ", percentVar[1],"% variance")) +
    theme_cowplot(font_size = 8, line_size = 0.25) +
    labs(subtitle = mysubtitle) +
    theme(legend.position = "none")
  
  mypca <- plot_grid(pca12, pca1, pca2, nrow = 1)
  return(mypca)
}

## plot candidate genes 
# e.g. plotcandidates(dds.female_hypothalamus, "female hypothalamus")

plotcandidates <- function(mydds, mysubtitle){
  
  vsd <- vst(mydds, blind=FALSE) # variance stabilized 
  
  # make dataframe counts
  DEGs <- assay(vsd)
  DEGs <- as.data.frame(DEGs)
  
  names(geneinfo)
  names(geneinfo)[4] <- "rownames"
  DEGs$rownames <- row.names(DEGs)
  
  # make dataframe with geneids and names and counts
  # how to gather: https://tidyr.tidyverse.org/reference/gather.html
  
  candidates <- full_join(geneinfo, DEGs)
  drop.cols <-colnames(candidates[,grep("padj|pval|pmin", colnames(candidates))])
  candidates <- candidates %>% dplyr::select(-one_of(drop.cols))
  candidates <- candidates %>%
    filter(Name %in% c("AR", "CYP19A1", "ESR1", "ESR2", "FSHR",
                       "GHRL", "GAL", "NPVF", "GNRH1", "LHCGR",
                       "PGR", "PRL", "PRLR", "VIP", "VIPR1")) 
  row.names(candidates) <- candidates$Name
  candidates <- candidates %>% select(-row.names, -rownames, -Name, -geneid)
  candidates <- candidates %>% drop_na()
  candidates <- as.data.frame(t(candidates))
  candidates$RNAseqID <- rownames(candidates)
  candidates <- candidates %>% gather(gene, value, -RNAseqID)  %>% 
    filter(RNAseqID != "gene")
  candidates$value <- as.numeric(candidates$value)
  candidates$V1  <- candidates$RNAseqID
  
  candidatecounts <- left_join(candidates, m.colData)
  candidatecounts$faketime <- as.numeric(candidatecounts$treatment)
  candidatecounts$gene <- as.factor(candidatecounts$gene)
  
  p1 <- ggplot(candidatecounts, aes(x = treatment, y = value, fill = treatment)) +
    geom_boxplot() +
    facet_wrap(~gene, scales = "free") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank(),
          legend.position = "bottom") +
    labs(x = NULL, subtitle = mysubtitle) +
    guides(fill= guide_legend(nrow=1)) 
  return(p1)
}


## make pheatmaps 
# e.g. makepheatmap(dds.female_hypothalamus, "female hypothalamus")

makepheatmap <- function(mydds, mysubtitle){
  
  vsd <- vst(mydds, blind=FALSE) # variance stabilized 
  
  # make dataframe counts
  DEGs <- assay(vsd)
  DEGs <- as.data.frame(DEGs)
  
  a <- levels(m.colData$treatment)
  b <- a
  
  for (i in a){
    for (j in b){
      if (i != j) {
        results <- returnpadj(i,j)
        DEGs <- cbind(DEGs,results)
      }
    }
    b <- b[-1]  # drop 1st element of second string to not recalculate DEGs
  }
  
  DEGsmatrix <- DEGs
  DEGsmatrix <- as.matrix(DEGs)
  padjmin <- rowMins(DEGsmatrix, na.rm = T) 
  padjmin <- as.data.frame(padjmin)
  
  sigDEGs <- cbind(DEGs,padjmin)
  sigDEGs <- sigDEGs %>% arrange(padjmin)
  sigDEGs <- head(sigDEGs,500)
  sigDEGs <- as.data.frame(sigDEGs)
  rownames(sigDEGs) <- sigDEGs$rownames
  drop.cols <-colnames(sigDEGs[,grep("padj|pval|pmin|rownames", colnames(sigDEGs))])
  sigDEGs <- sigDEGs %>% dplyr::select(-one_of(drop.cols))
  sigDEGs <- as.matrix(sigDEGs)
  sigDEGs <- sigDEGs - rowMeans(sigDEGs)
  
  paletteLength <- 30
  myBreaks <- c(seq(min(sigDEGs), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(sigDEGs)/paletteLength, max(sigDEGs), length.out=floor(paletteLength/2)))
  
  anndf <- m.colData %>% dplyr::select(treatment)
  rownames(anndf) <- m.colData$V1
  
  sigDEGs <- as.matrix(sigDEGs) 
  p1 <- pheatmap(sigDEGs, show_rownames = F, show_colnames = F,
           color = viridis(30),
           breaks=myBreaks,
           annotation_col=anndf,
           main = mysubtitle)
  return(p1)
  }



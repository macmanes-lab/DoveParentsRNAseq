# run DESeq on subset of data
# e.g. dds <- subsetDESeq("male_hypothalamus")
subsetDESeq <- function(colData, countData, eachgroup){
  
  # subset to look within one tissue in one sex
  colData <- colData %>%
    dplyr::filter(sextissue == eachgroup) %>%
    droplevels()
  row.names(colData) <- colData$V1
  
  # which counts to save
  savecols <- as.character(colData$V1) 
  savecols <- as.vector(savecols) 
  
  # save counts that match colData
  countData <- countData %>% dplyr::select(one_of(savecols)) 
  
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

returntotalDEGs <- function(dds){
  
  colData <- a.colData # %>%
   # filter(treatment %in% c("control", "bldg", "n9",
   #                         "inc.d17", "m.inc.d17",
   #                         "hatch",  "m.n2")) %>%
   # droplevels()
  
  group1 <- levels(colData$treatment)

  a <- group1
  b <- group1
  
  # comapre all contrasts, save to datafrmes
  totalDEGS=data.frame()
  for (i in a){
    for (j in b){
      if (i != j) {
        k <- paste(i,j, sep = ".") #assigns usique rownames
        print(k)
        totalDEGS[k,1]<-i               
        totalDEGS[k,2]<-j
        totalDEGS[k,3]<- numDEGs(dds, i,j) #caluculates number of DEGs
      }
    }
    b <- b[-1]  # drop 1st element of second string to not recalculate DEGs
  }
  
 print(totalDEGS)  
 return(totalDEGS)
  
}

plottotalDEGs <- function(myDEGS, mysubtitle){  
  totalDEGS <- myDEGS
  totalDEGS$V1 <- factor(totalDEGS$V1, levels =  c("control", "bldg", "lay",
                                                   "inc.d3", "inc.d9", "inc.d17",
                                                   "hatch", "n5", "n9",
                                                   
                                                   "m.inc.d3", "m.inc.d8", "m.inc.d9",
                                                   "m.inc.d17",  "m.n2"  ,
                                                   "prolong", "extend" ))
  totalDEGS$V2 <- factor(totalDEGS$V2, levels =  c("control", "bldg", "lay",
                                                   "inc.d3", "inc.d9", "inc.d17",
                                                   "hatch", "n5", "n9",
                                                   
                                                   "m.inc.d3", "m.inc.d8", "m.inc.d9",
                                                   "m.inc.d17",  "m.n2"  ,
                                                   "prolong", "extend" ))

  totalDEGS <- totalDEGS %>% dplyr::na_if(0)
  
  print(str(totalDEGS))
  
  allcontrasts <- totalDEGS %>%
    ggplot( aes(V1, V2)) +
    geom_tile(aes(fill = V3)) +
    theme_minimal(base_size = 12) + 
    geom_text(aes(label = round(V3, 1)), color = "black")+
    scale_fill_viridis(na.value="#bdbdbd", 
                       limits = c(0,7000)) +
    xlab(NULL) + ylab(NULL) +
    labs(fill = "# of DEGs",
         title = mysubtitle, subtitle = "  ", caption = "  ") +
    theme(axis.text.x = element_text(angle = 90))
  print(totalDEGS)
  plot(allcontrasts)
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
  
  totalDEGS$V1 <- as.factor(totalDEGS$V1)
  totalDEGS$V2 <- as.factor(totalDEGS$V2)
  
  allcontrasts <- totalDEGS %>%
    theme_minimal(base_size = 12) + 
    ggplot( aes(V1, V2)) +
    geom_tile(aes(fill = V3), size=6) +
    scale_fill_viridis(na.value="#440154") + 
    xlab(" ") + ylab("Treatment") +
    labs(fill = "# of DEGs",
         subtitle = mysubtitle) +
    
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
pcadataframe <- function (object, intgroup, ntop = 500, returnData = FALSE) 
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
returnPCAs <- function(dds){
  
  dds <- dds
  
  vsd <- vst(dds, blind=FALSE) # variance stabilized 
  
  # create the dataframe using my function pcadataframe
  pcadata <- pcadataframe(vsd, intgroup=c("lastday", "penultimate", "xlabel", "study"), returnData=TRUE)
  percentVar <- round(100 * attr(pcadata, "percentVar"))
  print(percentVar)
  
  print(summary(aov(PC1 ~ xlabel, data=pcadata)))
  #print(TukeyHSD(aov(PC1 ~ xlabel, data=pcadata), which = "xlabel"))
  print(summary(aov(PC2 ~ xlabel, data=pcadata))) 
  print(summary(aov(PC3 ~ xlabel, data=pcadata))) 
  print(summary(aov(PC4 ~ xlabel, data=pcadata))) 
  return(pcadata)
}

plotPC1 <- function(pcadata, mysubtitle, myxlab){ 
  
  pcadata <- pcadata
  
  #percentVar <- round(100 * attr(pcadata, "percentVar"))

  pca1 <- ggplot(pcadata, aes(xlabel, PC1, fill = lastday, color = penultimate)) + 
    geom_violin() +
    #geom_point() +
    theme_bw(base_size = 12) +
    ylab(myxlab) +
    xlab("Parental stages, with increasing days ->") +
    labs(subtitle = mysubtitle) +
    scale_fill_manual(values = colorlastday, 
                      guide = guide_legend(
                        direction = "horizontal",
                        label.position = "bottom")) +
    scale_color_manual(values = colorpenultimate,
                       guide = guide_legend(
                         direction = "horizontal",
                         label.position = "bottom")) +
    facet_wrap(~study, scales = "free_x")  +
    guides(col = guide_legend(nrow = 1),
           fill = guide_legend(nrow = 1))
  return(pca1)
}


plotPC2 <- function(pcadata, mysubtitle){ 
    
    pcadata <- pcadata
    
  pca2 <- ggplot(pcadata, aes(xlabel, PC2, color = penultimate, fill = lastday)) + 
    geom_violin() +
    #geom_point() +
    theme_bw(base_size = 12) +
    ylab(paste0("PC2")) +
    xlab("Parental stages, with increasing days ->") +
    labs(subtitle = mysubtitle) +
    scale_fill_manual(values = colorlastday, 
                      guide = guide_legend(
                        direction = "horizontal",
                        label.position = "bottom")) +
    scale_color_manual(values = colorpenultimate,
                       guide = guide_legend(
                         direction = "horizontal",
                         label.position = "bottom")) +
    facet_wrap(~study, scales = "free_x")  +
    guides(col = guide_legend(nrow = 1),
           fill = guide_legend(nrow = 1))
}  
 
plotPC12 <- function(pcadata, mysubtitle){ 
  
  pcadata <- pcadata
  
  pca12 <- ggplot(pcadata, aes(PC3, PC2, color = penultimate ,fill = lastday)) + 
    geom_point(pch=21, size = 3) +
    #stat_ellipse() +
    theme_bw(base_size = 12) +
    ylab(paste0("PC2")) +
    xlab(paste0("PC3")) +
    labs(subtitle = mysubtitle) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = colorpenultimate) +
    scale_fill_manual(values = colorlastday)  +
    guides(color = guide_legend(order = 1, ncol=2), 
           fill = guide_legend(order = 2, ncol=2)) +
    theme(legend.text = element_text(size=10)) 
  
}  
  #legend <- get_legend(pca1)
  #mypcatop <- plot_grid(pca1 + theme(legend.position = "none"), pca2, nrow = 1)  
  #mypca <- plot_grid(mypcatop, legend, ncol = 1, rel_heights = c(1,0.2))
  
## plot candidate genes 
# e.g. plotcandidates(dds.female_hypothalamus, "female hypothalamus")

plotcandidates <- function(mydds, colData, mysubtitle){
  
  vsd <- vst(mydds, blind=FALSE) # variance stabilized 
  
  # make dataframe counts
  DEGs <- assay(vsd)
  DEGs <- as.data.frame(DEGs)
  
  names(geneinfo)
  names(geneinfo)[4] <- "rownames"
  DEGs$rownames <- row.names(DEGs)
  print(head(DEGS))

  # make dataframe with geneids and names and counts
  # how to gather: https://tidyr.tidyverse.org/reference/gather.html
  
  candidates <- full_join(geneinfo, DEGs)
  candidates <- candidates %>%
    filter(pmin < 0.01) 
  drop.cols <-colnames(candidates[,grep("padj|pval|pmin", colnames(candidates))])
  candidates <- candidates %>% dplyr::select(-one_of(drop.cols))
  candidates <- candidates %>%
    filter(Name %in% c("AVP", "AVPR1A", "AVPR1B", "AVPR2",
                      "SERPINA4", "CRH", "CRHR1", "DRD1", "DRD2",
                      "DRD3",  "DRD4",  "DRD5",  "GABRQ",  "NR3C1",  "HSD11B1a", 
                      "HSD11B1b", "HSD11B2L",  "HSD11B1L",  "MC2R", "NR3C2",
                      "OXTR", "POMC",  "AR", "CYP19A1", "ESR1",  "ESR2", "FSHR", 
                      "FSHB", "GHRL", "GAL", "NPVF",  "NPFFR1", "GNRH1",  "GNRHR", 
                      "LEPR",  "LHCGR", "PGR", "PRL", "PRLR", "VIP", "VIPR1"))
  row.names(candidates) <- candidates$Name
  candidates <- candidates %>% select(-row.names, -rownames, -Name, -geneid)
  candidates <- candidates %>% drop_na()
  candidates <- as.data.frame(t(candidates))
  candidates$RNAseqID <- rownames(candidates)
  candidates <- candidates %>% gather(gene, value, -RNAseqID)  %>% 
    filter(RNAseqID != "gene")
  candidates$value <- as.numeric(candidates$value)
  candidates$V1  <- candidates$RNAseqID
  
  candidatecounts <- left_join(candidates, colData)
  candidatecounts$faketime <- as.numeric(candidatecounts$treatment)
  candidatecounts$gene <- as.factor(candidatecounts$gene)
  
  p1 <- ggplot(candidatecounts, aes(x = treatment, y = value, fill = treatment)) +
    geom_boxplot() +
    facet_wrap(~gene, scales = "free") +
    theme_bw(base_size = 8) +
    theme(axis.text.x = element_blank(),
          legend.position = "bottom") +
    labs(x = NULL, subtitle = mysubtitle) +
    guides(fill= guide_legend(nrow=2)) 
  return(p1)
}

## make pheatmaps 
# e.g. makepheatmap(dds.female_hypothalamus, "female hypothalamus")

makepheatmap <- function(mydds, colData, mysubtitle){
  dds <- mydds
  
  vsd <- vst(dds, blind=FALSE) # variance stabilized 
  
  # make dataframe counts
  DEGs <- assay(vsd)
  DEGs <- as.data.frame(DEGs)
  
  a <- levels(colData$treatment)
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
  
  anndf <- colData %>% dplyr::select(treatment)
  rownames(anndf) <- colData$V1
  
  sigDEGs <- as.matrix(sigDEGs) 
  p1 <- pheatmap(sigDEGs, show_rownames = F, show_colnames = F,
           color = viridis(30),
           breaks=myBreaks,
           annotation_col=anndf,
           main = mysubtitle)
  return(p1)
  }



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
         subtitle = mysubtitle)
  
  return(allcontrasts)
}

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
         subtitle = mysubtitle)
  
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



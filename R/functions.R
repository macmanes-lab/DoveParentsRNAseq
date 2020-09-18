

########## DESeq2 dds and vsd

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


wranglevsds <- function(pathtofile){
  df <- read_csv(pathtofile) %>%
    mutate(treatment = factor(treatment, levels = alllevels))
  return(df)
}



##### DEGs


createDEGdftreatmentvcontrols <- function(whichdds, whichgroup, downs, ups){
  
  for(down in downs){
    for(up in ups){
      if (up != down) {
        if (up != "prolong" | down != "early") {
          
          cat("\n\n")
          
          k <- paste(down, up, sep = " vs. ") #assigns usique rownames
          print(whichgroup)
          print(k)
          
          res <- results(whichdds, contrast = c("treatment", up, down), 
                         independentFiltering = T, alpha = 0.1,
                         lfcThreshold=log2(1.1))
          print(summary(res))
          
          DEGs <- data.frame(gene = row.names(res),
                             padj = res$padj, 
                             logpadj = -log10(res$padj),
                             lfc = res$log2FoldChange,
                             sextissue = whichgroup)
          DEGs <- na.omit(DEGs)
          
          DEGs <- DEGs %>%
            dplyr::mutate(direction = ifelse(lfc > 0.14 & padj < 0.1, yes = up, 
                                             ifelse(lfc < 0.14 & padj < 0.1, yes = down, 
                                                    no = "NS"))) %>% 
            dplyr::arrange(desc(lfc)) %>%
            dplyr::mutate(direction = factor(direction, levels = c(down, "NS", up)))
          
          # write DEGsframe of only significant genes
          DEGs <- DEGs %>% dplyr::filter(direction != "NS")
          
          top6 <- DEGs %>% head(.) %>% pull(gene) 
          print(top6)
          
          # return DEGs frome with all data, included NS genes
          partialfilename = paste("_", down, "_", up, sep = "")
          myfilename = paste0("results/DESeq2/treatment/", 
                              whichgroup, partialfilename, "_DEGs.csv")
          write.csv(DEGs, myfilename, row.names = F)
        }
      }  
    }
    up <- up[-1] 
  }
  down <- down[-1] 
}









## volcano plots

plot.volcano <- function(whichtissue, whichsex,  whichcomparison){
  
  volcano <- allDEG %>%
    dplyr::filter(tissue %in% whichtissue,
           comparison %in% whichcomparison,
           sex %in% whichsex) %>%
    dplyr::mutate(direction = factor(direction, levels = alllevels2)) %>%
    ggplot(aes(x = lfc, y = logpadj)) + 
    geom_point(aes(color = direction), size = 1, 
               alpha = 0.75, na.rm = T) + 
    theme_B3() +
    scale_color_manual(values = allcolors, 
                       name = "increased expression in:",
                       breaks = alllevels) +
    labs(y = expression(-log[10]("FDR")), 
         x = "Log-fold change (LFC)") +
    theme(legend.position = "bottom",
          legend.direction = "vertical",
          legend.key.height = unit(0, "cm"),      
          legend.spacing.y = unit(0, "cm")) +
    guides(color = guide_legend(nrow = 1))
  return(volcano)
}









## tsne fig 1 

subsetmaketsne <- function(whichtissue, whichtreatment, whichsex){
  
  colData <- colData %>%
    dplyr::filter(tissue %in% whichtissue,
                  treatment %in% whichtreatment,
                  sex %in% whichsex) 
  row.names(colData) <- colData$V1
  
  # save counts that match colData
  savecols <- as.character(colData$V1) 
  savecols <- as.vector(savecols) 
  
  countData <- as.data.frame(t(countData))
  countData <- countData %>% dplyr::select(one_of(savecols)) 
  countData <- as.data.frame(t(countData))
  
  euclidist <- dist(countData) # euclidean distances between the rows
  
  tsne_model <- Rtsne(euclidist, check_duplicates=FALSE)
  tsne_df = as.data.frame(tsne_model$Y) 
  
  # prep for adding columns
  colData2 <- colData 
  colData2$V1 <- NULL
  tsne_df_cols <- cbind(colData2, tsne_df)
  return(tsne_df_cols)
}


plottsne <- function(tsnedf, pointcolor, whichcolors){
  p <- ggplot(tsnedf, aes(x = V1, y = V2)) +
    
    geom_point(size = 1, aes(color = pointcolor, shape = sex)) +
    theme_B3() +
    labs(x = "tSNE 1", y = "tSNE 2", 
         title = " ") +
    scale_color_manual(values = whichcolors) +
    theme(legend.position = "none",
          axis.text = element_blank(), axis.ticks = element_blank()) 
  return(p)
}

## bar graphs fig 3

makebargraphv4 <- function(df, whichtissue, myylab, 
                           whichlevels, whichlabels, whichsex,
                           myymin, myymax){
  
  p <- df %>%
    dplyr::filter(tissue == whichtissue,
           sex == whichsex) %>%
    ggplot(aes(x = comparison, y = n, fill = direction)) +
    geom_bar(stat="identity") +
    theme_B3() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))  +
    scale_fill_manual(values = allcolors,
                      name = " ") +
    scale_color_manual(values = allcolors) +
    geom_text(stat='identity', aes(label= n), vjust =-0.5, 
              position = position_dodge(width = 1),
              size = 1.25, color = "black") +
    labs(x = NULL, y = myylab)  +
    scale_x_discrete(breaks = whichlevels,
                     labels = whichlabels,
                     drop = F,
                     position = "bottom") +
   # ylim(myymin, myymax) +
    scale_y_continuous(breaks = c(-1500,-1000, -500, 0, 500, 1000),
                       limits = c(myymin, myymax))
  return(p)
}

makebargraphv5 <- function(df, whichtissue, myylab, 
                           whichlevels, whichlabels, whichsex,
                           myymin, myymax){
  
  p <- df %>%
    dplyr::filter(tissue == whichtissue,
                  sex == whichsex) %>%
    ggplot(aes(x = comparison, y = n, fill = direction)) +
    geom_bar(stat="identity") +
    theme_B3() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))  +
    scale_fill_manual(values = allcolors,
                      name = " ") +
    scale_color_manual(values = allcolors) +
    #geom_text(stat='identity', aes(label= n), vjust =-0.5, 
    #          position = position_dodge(width = 1),
    #         size = 1.25, color = "black") +
    labs(x = NULL, y = myylab)  +
    scale_x_discrete(breaks = whichlevels,
                     labels = whichlabels,
                     drop = F,
                     position = "bottom") +
    scale_y_continuous(breaks = c(-4000,-3000, -2000, -1000, 0, 1000, 2000, 3000),
                       limits = c(myymin, myymax))
    
  return(p)
}












## candidate gene box plot 






## box plots fig 2

plotcandidatechar <- function(df, whichgene,  treatment1, treatment2){
  
  p <- df %>%
    dplyr::filter(gene  == whichgene,
                  treatment %in% charlevels)  %>%
    ggplot(aes(y =  counts, x = treatment, fill = treatment)) +
    geom_boxplot(lwd=0.5 , outlier.size = 0.1) +
    #geom_jitter(size = 0.25) +
    theme_B3() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_text(face = "italic")) +
    scale_fill_manual(values = allcolors)  +
    labs(y = whichgene, x = "Characterization") +
    geom_signif(comparisons = list(c(treatment1, treatment2)),
                map_signif_level=F,
                textsize = 1.5, family = 'Helvetica',
                vjust = 1.5, size = 0.5) 
  return(p)
  
}

plotcandidatemanip <- function(df, whichgene, treatment3, treatment4){
  
  p <- df %>%
    dplyr::filter(gene  == whichgene,
                  treatment %in% c(removallevels, timinglevels))  %>%
    ggplot(aes(y =  counts, x = treatment, fill = treatment)) +
    geom_boxplot(lwd=0.5, outlier.size = 0.1) +
   # geom_jitter(size = 0.25) +
    theme_B3() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_text(face = "italic")) +
    scale_fill_manual(values = allcolors)  +
    labs( y = whichgene,  x = "Removal and replacement") +
    geom_signif(comparisons = list(c(treatment3, treatment4)),
                map_signif_level= F,
                textsize = 1.5, family = 'Helvetica',
                vjust = 1.5,  size = 0.5) 
  return(p)
  
}

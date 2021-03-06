

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
    mutate(treatment = factor(treatment, levels = alllevels2))
  return(df)
}



##### DEGs


createDEGdfs <- function(whichdds, whichgroup, downs, ups){
  
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
                             pvalue = res$pvalue, 
                             padj = res$padj, 
                             logpadj = -log10(res$padj),
                             lfc = res$log2FoldChange,
                             sextissue = whichgroup)
          DEGs <- na.omit(DEGs)
          
          DEGs <- DEGs %>%
            dplyr::mutate(direction = ifelse(lfc > 0.14 & padj < 0.1, yes = up, 
                                             ifelse(lfc < -0.14 & padj < 0.1, yes = down, 
                                                    no = "NS"))) %>% 
            dplyr::mutate(direction2 = ifelse(lfc > 0.14 & pvalue < 0.05, yes = up, 
                                             ifelse(lfc < -0.14 & pvalue < 0.05, yes = down, 
                                                    no = "NS"))) %>% 
            dplyr::arrange(gene) %>%
            dplyr::mutate(direction = factor(direction, levels = c(down, "NS", up)),
                          direction2 = factor(direction2, levels = c(down, "NS", up)))
          

         
          
          # write DEGsframe of only significant genes
          DEGs <- DEGs %>% dplyr::filter(direction2 != "NS") %>%
            select(sextissue, gene, 
                   padj, direction, pvalue, direction2, 
                   lfc, logpadj)
          
          print("First 5 genes with un-adjusted p < 0.5")
          
          DEGs %>%
            select(-sextissue, -logpadj) %>%
            head() %>%
            print()
          
          print("Candidate genes with un-adjusted p < 0.05")
          listogenes <- DEGs %>% arrange(gene) %>%  
            mutate(gene = as.character(gene)) %>%
            filter(gene %in% candidategenes) %>%
            pull(gene) 
          print(listogenes)
          
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

plot.volcano <- function(whichtissue,  whichcomparison){
  
  volcano <- allDEG %>%
    dplyr::filter(tissue %in% whichtissue,
           comparison %in% whichcomparison) %>%
    dplyr::mutate(direction = factor(direction, levels = alllevels2)) %>%
    ggplot(aes(x = lfc, y = logpadj)) + 
    geom_point(aes(color = direction), size = 1, 
               alpha = 0.75, na.rm = T) + 
    theme_B3() +
    facet_wrap(~sex, nrow = 2, scales = "free_y") +
    scale_color_manual(values = allcolors, 
                       name = "higher in:",
                       breaks = alllevels) +
    labs(y = expression(-log[10]("FDR")), 
         x = "LFC",
         title = " ") +
    theme(legend.position = "top",
          legend.direction = "vertical",
          legend.key.height = unit(-0.2, "cm"),      
          #legend.spacing.y = unit(-0.2, "cm"),
          legend.spacing.x = unit(-0.2, "cm")) +
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
                           whichlevels, whichlabels){
  
  p <- df %>%
    dplyr::filter(tissue == whichtissue) %>%
    ggplot(aes(x = comparison, y = n, fill = direction)) +
    geom_hline(yintercept = 0) +
    geom_bar(stat="identity") +
    theme_B3() +
    facet_wrap(~sex, nrow =  2) +
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
                     position = "bottom")  
  return(p)
}

makebargraphv5 <- function(df, whichtissue, myylab, 
                           whichlevels, whichlabels){
  
  p <- df %>%
    dplyr::filter(tissue == whichtissue) %>%
    ggplot(aes(x = comparison, y = n, fill = direction)) +
    geom_hline(yintercept = 0) +
    geom_bar(stat="identity") +
    theme_B3() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))  +
    scale_fill_manual(values = allcolors,
                      name = " ") +
    scale_color_manual(values = allcolors) +
    facet_wrap(~sex, nrow = 2) +
    labs(x = NULL, y = myylab)  +
    scale_x_discrete(breaks = whichlevels,
                     labels = whichlabels,
                     drop = F,
                     position = "bottom")  
    
  return(p)
}



## candidate gene box plot 

plotcandidatechar <- function(df, whichgene){
  
  p <- df %>%
    dplyr::filter(gene  == whichgene,
                  treatment %in% charlevels)  %>%
    ggplot(aes(y =  counts, x = treatment, fill = treatment)) +
    geom_boxplot(lwd=0.5 , outlier.size = 0.1) +
    facet_wrap(~sex, nrow = 2, scales = "free_y") +
    theme_B3() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_text(face = "italic")) +
    scale_fill_manual(values = allcolors)  +
    labs(y = whichgene, x = "Stage", title =" ") +
    geom_signif(comparisons = list(c("bldg", "inc.d17"),
                                   c("inc.d9", "inc.d17")),
                map_signif_level=F, step_increase = 0.1,
                textsize = 1.5)
  return(p)
  
}

plotcandidatemanip <- function(df, whichgene){
  
  p <- df %>%
    dplyr::filter(gene  == whichgene,
                  treatment %in% c(removallevels, timinglevels))  %>%
    ggplot(aes(y =  counts, x = treatment, fill = treatment)) +
    geom_boxplot(lwd=0.5, outlier.size = 0.1) +
    theme_B3() +
    facet_wrap(~sex, nrow = 2, scales = "free_y") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_text(face = "italic")) +
    scale_fill_manual(values = allcolors)  +
    labs( y = whichgene,  x = "Treatment", title =" ")  +
    geom_signif(comparisons = list(c("inc.d17", "m.inc.d17"),
                                   c("inc.d17", "prolong")),
                map_signif_level=F, step_increase = 0.1,
                textsize = 1.5)
  return(p)
  
}

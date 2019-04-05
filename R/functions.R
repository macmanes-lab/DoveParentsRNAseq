# printplotcontrasts prints 
# 1) a summary of the DEGs for a given 2 way contrast
# 3) a MDS plot of the 2 way contrast

printplotcontrasts <- function(whichcontrast){
  cont <- whichcontrast
  print(summary(decideTestsDGE(
    glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
    adjust.method="fdr", p.value=0.01)))
  print(plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main=whichcontrast, frame.plot=F,
               ylim=c(-40,20)))
}

# this custom function returns a dataframe 
# with the pvalue,FDR, and LFC for all genes for a given contrast
# not currently used, but was useful for working out the plotVolcano function

returnDEGs <- function(whichcontrast){
  cont <- whichcontrast
  tt <- topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n= 14937)$table
  return(tt)
}


# plots volcano-ish plots of LFC versus FDR
plotVolcanos <- function(whichcontrast){
  cont <- whichcontrast
  tt <- topTags(glmTreat(fit, contrast=my.contrasts[,cont]),
                n = 14937)$table
  tt <- tt %>%
    mutate(DEGs = ifelse(tt$logFC > 1 & tt$FDR < 0.05, 
                         yes = "up", 
                         no = ifelse(tt$logFC < -1 & tt$FDR < 0.05, 
                                     yes = "down", 
                                     no = "NS")))
  tt$DEGs <- factor(tt$DEGs, levels = c("down", "NS","up"))

  myplot <- ggplot(tt, aes(x = logFC, y = FDR, color = DEGs)) +
    geom_point(size = 1.5, alpha = 0.8, na.rm = T) +
    scale_y_reverse(limits = c(1,0)) +
    scale_color_manual(values=c("down" = "#2166ac",
                                "NS" = "#bdbdbd", 
                                "up" = "#b2182b")) +
    xlim(-40, 40)  + 
    theme(legend.position = "bottom")
  print(myplot)
}

# this is an old version
totalDEGs <- function(contrastvector){
  res <- results(dds, contrast = c(contrastvector[1],contrastvector[2],contrastvector[3]), independentFiltering = T)
  sumpadj <- sum(res$padj < 0.1, na.rm = TRUE)
  print(sumpadj)
}

# this is a new version
numDEGs <- function(group1, group2){
  res <- results(dds, contrast = c("treatment", group1, group2), independentFiltering = T)
  sumpadj <- sum(res$padj < 0.1, na.rm = TRUE)
  print(sumpadj)
}


## I've taken the pca function from DESeq2 and elaborated it so that I could extract up to 6 PCs

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



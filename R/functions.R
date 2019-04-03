# printplotcontrasts prints 
# 1) a summary of the DEGs for a given 2 way contrast
# 3) a MDS plot of the 2 way contrast

printplotcontrasts <- function(whichcontrast){
  cont <- whichcontrast
  print(summary(decideTestsDGE(
    glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
    adjust.method="fdr", p.value=0.01)))
  print(plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main=whichcontrast, frame.plot=F))
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
    scale_y_reverse() +
    scale_color_manual(values=c("down" = "#2166ac",
                                "NS" = "#bdbdbd", 
                                "up" = "#b2182b")) +
    xlim(-40, 40)  +
    theme(legend.position = "bottom")
  print(myplot)
}

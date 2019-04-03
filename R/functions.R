# printplotcontrasts prints 
# 1) a summary of the DEGs for a given 2 way contrast
# 2) results of the top 5 DEGs
# 3) a volcano plot of the 2 way contrast

printplotcontrasts <- function(whichcontrast){
  cont <- whichcontrast
  print(summary(decideTestsDGE(
    glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
    adjust.method="fdr", p.value=0.01)))
  print(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)
  print(plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main=whichcontrast, frame.plot=F))
}
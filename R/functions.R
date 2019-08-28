# subset col for heatmap
subsetcolData <- function(colData, eachgroup){
  
  # subset to look within one tissue in one sex
  colData <- colData %>%
    dplyr::filter(sextissue == eachgroup) %>%
    droplevels()
  row.names(colData) <- colData$V1
  return(colData)
}

subsetcolData2 <- function(colData, eachgroup){
  
  # subset to look within one tissue in one sex
  colData <- colData %>%
    dplyr::filter(sextissue %in% eachgroup) %>%
    droplevels()
  row.names(colData) <- colData$V1
  return(colData)
}




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


subsetDESeq2 <- function(colData, countData, eachgroup){
  
  # subset to look within one tissue in one sex
  colData <- colData %>%
    dplyr::filter(sextissue %in% eachgroup) %>%
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
                                design = ~ treatment * sex )
  
  print(dds)
  dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
  print(dim(dds))
  
  dds <- DESeq(dds, parallel = TRUE) # Differential expression analysis
  return(dds)
}


# print total number of differntially expressed genes
# numDEGs('m.inc.d3', 'm.inc.d9')


numDEGs <- function(dds, group1, group2){
  res <- results(dds, contrast = c("treatment", group1, group2), independentFiltering = T)
  sumpadj <- sum(res$padj < 0.01, na.rm = TRUE)
  return(sumpadj)
}

# resturn pvalues for all genes

returnpadj <- function(group1, group2){
  res <- results(dds, contrast = c("treatment", group1, group2), independentFiltering = T)
  pvals <- as.data.frame(res$padj)
  padjcolname <- as.character(paste("padj", group1, group2, sep=""))
  colnames(pvals) <- c(padjcolname)
  return(pvals)
}


################## ALL DEG comparisons

returntotalDEGs <- function(dds){
  
  colData <- c.colData # %>%

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
  totalDEGS$V2 <- factor(totalDEGS$V2, levels =  c("control", "bldg", "lay",
                                                   "inc.d3", "inc.d9", "inc.d17",
                                                   "hatch", "n5", "n9"))
  
  totalDEGS$V1 <- factor(totalDEGS$V1, levels =  c("control", "bldg", "lay",
                                                   "inc.d3", "inc.d9", "inc.d17",
                                                   "hatch", "n5", "n9"))

  totalDEGS <- totalDEGS %>% dplyr::na_if(0)
  
  print(str(totalDEGS))
  
  
  allcontrasts <- totalDEGS %>%
    ggplot( aes(V1, V2)) +
    geom_tile(aes(fill = V3)) +
     theme_minimal(base_size = 12) + 
    geom_text(aes(label = round(V3, 1)), color = "black")+
    scale_fill_viridis(na.value="#bdbdbd", 
                      # limits = c(0,7000),
                       option = "C") +
    xlab(NULL) + ylab(NULL) +
    labs(fill = "# of DEGs",
         title = mysubtitle, caption = "  ") +
    theme(axis.text.x = element_text(angle = 90)) +
    coord_flip()
  print(totalDEGS)
  plot(allcontrasts)
}

################## characterization comparisons only

# bar plot with subset of characterization

serialtimepoints <- c("control.bldg" , "bldg.lay", "lay.inc.d3", "inc.d3.inc.d9", 
                      "inc.d9.inc.d17", "inc.d17.hatch", "hatch.n5", "n5.n9")

plotserialDEGs <- function(DEGs, mysubtitle, myfill){
  
  # subset to look within one tissue in one sex
  DEGs <- DEGs %>%
    dplyr::mutate(comparison = row.names(.)) %>%
    dplyr::filter(comparison %in% serialtimepoints) 
  
  DEGs$comparison <- factor(DEGs$comparison, levels = serialtimepoints)
  
  mybarplot <- ggplot(DEGs, aes(comparison)) +
    geom_bar(aes(weight = V3), fill = myfill) +
     theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()) +
    labs(subtitle = mysubtitle, y = "Number of DEGs", x = NULL) +
    ylim(0, 5105) +
    geom_line(aes(x=comparison, y=V3, group = 1))

  return(mybarplot)
}


## characterization line graphs


subsetDEGs <- function(DEGs, groupname){
  DEGs <- DEGs %>%
    dplyr::mutate(comparison = row.names(.),
                  sextissue = groupname) %>%
    dplyr::filter(comparison %in% serialtimepoints) 
  DEGs$comparison <- factor(DEGs$comparison, levels = serialtimepoints)
  return(DEGs)
}
  
  

########################### manipulations comparisons

manipVchar <- c("inc.d3.m.inc.d3",
                "inc.d9.m.inc.d8", "inc.d9.m.inc.d9",
                "inc.d17.m.inc.d17",
                "hatch.prolong", "hatch.extend", "hatch.m.n2")

plotmanipDEGs <- function(DEGs, mysubtitle, legendornot){
  
  # subset to look within one tissue in one sex
  DEGs <- DEGs %>%
    dplyr::mutate(comparison = row.names(.)) %>%
    dplyr::filter(comparison %in% manipVchar) 
  
  DEGs$description <- ifelse(grepl("inc.d9.m.inc.d8", DEGs$comparison),"early chicks", 
                             ifelse(grepl("d3|d17|d9|n2", DEGs$comparison),"remove offspring",
                                           ifelse(grepl("prolong", DEGs$comparison),"prolong inc.",
                                                  ifelse(grepl("extend", DEGs$comparison),"extend inc.", NA))))
  
  DEGs$comparison <- factor(DEGs$comparison, levels = manipVchar)
  DEGs$description <- factor(DEGs$description, levels = c("remove offspring", "early chicks", "prolong inc.", "extend inc."))
  
  mybarplot <- ggplot(DEGs, aes(comparison)) +
    geom_bar(aes(weight = V3, fill = description)) +
     theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(subtitle = mysubtitle, y = "Number of DEGs", x = NULL) +
    #ylim(0, 5105) +
    #geom_line(aes(x=comparison, y=V3, group = 1)) +
    scale_x_discrete(labels=c("Day 3\nRemove eggs", "Day 9\nAdd chicks", "Day 9\nRemove eggs",
                              "Day 17\nRemove eggs", "~Hatch\nProlong", "~Hatch\nExtend", "~Hatch\nRemove chicks")) +
    theme(legend.position = legendornot,
          legend.title = element_blank(),
          panel.grid = element_blank()) +
    scale_fill_manual(values = c("#FF0000", "#00A08A", "#F2AD00", "#5BBCD6"))
  
  return(mybarplot)
}



######### offspring removal


offspringremoval <- c("lay.inc.d3", "inc.d3.m.inc.d3", 
                      "inc.d3.inc.d9", "inc.d9.m.inc.d9",
                      "inc.d9.inc.d17", "inc.d17.m.inc.d17", 
                      "hatch.n5", "hatch.m.n2")

plotremoval <- function(DEGs, mysubtitle, legendornot){
  
  # subset to look within one tissue in one sex
  DEGs <- DEGs %>%
    dplyr::mutate(comparison = row.names(.)) %>%
    dplyr::filter(comparison %in% offspringremoval) 
  
  DEGs$description <- ifelse(grepl(".m.", DEGs$comparison),"offspring removal", "normal transition")
  
  DEGs$comparison <- factor(DEGs$comparison, levels = offspringremoval)

  mybarplot <- ggplot(DEGs, aes(comparison)) +
    geom_bar(aes(weight = V3, fill = description)) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(subtitle = mysubtitle, y = "Number of DEGs", x = NULL) +
    ylim(0, 2800) +
    theme(legend.position = legendornot,
          legend.title = element_blank(),
          panel.grid = element_blank()) 
  
  return(mybarplot)
}

######### prolong delay

prolongdelay <- c( "inc.d3.inc.d9", "inc.d9.inc.d17", 
                   "inc.d9.m.inc.d8",
                   "inc.d17.hatch", "hatch.n5",
                   "hatch.prolong", "hatch.extend")


plotprolongdelay <- function(DEGs, mysubtitle, legendornot){
  
  # subset to look within one tissue in one sex
  DEGs <- DEGs %>%
    dplyr::mutate(comparison = row.names(.)) %>%
    dplyr::filter(comparison %in% prolongdelay) 
  
  DEGs$description <- ifelse(grepl(".m.", DEGs$comparison),"chicks hatch early", 
                             ifelse(grepl("hatch.prolong", DEGs$comparison),"prolong inc",  
                                    ifelse(grepl("hatch.extend", DEGs$comparison),"extend hatch", "normal transition")))
  
  DEGs$comparison <- factor(DEGs$comparison, levels = prolongdelay)
  DEGs$description <- factor(DEGs$description, levels = c("normal transition", "chicks hatch early", "prolong inc", "extend hatch"))
  
  
  mybarplot <- ggplot(DEGs, aes(comparison)) +
    geom_bar(aes(weight = V3, fill = description)) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(subtitle = mysubtitle, y = "Number of DEGs", x = NULL) +
    ylim(0, 2800) +
    theme(legend.position = legendornot,
          legend.title = element_blank(),
          panel.grid = element_blank()) 
  
  return(mybarplot)
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
  pcadata <- pcadataframe(vsd, intgroup=c("treatment"), returnData=TRUE)
  percentVar <- round(100 * attr(pcadata, "percentVar"))
  print(percentVar)
  
  print(summary(aov(PC1 ~ treatment, data=pcadata)))
  print(TukeyHSD(aov(PC1 ~ treatment, data=pcadata), which = "treatment"))
  print(summary(aov(PC2 ~ treatment, data=pcadata))) 
  print(summary(aov(PC3 ~ treatment, data=pcadata))) 
  print(summary(aov(PC4 ~ treatment, data=pcadata))) 
  return(pcadata)
}


returnPCAs2 <- function(dds){
  
  dds <- dds
  
  vsd <- vst(dds, blind=FALSE) # variance stabilized 
  
  # create the dataframe using my function pcadataframe
  pcadata <- pcadataframe(vsd, intgroup=c("sex", "treatment"), returnData=TRUE)
  percentVar <- round(100 * attr(pcadata, "percentVar"))
  print(percentVar)
  
  print(summary(aov(PC1 ~ treatment * sex, data=pcadata)))
  print(summary(aov(PC2 ~ treatment * sex, data=pcadata))) 
  print(summary(aov(PC3 ~ treatment * sex, data=pcadata))) 
  print(summary(aov(PC4 ~ treatment * sex, data=pcadata))) 
  return(pcadata)
}


plotPC12 <- function(pcadata, mysubtitle){ 
  
  pcadata <- pcadata %>%
    mutate(hypothesis = fct_recode(treatment,
                            "anticipation" = "control",
                            "anticipation" = "bldg",
                            "incubation" = "lay",
                            "incubation" = "inc.d3",
                            "incubation" = "inc.d9",
                            "incubation" = "inc.d17",
                            "hatchling.care" = "hatch",
                            "hatchling.care" = "n5",
                            "hatchling.care" = "n9"))
  
  pca12 <- ggplot(pcadata, aes(PC1, PC2, color = treatment, shape = sex)) + 
    geom_point( size = 3) +
    #stat_ellipse() +
    #xlab(paste0("PC1")) +
    #ylab(paste0("PC2")) +
    labs(subtitle = mysubtitle) +
    #theme(legend.position = "bottom") +
    #guides( fill = guide_legend(order = 2, ncol=2)) +
    theme(legend.text = element_text(size=10)) 
  pca12
}  


## plot candidate genes 
# e.g. plotcandidates(dds.female_hypothalamus, "female hypothalamus")

plotcandidates <- function(mydds, colData, mysubtitle){
  
  mydds <- mydds
  
  vsd <- vst(mydds, blind=FALSE) # variance stabilized 
  
  # make dataframe counts
  DEGs <- assay(vsd)
  DEGs <- as.data.frame(DEGs)
  DEGs$entrezid <- row.names(DEGs)

  # make dataframe with geneids and names and counts
  # how to gather: https://tidyr.tidyverse.org/reference/gather.html
 
  candidates <- full_join(geneinfo, DEGs, by = "entrezid")
  head(candidates)
  
  candidates <- candidates %>%
    filter(Name %in% c("AVP", "AVPR1A", "AVPR1B", "AVPR2",
                      "SERPINA4", "CRH", "CRHR1", "DRD1", "DRD2",
                      "DRD3",  "DRD4",  "DRD5",  "GABRQ",  "NR3C1",  "HSD11B1a", 
                      "HSD11B1b", "HSD11B2L",  "HSD11B1L",  "MC2R", "NR3C2",
                      "OXTR", "POMC",  "AR", "CYP19A1", "ESR1",  "ESR2", "FSHR", 
                      "FSHB", "GHRL", "GAL", "NPVF",  "NPFFR1", "GNRH1",  "GNRHR", 
                      "LEPR",  "LHCGR", "PGR", "PRL", "PRLR", "VIP", "VIPR1"))
  row.names(candidates) <- candidates$Name
  candidates <- candidates %>% dplyr::select(-row.names, -entrezid, -Name, -geneid)
  candidates <- candidates %>% drop_na()
  candidates <- as.data.frame(t(candidates))
  candidates$RNAseqID <- rownames(candidates)
  candidates <- candidates %>% gather(gene, value, -RNAseqID)  %>% 
    filter(RNAseqID != "gene")
  candidates$value <- as.numeric(candidates$value)
  candidates$V1  <- candidates$RNAseqID
  
  candidatecounts <- left_join(candidates, colData, by = "V1")
  candidatecounts$faketime <- as.numeric(candidatecounts$treatment)
  candidatecounts$gene <- as.factor(candidatecounts$gene)
  
  p1 <- ggplot(candidatecounts, aes(x = treatment, y = value, fill = treatment)) +
    geom_boxplot() +
    facet_wrap(~gene, scales = "free") +
     theme_minimal(base_size = 8) +
    theme(axis.text.x = element_blank(),
          legend.position = "bottom") +
    labs(x = NULL, subtitle = mysubtitle) +
    guides(fill = guide_legend(nrow = 1))
  return(p1)
  
}


plotprolactin <- function(mydds, colData, mysubtitle){
  
  mydds <- mydds
  
  vsd <- vst(mydds, blind=FALSE) # variance stabilized 
  
  # make dataframe counts
  DEGs <- assay(vsd)
  DEGs <- as.data.frame(DEGs)
  DEGs$entrezid <- row.names(DEGs)
  
  # make dataframe with geneids and names and counts
  # how to gather: https://tidyr.tidyverse.org/reference/gather.html
  
  candidates <- full_join(geneinfo, DEGs, by = "entrezid")
  head(candidates)
  
  candidates <- candidates %>%
    filter(Name %in% c( "PRL"))
  row.names(candidates) <- candidates$Name
  candidates <- candidates %>% dplyr::select(-row.names, -entrezid, -Name, -geneid)
  candidates <- candidates %>% drop_na()
  candidates <- as.data.frame(t(candidates))
  candidates$RNAseqID <- rownames(candidates)
  candidates <- candidates %>% gather(gene, value, -RNAseqID)  %>% 
    filter(RNAseqID != "gene")
  candidates$value <- as.numeric(candidates$value)
  candidates$V1  <- candidates$RNAseqID
  
  candidatecounts <- left_join(candidates, colData, by = "V1")
  candidatecounts$faketime <- as.numeric(candidatecounts$treatment)
  candidatecounts$gene <- as.factor(candidatecounts$gene)
  
  p1 <- ggplot(candidatecounts, aes(x = treatment, y = value, fill = treatment)) +
    geom_boxplot() +
    facet_wrap(~gene, scales = "free") +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_blank(),
          legend.position = "bottom") +
    labs(x = NULL, subtitle = mysubtitle) +
    guides(fill = guide_legend(nrow = 1))
  return(p1)
  
}





## make pheatmaps 
# e.g. makepheatmap(dds.female_hypothalamus, "female hypothalamus")

makepheatmap <- function(mydds, colData, mysubtitle){

  vsd <- vst(mydds, blind=FALSE) # variance stabilized 
  
  # make dataframe counts
  DEGs <- assay(vsd)
  
  DEGs <- as.data.frame(DEGs)
  DEGs$entrezid <- row.names(DEGs)
  
  candidates <- full_join(geneinfo, DEGs, by = "entrezid")
  
  candidates <- candidates %>%
    filter(Name %in% c("AVP", "AVPR1A", "AVPR1B", "AVPR2",
                       "SERPINA4", "CRH", "CRHR1", "DRD1", "DRD2",
                       "DRD3",  "DRD4",  "DRD5",  "GABRQ",  "NR3C1",  "HSD11B1a", 
                       "HSD11B1b", "HSD11B2L",  "HSD11B1L",  "MC2R", "NR3C2",
                       "OXTR", "POMC",  "AR", "CYP19A1", "ESR1",  "ESR2", "FSHR", 
                       "FSHB", "GHRL", "GAL", "NPVF",  "NPFFR1", "GNRH1",  "GNRHR", 
                       "LEPR",  "LHCGR", "PGR", "PRL", "PRLR", "VIP", "VIPR1"))
  
  row.names(candidates) <- candidates$Name
  candidates <- candidates %>% dplyr::select(-row.names, -entrezid, -Name, -geneid)
  candidates <- candidates %>% drop_na()
  
  DEGsmatrix <- as.matrix(candidates)
  head(DEGsmatrix)
  
  DEGsmatrix <- DEGsmatrix - rowMeans(DEGsmatrix)

  paletteLength <- 30
  myBreaks <- c(seq(min(DEGsmatrix), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(DEGsmatrix)/paletteLength, max(DEGsmatrix), length.out=floor(paletteLength/2)))
  
  anndf <- colData %>% dplyr::select(treatment)
  rownames(anndf) <- colData$V1
  anncolors <- list(treatment = c("control" = "#F8766D",
                                          "bldg" = "#D39200" ,
                                          "lay"  = "#93AA00",
                                          "inc.d3" =  "#00BA38" ,
                                          "inc.d9" = "#00C19F",
                                          "inc.d17" = "#00B9E3",
                                          "hatch"  = "#619CFF",
                                          "n5"   = "#DB72FB",
                                          "n9" = "#FF61C3"))
  
  levels(c.colData$treatment)
  
  p1 <- pheatmap(DEGsmatrix, show_rownames = T, show_colnames = F,
           color = viridis(30),
           #breaks=myBreaks,
           annotation_col = anndf,
           annotation_colors = anncolors,
           main = mysubtitle)
  return(p1)
  }


## new correlation heatmap

plotcorrelationheatmaps <- function(mydds, mycoldata, mysubtitle){
  dds <- mydds
  vsd <- vst(dds, blind=FALSE) # variance stabilized 
  
  colnames(vsd) = mycoldata$treatment # set col names to group name
  
  vsdm <- assay(vsd) # create matrix
  
  vsdmmean <-sapply(unique(colnames(vsdm)), function(i)
    rowMeans(vsdm[,colnames(vsdm) == i]))
  
  myannotations <- mycoldata %>% 
    dplyr::distinct(lastday, penultimate, treatment) 
  rownames(myannotations) <- myannotations$treatment
  myannotations$treatment <- NULL
  
  mycor <- cor(vsdmmean)
  
  myBreaks = c(0.91,0.92,0.93,0.94,0.95, 0.96, 0.97,0.98, 0.99, 1.0)
  
  myBreaks <-seq(0.95, 1.0, length.out = 10)
  myBreaks
  
  pheatmap(cor(vsdmmean),
           annotation_row = myannotations,
           annotation_col = myannotations,
           annotation_colors = myannotationscolors,
           annotation_names_row = F,
           main= mysubtitle,
           color = inferno(10),
           show_rowname= F, show_colnames = F,
           breaks = myBreaks
  )
}


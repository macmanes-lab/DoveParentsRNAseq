######## subsetcolData #############

# subset colData to look one tissue in one sex
subsetcolData <- function(colData, eachgroup){
  
  colData <- colData %>%
    dplyr::filter(sextissue == eachgroup) %>%
    droplevels()
  row.names(colData) <- colData$V1
  return(colData)
}

######### subsetcolData2 ############

# subset to look within one tissue in two sexes
subsetcolData2 <- function(colData, eachgroup){
  
  colData <- colData %>%
    dplyr::filter(sextissue %in% eachgroup) %>%
    droplevels()
  row.names(colData) <- colData$V1
  return(colData)
}


############ subsetDESeq #########

# run DESeq on subset of data - treatment only
subsetDESeq <- function(colData, countData, eachgroup){
  
  colData <- subsetcolData(colData, eachgroup)
  
  # save counts that match colData
  savecols <- as.character(colData$V1) 
  savecols <- as.vector(savecols) 
  countData <- countData %>% dplyr::select(one_of(savecols)) 
  
  # assert that row and col lenghts are equal
  stopifnot(ncol(countData) == nrow(colData))
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ treatment )
  
  print(dds)
  dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
  print(dim(dds))
  
  dds <- DESeq(dds) # Differential expression analysis
  return(dds)
}


############ subsetDESeq2 #########

# run DESeq on subset of data sex * treatment 
subsetDESeq2 <- function(colData, countData, eachgroup){
  
  # subset to look within one tissue in one sex
  colData <- subsetcolData2(colData, eachgroup)
  
  # save counts that match colData
  savecols <- as.character(colData$V1) 
  savecols <- as.vector(savecols) 
  countData <- countData %>% dplyr::select(one_of(savecols)) 
  
  # check that row and col lenghts are equal
  stopifnot(ncol(countData) == nrow(colData))
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~ treatment * sex )
  
  print(dds)
  dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
  print(dim(dds))
  
  dds <- DESeq(dds, parallel = TRUE) # Differential expression analysis
  return(dds)
}

############ numDEGs #########

# print total number of differntially expressed genes
numDEGs <- function(dds, group1, group2){
  res <- results(dds, contrast = c("treatment", group1, group2), independentFiltering = T)
  sumpadj <- sum(res$padj < 0.01, na.rm = TRUE)
  return(sumpadj)
}


############ returnpadj #########

# return pvalues for all genes
returnpadj <- function(group1, group2){
  res <- results(dds, contrast = c("treatment", group1, group2), independentFiltering = T)
  pvals <- as.data.frame(res$padj)
  padjcolname <- as.character(paste("padj", group1, group2, sep=""))
  colnames(pvals) <- c(padjcolname)
  return(pvals)
}


############returntotalDEGs ############ 

returntotalDEGs <- function(dds){
  
  colData <- c.colData # %>%
  
  a <- levels(colData$treatment)
  b <- levels(colData$treatment)
  
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

############ plottotalDEGs ############ 

plottotalDEGs <- function(myDEGS, mysubtitle){  
  
  totalDEGS <- myDEGS
  
  totalDEGS$V2 <- factor(totalDEGS$V2, levels =  c("control", "bldg", "lay",
                                                   "inc.d3", "inc.d9", "inc.d17",
                                                   "hatch", "n5", "n9"))
  
  totalDEGS$V1 <- factor(totalDEGS$V1, levels =  c("control", "bldg", "lay",
                                                   "inc.d3", "inc.d9", "inc.d17",
                                                   "hatch", "n5", "n9"))

  #totalDEGS <- totalDEGS %>% dplyr::na_if(0)
  
  allcontrasts <- totalDEGS %>%
    ggplot( aes(V1, V2)) +
    geom_tile(aes(fill = V3)) +
     theme_minimal(base_size = 12) + 
    geom_text(aes(label = round(V3, 1)), color = "grey")+
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

############ plotserialDEGs ############ 

# make bar plot with subset of characterization
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


############ subsetDEGs ############ 

# used for making line graphs in 03_DESeq2_characterization.Rmd
subsetDEGs <- function(DEGs, groupname){
  DEGs <- DEGs %>%
    dplyr::mutate(comparison = row.names(.),
                  sextissue = groupname) %>%
    dplyr::filter(comparison %in% serialtimepoints) 
  DEGs$comparison <- factor(DEGs$comparison, levels = serialtimepoints)
  return(DEGs)
}
  


############ plotmanipDEGs ###############

# list of things to compare
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



######### plotremoval ######### 

# list of things to compare
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

######### plotprolongdelay ######### 

# list of things to compare
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

######### pcadataframe ######### 

# modified from DESeq2 and elaborated it so that I could extract up to 6 PCs

pcadataframe <- function(object, intgroup, ntop = 500, returnData = FALSE){
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
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], 
                  PC4 = pca$x[, 4],PC5 = pca$x[, 5],PC6 = pca$x[, 6],
                  group = group, intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:6]
    return(d)
  }
}

######### returnPCAs ######### 

# returns the stats for PC 1-4 for a single variable
returnPCAs <- function(vsd){
  
  # create the dataframe using my function pcadataframe
  pcadata <- pcadataframe(vsd, intgroup=c("treatment"), returnData=TRUE)
  percentVar <- round(100 * attr(pcadata, "percentVar"))
  
  print("Percent variance explained by PC 1-6")
  print(percentVar)
  
  print("PC1 ~ treatment, data=pcadata")
  print(summary(aov(PC1 ~ treatment, data=pcadata)))
  print(TukeyHSD(aov(PC1 ~ treatment, data=pcadata), which = "treatment"))
  
  print("PC2 ~ treatment, data=pcadata")
  print(summary(aov(PC2 ~ treatment, data=pcadata))) 
  
  print("PC3 ~ treatment, data=pcadata")
  print(summary(aov(PC3 ~ treatment, data=pcadata))) 
  
  print("PC4 ~ treatment, data=pcadata")
  print(summary(aov(PC4 ~ treatment, data=pcadata))) 
  
  print("PC5 ~ treatment, data=pcadata")
  print(summary(aov(PC5 ~ treatment, data=pcadata))) 
  
  print("PC6 ~ treatment, data=pcadata")
  print(summary(aov(PC6 ~ treatment, data=pcadata))) 
  return(pcadata)
}

######### returnPCA2 ######### 

# returns the stats for PC 1-4 for a two variables 

returnPCAs2 <- function(vsd){
  
  # create the dataframe using my function pcadataframe
  pcadata <- pcadataframe(vsd, intgroup=c("sex", "treatment"), returnData=TRUE)
  percentVar <- round(100 * attr(pcadata, "percentVar"))
  
  print("Percent variance explained by PC 1-6")
  print(percentVar)
  
  print("PC1 ~ treatment * sex, data=pcadata")
  print(summary(aov(PC1 ~ treatment * sex, data=pcadata)))
  
  print("PC2 ~ treatment * sex, data=pcadata")
  print(summary(aov(PC2 ~ treatment * sex, data=pcadata))) 
  
  print("PC3 ~ treatment * sex, data=pcadata")
  print(summary(aov(PC3 ~ treatment * sex, data=pcadata)))
  
  print("PC4 ~ treatment * sex, data=pcadata")
  print(summary(aov(PC4 ~ treatment * sex, data=pcadata)))
  
  
  print("PC5 ~ treatment * sex, data=pcadata")
  print(summary(aov(PC5 ~ treatment * sex, data=pcadata)))
  
  print("PC6 ~ treatment * sex, data=pcadata")
  print(summary(aov(PC6 ~ treatment * sex, data=pcadata))) 
  return(pcadata)
}

######### plotPC12 ######### 

# plot pc 1 and 2 for treatment and sex

plotPC12 <- function(pcadata, mysubtitle){ 

  pca12 <- ggplot(pcadata, aes(PC1, PC2, color = treatment, shape = sex)) + 
    geom_point() +
    #stat_ellipse() +
    #xlab(paste0("PC1")) +
    #ylab(paste0("PC2")) +
    labs(subtitle = mysubtitle) +
    #theme(legend.position = "bottom") +
    #guides( fill = guide_legend(order = 2, ncol=2)) +
    theme(legend.text = element_text(size=10)) 
  return(pca12)
}  


######### vsd.dataframe and corresponding col data ######### 

vsd.dataframe <- function(vsd){
  
  # extract varience stabilized and set rownames
  vsd.df <- assay(vsd)
  vsd.df <- as.data.frame(vsd.df)
  return(vsd.df)
}  



savevsdfiles <- function(myvsddf, mycolData, mytissue){
  myvsddf <- myvsddf
  mycolData <- mycolData
  
  vsdfilename <- paste0("../results/04_vsd_", mytissue, ".csv", sep = "" )
  colDatafilename <- paste0("../results/04_colData_", mytissue, ".csv", sep = "" )
  
  write.csv(myvsddf, vsdfilename, row.names = T)
  write.csv(mycolData, colDatafilename, row.names = F)
}


readvsd <- function(filename){
  vsd <- read_csv(filename)
  vsd <- as.data.frame(vsd)
  row.names(vsd) <- vsd$X1
  vsd$entrezid <- vsd$X1
  vsd$X1 <- NULL
  return(vsd)
}

readcolData <- function(filename){
  colData <- read_csv(filename)
  colData <- as.data.frame(colData)
  
  colData$treatment <- factor(colData$treatment, levels = 
                                  c("control",  "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", "hatch", "n5", "n9"))
  
  row.names(colData) <- colData$V1
  colData$sample <- colData$V1
  colData$V1 <- NULL
  return(colData)
}


### select candidate gene vsds
selectcandidatevsds <- function(listofgenes, vsd, colData){
  
  print(listofgenes)
  
  candidateentrezids <- geneinfo %>% 
    filter(Name %in% listofgenes) %>%  dplyr::select(entrezid) 
  
  print(candidateentrezids$entrezid)
  
  rowstofilter <- as.list(candidateentrezids$entrezid)
  
  candidatevsd <- vsd %>% 
    filter(entrezid %in% rowstofilter)
  candidatevsd <- as.data.frame(candidatevsd)
  row.names(candidatevsd) <- candidatevsd$entrezid
  
  candidatevsd <- right_join(geneinfo,candidatevsd) 
  candidatevsd <- as.data.frame(candidatevsd)
  row.names(candidatevsd) <- candidatevsd$Name
  candidatevsd <- candidatevsd[-c(1:3)]
  
  candidatevsd <- as.data.frame(t(candidatevsd))
  candidatevsd$sample <- row.names(candidatevsd)
  
  candidatevsd <- left_join(colData,candidatevsd)
  
  return(candidatevsd)
}



######### plotcandidates ######### 

plotcandidates <- function(vsd.df, colData, mysubtitle){
  
  # make dataframe with geneids and names and counts
  # how to gather: https://tidyr.tidyverse.org/reference/gather.html
 
  candidates <- full_join(geneinfo, vsd.df, by = "entrezid")
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

######### plotcandidates ######### 

plotprolactin <- function(vsd.df, colData, mysubtitle){

  candidates <- full_join(geneinfo, vsd.df, by = "entrezid")
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


###### plot wgcna candidates

plotWGCNAcandidates <- function(vsd, mygenelist, colData, mysubtitle){
  
  # make dataframe with geneids and names and counts
  # how to gather: https://tidyr.tidyverse.org/reference/gather.html
  
  vsd.df <- as.data.frame(assay(vsd))
  vsd.df$entrezid <- row.names(vsd.df)
  
  candidates <- full_join(geneinfo, vsd.df, by = "entrezid")
 
  # select genes of interst
  
  candidates <- candidates %>%
    filter(entrezid %in% mygenelist) %>% droplevels()
  
  candidates <- candidates %>% dplyr::select(-row.names, -entrezid, -geneid)
  
  candidates_long <- candidates %>% gather(-Name, key = "sample", value = "value")
  
  candidates_long$V1 <- candidates_long$sample
  
  candidatecounts <- left_join(candidates_long, colData, by = "V1")
  
  candidatecounts$Name <- as.factor(candidatecounts$Name)
  candidatecounts$treatment <- factor(candidatecounts$treatment, levels = charlevels)
  
  head(candidatecounts)
  
  bysextreatment <- group_by(candidatecounts, sex, treatment, Name)
  bysextreatment
  candidateST <- summarize(bysextreatment, expression = mean(value))
  
  
  p1 <- ggplot(candidateST, aes(x = as.numeric(treatment), y = expression, color = sex)) +
    geom_point() +
    geom_smooth(se = FALSE) +
    facet_wrap(~Name, scales = "free_y") +
    theme_rmh() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 60,  hjust=1),
          axis.text.y = element_blank()) +
    guides(fill = guide_legend(nrow = 1)) +
    labs(x = NULL, y = "Gene expression",
         subtitle = mysubtitle) +
    scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9),
                       labels=charlevels)
  return(p1)
  
}

plotWGCNAcandidatesManip <- function(vsd, mygenelist, colData, mysubtitle){
  
  # make dataframe with geneids and names and counts
  # how to gather: https://tidyr.tidyverse.org/reference/gather.html
  
  vsd.df <- as.data.frame(assay(vsd))
  vsd.df$entrezid <- row.names(vsd.df)
  
  candidates <- full_join(geneinfo, vsd.df, by = "entrezid")
  
  # select genes of interst
  
  candidates <- candidates %>%
    filter(entrezid %in% mygenelist) %>% droplevels()
  
  candidates <- candidates %>% dplyr::select(-row.names, -entrezid, -geneid)
  
  candidates_long <- candidates %>% gather(-Name, key = "sample", value = "value")
  
  candidates_long$V1 <- candidates_long$sample
  
  candidatecounts <- left_join(candidates_long, colData, by = "V1")
  
  candidatecounts$Name <- as.factor(candidatecounts$Name)
  candidatecounts$treatment <- factor(candidatecounts$treatment, levels = maniplevels)
  
  head(candidatecounts)
  
  bysextreatment <- group_by(candidatecounts, sex, treatment, Name)
  bysextreatment
  candidateST <- summarize(bysextreatment, expression = mean(value))
  
  
  p1 <- ggplot(candidateST, aes(x = as.numeric(treatment), y = expression, color = sex)) +
    geom_point() +
    geom_smooth(se = FALSE) +
    facet_wrap(~Name, scales = "free_y") +
    theme_rmh() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 60,  hjust=1),
          axis.text.y = element_blank()) +
    guides(fill = guide_legend(nrow = 1)) +
    labs(x = NULL, y = "Gene expression",
         subtitle = mysubtitle) +
    scale_x_continuous(breaks=c(1,2,3,4,5,6,7),
                       labels=maniplevels)
  return(p1)
  
}

######### makepheatmap ######### 

# makes a candidate heat map!!
makepheatmap <- function(vsd.df, colData, mysubtitle){

  candidates <- full_join(geneinfo, vsd.df, by = "entrezid")
  
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
  
  p1 <- pheatmap(DEGsmatrix, show_rownames = T, show_colnames = F,
           color = viridis(30),
           #breaks=myBreaks,
           annotation_col = anndf,
           annotation_colors = anncolors,
           main = mysubtitle)
  return(p1)
  }


######### plotcorrelationheatmaps ######### 

plotcorrelationheatmaps <- function(vsd, mycoldata, mysubtitle){

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



######### LDAdata.treatment ######### 

LDAdata.treatment <- function(vsd, mycolData, mypredictor){
  
  # data prep
  tvsd <- as.data.frame(t(assay(vsd)))   # get vsd counts and transform
  
  tvsd <- tvsd[, sample(ncol(tvsd), 20)] # select 20 random genes for testing
  tvsd$V1 <- row.names(tvsd)             # prep for joining
  
  mydata <- left_join(tvsd, mycolData, by = "V1")  # merge counts and variables
  mydata$V1 <- NULL
  mydata$sextissue <- NULL
  mydata$hypothesis <- NULL
  
  # Split the data into training (80%) and test set (20%)
  set.seed(175)
  
  training.samples <- mydata$treatment %>%
    createDataPartition(p = 0.8, list = FALSE)
  train.data <- mydata[training.samples, ]
  test.data <- mydata[-training.samples, ]
  
  # Normalize the data. Categorical variables are automatically ignored.
  # Estimate preprocessing parameters
  preproc.param <- train.data %>% 
    preProcess(method = c("center", "scale"))
  
  # Transform the data using the estimated parameters
  train.transformed <- preproc.param %>% predict(train.data)
  test.transformed <- preproc.param %>% predict(test.data)
  
  # LDA analysis
  # Fit the model
  model <- lda(treatment~ ., data = train.transformed)
  # Make predictions
  predictions <- model %>% predict(test.transformed)
  
  # Model accuracy
  print("model accuracy")
  print("predictions$class==test.transformed$treatment)")
  print(mean(predictions$class==test.transformed$treatment))
  
  # results
  print("the samples sizes")
  print(model$counts)
  
  print("the prior probabilities used")
  print(model$prior)
  
  print("svd: the singular values, which give the ratio of the between- and within-group standard deviations on the linear discriminant variables. Their squares are the canonical F-statistics.")
  print(model$svd)
  
  
  #  predictions
  predictions <- model %>% predict(test.transformed)
  
  # Predicted classes
  #print(predictions$class, 6)
  # Predicted probabilities of class memebership.
  #print(predictions$posterior, 6) 
  # Linear discriminants
  #print(predictions$x, 3)
  
  lda.data <- cbind(train.transformed, predict(model)$x)
  return(lda.data)
}  

######### LDAplot.treatment ######### 

LDAplot.treatment <- function(LDAdata, mytitle, mysubtitle, myxlab, myylab){
  p <- ggplot(data = LDAdata, aes(LD1, LD2, color = treatment, shape = sex)) +
    geom_point(size = 1.5) +
    labs(title = mytitle,
         subtitle = mysubtitle,
         x = myxlab,
         y = myylab) +
    scale_color_manual(name = "parental stage",
                       values = c("control" = "#F8766D",
                                  "bldg" = "#D39200" ,
                                  "lay"  = "#93AA00",
                                  "inc.d3" =  "#00BA38" ,
                                  "inc.d9" = "#00C19F",
                                  "inc.d17" = "#00B9E3",
                                  "hatch"  = "#619CFF",
                                  "n5"   = "#DB72FB",
                                  "n9" = "#FF61C3"),
                       labels = c("control" = "<span style='color:#F8766D'>control</span>",
                                  "bldg" = "<span style='color:#D39200'>nest bulding</span>",
                                  "lay" = "<span style='color:#93AA00'>egg laid</span>",
                                  "inc.d3" = "<span style='color:#00BA38'>incubation day 3</span>", 
                                  "inc.d9" = "<span style='color:#00B9E3'>incubation day 9</span>",
                                  "inc.d17" = "<span style='color:#00B9E3'>incubation day 17</span>",
                                  "hatch" = "<span style='color:#619CFF'>chicks hatch</span>", 
                                  "n5" = "<span style='color:#DB72FB'>nestling care day 5</span>",
                                  "n9" = "<span style='color:#FF61C3'>nestling care day 9</span>")) +
    theme(legend.text = element_markdown(size = 8))
  plot(p)
} 


LDAplot.manipulation <- function(LDAdata, mytitle, mysubtitle, myxlab, myylab){
  p <- ggplot(data = LDAdata, aes(LD1, LD2, color = treatment, shape = sex)) +
    geom_point(size = 1.5) +
    labs(title = mytitle,
         subtitle = mysubtitle,
         x = myxlab,
         y = myylab) +
    scale_color_manual(name = "parental stage",
                       values = c("remove.d03" = "#F8766D",
                                  "remove.d09" = "#D39200" ,
                                  "remove.d17"  = "#93AA00",
                                  "remove.d20" =  "#00BA38" ,
                                  "extend" = "#00C19F",
                                  "prolong" = "#00B9E3",
                                  "early"  = "#619CFF"),                 
                       labels = c("remove.d03" = "<span style='color:#F8766D'>remove.d03</span>",
                                  "remove.d09" = "<span style='color:#D39200'>remove.d09</span>",
                                  "remove.d17" = "<span style='color:#93AA00'>remove.d17</span>",
                                  "remove.d20" = "<span style='color:#00BA38'>remove.d20</span>", 
                                  "extend" = "<span style='color:#00B9E3'>extend</span>",
                                  "prolong" = "<span style='color:#00B9E3'>prolong</span>",
                                  "early" = "<span style='color:#619CFF'>early</span>")) +
    theme(legend.text = element_markdown(size = 8))
  plot(p)
} 


######### LDAdata.hypothesis ######### 

LDAdata.hypothesis <- function(vsd, mycolData, mypredictor){
  
  # data prep
  tvsd <- as.data.frame(t(assay(vsd)))   # get vsd counts and transform
  
  tvsd <- tvsd[, sample(ncol(tvsd), 20)] # select 20 random genes for testing
  tvsd$V1 <- row.names(tvsd)             # prep for joining
  
  mydata <- left_join(tvsd, mycolData, by = "V1")  # merge counts and variables
  mydata$V1 <- NULL
  mydata$sextissue <- NULL
  mydata$treatment <- NULL
  
  # Split the data into training (80%) and test set (20%)
  set.seed(175)
  
  training.samples <- mydata$hypothesis %>%
    createDataPartition(p = 0.8, list = FALSE)
  train.data <- mydata[training.samples, ]
  test.data <- mydata[-training.samples, ]
  
  # Normalize the data. Categorical variables are automatically ignored.
  # Estimate preprocessing parameters
  preproc.param <- train.data %>% 
    preProcess(method = c("center", "scale"))
  
  # Transform the data using the estimated parameters
  train.transformed <- preproc.param %>% predict(train.data)
  test.transformed <- preproc.param %>% predict(test.data)
  
  # LDA analysis
  # Fit the model
  model <- lda(hypothesis~ ., data = train.transformed)
  # Make predictions
  predictions <- model %>% predict(test.transformed)
  
  # Model accuracy
  print("model accuracy")
  print("predictions$class==test.transformed$hypothesis)")
  print(mean(predictions$class==test.transformed$hypothesis))
  
  # results
  print("the samples sizes")
  print(model$counts)
  
  print("the prior probabilities used")
  print(model$prior)
  print(model$terms)
  
  print("svd: the singular values, which give the ratio of the between- and within-group standard deviations on the linear discriminant variables. Their squares are the canonical F-statistics.")
  print(model$svd)
  

  #  predictions
  predictions <- model %>% predict(test.transformed)
  
  # Predicted classes
  #print(predictions$class, 6)
  # Predicted probabilities of class memebership.
  #print(predictions$posterior, 6) 
  # Linear discriminants
  #print(predictions$x, 3)
  
  lda.data <- cbind(train.transformed, predict(model)$x)
  return(lda.data)
}  

######### LDAplot.hypothesis ######### 

LDAplot.hypothesis <- function(LDAdata, mytitle, mysubtitle, myxlab, myylab){
  p <- ggplot(data = LDAdata, aes(LD1, LD2, color = hypothesis, shape = sex)) +
    geom_point(size = 1.5) +
    labs(title = mytitle,
         subtitle = mysubtitle,
         x = myxlab,
         y = myylab) +
    scale_color_manual(name = "parental state",
                       values = c("non-parental" = "#1b9e77", 
                                  "incubation" = "#d95f02",
                                  "hatchling.care" = "#7570b3"),
                       labels = c("non-parental" = "<span style='color:#1b9e77'>non-parental</span>",
                                  "incubation" = "<span style='color:#d95f02'>incubation</span>", 
                                  "hatchling.care" = "<span style='color:#7570b3'>hatchling care</span>")) +
    theme(legend.text = element_markdown(size = 8))
  plot(p)
}  

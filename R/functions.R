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
  
  totalDEGS$V2 <- factor(totalDEGS$V2, levels =  charlevels)
  totalDEGS$V1 <- factor(totalDEGS$V1, levels =  charlevels)

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
  
  colData$treatment <- factor(colData$treatment, levels = charlevels)
  
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
  
  vsd.df <- as.data.frame(assay(vsd))
  vsd.df$entrezid <- row.names(vsd.df)
  
  candidates <- full_join(geneinfo, vsd.df, by = "entrezid")
  candidates <- candidates %>%
    filter(entrezid %in% mygenelist) %>% droplevels()
  candidates <- candidates %>% dplyr::select(-row.names, -entrezid, -geneid)
  candidates_long <- candidates %>% gather(-Name, key = "sample", value = "value")
  candidates_long$V1 <- candidates_long$sample
  candidatecounts <- left_join(candidates_long, colData, by = "V1")
  candidatecounts$Name <- as.factor(candidatecounts$Name)
  candidatecounts$treatment <- factor(candidatecounts$treatment, levels = charlevels)
  bysextreatment <- group_by(candidatecounts, sex, treatment, Name)
  bysextreatment
  candidateST <- summarize(bysextreatment, expression = mean(value))
  
  p1 <- ggplot(candidateST, aes(x = as.numeric(treatment), y = expression, color = sex)) +
    geom_point() +
    geom_smooth(se = FALSE) +
    facet_wrap(~Name, scales = "free_y") +
    mytheme() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 60,  hjust=1)) +
    guides(fill = guide_legend(nrow = 1)) +
    labs(x = NULL, y = "Gene expression",
         subtitle = mysubtitle) +
    scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9),
                       labels=charlevels) +
    scale_color_manual(values = colorstreatmentsex)
  return(p1)
  
}

plotWGCNAcandidatesManip <- function(vsd, mygenelist, colData, mysubtitle){
  
  vsd.df <- as.data.frame(assay(vsd))
  vsd.df$entrezid <- row.names(vsd.df)
  
  candidates <- full_join(geneinfo, vsd.df, by = "entrezid")
  candidates <- candidates %>%
    filter(entrezid %in% mygenelist) %>% droplevels()
  candidates <- candidates %>% dplyr::select(-row.names, -entrezid, -geneid)
  candidates_long <- candidates %>% gather(-Name, key = "sample", value = "value")
  candidates_long$V1 <- candidates_long$sample
  candidatecounts <- left_join(candidates_long, colData, by = "V1")
  candidatecounts$Name <- as.factor(candidatecounts$Name)
  candidatecounts$treatment <- factor(candidatecounts$treatment, levels = maniplevels)
  bysextreatment <- group_by(candidatecounts, sex, treatment, Name)
  bysextreatment
  candidateST <- summarize(bysextreatment, expression = mean(value))
  
  p1 <- ggplot(candidateST, aes(x = as.numeric(treatment), y = expression)) +
    geom_point(aes(color = treatment)) +
    geom_smooth(se = FALSE, aes(color = sex)) +
    facet_wrap(~Name, scales = "free_y") +
    mytheme +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 60,  hjust=1)) +
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



######### LDAplot.treatment ######### 

LDAplot.treatment <- function(LDAdata, mytitle, mysubtitle, myxlab, myylab){
  p <- ggplot(data = LDAdata, aes(LD1, LD2, color = treatment, shape = sex)) +
    geom_point(size = 1.5) +
    labs(title = mytitle,
         subtitle = mysubtitle,
         x = myxlab,
         y = myylab) 
  plot(p)
} 


## volcano plots

plot.volcano <- function(data, whichfactor, up, down, mycolors){
  
  numbersup <- data %>% dplyr::filter(direction == up) %>%  
    group_by(direction) %>% summarize(n = n()) %>% pull(n)
  numbersdown <- data %>% dplyr::filter(direction == down) %>%  
    group_by(direction) %>% summarize(n = n()) %>% pull(n)
  
  volcano <- data %>%
    ggplot(aes(x = lfc, y = logpadj)) + 
    geom_point(aes(color = direction, shape = tissue), size = 1, 
               alpha = 0.75, na.rm = T) + 
    theme_B3() +
    scale_color_manual(values = mycolors,
                       name = " ",
                       drop = FALSE,
                       breaks=c(down, "NS", up)) +
    ylim(c(0,25)) +  
    xlim(c(-8,8)) +
    labs(y = "-log10(p)", x = " ")  +
    theme(legend.position = "none",
          legend.direction = "horizontal",
          legend.spacing.x = unit(-0.1, 'cm'),
          legend.margin=margin(t=-0, r=0, b=0, l=0, unit="cm"),
          panel.grid = element_blank()) +
    scale_shape_manual(values = myshapes) +
    guides(shape = F) +
    annotate("text", label = numbersdown, x = -4, y = 25, size = 2) +
    annotate("text", label = numbersup, x = 4, y = 25, size = 2) 
  return(volcano)
}


# made for beta testing
# usaage
# tempmanip <- subsetDESeq3(colData, countData, c("female_pituitary","male_pituitary"), 
#                                                  c("m.inc.d3", "inc.d3", "inc.d9", "m.inc.d9"))

subsetDESeq3 <- function(colData, countData, eachgroup, eachtreatment){
  
  colData <- colData %>%
    dplyr::filter(sextissue %in% eachgroup,
                  treatment %in% eachtreatment) %>%
    droplevels()
  row.names(colData) <- colData$V1
  
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
  
  dds <- DESeq(dds, parallel = TRUE) # Differential expression analysis
  return(dds)
}




selectcBRCA1vsds <- function(listofgenes, vsd, colData){
  
  print(listofgenes)
  
  candidateentrezids <- geneinfo %>% 
    filter(entrezid %in% listofgenes) %>%  dplyr::select(entrezid) 
  
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





### modified corrr function

network_plot.cor_df <- function(rdf,
                                min_cor = .30,
                                legend = TRUE,
                                colours = c("indianred2", "white", "skyblue1"),
                                repel = TRUE,
                                curved = TRUE,
                                colors) {
  
  if (min_cor < 0 || min_cor > 1) {
    stop ("min_cor must be a value ranging from zero to one.")
  }
  
  if (!missing(colors))
    colours <- colors
  
  rdf <-  as_matrix(rdf, diagonal = 1)
  distance <- sign(rdf) * (1 - abs(rdf))
  
  # Use multidimensional Scaling to obtain x and y coordinates for points.
  points <- data.frame(stats::cmdscale(abs(distance)))
  colnames(points) <-  c("x", "y")
  points$id <- rownames(points)
  
  # Create a proximity matrix of the paths to be plotted.
  proximity <- abs(rdf)
  proximity[upper.tri(proximity)] <- NA
  diag(proximity) <- NA
  proximity[proximity < min_cor] <- NA
  
  # Produce a data frame of data needed for plotting the paths.
  n_paths <- sum(!is.na(proximity))
  paths <- data.frame(matrix(nrow = n_paths, ncol = 6)) 
  colnames(paths) <- c("x", "y", "xend", "yend", "proximity", "sign")
  path <- 1
  for(row in 1:nrow(proximity)) {
    for(col in 1:ncol(proximity)) {
      path_proximity <- proximity[row, col]
      if (!is.na(path_proximity)) {
        path_sign <- sign(distance[row, col])
        x    <- points$x[row]
        y    <- points$y[row]
        xend <- points$x[col]
        yend <- points$y[col]
        paths[path, ] <- c(x, y, xend, yend, path_proximity, path_sign)
        path <- path + 1
      }
    }
  }
  
  plot_ <- list(
    # For plotting paths
    if (curved) geom_curve(data = paths,
                           aes(x = x, y = y, xend = xend, yend = yend,
                               alpha = proximity, 
                               colour = proximity*sign),
                           size = 2), 
    if (!curved) geom_segment(data = paths,
                              aes(x = x, y = y, xend = xend, yend = yend,
                                  alpha = proximity, 
                                  colour = proximity*sign),
                              size = 2), 
    scale_alpha(limits = c(0, 1)),
    scale_size(limits = c(0, 1)),
    scale_colour_gradientn(limits = c(-1, 1), colors = colours),
    # Plot the points
    geom_point(data = points,
               aes(x, y),
               size = 1, shape = 19, colour = "black"),
    # Plot variable labels
    if (repel) ggrepel::geom_text_repel(data = points,
                                        aes(x, y, label = id),
                                        size = 3,
                                        segment.size = 0.0,
                                        segment.color = "black"),
    if (!repel) geom_text(data = points,
                          aes(x, y, label = id),
                          fontface = 'italic', size = 3),
    # expand the axes to add space for curves
    expand_limits(x = c(min(points$x) - .1,
                        max(points$x) + .1),
                  y = c(min(points$y) - .1,
                        max(points$y) + .1)
    ),
    # Theme and legends
    theme_void(),
    guides(size = "none", alpha = "none"),
    if (legend)  labs(colour = NULL),
    if (!legend) theme(legend.position = "none")
  )
  
  ggplot() + plot_
  
}

## PCA and FVIZ plots

subsetmakepca <- function(whichtissue, whichtreatment, whichsex){	
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
  
  mypca <- prcomp(countData)	
  mypcadf <- data.frame(PC1 = mypca$x[, 1], PC2 = mypca$x[, 2], PC3 = mypca$x[, 3], 	
                        PC4 = mypca$x[, 4],PC5 = mypca$x[, 5],PC6 = mypca$x[, 6],	
                        ID = row.names(countData))	
  mypcadf$V1 <- row.names(mypcadf)	
  mypcadf <- left_join(colData, mypcadf)	
  mypcadf <- mypcadf %>% dplyr::select(bird,sex,tissue,treatment,PC1:PC6)	
  return(mypcadf)	
}	


plotcolorfulpcs <- function(mypcadf,  whichfactor, whichcolors){	
  p <- mypcadf %>%	
    ggplot(aes(x = PC1, y = PC2, shape = tissue, color = whichfactor )) +	
    geom_point(size = 1)  +	
    theme_B3() +	
    theme(legend.title = element_blank(),	
          axis.text = element_blank(),	
          legend.position = "none") +	
    labs(x = "PC1", y = "PC2")  +	
    scale_color_manual(values = whichcolors) +	
    scale_shape_manual(values = myshapes)	+
    stat_ellipse( )
  
  return(p)	
}	


makefvizdf <-  function(whichtissue, whichtreatment, whichsex){	
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
  mypca <- prcomp(countData)	
  return(mypca)	
}	


plotfriz <- function(frizdf){
  
  p <- fviz_pca_var(frizdf,  
                    axes.linetype = "blank", 
                    repel = T , 
                    select.var= list(contrib = 3))  + 
    labs(title = NULL) + 
    theme_B3() 
  return(p)
}




## create DEGs

createDEGdfsave <- function(up, down, mytissue){
  
  res <- results(dds, contrast = c("treatment", up, down), independentFiltering = T, alpha = 0.1)
  
  DEGs <- data.frame(gene = row.names(res),
                        padj = res$padj, 
                        logpadj = -log10(res$padj),
                        lfc = res$log2FoldChange,
                        sextissue = mytissue)
  DEGs <- na.omit(DEGs)
  DEGs <- DEGs %>%
    dplyr::mutate(direction = ifelse(DEGs$lfc > 0 & DEGs$padj < 0.1, 
                                     yes = up, no = ifelse(DEGs$lfc < 0 & DEGs$padj < 0.1, 
                                                           yes = down, no = "NS"))) %>% 
    dplyr::arrange(desc(lfc)) 
  
  DEGs$direction <- factor(DEGs$direction, levels = c(down, "NS", up)) 
  
  # write DEGsframe of only significant genes
  DEGs <- DEGs %>% dplyr::filter(direction != "NS")
  print(str(DEGs))
  
  partialfilename = paste("_", down, "_", up, sep = "")
  myfilename = paste0("../results/DESeq2/", mytissue, partialfilename, "_DEGs.csv")
  
  write.csv(DEGs, myfilename, row.names = F)
  # return DEGs frome with all data, included NS genes
  print(head(DEGs))
}  


## prolactin plots 

plotprolactin <- function(df, myy, myylab, mysubtitle ){
  
  p <-  ggplot(df, aes(x = treatment, y = myy)) +
    geom_boxplot(aes(fill = treatment, color = sex)) +
    theme_B3() +
    scale_fill_manual(values = allcolors) +
    scale_color_manual(values = sexcolors) +
    labs(y = "prolactin (ng/mL)", x = NULL) +
    theme(legend.position = c(0.85,0.15), legend.direction = "horizontal") + 
    labs(x = "parental stage", subtitle = mysubtitle, y= myylab) +
    guides(fill = guide_legend(nrow = 1)) 
  return(p)
}


## tsne  figure 1 and 5

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
  
  tsne_model <- Rtsne(euclidist, check_duplicates=FALSE, pca=TRUE, perplexity=10, theta=0.5, dims=2)
  tsne_df = as.data.frame(tsne_model$Y) 
  
  # prep for adding columns
  colData2 <- colData 
  colData2$V1 <- NULL
  tsne_df_cols <- cbind(colData2, tsne_df)
  return(tsne_df_cols)
}

plottsneelipse <- function(tsnedf, pointcolor, whichcolors){
  p <- ggplot(tsnedf, aes(x = V1, y = V2)) +
    geom_point(size = 1, aes(color = pointcolor)) +
    theme_B3() +
    labs(x = "tSNE 1", y = "tSNE 2") +
    scale_color_manual(values = whichcolors) +
    theme(legend.position = "none",
          axis.text = element_blank()) +
    stat_ellipse(linetype = 1, aes(color = tissue )) 
  return(p)
}

plottsneelipsev2 <- function(tsnedf, pointcolor, whichcolors){
  p <- ggplot(tsnedf, aes(x = V1, y = V2)) +
    geom_point(size = 1, aes(color = pointcolor)) +
    theme_B3() +
    labs(x = "tSNE 1", y = "tSNE 2") +
    scale_color_manual(values = whichcolors) +
    theme(legend.position = "none",
          axis.text = element_blank()) +
    stat_ellipse(linetype = 1, aes(color = pointcolor)) 
  return(p)
}



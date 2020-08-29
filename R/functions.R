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

subsetcolData3 <- function(colData, eachgroup){
  
  colData <- colData %>%
    dplyr::filter(sextissue %in% eachgroup) %>%
    droplevels()
  colData <- as.data.frame(colData)
  row.names(colData) <- colData$sample
  return(colData)
}

subsetcolData4 <- function(colData, eachgroup){
  
  colData <- colData %>%
    dplyr::filter(tissue %in% eachgroup) %>%
    droplevels()
  colData <- as.data.frame(colData)
  row.names(colData) <- colData$sample
  return(colData)
}

subsetcountData3 <- function(df){
  savecols <- as.character(df$sample) 
  savecols <- as.vector(savecols)  
  df2 <- countData %>% dplyr::select(one_of(savecols)) 
  return(df2)
}


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

plot.volcano <- function(whichtissue, whichsex,  whichcomparison){
  
  volcano <- allDEG %>%
    filter(tissue == whichtissue,
           comparison == whichcomparison,
           sex %in% whichsex) %>%
    ggplot(aes(x = lfc, y = logpadj)) + 
    geom_point(aes(color = direction), size = 1, 
               alpha = 0.75, na.rm = T) + 
    theme_B3() +
    scale_color_manual(values = allcolors, 
                       name = "increased expression in:") +
    labs(y = "-log10(adj. p-value)", 
         x = "Log-fold change (LFC)") +
    theme(legend.position = "bottom",
          legend.direction = "vertical",
          legend.key.height = unit(0, "cm"),      
          legend.spacing.y = unit(0, "cm")) +
    guides(color = guide_legend(nrow = 1))
  return(volcano)
}

plot.volcano.sex <- function(df){
  
  volcano <- df %>%
    ggplot(aes(x = lfc, y = logpadj)) + 
    geom_point(aes(color = direction), size = 1, 
               alpha = 0.75, na.rm = T) + 
    theme_B3() +
    scale_color_manual(values = allcolors, 
                       name = NULL) +
    labs(y = "-log10(adj. p-value)", x = "Log-fold change (LFC)") +
    theme(legend.position = "top") 
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

# plot colorful pca and vector based pcas

plotpc12 <- function(df1color, df2fviz, whichfactor, whichcolors, mysubtitle, mytitle){
  
  pc12color <- df1color %>%	
    ggplot(aes(x = PC1, y = PC2, color = whichfactor )) +	
    geom_point(size = 1)  +	
    theme_B3() +	
    theme(legend.title = element_blank(),	
          axis.text = element_blank(),	
          legend.position = "none",
          axis.ticks = element_blank()) +	
    labs(x = "PC1", y = "PC2", subtitle = mysubtitle, title = mytitle)  +	
    scale_color_manual(values = whichcolors)  
  
  
  pc12vector <- fviz_pca_var(df2fviz,  
                             axes.linetype = "blank", 
                             repel = T , 
                             select.var= list(contrib = 5),
                             labelsize = 3)  + 
    labs(title = mytitle, subtitle = " ") + 
    theme_B3() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank()) 
  
  p <- plot_grid(pc12color, pc12vector, nrow = 1, rel_widths = c(1,1))
  return(p)
  
}


## bar graphs for PRL fig 4

makenewbargraph <- function(whichtissue, whichsex,  whichcomparison, lowlim, higherlim){
  p <- allDEG %>%
    filter(tissue == whichtissue,
           comparison == whichcomparison,
           sex == whichsex) %>%
    ggplot(aes(x = comparison,  fill = direction)) +
    geom_bar(position = "dodge", drop = FALSE) +
    theme_B3() +
    theme(legend.position = "none")  +
    guides(fill = guide_legend(nrow = 1)) +
    labs( y = "DEGs w/ + LFC", x = NULL) +
    geom_text(stat='count', aes(label=..count..), 
              vjust =-0.5, na.rm = F,
              position = position_dodge(width = 1),
              size = 2, color = "black")  + 
    ylim(lowlim, higherlim) +
    scale_fill_manual(values = allcolors, name = "higher in")    
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


createDEGdfsavesex <- function(up, down, mytissue){
  
  res <- results(dds, contrast = c("sex", up, down), independentFiltering = T, alpha = 0.1)
  
  DEGs <- data.frame(gene = row.names(res),
                     padj = res$padj, 
                     logpadj = -log10(res$padj),
                     lfc = res$log2FoldChange,
                     tissue = mytissue)
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
  myfilename = paste0("../results/DESeq2/sex/", mytissue, partialfilename, "_DEGs.csv")
  
  write.csv(DEGs, myfilename, row.names = F)
  # return DEGs frome with all data, included NS genes
  print(head(DEGs))
} 



createDEGdfsavestissue <- function(up, down){
  
  res <- results(dds, contrast = c("tissue", up, down), independentFiltering = T, alpha = 0.1)
  
  DEGs <- data.frame(gene = row.names(res),
                     padj = res$padj, 
                     logpadj = -log10(res$padj),
                     lfc = res$log2FoldChange)
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
  
  partialfilename = paste(down, "_", up, sep = "")
  myfilename = paste0("../results/DESeq2/tissue/", partialfilename, "_DEGs.csv")
  
  write.csv(DEGs, myfilename, row.names = F)
  # return DEGs frome with all data, included NS genes
  print(head(DEGs))
} 



## prolactin plots 

plotprolactin <- function(df, myy, myylab, mysubtitle ){
  
  p <-  ggplot(df, aes(x = treatment, y = myy)) +
    geom_boxplot(aes(fill = treatment, color = sex), outlier.shape = NA,  size = 0.25) +
    geom_jitter(size = 0.25, aes(color = sex)) +
    theme_B3() +
    scale_fill_manual(values = allcolors) +
    scale_color_manual(values = sexcolors) +
    theme(legend.position = c(0.85,0.15), legend.direction = "horizontal") + 
     theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none",
          axis.ticks = element_blank(),
          axis.title.y = element_text(face = "italic")) +
    labs(subtitle = mysubtitle, y= myylab) +
    guides(fill = guide_legend(nrow = 1)) +
    geom_signif(comparisons = list(c( "bldg", "lay"),
                                   c( "lay", "inc.d3"),
                                   c("inc.d3", "inc.d9"),
                                   c( "inc.d9", "inc.d17"),
                                   c( "inc.d17", "hatch"),
                                   c("hatch", "n5"),
                                   c( "n5", "n9")),  
                map_signif_level=TRUE,
                textsize = 2, family = 'Helvetica',
                vjust = 1.5, size = 0.2) 
  return(p)
}

## subsetmaketsne 

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
  
  tsne_model <- Rtsne(euclidist, check_duplicates=FALSE,
                      dims = 3, perplexity = 30 )
  tsne_df = as.data.frame(tsne_model$Y) 
  
  # prep for adding columns
  colData2 <- colData 
  colData2$V1 <- NULL
  tsne_df_cols <- cbind(colData2, tsne_df)
  return(tsne_df_cols)
}


plottsne <- function(tsnedf, pointcolor, whichcolors){
  p <- ggplot(tsnedf, aes(x = V1, y = V2)) +
    #stat_ellipse(linetype = 1, aes(color = treatment)) +
    geom_point(size = 1, aes(color = pointcolor, shape = sex)) +
    theme_B3() +
    labs(x = "tSNE 1", y = "tSNE 2", 
         title = " ") +
    scale_color_manual(values = whichcolors) +
    theme(legend.position = "none",
          axis.text = element_blank(), axis.ticks = element_blank()) 
  return(p)
}

makebargraph <- function(df, whichtissue, myylab, lowlim, higherlim, mylabels){
  p <- df %>%
    filter(tissue == whichtissue) %>%
    ggplot(aes(x = comparison,  fill = direction)) +
    geom_bar(position = "dodge") +
    facet_wrap(~sex) +
    theme_B3() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))  +
    labs(x = NULL, y = myylab) +
    scale_fill_manual(values = allcolors,
                      name = " ",
                      drop = FALSE) +
    scale_color_manual(values = allcolors) +
    geom_text(stat='count', aes(label=..count..), vjust =-0.5, 
              position = position_dodge(width = 1),
              size = 1.75, color = "black")  +
    ylim(lowlim, higherlim) +
    scale_x_discrete(labels = mylabels)
  return(p)
}



makebargraphv4 <- function(df, whichtissue, myylab, mybreaks, mylabels, whichsex){
  
  p <- df %>%
    filter(tissue == whichtissue,
           sex == whichsex) %>%
    ggplot(aes(x = comparison, y = n, fill = direction)) +
    geom_bar(stat="identity") +
    facet_wrap(~sex) +
    theme_B3() +
    theme(legend.position = "none")  +
    scale_fill_manual(values = allcolors,
                      name = " ",
                      drop = FALSE) +
    scale_color_manual(values = allcolors) +
    geom_text(stat='identity', aes(label= n), vjust =-0.5, 
              position = position_dodge(width = 1),
              size = 1.75, color = "black") +
    labs(x = NULL, y = myylab)  +
    scale_x_discrete(breaks = mybreaks,
                     labels = mylabels)
  return(p)
}




getcandidatevsd <- function(whichgenes, whichtissue, whichsex){
  candidates  <- allvsd %>%
    filter(gene %in% whichgenes) %>%
    dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
    dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/treatment/'), "[", 2)) %>%
    dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                  tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                  treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
    dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
    dplyr::select(sex, tissue, treatment, gene, samples, counts) %>%
    filter(tissue == whichtissue, sex %in% whichsex)  %>%
    drop_na()
  candidates$treatment <- factor(candidates$treatment, levels = alllevels)
  return(candidates)
}




#### modified rplot
# from https://github.com/tidymodels/corrr/blob/207ba60b3d2cd771643e97b5844369f05792c411/R/cor_df.R#L119 

rplot2 <- function(rdf,
                         legend = TRUE,
                         shape = 16,
                         colours = c("#2c7bb6", "#ffffbf",  "#d7191c"),
                         print_cor = FALSE,
                         colors) {
  
  if (!missing(colors))
    colours <- colors
  
  # Store order for factoring the variables
  row_order <- rdf$rowname
  
  # Convert data to relevant format for plotting
  pd <- stretch(rdf, na.rm = TRUE) 
  pd$size = abs(pd$r)
  pd$label = fashion(pd$r)
  
  plot_ <- list(
    # Geoms
    geom_point(shape = shape, size = 2),
    if (print_cor) geom_text(color = "black", size = 1, show.legend = FALSE),
    scale_colour_gradientn(limits = c(-1, 1), colors = colours),
    # Theme, labels, and legends
    theme_classic(),
    labs(x = "", y =""),
    guides(size = "none", alpha = "none"),
    if (!legend) theme(legend.position = "none")
  )
  
  ggplot(pd, aes_string(x = "x", y = "y", color = "r",
                        size = "size", alpha = "size",
                        label = "label")) +
    plot_
}

##### 3 funcitons  for correlation plots
makecorrdf <- function(whichsex, whichtissue, whichgenes){
  corrrdf <- candidatevsd %>%
    filter(sex == whichsex, tissue == whichtissue,
           gene %in% whichgenes) %>%
    drop_na() %>%
    pivot_wider(names_from = gene, values_from = counts) %>%
    select(-sex, -tissue, -treatment, -samples) %>%
    correlate() 
  # print(head(corrrdf))
  return(corrrdf)
}

subsetcandidatevsdwide <- function(whichsex, whichtissue){
  df <- candidatevsdwide %>%
    filter(sex == whichsex, tissue == whichtissue,
           treatment != "control") 
  return(df)
}

plotcorrplot <- function(df, mysubtitle){  
  p <- df %>%
    rplot2() +
    theme_B3() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
           axis.text = element_text(face = "italic"),
          axis.ticks = element_blank(),
          axis.line = element_blank()) +
    labs(subtitle = mysubtitle, y = NULL) + 
    theme(legend.position = "none")
  
  return(p)
}


## fig 3 top degs

plottopDEGs <- function(df, whichsex, myylab, mysubtitle){
  
  p <- df %>%
    arrange(desc(logpadj)) %>%
    filter(sex == whichsex) %>% head(10) %>% 
    ggplot(aes(x = reorder(gene, lfc), y = lfc)) + 
    geom_bar(stat = "identity",  aes(fill = direction)) +
    theme_B3() +
    theme(axis.text.y = element_text(face = "italic"),
          legend.position = "none") +
    labs(x = " ", y = myylab, subtitle = mysubtitle) +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = colorhypothesis) 
    
  return(p)
}

candidateboxplot <- function(whichtissue, whichgenes, whichsex){
  
  p <- candidatevsd %>%
    filter(tissue %in% whichtissue,
           gene %in% whichgenes,
           sex %in% whichsex) %>%
    ggplot(aes(x = treatment, y = counts)) +
    geom_boxplot(aes(fill = treatment, color = sex), outlier.shape = NA, size = 0.25) +
    geom_jitter(size = 0.25, aes(color = sex)) +
    facet_wrap(~sex,  nrow = 1) +
    scale_fill_manual(values = allcolors) +
    scale_color_manual(values = allcolors) +
    theme_B3() + 
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_text(face = "italic"),
          axis.ticks = element_blank()) +
    labs(y = whichgenes,
         x = NULL) +
    geom_signif(comparisons = list(c( "control", "bldg"),
                                   c( "bldg", "lay"),
                                   c( "lay", "inc.d3"),
                                   c("inc.d3", "inc.d9"),
                                   c( "inc.d9", "inc.d17"),
                                   c( "inc.d17", "hatch"),
                                   c("hatch", "n5"),
                                   c( "n5", "n9")),  
                map_signif_level=TRUE,
                textsize = 1.5, family = 'Helvetica',
                vjust = 1.5, size = 0) 
  
  return(p)
}

scattercorrelations <- function(df, gene1, myylab, gene2, myxlab,  mylinecolor){
  p <- ggplot(df, aes(x = gene2, y = gene1)) +
    labs(x = myxlab, y = myylab) + 
    scale_color_manual(values = allcolors) +
    geom_smooth(method = "glm",  aes(color = sex)) +
    geom_point(aes(color = treatment))  +
    theme_B3() +  
    theme(axis.title = element_text(face = "italic"), 
          legend.position = "none") +
    stat_cor( size = 2)
  return(p)
}



plotcandidatemanipquad <- function(df, whichtissue, whichgenes){
  
  plotcandidatemanip <- function(whichsex, whichlevels){
  
  p <- df %>%
    filter(tissue == whichtissue, sex == whichsex) %>%
    filter(gene %in% whichgenes) %>%
    filter(treatment %in% whichlevels) %>%
    mutate(treatment = factor(treatment, levels = alllevels)) %>%
    ggplot(aes(y =  counts, x = treatment, fill = treatment, color = sex)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.25, aes(color = sex)) +
    theme_B3() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          strip.text = element_text(face = "italic")) +
    scale_color_manual(values = allcolors) +
    scale_fill_manual(values = allcolors)  +
    labs(subtitle = paste(whichsex, whichtissue, sep = " "), y = whichgenes)
  
  return(p)
  
  }  
  
  p1 <- plotcandidatemanip("female", c(timinglevels, removallevels))  +
    geom_signif(comparisons = list(c( "inc.d9", "early"),
                                   c( "inc.d17", "prolong"),
                                   c( "hatch", "prolong"),
                                   c( "hatch", "extend")),
                map_signif_level=TRUE,
                textsize = 3, family = 'Helvetica',
                vjust = 1.5, size = 0)
  p2 <- plotcandidatemanip("male", c(timinglevels, removallevels)) +
    geom_signif(comparisons = list(c( "inc.d9", "early"),
                                   c( "inc.d17", "prolong"),
                                   c( "hatch", "prolong"),
                                   c( "hatch", "extend")),
                map_signif_level=TRUE,
                textsize = 3, family = 'Helvetica',
                vjust = 1.5, size = 0)
  
  p3 <- plotcandidatemanip("female", c(timinglevels, removallevels))  +
    geom_signif(comparisons = list(c( "inc.d3", "m.inc.d3"),
                                   c( "inc.d9", "m.inc.d9"),
                                   c( "inc.d17", "m.inc.d17"),
                                   c( "hatch", "m.n2")),
                map_signif_level=TRUE,
                textsize = 3, family = 'Helvetica',
                vjust = 1.5, size = 0)
  p4 <- plotcandidatemanip("male", c(timinglevels, removallevels))  +
    geom_signif(comparisons = list(c( "inc.d3", "m.inc.d3"),
                                   c( "inc.d9", "m.inc.d9"),
                                   c( "inc.d17", "m.inc.d17"),
                                   c( "hatch", "m.n2")),
                map_signif_level= TRUE,
                textsize = 3, family = 'Helvetica',
                vjust = 1.5, size =0)
  
  p <- plot_grid(p1,p3,p2,p4, nrow = 2, rel_widths = c(1,1.2))
  
  
  return(p)
}

## candidate gene box plot 

newcandidateboxplot <- function(whichtissue, whichgenes, whichsex, mytitle, mysubtitle){
  
  p <- candidatevsd  %>%
    dplyr::filter(sex %in% whichsex, 
                  gene %in% whichgenes, 
                  tissue %in% whichtissue) %>%
    DT::datatable() %>%
    ggplot(aes(x = treatment, y = counts)) +
    geom_boxplot(aes(fill = treatment, color = sex), outlier.shape = NA) +
    geom_jitter(size = 0.25, aes(color = sex)) +
    #facet_wrap(~sex,  nrow = 1) +
    scale_fill_manual(values = allcolors) +
    scale_color_manual(values = allcolors) +
    theme_B3() + 
    theme(legend.position = "none",
          axis.title.y = element_text(face = "italic"),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line.x = element_blank()) +
    labs(y = whichgenes,
         x = NULL,
         title = mytitle, 
         subtitle = mysubtitle) +
    geom_signif(comparisons = list( c( "control", "bldg"),
                                    c( "bldg", "lay"),
                                    c( "lay", "inc.d3"),
                                    c("inc.d3", "inc.d9"),
                                    c( "inc.d9", "inc.d17"),
                                    c( "inc.d17", "hatch"),
                                    c("hatch", "n5"),
                                    c( "n5", "n9")),  
                map_signif_level=TRUE,
                textsize = 1.5, family = 'Helvetica',
                vjust = 1.5, size = 0.2) 
  
  return(p)
}


## candidate gene box plot 

externalboxplots <- function(whichtissue, whichgenes, whichsex, mytitle, mysubtitle){
  
  p <- candidatevsd %>%
    dplyr::filter(sex %in% whichsex, gene %in% whichgenes, tissue %in% whichtissue) %>%
    DT::datatable() %>%
    ggplot(aes(x = external, y = counts)) +
    geom_boxplot(aes(fill = external, color = sex), outlier.shape = NA) +
    geom_jitter(size = 0.25, aes(color = treatment)) +
    #facet_wrap(~sex,  nrow = 1) +
    scale_fill_manual(values = allcolors) +
    scale_color_manual(values = allcolors) +
    theme_B3() + 
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.line.x = element_blank()) +
    labs(y = whichgenes,
         x = NULL,
         title = mytitle, 
         subtitle = mysubtitle) +
    geom_signif(comparisons = list( c( "none", "eggs"),
                                    c( "eggs", "chicks")),  
                map_signif_level=TRUE,
                textsize = 1.5, family = 'Helvetica',
                vjust = 1.5, size = 0.2) 
  
  return(p)
}


# for manipulation tsne 

addgroupings <- function(df){
  
  df <- df %>% mutate(earlylate = fct_collapse(treatment, 
                                               early = c("early", "m.inc.d3", "m.inc.d9", "inc.d3", 
                                                         "inc.d9",  "bldg"),
                                               late = c("inc.d17",   "m.inc.d17", "prolong" ,  "hatch" , 
                                                        "m.n2", "extend", "n5"),
                                               reference = "control"),
                      extint = fct_collapse(treatment, 
                                            nest = c("m.inc.d3", "m.inc.d9", "m.n2", "m.inc.d17", "bldg"),
                                            eggs = c("inc.d3", "inc.d9", "inc.d17", "prolong"),
                                            chicks = c("early", "hatch", "extend", "n5"),
                                            reference = "control")) 
  return(df)
}


## box plots 

plotcandidatechar <- function(df, whichgene){
  
  p <- df %>%
    dplyr::filter(gene  == whichgene,
                  treatment %in% charlevels)  %>%
    ggplot(aes(y =  counts, x = treatment, fill = treatment, color = sex)) +
    geom_boxplot(outlier.shape = NA, lwd=0.5) +
    facet_wrap(~sex) +
    geom_jitter(size = 0.25, aes(color = sex)) +
    theme_B3() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_text(face = "italic")) +
    scale_color_manual(values = allcolors) +
    scale_fill_manual(values = allcolors)  +
    labs(y = whichgene, x = "Sequential stages") +
    geom_signif(comparisons = list(c("control", "bldg"), 
                                   c( "bldg", "lay"), 
                                   c( "lay", "inc.d3")),
                map_signif_level=TRUE,
                textsize = 2, family = 'Helvetica',
                vjust = 0, size = 0.5, step_increase = 0.075) +
    geom_signif(comparisons = list(c( "inc.d3", "inc.d9"),
                                   c( "inc.d9", "inc.d17"),
                                   c( "inc.d17", "hatch")),
                map_signif_level=TRUE,
                textsize = 2, family = 'Helvetica',
                vjust = 0,  size = 0.5, step_increase = 0.075) +
    geom_signif(comparisons = list(c( "hatch", "n5"), 
                                   c( "n5", "n9")),
                map_signif_level=TRUE,
                textsize = 2, family = 'Helvetica',
                vjust = 0, size = 0.5, step_increase = 0.075)
  
  return(p)
  
}

plotremoval <- function(df, whichgene){
  
  p <- df %>%
    dplyr::filter(gene  == whichgene,
                  treatment %in% removallevels)  %>%
    ggplot(aes(y =  counts, x = treatment, fill = treatment, color = sex)) +
    geom_boxplot(outlier.shape = NA, lwd=0.5) +
    facet_wrap(~sex) +
    geom_jitter(size = 0.25, aes(color = sex)) +
    theme_B3() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_text(face = "italic")) +
    scale_color_manual(values = allcolors) +
    scale_fill_manual(values = allcolors)  +
    labs( y = whichgene,  x = "Offspring removal") +
    
    geom_signif(comparisons = list(c( "inc.d3", "m.inc.d3"),
                                   c( "inc.d9", "m.inc.d9"),
                                   c( "inc.d17", "m.inc.d17")),
                map_signif_level=TRUE,
                textsize = 2, family = 'Helvetica',
                vjust = 0,  size = 0.5, step_increase = 0.075) +
    geom_signif(comparisons = list(c( "hatch", "m.n2")),
                map_signif_level=TRUE,
                textsize = 2, family = 'Helvetica',
                vjust = 0,  size = 0.5, step_increase = 0.075) 
  return(p)
  
}


plotreplacement <- function(df, whichgene){
  
  p <- df %>%
    dplyr::filter(gene  == whichgene,
                  treatment %in% timinglevels)  %>%
    ggplot(aes(y =  counts, x = treatment, fill = treatment, color = sex)) +
    geom_boxplot(outlier.shape = NA, lwd=0.5) +
    facet_wrap(~sex) +
    geom_jitter(size = 0.25, aes(color = sex)) +
    theme_B3() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.y = element_text(face = "italic")) +
    scale_color_manual(values = allcolors) +
    scale_fill_manual(values = allcolors)  +
    labs(y = whichgene,  x = "Offspring replacement") +

    geom_signif(comparisons = list(c("inc.d9", "early"),
                                   c("inc.d17", "prolong"),
                                   c( "hatch", "extend")),
                map_signif_level=TRUE,
                textsize = 2, family = 'Helvetica',
                vjust = 0, size = 0.5, step_increase = 0.075) +
    geom_signif(comparisons = list(c("hatch", "early"),
                                   c("hatch", "prolong"),
                                   c( "n5", "extend")),
                map_signif_level=TRUE,
                textsize = 2, family = 'Helvetica',
                vjust = 0, size = 0.5, step_increase = 0.075) +
    geom_signif(comparisons = list(),
                map_signif_level=TRUE,
                textsize = 2, family = 'Helvetica',
                vjust = 0, size = 0.5, step_increase = 0.075) +
    geom_signif(comparisons = list(),
                map_signif_level=TRUE,
                textsize = 2, family = 'Helvetica',
                vjust = 0, size = 0.5, step_increase = 0.075) 
  return(p)
  
}


boxplotextint <- function(df, whichgene, whichtissue){
  
  p1 <- df %>%
    dplyr::filter(gene  == whichgene)  %>%
    mutate(earlylate = factor(earlylate, levels = levelhypothesis)) %>%
    ggplot(aes(y =  counts, x = earlylate, fill = earlylate, color = sex)) +
    geom_boxplot(outlier.shape = NA, lwd=0.5) +
    #facet_wrap(~sex) +
    geom_jitter(size = 0.25, aes(color = treatment)) +
    theme_B3() +
    scale_color_manual(values = allcolors) +
    scale_fill_manual(values = allcolors)  +
    geom_signif(comparisons = list(c( "reference", "early"), 
                                   c( "early", "late"),
                                   c( "reference", "late")),
                map_signif_level=TRUE,
                textsize = 2, family = 'Helvetica',
                vjust = 0, size = 0.5, step_increase = 0.075) +
    labs(subtitle = " ", y = NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- df %>%
    dplyr::filter(gene  == whichgene)  %>%
    mutate(extint = factor(extint, levels = levelhypothesis)) %>%
    ggplot(aes(y =  counts, x = extint, fill = extint, color = sex)) +
    geom_boxplot(outlier.shape = NA, lwd=0.5) +
    #facet_wrap(~sex) +
    geom_jitter(size = 0.25, aes(color = treatment)) +
    theme_B3() +
    scale_color_manual(values = allcolors) +
    scale_fill_manual(values = allcolors)  +
    geom_signif(comparisons = list(c( "reference", "eggs"), 
                                   c( "eggs", "chicks"),
                                   c( "chicks", "loss"),
                                   c( "reference", "chicks"),
                                   c( "eggs", "loss"),
                                   c( "reference", "loss")),
                map_signif_level=TRUE,
                textsize = 2, family = 'Helvetica',
                vjust = 0, size = 0.5, step_increase = 0.075) +
    labs(subtitle = " " , y = NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  p <- plot_grid(p1,p2, rel_widths = c(3.5,4))
  return(p)
}
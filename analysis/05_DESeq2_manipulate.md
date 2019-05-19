    library(tidyverse)
    library(DESeq2)
    library(cowplot)
    library(RColorBrewer)
    library(pheatmap)
    library(kableExtra)
    library(viridis)
    library(edgeR)

    # load custom functions  
    source("../R/functions.R") 

    knitr::opts_chunk$set(fig.path = '../figures/manipulation/', cache = TRUE)

Manipulation data
-----------------

    # import "colData" which contains sample information and "countData" which contains read counts
    m.colData <- read.csv("../metadata/00_colData_manipluation.csv", header = T, row.names = 1)
    m.countData <- read.csv("../results/00_countData_manipluation.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)

    # set levels
    m.colData$treatment <- factor(m.colData$treatment, levels = 
                                  c("m.inc.d3",  "m.inc.d8",
                                    "m.inc.d9", "m.inc.d17",
                                    "prolong", "extend", "m.n2"))

    m.colData$sextissue <- as.factor(paste(m.colData$sex, m.colData$tissue, sep = "_"))

    m.colData$outcome <- ifelse(grepl("d3|d9|d17", m.colData$treatment), "end inc", 
                         ifelse(grepl("d8|n2", m.colData$treatment),"end hatch",
                         ifelse(grepl("prolong", m.colData$treatment),"prolong inc",
                         ifelse(grepl("extend", m.colData$treatment),"delay hatch", NA))))

    m.colData$outcome <- factor(m.colData$outcome, levels = 
                                  c("end inc",  "end hatch",
                                    "prolong inc", "delay hatch"))
    summary(m.colData[c(7,3,4,5,8,9)])

    ##           study         sex               tissue        treatment 
    ##  manipulation:411   female:208   gonad       :136   m.inc.d3 :60  
    ##                     male  :203   hypothalamus:138   m.inc.d8 :60  
    ##                                  pituitary   :137   m.inc.d9 :49  
    ##                                                     m.inc.d17:63  
    ##                                                     prolong  :60  
    ##                                                     extend   :60  
    ##                                                     m.n2     :59  
    ##                sextissue         outcome   
    ##  female_gonad       :69   end inc    :172  
    ##  female_hypothalamus:70   end hatch  :119  
    ##  female_pituitary   :69   prolong inc: 60  
    ##  male_gonad         :67   delay hatch: 60  
    ##  male_hypothalamus  :68                    
    ##  male_pituitary     :68                    
    ## 

Make a PCA plot, heat map, and box plots for each tissue for each sex
---------------------------------------------------------------------

    ## subset for test 
    m.colData <- m.colData %>%
      filter(sextissue == "male_hypothalamus") %>%
     droplevels()

    for (eachgroup in levels(m.colData$sextissue)){
      
      print(eachgroup)
      
      colData <- m.colData %>%
          dplyr::filter(sextissue == eachgroup) %>%
          droplevels()
      row.names(colData) <- colData$V1
      
      savecols <- as.character(colData$V1) 
      savecols <- as.vector(savecols) 

      countData <- m.countData %>% dplyr::select(one_of(savecols)) 

      # check that row and col lenghts are equal
      print(ncol(countData) == nrow(colData))

      dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ treatment )
      
      print(dds)
      nrow(colData) > 50
      
      dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
      print(dds)
      
      dds <- DESeq(dds) # Differential expression analysis
      
      vsd <- vst(dds, blind=FALSE) # variance stabilized 


    #create list of groups
    a <- levels(colData$treatment)
    b <- levels(colData$treatment)


    # comapre all contrasts, save to datafrmes
    totalDEGS=data.frame()
    for (i in a){
      for (j in b){
        if (i != j) {
          k <- paste(i,j, sep = "") #assigns usique rownames
          print(k)
          totalDEGS[k,1]<-i               
          totalDEGS[k,2]<-j
          totalDEGS[k,3]<- numDEGs(i,j) #caluculates number of DEGs
        }
      }
    }

    head(totalDEGS)

    # widen data to create table of degs
    rownames(totalDEGS) <- NULL #remove row names
    totalDEGS_wide <- spread(totalDEGS, V2, V3)
    print(totalDEGS_wide) 

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
        scale_fill_viridis(na.value="#FFFFFF00", 
                         limits = c(0, 3000),
                         breaks = c(0, 1000, 2000, 3000)) + 
        xlab(" ") + ylab("Timepoint") +
        labs(fill = "# of DEGs",
             subtitle = "eachgroup")
    plot(allcontrasts)


    # create the dataframe using my function pcadataframe
    pcadata <- pcadataframe(vsd, intgroup=c("treatment"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))
    print(percentVar)


    pca1 <- ggplot(pcadata, aes(treatment, PC1,color = treatment)) + 
      geom_boxplot() +
      ylab(paste0("PC1: ", percentVar[1],"% variance")) +
      xlab(NULL) +
      theme_cowplot(font_size = 8, line_size = 0.25) +
      labs(subtitle = "eachgroup") +
      theme(legend.position = "none")


    pca2 <- ggplot(pcadata, aes(treatment, PC2,color = treatment)) + 
      geom_boxplot() +
      ylab(paste0("PC2: ", percentVar[2],"% variance")) +
      xlab(NULL) +
      theme_cowplot(font_size = 8, line_size = 0.25) +
      labs(subtitle = "eachgroup") +
      theme(legend.position = "none")

    plot_grid(pca1, pca2)



    print(summary(aov(PC1 ~ treatment, data=pcadata)))
    print(TukeyHSD(aov(PC1 ~ treatment, data=pcadata), which = "treatment"))

    print(summary(aov(PC2 ~ treatment, data=pcadata))) 
    print(summary(aov(PC3 ~ treatment, data=pcadata))) 
    print(summary(aov(PC4 ~ treatment, data=pcadata))) 


    ## heamap with minimum pvalue

    # make dataframe counts
    DEGs <- assay(vsd)
    DEGs <- as.data.frame(DEGs)

    # add pvalues

    for (i in a){
      for (j in b){
        if (i != j) {
          k <- paste(i,j, sep = " vs ")
          print(k)
          results <- resvals2(i,j)
          DEGs <- cbind(DEGs,results)
        }
      }
    }



    DEGsmatrix <- DEGs
    DEGsmatrix <- as.matrix(DEGs)
    padjmin <- rowMins(DEGsmatrix, na.rm = T) 
    padjmin <- as.data.frame(padjmin)

    sigDEGs <- cbind(DEGs,padjmin)

    sigDEGs <- sigDEGs %>% arrange(padjmin)
    sigDEGs <- head(sigDEGs,500)
    sigDEGs <- as.data.frame(sigDEGs)
    rownames(sigDEGs) <- sigDEGs$rownames
    drop.cols <-colnames(sigDEGs[,grep("padj|pval|pmin|rownames", colnames(sigDEGs))])
    sigDEGs <- sigDEGs %>% dplyr::select(-one_of(drop.cols))
    sigDEGs <- as.matrix(sigDEGs)
    sigDEGs <- sigDEGs - rowMeans(sigDEGs)

    paletteLength <- 30
    myBreaks <- c(seq(min(sigDEGs), 0, length.out=ceiling(paletteLength/2) + 1), 
                    seq(max(sigDEGs)/paletteLength, max(sigDEGs), length.out=floor(paletteLength/2)))

    anndf <- m.colData %>% dplyr::select(treatment)
    rownames(anndf) <- m.colData$V1

    sigDEGs <- as.matrix(sigDEGs) 
    pheatmap(sigDEGs, show_rownames = F, show_colnames = F,
             color = viridis(30),
             breaks=myBreaks,
             annotation_col=anndf,
             main = "eachgroup")

    pheatmap(sigDEGs, kmeans_k = 5,
             show_rownames = F, show_colnames = F,
             color = viridis(30),
             breaks=myBreaks,
             annotation_col=anndf,
             main = "eachgroup")

    ## candidate genes
    names(geneinfo)
    names(geneinfo)[4] <- "rownames"
    DEGs$rownames <- row.names(DEGs)

    # make dataframe with geneids and names and counts

    candidates <- full_join(geneinfo, DEGs)
    drop.cols <-colnames(candidates[,grep("padj|pval|pmin", colnames(candidates))])
    candidates <- candidates %>% dplyr::select(-one_of(drop.cols))
    candidates <- candidates %>%
      filter(Name %in% c( #"JUN", "JUND", "EGR",  "AVP", "AVPR1A", "AVPR1B", "AVPR2", "OXT",  
                         "AR", "CYP19A1", "ESR1", "ESR2", "FSHR",
                         "GHRL", "GAL", "NPVF", "GNRH1", "LHCGR",
                         "PGR", "PRL", "PRLR", "VIP", "VIPR1")) 
    row.names(candidates) <- candidates$Name
    candidates <- candidates %>% select(-row.names, -rownames, -Name, -geneid)
    candidates <- candidates %>% drop_na()

    head(candidates)

    candidates <- as.data.frame(t(candidates))
    candidates$RNAseqID <- rownames(candidates)
    head(candidates)

    candidates <- candidates %>% gather(gene, value, -RNAseqID)  %>% # https://tidyr.tidyverse.org/reference/gather.html
      filter(RNAseqID != "gene")
    candidates$value <- as.numeric(candidates$value)
    candidates$V1  <- candidates$RNAseqID
    head(candidates)

    candidatecounts <- left_join(candidates, colData)
    #head(candidatecounts)

    candidatecounts$faketime <- as.numeric(candidatecounts$treatment)
    #head(candidatecounts)

    candidatecounts$gene <- as.factor(candidatecounts$gene)
    levels(candidatecounts$gene)

    p1 <- candidatecounts %>%
      filter(gene %in% c( "AR", "CYP19A1", "ESR1", "ESR2", "FSHR")) %>%
      ggplot(aes(x = treatment, y = value, fill = treatment)) +
      geom_boxplot() +
      facet_wrap(~gene, scales = "free", nrow = 1) +
      theme(axis.text.x = element_blank(),
            legend.position = "bottom") +
      labs(subtitle = "eachgroup")
    p1

    print(p1)

    cat("\n")

    p2 <- candidatecounts %>%
      filter(gene %in% c("GHRL", "GAL", "NPVF", "GNRH1", "LHCGR")) %>%
      ggplot(aes(x = treatment, y = value, fill = treatment)) +
      geom_boxplot() +
      facet_wrap(~gene, scales = "free", nrow = 1) +
      theme(axis.text.x = element_blank(),
            legend.position = "bottom") +
      labs(subtitle = "eachgroup")
    p2

    print(p2)

    cat("\n")

    p3 <- candidatecounts %>%
      filter(gene %in% c("PGR", "PRL", "PRLR", "VIP", "VIPR1")) %>%
      ggplot(aes(x = treatment, y = value, fill = treatment)) +
      geom_boxplot() +
      facet_wrap(~gene, scales = "free", nrow = 1) +
      theme(axis.text.x = element_blank(),
            legend.position = "bottom") +
      labs(subtitle = "eachgroup")
    p3


    print(p3)

    cat("\n")

    cat("Done with eachgroup")

    cat("\n")

    }

    ## [1] "male_hypothalamus"
    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 68 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(68): blk.s030.o.g_male_hypothalamus_prolong
    ##   blk5.x_male_hypothalamus_m.inc.d3 ...
    ##   y55.x_male_hypothalamus_m.inc.d8
    ##   y63.x_male_hypothalamus_m.inc.d9
    ## colData names(9): V1 bird ... sextissue outcome
    ## class: DESeqDataSet 
    ## dim: 14339 68 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14339): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(68): blk.s030.o.g_male_hypothalamus_prolong
    ##   blk5.x_male_hypothalamus_m.inc.d3 ...
    ##   y55.x_male_hypothalamus_m.inc.d8
    ##   y63.x_male_hypothalamus_m.inc.d9
    ## colData names(9): V1 bird ... sextissue outcome

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 7 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    ## [1] "m.inc.d3m.inc.d8"
    ## [1] "m.inc.d3m.inc.d9"
    ## [1] "m.inc.d3m.inc.d17"
    ## [1] "m.inc.d3prolong"
    ## [1] "m.inc.d3extend"
    ## [1] "m.inc.d3m.n2"
    ## [1] "m.inc.d8m.inc.d3"
    ## [1] "m.inc.d8m.inc.d9"
    ## [1] "m.inc.d8m.inc.d17"
    ## [1] "m.inc.d8prolong"
    ## [1] "m.inc.d8extend"
    ## [1] "m.inc.d8m.n2"
    ## [1] "m.inc.d9m.inc.d3"
    ## [1] "m.inc.d9m.inc.d8"
    ## [1] "m.inc.d9m.inc.d17"
    ## [1] "m.inc.d9prolong"
    ## [1] "m.inc.d9extend"
    ## [1] "m.inc.d9m.n2"
    ## [1] "m.inc.d17m.inc.d3"
    ## [1] "m.inc.d17m.inc.d8"
    ## [1] "m.inc.d17m.inc.d9"
    ## [1] "m.inc.d17prolong"
    ## [1] "m.inc.d17extend"
    ## [1] "m.inc.d17m.n2"
    ## [1] "prolongm.inc.d3"
    ## [1] "prolongm.inc.d8"
    ## [1] "prolongm.inc.d9"
    ## [1] "prolongm.inc.d17"
    ## [1] "prolongextend"
    ## [1] "prolongm.n2"
    ## [1] "extendm.inc.d3"
    ## [1] "extendm.inc.d8"
    ## [1] "extendm.inc.d9"
    ## [1] "extendm.inc.d17"
    ## [1] "extendprolong"
    ## [1] "extendm.n2"
    ## [1] "m.n2m.inc.d3"
    ## [1] "m.n2m.inc.d8"
    ## [1] "m.n2m.inc.d9"
    ## [1] "m.n2m.inc.d17"
    ## [1] "m.n2prolong"
    ## [1] "m.n2extend"
    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA         0        0        0     1099    1       0
    ## 2 m.inc.d17      0        NA        0        1      658    0       0
    ## 3  m.inc.d3      0         0       NA        0      832    0       0
    ## 4  m.inc.d8      0         1        0       NA     2934    0       1
    ## 5  m.inc.d9   1099       658      832     2934       NA  853     550
    ## 6      m.n2      1         0        0        0      853   NA       0
    ## 7   prolong      0         0        0        1      550    0      NA

![](../figures/manipulation/althings-1.png)

    ## [1] 21  9  6  5  3  3
    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment    6  672.9  112.15   5.567 0.000118 ***
    ## Residuals   61 1228.8   20.14                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                          diff        lwr        upr     p adj
    ## m.inc.d8-m.inc.d3  -2.9425407  -9.062077  3.1769957 0.7633728
    ## m.inc.d9-m.inc.d3   8.4818829   1.991134 14.9726314 0.0033023
    ## m.inc.d17-m.inc.d3  0.2285568  -5.890980  6.3480931 0.9999998
    ## prolong-m.inc.d3    1.0518632  -5.067673  7.1713995 0.9983936
    ## extend-m.inc.d3    -0.3208496  -6.440386  5.7986867 0.9999985
    ## m.n2-m.inc.d3      -1.2802542  -7.399791  4.8392822 0.9952286
    ## m.inc.d9-m.inc.d8  11.4244236   4.933675 17.9151720 0.0000263
    ## m.inc.d17-m.inc.d8  3.1710975  -2.948439  9.2906338 0.6953497
    ## prolong-m.inc.d8    3.9944038  -2.125133 10.1139402 0.4314036
    ## extend-m.inc.d8     2.6216910  -3.497845  8.7412274 0.8464640
    ## m.n2-m.inc.d8       1.6622865  -4.457250  7.7818228 0.9810358
    ## m.inc.d17-m.inc.d9 -8.2533261 -14.744075 -1.7625776 0.0046418
    ## prolong-m.inc.d9   -7.4300197 -13.920768 -0.9392713 0.0149496
    ## extend-m.inc.d9    -8.8027325 -15.293481 -2.3119841 0.0020265
    ## m.n2-m.inc.d9      -9.7621371 -16.252886 -3.2713886 0.0004413
    ## prolong-m.inc.d17   0.8233064  -5.296230  6.9428427 0.9996024
    ## extend-m.inc.d17   -0.5494065  -6.668943  5.5701299 0.9999625
    ## m.n2-m.inc.d17     -1.5088110  -7.628347  4.6107254 0.9885084
    ## extend-prolong     -1.3727128  -7.492249  4.7468235 0.9930426
    ## m.n2-prolong       -2.3321173  -8.451654  3.7874190 0.9054801
    ## m.n2-extend        -0.9594045  -7.078941  5.1601318 0.9990452
    ## 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    6   73.2   12.21   1.019  0.422
    ## Residuals   61  730.9   11.98               
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    6   67.7  11.292   1.577  0.169
    ## Residuals   61  436.7   7.159               
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    6   28.6   4.768   0.745  0.616
    ## Residuals   61  390.4   6.400               
    ## [1] "m.inc.d3 vs m.inc.d8"
    ## [1] 0
    ## [1] "m.inc.d3 vs m.inc.d9"
    ## [1] 832
    ## [1] "m.inc.d3 vs m.inc.d17"
    ## [1] 0
    ## [1] "m.inc.d3 vs prolong"
    ## [1] 0
    ## [1] "m.inc.d3 vs extend"
    ## [1] 0
    ## [1] "m.inc.d3 vs m.n2"
    ## [1] 0
    ## [1] "m.inc.d8 vs m.inc.d3"
    ## [1] 0
    ## [1] "m.inc.d8 vs m.inc.d9"
    ## [1] 2934
    ## [1] "m.inc.d8 vs m.inc.d17"
    ## [1] 1
    ## [1] "m.inc.d8 vs prolong"
    ## [1] 1
    ## [1] "m.inc.d8 vs extend"
    ## [1] 0
    ## [1] "m.inc.d8 vs m.n2"
    ## [1] 0
    ## [1] "m.inc.d9 vs m.inc.d3"
    ## [1] 832
    ## [1] "m.inc.d9 vs m.inc.d8"
    ## [1] 2934
    ## [1] "m.inc.d9 vs m.inc.d17"
    ## [1] 658
    ## [1] "m.inc.d9 vs prolong"
    ## [1] 550
    ## [1] "m.inc.d9 vs extend"
    ## [1] 1099
    ## [1] "m.inc.d9 vs m.n2"
    ## [1] 853
    ## [1] "m.inc.d17 vs m.inc.d3"
    ## [1] 0
    ## [1] "m.inc.d17 vs m.inc.d8"
    ## [1] 1
    ## [1] "m.inc.d17 vs m.inc.d9"
    ## [1] 658
    ## [1] "m.inc.d17 vs prolong"
    ## [1] 0
    ## [1] "m.inc.d17 vs extend"
    ## [1] 0
    ## [1] "m.inc.d17 vs m.n2"
    ## [1] 0
    ## [1] "prolong vs m.inc.d3"
    ## [1] 0
    ## [1] "prolong vs m.inc.d8"
    ## [1] 1
    ## [1] "prolong vs m.inc.d9"
    ## [1] 550
    ## [1] "prolong vs m.inc.d17"
    ## [1] 0
    ## [1] "prolong vs extend"
    ## [1] 0
    ## [1] "prolong vs m.n2"
    ## [1] 0
    ## [1] "extend vs m.inc.d3"
    ## [1] 0
    ## [1] "extend vs m.inc.d8"
    ## [1] 0
    ## [1] "extend vs m.inc.d9"
    ## [1] 1099
    ## [1] "extend vs m.inc.d17"
    ## [1] 0
    ## [1] "extend vs prolong"
    ## [1] 0
    ## [1] "extend vs m.n2"
    ## [1] 1
    ## [1] "m.n2 vs m.inc.d3"
    ## [1] 0
    ## [1] "m.n2 vs m.inc.d8"
    ## [1] 0
    ## [1] "m.n2 vs m.inc.d9"
    ## [1] 853
    ## [1] "m.n2 vs m.inc.d17"
    ## [1] 0
    ## [1] "m.n2 vs prolong"
    ## [1] 0
    ## [1] "m.n2 vs extend"
    ## [1] 1

![](../figures/manipulation/althings-2.png)

    ## Joining, by = "rownames"

    ## Warning: Column `rownames` joining factor and character vector, coercing
    ## into character vector

    ## Joining, by = "V1"

    ## Warning: Column `V1` joining character vector and factor, coercing into
    ## character vector

![](../figures/manipulation/althings-3.png)![](../figures/manipulation/althings-4.png)
![](../figures/manipulation/althings-5.png)
![](../figures/manipulation/althings-6.png)

    ## 
    ## Done with eachgroup

    library(tidyverse)
    library(DESeq2)
    library(cowplot)
    library(RColorBrewer)
    library(pheatmap)
    library(kableExtra)
    library(viridis)

    # load custom functions  
    source("../R/functions.R") 

    knitr::opts_chunk$set(fig.path = '../figures/manipulation/',cache=TRUE)

Manipulation data
=================

    # import "colData" which contains sample information and "countData" which contains read counts
    m.colData <- read.csv("../results/00_colData_manipluation.csv", header = T, row.names = 1)
    m.countData <- read.csv("../results/00_countData_manipluation.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../results/00_geneinfo.csv", row.names = 1)

    # set levels
    m.colData$treatment <- factor(m.colData$treatment, levels = 
                                  c("m.inc.d3", "m.inc.d8", "m.inc.d9",
                                    "m.inc.d17",  "prolong", "extend", "m.n2"))

    m.colData$sextissue <- as.factor(paste(m.colData$sex, m.colData$tissue, sep = "_"))

Write for loop to do this one for every tissue and for every treatment
======================================================================

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
      dds <- dds[ rowSums(counts(dds)) > 2, ] ## pre-filter genes 
      dds <- DESeq(dds) # Differential expression analysis
      vsd <- vst(dds, blind=FALSE) # variance stabilized 

    #create list of groups
    a <- levels(colData$treatment)
    b <- levels(colData$treatment)

    # comapre all contrasts, save to datafrmes
    dat=data.frame()
    for (i in a){
      for (j in b){
        if (i != j) {
          k <- paste(i,j, sep = "") #assigns usique rownames
          dat[k,1]<-i               
          dat[k,2]<-j
          dat[k,3]<- numDEGs(i,j) #caluculates number of DEGs
        }
      }
    }

    head(dat)

    # widen data to create table of degs
    rownames(dat) <- NULL #remove row names
    data_wide <- spread(dat, V2, V3)
    print(data_wide) 

    dat$V1 <- factor(dat$V1, levels = 
                                  c("m.inc.d3", "m.inc.d8", "m.inc.d9",
                                    "m.inc.d17",  "prolong", "extend", "m.n2"))
    dat$V2 <- factor(dat$V2, levels = 
                                  c("m.inc.d3", "m.inc.d8", "m.inc.d9",
                                    "m.inc.d17",  "prolong", "extend", "m.n2"))

    allcontrasts <- dat %>%
      ggplot( aes(V1, V2)) +
        geom_tile(aes(fill = V3)) +
        scale_fill_viridis(na.value="#FFFFFF00", 
                         limits = c(0, 5025),
                         breaks = c(0, 1000, 2000, 3000, 4000, 5000)) + 
        xlab(" ") + ylab("Timepoint") +
        labs(fill = "# of DEGs",
             subtitle = eachgroup)
    plot(allcontrasts)

    # create the dataframe using my function pcadataframe
    pcadata <- pcadataframe(vsd, intgroup=c("treatment"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))
    percentVar

    pca12 <- ggplot(pcadata, aes(PC1, PC2,color = treatment)) + 
      geom_point(size = 2, alpha = 1) +
      stat_ellipse(type = "t") +
      xlab(paste0("PC1: ", percentVar[1],"% variance")) +
      ylab(paste0("PC2: ", percentVar[2],"% variance")) +
      theme_cowplot(font_size = 8, line_size = 0.25) +
      labs(subtitle = eachgroup)
    print(pca12)

    summary(aov(PC1 ~ treatment, data=pcadata)) 
    TukeyHSD(aov(PC1 ~ treatment, data=pcadata), which = "treatment") 

    summary(aov(PC2 ~ treatment, data=pcadata)) 
    TukeyHSD(aov(PC2 ~ treatment, data=pcadata), which = "treatment") 

    pca34 <- ggplot(pcadata, aes(PC3, PC4,color = treatment)) + 
      geom_point(size = 2, alpha = 1) +
      stat_ellipse(type = "t") +
      xlab(paste0("PC3: ", percentVar[3],"% variance")) +
      ylab(paste0("PC4: ", percentVar[4],"% variance")) +
      theme_cowplot(font_size = 8, line_size = 0.25) +
      labs(subtitle = eachgroup)
    print(pca34)

    summary(aov(PC3 ~ treatment, data=pcadata)) 

    summary(aov(PC4 ~ treatment, data=pcadata)) 

    summary(aov(PC5 ~ treatment, data=pcadata)) 

    summary(aov(PC6 ~ treatment, data=pcadata)) 

    # see http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix
    sampleDists <- dist(t(assay(vsd)))

    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- NULL
    colnames(sampleDistMatrix) <- colData$treatment
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors,
             fontsize = 6, 
             main = eachgroup)
    }

    ## [1] "female_gonad"
    ## [1] TRUE

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 507 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA       810      102      170      301  777     300
    ## 2 m.inc.d17    810        NA      860      589        1    0     610
    ## 3  m.inc.d3    102       860       NA        0       41  678     640
    ## 4  m.inc.d8    170       589        0       NA       51  734    1440
    ## 5  m.inc.d9    301         1       41       51       NA    3     129
    ## 6      m.n2    777         0      678      734        3   NA     425
    ## 7   prolong    300       610      640     1440      129  425      NA

![](../figures/manipulation/althings-1.png)

    ## Warning in MASS::cov.trob(data[, vars]): Probable convergence failure

![](../figures/manipulation/althings-2.png)![](../figures/manipulation/althings-3.png)

    ## [1] "female_hypothalamus"
    ## [1] TRUE

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 26 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

![](../figures/manipulation/althings-4.png)

    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA       336     1051        0      609   16       6
    ## 2 m.inc.d17    336        NA        0      697      448   15      45
    ## 3  m.inc.d3   1051         0       NA     1630     1276  239     116
    ## 4  m.inc.d8      0       697     1630       NA      192    5      20
    ## 5  m.inc.d9    609       448     1276      192       NA   88     935
    ## 6      m.n2     16        15      239        5       88   NA     142
    ## 7   prolong      6        45      116       20      935  142      NA

![](../figures/manipulation/althings-5.png)![](../figures/manipulation/althings-6.png)![](../figures/manipulation/althings-7.png)

    ## [1] "female_pituitary"
    ## [1] TRUE

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 65 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

![](../figures/manipulation/althings-8.png)

    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA      2014     2670     2837     2842 1381     284
    ## 2 m.inc.d17   2014        NA      317     1978     1525   13    1272
    ## 3  m.inc.d3   2670       317       NA      187      522 1249    1700
    ## 4  m.inc.d8   2837      1978      187       NA     1156 3031    2380
    ## 5  m.inc.d9   2842      1525      522     1156       NA 1621    1547
    ## 6      m.n2   1381        13     1249     3031     1621   NA     943
    ## 7   prolong    284      1272     1700     2380     1547  943      NA

![](../figures/manipulation/althings-9.png)![](../figures/manipulation/althings-10.png)

    ## Warning in MASS::cov.trob(data[, vars]): Probable convergence failure

![](../figures/manipulation/althings-11.png)

    ## [1] "male_gonad"
    ## [1] TRUE

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 493 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

![](../figures/manipulation/althings-12.png)

    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA         2       39      214       45    3     650
    ## 2 m.inc.d17      2        NA        1        0        0    1     577
    ## 3  m.inc.d3     39         1       NA        0        6    0     556
    ## 4  m.inc.d8    214         0        0       NA        7   19     679
    ## 5  m.inc.d9     45         0        6        7       NA    0     378
    ## 6      m.n2      3         1        0       19        0   NA     157
    ## 7   prolong    650       577      556      679      378  157      NA

![](../figures/manipulation/althings-13.png)

    ## Warning in MASS::cov.trob(data[, vars]): Probable convergence failure

    ## Warning in MASS::cov.trob(data[, vars]): Probable convergence failure

![](../figures/manipulation/althings-14.png)![](../figures/manipulation/althings-15.png)

    ## [1] "male_hypothalamus"
    ## [1] TRUE

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 13 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

![](../figures/manipulation/althings-16.png)

    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA         0        4        1     3575    1       0
    ## 2 m.inc.d17      0        NA        0      130     2721    0       0
    ## 3  m.inc.d3      4         0       NA      254     2974    0       0
    ## 4  m.inc.d8      1       130      254       NA     5669    0     129
    ## 5  m.inc.d9   3575      2721     2974     5669       NA 3405    2653
    ## 6      m.n2      1         0        0        0     3405   NA       0
    ## 7   prolong      0         0        0      129     2653    0      NA

![](../figures/manipulation/althings-17.png)![](../figures/manipulation/althings-18.png)![](../figures/manipulation/althings-19.png)

    ## [1] "male_pituitary"
    ## [1] TRUE

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 233 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

![](../figures/manipulation/althings-20.png)

    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA        87     1070     2153     1605  292       6
    ## 2 m.inc.d17     87        NA      250     1165      865    3      10
    ## 3  m.inc.d3   1070       250       NA      250      400  387     396
    ## 4  m.inc.d8   2153      1165      250       NA      625 1569    1140
    ## 5  m.inc.d9   1605       865      400      625       NA  579     570
    ## 6      m.n2    292         3      387     1569      579   NA     100
    ## 7   prolong      6        10      396     1140      570  100      NA

![](../figures/manipulation/althings-21.png)![](../figures/manipulation/althings-22.png)![](../figures/manipulation/althings-23.png)![](../figures/manipulation/althings-24.png)

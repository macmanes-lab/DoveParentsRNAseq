    library(tidyverse)
    library(DESeq2)
    library(cowplot)
    library(RColorBrewer)
    library(pheatmap)
    library(kableExtra)
    library(viridis)

    # load custom functions  
    source("../R/functions.R") 

    knitr::opts_chunk$set(fig.path = '../figures/manipulation/')

Manipulation data
=================

    # import "colData" which contains sample information and "countData" which contains read counts
    m.colData <- read.csv("../metadata/00_colData_manipluation.csv", header = T, row.names = 1)
    m.countData <- read.csv("../results/00_countData_manipluation.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)

    # set levels
    m.colData$treatment <- factor(m.colData$treatment, levels = 
                                  c("m.inc.d3",  "m.inc.d9",
                                    "m.inc.d17", "m.n2",
                                    "m.inc.d8", "prolong", "extend"))

    m.colData$sextissue <- as.factor(paste(m.colData$sex, m.colData$tissue, sep = "_"))

    m.colData$outcome <- ifelse(grepl("d3|d9|d17|n2", m.colData$treatment), "remove", 
                         ifelse(grepl("d8", m.colData$treatment),"early",
                         ifelse(grepl("prolong", m.colData$treatment),"prolong",
                         ifelse(grepl("extend", m.colData$treatment),"extend", NA))))
    m.colData$outcome <- factor(m.colData$outcome, levels = 
                                  c("remove",  "prolong",
                                    "early", "extend"))
    summary(m.colData[c(7,3,4,5,8,9)])

    ##           study         sex               tissue        treatment 
    ##  manipulation:411   female:208   gonad       :136   m.inc.d3 :60  
    ##                     male  :203   hypothalamus:138   m.inc.d9 :49  
    ##                                  pituitary   :137   m.inc.d17:63  
    ##                                                     m.n2     :59  
    ##                                                     m.inc.d8 :60  
    ##                                                     prolong  :60  
    ##                                                     extend   :60  
    ##                sextissue     outcome   
    ##  female_gonad       :69   remove :231  
    ##  female_hypothalamus:70   prolong: 60  
    ##  female_pituitary   :69   early  : 60  
    ##  male_gonad         :67   extend : 60  
    ##  male_hypothalamus  :68                
    ##  male_pituitary     :68                
    ## 

Write for loop to do this one for every tissue and for every treatment
======================================================================

    ## subset for test 

    m.colData <- m.colData %>%
      filter(sextissue == "female_pituitary") %>%
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
                                  c("m.inc.d3",  "m.inc.d9",
                                    "m.inc.d17", "m.n2",
                                    "m.inc.d8", "prolong", "extend"))
    dat$V2 <- factor(dat$V2, levels = 
                                  c("m.inc.d3",  "m.inc.d9",
                                    "m.inc.d17", "m.n2",
                                    "m.inc.d8", "prolong", "extend"))

    allcontrasts <- dat %>%
      ggplot( aes(V1, V2)) +
        geom_tile(aes(fill = V3)) +
        scale_fill_viridis(na.value="#FFFFFF00", 
                         limits = c(0, 6000),
                         breaks = c(0, 1000, 2000, 3000, 4000, 5000, 6000)) + 
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

    print(summary(aov(PC1 ~ treatment, data=pcadata)))
    print(TukeyHSD(aov(PC1 ~ treatment, data=pcadata), which = "treatment"))

    print(summary(aov(PC2 ~ treatment, data=pcadata))) 
    print(TukeyHSD(aov(PC2 ~ treatment, data=pcadata), which = "treatment")) 

    pca34 <- ggplot(pcadata, aes(PC3, PC4,color = treatment)) + 
      geom_point(size = 2, alpha = 1) +
      stat_ellipse(type = "t") +
      xlab(paste0("PC3: ", percentVar[3],"% variance")) +
      ylab(paste0("PC4: ", percentVar[4],"% variance")) +
      theme_cowplot(font_size = 8, line_size = 0.25) +
      labs(subtitle = eachgroup)
    #print(pca34)

    # see http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix
    sampleDists <- dist(t(assay(vsd)))

    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- NULL
    colnames(sampleDistMatrix) <- colData$treatment
    colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)

    #pheatmap(sampleDistMatrix,
    #         clustering_distance_rows=sampleDists,
    #         clustering_distance_cols=sampleDists,
    #         col=colors,
     #        fontsize = 6, 
     #        main = eachgroup)

    # Multi-factor designs
    # revals is my function easily calculate multiple multi factor designs
    # it is built around the documentaion from DESeq which looks more simply like
    # res <- results(dds,contrast=c("columnname", "factor1", "factor2"))
    # then I extract the p-values from within the results for showing number of DEGes
    resvals <- function(contrastvector, mypval){
      res <- results(dds, contrast = c(contrastvector[1],contrastvector[2],contrastvector[3]), independentFiltering = T)
      sumpvalue <- sum(res$pvalue < mypval, na.rm = TRUE)
      #print(sumpvalue)
      sumpadj <- sum(res$padj < mypval, na.rm = TRUE)
      print(sumpadj)
      vals <- cbind(res$pvalue, res$padj)
      pvalcolname <- as.character(paste("pval",contrastvector[1],contrastvector[2],contrastvector[3], sep=""))
      padjcolname <- as.character(paste("padj",contrastvector[1],contrastvector[2],contrastvector[3], sep=""))
      colnames(vals) <- c(pvalcolname, padjcolname)
      return(vals)
    }
    ## DEG by contrasts at 0.1 pvalue
    contrast1 <- resvals(contrastvector = c('treatment', 'm.inc.d3', 'm.inc.d9'), mypval = 0.01) #12
    contrast2 <- resvals(contrastvector = c('treatment', 'm.inc.d3', 'm.inc.d17'), mypval = 0.01) #73
    contrast3 <- resvals(contrastvector = c('treatment', 'm.inc.d3', 'm.n2'), mypval = 0.01) #284
    contrast4 <- resvals(contrastvector = c('treatment', 'm.inc.d3', 'm.inc.d8'), mypval = 0.01) #15
    contrast5 <- resvals(contrastvector = c('treatment', 'm.inc.d3', 'prolong'), mypval = 0.01) #555
    contrast6 <- resvals(contrastvector = c('treatment', 'm.inc.d3', 'extend'), mypval = 0.01) #1109
    contrast7 <- resvals(contrastvector = c('treatment', 'm.inc.d9', 'm.inc.d17'), mypval = 0.01) #280
    contrast8 <- resvals(contrastvector = c('treatment', 'm.inc.d9', 'm.n2'), mypval = 0.01) #440
    contrast9 <- resvals(contrastvector = c('treatment', 'm.inc.d9', 'm.inc.d8'), mypval = 0.01) #484
    contrast10 <- resvals(contrastvector = c('treatment', 'm.inc.d9', 'prolong'), mypval = 0.01) #529
    contrast11 <- resvals(contrastvector = c('treatment', 'm.inc.d9', 'extend'), mypval = 0.01) #1177

    contrast12 <- resvals(contrastvector = c('treatment', 'm.inc.d17', 'm.n2'), mypval = 0.01) #5
    contrast13 <- resvals(contrastvector = c('treatment', 'm.inc.d17', 'm.inc.d8'), mypval = 0.01) #553
    contrast14 <- resvals(contrastvector = c('treatment', 'm.inc.d17', 'prolong'), mypval = 0.01) #339
    contrast15 <- resvals(contrastvector = c('treatment', 'm.inc.d17', 'extend'), mypval = 0.01) #782
    contrast16 <- resvals(contrastvector = c('treatment', 'm.n2', 'm.inc.d8'), mypval = 0.01) #1123
    contrast17 <- resvals(contrastvector = c('treatment', 'm.n2', 'prolong'), mypval = 0.01) #265
    contrast18 <- resvals(contrastvector = c('treatment', 'm.n2', 'extend'), mypval = 0.01) #552
    contrast19 <- resvals(contrastvector = c('treatment', 'm.inc.d8', 'prolong'), mypval = 0.01) #828
    contrast20 <- resvals(contrastvector = c('treatment', 'm.inc.d8', 'extend'), mypval = 0.01) #1110
    contrast21 <- resvals(contrastvector = c('treatment', 'prolong', 'extend'), mypval = 0.01) #64


    DEGs <- assay(vsd)
    DEGs <- cbind(DEGs, contrast1, contrast2, contrast3, contrast4, contrast5,
                  contrast6, contrast7, contrast8, contrast9, contrast10,
                  contrast11, contrast12, contrast13, contrast14, contrast15,
                  contrast16, contrast17, contrast18, contrast19, contrast20,contrast21)
    DEGs <- as.data.frame(DEGs) # convert matrix to dataframe
    DEGs$rownames <- rownames(DEGs)  # add the rownames to the dataframe

    #DEGs[is.na(DEGs)] <-  1


    DEGs$padjmin <- with(DEGs, pmin(pvaltreatmentm.inc.d3m.inc.d9, 
      padjtreatmentm.inc.d3m.inc.d9,               
      pvaltreatmentm.inc.d3m.inc.d17, padjtreatmentm.inc.d3m.inc.d17,              
      pvaltreatmentm.inc.d3m.n2, padjtreatmentm.inc.d3m.n2,                  
      pvaltreatmentm.inc.d3m.inc.d8, padjtreatmentm.inc.d3m.inc.d8,               
      pvaltreatmentm.inc.d3prolong, padjtreatmentm.inc.d3prolong,               
      pvaltreatmentm.inc.d3extend, padjtreatmentm.inc.d3extend,                
      pvaltreatmentm.inc.d9m.inc.d17, padjtreatmentm.inc.d9m.inc.d17,             
      pvaltreatmentm.inc.d9m.n2, padjtreatmentm.inc.d9m.n2,                  
      pvaltreatmentm.inc.d9m.inc.d8, padjtreatmentm.inc.d9m.inc.d8,              
      pvaltreatmentm.inc.d9prolong, padjtreatmentm.inc.d9prolong,               
      pvaltreatmentm.inc.d9extend, padjtreatmentm.inc.d9extend,                
      pvaltreatmentm.inc.d17m.n2, padjtreatmentm.inc.d17m.n2,                 
      pvaltreatmentm.inc.d17m.inc.d8, padjtreatmentm.inc.d17m.inc.d8,             
      pvaltreatmentm.inc.d17prolong, padjtreatmentm.inc.d17prolong,
      pvaltreatmentm.inc.d17extend,padjtreatmentm.inc.d17extend, 
      pvaltreatmentm.n2m.inc.d8, padjtreatmentm.n2m.inc.d8,                  
      pvaltreatmentm.n2prolong, padjtreatmentm.n2prolong,                   
      pvaltreatmentm.n2extend, padjtreatmentm.n2extend,                    
      pvaltreatmentm.inc.d8prolong, padjtreatmentm.inc.d8prolong,               
      pvaltreatmentm.inc.d8extend, padjtreatmentm.inc.d8extend,                
      pvaltreatmentprolongextend, padjtreatmentprolongextend)) 

    sigDEGs <- DEGs %>% arrange(padjmin)
    sigDEGs <- head(sigDEGs,30)
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
             annotation_col=anndf)

    pheatmap(sigDEGs, kmeans_k = 5,
             show_rownames = F, show_colnames = F,
             color = viridis(30),
             breaks=myBreaks,
             annotation_col=anndf)
    }

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

    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA      2014     2670     2837     2842 1381     284
    ## 2 m.inc.d17   2014        NA      317     1978     1525   13    1272
    ## 3  m.inc.d3   2670       317       NA      187      522 1249    1700
    ## 4  m.inc.d8   2837      1978      187       NA     1156 3031    2380
    ## 5  m.inc.d9   2842      1525      522     1156       NA 1621    1547
    ## 6      m.n2   1381        13     1249     3031     1621   NA     943
    ## 7   prolong    284      1272     1700     2380     1547  943      NA

![](../figures/manipulation/althings-1.png)![](../figures/manipulation/althings-2.png)

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment    6 1116.2  186.04   20.54 4.41e-13 ***
    ## Residuals   62  561.4    9.06                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                           diff        lwr        upr     p adj
    ## m.inc.d9-m.inc.d3    1.8634057  -2.486171  6.2129819 0.8468409
    ## m.inc.d17-m.inc.d3   1.4336881  -2.572847  5.4402236 0.9286328
    ## m.n2-m.inc.d3        0.7097154  -3.391104  4.8105352 0.9983377
    ## m.inc.d8-m.inc.d3   -3.8655271  -7.966347  0.2352927 0.0772733
    ## prolong-m.inc.d3    -6.3960176 -10.496837 -2.2951978 0.0002398
    ## extend-m.inc.d3     -9.2890257 -13.389846 -5.1882059 0.0000001
    ## m.inc.d17-m.inc.d9  -0.4297175  -4.690518  3.8310825 0.9999261
    ## m.n2-m.inc.d9       -1.1536903  -5.503267  3.1958859 0.9832727
    ## m.inc.d8-m.inc.d9   -5.7289328 -10.078509 -1.3793565 0.0029668
    ## prolong-m.inc.d9    -8.2594233 -12.609000 -3.9098470 0.0000052
    ## extend-m.inc.d9    -11.1524314 -15.502008 -6.8028551 0.0000000
    ## m.n2-m.inc.d17      -0.7239728  -4.730508  3.2825627 0.9978832
    ## m.inc.d8-m.inc.d17  -5.2992152  -9.305751 -1.2926797 0.0028098
    ## prolong-m.inc.d17   -7.8297058 -11.836241 -3.8231702 0.0000027
    ## extend-m.inc.d17   -10.7227138 -14.729249 -6.7161783 0.0000000
    ## m.inc.d8-m.n2       -4.5752424  -8.676062 -0.4744226 0.0192398
    ## prolong-m.n2        -7.1057330 -11.206553 -3.0049131 0.0000351
    ## extend-m.n2         -9.9987411 -14.099561 -5.8979212 0.0000000
    ## prolong-m.inc.d8    -2.5304905  -6.631310  1.5703293 0.5007634
    ## extend-m.inc.d8     -5.4234986  -9.524318 -1.3226788 0.0028127
    ## extend-prolong      -2.8930081  -6.993828  1.2078117 0.3376042
    ## 
    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment    6  526.8   87.79   7.915 2.38e-06 ***
    ## Residuals   62  687.7   11.09                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                           diff          lwr       upr     p adj
    ## m.inc.d9-m.inc.d3  -1.95737666  -6.77119833  2.856445 0.8760019
    ## m.inc.d17-m.inc.d3  4.34979304  -0.08437393  8.783960 0.0579389
    ## m.n2-m.inc.d3       5.09556144   0.55704685  9.634076 0.0180968
    ## m.inc.d8-m.inc.d3  -2.04335880  -6.58187340  2.495156 0.8143743
    ## prolong-m.inc.d3    3.24139278  -1.29712181  7.779907 0.3231200
    ## extend-m.inc.d3     3.63231953  -0.90619507  8.170834 0.2003816
    ## m.inc.d17-m.inc.d9  6.30716970   1.59159963 11.022740 0.0024250
    ## m.n2-m.inc.d9       7.05293810   2.23911643 11.866760 0.0006577
    ## m.inc.d8-m.inc.d9  -0.08598215  -4.89980382  4.727840 1.0000000
    ## prolong-m.inc.d9    5.19876944   0.38494777 10.012591 0.0261070
    ## extend-m.inc.d9     5.58969618   0.77587451 10.403518 0.0128838
    ## m.n2-m.inc.d17      0.74576840  -3.68839858  5.179935 0.9985849
    ## m.inc.d8-m.inc.d17 -6.39315185 -10.82731883 -1.958985 0.0008394
    ## prolong-m.inc.d17  -1.10840026  -5.54256724  3.325767 0.9876982
    ## extend-m.inc.d17   -0.71747352  -5.15164049  3.716693 0.9988628
    ## m.inc.d8-m.n2      -7.13892025 -11.67743484 -2.600406 0.0002076
    ## prolong-m.n2       -1.85416866  -6.39268326  2.684346 0.8735637
    ## extend-m.n2        -1.46324191  -6.00175651  3.075273 0.9559841
    ## prolong-m.inc.d8    5.28475159   0.74623699  9.823266 0.0125139
    ## extend-m.inc.d8     5.67567833   1.13716374 10.214193 0.0056456
    ## extend-prolong      0.39092675  -4.14758785  4.929441 0.9999708
    ## 
    ## [1] 12
    ## [1] 73
    ## [1] 348
    ## [1] 15
    ## [1] 555
    ## [1] 1109
    ## [1] 280
    ## [1] 440
    ## [1] 241
    ## [1] 529
    ## [1] 1177
    ## [1] 5
    ## [1] 553
    ## [1] 339
    ## [1] 782
    ## [1] 1123
    ## [1] 265
    ## [1] 552
    ## [1] 828
    ## [1] 1010
    ## [1] 62

![](../figures/manipulation/althings-3.png)![](../figures/manipulation/althings-4.png)

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
                                  design = ~ outcome )
      dds <- dds[ rowSums(counts(dds)) > 2, ] ## pre-filter genes 
      dds <- DESeq(dds) # Differential expression analysis
      vsd <- vst(dds, blind=FALSE) # variance stabilized 
    head(vsd)

    numDEGs <- function(group1, group2){
      res <- results(dds, contrast = c("outcome", group1, group2), independentFiltering = T)
      sumpadj <- sum(res$padj < 0.1, na.rm = TRUE)
      return(sumpadj)}

    #create list of groups
    a <- levels(colData$outcome)
    b <- levels(colData$outcome)

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
                                  c("remove",  "prolong",
                                    "early", "extend"))
    dat$V2 <- factor(dat$V2, levels = 
                                  c("remove",  "prolong",
                                    "early", "extend"))

    allcontrasts <- dat %>%
      ggplot( aes(V1, V2)) +
        geom_tile(aes(fill = V3)) +
        scale_fill_viridis(na.value="#FFFFFF00", 
                         limits = c(0, 6000),
                         breaks = c(0, 1000, 2000, 3000, 4000, 5000, 6000)) + 
        xlab(" ") + ylab("Timepoint") +
        labs(fill = "# of DEGs",
             subtitle = eachgroup)
    plot(allcontrasts)

    # create the dataframe using my function pcadataframe
    pcadata <- pcadataframe(vsd, intgroup=c("outcome"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))
    percentVar

    pca12 <- ggplot(pcadata, aes(PC1, PC2,color = outcome)) + 
      geom_point(size = 2, alpha = 1) +
      stat_ellipse(type = "t") +
      xlab(paste0("PC1: ", percentVar[1],"% variance")) +
      ylab(paste0("PC2: ", percentVar[2],"% variance")) +
      theme_cowplot(font_size = 8, line_size = 0.25) +
      labs(subtitle = eachgroup)
    print(pca12)

    print(summary(aov(PC1 ~ outcome, data=pcadata)))
    print(TukeyHSD(aov(PC1 ~ outcome, data=pcadata), which = "outcome"))

    print(summary(aov(PC2 ~ outcome, data=pcadata))) 
    print(TukeyHSD(aov(PC2 ~ outcome, data=pcadata), which = "outcome")) 

    pca34 <- ggplot(pcadata, aes(PC3, PC4,color = outcome)) + 
      geom_point(size = 2, alpha = 1) +
      stat_ellipse(type = "t") +
      xlab(paste0("PC3: ", percentVar[3],"% variance")) +
      ylab(paste0("PC4: ", percentVar[4],"% variance")) +
      theme_cowplot(font_size = 8, line_size = 0.25) +
      labs(subtitle = eachgroup)
    #print(pca34)

    }

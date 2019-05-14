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
    summary(m.colData[c(7,3,4,5,9)])

    ##           study         sex               tissue        treatment 
    ##  manipulation:411   female:208   gonad       :136   m.inc.d3 :60  
    ##                     male  :203   hypothalamus:138   m.inc.d9 :49  
    ##                                  pituitary   :137   m.inc.d17:63  
    ##                                                     m.n2     :59  
    ##                                                     m.inc.d8 :60  
    ##                                                     prolong  :60  
    ##                                                     extend   :60  
    ##     outcome   
    ##  remove :231  
    ##  prolong: 60  
    ##  early  : 60  
    ##  extend : 60  
    ##               
    ##               
    ## 

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

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    6    727   121.2    0.61  0.721
    ## Residuals   62  12320   198.7               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                            diff        lwr      upr     p adj
    ## m.inc.d9-m.inc.d3   8.288419936 -12.086705 28.66355 0.8757811
    ## m.inc.d17-m.inc.d3  9.170693472  -9.597494 27.93888 0.7501910
    ## m.n2-m.inc.d3       6.181970794 -13.027882 25.39182 0.9563682
    ## m.inc.d8-m.inc.d3   3.985928139 -15.223924 23.19578 0.9954565
    ## prolong-m.inc.d3    3.172718959 -16.037133 22.38257 0.9987228
    ## extend-m.inc.d3     9.176996148 -10.032856 28.38685 0.7692394
    ## m.inc.d17-m.inc.d9  0.882273535 -19.076989 20.84154 0.9999994
    ## m.n2-m.inc.d9      -2.106449142 -22.481574 18.26868 0.9999145
    ## m.inc.d8-m.inc.d9  -4.302491797 -24.677617 16.07263 0.9950010
    ## prolong-m.inc.d9   -5.115700977 -25.490826 15.25942 0.9874108
    ## extend-m.inc.d9     0.888576212 -19.486549 21.26370 0.9999995
    ## m.n2-m.inc.d17     -2.988722677 -21.756910 15.77946 0.9989612
    ## m.inc.d8-m.inc.d17 -5.184765333 -23.952953 13.58342 0.9794158
    ## prolong-m.inc.d17  -5.997974513 -24.766162 12.77021 0.9577923
    ## extend-m.inc.d17    0.006302676 -18.761885 18.77449 1.0000000
    ## m.inc.d8-m.n2      -2.196042656 -21.405895 17.01381 0.9998461
    ## prolong-m.n2       -3.009251835 -22.219104 16.20060 0.9990537
    ## extend-m.n2         2.995025353 -16.214827 22.20488 0.9990788
    ## prolong-m.inc.d8   -0.813209180 -20.023062 18.39664 0.9999996
    ## extend-m.inc.d8     5.191068009 -14.018784 24.40092 0.9815988
    ## extend-prolong      6.004277189 -13.205575 25.21413 0.9620632
    ## 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    6    358   59.62   0.747  0.614
    ## Residuals   62   4951   79.85               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                          diff        lwr       upr     p adj
    ## m.inc.d9-m.inc.d3   3.2502224  -9.666037 16.166482 0.9872630
    ## m.inc.d17-m.inc.d3  4.6824297  -7.215155 16.580014 0.8918193
    ## m.n2-m.inc.d3       2.9309140  -9.246652 15.108480 0.9899169
    ## m.inc.d8-m.inc.d3   2.2245065  -9.953059 14.402072 0.9977509
    ## prolong-m.inc.d3   -2.6957459 -14.873312  9.481820 0.9935459
    ## extend-m.inc.d3     1.4553941 -10.722172 13.632960 0.9998007
    ## m.inc.d17-m.inc.d9  1.4322073 -11.220427 14.084842 0.9998547
    ## m.n2-m.inc.d9      -0.3193084 -13.235568 12.596951 1.0000000
    ## m.inc.d8-m.inc.d9  -1.0257159 -13.941975 11.890543 0.9999819
    ## prolong-m.inc.d9   -5.9459683 -18.862227  6.970291 0.7982752
    ## extend-m.inc.d9    -1.7948282 -14.711087 11.121431 0.9995232
    ## m.n2-m.inc.d17     -1.7515157 -13.649100 10.146069 0.9993360
    ## m.inc.d8-m.inc.d17 -2.4579232 -14.355508  9.439661 0.9955635
    ## prolong-m.inc.d17  -7.3781756 -19.275760  4.519409 0.4947226
    ## extend-m.inc.d17   -3.2270355 -15.124620  8.670549 0.9812470
    ## m.inc.d8-m.n2      -0.7064075 -12.883973 11.471158 0.9999972
    ## prolong-m.n2       -5.6266599 -17.804226  6.550906 0.7955043
    ## extend-m.n2        -1.4755198 -13.653086 10.702046 0.9997841
    ## prolong-m.inc.d8   -4.9202524 -17.097818  7.257314 0.8792194
    ## extend-m.inc.d8    -0.7691124 -12.946678 11.408454 0.9999954
    ## extend-prolong      4.1511401  -8.026426 16.328706 0.9428206
    ## 
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

![](../figures/manipulation/althings-2.png)

    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA       336     1051        0      609   16       6
    ## 2 m.inc.d17    336        NA        0      697      448   15      45
    ## 3  m.inc.d3   1051         0       NA     1630     1276  239     116
    ## 4  m.inc.d8      0       697     1630       NA      192    5      20
    ## 5  m.inc.d9    609       448     1276      192       NA   88     935
    ## 6      m.n2     16        15      239        5       88   NA     142
    ## 7   prolong      6        45      116       20      935  142      NA

![](../figures/manipulation/althings-3.png)

    ##             Df Sum Sq Mean Sq F value Pr(>F)   
    ## treatment    6  366.6   61.10    3.71 0.0032 **
    ## Residuals   63 1037.6   16.47                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                          diff         lwr        upr     p adj
    ## m.inc.d9-m.inc.d3  -6.5670766 -12.2460429 -0.8881103 0.0134194
    ## m.inc.d17-m.inc.d3 -0.1061703  -5.5065842  5.2942436 1.0000000
    ## m.n2-m.inc.d3      -3.7030360  -9.2305358  1.8244639 0.4004577
    ## m.inc.d8-m.inc.d3  -5.0656415 -10.5931413  0.4618584 0.0934475
    ## prolong-m.inc.d3   -1.8940812  -7.4215810  3.6334187 0.9416039
    ## extend-m.inc.d3    -4.0127781  -9.5402780  1.5147217 0.3045198
    ## m.inc.d17-m.inc.d9  6.4609063   0.9055599 12.0162527 0.0126456
    ## m.n2-m.inc.d9       2.8640406  -2.8149256  8.5430069 0.7222966
    ## m.inc.d8-m.inc.d9   1.5014351  -4.1775311  7.1804014 0.9836074
    ## prolong-m.inc.d9    4.6729954  -1.0059708 10.3519617 0.1747660
    ## extend-m.inc.d9     2.5542985  -3.1246678  8.2332647 0.8154718
    ## m.n2-m.inc.d17     -3.5968656  -8.9972795  1.8035483 0.4076012
    ## m.inc.d8-m.inc.d17 -4.9594711 -10.3598851  0.4409428 0.0921851
    ## prolong-m.inc.d17  -1.7879108  -7.1883248  3.6125031 0.9502585
    ## extend-m.inc.d17   -3.9066078  -9.3070217  1.4938061 0.3086197
    ## m.inc.d8-m.n2      -1.3626055  -6.8901054  4.1648944 0.9886051
    ## prolong-m.n2        1.8089548  -3.7185451  7.3364547 0.9528964
    ## extend-m.n2        -0.3097422  -5.8372420  5.2177577 0.9999977
    ## prolong-m.inc.d8    3.1715603  -2.3559396  8.6990602 0.5875109
    ## extend-m.inc.d8     1.0528633  -4.4746365  6.5803632 0.9971735
    ## extend-prolong     -2.1186970  -7.6461968  3.4088029 0.9036262
    ## 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    6   34.8   5.796   0.325  0.921
    ## Residuals   63 1122.0  17.809               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                           diff       lwr      upr     p adj
    ## m.inc.d9-m.inc.d3  -1.69811389 -7.603552 4.207324 0.9749032
    ## m.inc.d17-m.inc.d3 -1.85827204 -7.474049 3.757505 0.9503765
    ## m.n2-m.inc.d3      -0.82437632 -6.572307 4.923555 0.9994307
    ## m.inc.d8-m.inc.d3  -2.05879927 -7.806730 3.689132 0.9285318
    ## prolong-m.inc.d3   -0.61830736 -6.366238 5.129624 0.9998926
    ## extend-m.inc.d3    -0.79932757 -6.547259 4.948603 0.9995232
    ## m.inc.d17-m.inc.d9 -0.16015815 -5.937046 5.616730 1.0000000
    ## m.n2-m.inc.d9       0.87373757 -5.031700 6.779175 0.9993196
    ## m.inc.d8-m.inc.d9  -0.36068538 -6.266123 5.544752 0.9999962
    ## prolong-m.inc.d9    1.07980653 -4.825631 6.985244 0.9977479
    ## extend-m.inc.d9     0.89878632 -5.006651 6.804224 0.9992004
    ## m.n2-m.inc.d17      1.03389573 -4.581881 6.649673 0.9976604
    ## m.inc.d8-m.inc.d17 -0.20052723 -5.816304 5.415250 0.9999998
    ## prolong-m.inc.d17   1.23996468 -4.375812 6.855742 0.9936595
    ## extend-m.inc.d17    1.05894447 -4.556833 6.674721 0.9973273
    ## m.inc.d8-m.n2      -1.23442295 -6.982354 4.513508 0.9945425
    ## prolong-m.n2        0.20606896 -5.541862 5.954000 0.9999998
    ## extend-m.n2         0.02504874 -5.722882 5.772980 1.0000000
    ## prolong-m.inc.d8    1.44049191 -4.307439 7.188423 0.9875781
    ## extend-m.inc.d8     1.25947170 -4.488459 7.007403 0.9939159
    ## extend-prolong     -0.18102021 -5.928951 5.566911 0.9999999
    ## 
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

![](../figures/manipulation/althings-4.png)

    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA      2014     2670     2837     2842 1381     284
    ## 2 m.inc.d17   2014        NA      317     1978     1525   13    1272
    ## 3  m.inc.d3   2670       317       NA      187      522 1249    1700
    ## 4  m.inc.d8   2837      1978      187       NA     1156 3031    2380
    ## 5  m.inc.d9   2842      1525      522     1156       NA 1621    1547
    ## 6      m.n2   1381        13     1249     3031     1621   NA     943
    ## 7   prolong    284      1272     1700     2380     1547  943      NA

![](../figures/manipulation/althings-5.png)

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

![](../figures/manipulation/althings-6.png)

    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA         2       39      214       45    3     650
    ## 2 m.inc.d17      2        NA        1        0        0    1     577
    ## 3  m.inc.d3     39         1       NA        0        6    0     556
    ## 4  m.inc.d8    214         0        0       NA        7   19     679
    ## 5  m.inc.d9     45         0        6        7       NA    0     378
    ## 6      m.n2      3         1        0       19        0   NA     157
    ## 7   prolong    650       577      556      679      378  157      NA

![](../figures/manipulation/althings-7.png)

    ## Warning in MASS::cov.trob(data[, vars]): Probable convergence failure

    ## Warning in MASS::cov.trob(data[, vars]): Probable convergence failure

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    6    521   86.86   1.011  0.427
    ## Residuals   60   5155   85.91               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                          diff        lwr       upr     p adj
    ## m.inc.d9-m.inc.d3  -0.4910682 -13.902995 12.920858 0.9999998
    ## m.inc.d17-m.inc.d3  0.6543716 -11.990514 13.299257 0.9999986
    ## m.n2-m.inc.d3       0.8730649 -12.118320 13.864450 0.9999932
    ## m.inc.d8-m.inc.d3   1.1960271 -11.448858 13.840912 0.9999488
    ## prolong-m.inc.d3    8.2976615  -4.347224 20.942547 0.4243613
    ## extend-m.inc.d3     2.0009188 -10.643967 14.645804 0.9989889
    ## m.inc.d17-m.inc.d9  1.1454399 -12.266486 14.557366 0.9999720
    ## m.n2-m.inc.d9       1.3641331 -12.374961 15.103228 0.9999319
    ## m.inc.d8-m.inc.d9   1.6870953 -11.724831 15.099022 0.9997292
    ## prolong-m.inc.d9    8.7887298  -4.623197 22.200656 0.4260734
    ## extend-m.inc.d9     2.4919870 -10.919939 15.903913 0.9975063
    ## m.n2-m.inc.d17      0.2186933 -12.772691 13.210078 1.0000000
    ## m.inc.d8-m.inc.d17  0.5416554 -12.103230 13.186541 0.9999995
    ## prolong-m.inc.d17   7.6432899  -5.001595 20.288175 0.5246010
    ## extend-m.inc.d17    1.3465471 -11.298338 13.991433 0.9998974
    ## m.inc.d8-m.n2       0.3229622 -12.668423 13.314347 1.0000000
    ## prolong-m.n2        7.4245966  -5.566788 20.415981 0.5903757
    ## extend-m.n2         1.1278539 -11.863531 14.119239 0.9999691
    ## prolong-m.inc.d8    7.1016345  -5.543251 19.746520 0.6100735
    ## extend-m.inc.d8     0.8048917 -11.839994 13.449777 0.9999951
    ## extend-prolong     -6.2967428 -18.941628  6.348143 0.7324015
    ## 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    6  296.6   49.44   1.543   0.18
    ## Residuals   60 1922.2   32.04               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                           diff        lwr       upr     p adj
    ## m.inc.d9-m.inc.d3  -2.36319936 -10.553441  5.827042 0.9741974
    ## m.inc.d17-m.inc.d3  3.39787109  -4.323962 11.119705 0.8291233
    ## m.n2-m.inc.d3       2.07288480  -5.860545 10.006315 0.9844001
    ## m.inc.d8-m.inc.d3  -2.28577754 -10.007611  5.436056 0.9707073
    ## prolong-m.inc.d3    1.03079949  -6.691034  8.752633 0.9996182
    ## extend-m.inc.d3     2.61555489  -5.106279 10.337388 0.9441258
    ## m.inc.d17-m.inc.d9  5.76107045  -2.429171 13.951312 0.3402203
    ## m.n2-m.inc.d9       4.43608416  -3.953948 12.826117 0.6746456
    ## m.inc.d8-m.inc.d9   0.07742183  -8.112819  8.267663 1.0000000
    ## prolong-m.inc.d9    3.39399885  -4.796242 11.584240 0.8652873
    ## extend-m.inc.d9     4.97875425  -3.211487 13.168995 0.5177934
    ## m.n2-m.inc.d17     -1.32498629  -9.258416  6.608443 0.9986278
    ## m.inc.d8-m.inc.d17 -5.68364862 -13.405482  2.038185 0.2876953
    ## prolong-m.inc.d17  -2.36707160 -10.088905  5.354762 0.9652369
    ## extend-m.inc.d17   -0.78231620  -8.504150  6.939517 0.9999234
    ## m.inc.d8-m.n2      -4.35866234 -12.292092  3.574767 0.6342868
    ## prolong-m.n2       -1.04208531  -8.975515  6.891344 0.9996522
    ## extend-m.n2         0.54267009  -7.390760  8.476100 0.9999925
    ## prolong-m.inc.d8    3.31657702  -4.405256 11.038411 0.8445221
    ## extend-m.inc.d8     4.90133242  -2.820501 12.623166 0.4651975
    ## extend-prolong      1.58475540  -6.137078  9.306589 0.9956849
    ## 
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

![](../figures/manipulation/althings-8.png)

    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA         0        4        1     3575    1       0
    ## 2 m.inc.d17      0        NA        0      130     2721    0       0
    ## 3  m.inc.d3      4         0       NA      254     2974    0       0
    ## 4  m.inc.d8      1       130      254       NA     5669    0     129
    ## 5  m.inc.d9   3575      2721     2974     5669       NA 3405    2653
    ## 6      m.n2      1         0        0        0     3405   NA       0
    ## 7   prolong      0         0        0      129     2653    0      NA

![](../figures/manipulation/althings-9.png)

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
    ##                           diff        lwr        upr     p adj
    ## m.inc.d9-m.inc.d3    8.4818829   1.991134 14.9726314 0.0033023
    ## m.inc.d17-m.inc.d3   0.2285568  -5.890980  6.3480931 0.9999998
    ## m.n2-m.inc.d3       -1.2802542  -7.399791  4.8392822 0.9952286
    ## m.inc.d8-m.inc.d3   -2.9425407  -9.062077  3.1769957 0.7633728
    ## prolong-m.inc.d3     1.0518632  -5.067673  7.1713995 0.9983936
    ## extend-m.inc.d3     -0.3208496  -6.440386  5.7986867 0.9999985
    ## m.inc.d17-m.inc.d9  -8.2533261 -14.744075 -1.7625776 0.0046418
    ## m.n2-m.inc.d9       -9.7621371 -16.252886 -3.2713886 0.0004413
    ## m.inc.d8-m.inc.d9  -11.4244236 -17.915172 -4.9336751 0.0000263
    ## prolong-m.inc.d9    -7.4300197 -13.920768 -0.9392713 0.0149496
    ## extend-m.inc.d9     -8.8027325 -15.293481 -2.3119841 0.0020265
    ## m.n2-m.inc.d17      -1.5088110  -7.628347  4.6107254 0.9885084
    ## m.inc.d8-m.inc.d17  -3.1710975  -9.290634  2.9484389 0.6953497
    ## prolong-m.inc.d17    0.8233064  -5.296230  6.9428427 0.9996024
    ## extend-m.inc.d17    -0.5494065  -6.668943  5.5701299 0.9999625
    ## m.inc.d8-m.n2       -1.6622865  -7.781823  4.4572498 0.9810358
    ## prolong-m.n2         2.3321173  -3.787419  8.4516537 0.9054801
    ## extend-m.n2          0.9594045  -5.160132  7.0789409 0.9990452
    ## prolong-m.inc.d8     3.9944038  -2.125133 10.1139402 0.4314036
    ## extend-m.inc.d8      2.6216910  -3.497845  8.7412274 0.8464640
    ## extend-prolong      -1.3727128  -7.492249  4.7468235 0.9930426
    ## 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    6   73.2   12.21   1.019  0.422
    ## Residuals   61  730.9   11.98               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                          diff       lwr      upr     p adj
    ## m.inc.d9-m.inc.d3  -0.4538550 -5.460023 4.552313 0.9999603
    ## m.inc.d17-m.inc.d3 -0.3206954 -5.040556 4.399165 0.9999928
    ## m.n2-m.inc.d3      -1.5182015 -6.238062 3.201659 0.9563218
    ## m.inc.d8-m.inc.d3   0.2628505 -4.457010 4.982711 0.9999978
    ## prolong-m.inc.d3   -2.6092873 -7.329148 2.110573 0.6280903
    ## extend-m.inc.d3    -2.0837008 -6.803561 2.636160 0.8273739
    ## m.inc.d17-m.inc.d9  0.1331596 -4.873008 5.139328 1.0000000
    ## m.n2-m.inc.d9      -1.0643465 -6.070515 3.941822 0.9947914
    ## m.inc.d8-m.inc.d9   0.7167055 -4.289463 5.722874 0.9994316
    ## prolong-m.inc.d9   -2.1554324 -7.161600 2.850736 0.8434127
    ## extend-m.inc.d9    -1.6298458 -6.636014 3.376322 0.9537459
    ## m.n2-m.inc.d17     -1.1975061 -5.917367 3.522354 0.9866549
    ## m.inc.d8-m.inc.d17  0.5835459 -4.136315 5.303406 0.9997562
    ## prolong-m.inc.d17  -2.2885920 -7.008452 2.431269 0.7563124
    ## extend-m.inc.d17   -1.7630054 -6.482866 2.956855 0.9133821
    ## m.inc.d8-m.n2       1.7810520 -2.938809 6.500912 0.9094414
    ## prolong-m.n2       -1.0910859 -5.810946 3.628775 0.9918272
    ## extend-m.n2        -0.5654993 -5.285360 4.154361 0.9997968
    ## prolong-m.inc.d8   -2.8721378 -7.591998 1.847723 0.5171058
    ## extend-m.inc.d8    -2.3465512 -7.066412 2.373309 0.7343801
    ## extend-prolong      0.5255866 -4.194274 5.245447 0.9998674
    ## 
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

![](../figures/manipulation/althings-10.png)

    ##          V1 extend m.inc.d17 m.inc.d3 m.inc.d8 m.inc.d9 m.n2 prolong
    ## 1    extend     NA        87     1070     2153     1605  292       6
    ## 2 m.inc.d17     87        NA      250     1165      865    3      10
    ## 3  m.inc.d3   1070       250       NA      250      400  387     396
    ## 4  m.inc.d8   2153      1165      250       NA      625 1569    1140
    ## 5  m.inc.d9   1605       865      400      625       NA  579     570
    ## 6      m.n2    292         3      387     1569      579   NA     100
    ## 7   prolong      6        10      396     1140      570  100      NA

![](../figures/manipulation/althings-11.png)![](../figures/manipulation/althings-12.png)

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    6  546.2   91.03   1.873    0.1
    ## Residuals   61 2964.7   48.60               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                          diff        lwr       upr     p adj
    ## m.inc.d9-m.inc.d3  -5.0571546 -15.139226  5.024917 0.7262959
    ## m.inc.d17-m.inc.d3  1.1121817  -8.393286 10.617649 0.9998230
    ## m.n2-m.inc.d3       0.1740609  -9.331407  9.679529 1.0000000
    ## m.inc.d8-m.inc.d3  -6.1437379 -15.649206  3.361730 0.4434990
    ## prolong-m.inc.d3    0.4604334  -9.045034  9.965901 0.9999990
    ## extend-m.inc.d3     1.4539679  -8.051500 10.959436 0.9991701
    ## m.inc.d17-m.inc.d9  6.1693362  -3.912735 16.251407 0.5103964
    ## m.n2-m.inc.d9       5.2312154  -4.850856 15.313287 0.6940780
    ## m.inc.d8-m.inc.d9  -1.0865834 -11.168654  8.995488 0.9998904
    ## prolong-m.inc.d9    5.5175879  -4.564483 15.599659 0.6390724
    ## extend-m.inc.d9     6.5111225  -3.570949 16.593194 0.4444981
    ## m.n2-m.inc.d17     -0.9381208 -10.443589  8.567347 0.9999346
    ## m.inc.d8-m.inc.d17 -7.2559196 -16.761387  2.249548 0.2481345
    ## prolong-m.inc.d17  -0.6517483 -10.157216  8.853719 0.9999924
    ## extend-m.inc.d17    0.3417862  -9.163682  9.847254 0.9999998
    ## m.inc.d8-m.n2      -6.3177988 -15.823267  3.187669 0.4091962
    ## prolong-m.n2        0.2863725  -9.219095  9.791840 0.9999999
    ## extend-m.n2         1.2799070  -8.225561 10.785375 0.9996005
    ## prolong-m.inc.d8    6.6041713  -2.901297 16.109639 0.3553840
    ## extend-m.inc.d8     7.5977058  -1.907762 17.103174 0.2013881
    ## extend-prolong      0.9935345  -8.511933 10.499002 0.9999084
    ## 
    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment    6  620.3  103.38   14.13 5.23e-10 ***
    ## Residuals   61  446.3    7.32                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                          diff         lwr         upr     p adj
    ## m.inc.d9-m.inc.d3  -0.7454222  -4.6572913  3.16644693 0.9971438
    ## m.inc.d17-m.inc.d3 -5.5922585  -9.2804040 -1.90411292 0.0003874
    ## m.n2-m.inc.d3      -4.6194972  -8.3076427 -0.93135163 0.0055643
    ## m.inc.d8-m.inc.d3  -2.1456390  -5.8337846  1.54250651 0.5704396
    ## prolong-m.inc.d3   -6.8333494 -10.5214949 -3.14520383 0.0000091
    ## extend-m.inc.d3    -8.8443661 -12.5325116 -5.15622054 0.0000000
    ## m.inc.d17-m.inc.d9 -4.8468363  -8.7587054 -0.93496721 0.0063216
    ## m.n2-m.inc.d9      -3.8740750  -7.7859441  0.03779408 0.0538991
    ## m.inc.d8-m.inc.d9  -1.4002169  -5.3120860  2.51165222 0.9283253
    ## prolong-m.inc.d9   -6.0879272  -9.9997963 -2.17605812 0.0002523
    ## extend-m.inc.d9    -8.0989439 -12.0108130 -4.18707483 0.0000007
    ## m.n2-m.inc.d17      0.9727613  -2.7153843  4.66090684 0.9836890
    ## m.inc.d8-m.inc.d17  3.4466194  -0.2415261  7.13476498 0.0819476
    ## prolong-m.inc.d17  -1.2410909  -4.9292365  2.44705464 0.9459795
    ## extend-m.inc.d17   -3.2521076  -6.9402532  0.43603793 0.1187718
    ## m.inc.d8-m.n2       2.4738581  -1.2142874  6.16200369 0.3980001
    ## prolong-m.n2       -2.2138522  -5.9019978  1.47429335 0.5335139
    ## extend-m.n2        -4.2248689  -7.9130145 -0.53672336 0.0148419
    ## prolong-m.inc.d8   -4.6877103  -8.3758559 -0.99956479 0.0046655
    ## extend-m.inc.d8    -6.6987271 -10.3868726 -3.01058150 0.0000139
    ## extend-prolong     -2.0110167  -5.6991623  1.67712884 0.6430113

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
        scale_fill_viridis(na.value="#FFFFFF00") + 
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

    ## [1] "female_gonad"
    ## [1] TRUE

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 824 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    ##        V1 early extend prolong remove
    ## 1   early    NA    150    1008    227
    ## 2  extend   150     NA      63    762
    ## 3 prolong  1008     63      NA    557
    ## 4  remove   227    762     557     NA

![](../figures/manipulation/allthings-outcome-1.png)

    ## Warning in MASS::cov.trob(data[, vars]): Probable convergence failure

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## outcome      3    215   71.79   0.362  0.781
    ## Residuals   65  12890  198.30               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ outcome, data = pcadata)
    ## 
    ## $outcome
    ##                      diff        lwr      upr     p adj
    ## prolong-remove -2.6840623 -15.845492 10.47737 0.9494900
    ## early-remove   -1.8680995 -15.029530 11.29333 0.9819670
    ## extend-remove   3.3260144  -9.835416 16.48744 0.9093582
    ## early-prolong   0.8159628 -15.789552 17.42148 0.9992150
    ## extend-prolong  6.0100768 -10.595438 22.61559 0.7756094
    ## extend-early    5.1941140 -11.411401 21.79963 0.8425371
    ## 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## outcome      3    240   79.85   1.019   0.39
    ## Residuals   65   5093   78.35               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ outcome, data = pcadata)
    ## 
    ## $outcome
    ##                      diff        lwr       upr     p adj
    ## prolong-remove -5.4486918 -13.721663  2.824279 0.3134260
    ## early-remove   -0.5465411  -8.819512  7.726430 0.9981024
    ## extend-remove  -1.2791583  -9.552129  6.993812 0.9769231
    ## early-prolong   4.9021507  -5.535692 15.339994 0.6050126
    ## extend-prolong  4.1695335  -6.268309 14.607376 0.7188381
    ## extend-early   -0.7326172 -11.170460  9.705226 0.9977278
    ## 
    ## [1] "female_hypothalamus"
    ## [1] TRUE

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 31 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

![](../figures/manipulation/allthings-outcome-2.png)

    ##        V1 early extend prolong remove
    ## 1   early    NA      0       3    326
    ## 2  extend     0     NA       8    429
    ## 3 prolong     3      8      NA    207
    ## 4  remove   326    429     207     NA

![](../figures/manipulation/allthings-outcome-3.png)

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## outcome      3   79.1   26.36   1.298  0.283
    ## Residuals   66 1340.2   20.31               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ outcome, data = pcadata)
    ## 
    ## $outcome
    ##                      diff       lwr      upr     p adj
    ## prolong-remove  0.5497012 -3.649518 4.748921 0.9857567
    ## early-remove   -2.6403307 -6.839550 1.558889 0.3543892
    ## extend-remove  -1.6089119 -5.808131 2.590308 0.7442375
    ## early-prolong  -3.1900320 -8.501671 2.121607 0.3952817
    ## extend-prolong -2.1586131 -7.470252 3.153026 0.7081884
    ## extend-early    1.0314189 -4.280220 6.343058 0.9559977
    ## 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## outcome      3   13.2   4.413   0.251   0.86
    ## Residuals   66 1160.5  17.583               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ outcome, data = pcadata)
    ## 
    ## $outcome
    ##                      diff       lwr      upr     p adj
    ## prolong-remove  0.4933244 -3.414170 4.400819 0.9871815
    ## early-remove   -1.0082429 -4.915738 2.899252 0.9043016
    ## extend-remove   0.2694912 -3.638004 4.176986 0.9978462
    ## early-prolong  -1.5015672 -6.444201 3.441066 0.8538279
    ## extend-prolong -0.2238332 -5.166467 4.718800 0.9993857
    ## extend-early    1.2777341 -3.664899 6.220368 0.9038238
    ## 
    ## [1] "female_pituitary"
    ## [1] TRUE

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 121 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

![](../figures/manipulation/allthings-outcome-4.png)

    ##        V1 early extend prolong remove
    ## 1   early    NA   2624    2167   2071
    ## 2  extend  2624     NA     204   2564
    ## 3 prolong  2167    204      NA   1489
    ## 4  remove  2071   2564    1489     NA

![](../figures/manipulation/allthings-outcome-5.png)

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## outcome      3 1107.2   369.1   40.77 6.09e-15 ***
    ## Residuals   65  588.4     9.1                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ outcome, data = pcadata)
    ## 
    ## $outcome
    ##                      diff        lwr        upr     p adj
    ## prolong-remove  -7.390390 -10.202287 -4.5784936 0.0000000
    ## early-remove    -4.908275  -7.720172 -2.0963783 0.0001151
    ## extend-remove  -10.295604 -13.107501 -7.4837074 0.0000000
    ## early-prolong    2.482115  -1.065599  6.0298296 0.2621399
    ## extend-prolong  -2.905214  -6.452928  0.6425005 0.1456989
    ## extend-early    -5.387329  -8.935043 -1.8396148 0.0009152
    ## 
    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## outcome      3  203.1   67.71   4.327 0.00765 **
    ## Residuals   65 1017.2   15.65                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ outcome, data = pcadata)
    ## 
    ## $outcome
    ##                      diff        lwr        upr     p adj
    ## prolong-remove  1.1410270 -2.5562941  4.8383482 0.8477669
    ## early-remove   -4.1144041 -7.8117253 -0.4170830 0.0233190
    ## extend-remove   1.5424843 -2.1548368  5.2398054 0.6907838
    ## early-prolong  -5.2554312 -9.9202677 -0.5905946 0.0211434
    ## extend-prolong  0.4014572 -4.2633793  5.0662938 0.9958396
    ## extend-early    5.6568884  0.9920518 10.3217250 0.0112364
    ## 
    ## [1] "male_gonad"
    ## [1] TRUE

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 792 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

![](../figures/manipulation/allthings-outcome-6.png)

    ##        V1 early extend prolong remove
    ## 1   early    NA    202     379     20
    ## 2  extend   202     NA     420      8
    ## 3 prolong   379    420      NA   1092
    ## 4  remove    20      8    1092     NA

![](../figures/manipulation/allthings-outcome-7.png)

    ## Warning in MASS::cov.trob(data[, vars]): Probable convergence failure

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## outcome      3    512  170.74   2.077  0.112
    ## Residuals   63   5178   82.19               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ outcome, data = pcadata)
    ## 
    ## $outcome
    ##                      diff         lwr       upr     p adj
    ## prolong-remove  8.0234518  -0.5035053 16.550409 0.0723999
    ## early-remove    0.9132040  -7.6137531  9.440161 0.9920482
    ## extend-remove   1.7327077  -6.7942495 10.259665 0.9498650
    ## early-prolong  -7.1102478 -17.8096869  3.589191 0.3052289
    ## extend-prolong -6.2907442 -16.9901833  4.408695 0.4134790
    ## extend-early    0.8195037  -9.8799355 11.518943 0.9970466
    ## 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## outcome      3  128.1   42.69   1.284  0.288
    ## Residuals   63 2094.5   33.25               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ outcome, data = pcadata)
    ## 
    ## $outcome
    ##                      diff        lwr       upr     p adj
    ## prolong-remove  0.1004901  -5.322591  5.523571 0.9999575
    ## early-remove   -3.1921135  -8.615194  2.230967 0.4124692
    ## extend-remove   1.7015944  -3.721487  7.124675 0.8409715
    ## early-prolong  -3.2926036 -10.097367  3.512160 0.5808804
    ## extend-prolong  1.6011043  -5.203659  8.405868 0.9250108
    ## extend-early    4.8937079  -1.911055 11.698471 0.2393912
    ## 
    ## [1] "male_hypothalamus"
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

![](../figures/manipulation/allthings-outcome-8.png)

    ##        V1 early extend prolong remove
    ## 1   early    NA      0      12   2246
    ## 2  extend     0     NA       0      6
    ## 3 prolong    12      0      NA      0
    ## 4  remove  2246      6       0     NA

![](../figures/manipulation/allthings-outcome-9.png)

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## outcome      3  173.1   57.69   2.062  0.114
    ## Residuals   64 1790.3   27.97               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ outcome, data = pcadata)
    ## 
    ## $outcome
    ##                      diff        lwr       upr     p adj
    ## prolong-remove -0.5072784  -5.465730 4.4511735 0.9930591
    ## early-remove   -4.5307231  -9.489175 0.4277288 0.0853089
    ## extend-remove  -1.8953664  -6.853818 3.0630855 0.7451485
    ## early-prolong  -4.0234447 -10.262693 2.2158038 0.3316212
    ## extend-prolong -1.3880880  -7.627337 4.8511604 0.9356999
    ## extend-early    2.6353566  -3.603892 8.8746051 0.6821910
    ## 
    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## outcome      3   66.4   22.13   1.803  0.156
    ## Residuals   64  785.8   12.28               
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ outcome, data = pcadata)
    ## 
    ## $outcome
    ##                      diff       lwr      upr     p adj
    ## prolong-remove -2.0802474 -5.365358 1.204863 0.3476464
    ## early-remove    1.0177474 -2.267363 4.302858 0.8461178
    ## extend-remove  -1.4965806 -4.781691 1.788530 0.6280266
    ## early-prolong   3.0979948 -1.035679 7.231668 0.2074453
    ## extend-prolong  0.5836668 -3.550007 4.717340 0.9822138
    ## extend-early   -2.5143279 -6.648002 1.619346 0.3834251
    ## 
    ## [1] "male_pituitary"
    ## [1] TRUE

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 477 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

![](../figures/manipulation/allthings-outcome-10.png)

    ##        V1 early extend prolong remove
    ## 1   early    NA   1783     848   1048
    ## 2  extend  1783     NA       5    683
    ## 3 prolong   848      5      NA    146
    ## 4  remove  1048    683     146     NA

![](../figures/manipulation/allthings-outcome-11.png)![](../figures/manipulation/allthings-outcome-12.png)

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## outcome      3  338.7  112.92   2.373 0.0785 .
    ## Residuals   64 3045.3   47.58                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ outcome, data = pcadata)
    ## 
    ## $outcome
    ##                      diff        lwr       upr     p adj
    ## prolong-remove  1.1357123  -5.331271  7.602695 0.9667659
    ## early-remove   -5.3561456 -11.823129  1.110838 0.1385003
    ## extend-remove   2.1304584  -4.336525  8.597442 0.8207883
    ## early-prolong  -6.4918578 -14.629300  1.645584 0.1627185
    ## extend-prolong  0.9947461  -7.142696  9.132188 0.9883017
    ## extend-early    7.4866040  -0.650838 15.624046 0.0821767
    ## 
    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## outcome      3  389.3  129.78   12.88 1.11e-06 ***
    ## Residuals   64  645.0   10.08                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ outcome, data = pcadata)
    ## 
    ## $outcome
    ##                      diff         lwr        upr     p adj
    ## prolong-remove -3.9569181  -6.9330940 -0.9807421 0.0045231
    ## early-remove    0.6389908  -2.3371852  3.6151667 0.9416867
    ## extend-remove  -5.9604883  -8.9366643 -2.9843124 0.0000096
    ## early-prolong   4.5959088   0.8509696  8.3408481 0.0100807
    ## extend-prolong -2.0035703  -5.7485095  1.7413690 0.4970789
    ## extend-early   -6.5994791 -10.3444184 -2.8545398 0.0000996

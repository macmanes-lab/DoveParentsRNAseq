    library(tidyverse)
    library(DESeq2)
    library(cowplot)
    library(RColorBrewer)
    library(pheatmap)
    library(kableExtra)
    library(viridis)

    # load custom functions  
    source("../R/functions.R")  

    knitr::opts_chunk$set(fig.path = '../figures/pit/',cache=TRUE)

Females only
============

    # import "colData" which contains sample information and "countData" which contains read counts
    colData <- read.csv("../results/00_colData_characterization.csv", header = T, row.names = 1)
    countData <- read.csv("../results/00_countData_characterization.csv", header = T, row.names = 1)  
    geneinfo <- read.csv("../results/00_geneinfo.csv", row.names = 1)

    colData$treatment <- factor(colData$treatment, levels = 
                                  c("control", "bldg", "lay", "inc.d3", "inc.d9", 
                                    "inc.d17", "hatch", "n5", "n9"))  

    colData <- colData %>%
      dplyr::filter(grepl('pituitary', tissue)) %>%
      dplyr::filter(sex == "female") %>%
      droplevels()
    row.names(colData) <- colData$V1

    # print sample sizes
    colData %>% select(group, tissue)  %>%  summary()

    ##                       group          tissue  
    ##  female.pituitary.inc.d9 :13   pituitary:96  
    ##  female.pituitary.control:11                 
    ##  female.pituitary.inc.d17:11                 
    ##  female.pituitary.n9     :11                 
    ##  female.pituitary.bldg   :10                 
    ##  female.pituitary.hatch  :10                 
    ##  (Other)                 :30

    savecols <- as.character(colData$V1) 
    savecols <- as.vector(savecols) 
    countData <- countData %>% dplyr::select(one_of(savecols)) 

    # check that row and col lenghts are equal
    ncol(countData) == nrow(colData)     

    ## [1] TRUE

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ treatment ) 
    dds <- dds[ rowSums(counts(dds)) > 2, ] ## pre-filter genes 
    dds <- DESeq(dds) # Differential expression analysis

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 80 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    vsd <- vst(dds, blind=FALSE) # variance stabilized   

    #create list of groups
    a <- levels(colData$treatment)
    b <- levels(colData$treatment)

    # slim for testing
    #a <- c("n9", "bldg" , "lay" )
    #b <- c("n9", "bldg" , "lay" )

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

    ##                     V1      V2   V3
    ## controlbldg    control    bldg 6448
    ## controllay     control     lay 6315
    ## controlinc.d3  control  inc.d3 5550
    ## controlinc.d9  control  inc.d9 6244
    ## controlinc.d17 control inc.d17 5341
    ## controlhatch   control   hatch 5966

    # widen data to create table of degs
    rownames(dat) <- NULL #remove row names
    data_wide <- spread(dat, V2, V3)
    data_wide

    ##        V1 bldg control hatch inc.d17 inc.d3 inc.d9  lay   n5   n9
    ## 1    bldg   NA    6448  2682    3427     37    103  280  494  294
    ## 2 control 6448      NA  5966    5341   5550   6244 6315 6140 6141
    ## 3   hatch 2682    5966    NA       6   1083   1671 3348 1089 2319
    ## 4 inc.d17 3427    5341     6      NA   1132   2017 3995 1711 2883
    ## 5  inc.d3   37    5550  1083    1132     NA      1  333  198  377
    ## 6  inc.d9  103    6244  1671    2017      1     NA  107   87  277
    ## 7     lay  280    6315  3348    3995    333    107   NA  460  570
    ## 8      n5  494    6140  1089    1711    198     87  460   NA   23
    ## 9      n9  294    6141  2319    2883    377    277  570   23   NA

    kable(data_wide) 

<table>
<thead>
<tr>
<th style="text-align:left;">
V1
</th>
<th style="text-align:right;">
bldg
</th>
<th style="text-align:right;">
control
</th>
<th style="text-align:right;">
hatch
</th>
<th style="text-align:right;">
inc.d17
</th>
<th style="text-align:right;">
inc.d3
</th>
<th style="text-align:right;">
inc.d9
</th>
<th style="text-align:right;">
lay
</th>
<th style="text-align:right;">
n5
</th>
<th style="text-align:right;">
n9
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
bldg
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
6448
</td>
<td style="text-align:right;">
2682
</td>
<td style="text-align:right;">
3427
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
103
</td>
<td style="text-align:right;">
280
</td>
<td style="text-align:right;">
494
</td>
<td style="text-align:right;">
294
</td>
</tr>
<tr>
<td style="text-align:left;">
control
</td>
<td style="text-align:right;">
6448
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
5966
</td>
<td style="text-align:right;">
5341
</td>
<td style="text-align:right;">
5550
</td>
<td style="text-align:right;">
6244
</td>
<td style="text-align:right;">
6315
</td>
<td style="text-align:right;">
6140
</td>
<td style="text-align:right;">
6141
</td>
</tr>
<tr>
<td style="text-align:left;">
hatch
</td>
<td style="text-align:right;">
2682
</td>
<td style="text-align:right;">
5966
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
1083
</td>
<td style="text-align:right;">
1671
</td>
<td style="text-align:right;">
3348
</td>
<td style="text-align:right;">
1089
</td>
<td style="text-align:right;">
2319
</td>
</tr>
<tr>
<td style="text-align:left;">
inc.d17
</td>
<td style="text-align:right;">
3427
</td>
<td style="text-align:right;">
5341
</td>
<td style="text-align:right;">
6
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1132
</td>
<td style="text-align:right;">
2017
</td>
<td style="text-align:right;">
3995
</td>
<td style="text-align:right;">
1711
</td>
<td style="text-align:right;">
2883
</td>
</tr>
<tr>
<td style="text-align:left;">
inc.d3
</td>
<td style="text-align:right;">
37
</td>
<td style="text-align:right;">
5550
</td>
<td style="text-align:right;">
1083
</td>
<td style="text-align:right;">
1132
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
333
</td>
<td style="text-align:right;">
198
</td>
<td style="text-align:right;">
377
</td>
</tr>
<tr>
<td style="text-align:left;">
inc.d9
</td>
<td style="text-align:right;">
103
</td>
<td style="text-align:right;">
6244
</td>
<td style="text-align:right;">
1671
</td>
<td style="text-align:right;">
2017
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
107
</td>
<td style="text-align:right;">
87
</td>
<td style="text-align:right;">
277
</td>
</tr>
<tr>
<td style="text-align:left;">
lay
</td>
<td style="text-align:right;">
280
</td>
<td style="text-align:right;">
6315
</td>
<td style="text-align:right;">
3348
</td>
<td style="text-align:right;">
3995
</td>
<td style="text-align:right;">
333
</td>
<td style="text-align:right;">
107
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
460
</td>
<td style="text-align:right;">
570
</td>
</tr>
<tr>
<td style="text-align:left;">
n5
</td>
<td style="text-align:right;">
494
</td>
<td style="text-align:right;">
6140
</td>
<td style="text-align:right;">
1089
</td>
<td style="text-align:right;">
1711
</td>
<td style="text-align:right;">
198
</td>
<td style="text-align:right;">
87
</td>
<td style="text-align:right;">
460
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
23
</td>
</tr>
<tr>
<td style="text-align:left;">
n9
</td>
<td style="text-align:right;">
294
</td>
<td style="text-align:right;">
6141
</td>
<td style="text-align:right;">
2319
</td>
<td style="text-align:right;">
2883
</td>
<td style="text-align:right;">
377
</td>
<td style="text-align:right;">
277
</td>
<td style="text-align:right;">
570
</td>
<td style="text-align:right;">
23
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>

    dat$V1 <- factor(dat$V1, levels = 
                                  c("control", "bldg", "lay", "inc.d3", "inc.d9", 
                                    "inc.d17", "hatch", "n5", "n9"))
    dat$V2 <- factor(dat$V2, levels = 
                                  c("control", "bldg", "lay", "inc.d3", "inc.d9", 
                                    "inc.d17", "hatch", "n5", "n9"))

    dat %>%
      filter(V1 != "control", V2 != "control")  %>%
      ggplot( aes(V1, V2)) +
        geom_tile(aes(fill = V3)) +
        scale_fill_viridis(na.value="#FFFFFF00", 
                         limits = c(0, 5025),
                         breaks = c(0, 1000, 2000, 3000, 4000, 5000)) + 
        xlab(" ") + ylab("Timepoint") +
        labs(fill = "# of DEGs") 

![](../figures/pit/restable-1.png)

    # create the dataframe using my function pcadataframe
    pcadata <- pcadataframe(vsd, intgroup=c("treatment"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))
    percentVar

    ## [1] 12 10  8  6  4  3

    ggplot(pcadata, aes(PC1, PC2,color = treatment)) + 
      geom_point(size = 2, alpha = 1) +
      stat_ellipse(type = "t") +
      xlab(paste0("PC1: ", percentVar[1],"% variance")) +
      ylab(paste0("PC2: ", percentVar[2],"% variance")) +
      theme_cowplot(font_size = 8, line_size = 0.25) 

![](../figures/pit/PCA-1.png)

    summary(aov(PC1 ~ treatment, data=pcadata)) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment    8   1596  199.55   13.76 1.04e-12 ***
    ## Residuals   87   1262   14.51                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov(PC1 ~ treatment, data=pcadata), which = "treatment") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                       diff          lwr        upr     p adj
    ## bldg-control     0.8726313  -4.42145379  6.1667163 0.9998419
    ## lay-control      1.5629395  -3.73114553  6.8570246 0.9899702
    ## inc.d3-control   3.7211880  -1.57289707  9.0152730 0.3923453
    ## inc.d9-control   4.5567549  -0.40705584  9.5205656 0.0979859
    ## inc.d17-control 11.5434875   6.37698955 16.7099855 0.0000000
    ## hatch-control   11.7328614   6.43877632 17.0269464 0.0000000
    ## n5-control       7.0487637   1.75467865 12.3428488 0.0017788
    ## n9-control       2.8802424  -2.28625556  8.0467404 0.6991175
    ## lay-bldg         0.6903083  -4.72836055  6.1089771 0.9999778
    ## inc.d3-bldg      2.8485567  -2.57011208  8.2672255 0.7613604
    ## inc.d9-bldg      3.6841236  -1.41235085  8.7805981 0.3538560
    ## inc.d17-bldg    10.6708563   5.37677122 15.9649413 0.0000003
    ## hatch-bldg      10.8602301   5.44156130 16.2788989 0.0000003
    ## n5-bldg          6.1761324   0.75746363 11.5948012 0.0136505
    ## n9-bldg          2.0076112  -3.28647389  7.3016962 0.9530338
    ## inc.d3-lay       2.1582485  -3.26042034  7.5769173 0.9380866
    ## inc.d9-lay       2.9938153  -2.10265911  8.0902898 0.6366080
    ## inc.d17-lay      9.9805480   4.68646296 15.2746331 0.0000016
    ## hatch-lay       10.1699218   4.75125304 15.5885907 0.0000018
    ## n5-lay           5.4858242   0.06715537 10.9044930 0.0448940
    ## n9-lay           1.3173029  -3.97678215  6.6113880 0.9968417
    ## inc.d9-inc.d3    0.8355669  -4.26090757  5.9320413 0.9998481
    ## inc.d17-inc.d3   7.8222996   2.52821450 13.1163846 0.0003200
    ## hatch-inc.d3     8.0116734   2.59300458 13.4303422 0.0003163
    ## n5-inc.d3        3.3275757  -2.09109309  8.7462445 0.5792243
    ## n9-inc.d3       -0.8409456  -6.13503061  4.4531395 0.9998803
    ## inc.d17-inc.d9   6.9867327   2.02292197 11.9505434 0.0007394
    ## hatch-inc.d9     7.1761065   2.07963206 12.2725810 0.0007349
    ## n5-inc.d9        2.4920088  -2.60446561  7.5884833 0.8255836
    ## n9-inc.d9       -1.6765124  -6.64032314  3.2872983 0.9764062
    ## hatch-inc.d17    0.1893738  -5.10471122  5.4834589 1.0000000
    ## n5-inc.d17      -4.4947238  -9.78880889  0.7993612 0.1628792
    ## n9-inc.d17      -8.6632451 -13.82974310 -3.4967471 0.0000259
    ## n5-hatch        -4.6840977 -10.10276648  0.7345711 0.1460890
    ## n9-hatch        -8.8526189 -14.14670400 -3.5585339 0.0000275
    ## n9-n5           -4.1685213  -9.46260633  1.1255638 0.2442148

    summary(aov(PC2 ~ treatment, data=pcadata)) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment    8 1530.7  191.33   17.62 2.44e-15 ***
    ## Residuals   87  944.5   10.86                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov(PC2 ~ treatment, data=pcadata), which = "treatment") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                        diff         lwr       upr     p adj
    ## bldg-control    -13.6660334 -18.2458104 -9.086256 0.0000000
    ## lay-control     -12.2955919 -16.8753689 -7.715815 0.0000000
    ## inc.d3-control  -11.4844565 -16.0642336 -6.904679 0.0000000
    ## inc.d9-control  -12.7368811 -17.0309463 -8.442816 0.0000000
    ## inc.d17-control  -7.8891232 -12.3585279 -3.419718 0.0000080
    ## hatch-control    -7.4285976 -12.0083746 -2.848821 0.0000526
    ## n5-control      -11.0479415 -15.6277186 -6.468165 0.0000000
    ## n9-control      -10.8056969 -15.2751016 -6.336292 0.0000000
    ## lay-bldg          1.3704415  -3.3171097  6.057993 0.9905922
    ## inc.d3-bldg       2.1815769  -2.5059743  6.869128 0.8614767
    ## inc.d9-bldg       0.9291523  -3.4796769  5.337981 0.9990282
    ## inc.d17-bldg      5.7769102   1.1971332 10.356687 0.0038627
    ## hatch-bldg        6.2374359   1.5498846 10.924987 0.0017950
    ## n5-bldg           2.6180919  -2.0694594  7.305643 0.6970080
    ## n9-bldg           2.8603366  -1.7194405  7.440114 0.5565557
    ## inc.d3-lay        0.8111354  -3.8764159  5.498687 0.9997723
    ## inc.d9-lay       -0.4412892  -4.8501184  3.967540 0.9999966
    ## inc.d17-lay       4.4064687  -0.1733083  8.986246 0.0688078
    ## hatch-lay         4.8669943   0.1794431  9.554546 0.0356750
    ## n5-lay            1.2476504  -3.4399009  5.935202 0.9949811
    ## n9-lay            1.4898950  -3.0898820  6.069672 0.9813060
    ## inc.d9-inc.d3    -1.2524246  -5.6612537  3.156405 0.9922213
    ## inc.d17-inc.d3    3.5953334  -0.9844437  8.175110 0.2477733
    ## hatch-inc.d3      4.0558590  -0.6316923  8.743410 0.1452522
    ## n5-inc.d3         0.4365150  -4.2510363  5.124066 0.9999981
    ## n9-inc.d3         0.6787597  -3.9010174  5.258537 0.9999293
    ## inc.d17-inc.d9    4.8477579   0.5536928  9.141823 0.0151908
    ## hatch-inc.d9      5.3082835   0.8994544  9.717113 0.0071072
    ## n5-inc.d9         1.6889395  -2.7198896  6.097769 0.9502364
    ## n9-inc.d9         1.9311842  -2.3628809  6.225249 0.8827968
    ## hatch-inc.d17     0.4605256  -4.1192514  5.040303 0.9999965
    ## n5-inc.d17       -3.1588184  -7.7385954  1.420959 0.4186688
    ## n9-inc.d17       -2.9165737  -7.3859784  1.552831 0.4962133
    ## n5-hatch         -3.6193440  -8.3068952  1.068207 0.2679772
    ## n9-hatch         -3.3770993  -7.9568763  1.202678 0.3271377
    ## n9-n5             0.2422447  -4.3375323  4.822022 1.0000000

    ggplot(pcadata, aes(PC3, PC4,color = treatment)) + 
      geom_point(size = 2, alpha = 1) +
      stat_ellipse(type = "t") +
      xlab(paste0("PC3: ", percentVar[3],"% variance")) +
      ylab(paste0("PC4: ", percentVar[4],"% variance")) +
      theme_cowplot(font_size = 8, line_size = 0.25) 

![](../figures/pit/PCA-2.png)

    summary(aov(PC3 ~ treatment, data=pcadata)) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment    8  597.7   74.71   4.702 8.25e-05 ***
    ## Residuals   87 1382.5   15.89                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    summary(aov(PC4 ~ treatment, data=pcadata)) 

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment    8  515.8   64.47   5.312 1.95e-05 ***
    ## Residuals   87 1056.0   12.14                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    summary(aov(PC5 ~ treatment, data=pcadata)) 

    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## treatment    8  253.8   31.72   3.399 0.00192 **
    ## Residuals   87  812.0    9.33                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    summary(aov(PC6 ~ treatment, data=pcadata))   

    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## treatment    8  178.5  22.311   2.856 0.00724 **
    ## Residuals   87  679.7   7.813                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Heatmaps

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
             fontsize = 6)  

![](../figures/pit/heatmap-1.png)

Males only
==========

    # import "colData" which contains sample information and "countData" which contains read counts
    colData <- read.csv("../results/00_colData_characterization.csv", header = T, row.names = 1)
    countData <- read.csv("../results/00_countData_characterization.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../results/00_geneinfo.csv", row.names = 1)  

    colData$treatment <- factor(colData$treatment, levels = 
                                  c("control", "bldg", "lay", "inc.d3", "inc.d9", 
                                    "inc.d17", "hatch", "n5", "n9"))

    colData <- colData %>%
      dplyr::filter(grepl('pituitary', tissue)) %>%
      dplyr::filter(sex == "male") %>%
      droplevels()
    row.names(colData) <- colData$V1

    # print sample sizes
    colData %>% select(group, tissue)  %>%  summary()

    ##                     group          tissue  
    ##  male.pituitary.control:14   pituitary:97  
    ##  male.pituitary.inc.d17:11                 
    ##  male.pituitary.inc.d9 :11                 
    ##  male.pituitary.n9     :11                 
    ##  male.pituitary.bldg   :10                 
    ##  male.pituitary.hatch  :10                 
    ##  (Other)               :30

    savecols <- as.character(colData$V1) 
    savecols <- as.vector(savecols) 
    countData <- countData %>% dplyr::select(one_of(savecols)) 

    # check that row and col lenghts are equal
    ncol(countData) == nrow(colData) 

    ## [1] TRUE

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ treatment )
    dds <- dds[ rowSums(counts(dds)) > 2, ] ## pre-filter genes 
    dds <- DESeq(dds) # Differential expression analysis

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

    vsd <- vst(dds, blind=FALSE) # variance stabilized   

    #create list of groups
    a <- levels(colData$treatment)
    b <- levels(colData$treatment)

    # slim for testing
    #a <- c("n9", "bldg" , "lay" )
    #b <- c("n9", "bldg" , "lay" )

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

    ##                     V1      V2   V3
    ## controlbldg    control    bldg 7026
    ## controllay     control     lay 7304
    ## controlinc.d3  control  inc.d3 7028
    ## controlinc.d9  control  inc.d9 6997
    ## controlinc.d17 control inc.d17 6724
    ## controlhatch   control   hatch 7026

    # widen data to create table of degs
    rownames(dat) <- NULL #remove row names
    data_wide <- spread(dat, V2, V3)
    data_wide

    ##        V1 bldg control hatch inc.d17 inc.d3 inc.d9  lay   n5   n9
    ## 1    bldg   NA    7026  1953     880     14      0    0  612   99
    ## 2 control 7026      NA  7026    6724   7028   6997 7304 7304 7114
    ## 3   hatch 1953    7026    NA      56   1880   1157 2231  155 1207
    ## 4 inc.d17  880    6724    56      NA    660    433  951  142  486
    ## 5  inc.d3   14    7028  1880     660     NA      1    5  703  183
    ## 6  inc.d9    0    6997  1157     433      1     NA    5  541   35
    ## 7     lay    0    7304  2231     951      5      5   NA  677   75
    ## 8      n5  612    7304   155     142    703    541  677   NA   63
    ## 9      n9   99    7114  1207     486    183     35   75   63   NA

    kable(data_wide) 

<table>
<thead>
<tr>
<th style="text-align:left;">
V1
</th>
<th style="text-align:right;">
bldg
</th>
<th style="text-align:right;">
control
</th>
<th style="text-align:right;">
hatch
</th>
<th style="text-align:right;">
inc.d17
</th>
<th style="text-align:right;">
inc.d3
</th>
<th style="text-align:right;">
inc.d9
</th>
<th style="text-align:right;">
lay
</th>
<th style="text-align:right;">
n5
</th>
<th style="text-align:right;">
n9
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
bldg
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
7026
</td>
<td style="text-align:right;">
1953
</td>
<td style="text-align:right;">
880
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
612
</td>
<td style="text-align:right;">
99
</td>
</tr>
<tr>
<td style="text-align:left;">
control
</td>
<td style="text-align:right;">
7026
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
7026
</td>
<td style="text-align:right;">
6724
</td>
<td style="text-align:right;">
7028
</td>
<td style="text-align:right;">
6997
</td>
<td style="text-align:right;">
7304
</td>
<td style="text-align:right;">
7304
</td>
<td style="text-align:right;">
7114
</td>
</tr>
<tr>
<td style="text-align:left;">
hatch
</td>
<td style="text-align:right;">
1953
</td>
<td style="text-align:right;">
7026
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
56
</td>
<td style="text-align:right;">
1880
</td>
<td style="text-align:right;">
1157
</td>
<td style="text-align:right;">
2231
</td>
<td style="text-align:right;">
155
</td>
<td style="text-align:right;">
1207
</td>
</tr>
<tr>
<td style="text-align:left;">
inc.d17
</td>
<td style="text-align:right;">
880
</td>
<td style="text-align:right;">
6724
</td>
<td style="text-align:right;">
56
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
660
</td>
<td style="text-align:right;">
433
</td>
<td style="text-align:right;">
951
</td>
<td style="text-align:right;">
142
</td>
<td style="text-align:right;">
486
</td>
</tr>
<tr>
<td style="text-align:left;">
inc.d3
</td>
<td style="text-align:right;">
14
</td>
<td style="text-align:right;">
7028
</td>
<td style="text-align:right;">
1880
</td>
<td style="text-align:right;">
660
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
703
</td>
<td style="text-align:right;">
183
</td>
</tr>
<tr>
<td style="text-align:left;">
inc.d9
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
6997
</td>
<td style="text-align:right;">
1157
</td>
<td style="text-align:right;">
433
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
35
</td>
</tr>
<tr>
<td style="text-align:left;">
lay
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
7304
</td>
<td style="text-align:right;">
2231
</td>
<td style="text-align:right;">
951
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
5
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
677
</td>
<td style="text-align:right;">
75
</td>
</tr>
<tr>
<td style="text-align:left;">
n5
</td>
<td style="text-align:right;">
612
</td>
<td style="text-align:right;">
7304
</td>
<td style="text-align:right;">
155
</td>
<td style="text-align:right;">
142
</td>
<td style="text-align:right;">
703
</td>
<td style="text-align:right;">
541
</td>
<td style="text-align:right;">
677
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
63
</td>
</tr>
<tr>
<td style="text-align:left;">
n9
</td>
<td style="text-align:right;">
99
</td>
<td style="text-align:right;">
7114
</td>
<td style="text-align:right;">
1207
</td>
<td style="text-align:right;">
486
</td>
<td style="text-align:right;">
183
</td>
<td style="text-align:right;">
35
</td>
<td style="text-align:right;">
75
</td>
<td style="text-align:right;">
63
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>

    dat$V1 <- factor(dat$V1, levels = 
                                  c("control", "bldg", "lay", "inc.d3", "inc.d9", 
                                    "inc.d17", "hatch", "n5", "n9"))
    dat$V2 <- factor(dat$V2, levels = 
                                  c("control", "bldg", "lay", "inc.d3", "inc.d9", 
                                    "inc.d17", "hatch", "n5", "n9"))


    dat %>%
      filter(V1 != "control", V2 != "control")  %>%
      ggplot( aes(V1, V2)) +
        geom_tile(aes(fill = V3)) +
        scale_fill_viridis(na.value="#FFFFFF00", 
                         limits = c(0, 5025),
                         breaks = c(0, 1000, 2000, 3000, 4000, 5000)) + 
        xlab(" ") + ylab("Timepoint") +
        labs(fill = "# of DEGs")

![](../figures/pit/restableMale-1.png)

    # create the dataframe using my function pcadataframe
    pcadata <- pcadataframe(vsd, intgroup=c("treatment"), returnData=TRUE)
    percentVar <- round(100 * attr(pcadata, "percentVar"))
    percentVar

    ## [1] 19  9  6  5  4  3

    ggplot(pcadata, aes(PC1, PC2,color = treatment)) + 
      geom_point(size = 2, alpha = 1) +
      stat_ellipse(type = "t") +
      xlab(paste0("PC1: ", percentVar[1],"% variance")) +
      ylab(paste0("PC2: ", percentVar[2],"% variance")) +
      theme_cowplot(font_size = 8, line_size = 0.25) 

![](../figures/pit/PCAmale-1.png)

    summary(aov(PC1 ~ treatment, data=pcadata)) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)    
    ## treatment    8   3357   419.6   30.11 <2e-16 ***
    ## Residuals   88   1227    13.9                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov(PC1 ~ treatment, data=pcadata), which = "treatment") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC1 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                        diff       lwr       upr     p adj
    ## bldg-control    16.46575142 11.549636 21.381866 0.0000000
    ## lay-control     15.75436897 10.838254 20.670484 0.0000000
    ## inc.d3-control  16.06123762 11.145123 20.977353 0.0000000
    ## inc.d9-control  16.62790213 11.843914 21.411890 0.0000000
    ## inc.d17-control 16.00898345 11.224995 20.792971 0.0000000
    ## hatch-control   17.56612897 12.650014 22.482244 0.0000000
    ## n5-control      18.16332497 13.247210 23.079440 0.0000000
    ## n9-control      16.38950010 11.605512 21.173488 0.0000000
    ## lay-bldg        -0.71138245 -6.021394  4.598629 0.9999674
    ## inc.d3-bldg     -0.40451380 -5.714525  4.905497 0.9999996
    ## inc.d9-bldg      0.16215071 -5.025775  5.350076 1.0000000
    ## inc.d17-bldg    -0.45676797 -5.644694  4.731158 0.9999988
    ## hatch-bldg       1.10037755 -4.209634  6.410389 0.9991426
    ## n5-bldg          1.69757355 -3.612438  7.007585 0.9832944
    ## n9-bldg         -0.07625132 -5.264177  5.111674 1.0000000
    ## inc.d3-lay       0.30686865 -5.003143  5.616880 1.0000000
    ## inc.d9-lay       0.87353316 -4.314392  6.061459 0.9998151
    ## inc.d17-lay      0.25461448 -4.933311  5.442540 1.0000000
    ## hatch-lay        1.81176000 -3.498251  7.121771 0.9749506
    ## n5-lay           2.40895600 -2.901055  7.718967 0.8778542
    ## n9-lay           0.63513114 -4.552794  5.823057 0.9999838
    ## inc.d9-inc.d3    0.56666451 -4.621261  5.754590 0.9999933
    ## inc.d17-inc.d3  -0.05225417 -5.240180  5.135671 1.0000000
    ## hatch-inc.d3     1.50489135 -3.805120  6.814903 0.9923632
    ## n5-inc.d3        2.10208735 -3.207924  7.412099 0.9402724
    ## n9-inc.d3        0.32826249 -4.859663  5.516188 0.9999999
    ## inc.d17-inc.d9  -0.61891868 -5.681816  4.443978 0.9999839
    ## hatch-inc.d9     0.93822684 -4.249699  6.126152 0.9996847
    ## n5-inc.d9        1.53542284 -3.652503  6.723348 0.9898329
    ## n9-inc.d9       -0.23840202 -5.301299  4.824495 1.0000000
    ## hatch-inc.d17    1.55714552 -3.630780  6.745071 0.9888566
    ## n5-inc.d17       2.15434151 -3.033584  7.342267 0.9224832
    ## n9-inc.d17       0.38051665 -4.682380  5.443414 0.9999996
    ## n5-hatch         0.59719600 -4.712815  5.907207 0.9999916
    ## n9-hatch        -1.17662886 -6.364554  4.011297 0.9983595
    ## n9-n5           -1.77382486 -6.961750  3.414101 0.9746271

    summary(aov(PC2 ~ treatment, data=pcadata)) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)    
    ## treatment    8 1542.7  192.84   27.67 <2e-16 ***
    ## Residuals   88  613.3    6.97                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov(PC2 ~ treatment, data=pcadata), which = "treatment") 

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = PC2 ~ treatment, data = pcadata)
    ## 
    ## $treatment
    ##                       diff         lwr        upr     p adj
    ## bldg-control    -4.7712123  -8.2474830 -1.2949417 0.0011071
    ## lay-control     -3.9981196  -7.4743902 -0.5218490 0.0123013
    ## inc.d3-control  -4.5373430  -8.0136136 -1.0610724 0.0023783
    ## inc.d9-control  -2.9718195  -6.3546608  0.4110217 0.1320462
    ## inc.d17-control  4.2592271   0.8763858  7.6420683 0.0039503
    ## hatch-control    6.9493928   3.4731222 10.4256634 0.0000003
    ## n5-control       4.1288817   0.6526111  7.6051523 0.0083966
    ## n9-control      -0.1294669  -3.5123082  3.2533743 1.0000000
    ## lay-bldg         0.7730928  -2.9817087  4.5278942 0.9991819
    ## inc.d3-bldg      0.2338694  -3.5209320  3.9886708 0.9999999
    ## inc.d9-bldg      1.7993928  -1.8690798  5.4678654 0.8233582
    ## inc.d17-bldg     9.0304394   5.3619668 12.6989120 0.0000000
    ## hatch-bldg      11.7206051   7.9658037 15.4754065 0.0000000
    ## n5-bldg          8.9000940   5.1452926 12.6548954 0.0000000
    ## n9-bldg          4.6417454   0.9732728  8.3102180 0.0036917
    ## inc.d3-lay      -0.5392234  -4.2940248  3.2155780 0.9999446
    ## inc.d9-lay       1.0263001  -2.6421725  4.6947727 0.9929952
    ## inc.d17-lay      8.2573467   4.5888741 11.9258193 0.0000000
    ## hatch-lay       10.9475124   7.1927110 14.7023138 0.0000000
    ## n5-lay           8.1270013   4.3721999 11.8818027 0.0000000
    ## n9-lay           3.8686527   0.2001801  7.5371253 0.0307585
    ## inc.d9-inc.d3    1.5655235  -2.1029491  5.2339961 0.9104541
    ## inc.d17-inc.d3   8.7965701   5.1280975 12.4650426 0.0000000
    ## hatch-inc.d3    11.4867358   7.7319343 15.2415372 0.0000000
    ## n5-inc.d3        8.6662246   4.9114232 12.4210260 0.0000000
    ## n9-inc.d3        4.4078761   0.7394035  8.0763487 0.0072768
    ## inc.d17-inc.d9   7.2310466   3.6509839 10.8111093 0.0000002
    ## hatch-inc.d9     9.9212123   6.2527397 13.5896849 0.0000000
    ## n5-inc.d9        7.1007012   3.4322286 10.7691738 0.0000008
    ## n9-inc.d9        2.8423526  -0.7377101  6.4224153 0.2346367
    ## hatch-inc.d17    2.6901657  -0.9783069  6.3586383 0.3347728
    ## n5-inc.d17      -0.1303454  -3.7988180  3.5381272 1.0000000
    ## n9-inc.d17      -4.3886940  -7.9687567 -0.8086313 0.0056372
    ## n5-hatch        -2.8205111  -6.5753125  0.9342903 0.3030569
    ## n9-hatch        -7.0788597 -10.7473323 -3.4103871 0.0000008
    ## n9-n5           -4.2583486  -7.9268212 -0.5898760 0.0110502

    ggplot(pcadata, aes(PC3, PC4,color = treatment)) + 
      geom_point(size = 2, alpha = 1) +
      stat_ellipse(type = "t") +
      xlab(paste0("PC3: ", percentVar[3],"% variance")) +
      ylab(paste0("PC4: ", percentVar[4],"% variance")) +
      theme_cowplot(font_size = 8, line_size = 0.25) 

![](../figures/pit/PCAmale-2.png)

    summary(aov(PC3 ~ treatment, data=pcadata)) 

    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## treatment    8  299.4   37.43   2.883 0.00673 **
    ## Residuals   88 1142.3   12.98                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    summary(aov(PC4 ~ treatment, data=pcadata)) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## treatment    8  175.9   21.98   2.089 0.0452 *
    ## Residuals   88  926.0   10.52                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    summary(aov(PC5 ~ treatment, data=pcadata)) 

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    8  127.2  15.898   1.731  0.102
    ## Residuals   88  808.1   9.183

    summary(aov(PC6 ~ treatment, data=pcadata))   

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    8   96.0  12.001   1.594  0.138
    ## Residuals   88  662.4   7.527

Heatmaps

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
             fontsize = 6)  

![](../figures/pit/heatmapMale-1.png)

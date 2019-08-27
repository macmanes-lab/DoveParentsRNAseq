    library(tidyverse)
    library(DESeq2)
    library(cowplot)
    library(RColorBrewer)
    library(pheatmap)
    library(kableExtra)
    library(viridis)
    library(forcats)
    library("wesanderson")

    library("BiocParallel")
    register(MulticoreParam(4))

    source("../R/functions.R")  # load custom functions 
    source("../R/themes.R")  # load custom themes and color palletes

    knitr::opts_chunk$set(fig.path = '../figures/sexes/', cache = TRUE)

Starting with all the data
--------------------------

    # import "colData" which contains sample information and "countData" which contains read counts
    c.colData <- read.csv("../metadata/00_colData_characterization.csv", header = T, row.names = 1)
    c.countData <- read.csv("../results/00_countData_characterization.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)

    # set levels
    c.colData$treatment <- factor(c.colData$treatment, levels = 
                                  c("control",  "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", "hatch", "n5", "n9"))
    levels(c.colData$treatment)

    ## [1] "control" "bldg"    "lay"     "inc.d3"  "inc.d9"  "inc.d17" "hatch"  
    ## [8] "n5"      "n9"

    c.colData$sextissue <- as.factor(paste(c.colData$sex, c.colData$tissue, sep = "_"))
    summary(c.colData[c(7,3,4,5,8)])

    ##              study         sex               tissue      treatment  
    ##  charcterization:576   female:289   gonad       :194   control: 73  
    ##                        male  :287   hypothalamus:189   inc.d9 : 71  
    ##                                     pituitary   :193   inc.d17: 66  
    ##                                                        n9     : 66  
    ##                                                        bldg   : 60  
    ##                                                        lay    : 60  
    ##                                                        (Other):180  
    ##                sextissue 
    ##  female_gonad       :98  
    ##  female_hypothalamus:95  
    ##  female_pituitary   :96  
    ##  male_gonad         :96  
    ##  male_hypothalamus  :94  
    ##  male_pituitary     :97  
    ## 

Run DESeq on all subsets of the data
------------------------------------

    dds.hypothalamus <- subsetDESeq2(c.colData,  c.countData, c("female_hypothalamus","male_hypothalamus") )

    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 189 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(189): L.Blu13_male_hypothalamus_control.NYNO
    ##   L.G107_male_hypothalamus_control ...
    ##   y97.x_female_hypothalamus_n9 y98.o50.x_male_hypothalamus_inc.d3
    ## colData names(8): V1 bird ... study sextissue
    ## [1] 14597   189

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 4 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 4 workers

    ## -- replacing outliers and refitting for 6 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    dds.pituitary <- subsetDESeq2(c.colData,  c.countData, c("female_pituitary", "male_pituitary" ) )

    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 193 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(193): L.Blu13_male_pituitary_control.NYNO
    ##   L.G107_male_pituitary_control ... y97.x_female_pituitary_n9
    ##   y98.o50.x_male_pituitary_inc.d3
    ## colData names(8): V1 bird ... study sextissue
    ## [1] 14484   193

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 4 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 4 workers

    ## -- replacing outliers and refitting for 28 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    dds.gonad <- subsetDESeq2(c.colData,  c.countData, c("female_gonad", "male_gonad") )

    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 194 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(194): L.Blu13_male_gonad_control.NYNO
    ##   L.G107_male_gonad_control ... y97.x_female_gonad_n9
    ##   y98.o50.x_male_gonad_inc.d3
    ## colData names(8): V1 bird ... study sextissue
    ## [1] 14843   194

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 4 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 4 workers

    ## -- replacing outliers and refitting for 8 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    hypPCA <- returnPCAs2(dds.hypothalamus)

    ## [1] 20 10  7  6  3  3
    ##                Df Sum Sq Mean Sq F value Pr(>F)    
    ## treatment       8   4783   597.8  17.824 <2e-16 ***
    ## sex             1     45    44.7   1.333  0.250    
    ## treatment:sex   8     53     6.6   0.197  0.991    
    ## Residuals     171   5735    33.5                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment       8   1657  207.09  10.524 5.86e-12 ***
    ## sex             1     33   33.14   1.684    0.196    
    ## treatment:sex   8    119   14.93   0.759    0.640    
    ## Residuals     171   3365   19.68                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment       8    185   23.16   1.128  0.347
    ## sex             1     24   23.82   1.160  0.283
    ## treatment:sex   8    147   18.38   0.895  0.522
    ## Residuals     171   3510   20.53               
    ##                Df Sum Sq Mean Sq F value Pr(>F)    
    ## treatment       8  137.2    17.2   1.756 0.0890 .  
    ## sex             1 1064.2  1064.2 108.914 <2e-16 ***
    ## treatment:sex   8  154.8    19.4   1.981 0.0515 .  
    ## Residuals     171 1670.8     9.8                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    plotPC12(hypPCA, "female and male hypothalamus")

![](../figures/sexes/pca-1.png)

    pitPCA <- returnPCAs2(dds.pituitary)

    ## [1] 13  9  7  5  4  4
    ##                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment       8   3851   481.4  49.296  < 2e-16 ***
    ## sex             1    627   626.6  64.172 1.54e-13 ***
    ## treatment:sex   8     60     7.5   0.771    0.628    
    ## Residuals     175   1709     9.8                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value Pr(>F)    
    ## treatment       8 2276.0   284.5  45.986 <2e-16 ***
    ## sex             1 1007.5  1007.5 162.847 <2e-16 ***
    ## treatment:sex   8   91.1    11.4   1.842 0.0723 .  
    ## Residuals     175 1082.7     6.2                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value Pr(>F)    
    ## treatment       8 1123.8   140.5   47.77 <2e-16 ***
    ## sex             1 1979.9  1979.9  673.36 <2e-16 ***
    ## treatment:sex   8   36.0     4.5    1.53   0.15    
    ## Residuals     175  514.6     2.9                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value Pr(>F)  
    ## treatment       8  219.9   27.49   2.278 0.0242 *
    ## sex             1   72.6   72.56   6.013 0.0152 *
    ## treatment:sex   8   86.7   10.84   0.898 0.5192  
    ## Residuals     175 2111.7   12.07                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    plotPC12(pitPCA, "female and male pituitary")

![](../figures/sexes/pca-2.png)

    gonPCA <- returnPCAs2(dds.gonad)

    ## [1] 91  3  1  1  0  0
    ##                Df Sum Sq Mean Sq   F value   Pr(>F)    
    ## treatment       8    314      39     3.398 0.001182 ** 
    ## sex             1 433334  433334 37453.286  < 2e-16 ***
    ## treatment:sex   8    334      42     3.606 0.000662 ***
    ## Residuals     176   2036      12                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment       8   2349  293.60   5.181 8.06e-06 ***
    ## sex             1      7    7.32   0.129     0.72    
    ## treatment:sex   8   2340  292.53   5.162 8.50e-06 ***
    ## Residuals     176   9974   56.67                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value  Pr(>F)    
    ## treatment       8    914  114.31   3.879 0.00031 ***
    ## sex             1      1    0.81   0.028 0.86842    
    ## treatment:sex   8    655   81.88   2.779 0.00643 ** 
    ## Residuals     176   5187   29.47                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value Pr(>F)  
    ## treatment       8  208.4  26.047   1.994 0.0496 *
    ## sex             1    1.5   1.526   0.117 0.7329  
    ## treatment:sex   8  173.2  21.655   1.658 0.1117  
    ## Residuals     176 2298.6  13.060                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    plotPC12(gonPCA, "female and male gonad")

![](../figures/sexes/pca-3.png)

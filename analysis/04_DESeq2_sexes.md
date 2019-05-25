    library(tidyverse)
    library(DESeq2)
    library(cowplot)
    library(RColorBrewer)
    library(pheatmap)
    library(kableExtra)
    library(viridis)

    library("BiocParallel")
    register(MulticoreParam(4))

    source("../R/functions.R")  # load custom functions 
    source("../R/themes.R")  # load custom themes and color palletes

    knitr::opts_chunk$set(fig.path = '../figures/sexes/', cache = TRUE)

Starting with all the data
--------------------------

    # import "colData" which contains sample information and "countData" which contains read counts
    a.colData <- read.csv("../metadata/00_samples.csv", header = T, row.names = 1)
    a.countData <- read.csv("../results/00_counts.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)

    # set levels
    levels(a.colData$treatment)

    ##  [1] "bldg"      "control"   "extend"    "hatch"     "inc.d17"  
    ##  [6] "inc.d3"    "inc.d9"    "lay"       "m.inc.d17" "m.inc.d3" 
    ## [11] "m.inc.d8"  "m.inc.d9"  "m.n2"      "n5"        "n9"       
    ## [16] "prolong"

    a.colData$treatment <- factor(a.colData$treatment, levels = 
                                 c("control", "bldg", "lay",
                                   "inc.d3", "m.inc.d3", 
                                   "inc.d9", "m.inc.d8", "m.inc.d9",
                                   "inc.d17", "m.inc.d17",
                                   "hatch",  "m.n2"  ,
                                   "n5", "prolong", "extend", "n9" ))
                                   
                                   
    a.colData$sextissue <- as.factor(paste(a.colData$sex, a.colData$tissue, sep = "_"))

    a.colData$lastday <- ifelse(grepl("m.inc.d3|m.inc.d9|m.inc.d17|m.n2", a.colData$treatment), "empty nest", 
                        ifelse(grepl("m.inc.d8|hatch|extend", a.colData$treatment),"chicks hatch",      
                         ifelse(grepl("n5", a.colData$treatment),"chicks early",
                         ifelse(grepl("n9", a.colData$treatment),"chicks later",
                         ifelse(grepl("control", a.colData$treatment),"control",
                         ifelse(grepl("bldg", a.colData$treatment),"nest building",
                         ifelse(grepl("prolong", a.colData$treatment),"eggs delay",
                          ifelse(grepl("lay", a.colData$treatment),"eggs lay",
                          ifelse(grepl("inc.d3", a.colData$treatment),"eggs early",
                         ifelse(grepl("inc.d9", a.colData$treatment),"eggs middle",
                        ifelse(grepl("inc.d17", a.colData$treatment),"eggs later", NA)))))))))))

    a.colData$penultimate <-  ifelse(grepl("extend", a.colData$treatment),"eggs delay",  
                              ifelse(grepl("m.n2", a.colData$treatment),"chicks hatch",
                         ifelse(grepl("n5", a.colData$treatment),"chicks early",
                         ifelse(grepl("n9", a.colData$treatment),"chicks later",
                         ifelse(grepl("control", a.colData$treatment),"control",
                         ifelse(grepl("bldg|lay", a.colData$treatment),"nest building",
                         ifelse(grepl("prolong", a.colData$treatment),"eggs delay",
                          ifelse(grepl("inc.d3|m.inc.d3", a.colData$treatment),"eggs early",
                         ifelse(grepl("inc.d9|m.inc.d9|m.inc.d8", a.colData$treatment),"eggs middle",
                        ifelse(grepl("inc.d17|hatch", a.colData$treatment),"eggs later", NA))))))))))

    a.colData$xlabel <- a.colData$treatment

    levels(a.colData$xlabel) <-  c("control", "nest.building", "egg.lay",
                                   "eggs.early", "eggs.early.remove", 
                                   "eggs.mid", "eggs.mid.hatch", "eggs.mid.remove",
                                   "eggs.end", "eggs.end.remove",
                                   "chicks.hatch",  "chicks.hatch.remove"  ,
                                   "chicks.mid", "eggs.delay", "eggs.delay.hatch", "chicks.end")


    a.colData$lastday <- factor(a.colData$lastday, levels =  c("control", "nest building", "eggs lay",
                                   "eggs early", "eggs middle",  "eggs later","eggs delay",
                                   "chicks hatch", "chicks early",  "chicks later", "empty nest"))

    a.colData$penultimate <- factor(a.colData$penultimate, levels =  c("control", "nest building",
                                   "eggs early",  "eggs middle",  "eggs later",  "eggs delay",
                                   "chicks hatch", "chicks early",  "chicks later"))

    summary(a.colData[c(7,3,4,5,8,9, 10,11)])

    ##              study         sex               tissue        treatment  
    ##  charcterization:576   female:497   gonad       :330   control  : 73  
    ##  manipulation   :411   male  :490   hypothalamus:327   inc.d9   : 71  
    ##                                     pituitary   :330   inc.d17  : 66  
    ##                                                        n9       : 66  
    ##                                                        m.inc.d17: 63  
    ##                                                        bldg     : 60  
    ##                                                        (Other)  :588  
    ##                sextissue           lastday           penultimate 
    ##  female_gonad       :167   empty nest  :231   eggs later   :189  
    ##  female_hypothalamus:165   chicks hatch:180   eggs middle  :180  
    ##  female_pituitary   :165   control     : 73   nest building:120  
    ##  male_gonad         :163   eggs middle : 71   eggs early   :120  
    ##  male_hypothalamus  :162   eggs later  : 66   eggs delay   :120  
    ##  male_pituitary     :165   chicks later: 66   control      : 73  
    ##                            (Other)     :300   (Other)      :185  
    ##              xlabel   
    ##  control        : 73  
    ##  eggs.mid       : 71  
    ##  eggs.end       : 66  
    ##  chicks.end     : 66  
    ##  eggs.end.remove: 63  
    ##  nest.building  : 60  
    ##  (Other)        :588

Run DESeq on all subsets of the data
------------------------------------

    dds.female_hypothalamus <- subsetDESeq(a.colData, a.countData, "female_hypothalamus")

    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 165 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(165): L.G118_female_hypothalamus_control.NYNO
    ##   R.G106_female_hypothalamus_control ...
    ##   y97.x_female_hypothalamus_n9 y98.g54_female_hypothalamus_m.hatch
    ## colData names(11): V1 bird ... penultimate xlabel
    ## [1] 14576   165

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 9 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    dds.female_pituitary <- subsetDESeq(a.colData, a.countData, "female_pituitary" )

    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 165 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(165): L.G118_female_pituitary_control.NYNO
    ##   R.G106_female_pituitary_control ... y97.x_female_pituitary_n9
    ##   y98.g54_female_pituitary_m.hatch
    ## colData names(11): V1 bird ... penultimate xlabel
    ## [1] 14496   165

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 49 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    dds.female_gonad <- subsetDESeq(a.colData, a.countData, "female_gonad" )

    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 167 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(167): L.G118_female_gonad_control
    ##   R.G106_female_gonad_control ... y97.x_female_gonad_n9
    ##   y98.g54_female_gonad_m.hatch
    ## colData names(11): V1 bird ... penultimate xlabel
    ## [1] 14746   167

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 156 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    dds.male_hypothalamus <- subsetDESeq(a.colData, a.countData, "male_hypothalamus" )

    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 162 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(162): L.Blu13_male_hypothalamus_control.NYNO
    ##   L.G107_male_hypothalamus_control ...
    ##   y95.g131.x_male_hypothalamus_inc.d9
    ##   y98.o50.x_male_hypothalamus_inc.d3
    ## colData names(11): V1 bird ... penultimate xlabel
    ## [1] 14536   162

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

    dds.male_pituitary <- subsetDESeq(a.colData, a.countData, "male_pituitary"  )

    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 165 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(165): L.Blu13_male_pituitary_control.NYNO
    ##   L.G107_male_pituitary_control ...
    ##   y95.g131.x_male_pituitary_inc.d9 y98.o50.x_male_pituitary_inc.d3
    ## colData names(11): V1 bird ... penultimate xlabel
    ## [1] 14480   165

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 51 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    dds.male_gondad <- subsetDESeq(a.colData, a.countData, "male_gonad")

    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 163 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(163): L.Blu13_male_gonad_control.NYNO
    ##   L.G107_male_gonad_control ... y95.g131.x_male_gonad_inc.d9
    ##   y98.o50.x_male_gonad_inc.d3
    ## colData names(11): V1 bird ... penultimate xlabel
    ## [1] 14765   163

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 123 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

Calculate and plot total DEGs
-----------------------------

    #DEGs.female_hypothalamus <- returntotalDEGs(dds.female_hypothalamus)
    #DEGs.female_pituitary <- returntotalDEGs(dds.female_pituitary)
    #DEGs.female_gonad <- returntotalDEGs(dds.female_gonad)
    #DEGs.male_hypothalamus <- returntotalDEGs(dds.male_hypothalamus)
    DEGs.male_pituitary <- returntotalDEGs(dds.male_pituitary)

    ## [1] "control.bldg"
    ## [1] "control.inc.d17"
    ## [1] "control.m.inc.d17"
    ## [1] "control.hatch"
    ## [1] "control.m.n2"
    ## [1] "control.n9"
    ## [1] "bldg.inc.d17"
    ## [1] "bldg.m.inc.d17"
    ## [1] "bldg.hatch"
    ## [1] "bldg.m.n2"
    ## [1] "bldg.n9"
    ## [1] "inc.d17.m.inc.d17"
    ## [1] "inc.d17.hatch"
    ## [1] "inc.d17.m.n2"
    ## [1] "inc.d17.n9"
    ## [1] "m.inc.d17.hatch"
    ## [1] "m.inc.d17.m.n2"
    ## [1] "m.inc.d17.n9"
    ## [1] "hatch.m.n2"
    ## [1] "hatch.n9"
    ## [1] "m.n2.n9"
    ##                          V1        V2   V3
    ## control.bldg        control      bldg 4969
    ## control.inc.d17     control   inc.d17 4567
    ## control.m.inc.d17   control m.inc.d17 4893
    ## control.hatch       control     hatch 4877
    ## control.m.n2        control      m.n2 4815
    ## control.n9          control        n9 5020
    ## bldg.inc.d17           bldg   inc.d17  342
    ## bldg.m.inc.d17         bldg m.inc.d17  163
    ## bldg.hatch             bldg     hatch  840
    ## bldg.m.n2              bldg      m.n2  171
    ## bldg.n9                bldg        n9   38
    ## inc.d17.m.inc.d17   inc.d17 m.inc.d17  278
    ## inc.d17.hatch       inc.d17     hatch    9
    ## inc.d17.m.n2        inc.d17      m.n2  201
    ## inc.d17.n9          inc.d17        n9  166
    ## m.inc.d17.hatch   m.inc.d17     hatch  737
    ## m.inc.d17.m.n2    m.inc.d17      m.n2    4
    ## m.inc.d17.n9      m.inc.d17        n9    8
    ## hatch.m.n2            hatch      m.n2  688
    ## hatch.n9              hatch        n9  433
    ## m.n2.n9                m.n2        n9   19

    #DEGs.male_gondad <- returntotalDEGs(dds.male_gondad)

    #a <- plottotalDEGs(DEGs.female_hypothalamus, "female hypothalamus")
    #b <- plottotalDEGs(DEGs.female_pituitary, "female pituitary")
    #c <- plottotalDEGs(DEGs.female_gonad, "female gonad")
    #d <- plottotalDEGs(DEGs.male_hypothalamus, "male hypothalamus")
    e <- plottotalDEGs(DEGs.male_pituitary, "male pituitary")

    ##                          V1        V2   V3
    ## control.bldg        control      bldg 4969
    ## control.inc.d17     control   inc.d17 4567
    ## control.m.inc.d17   control m.inc.d17 4893
    ## control.hatch       control     hatch 4877
    ## control.m.n2        control      m.n2 4815
    ## control.n9          control        n9 5020
    ## bldg.inc.d17           bldg   inc.d17  342
    ## bldg.m.inc.d17         bldg m.inc.d17  163
    ## bldg.hatch             bldg     hatch  840
    ## bldg.m.n2              bldg      m.n2  171
    ## bldg.n9                bldg        n9   38
    ## inc.d17.m.inc.d17   inc.d17 m.inc.d17  278
    ## inc.d17.hatch       inc.d17     hatch    9
    ## inc.d17.m.n2        inc.d17      m.n2  201
    ## inc.d17.n9          inc.d17        n9  166
    ## m.inc.d17.hatch   m.inc.d17     hatch  737
    ## m.inc.d17.m.n2    m.inc.d17      m.n2    4
    ## m.inc.d17.n9      m.inc.d17        n9    8
    ## hatch.m.n2            hatch      m.n2  688
    ## hatch.n9              hatch        n9  433
    ## m.n2.n9                m.n2        n9   19

![](../figures/sexes/plotDEGs-1.png)

    #f <- plottotalDEGs(DEGs.male_gondad, "male gonad")

    plot_grid(#a + theme(legend.position = "none"),
              #b + theme(legend.position = "none"),
              #c + theme(legend.position = "none"),
              #d + theme(legend.position = "none"),
              e + theme(legend.position = "none"),
              #f + theme(legend.position = "none"),
              nrow = 2) 

![](../figures/sexes/plottotalDEGs-1.png)

Calculate and plot principal components
---------------------------------------

    pcas.female_hypothalamus <- returnPCAs(dds.female_hypothalamus)

    ## [1] 26  9  8  4  3  3
    ##              Df Sum Sq Mean Sq F value Pr(>F)    
    ## xlabel       15   6078   405.2   18.67 <2e-16 ***
    ## Residuals   149   3235    21.7                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value Pr(>F)  
    ## xlabel       15  490.5   32.70   1.684 0.0596 .
    ## Residuals   149 2892.8   19.41                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## xlabel       15  749.2   49.95   3.791 1.12e-05 ***
    ## Residuals   149 1963.1   13.18                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## xlabel       15  294.6   19.64   3.012 0.000306 ***
    ## Residuals   149  971.5    6.52                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    pcas.female_pituitary <- returnPCAs(dds.female_pituitary)      

    ## [1] 11  8  7  6  4  3
    ##              Df Sum Sq Mean Sq F value Pr(>F)    
    ## xlabel       15   2480  165.31   13.99 <2e-16 ***
    ## Residuals   149   1761   11.82                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## xlabel       15   1076   71.70   5.318 1.83e-08 ***
    ## Residuals   149   2009   13.48                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value Pr(>F)    
    ## xlabel       15   1342   89.46   10.26 <2e-16 ***
    ## Residuals   149   1300    8.72                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value Pr(>F)    
    ## xlabel       15 1195.8   79.72   12.02 <2e-16 ***
    ## Residuals   149  988.5    6.63                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    pcas.female_gonad <- returnPCAs(dds.female_gonad)

    ## [1] 37 12  7  4  4  3
    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## xlabel       15  11694   779.6   3.644 2.04e-05 ***
    ## Residuals   151  32303   213.9                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## xlabel       15   3103  206.83   2.799 0.000738 ***
    ## Residuals   151  11160   73.91                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value Pr(>F)
    ## xlabel       15    803   53.56   1.109  0.353
    ## Residuals   151   7294   48.31               
    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## xlabel       15   1487   99.11   4.148 2.37e-06 ***
    ## Residuals   151   3607   23.89                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    pcas.male_hypothalamus <- returnPCAs(dds.male_hypothalamus)

    ## [1] 28  9  6  4  3  3
    ##              Df Sum Sq Mean Sq F value Pr(>F)    
    ## xlabel       15   5836   389.1   10.54 <2e-16 ***
    ## Residuals   146   5388    36.9                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## xlabel       15   1230   82.00   4.879 1.21e-07 ***
    ## Residuals   146   2454   16.81                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value Pr(>F)
    ## xlabel       15  222.1   14.81   0.959  0.501
    ## Residuals   146 2253.7   15.44               
    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## xlabel       15  419.1  27.938    3.23 0.000125 ***
    ## Residuals   146 1263.0   8.651                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    pcas.male_pituitary <- returnPCAs(dds.male_pituitary)

    ## [1] 11  9  7  5  4  3
    ##              Df Sum Sq Mean Sq F value Pr(>F)    
    ## xlabel       15 2676.7  178.44   30.47 <2e-16 ***
    ## Residuals   149  872.6    5.86                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value  Pr(>F)    
    ## xlabel       15   1412   94.15   9.431 2.8e-15 ***
    ## Residuals   149   1488    9.98                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## xlabel       15  762.4   50.82   4.975 7.58e-08 ***
    ## Residuals   149 1522.3   10.22                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## xlabel       15  549.6   36.64   4.769 1.79e-07 ***
    ## Residuals   149 1144.8    7.68                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    pcas.male_gondad <- returnPCAs(dds.male_gondad)

    ## [1] 14 10  5  3  3  2
    ##              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## xlabel       15   1476   98.37   3.674 1.89e-05 ***
    ## Residuals   147   3936   26.77                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value  Pr(>F)   
    ## xlabel       15  830.7   55.38   2.679 0.00124 **
    ## Residuals   147 3038.4   20.67                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value Pr(>F)    
    ## xlabel       15 1459.0   97.27   22.58 <2e-16 ***
    ## Residuals   147  633.3    4.31                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##              Df Sum Sq Mean Sq F value  Pr(>F)   
    ## xlabel       15  235.5   15.70   2.488 0.00273 **
    ## Residuals   147  927.6    6.31                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    p1 <- plotPCAs(pcas.female_hypothalamus, "female hypothalamus") 

![](../figures/sexes/plotpca-1.png)![](../figures/sexes/plotpca-2.png)![](../figures/sexes/plotpca-3.png)

    p2 <- plotPCAs(pcas.female_pituitary, "female pituitary")   

![](../figures/sexes/plotpca-4.png)![](../figures/sexes/plotpca-5.png)![](../figures/sexes/plotpca-6.png)

    p3 <- plotPCAs(pcas.female_gonad, "female gonad") 

![](../figures/sexes/plotpca-7.png)![](../figures/sexes/plotpca-8.png)![](../figures/sexes/plotpca-9.png)

    p4 <- plotPCAs(pcas.male_hypothalamus, "male hypothalamus") 

![](../figures/sexes/plotpca-10.png)![](../figures/sexes/plotpca-11.png)![](../figures/sexes/plotpca-12.png)

    p5 <- plotPCAs(pcas.male_pituitary, "male pituitary") 

![](../figures/sexes/plotpca-13.png)![](../figures/sexes/plotpca-14.png)![](../figures/sexes/plotpca-15.png)

    p6 <- plotPCAs(pcas.male_gondad, "male gonad") 

![](../figures/sexes/plotpca-16.png)![](../figures/sexes/plotpca-17.png)![](../figures/sexes/plotpca-18.png)

    p1 + theme_noaxislabels

![](../figures/sexes/plotpca-19.png)

    p2 + theme_noaxislabels

![](../figures/sexes/plotpca-20.png)

    p3 + theme_noaxislabels

![](../figures/sexes/plotpca-21.png)

    p4 + theme_noaxislabels

![](../figures/sexes/plotpca-22.png)

    p5 + theme_noaxislabels

![](../figures/sexes/plotpca-23.png)

    p6 + theme_noaxislabels

![](../figures/sexes/plotpca-24.png)

    mylegend <- get_legend(p1)

    allPC1s <- plot_grid(p1 + theme_noaxislabels, 
                         p2 + theme_noaxislabels, 
                         p3 + theme_noaxislabels,
                         p4 + theme_noaxislabels,
                         p5 + theme_noaxislabels,
                         p6 + theme_noaxislabels,
                         nrow = 3, rel_heights = c(0.3, 0.3, 0.4))

    pc1 <- plot_grid(allPC1s,mylegend,  nrow = 2, rel_heights = c(1.0, 0.15))
    pc1

    pdf("../figures/sexes/pca-1.pdf", width = 12, height = 10)
    plot(pc1)
    dev.off()

heamap with minimum pvalue
--------------------------

    makepheatmap(dds.female_hypothalamus, a.colData, "female hypothalamus")
    makepheatmap(dds.female_pituitary, "female pituitary")
    makepheatmap(dds.female_gonad, "female gonad")
    makepheatmap(dds.male_hypothalamus, "male hypothalamus")
    makepheatmap(dds.male_pituitary, "male pituitary")
    makepheatmap(dds.male_gondad, "male gonad")        

candidate genes
---------------

    plotcandidates(dds.female_hypothalamus, a.colData, "female hypothalamus")
    plotcandidates(dds.female_pituitary, "female pituitary")
    plotcandidates(dds.female_gonad, "female gonad")
    plotcandidates(dds.male_hypothalamus, "male hypothalamus")
    plotcandidates(dds.male_pituitary, "male pituitary")
    plotcandidates(dds.male_gondad, "male gonad")

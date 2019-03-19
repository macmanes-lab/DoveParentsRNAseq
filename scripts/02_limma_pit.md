DESeq2 is *not* recommended for experiments with more than 100 samples
([see Mike Love’s
post](https://mikelove.wordpress.com/2016/09/28/deseq2-or-edger/)), so I
decided to try the limma package. I followed [this
tutorial](https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html).

    library(tidyverse)

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.0       ✔ purrr   0.3.1  
    ## ✔ tibble  2.0.1       ✔ dplyr   0.8.0.1
    ## ✔ tidyr   0.8.3       ✔ stringr 1.4.0  
    ## ✔ readr   1.3.1       ✔ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

    library(limma)
    library(Glimma)
    library(edgeR)
    library(kableExtra)

    ## 
    ## Attaching package: 'kableExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    library(ggplot2)

    knitr::opts_chunk$set(fig.path = '../figures/pit/',cache=TRUE)

First, I read in the data I processed in 00\_datawrangling.Rmd.

    # import "colData" which contains sample information and "countData" which contains read counts
    colData <- read.csv("../results/00_colData_characterization.csv", header = T, row.names = 1)
    countData <- read.csv("../results/00_countData_characterization.csv", header = T, row.names = 1)
    head(colData)

    ##                                                                            V1
    ## L.Blu13_male_gonad_control.NYNO               L.Blu13_male_gonad_control.NYNO
    ## L.Blu13_male_hypothalamus_control.NYNO L.Blu13_male_hypothalamus_control.NYNO
    ## L.Blu13_male_pituitary_control.NYNO       L.Blu13_male_pituitary_control.NYNO
    ## L.G107_male_gonad_control                           L.G107_male_gonad_control
    ## L.G107_male_hypothalamus_control             L.G107_male_hypothalamus_control
    ## L.G107_male_pituitary_control                   L.G107_male_pituitary_control
    ##                                           bird  sex       tissue NYNO
    ## L.Blu13_male_gonad_control.NYNO        L.Blu13 male        gonad NYNO
    ## L.Blu13_male_hypothalamus_control.NYNO L.Blu13 male hypothalamus NYNO
    ## L.Blu13_male_pituitary_control.NYNO    L.Blu13 male    pituitary NYNO
    ## L.G107_male_gonad_control               L.G107 male        gonad <NA>
    ## L.G107_male_hypothalamus_control        L.G107 male hypothalamus <NA>
    ## L.G107_male_pituitary_control           L.G107 male    pituitary <NA>
    ##                                        treatment                     group
    ## L.Blu13_male_gonad_control.NYNO          control        male.gonad.control
    ## L.Blu13_male_hypothalamus_control.NYNO   control male.hypothalamus.control
    ## L.Blu13_male_pituitary_control.NYNO      control    male.pituitary.control
    ## L.G107_male_gonad_control                control        male.gonad.control
    ## L.G107_male_hypothalamus_control         control male.hypothalamus.control
    ## L.G107_male_pituitary_control            control    male.pituitary.control

    geneinfo <- read.csv("../results/00_geneinfo.csv", row.names = 1)

    colData <- colData %>%
      dplyr::filter(grepl('pituitary', tissue)) %>%
      droplevels()
    row.names(colData) <- colData$V1

    # print sample sizes
    colData %>% select(sex,treatment, tissue)  %>%  summary()

    ##      sex       treatment        tissue   
    ##  female:96   control:25   pituitary:193  
    ##  male  :97   inc.d9 :24                  
    ##              inc.d17:22                  
    ##              n9     :22                  
    ##              bldg   :20                  
    ##              hatch  :20                  
    ##              (Other):60

    savecols <- as.character(colData$V1) 
    savecols <- as.vector(savecols) 
    countData <- countData %>% dplyr::select(one_of(savecols)) 

    # check that row and col lenghts are equal
    ncol(countData) == nrow(colData)

    ## [1] TRUE

Then, I followed the steps from
<a href="https://github.com/macmanes-lab/RockDove/blob/master/parental_care/parental_analysis.Rmd" class="uri">https://github.com/macmanes-lab/RockDove/blob/master/parental_care/parental_analysis.Rmd</a>.

    # create a large DGEList with 3 elements
    parentalobject <- DGEList(counts=countData, genes=geneinfo, group=colData$group)

    # transform raw counts to countspermillion
    cpms <- cpm(parentalobject)

    # calculate number of lowly lowly expressed genes and remove them
    table(rowSums(parentalobject$counts==0)==10)

    ## 
    ## FALSE  TRUE 
    ## 14877    60

    keep_genes <- rowSums(cpms >= 1) >= 10
    dge <- parentalobject[keep_genes, ]

    # specific the design
    parentaldesign <- model.matrix(~ colData$group )
    colnames(parentaldesign) <- levels(colData$group)

    # The TMM normalization
    parentalobject <- calcNormFactors(parentalobject)
    parentalobject <- estimateCommonDisp(parentalobject)
    parentalobject <- estimateTagwiseDisp(parentalobject)
    parentalobject <- estimateDisp(parentalobject, parentaldesign)
    parentalobject <- estimateGLMCommonDisp(parentalobject, parentaldesign, verbose=TRUE)

    ## Disp = 0.08405 , BCV = 0.2899

    parentalobject <- estimateGLMTrendedDisp(parentalobject, parentaldesign)
    parentalobject <- estimateGLMTagwiseDisp(parentalobject, parentaldesign)

    #  perform likelihood ratio test and thresholded testing
    fit <- glmFit( parentalobject, parentaldesign, robust=T)
    tr <- glmTreat(fit, lfc = 1)
    topTags(tr)

    ## Coefficient:  male.pituitary.n9 
    ##                row.names         Name    geneid       entrezid     logFC
    ## XP_015135776.1    427254       ZFAND5    427254 XP_015135776.1 -6.197065
    ## XP_001234565.2    771273       MRPS36    771273 XP_001234565.2 -7.670507
    ## XP_004944394.1 101750188 LOC101750188 101750188 XP_004944394.1 -7.828419
    ## NP_001264524.1    415770  C11H19ORF40    415770 NP_001264524.1 -8.090663
    ## XP_015135862.1    427458       HNRNPK    427458 XP_015135862.1 -6.744749
    ## XP_015132890.1    770140         CTIF    770140 XP_015132890.1 -5.270510
    ## XP_015131806.1    448833          BTC    448833 XP_015131806.1 -3.611849
    ## NP_990826.1       396491       LGALS1    396491    NP_990826.1 -5.398852
    ## XP_004937386.1    427645        PDE8B    427645 XP_004937386.1  2.521983
    ## NP_001278710.1    417513        RASA4    417513 NP_001278710.1 -3.107554
    ##                unshrunk.logFC     logCPM       PValue          FDR
    ## XP_015135776.1      -6.232282  4.3516310 1.060502e-30 1.584072e-26
    ## XP_001234565.2      -8.453837  1.8913376 2.841630e-30 2.122271e-26
    ## XP_004944394.1      -8.194004  2.7337818 1.031776e-29 5.137213e-26
    ## NP_001264524.1      -8.535391  2.8561678 1.460167e-26 5.452628e-23
    ## XP_015135862.1      -7.525839  1.0696834 1.383062e-22 4.131760e-19
    ## XP_015132890.1      -6.540394 -0.8009548 1.455972e-12 3.624643e-09
    ## XP_015131806.1      -3.621634  3.0316585 2.754836e-11 5.878427e-08
    ## NP_990826.1         -5.428065  2.9709474 6.272613e-11 1.171175e-07
    ## XP_004937386.1       2.524241  4.7704921 4.020179e-09 6.672157e-06
    ## NP_001278710.1      -3.109632  4.3810017 1.642355e-08 2.453186e-05

plotMDS (multidimential scaling)
================================

    plotMDS(parentalobject, cex = 0.5)

![](../figures/pit/plotMDS-1.png)

    levels(colData$treatment)

    ## [1] "bldg"    "control" "hatch"   "inc.d17" "inc.d3"  "inc.d9"  "lay"    
    ## [8] "n5"      "n9"

    col.treatment <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6")[colData$treatment]

    plotMDS(parentalobject,col=col.treatment, labels = colData$sex)
    legend("bottomright",fill=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6"),legend=levels(colData$treatment))
    title("Pituitary Colored by Treatment")

![](../figures/pit/plotMDS-lables-1.png)

    plotMDS(parentalobject,dim=c(3,4), col=col.treatment, labels = colData$sex)
    legend("bottomright",fill=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6"),legend=levels(colData$treatment))
    title("Pituitary Colored by Treatment")

![](../figures/pit/plotMDS-lables-2.png)

specify contrasts and make MA plots (currently only a subset)
=============================================================

    # view all levels
    levels(colData$group)

    ##  [1] "female.pituitary.bldg"    "female.pituitary.control"
    ##  [3] "female.pituitary.hatch"   "female.pituitary.inc.d17"
    ##  [5] "female.pituitary.inc.d3"  "female.pituitary.inc.d9" 
    ##  [7] "female.pituitary.lay"     "female.pituitary.n5"     
    ##  [9] "female.pituitary.n9"      "male.pituitary.bldg"     
    ## [11] "male.pituitary.control"   "male.pituitary.hatch"    
    ## [13] "male.pituitary.inc.d17"   "male.pituitary.inc.d3"   
    ## [15] "male.pituitary.inc.d9"    "male.pituitary.lay"      
    ## [17] "male.pituitary.n5"        "male.pituitary.n9"

    # subset of conrasts - sex specific comparing hatch to lay
    my.contrasts <- makeContrasts(
                 FP_CB = female.pituitary.control - female.pituitary.bldg,
                 FP_BL = female.pituitary.bldg - female.pituitary.lay,
                 FP_Li3 = female.pituitary.lay - female.pituitary.inc.d3,
                 FP_i39 = female.pituitary.inc.d3 - female.pituitary.inc.d9,
                 FP_i917 = female.pituitary.inc.d9 - female.pituitary.inc.d17,
                 FP_i17H = female.pituitary.inc.d17 - female.pituitary.hatch,
                 FP_H5 = female.pituitary.hatch -  female.pituitary.n5,
                 FP_n59 = female.pituitary.n5 - female.pituitary.n9,
                 FP_n9C = female.pituitary.n9 - female.pituitary.control,
    levels=parentaldesign)

    cont <- "FP_CB"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*female.pituitary.bldg 1*female.pituitary.control
    ## Down                                                     0
    ## NotSig                                                 344
    ## Up                                                   14593

    kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)

<table class="kable_wrapper">
<tbody>
<tr>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
row.names
</th>
<th style="text-align:left;">
Name
</th>
<th style="text-align:right;">
geneid
</th>
<th style="text-align:left;">
entrezid
</th>
<th style="text-align:right;">
logFC
</th>
<th style="text-align:right;">
unshrunk.logFC
</th>
<th style="text-align:right;">
logCPM
</th>
<th style="text-align:right;">
PValue
</th>
<th style="text-align:right;">
FDR
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
XP\_004951164.1
</td>
<td style="text-align:right;">
101751744
</td>
<td style="text-align:left;">
SLC35D2
</td>
<td style="text-align:right;">
101751744
</td>
<td style="text-align:left;">
XP\_004951164.1
</td>
<td style="text-align:right;">
20.39
</td>
<td style="text-align:right;">
20.41
</td>
<td style="text-align:right;">
1.15
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_004942079.1
</td>
<td style="text-align:right;">
423703
</td>
<td style="text-align:left;">
NUDT13
</td>
<td style="text-align:right;">
423703
</td>
<td style="text-align:left;">
XP\_004942079.1
</td>
<td style="text-align:right;">
19.69
</td>
<td style="text-align:right;">
19.70
</td>
<td style="text-align:right;">
1.23
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001186323.1
</td>
<td style="text-align:right;">
416716
</td>
<td style="text-align:left;">
CST7
</td>
<td style="text-align:right;">
416716
</td>
<td style="text-align:left;">
NP\_001186323.1
</td>
<td style="text-align:right;">
19.53
</td>
<td style="text-align:right;">
19.54
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_001233261.2
</td>
<td style="text-align:right;">
107049048
</td>
<td style="text-align:left;">
CMC4
</td>
<td style="text-align:right;">
107049048
</td>
<td style="text-align:left;">
XP\_001233261.2
</td>
<td style="text-align:right;">
19.48
</td>
<td style="text-align:right;">
19.49
</td>
<td style="text-align:right;">
1.13
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015131033.1
</td>
<td style="text-align:right;">
416881
</td>
<td style="text-align:left;">
TMEM116
</td>
<td style="text-align:right;">
416881
</td>
<td style="text-align:left;">
XP\_015131033.1
</td>
<td style="text-align:right;">
19.35
</td>
<td style="text-align:right;">
19.37
</td>
<td style="text-align:right;">
1.33
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BH
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
-1*female.pituitary.bldg 1*female.pituitary.control
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
glm
</td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FP_CB', frame.plot=F)

![](../figures/pit/01-contrasts-1.png)

    cont <- "FP_BL"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*female.pituitary.bldg -1*female.pituitary.lay
    ## Down                                             14678
    ## NotSig                                             259
    ## Up                                                   0

    kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)

<table class="kable_wrapper">
<tbody>
<tr>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
row.names
</th>
<th style="text-align:left;">
Name
</th>
<th style="text-align:right;">
geneid
</th>
<th style="text-align:left;">
entrezid
</th>
<th style="text-align:right;">
logFC
</th>
<th style="text-align:right;">
unshrunk.logFC
</th>
<th style="text-align:right;">
logCPM
</th>
<th style="text-align:right;">
PValue
</th>
<th style="text-align:right;">
FDR
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
NP\_001001755.1
</td>
<td style="text-align:right;">
414837
</td>
<td style="text-align:left;">
THBS2
</td>
<td style="text-align:right;">
414837
</td>
<td style="text-align:left;">
NP\_001001755.1
</td>
<td style="text-align:right;">
-20.52
</td>
<td style="text-align:right;">
-20.54
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_001233219.1
</td>
<td style="text-align:right;">
769909
</td>
<td style="text-align:left;">
TNFRSF11A
</td>
<td style="text-align:right;">
769909
</td>
<td style="text-align:left;">
XP\_001233219.1
</td>
<td style="text-align:right;">
-20.15
</td>
<td style="text-align:right;">
-20.17
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_004951164.1
</td>
<td style="text-align:right;">
101751744
</td>
<td style="text-align:left;">
SLC35D2
</td>
<td style="text-align:right;">
101751744
</td>
<td style="text-align:left;">
XP\_004951164.1
</td>
<td style="text-align:right;">
-19.98
</td>
<td style="text-align:right;">
-20.00
</td>
<td style="text-align:right;">
1.15
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015155022.1
</td>
<td style="text-align:right;">
771869
</td>
<td style="text-align:left;">
FBXL20
</td>
<td style="text-align:right;">
771869
</td>
<td style="text-align:left;">
XP\_015155022.1
</td>
<td style="text-align:right;">
-19.83
</td>
<td style="text-align:right;">
-19.84
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001025839.2
</td>
<td style="text-align:right;">
416967
</td>
<td style="text-align:left;">
HNF1A
</td>
<td style="text-align:right;">
416967
</td>
<td style="text-align:left;">
NP\_001025839.2
</td>
<td style="text-align:right;">
-19.79
</td>
<td style="text-align:right;">
-19.81
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BH
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1*female.pituitary.bldg -1*female.pituitary.lay
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
glm
</td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FP_BL', frame.plot=F)

![](../figures/pit/01-contrasts-2.png)

    cont <- "FP_Li3"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*female.pituitary.inc.d3 1*female.pituitary.lay
    ## Down                                                   0
    ## NotSig                                             14937
    ## Up                                                     0

    kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)

<table class="kable_wrapper">
<tbody>
<tr>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
row.names
</th>
<th style="text-align:left;">
Name
</th>
<th style="text-align:right;">
geneid
</th>
<th style="text-align:left;">
entrezid
</th>
<th style="text-align:right;">
logFC
</th>
<th style="text-align:right;">
unshrunk.logFC
</th>
<th style="text-align:right;">
logCPM
</th>
<th style="text-align:right;">
PValue
</th>
<th style="text-align:right;">
FDR
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
XP\_015131806.1
</td>
<td style="text-align:right;">
448833
</td>
<td style="text-align:left;">
BTC
</td>
<td style="text-align:right;">
448833
</td>
<td style="text-align:left;">
XP\_015131806.1
</td>
<td style="text-align:right;">
-2.34
</td>
<td style="text-align:right;">
-2.35
</td>
<td style="text-align:right;">
3.03
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.00
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015133687.1
</td>
<td style="text-align:right;">
396214
</td>
<td style="text-align:left;">
PLP1
</td>
<td style="text-align:right;">
396214
</td>
<td style="text-align:left;">
XP\_015133687.1
</td>
<td style="text-align:right;">
7.06
</td>
<td style="text-align:right;">
7.07
</td>
<td style="text-align:right;">
5.46
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.00
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015152637.1
</td>
<td style="text-align:right;">
419443
</td>
<td style="text-align:left;">
TMEM201
</td>
<td style="text-align:right;">
419443
</td>
<td style="text-align:left;">
XP\_015152637.1
</td>
<td style="text-align:right;">
1.05
</td>
<td style="text-align:right;">
1.05
</td>
<td style="text-align:right;">
5.77
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.00
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015153698.1
</td>
<td style="text-align:right;">
419759
</td>
<td style="text-align:left;">
ZBTB16
</td>
<td style="text-align:right;">
419759
</td>
<td style="text-align:left;">
XP\_015153698.1
</td>
<td style="text-align:right;">
1.57
</td>
<td style="text-align:right;">
1.57
</td>
<td style="text-align:right;">
3.72
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001025732.1
</td>
<td style="text-align:right;">
415650
</td>
<td style="text-align:left;">
PLLP
</td>
<td style="text-align:right;">
415650
</td>
<td style="text-align:left;">
NP\_001025732.1
</td>
<td style="text-align:right;">
3.45
</td>
<td style="text-align:right;">
3.46
</td>
<td style="text-align:right;">
2.41
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BH
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
-1*female.pituitary.inc.d3 1*female.pituitary.lay
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
glm
</td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FP_Li3', frame.plot=F)

![](../figures/pit/01-contrasts-3.png)

    cont <- "FP_i39"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*female.pituitary.inc.d3 -1*female.pituitary.inc.d9
    ## Down                                                      0
    ## NotSig                                                14937
    ## Up                                                        0

    kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)

<table class="kable_wrapper">
<tbody>
<tr>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
row.names
</th>
<th style="text-align:left;">
Name
</th>
<th style="text-align:right;">
geneid
</th>
<th style="text-align:left;">
entrezid
</th>
<th style="text-align:right;">
logFC
</th>
<th style="text-align:right;">
unshrunk.logFC
</th>
<th style="text-align:right;">
logCPM
</th>
<th style="text-align:right;">
PValue
</th>
<th style="text-align:right;">
FDR
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
XP\_015152100.1
</td>
<td style="text-align:right;">
107054855
</td>
<td style="text-align:left;">
LOC107054855
</td>
<td style="text-align:right;">
107054855
</td>
<td style="text-align:left;">
XP\_015152100.1
</td>
<td style="text-align:right;">
3.14
</td>
<td style="text-align:right;">
3.19
</td>
<td style="text-align:right;">
1.24
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015128155.1
</td>
<td style="text-align:right;">
107049005
</td>
<td style="text-align:left;">
LOC107049005
</td>
<td style="text-align:right;">
107049005
</td>
<td style="text-align:left;">
XP\_015128155.1
</td>
<td style="text-align:right;">
7.24
</td>
<td style="text-align:right;">
7.36
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015151573.1
</td>
<td style="text-align:right;">
100858782
</td>
<td style="text-align:left;">
PROCA1
</td>
<td style="text-align:right;">
100858782
</td>
<td style="text-align:left;">
XP\_015151573.1
</td>
<td style="text-align:right;">
3.73
</td>
<td style="text-align:right;">
3.77
</td>
<td style="text-align:right;">
1.65
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015141674.1
</td>
<td style="text-align:right;">
422973
</td>
<td style="text-align:left;">
ANO5
</td>
<td style="text-align:right;">
422973
</td>
<td style="text-align:left;">
XP\_015141674.1
</td>
<td style="text-align:right;">
-5.12
</td>
<td style="text-align:right;">
-144269482.98
</td>
<td style="text-align:right;">
-1.03
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001264456.1
</td>
<td style="text-align:right;">
424379
</td>
<td style="text-align:left;">
REG4
</td>
<td style="text-align:right;">
424379
</td>
<td style="text-align:left;">
NP\_001264456.1
</td>
<td style="text-align:right;">
3.57
</td>
<td style="text-align:right;">
3.58
</td>
<td style="text-align:right;">
2.10
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BH
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1*female.pituitary.inc.d3 -1*female.pituitary.inc.d9
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
glm
</td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FP_i39', frame.plot=F)

![](../figures/pit/01-contrasts-4.png)

    cont <- "FP_i917"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*female.pituitary.inc.d17 1*female.pituitary.inc.d9
    ## Down                                                      57
    ## NotSig                                                 14880
    ## Up                                                         0

    kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)

<table class="kable_wrapper">
<tbody>
<tr>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
row.names
</th>
<th style="text-align:left;">
Name
</th>
<th style="text-align:right;">
geneid
</th>
<th style="text-align:left;">
entrezid
</th>
<th style="text-align:right;">
logFC
</th>
<th style="text-align:right;">
unshrunk.logFC
</th>
<th style="text-align:right;">
logCPM
</th>
<th style="text-align:right;">
PValue
</th>
<th style="text-align:right;">
FDR
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
XP\_015143549.1
</td>
<td style="text-align:right;">
396252
</td>
<td style="text-align:left;">
CDK1
</td>
<td style="text-align:right;">
396252
</td>
<td style="text-align:left;">
XP\_015143549.1
</td>
<td style="text-align:right;">
-4.47
</td>
<td style="text-align:right;">
-4.49
</td>
<td style="text-align:right;">
2.41
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001012888.1
</td>
<td style="text-align:right;">
421226
</td>
<td style="text-align:left;">
BUB1
</td>
<td style="text-align:right;">
421226
</td>
<td style="text-align:left;">
NP\_001012888.1
</td>
<td style="text-align:right;">
-3.79
</td>
<td style="text-align:right;">
-3.82
</td>
<td style="text-align:right;">
1.05
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015135377.1
</td>
<td style="text-align:right;">
417420
</td>
<td style="text-align:left;">
KPNA2
</td>
<td style="text-align:right;">
417420
</td>
<td style="text-align:left;">
XP\_015135377.1
</td>
<td style="text-align:right;">
-4.86
</td>
<td style="text-align:right;">
-4.88
</td>
<td style="text-align:right;">
2.52
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001006274.1
</td>
<td style="text-align:right;">
418882
</td>
<td style="text-align:left;">
CKAP2
</td>
<td style="text-align:right;">
418882
</td>
<td style="text-align:left;">
NP\_001006274.1
</td>
<td style="text-align:right;">
-4.39
</td>
<td style="text-align:right;">
-4.41
</td>
<td style="text-align:right;">
2.60
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_004949789.2
</td>
<td style="text-align:right;">
426884
</td>
<td style="text-align:left;">
RACGAP1
</td>
<td style="text-align:right;">
426884
</td>
<td style="text-align:left;">
XP\_004949789.2
</td>
<td style="text-align:right;">
-3.69
</td>
<td style="text-align:right;">
-3.71
</td>
<td style="text-align:right;">
1.67
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BH
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
-1*female.pituitary.inc.d17 1*female.pituitary.inc.d9
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
glm
</td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FP_i917', frame.plot=F)

![](../figures/pit/01-contrasts-5.png)

    cont <- "FP_i17H"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*female.pituitary.hatch 1*female.pituitary.inc.d17
    ## Down                                                      3
    ## NotSig                                                14934
    ## Up                                                        0

    kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)

<table class="kable_wrapper">
<tbody>
<tr>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
row.names
</th>
<th style="text-align:left;">
Name
</th>
<th style="text-align:right;">
geneid
</th>
<th style="text-align:left;">
entrezid
</th>
<th style="text-align:right;">
logFC
</th>
<th style="text-align:right;">
unshrunk.logFC
</th>
<th style="text-align:right;">
logCPM
</th>
<th style="text-align:right;">
PValue
</th>
<th style="text-align:right;">
FDR
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
XP\_015135755.1
</td>
<td style="text-align:right;">
404271
</td>
<td style="text-align:left;">
ANXA1
</td>
<td style="text-align:right;">
404271
</td>
<td style="text-align:left;">
XP\_015135755.1
</td>
<td style="text-align:right;">
-5.19
</td>
<td style="text-align:right;">
-5.21
</td>
<td style="text-align:right;">
1.87
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001006368.1
</td>
<td style="text-align:right;">
420706
</td>
<td style="text-align:left;">
TGM4
</td>
<td style="text-align:right;">
420706
</td>
<td style="text-align:left;">
NP\_001006368.1
</td>
<td style="text-align:right;">
-5.70
</td>
<td style="text-align:right;">
-5.73
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015133687.1
</td>
<td style="text-align:right;">
396214
</td>
<td style="text-align:left;">
PLP1
</td>
<td style="text-align:right;">
396214
</td>
<td style="text-align:left;">
XP\_015133687.1
</td>
<td style="text-align:right;">
-8.40
</td>
<td style="text-align:right;">
-8.47
</td>
<td style="text-align:right;">
5.46
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_990228.1
</td>
<td style="text-align:right;">
395715
</td>
<td style="text-align:left;">
SERPINB10
</td>
<td style="text-align:right;">
395715
</td>
<td style="text-align:left;">
NP\_990228.1
</td>
<td style="text-align:right;">
-2.89
</td>
<td style="text-align:right;">
-2.90
</td>
<td style="text-align:right;">
1.73
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_001231917.1
</td>
<td style="text-align:right;">
769726
</td>
<td style="text-align:left;">
LOC769726
</td>
<td style="text-align:right;">
769726
</td>
<td style="text-align:left;">
XP\_001231917.1
</td>
<td style="text-align:right;">
-4.85
</td>
<td style="text-align:right;">
-4.88
</td>
<td style="text-align:right;">
3.59
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BH
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
-1*female.pituitary.hatch 1*female.pituitary.inc.d17
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
glm
</td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FP_i17H', frame.plot=F)

![](../figures/pit/01-contrasts-6.png)

    cont <- "FP_H5"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*female.pituitary.hatch -1*female.pituitary.n5
    ## Down                                                 0
    ## NotSig                                           14934
    ## Up                                                   3

    kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)

<table class="kable_wrapper">
<tbody>
<tr>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
row.names
</th>
<th style="text-align:left;">
Name
</th>
<th style="text-align:right;">
geneid
</th>
<th style="text-align:left;">
entrezid
</th>
<th style="text-align:right;">
logFC
</th>
<th style="text-align:right;">
unshrunk.logFC
</th>
<th style="text-align:right;">
logCPM
</th>
<th style="text-align:right;">
PValue
</th>
<th style="text-align:right;">
FDR
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
XP\_015135755.1
</td>
<td style="text-align:right;">
404271
</td>
<td style="text-align:left;">
ANXA1
</td>
<td style="text-align:right;">
404271
</td>
<td style="text-align:left;">
XP\_015135755.1
</td>
<td style="text-align:right;">
5.89
</td>
<td style="text-align:right;">
5.91
</td>
<td style="text-align:right;">
1.87
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001006368.1
</td>
<td style="text-align:right;">
420706
</td>
<td style="text-align:left;">
TGM4
</td>
<td style="text-align:right;">
420706
</td>
<td style="text-align:left;">
NP\_001006368.1
</td>
<td style="text-align:right;">
6.85
</td>
<td style="text-align:right;">
6.93
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001239091.1
</td>
<td style="text-align:right;">
423562
</td>
<td style="text-align:left;">
CDKN3
</td>
<td style="text-align:right;">
423562
</td>
<td style="text-align:left;">
NP\_001239091.1
</td>
<td style="text-align:right;">
1.64
</td>
<td style="text-align:right;">
1.64
</td>
<td style="text-align:right;">
2.69
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_004941667.2
</td>
<td style="text-align:right;">
100859610
</td>
<td style="text-align:left;">
CREB3L1
</td>
<td style="text-align:right;">
100859610
</td>
<td style="text-align:left;">
XP\_004941667.2
</td>
<td style="text-align:right;">
1.34
</td>
<td style="text-align:right;">
1.34
</td>
<td style="text-align:right;">
6.34
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001264404.1
</td>
<td style="text-align:right;">
425027
</td>
<td style="text-align:left;">
SSR3
</td>
<td style="text-align:right;">
425027
</td>
<td style="text-align:left;">
NP\_001264404.1
</td>
<td style="text-align:right;">
1.08
</td>
<td style="text-align:right;">
1.08
</td>
<td style="text-align:right;">
6.85
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BH
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1*female.pituitary.hatch -1*female.pituitary.n5
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
glm
</td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FP_H5', frame.plot=F)

![](../figures/pit/01-contrasts-7.png)

    cont <- "FP_n59"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*female.pituitary.n5 -1*female.pituitary.n9
    ## Down                                              0
    ## NotSig                                        14937
    ## Up                                                0

    kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)

<table class="kable_wrapper">
<tbody>
<tr>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
row.names
</th>
<th style="text-align:left;">
Name
</th>
<th style="text-align:right;">
geneid
</th>
<th style="text-align:left;">
entrezid
</th>
<th style="text-align:right;">
logFC
</th>
<th style="text-align:right;">
unshrunk.logFC
</th>
<th style="text-align:right;">
logCPM
</th>
<th style="text-align:right;">
PValue
</th>
<th style="text-align:right;">
FDR
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
NP\_990049.1
</td>
<td style="text-align:right;">
395465
</td>
<td style="text-align:left;">
CITED4
</td>
<td style="text-align:right;">
395465
</td>
<td style="text-align:left;">
NP\_990049.1
</td>
<td style="text-align:right;">
-1.98
</td>
<td style="text-align:right;">
-1.98
</td>
<td style="text-align:right;">
4.24
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.00
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001278710.1
</td>
<td style="text-align:right;">
417513
</td>
<td style="text-align:left;">
RASA4
</td>
<td style="text-align:right;">
417513
</td>
<td style="text-align:left;">
NP\_001278710.1
</td>
<td style="text-align:right;">
-2.24
</td>
<td style="text-align:right;">
-2.24
</td>
<td style="text-align:right;">
4.38
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.00
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_004947241.2
</td>
<td style="text-align:right;">
101750367
</td>
<td style="text-align:left;">
LOC101750367
</td>
<td style="text-align:right;">
101750367
</td>
<td style="text-align:left;">
XP\_004947241.2
</td>
<td style="text-align:right;">
2.68
</td>
<td style="text-align:right;">
2.74
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_990392.1
</td>
<td style="text-align:right;">
395935
</td>
<td style="text-align:left;">
VTN
</td>
<td style="text-align:right;">
395935
</td>
<td style="text-align:left;">
NP\_990392.1
</td>
<td style="text-align:right;">
-2.34
</td>
<td style="text-align:right;">
-2.34
</td>
<td style="text-align:right;">
4.11
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_415813.4
</td>
<td style="text-align:right;">
417567
</td>
<td style="text-align:left;">
SEBOX
</td>
<td style="text-align:right;">
417567
</td>
<td style="text-align:left;">
XP\_415813.4
</td>
<td style="text-align:right;">
-2.86
</td>
<td style="text-align:right;">
-2.89
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BH
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
1*female.pituitary.n5 -1*female.pituitary.n9
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
glm
</td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FP_n59', frame.plot=F)

![](../figures/pit/01-contrasts-8.png)

    cont <- "FP_n9C"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*female.pituitary.control 1*female.pituitary.n9
    ## Down                                                  34
    ## NotSig                                             14893
    ## Up                                                    10

    kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)

<table class="kable_wrapper">
<tbody>
<tr>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
row.names
</th>
<th style="text-align:left;">
Name
</th>
<th style="text-align:right;">
geneid
</th>
<th style="text-align:left;">
entrezid
</th>
<th style="text-align:right;">
logFC
</th>
<th style="text-align:right;">
unshrunk.logFC
</th>
<th style="text-align:right;">
logCPM
</th>
<th style="text-align:right;">
PValue
</th>
<th style="text-align:right;">
FDR
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
XP\_015128865.1
</td>
<td style="text-align:right;">
101748081
</td>
<td style="text-align:left;">
LOC101748081
</td>
<td style="text-align:right;">
101748081
</td>
<td style="text-align:left;">
XP\_015128865.1
</td>
<td style="text-align:right;">
2.61
</td>
<td style="text-align:right;">
2.62
</td>
<td style="text-align:right;">
2.80
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_003642711.2
</td>
<td style="text-align:right;">
100858704
</td>
<td style="text-align:left;">
ECM1
</td>
<td style="text-align:right;">
100858704
</td>
<td style="text-align:left;">
XP\_003642711.2
</td>
<td style="text-align:right;">
2.73
</td>
<td style="text-align:right;">
2.74
</td>
<td style="text-align:right;">
3.63
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_989881.1
</td>
<td style="text-align:right;">
395233
</td>
<td style="text-align:left;">
DYRK1A
</td>
<td style="text-align:right;">
395233
</td>
<td style="text-align:left;">
NP\_989881.1
</td>
<td style="text-align:right;">
-1.61
</td>
<td style="text-align:right;">
-1.61
</td>
<td style="text-align:right;">
7.88
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_003642170.2
</td>
<td style="text-align:right;">
427675
</td>
<td style="text-align:left;">
RPS15A
</td>
<td style="text-align:right;">
427675
</td>
<td style="text-align:left;">
XP\_003642170.2
</td>
<td style="text-align:right;">
-2.08
</td>
<td style="text-align:right;">
-2.08
</td>
<td style="text-align:right;">
7.98
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015158171.1
</td>
<td style="text-align:right;">
426166
</td>
<td style="text-align:left;">
RASAL3
</td>
<td style="text-align:right;">
426166
</td>
<td style="text-align:left;">
XP\_015158171.1
</td>
<td style="text-align:right;">
1.40
</td>
<td style="text-align:right;">
1.40
</td>
<td style="text-align:right;">
6.08
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
BH
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
-1*female.pituitary.control 1*female.pituitary.n9
</td>
</tr>
</tbody>
</table>
</td>
<td>
<table>
<thead>
<tr>
<th style="text-align:left;">
x
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
glm
</td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FP_n9C', frame.plot=F)

![](../figures/pit/01-contrasts-9.png)

volcano plots
=============

    # from http://www.compbio.dundee.ac.uk/user/pschofield/Projects/teaching_pg/workshops/biocDGE.html#maplots

    lrt <- glmLRT(fit,coef=2)
    topTags(lrt)

    ## Coefficient:  female.pituitary.control 
    ##                row.names         Name    geneid       entrezid      logFC
    ## XP_003642240.1 100857736       EIF2B1 100857736 XP_003642240.1  1.0104389
    ## NP_989881.1       395233       DYRK1A    395233    NP_989881.1  1.7904529
    ## NP_001025818.1    416678        UBE2H    416678 NP_001025818.1 -0.6072189
    ## XP_015158171.1    426166       RASAL3    426166 XP_015158171.1 -1.4583484
    ## XP_015139992.1    396473       MARCKS    396473 XP_015139992.1 -1.3127748
    ## NP_001263232.1    396544         NACA    396544 NP_001263232.1  1.1069333
    ## XP_015143203.1    429422      C14ORF4    429422 XP_015143203.1 -1.2529836
    ## XP_003642711.2 100858704         ECM1 100858704 XP_003642711.2 -2.8349924
    ## XP_015128865.1 101748081 LOC101748081 101748081 XP_015128865.1 -2.5727037
    ## NP_990604.1       396210         NFIA    396210    NP_990604.1 -1.4058106
    ##                  logCPM        LR       PValue          FDR
    ## XP_003642240.1 5.085260 120.91241 3.993812e-28 5.965557e-24
    ## NP_989881.1    7.880155 109.54561 1.232359e-25 9.203874e-22
    ## NP_001025818.1 6.843615  98.07686 4.024535e-23 2.003816e-19
    ## XP_015158171.1 6.077579  97.14707 6.436450e-23 2.403531e-19
    ## XP_015139992.1 5.971584  96.38034 9.480395e-23 2.832173e-19
    ## NP_001263232.1 7.066068  94.66701 2.252626e-22 5.607912e-19
    ## XP_015143203.1 6.642067  91.34341 1.207818e-21 2.577312e-18
    ## XP_003642711.2 3.634435  88.64130 4.733043e-21 8.837183e-18
    ## XP_015128865.1 2.797257  87.15428 1.003777e-20 1.665935e-17
    ## NP_990604.1    6.243470  85.91467 1.878697e-20 2.806209e-17

    tt <- topTags(lrt,n=10000)$table

    ggplot(data=tt) + geom_point(aes(x=logFC,y=-log(FDR),color=logCPM)) +
      scale_colour_gradientn(colours=c("#000000" ,"#FF0000" ))

![](../figures/pit/volcanoplots-1.png)

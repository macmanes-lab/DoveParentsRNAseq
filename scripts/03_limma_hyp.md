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

    knitr::opts_chunk$set(fig.path = '../figures/hyp/',cache=TRUE)

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
      dplyr::filter(grepl('hypothalamus', tissue)) %>%
      droplevels()
    row.names(colData) <- colData$V1

    # print sample sizes
    colData %>% select(sex,treatment, tissue)  %>%  summary()

    ##      sex       treatment           tissue   
    ##  female:95   inc.d9 :23   hypothalamus:189  
    ##  male  :94   control:22                     
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
    ## 14862    75

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

    ## Disp = 0.07817 , BCV = 0.2796

    parentalobject <- estimateGLMTrendedDisp(parentalobject, parentaldesign)
    parentalobject <- estimateGLMTagwiseDisp(parentalobject, parentaldesign)

    #  perform likelihood ratio test and thresholded testing
    fit <- glmFit( parentalobject, parentaldesign, robust=T)
    tr <- glmTreat(fit, lfc = 1)
    topTags(tr)

    ## Coefficient:  male.hypothalamus.n9 
    ##                row.names         Name    geneid       entrezid     logFC
    ## XP_004944394.1 101750188 LOC101750188 101750188 XP_004944394.1 -6.348326
    ## XP_015135862.1    427458       HNRNPK    427458 XP_015135862.1 -7.822022
    ## XP_015135776.1    427254       ZFAND5    427254 XP_015135776.1 -5.949506
    ## XP_001234565.2    771273       MRPS36    771273 XP_001234565.2 -8.187554
    ## XP_015132890.1    770140         CTIF    770140 XP_015132890.1 -5.656195
    ## NP_001264524.1    415770  C11H19ORF40    415770 NP_001264524.1 -6.157800
    ## XP_015128195.1 107049327 LOC107049327 107049327 XP_015128195.1 -5.001708
    ## XP_015128135.1 107049275 LOC107049275 107049275 XP_015128135.1 -3.426814
    ## XP_015151999.1    418412        TRAT1    418412 XP_015151999.1 -3.287743
    ## XP_004937352.1    427168        SREK1    427168 XP_004937352.1 -2.042345
    ##                unshrunk.logFC      logCPM       PValue          FDR
    ## XP_004944394.1  -6.622678e+00  1.34279932 1.623913e-21 2.425639e-17
    ## XP_015135862.1  -1.442695e+08  0.61756011 3.277354e-20 2.447691e-16
    ## XP_015135776.1  -5.998842e+00  3.16280330 2.977696e-19 1.482595e-15
    ## XP_001234565.2  -8.581521e+00  2.66680542 5.524717e-19 2.063067e-15
    ## XP_015132890.1  -6.160934e+00  0.05238573 1.846679e-16 5.516768e-13
    ## NP_001264524.1  -6.315726e+00  1.81575906 5.245826e-13 1.305948e-09
    ## XP_015128195.1  -5.698986e+00 -1.15246178 1.062948e-08 2.268180e-05
    ## XP_015128135.1  -3.532485e+00 -0.21607745 1.158088e-07 2.162295e-04
    ## XP_015151999.1  -3.415411e+00 -0.52286911 9.214142e-06 1.529240e-02
    ## XP_004937352.1  -2.043454e+00  4.68201298 1.117485e-05 1.669187e-02

plotMDS (multidimential scaling)
================================

    plotMDS(parentalobject, cex = 0.5)

![](../figures/hyp/plotMDS-1.png)

    levels(colData$treatment)

    ## [1] "bldg"    "control" "hatch"   "inc.d17" "inc.d3"  "inc.d9"  "lay"    
    ## [8] "n5"      "n9"

    col.treatment <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6")[colData$treatment]

    plotMDS(parentalobject,col=col.treatment, labels = colData$sex)
    legend("bottomright",fill=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6"),legend=levels(colData$treatment))
    title("Hypothalamus Colored by Treatment")

![](../figures/hyp/plotMDS-lables-1.png)

    plotMDS(parentalobject,dim=c(3,4), col=col.treatment, labels = colData$sex)
    legend("bottomright",fill=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6"),legend=levels(colData$treatment))
    title("Hypothalamus Colored by Treatment")

![](../figures/hyp/plotMDS-lables-2.png)

specify contrasts and make MA plots (currently only a subset)
=============================================================

    # view all levels
    levels(colData$group)

    ##  [1] "female.hypothalamus.bldg"    "female.hypothalamus.control"
    ##  [3] "female.hypothalamus.hatch"   "female.hypothalamus.inc.d17"
    ##  [5] "female.hypothalamus.inc.d3"  "female.hypothalamus.inc.d9" 
    ##  [7] "female.hypothalamus.lay"     "female.hypothalamus.n5"     
    ##  [9] "female.hypothalamus.n9"      "male.hypothalamus.bldg"     
    ## [11] "male.hypothalamus.control"   "male.hypothalamus.hatch"    
    ## [13] "male.hypothalamus.inc.d17"   "male.hypothalamus.inc.d3"   
    ## [15] "male.hypothalamus.inc.d9"    "male.hypothalamus.lay"      
    ## [17] "male.hypothalamus.n5"        "male.hypothalamus.n9"

    # subset of conrasts - sex specific comparing hatch to lay
    my.contrasts <- makeContrasts(
                 FH_CB = female.hypothalamus.control - female.hypothalamus.bldg,
                 FH_BL = female.hypothalamus.bldg - female.hypothalamus.lay,
                 FH_Li3 = female.hypothalamus.lay - female.hypothalamus.inc.d3,
                 FH_i39 = female.hypothalamus.inc.d3 - female.hypothalamus.inc.d9,
                 FH_i917 = female.hypothalamus.inc.d9 - female.hypothalamus.inc.d17,
                 FH_i17H = female.hypothalamus.inc.d17 - female.hypothalamus.hatch,
                 FH_H5 = female.hypothalamus.hatch -  female.hypothalamus.n5,
                 FH_n59 = female.hypothalamus.n5 - female.hypothalamus.n9,
                 FH_n9C = female.hypothalamus.n9 - female.hypothalamus.control,
                 
                 MH_CB = male.hypothalamus.control - male.hypothalamus.bldg,
                 MH_BL = male.hypothalamus.bldg - male.hypothalamus.lay,
                 MH_Li3 = male.hypothalamus.lay - male.hypothalamus.inc.d3,
                 MH_i39 = male.hypothalamus.inc.d3 - male.hypothalamus.inc.d9,
                 MH_i917 = male.hypothalamus.inc.d9 - male.hypothalamus.inc.d17,
                 MH_i17H = male.hypothalamus.inc.d17 - male.hypothalamus.hatch,
                 MH_H5 = male.hypothalamus.hatch -  male.hypothalamus.n5,
                 MH_n59 = male.hypothalamus.n5 - male.hypothalamus.n9,
                 MH_n9C = male.hypothalamus.n9 - male.hypothalamus.control,
    levels=parentaldesign)

    # female comparisons
    cont <- "FH_CB"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*female.hypothalamus.bldg 1*female.hypothalamus.control
    ## Down                                                           0
    ## NotSig                                                       263
    ## Up                                                         14674

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
XP\_015140781.1
</td>
<td style="text-align:right;">
422420
</td>
<td style="text-align:left;">
MARCH1
</td>
<td style="text-align:right;">
422420
</td>
<td style="text-align:left;">
XP\_015140781.1
</td>
<td style="text-align:right;">
20.72
</td>
<td style="text-align:right;">
20.74
</td>
<td style="text-align:right;">
0.08
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
XP\_015155805.1
</td>
<td style="text-align:right;">
107055406
</td>
<td style="text-align:left;">
LMBR1L
</td>
<td style="text-align:right;">
107055406
</td>
<td style="text-align:left;">
XP\_015155805.1
</td>
<td style="text-align:right;">
20.49
</td>
<td style="text-align:right;">
20.51
</td>
<td style="text-align:right;">
0.32
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
XP\_015151240.1
</td>
<td style="text-align:right;">
417564
</td>
<td style="text-align:left;">
DPH1
</td>
<td style="text-align:right;">
417564
</td>
<td style="text-align:left;">
XP\_015151240.1
</td>
<td style="text-align:right;">
20.17
</td>
<td style="text-align:right;">
20.18
</td>
<td style="text-align:right;">
1.01
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
NP\_001072968.1
</td>
<td style="text-align:right;">
770922
</td>
<td style="text-align:left;">
DGUOK
</td>
<td style="text-align:right;">
770922
</td>
<td style="text-align:left;">
NP\_001072968.1
</td>
<td style="text-align:right;">
20.09
</td>
<td style="text-align:right;">
20.10
</td>
<td style="text-align:right;">
0.65
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
XP\_015128664.1
</td>
<td style="text-align:right;">
100857346
</td>
<td style="text-align:left;">
CILD27
</td>
<td style="text-align:right;">
100857346
</td>
<td style="text-align:left;">
XP\_015128664.1
</td>
<td style="text-align:right;">
20.06
</td>
<td style="text-align:right;">
20.07
</td>
<td style="text-align:right;">
0.94
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
-1*female.hypothalamus.bldg 1*female.hypothalamus.control
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FH_CB', frame.plot=F)

![](../figures/hyp/01-contrasts-1.png)

    cont <- "FH_BL"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*female.hypothalamus.bldg -1*female.hypothalamus.lay
    ## Down                                                   14742
    ## NotSig                                                   195
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
NP\_989878.1
</td>
<td style="text-align:right;">
395229
</td>
<td style="text-align:left;">
TCF7
</td>
<td style="text-align:right;">
395229
</td>
<td style="text-align:left;">
NP\_989878.1
</td>
<td style="text-align:right;">
-21.00
</td>
<td style="text-align:right;">
-21.03
</td>
<td style="text-align:right;">
-0.59
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
XP\_015144243.1
</td>
<td style="text-align:right;">
428973
</td>
<td style="text-align:left;">
LZTS2
</td>
<td style="text-align:right;">
428973
</td>
<td style="text-align:left;">
XP\_015144243.1
</td>
<td style="text-align:right;">
-20.53
</td>
<td style="text-align:right;">
-20.55
</td>
<td style="text-align:right;">
0.02
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
XP\_015155805.1
</td>
<td style="text-align:right;">
107055406
</td>
<td style="text-align:left;">
LMBR1L
</td>
<td style="text-align:right;">
107055406
</td>
<td style="text-align:left;">
XP\_015155805.1
</td>
<td style="text-align:right;">
-20.41
</td>
<td style="text-align:right;">
-20.43
</td>
<td style="text-align:right;">
0.32
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
XP\_015130583.1
</td>
<td style="text-align:right;">
107051524
</td>
<td style="text-align:left;">
LOC107051524
</td>
<td style="text-align:right;">
107051524
</td>
<td style="text-align:left;">
XP\_015130583.1
</td>
<td style="text-align:right;">
-20.13
</td>
<td style="text-align:right;">
-20.14
</td>
<td style="text-align:right;">
0.90
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
XP\_015140781.1
</td>
<td style="text-align:right;">
422420
</td>
<td style="text-align:left;">
MARCH1
</td>
<td style="text-align:right;">
422420
</td>
<td style="text-align:left;">
XP\_015140781.1
</td>
<td style="text-align:right;">
-20.06
</td>
<td style="text-align:right;">
-20.08
</td>
<td style="text-align:right;">
0.08
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
1*female.hypothalamus.bldg -1*female.hypothalamus.lay
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FH_BL', frame.plot=F)

![](../figures/hyp/01-contrasts-2.png)

    cont <- "FH_Li3"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*female.hypothalamus.inc.d3 1*female.hypothalamus.lay
    ## Down                                                         0
    ## NotSig                                                   14937
    ## Up                                                           0

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
NP\_001231522.1
</td>
<td style="text-align:right;">
420588
</td>
<td style="text-align:left;">
SCIN
</td>
<td style="text-align:right;">
420588
</td>
<td style="text-align:left;">
NP\_001231522.1
</td>
<td style="text-align:right;">
1.16
</td>
<td style="text-align:right;">
1.16
</td>
<td style="text-align:right;">
3.65
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.25
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001039297.1
</td>
<td style="text-align:right;">
417465
</td>
<td style="text-align:left;">
CCL5
</td>
<td style="text-align:right;">
417465
</td>
<td style="text-align:left;">
NP\_001039297.1
</td>
<td style="text-align:right;">
3.74
</td>
<td style="text-align:right;">
4.15
</td>
<td style="text-align:right;">
-1.87
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.98
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001186559.1
</td>
<td style="text-align:right;">
417852
</td>
<td style="text-align:left;">
KCNMB4
</td>
<td style="text-align:right;">
417852
</td>
<td style="text-align:left;">
NP\_001186559.1
</td>
<td style="text-align:right;">
2.18
</td>
<td style="text-align:right;">
2.56
</td>
<td style="text-align:right;">
-1.99
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.98
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015149500.1
</td>
<td style="text-align:right;">
396384
</td>
<td style="text-align:left;">
IRF1
</td>
<td style="text-align:right;">
396384
</td>
<td style="text-align:left;">
XP\_015149500.1
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:right;">
1.04
</td>
<td style="text-align:right;">
2.54
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.98
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001231990.1
</td>
<td style="text-align:right;">
430189
</td>
<td style="text-align:left;">
HLA-DRA
</td>
<td style="text-align:right;">
430189
</td>
<td style="text-align:left;">
NP\_001231990.1
</td>
<td style="text-align:right;">
1.48
</td>
<td style="text-align:right;">
1.48
</td>
<td style="text-align:right;">
3.13
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.98
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
-1*female.hypothalamus.inc.d3 1*female.hypothalamus.lay
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FH_Li3', frame.plot=F)

![](../figures/hyp/01-contrasts-3.png)

    cont <- "FH_i39"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*female.hypothalamus.inc.d3 -1*female.hypothalamus.inc.d9
    ## Down                                                            0
    ## NotSig                                                      14937
    ## Up                                                              0

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
XP\_419585.2
</td>
<td style="text-align:right;">
421544
</td>
<td style="text-align:left;">
CAPN9
</td>
<td style="text-align:right;">
421544
</td>
<td style="text-align:left;">
XP\_419585.2
</td>
<td style="text-align:right;">
1.85
</td>
<td style="text-align:right;">
2.02
</td>
<td style="text-align:right;">
-1.86
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_004934606.1
</td>
<td style="text-align:right;">
427984
</td>
<td style="text-align:left;">
KCNJ15
</td>
<td style="text-align:right;">
427984
</td>
<td style="text-align:left;">
XP\_004934606.1
</td>
<td style="text-align:right;">
1.72
</td>
<td style="text-align:right;">
1.87
</td>
<td style="text-align:right;">
-1.69
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001186559.1
</td>
<td style="text-align:right;">
417852
</td>
<td style="text-align:left;">
KCNMB4
</td>
<td style="text-align:right;">
417852
</td>
<td style="text-align:left;">
NP\_001186559.1
</td>
<td style="text-align:right;">
-1.93
</td>
<td style="text-align:right;">
-2.28
</td>
<td style="text-align:right;">
-1.99
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_415429.5
</td>
<td style="text-align:right;">
395549
</td>
<td style="text-align:left;">
DBH
</td>
<td style="text-align:right;">
395549
</td>
<td style="text-align:left;">
XP\_415429.5
</td>
<td style="text-align:right;">
2.22
</td>
<td style="text-align:right;">
2.61
</td>
<td style="text-align:right;">
-1.81
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
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
4.90
</td>
<td style="text-align:right;">
5.46
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
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
1*female.hypothalamus.inc.d3 -1*female.hypothalamus.inc.d9
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FH_i39', frame.plot=F)

![](../figures/hyp/01-contrasts-4.png)

    cont <- "FH_i917"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*female.hypothalamus.inc.d17 1*female.hypothalamus.inc.d9
    ## Down                                                             0
    ## NotSig                                                       14937
    ## Up                                                               0

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
XP\_015130802.1
</td>
<td style="text-align:right;">
416928
</td>
<td style="text-align:left;">
IGLL1
</td>
<td style="text-align:right;">
416928
</td>
<td style="text-align:left;">
XP\_015130802.1
</td>
<td style="text-align:right;">
-4.19
</td>
<td style="text-align:right;">
-4.20
</td>
<td style="text-align:right;">
2.97
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.37
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_004938467.2
</td>
<td style="text-align:right;">
418651
</td>
<td style="text-align:left;">
SHROOM2
</td>
<td style="text-align:right;">
418651
</td>
<td style="text-align:left;">
XP\_004938467.2
</td>
<td style="text-align:right;">
-0.90
</td>
<td style="text-align:right;">
-0.90
</td>
<td style="text-align:right;">
7.82
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.43
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_003642428.2
</td>
<td style="text-align:right;">
417471
</td>
<td style="text-align:left;">
HSF5
</td>
<td style="text-align:right;">
417471
</td>
<td style="text-align:left;">
XP\_003642428.2
</td>
<td style="text-align:right;">
3.84
</td>
<td style="text-align:right;">
144269481.23
</td>
<td style="text-align:right;">
-1.80
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_989594.1
</td>
<td style="text-align:right;">
374117
</td>
<td style="text-align:left;">
IGJ
</td>
<td style="text-align:right;">
374117
</td>
<td style="text-align:left;">
NP\_989594.1
</td>
<td style="text-align:right;">
-3.10
</td>
<td style="text-align:right;">
-3.11
</td>
<td style="text-align:right;">
1.95
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015128794.1
</td>
<td style="text-align:right;">
107049579
</td>
<td style="text-align:left;">
LOC107049579
</td>
<td style="text-align:right;">
107049579
</td>
<td style="text-align:left;">
XP\_015128794.1
</td>
<td style="text-align:right;">
-4.32
</td>
<td style="text-align:right;">
-144269481.74
</td>
<td style="text-align:right;">
-1.95
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.76
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
-1*female.hypothalamus.inc.d17 1*female.hypothalamus.inc.d9
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FH_i917', frame.plot=F)

![](../figures/hyp/01-contrasts-5.png)

    cont <- "FH_i17H"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*female.hypothalamus.hatch 1*female.hypothalamus.inc.d17
    ## Down                                                            0
    ## NotSig                                                      14937
    ## Up                                                              0

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
NP\_990411.1
</td>
<td style="text-align:right;">
395963
</td>
<td style="text-align:left;">
CAPN2
</td>
<td style="text-align:right;">
395963
</td>
<td style="text-align:left;">
NP\_990411.1
</td>
<td style="text-align:right;">
-2.17
</td>
<td style="text-align:right;">
-2.17
</td>
<td style="text-align:right;">
7.82
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
XP\_015137344.1
</td>
<td style="text-align:right;">
420606
</td>
<td style="text-align:left;">
ABCB5
</td>
<td style="text-align:right;">
420606
</td>
<td style="text-align:left;">
XP\_015137344.1
</td>
<td style="text-align:right;">
-4.29
</td>
<td style="text-align:right;">
-4.43
</td>
<td style="text-align:right;">
-1.45
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001026427.2
</td>
<td style="text-align:right;">
424189
</td>
<td style="text-align:left;">
PSMD14
</td>
<td style="text-align:right;">
424189
</td>
<td style="text-align:left;">
NP\_001026427.2
</td>
<td style="text-align:right;">
-1.05
</td>
<td style="text-align:right;">
-1.05
</td>
<td style="text-align:right;">
4.41
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.26
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001072943.1
</td>
<td style="text-align:right;">
404299
</td>
<td style="text-align:left;">
HIST1H2A4
</td>
<td style="text-align:right;">
404299
</td>
<td style="text-align:left;">
NP\_001072943.1
</td>
<td style="text-align:right;">
-1.19
</td>
<td style="text-align:right;">
-1.19
</td>
<td style="text-align:right;">
2.86
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.40
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001075178.3
</td>
<td style="text-align:right;">
771173
</td>
<td style="text-align:left;">
TLR1B
</td>
<td style="text-align:right;">
771173
</td>
<td style="text-align:left;">
NP\_001075178.3
</td>
<td style="text-align:right;">
-1.81
</td>
<td style="text-align:right;">
-1.86
</td>
<td style="text-align:right;">
-1.13
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.44
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
-1*female.hypothalamus.hatch 1*female.hypothalamus.inc.d17
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FH_i17H', frame.plot=F)

![](../figures/hyp/01-contrasts-6.png)

    cont <- "FH_H5"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*female.hypothalamus.hatch -1*female.hypothalamus.n5
    ## Down                                                       0
    ## NotSig                                                 14936
    ## Up                                                         1

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
NP\_990411.1
</td>
<td style="text-align:right;">
395963
</td>
<td style="text-align:left;">
CAPN2
</td>
<td style="text-align:right;">
395963
</td>
<td style="text-align:left;">
NP\_990411.1
</td>
<td style="text-align:right;">
2.68
</td>
<td style="text-align:right;">
2.69
</td>
<td style="text-align:right;">
7.82
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
XP\_015137344.1
</td>
<td style="text-align:right;">
420606
</td>
<td style="text-align:left;">
ABCB5
</td>
<td style="text-align:right;">
420606
</td>
<td style="text-align:left;">
XP\_015137344.1
</td>
<td style="text-align:right;">
5.73
</td>
<td style="text-align:right;">
6.16
</td>
<td style="text-align:right;">
-1.45
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
2.53
</td>
<td style="text-align:right;">
2.53
</td>
<td style="text-align:right;">
2.62
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
NP\_990229.1
</td>
<td style="text-align:right;">
395719
</td>
<td style="text-align:left;">
TK1
</td>
<td style="text-align:right;">
395719
</td>
<td style="text-align:left;">
NP\_990229.1
</td>
<td style="text-align:right;">
1.01
</td>
<td style="text-align:right;">
1.02
</td>
<td style="text-align:right;">
1.68
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015155605.1
</td>
<td style="text-align:right;">
100858365
</td>
<td style="text-align:left;">
HAPLN4
</td>
<td style="text-align:right;">
100858365
</td>
<td style="text-align:left;">
XP\_015155605.1
</td>
<td style="text-align:right;">
1.69
</td>
<td style="text-align:right;">
1.69
</td>
<td style="text-align:right;">
3.55
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.09
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
1*female.hypothalamus.hatch -1*female.hypothalamus.n5
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FH_H5', frame.plot=F)

![](../figures/hyp/01-contrasts-7.png)

    cont <- "FH_n59"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*female.hypothalamus.n5 -1*female.hypothalamus.n9
    ## Down                                                    0
    ## NotSig                                              14937
    ## Up                                                      0

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
XP\_416510.2
</td>
<td style="text-align:right;">
418287
</td>
<td style="text-align:left;">
LAG3
</td>
<td style="text-align:right;">
418287
</td>
<td style="text-align:left;">
XP\_416510.2
</td>
<td style="text-align:right;">
3.35
</td>
<td style="text-align:right;">
4.50
</td>
<td style="text-align:right;">
-2.17
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015146387.1
</td>
<td style="text-align:right;">
424552
</td>
<td style="text-align:left;">
NEXN
</td>
<td style="text-align:right;">
424552
</td>
<td style="text-align:left;">
XP\_015146387.1
</td>
<td style="text-align:right;">
3.93
</td>
<td style="text-align:right;">
144269481.33
</td>
<td style="text-align:right;">
-2.02
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_990847.1
</td>
<td style="text-align:right;">
396522
</td>
<td style="text-align:left;">
CNN1
</td>
<td style="text-align:right;">
396522
</td>
<td style="text-align:left;">
NP\_990847.1
</td>
<td style="text-align:right;">
-1.15
</td>
<td style="text-align:right;">
-1.18
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015134716.1
</td>
<td style="text-align:right;">
415686
</td>
<td style="text-align:left;">
MLKL
</td>
<td style="text-align:right;">
415686
</td>
<td style="text-align:left;">
XP\_015134716.1
</td>
<td style="text-align:right;">
-2.68
</td>
<td style="text-align:right;">
-3.34
</td>
<td style="text-align:right;">
-2.08
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_419781.2
</td>
<td style="text-align:right;">
421750
</td>
<td style="text-align:left;">
TUBE1
</td>
<td style="text-align:right;">
421750
</td>
<td style="text-align:left;">
XP\_419781.2
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
2.44
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
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
1*female.hypothalamus.n5 -1*female.hypothalamus.n9
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FH_n59', frame.plot=F)

![](../figures/hyp/01-contrasts-8.png)

    cont <- "FH_n9C"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*female.hypothalamus.control 1*female.hypothalamus.n9
    ## Down                                                        70
    ## NotSig                                                   14858
    ## Up                                                           9

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
XP\_004950299.2
</td>
<td style="text-align:right;">
107049137
</td>
<td style="text-align:left;">
LOC107049137
</td>
<td style="text-align:right;">
107049137
</td>
<td style="text-align:left;">
XP\_004950299.2
</td>
<td style="text-align:right;">
-4.16
</td>
<td style="text-align:right;">
-4.16
</td>
<td style="text-align:right;">
3.09
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
XP\_015129796.1
</td>
<td style="text-align:right;">
100859555
</td>
<td style="text-align:left;">
CLAPS2
</td>
<td style="text-align:right;">
100859555
</td>
<td style="text-align:left;">
XP\_015129796.1
</td>
<td style="text-align:right;">
-2.38
</td>
<td style="text-align:right;">
-2.38
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
XP\_015155850.1
</td>
<td style="text-align:right;">
101749643
</td>
<td style="text-align:left;">
LOC101749643
</td>
<td style="text-align:right;">
101749643
</td>
<td style="text-align:left;">
XP\_015155850.1
</td>
<td style="text-align:right;">
-3.29
</td>
<td style="text-align:right;">
-3.29
</td>
<td style="text-align:right;">
1.53
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
XP\_015147109.1
</td>
<td style="text-align:right;">
424884
</td>
<td style="text-align:left;">
PLSCR5
</td>
<td style="text-align:right;">
424884
</td>
<td style="text-align:left;">
XP\_015147109.1
</td>
<td style="text-align:right;">
-2.46
</td>
<td style="text-align:right;">
-2.46
</td>
<td style="text-align:right;">
2.48
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
NP\_006917.1
</td>
<td style="text-align:right;">
807639
</td>
<td style="text-align:left;">
COX1
</td>
<td style="text-align:right;">
807639
</td>
<td style="text-align:left;">
NP\_006917.1
</td>
<td style="text-align:right;">
-1.57
</td>
<td style="text-align:right;">
-1.57
</td>
<td style="text-align:right;">
16.11
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
-1*female.hypothalamus.control 1*female.hypothalamus.n9
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='FH_n9C', frame.plot=F)

![](../figures/hyp/01-contrasts-9.png)

    ## male comparisons

    cont <- "MH_CB"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*male.hypothalamus.bldg 1*male.hypothalamus.control
    ## Down                                                      11
    ## NotSig                                                 14820
    ## Up                                                       106

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
NP\_006917.1
</td>
<td style="text-align:right;">
807639
</td>
<td style="text-align:left;">
COX1
</td>
<td style="text-align:right;">
807639
</td>
<td style="text-align:left;">
NP\_006917.1
</td>
<td style="text-align:right;">
1.96
</td>
<td style="text-align:right;">
1.96
</td>
<td style="text-align:right;">
16.11
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
NP\_001264272.1
</td>
<td style="text-align:right;">
423726
</td>
<td style="text-align:left;">
RPS24
</td>
<td style="text-align:right;">
423726
</td>
<td style="text-align:left;">
NP\_001264272.1
</td>
<td style="text-align:right;">
2.13
</td>
<td style="text-align:right;">
2.13
</td>
<td style="text-align:right;">
7.14
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
XP\_003640629.1
</td>
<td style="text-align:right;">
771615
</td>
<td style="text-align:left;">
MRP63
</td>
<td style="text-align:right;">
771615
</td>
<td style="text-align:left;">
XP\_003640629.1
</td>
<td style="text-align:right;">
2.24
</td>
<td style="text-align:right;">
2.24
</td>
<td style="text-align:right;">
4.10
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
NP\_001007968.1
</td>
<td style="text-align:right;">
425416
</td>
<td style="text-align:left;">
RPL30
</td>
<td style="text-align:right;">
425416
</td>
<td style="text-align:left;">
NP\_001007968.1
</td>
<td style="text-align:right;">
2.36
</td>
<td style="text-align:right;">
2.36
</td>
<td style="text-align:right;">
6.58
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
NP\_001279014.1
</td>
<td style="text-align:right;">
768911
</td>
<td style="text-align:left;">
RSL24D1
</td>
<td style="text-align:right;">
768911
</td>
<td style="text-align:left;">
NP\_001279014.1
</td>
<td style="text-align:right;">
1.62
</td>
<td style="text-align:right;">
1.62
</td>
<td style="text-align:right;">
5.17
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
-1*male.hypothalamus.bldg 1*male.hypothalamus.control
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='MH_CB', frame.plot=F)

![](../figures/hyp/01-contrasts-10.png)

    cont <- "MH_BL"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*male.hypothalamus.bldg -1*male.hypothalamus.lay
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
NP\_001001301.1
</td>
<td style="text-align:right;">
408026
</td>
<td style="text-align:left;">
TPH2
</td>
<td style="text-align:right;">
408026
</td>
<td style="text-align:left;">
NP\_001001301.1
</td>
<td style="text-align:right;">
3.20
</td>
<td style="text-align:right;">
3.22
</td>
<td style="text-align:right;">
0.52
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
NP\_990136.1
</td>
<td style="text-align:right;">
395592
</td>
<td style="text-align:left;">
TH
</td>
<td style="text-align:right;">
395592
</td>
<td style="text-align:left;">
NP\_990136.1
</td>
<td style="text-align:right;">
2.08
</td>
<td style="text-align:right;">
2.08
</td>
<td style="text-align:right;">
2.85
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
XP\_015154090.1
</td>
<td style="text-align:right;">
107055138
</td>
<td style="text-align:left;">
HAX1
</td>
<td style="text-align:right;">
107055138
</td>
<td style="text-align:left;">
XP\_015154090.1
</td>
<td style="text-align:right;">
-1.44
</td>
<td style="text-align:right;">
-1.44
</td>
<td style="text-align:right;">
4.01
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_004935201.1
</td>
<td style="text-align:right;">
420947
</td>
<td style="text-align:left;">
DDC
</td>
<td style="text-align:right;">
420947
</td>
<td style="text-align:left;">
XP\_004935201.1
</td>
<td style="text-align:right;">
1.23
</td>
<td style="text-align:right;">
1.23
</td>
<td style="text-align:right;">
3.54
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_004941024.1
</td>
<td style="text-align:right;">
101749269
</td>
<td style="text-align:left;">
LOC101749269
</td>
<td style="text-align:right;">
101749269
</td>
<td style="text-align:left;">
XP\_004941024.1
</td>
<td style="text-align:right;">
3.23
</td>
<td style="text-align:right;">
3.39
</td>
<td style="text-align:right;">
-1.25
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.21
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
1*male.hypothalamus.bldg -1*male.hypothalamus.lay
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='MH_BL', frame.plot=F)

![](../figures/hyp/01-contrasts-11.png)

    cont <- "MH_Li3"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*male.hypothalamus.inc.d3 1*male.hypothalamus.lay
    ## Down                                                     0
    ## NotSig                                               14937
    ## Up                                                       0

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
XP\_015154090.1
</td>
<td style="text-align:right;">
107055138
</td>
<td style="text-align:left;">
HAX1
</td>
<td style="text-align:right;">
107055138
</td>
<td style="text-align:left;">
XP\_015154090.1
</td>
<td style="text-align:right;">
1.54
</td>
<td style="text-align:right;">
1.54
</td>
<td style="text-align:right;">
4.01
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
NP\_001289086.1
</td>
<td style="text-align:right;">
426880
</td>
<td style="text-align:left;">
COX14
</td>
<td style="text-align:right;">
426880
</td>
<td style="text-align:left;">
NP\_001289086.1
</td>
<td style="text-align:right;">
1.44
</td>
<td style="text-align:right;">
1.44
</td>
<td style="text-align:right;">
2.81
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.34
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_424461.2
</td>
<td style="text-align:right;">
426853
</td>
<td style="text-align:left;">
ALPK2
</td>
<td style="text-align:right;">
426853
</td>
<td style="text-align:left;">
XP\_424461.2
</td>
<td style="text-align:right;">
-1.34
</td>
<td style="text-align:right;">
-1.36
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.60
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_003640374.1
</td>
<td style="text-align:right;">
429671
</td>
<td style="text-align:left;">
SCO2
</td>
<td style="text-align:right;">
429671
</td>
<td style="text-align:left;">
XP\_003640374.1
</td>
<td style="text-align:right;">
1.09
</td>
<td style="text-align:right;">
1.10
</td>
<td style="text-align:right;">
2.28
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015145655.1
</td>
<td style="text-align:right;">
107053895
</td>
<td style="text-align:left;">
RPRM
</td>
<td style="text-align:right;">
107053895
</td>
<td style="text-align:left;">
XP\_015145655.1
</td>
<td style="text-align:right;">
1.30
</td>
<td style="text-align:right;">
1.30
</td>
<td style="text-align:right;">
2.52
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.00
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
-1*male.hypothalamus.inc.d3 1*male.hypothalamus.lay
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='MH_Li3', frame.plot=F)

![](../figures/hyp/01-contrasts-12.png)

    cont <- "MH_i39"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*male.hypothalamus.inc.d3 -1*male.hypothalamus.inc.d9
    ## Down                                                        1
    ## NotSig                                                  14936
    ## Up                                                          0

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
XP\_001234565.2
</td>
<td style="text-align:right;">
771273
</td>
<td style="text-align:left;">
MRPS36
</td>
<td style="text-align:right;">
771273
</td>
<td style="text-align:left;">
XP\_001234565.2
</td>
<td style="text-align:right;">
-6.29
</td>
<td style="text-align:right;">
-7.59
</td>
<td style="text-align:right;">
2.67
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
NP\_001264524.1
</td>
<td style="text-align:right;">
415770
</td>
<td style="text-align:left;">
C11H19ORF40
</td>
<td style="text-align:right;">
415770
</td>
<td style="text-align:left;">
NP\_001264524.1
</td>
<td style="text-align:right;">
-3.65
</td>
<td style="text-align:right;">
-3.95
</td>
<td style="text-align:right;">
1.82
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_004944394.1
</td>
<td style="text-align:right;">
101750188
</td>
<td style="text-align:left;">
LOC101750188
</td>
<td style="text-align:right;">
101750188
</td>
<td style="text-align:left;">
XP\_004944394.1
</td>
<td style="text-align:right;">
-2.86
</td>
<td style="text-align:right;">
-3.13
</td>
<td style="text-align:right;">
1.34
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015135776.1
</td>
<td style="text-align:right;">
427254
</td>
<td style="text-align:left;">
ZFAND5
</td>
<td style="text-align:right;">
427254
</td>
<td style="text-align:left;">
XP\_015135776.1
</td>
<td style="text-align:right;">
-2.44
</td>
<td style="text-align:right;">
-2.49
</td>
<td style="text-align:right;">
3.16
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015130802.1
</td>
<td style="text-align:right;">
416928
</td>
<td style="text-align:left;">
IGLL1
</td>
<td style="text-align:right;">
416928
</td>
<td style="text-align:left;">
XP\_015130802.1
</td>
<td style="text-align:right;">
4.16
</td>
<td style="text-align:right;">
4.18
</td>
<td style="text-align:right;">
2.97
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.19
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
1*male.hypothalamus.inc.d3 -1*male.hypothalamus.inc.d9
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='MH_i39', frame.plot=F)

![](../figures/hyp/01-contrasts-13.png)

    cont <- "MH_i917"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*male.hypothalamus.inc.d17 1*male.hypothalamus.inc.d9
    ## Down                                                         0
    ## NotSig                                                   14935
    ## Up                                                           2

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
XP\_001234565.2
</td>
<td style="text-align:right;">
771273
</td>
<td style="text-align:left;">
MRPS36
</td>
<td style="text-align:right;">
771273
</td>
<td style="text-align:left;">
XP\_001234565.2
</td>
<td style="text-align:right;">
5.81
</td>
<td style="text-align:right;">
6.62
</td>
<td style="text-align:right;">
2.67
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
NP\_001264524.1
</td>
<td style="text-align:right;">
415770
</td>
<td style="text-align:left;">
C11H19ORF40
</td>
<td style="text-align:right;">
415770
</td>
<td style="text-align:left;">
NP\_001264524.1
</td>
<td style="text-align:right;">
4.64
</td>
<td style="text-align:right;">
5.44
</td>
<td style="text-align:right;">
1.82
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
XP\_004944394.1
</td>
<td style="text-align:right;">
101750188
</td>
<td style="text-align:left;">
LOC101750188
</td>
<td style="text-align:right;">
101750188
</td>
<td style="text-align:left;">
XP\_004944394.1
</td>
<td style="text-align:right;">
3.62
</td>
<td style="text-align:right;">
4.18
</td>
<td style="text-align:right;">
1.34
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
XP\_015128714.1
</td>
<td style="text-align:right;">
395532
</td>
<td style="text-align:left;">
COL1A1
</td>
<td style="text-align:right;">
395532
</td>
<td style="text-align:left;">
XP\_015128714.1
</td>
<td style="text-align:right;">
-1.28
</td>
<td style="text-align:right;">
-1.29
</td>
<td style="text-align:right;">
2.29
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015132890.1
</td>
<td style="text-align:right;">
770140
</td>
<td style="text-align:left;">
CTIF
</td>
<td style="text-align:right;">
770140
</td>
<td style="text-align:left;">
XP\_015132890.1
</td>
<td style="text-align:right;">
3.47
</td>
<td style="text-align:right;">
144269480.80
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.11
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
-1*male.hypothalamus.inc.d17 1*male.hypothalamus.inc.d9
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='MH_i917', frame.plot=F)

![](../figures/hyp/01-contrasts-14.png)

    cont <- "MH_i17H"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*male.hypothalamus.hatch 1*male.hypothalamus.inc.d17
    ## Down                                                        0
    ## NotSig                                                  14937
    ## Up                                                          0

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
XP\_004942140.1
</td>
<td style="text-align:right;">
423744
</td>
<td style="text-align:left;">
MYOZ1
</td>
<td style="text-align:right;">
423744
</td>
<td style="text-align:left;">
XP\_004942140.1
</td>
<td style="text-align:right;">
-2.82
</td>
<td style="text-align:right;">
-2.98
</td>
<td style="text-align:right;">
-1.38
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
NP\_001090995.1
</td>
<td style="text-align:right;">
417506
</td>
<td style="text-align:left;">
MYL10
</td>
<td style="text-align:right;">
417506
</td>
<td style="text-align:left;">
NP\_001090995.1
</td>
<td style="text-align:right;">
-2.45
</td>
<td style="text-align:right;">
-2.49
</td>
<td style="text-align:right;">
-0.69
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
NP\_990781.1
</td>
<td style="text-align:right;">
396434
</td>
<td style="text-align:left;">
TNNC2
</td>
<td style="text-align:right;">
396434
</td>
<td style="text-align:left;">
NP\_990781.1
</td>
<td style="text-align:right;">
-5.50
</td>
<td style="text-align:right;">
-5.75
</td>
<td style="text-align:right;">
-1.37
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
XP\_015157102.1
</td>
<td style="text-align:right;">
101748402
</td>
<td style="text-align:left;">
LOC101748402
</td>
<td style="text-align:right;">
101748402
</td>
<td style="text-align:left;">
XP\_015157102.1
</td>
<td style="text-align:right;">
-2.43
</td>
<td style="text-align:right;">
-2.45
</td>
<td style="text-align:right;">
0.14
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
XP\_015143153.1
</td>
<td style="text-align:right;">
418099
</td>
<td style="text-align:left;">
MYBPC1
</td>
<td style="text-align:right;">
418099
</td>
<td style="text-align:left;">
XP\_015143153.1
</td>
<td style="text-align:right;">
-3.46
</td>
<td style="text-align:right;">
-3.55
</td>
<td style="text-align:right;">
-1.41
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.04
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
-1*male.hypothalamus.hatch 1*male.hypothalamus.inc.d17
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='MH_i17H', frame.plot=F)

![](../figures/hyp/01-contrasts-15.png)

    cont <- "MH_H5"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*male.hypothalamus.hatch -1*male.hypothalamus.n5
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
XP\_015143153.1
</td>
<td style="text-align:right;">
418099
</td>
<td style="text-align:left;">
MYBPC1
</td>
<td style="text-align:right;">
418099
</td>
<td style="text-align:left;">
XP\_015143153.1
</td>
<td style="text-align:right;">
4.23
</td>
<td style="text-align:right;">
4.38
</td>
<td style="text-align:right;">
-1.41
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
XP\_015145352.1
</td>
<td style="text-align:right;">
374027
</td>
<td style="text-align:left;">
NEB
</td>
<td style="text-align:right;">
374027
</td>
<td style="text-align:left;">
XP\_015145352.1
</td>
<td style="text-align:right;">
1.61
</td>
<td style="text-align:right;">
1.61
</td>
<td style="text-align:right;">
1.56
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015141730.1
</td>
<td style="text-align:right;">
422974
</td>
<td style="text-align:left;">
SLC6A5
</td>
<td style="text-align:right;">
422974
</td>
<td style="text-align:left;">
XP\_015141730.1
</td>
<td style="text-align:right;">
-3.03
</td>
<td style="text-align:right;">
-3.14
</td>
<td style="text-align:right;">
-1.29
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_990748.1
</td>
<td style="text-align:right;">
396386
</td>
<td style="text-align:left;">
TNNI2
</td>
<td style="text-align:right;">
396386
</td>
<td style="text-align:left;">
NP\_990748.1
</td>
<td style="text-align:right;">
3.42
</td>
<td style="text-align:right;">
3.47
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001001301.1
</td>
<td style="text-align:right;">
408026
</td>
<td style="text-align:left;">
TPH2
</td>
<td style="text-align:right;">
408026
</td>
<td style="text-align:left;">
NP\_001001301.1
</td>
<td style="text-align:right;">
-2.39
</td>
<td style="text-align:right;">
-2.41
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.09
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
1*male.hypothalamus.hatch -1*male.hypothalamus.n5
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='MH_H5', frame.plot=F)

![](../figures/hyp/01-contrasts-16.png)

    cont <- "MH_n59"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*male.hypothalamus.n5 -1*male.hypothalamus.n9
    ## Down                                                0
    ## NotSig                                          14937
    ## Up                                                  0

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
NP\_001001296.1
</td>
<td style="text-align:right;">
403120
</td>
<td style="text-align:left;">
IFI6
</td>
<td style="text-align:right;">
403120
</td>
<td style="text-align:left;">
NP\_001001296.1
</td>
<td style="text-align:right;">
3.23
</td>
<td style="text-align:right;">
3.28
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
NP\_001001301.1
</td>
<td style="text-align:right;">
408026
</td>
<td style="text-align:left;">
TPH2
</td>
<td style="text-align:right;">
408026
</td>
<td style="text-align:left;">
NP\_001001301.1
</td>
<td style="text-align:right;">
2.57
</td>
<td style="text-align:right;">
2.58
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015142180.1
</td>
<td style="text-align:right;">
107053409
</td>
<td style="text-align:left;">
LOC107053409
</td>
<td style="text-align:right;">
107053409
</td>
<td style="text-align:left;">
XP\_015142180.1
</td>
<td style="text-align:right;">
-2.73
</td>
<td style="text-align:right;">
-3.39
</td>
<td style="text-align:right;">
-1.97
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_015141730.1
</td>
<td style="text-align:right;">
422974
</td>
<td style="text-align:left;">
SLC6A5
</td>
<td style="text-align:right;">
422974
</td>
<td style="text-align:left;">
XP\_015141730.1
</td>
<td style="text-align:right;">
2.32
</td>
<td style="text-align:right;">
2.38
</td>
<td style="text-align:right;">
-1.29
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
XP\_004936174.2
</td>
<td style="text-align:right;">
426317
</td>
<td style="text-align:left;">
TMEM156
</td>
<td style="text-align:right;">
426317
</td>
<td style="text-align:left;">
XP\_004936174.2
</td>
<td style="text-align:right;">
3.10
</td>
<td style="text-align:right;">
144269480.41
</td>
<td style="text-align:right;">
-2.18
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1.00
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
1*male.hypothalamus.n5 -1*male.hypothalamus.n9
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='MH_n59', frame.plot=F)

![](../figures/hyp/01-contrasts-17.png)

    cont <- "MH_n9C"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        -1*male.hypothalamus.control 1*male.hypothalamus.n9
    ## Down                                                   199
    ## NotSig                                               14717
    ## Up                                                      21

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
NP\_006917.1
</td>
<td style="text-align:right;">
807639
</td>
<td style="text-align:left;">
COX1
</td>
<td style="text-align:right;">
807639
</td>
<td style="text-align:left;">
NP\_006917.1
</td>
<td style="text-align:right;">
-2.03
</td>
<td style="text-align:right;">
-2.03
</td>
<td style="text-align:right;">
16.11
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
XP\_003640629.1
</td>
<td style="text-align:right;">
771615
</td>
<td style="text-align:left;">
MRP63
</td>
<td style="text-align:right;">
771615
</td>
<td style="text-align:left;">
XP\_003640629.1
</td>
<td style="text-align:right;">
-2.37
</td>
<td style="text-align:right;">
-2.37
</td>
<td style="text-align:right;">
4.10
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
NP\_001264272.1
</td>
<td style="text-align:right;">
423726
</td>
<td style="text-align:left;">
RPS24
</td>
<td style="text-align:right;">
423726
</td>
<td style="text-align:left;">
NP\_001264272.1
</td>
<td style="text-align:right;">
-2.17
</td>
<td style="text-align:right;">
-2.17
</td>
<td style="text-align:right;">
7.14
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
NP\_001007968.1
</td>
<td style="text-align:right;">
425416
</td>
<td style="text-align:left;">
RPL30
</td>
<td style="text-align:right;">
425416
</td>
<td style="text-align:left;">
NP\_001007968.1
</td>
<td style="text-align:right;">
-2.44
</td>
<td style="text-align:right;">
-2.44
</td>
<td style="text-align:right;">
6.58
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
NP\_001171687.1
</td>
<td style="text-align:right;">
417907
</td>
<td style="text-align:left;">
NDUFA12
</td>
<td style="text-align:right;">
417907
</td>
<td style="text-align:left;">
NP\_001171687.1
</td>
<td style="text-align:right;">
-2.87
</td>
<td style="text-align:right;">
-2.87
</td>
<td style="text-align:right;">
6.02
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
-1*male.hypothalamus.control 1*male.hypothalamus.n9
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

    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='MH_n9C', frame.plot=F)

![](../figures/hyp/01-contrasts-18.png)

volcano plots
=============

    # from http://www.compbio.dundee.ac.uk/user/pschofield/Projects/teaching_pg/workshops/biocDGE.html#maplots

    lrt <- glmLRT(fit,coef=2)
    topTags(lrt)

    ## Coefficient:  female.hypothalamus.control 
    ##                row.names         Name    geneid       entrezid      logFC
    ## XP_004950299.2 107049137 LOC107049137 107049137 XP_004950299.2  4.4008201
    ## XP_015146711.1    424629         FAF1    424629 XP_015146711.1 -0.7307621
    ## XP_015129796.1 100859555       CLAPS2 100859555 XP_015129796.1  2.3909235
    ## NP_001263232.1    396544         NACA    396544 NP_001263232.1  1.4639432
    ## XP_001234052.1    770722        RPS25    770722 XP_001234052.1  1.8754281
    ## XP_004936123.1    422772       GNPDA2    422772 XP_004936123.1  1.6739571
    ## XP_015155333.1    420089        MLLT1    420089 XP_015155333.1 -0.8824521
    ## NP_001157867.1    419653        EIF3I    419653 NP_001157867.1  1.4749709
    ## XP_004939591.1    420742      HERPUD2    420742 XP_004939591.1  0.7906706
    ## NP_006917.1       807639         COX1    807639    NP_006917.1  1.5368638
    ##                   logCPM        LR       PValue          FDR
    ## XP_004950299.2  3.092408 136.93628 1.245174e-31 1.859916e-27
    ## XP_015146711.1  6.565959 115.29650 6.776817e-27 5.061266e-23
    ## XP_015129796.1  6.343011 103.14518 3.114660e-24 1.550789e-20
    ## NP_001263232.1  6.546654 102.23970 4.919396e-24 1.837026e-20
    ## XP_001234052.1  7.923684  98.91595 2.634486e-23 7.202575e-20
    ## XP_004936123.1  5.224113  98.73046 2.893181e-23 7.202575e-20
    ## XP_015155333.1  6.134614  95.20223 1.718963e-22 3.668021e-19
    ## NP_001157867.1  6.735818  93.52854 4.003860e-22 7.475706e-19
    ## XP_004939591.1  5.528968  92.61971 6.337311e-22 1.051782e-18
    ## NP_006917.1    16.110350  92.40566 7.061200e-22 1.054731e-18

    tt <- topTags(lrt,n=10000)$table

    ggplot(data=tt) + geom_point(aes(x=logFC,y=-log(FDR),color=logCPM)) +
      scale_colour_gradientn(colours=c("#000000" ,"#FF0000" ))

![](../figures/hyp/volcanoplots-1.png)

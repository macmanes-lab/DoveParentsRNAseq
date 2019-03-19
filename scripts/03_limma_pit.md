DESeq2 is *not* recommended for experiments with more than 100 samples
([see Mike Love’s
post](https://mikelove.wordpress.com/2016/09/28/deseq2-or-edger/)), so I
decided to try the limma package. I followed [this
tutorial](https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html).

    library(tidyverse)

    ## ── Attaching packages ─────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.0       ✔ purrr   0.3.1  
    ## ✔ tibble  2.0.1       ✔ dplyr   0.8.0.1
    ## ✔ tidyr   0.8.3       ✔ stringr 1.4.0  
    ## ✔ readr   1.3.1       ✔ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────── tidyverse_conflicts() ──
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

For color coding, I used this tutorial for guidance
<a href="https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html" class="uri">https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html</a>.

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

specify contrasts and make MA plots
===================================

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
                 
                 MP_CB = male.pituitary.control - male.pituitary.bldg,
                 MP_BL = male.pituitary.bldg - male.pituitary.lay,
                 MP_Li3 = male.pituitary.lay - male.pituitary.inc.d3,
                 MP_i39 = male.pituitary.inc.d3 - male.pituitary.inc.d9,
                 MP_i917 = male.pituitary.inc.d9 - male.pituitary.inc.d17,
                 MP_i17H = male.pituitary.inc.d17 - male.pituitary.hatch,
                 MP_H5 = male.pituitary.hatch -  male.pituitary.n5,
                 MP_n59 = male.pituitary.n5 - male.pituitary.n9,
                 MP_n9C = male.pituitary.n9 - male.pituitary.control,
                 
                 FP_n9B = female.pituitary.n9 - female.pituitary.bldg,
                 MP_n9B = male.pituitary.n9 - male.pituitary.bldg,
    levels=parentaldesign)


    mycontrasts <- c("FP_CB", "FP_BL", "FP_Li3", "FP_i39", "FP_i917", "FP_i17H", "FP_H5", "FP_n59", "FP_n9C",
                     "MP_CB", "MP_BL", "MP_Li3", "MP_i39", "MP_i917", "MP_i17H", "MP_H5", "MP_n59", "MP_n9C",
                     "FP_n9B", "MP_n9B")


    printplotcontrasts <- function(whichcontrast){
      cont <- whichcontrast
      summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))
      return(kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1))
      return(plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main=whichcontrast, frame.plot=F))
    }

    for(i in mycontrasts){
      printplotcontrasts(i)
    }

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

This anlaysis will *exclude* the control timepoint and will example
males and females separately.

    library(tidyverse)

    ## ── Attaching packages ───────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.0       ✔ purrr   0.3.1  
    ## ✔ tibble  2.0.1       ✔ dplyr   0.8.0.1
    ## ✔ tidyr   0.8.3       ✔ stringr 1.4.0  
    ## ✔ readr   1.3.1       ✔ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────── tidyverse_conflicts() ──
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

    # load custom functions
    source("../R/functions.R") 

    knitr::opts_chunk$set(fig.path = '../figures/gon_F/',cache=TRUE)

Females first
=============

    # import "colData" which contains sample information and "countData" which contains read counts
    colData <- read.csv("../results/00_colData_characterization.csv", header = T, row.names = 1)
    countData <- read.csv("../results/00_countData_characterization.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../results/00_geneinfo.csv", row.names = 1)

    colData <- colData %>%
      dplyr::filter(grepl('gonad', tissue)) %>%
      dplyr::filter(treatment != "control") %>%
      dplyr::filter(sex == "female") %>%
      droplevels()
    row.names(colData) <- colData$V1

    # print sample sizes
    colData %>% select(treatment, tissue)  %>%  summary()

    ##    treatment    tissue  
    ##  inc.d9 :13   gonad:85  
    ##  inc.d17:11             
    ##  n9     :11             
    ##  bldg   :10             
    ##  hatch  :10             
    ##  inc.d3 :10             
    ##  (Other):20

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
    ## 14878    59

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

    ## Disp = 0.16993 , BCV = 0.4122

    parentalobject <- estimateGLMTrendedDisp(parentalobject, parentaldesign)
    parentalobject <- estimateGLMTagwiseDisp(parentalobject, parentaldesign)

    #  perform likelihood ratio test and thresholded testing
    fit <- glmFit( parentalobject, parentaldesign, robust=T)
    tr <- glmTreat(fit, lfc = 1)
    topTags(tr)

    ## Coefficient:  female.gonad.n9 
    ##                row.names         Name    geneid       entrezid     logFC
    ## XP_004937800.2 101747255 LOC101747255 101747255 XP_004937800.2  6.013288
    ## XP_420169.3       422172        FDX1L    422172    XP_420169.3 -3.034067
    ## XP_425958.3       428397      GHRH-LR    428397    XP_425958.3  3.508344
    ## XP_004939031.1 101748683 LOC101748683 101748683 XP_004939031.1  5.496202
    ## XP_004949005.1    420152        ZNRF4    420152 XP_004949005.1  5.119274
    ## XP_015135864.1    427459      SLC28A3    427459 XP_015135864.1  3.727855
    ## XP_414866.2       416566         AQP8    416566    XP_414866.2 -4.115631
    ## NP_989747.1       386578       CHRNA3    386578    NP_989747.1 -4.157425
    ## NP_001001194.1    407777        AvBD7    407777 NP_001001194.1 -3.485892
    ## XP_015153171.1    418420       TAGLN3    418420 XP_015153171.1 -2.771113
    ##                unshrunk.logFC     logCPM       PValue        FDR
    ## XP_004937800.2       6.613134 -0.5948276 3.674670e-06 0.03451517
    ## XP_420169.3         -3.035814  4.8806343 4.621433e-06 0.03451517
    ## XP_425958.3          3.528872  1.4109742 7.983999e-06 0.03688867
    ## XP_004939031.1       5.831383 -0.4384580 1.017902e-05 0.03688867
    ## XP_004949005.1       5.424871 -0.4533770 1.297026e-05 0.03688867
    ## XP_015135864.1       3.813099 -0.1112861 1.481770e-05 0.03688867
    ## XP_414866.2         -4.420396 -0.7032106 2.246445e-05 0.04246466
    ## NP_989747.1         -4.206947  0.7086444 2.274334e-05 0.04246466
    ## NP_001001194.1      -3.505514  2.6512374 6.508329e-05 0.10801656
    ## XP_015153171.1      -2.798100  0.5280188 9.768677e-05 0.14306495

    head(tr$table)

    ##                      logFC unshrunk.logFC    logCPM       PValue
    ## NP_001001127.1  0.48390412     0.48472128  3.431749 7.805434e-01
    ## NP_001001129.1 -0.01647735    -0.01651362  3.535450 9.877572e-01
    ## NP_001001189.1 -0.18074725    -0.18079334  5.730961 9.999990e-01
    ## NP_001001194.1 -3.48589224    -3.50551447  2.651237 6.508329e-05
    ## NP_001001195.1  0.89056964     0.94387531 -1.059992 2.980430e-01
    ## NP_001001201.1  0.38736800     0.38783239  2.963980 9.031657e-01

plotMDS (multidimential scaling)
================================

    levels(colData$treatment)

    ## [1] "bldg"    "hatch"   "inc.d17" "inc.d3"  "inc.d9"  "lay"     "n5"     
    ## [8] "n9"

    col.treatment <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6")[colData$treatment]
    myfill <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6")

    plotMDS(parentalobject, cex = 0.5, labels = colData$treatment, col=col.treatment)
    title("Female Gonad MDS 1v2")

![](../figures/gon_F/plotMDS-1.png)

    plotMDS(parentalobject, cex = 0.5, labels = colData$treatment, col=col.treatment, dim=c(3,4))
    title("Female Gonad MDS 3v4")

![](../figures/gon_F/plotMDS-2.png)

    plotMDS(parentalobject, cex = 0.5, labels = colData$treatment, col=col.treatment, dim=c(5,6))
    title("Female Gonad MDS 5v6")

![](../figures/gon_F/plotMDS-3.png)

For color coding, I used this tutorial for guidance
<a href="https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html" class="uri">https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html</a>.

specify contrasts and make MA plots
===================================

    # view all levels
    levels(colData$group)

    ## [1] "female.gonad.bldg"    "female.gonad.hatch"   "female.gonad.inc.d17"
    ## [4] "female.gonad.inc.d3"  "female.gonad.inc.d9"  "female.gonad.lay"    
    ## [7] "female.gonad.n5"      "female.gonad.n9"

    # subset of conrasts - sex specific comparing hatch to lay
    my.contrasts <- makeContrasts(
                 FG_BL = female.gonad.bldg - female.gonad.lay,
                 FG_Li3 = female.gonad.lay - female.gonad.inc.d3,
                 FG_i3i9 = female.gonad.inc.d3 - female.gonad.inc.d9,
                 FG_i9i17 = female.gonad.inc.d9 - female.gonad.inc.d17,
                 FG_i17H = female.gonad.inc.d17 - female.gonad.hatch,
                 FG_Hn5 = female.gonad.hatch -  female.gonad.n5,
                 FG_n5n9 = female.gonad.n5 - female.gonad.n9,
                 
                 FG_Bi3 = female.gonad.bldg - female.gonad.inc.d3,
                 FG_Bi9 = female.gonad.bldg - female.gonad.inc.d9,
                 FG_Bi17 = female.gonad.bldg - female.gonad.inc.d17,
                 FG_BH = female.gonad.bldg - female.gonad.hatch,
                 FG_Bn5 = female.gonad.bldg - female.gonad.n5,
                 FG_Bn9 = female.gonad.bldg - female.gonad.n9,
                 
                 FG_Li9 = female.gonad.lay - female.gonad.inc.d9,
               FG_Li17 = female.gonad.lay - female.gonad.inc.d17, 
                 FG_LH = female.gonad.lay - female.gonad.hatch,
                 FG_Ln5 = female.gonad.lay - female.gonad.n5,
                 FG_Ln9 = female.gonad.lay - female.gonad.n9,

                 FG_i3i17 = female.gonad.inc.d3 - female.gonad.inc.d17,
                 FG_i3H = female.gonad.inc.d3 - female.gonad.hatch,
                 FG_i3n5 = female.gonad.inc.d3 - female.gonad.n5,
                 FG_i3n9 = female.gonad.inc.d3 - female.gonad.n9,
                
                 FG_i9H = female.gonad.inc.d9 - female.gonad.hatch,
                 FG_i9n5 = female.gonad.inc.d9 - female.gonad.n5,
                 FG_i9n9 = female.gonad.inc.d9 - female.gonad.n9,
                 
                 FG_i17n5 = female.gonad.inc.d17 - female.gonad.n5,
                 FG_i17n9 = female.gonad.inc.d17 - female.gonad.n9,
                 
                 FG_Hn9 = female.gonad.hatch -  female.gonad.n9,
                 
    levels=parentaldesign)

    colnames(my.contrasts)

    ##  [1] "FG_BL"    "FG_Li3"   "FG_i3i9"  "FG_i9i17" "FG_i17H"  "FG_Hn5"  
    ##  [7] "FG_n5n9"  "FG_Bi3"   "FG_Bi9"   "FG_Bi17"  "FG_BH"    "FG_Bn5"  
    ## [13] "FG_Bn9"   "FG_Li9"   "FG_Li17"  "FG_LH"    "FG_Ln5"   "FG_Ln9"  
    ## [19] "FG_i3i17" "FG_i3H"   "FG_i3n5"  "FG_i3n9"  "FG_i9H"   "FG_i9n5" 
    ## [25] "FG_i9n9"  "FG_i17n5" "FG_i17n9" "FG_Hn9"

    # create a list with all the two way contrasts

    mycontrasts <- colnames(my.contrasts)
    mycontrasts

    ##  [1] "FG_BL"    "FG_Li3"   "FG_i3i9"  "FG_i9i17" "FG_i17H"  "FG_Hn5"  
    ##  [7] "FG_n5n9"  "FG_Bi3"   "FG_Bi9"   "FG_Bi17"  "FG_BH"    "FG_Bn5"  
    ## [13] "FG_Bn9"   "FG_Li9"   "FG_Li17"  "FG_LH"    "FG_Ln5"   "FG_Ln9"  
    ## [19] "FG_i3i17" "FG_i3H"   "FG_i3n5"  "FG_i3n9"  "FG_i9H"   "FG_i9n5" 
    ## [25] "FG_i9n9"  "FG_i17n5" "FG_i17n9" "FG_Hn9"

    # use the printplotcontrasts function to print summary stats and a volcano plot

    for(i in mycontrasts){
      printplotcontrasts(i)
    }

    ##        1*female.gonad.bldg -1*female.gonad.lay
    ## Down                                     14833
    ## NotSig                                     104
    ## Up                                           0

![](../figures/gon_F/MDSplots-1.png)

    ## NULL
    ##        -1*female.gonad.inc.d3 1*female.gonad.lay
    ## Down                                           1
    ## NotSig                                     14929
    ## Up                                             7

![](../figures/gon_F/MDSplots-2.png)

    ## NULL
    ##        1*female.gonad.inc.d3 -1*female.gonad.inc.d9
    ## Down                                              0
    ## NotSig                                        14936
    ## Up                                                1

![](../figures/gon_F/MDSplots-3.png)

    ## NULL
    ##        -1*female.gonad.inc.d17 1*female.gonad.inc.d9
    ## Down                                               0
    ## NotSig                                         14937
    ## Up                                                 0

![](../figures/gon_F/MDSplots-4.png)

    ## NULL
    ##        -1*female.gonad.hatch 1*female.gonad.inc.d17
    ## Down                                              0
    ## NotSig                                        14937
    ## Up                                                0

![](../figures/gon_F/MDSplots-5.png)

    ## NULL
    ##        1*female.gonad.hatch -1*female.gonad.n5
    ## Down                                         0
    ## NotSig                                   14937
    ## Up                                           0

![](../figures/gon_F/MDSplots-6.png)

    ## NULL
    ##        1*female.gonad.n5 -1*female.gonad.n9
    ## Down                                      2
    ## NotSig                                14935
    ## Up                                        0

![](../figures/gon_F/MDSplots-7.png)

    ## NULL
    ##        1*female.gonad.bldg -1*female.gonad.inc.d3
    ## Down                                        14808
    ## NotSig                                        129
    ## Up                                              0

![](../figures/gon_F/MDSplots-8.png)

    ## NULL
    ##        1*female.gonad.bldg -1*female.gonad.inc.d9
    ## Down                                        14816
    ## NotSig                                        121
    ## Up                                              0

![](../figures/gon_F/MDSplots-9.png)

    ## NULL
    ##        1*female.gonad.bldg -1*female.gonad.inc.d17
    ## Down                                         14784
    ## NotSig                                         153
    ## Up                                               0

![](../figures/gon_F/MDSplots-10.png)

    ## NULL
    ##        1*female.gonad.bldg -1*female.gonad.hatch
    ## Down                                       14773
    ## NotSig                                       164
    ## Up                                             0

![](../figures/gon_F/MDSplots-11.png)

    ## NULL
    ##        1*female.gonad.bldg -1*female.gonad.n5
    ## Down                                    14781
    ## NotSig                                    156
    ## Up                                          0

![](../figures/gon_F/MDSplots-12.png)

    ## NULL
    ##        1*female.gonad.bldg -1*female.gonad.n9
    ## Down                                    14840
    ## NotSig                                     97
    ## Up                                          0

![](../figures/gon_F/MDSplots-13.png)

    ## NULL
    ##        -1*female.gonad.inc.d9 1*female.gonad.lay
    ## Down                                           0
    ## NotSig                                     14913
    ## Up                                            24

![](../figures/gon_F/MDSplots-14.png)

    ## NULL
    ##        -1*female.gonad.inc.d17 1*female.gonad.lay
    ## Down                                            0
    ## NotSig                                      14913
    ## Up                                             24

![](../figures/gon_F/MDSplots-15.png)

    ## NULL
    ##        -1*female.gonad.hatch 1*female.gonad.lay
    ## Down                                          0
    ## NotSig                                    14929
    ## Up                                            8

![](../figures/gon_F/MDSplots-16.png)

    ## NULL
    ##        1*female.gonad.lay -1*female.gonad.n5
    ## Down                                       0
    ## NotSig                                 14929
    ## Up                                         8

![](../figures/gon_F/MDSplots-17.png)

    ## NULL
    ##        1*female.gonad.lay -1*female.gonad.n9
    ## Down                                       1
    ## NotSig                                 14934
    ## Up                                         2

![](../figures/gon_F/MDSplots-18.png)

    ## NULL
    ##        -1*female.gonad.inc.d17 1*female.gonad.inc.d3
    ## Down                                               0
    ## NotSig                                         14937
    ## Up                                                 0

![](../figures/gon_F/MDSplots-19.png)

    ## NULL
    ##        -1*female.gonad.hatch 1*female.gonad.inc.d3
    ## Down                                             0
    ## NotSig                                       14936
    ## Up                                               1

![](../figures/gon_F/MDSplots-20.png)

    ## NULL
    ##        1*female.gonad.inc.d3 -1*female.gonad.n5
    ## Down                                          1
    ## NotSig                                    14934
    ## Up                                            2

![](../figures/gon_F/MDSplots-21.png)

    ## NULL
    ##        1*female.gonad.inc.d3 -1*female.gonad.n9
    ## Down                                         11
    ## NotSig                                    14924
    ## Up                                            2

![](../figures/gon_F/MDSplots-22.png)

    ## NULL
    ##        -1*female.gonad.hatch 1*female.gonad.inc.d9
    ## Down                                             0
    ## NotSig                                       14937
    ## Up                                               0

![](../figures/gon_F/MDSplots-23.png)

    ## NULL
    ##        1*female.gonad.inc.d9 -1*female.gonad.n5
    ## Down                                          0
    ## NotSig                                    14937
    ## Up                                            0

![](../figures/gon_F/MDSplots-24.png)

    ## NULL
    ##        1*female.gonad.inc.d9 -1*female.gonad.n9
    ## Down                                         14
    ## NotSig                                    14923
    ## Up                                            0

![](../figures/gon_F/MDSplots-25.png)

    ## NULL
    ##        1*female.gonad.inc.d17 -1*female.gonad.n5
    ## Down                                           0
    ## NotSig                                     14937
    ## Up                                             0

![](../figures/gon_F/MDSplots-26.png)

    ## NULL
    ##        1*female.gonad.inc.d17 -1*female.gonad.n9
    ## Down                                          19
    ## NotSig                                     14918
    ## Up                                             0

![](../figures/gon_F/MDSplots-27.png)

    ## NULL
    ##        1*female.gonad.hatch -1*female.gonad.n9
    ## Down                                         1
    ## NotSig                                   14936
    ## Up                                           0

![](../figures/gon_F/MDSplots-28.png)

    ## NULL

    for(i in mycontrasts){
      plotVolcanos(i)
    }

![](../figures/gon_F/volcanoplots-1.png)![](../figures/gon_F/volcanoplots-2.png)![](../figures/gon_F/volcanoplots-3.png)![](../figures/gon_F/volcanoplots-4.png)![](../figures/gon_F/volcanoplots-5.png)![](../figures/gon_F/volcanoplots-6.png)![](../figures/gon_F/volcanoplots-7.png)![](../figures/gon_F/volcanoplots-8.png)![](../figures/gon_F/volcanoplots-9.png)![](../figures/gon_F/volcanoplots-10.png)![](../figures/gon_F/volcanoplots-11.png)![](../figures/gon_F/volcanoplots-12.png)![](../figures/gon_F/volcanoplots-13.png)![](../figures/gon_F/volcanoplots-14.png)![](../figures/gon_F/volcanoplots-15.png)![](../figures/gon_F/volcanoplots-16.png)![](../figures/gon_F/volcanoplots-17.png)![](../figures/gon_F/volcanoplots-18.png)![](../figures/gon_F/volcanoplots-19.png)![](../figures/gon_F/volcanoplots-20.png)![](../figures/gon_F/volcanoplots-21.png)![](../figures/gon_F/volcanoplots-22.png)![](../figures/gon_F/volcanoplots-23.png)![](../figures/gon_F/volcanoplots-24.png)![](../figures/gon_F/volcanoplots-25.png)![](../figures/gon_F/volcanoplots-26.png)![](../figures/gon_F/volcanoplots-27.png)![](../figures/gon_F/volcanoplots-28.png)

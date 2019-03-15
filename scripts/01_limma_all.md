DESeq2 is *not* recommended for experiments with more than 100 samples
([see Mike Love’s
post](https://mikelove.wordpress.com/2016/09/28/deseq2-or-edger/)), so I
decided to try the limma package. I followed [this
tutorial](https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html).

    library(limma)
    library(Glimma)
    library(edgeR)
    library(kableExtra)

    knitr::opts_chunk$set(fig.path = '../figures/',cache=TRUE)

First, I read in the data I processed in 00\_datawrangling.Rmd.

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
    ## 14867    70

    keep_genes <- rowSums(cpms >= 1) >= 10
    dge <- parentalobject[keep_genes, ]

    # magic
    parentaldesign <- model.matrix(~ colData$group )
    colnames(parentaldesign) <- levels(colData$group)
    parentalobject <- calcNormFactors(parentalobject)
    parentalobject <- estimateCommonDisp(parentalobject)
    parentalobject <- estimateTagwiseDisp(parentalobject)
    parentalobject <- estimateDisp(parentalobject, parentaldesign)
    parentalobject <- estimateGLMCommonDisp(parentalobject, parentaldesign, verbose=TRUE)

    ## Disp = 0.09128 , BCV = 0.3021

    parentalobject <- estimateGLMTrendedDisp(parentalobject, parentaldesign)
    parentalobject <- estimateGLMTagwiseDisp(parentalobject, parentaldesign)
    fit <- glmFit( parentalobject, parentaldesign, robust=T)
    treatres <- glmTreat(fit, lfc = 1)

specify contrasts
=================

    levels(colData$group)

    ##  [1] "female.gonad.bldg"             "female.gonad.control"         
    ##  [3] "female.gonad.hatch"            "female.gonad.inc.d17"         
    ##  [5] "female.gonad.inc.d3"           "female.gonad.inc.d9"          
    ##  [7] "female.gonad.lay"              "female.gonad.m.hatch"         
    ##  [9] "female.gonad.m.inc.d17"        "female.gonad.m.inc.d3"        
    ## [11] "female.gonad.m.inc.d8"         "female.gonad.m.inc.d9"        
    ## [13] "female.gonad.m.n2"             "female.gonad.n5"              
    ## [15] "female.gonad.n9"               "female.hypothalamus.bldg"     
    ## [17] "female.hypothalamus.control"   "female.hypothalamus.hatch"    
    ## [19] "female.hypothalamus.inc.d17"   "female.hypothalamus.inc.d3"   
    ## [21] "female.hypothalamus.inc.d9"    "female.hypothalamus.lay"      
    ## [23] "female.hypothalamus.m.hatch"   "female.hypothalamus.m.inc.d17"
    ## [25] "female.hypothalamus.m.inc.d3"  "female.hypothalamus.m.inc.d8" 
    ## [27] "female.hypothalamus.m.inc.d9"  "female.hypothalamus.m.n2"     
    ## [29] "female.hypothalamus.n5"        "female.hypothalamus.n9"       
    ## [31] "female.pituitary.bldg"         "female.pituitary.control"     
    ## [33] "female.pituitary.hatch"        "female.pituitary.inc.d17"     
    ## [35] "female.pituitary.inc.d3"       "female.pituitary.inc.d9"      
    ## [37] "female.pituitary.lay"          "female.pituitary.m.hatch"     
    ## [39] "female.pituitary.m.inc.d17"    "female.pituitary.m.inc.d3"    
    ## [41] "female.pituitary.m.inc.d8"     "female.pituitary.m.inc.d9"    
    ## [43] "female.pituitary.m.n2"         "female.pituitary.n5"          
    ## [45] "female.pituitary.n9"           "male.gonad.bldg"              
    ## [47] "male.gonad.control"            "male.gonad.hatch"             
    ## [49] "male.gonad.inc.d17"            "male.gonad.inc.d3"            
    ## [51] "male.gonad.inc.d9"             "male.gonad.lay"               
    ## [53] "male.gonad.m.hatch"            "male.gonad.m.inc.d17"         
    ## [55] "male.gonad.m.inc.d3"           "male.gonad.m.inc.d8"          
    ## [57] "male.gonad.m.inc.d9"           "male.gonad.m.n2"              
    ## [59] "male.gonad.n5"                 "male.gonad.n9"                
    ## [61] "male.hypothalamus.bldg"        "male.hypothalamus.control"    
    ## [63] "male.hypothalamus.hatch"       "male.hypothalamus.inc.d17"    
    ## [65] "male.hypothalamus.inc.d3"      "male.hypothalamus.inc.d9"     
    ## [67] "male.hypothalamus.lay"         "male.hypothalamus.m.hatch"    
    ## [69] "male.hypothalamus.m.inc.d17"   "male.hypothalamus.m.inc.d3"   
    ## [71] "male.hypothalamus.m.inc.d8"    "male.hypothalamus.m.inc.d9"   
    ## [73] "male.hypothalamus.m.n2"        "male.hypothalamus.n5"         
    ## [75] "male.hypothalamus.n9"          "male.pituitary.bldg"          
    ## [77] "male.pituitary.control"        "male.pituitary.hatch"         
    ## [79] "male.pituitary.inc.d17"        "male.pituitary.inc.d3"        
    ## [81] "male.pituitary.inc.d9"         "male.pituitary.lay"           
    ## [83] "male.pituitary.m.hatch"        "male.pituitary.m.inc.d17"     
    ## [85] "male.pituitary.m.inc.d3"       "male.pituitary.m.inc.d8"      
    ## [87] "male.pituitary.m.inc.d9"       "male.pituitary.m.n2"          
    ## [89] "male.pituitary.n5"             "male.pituitary.n9"

    my.contrasts <- makeContrasts(
                 FG_HL = female.gonad.hatch - female.gonad.lay,
                 FH_HL = female.hypothalamus.hatch - female.hypothalamus.lay,
                 FP_HL = male.pituitary.hatch - male.pituitary.lay,
                 MP_HL = male.pituitary.hatch - male.pituitary.lay,
                 MH_HL = male.hypothalamus.hatch - male.hypothalamus.lay,
                 MG_HL = male.gonad.hatch - male.gonad.lay,          
    levels=parentaldesign)

    cont <- "FG_HL"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*female.gonad.hatch -1*female.gonad.lay
    ## Down                                         40
    ## NotSig                                    14895
    ## Up                                            2

    #kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)
    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='Female gonad lay to hatch', frame.plot=F)

![](../figures/01-contrasts-1.png)

    cont <- "FH_HL"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*female.hypothalamus.hatch -1*female.hypothalamus.lay
    ## Down                                                        0
    ## NotSig                                                  14935
    ## Up                                                          2

    #kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)
    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='Female hypothalamus lay to hatch', frame.plot=F)

![](../figures/01-contrasts-2.png)

    cont <- "FP_HL"
    summary(decideTestsDGE(
        glmTreat(fit, contrast=my.contrasts[,cont], lfc = 1), 
        adjust.method="fdr", p.value=0.01))

    ##        1*male.pituitary.hatch -1*male.pituitary.lay
    ## Down                                              0
    ## NotSig                                        14882
    ## Up                                               55

    #kable(topTags(glmTreat(fit, contrast=my.contrasts[,cont]), n=5), digits=2, lfc = 1)
    plotMD(glmTreat(fit, contrast=my.contrasts[,cont], lfc=1), main='Female pituitary lay to hatch', frame.plot=F)

![](../figures/01-contrasts-3.png)

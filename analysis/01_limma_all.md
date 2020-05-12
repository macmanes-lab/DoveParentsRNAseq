DESeq2 is *not* recommended for experiments with more than 100 samples
([see Mike Love’s
post](https://mikelove.wordpress.com/2016/09/28/deseq2-or-edger/)), so I
decided to try the limma package. I followed [this
tutorial](https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html).

    library(limma)
    library(Glimma)
    library(edgeR)
    library(kableExtra)
    library(cowplot)

    ## Loading required package: ggplot2

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    knitr::opts_chunk$set(fig.path = '../figures/01_limma/')

First, I read in the data I processed in 00\_datawrangling.Rmd.

    # import "colData" which contains sample information and "countData" which contains read counts
    colData <- read.csv("../metadata/00_colData.csv", header = T, row.names = 1)
    countData <- read.csv("../results/00_counts.csv", header = T, row.names = 1)

    #countData  <- head(countData, 1000) # suset for quick analysis

    geneinfo <- row.names(countData)

Then, I followed the steps from
<a href="https://github.com/macmanes-lab/RockDove/blob/master/parental_care/parental_analysis.Rmd" class="uri">https://github.com/macmanes-lab/RockDove/blob/master/parental_care/parental_analysis.Rmd</a>.

    head(colData$group)

    ## [1] male.gonads.control       male.hypothalamus.control
    ## [3] male.pituitary.control    male.gonads.control      
    ## [5] male.hypothalamus.control male.pituitary.control   
    ## 96 Levels: female.gonads.bldg ... male.pituitary.prolong

    # create a large DGEList with 3 elements
    parentalobject <- DGEList(counts=countData, genes=geneinfo, group=colData$group)

    # transform raw counts to countspermillion
    cpms <- cpm(parentalobject)

    # calculate number of lowly lowly expressed genes and remove them
    table(rowSums(parentalobject$counts==0)==10)

    ## 
    ## FALSE  TRUE 
    ## 13908    58

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

    ## Disp = 0.08632 , BCV = 0.2938

    parentalobject <- estimateGLMTrendedDisp(parentalobject, parentaldesign)
    parentalobject <- estimateGLMTagwiseDisp(parentalobject, parentaldesign)

    # find and print data
    names(parentalobject)

    ##  [1] "counts"             "samples"            "genes"             
    ##  [4] "common.dispersion"  "pseudo.counts"      "pseudo.lib.size"   
    ##  [7] "AveLogCPM"          "prior.df"           "prior.n"           
    ## [10] "tagwise.dispersion" "span"               "design"            
    ## [13] "trended.dispersion" "trend.method"

    #head(countData)
    #head(parentalobject$counts)
    #head(parentalobject$pseudo.counts)


    write.csv(parentalobject$pseudo.counts, "../results/01_limma.csv")

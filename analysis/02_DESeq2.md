    library(tidyverse)

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(DESeq2)

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
    ##     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
    ##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
    ##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
    ##     setdiff, sort, table, tapply, union, unique, unsplit, which,
    ##     which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## Loading required package: BiocParallel

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

    library(BiocParallel)
    register(MulticoreParam(6))
    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    source("../R/functions.R")  # load custom functions 
    source("../R/themes.R")  # load custom themes and color palletes

    knitr::opts_chunk$set(fig.path = '../figures/DESeq2/',  message=F, comment=F, warning=F)

DESeq2
======

-   DESeq2 was not designed to run on 100+ samples. But, I really like
    it, so I do it anyways. Some of these commands take like 15 min to
    run using 6 cores. \_ Also, I don’t run every chuck every time. When
    I want new analysese, I add them to the bottom and set `eval = F`
    for chunks I don’t want to rerurn

<!-- -->

    countData <- read.csv("../results/00_counts.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)

    colData <- read.csv("../metadata/00_colData.csv", header = T, row.names = 1) %>%
      mutate(tissue = factor(tissue, levels = tissuelevels),
             treatment <- factor(treatment, levels = alllevels)) %>%
      mutate(sextissue = as.factor(paste(sex, tissue, sep = "_"))) 
    colData <- colData %>% drop_na()
    # 987 samples, 8 variables

    ncol(countData) == nrow(colData)

    FALSE [1] TRUE

tissue and sex differences
--------------------------

    ## first or sex

    for(i in levels(colData$tissue)){
     
      newcolData <- subsetcolData4(colData, i)
      print(head(newcolData))
      # save counts that match colData
      savecols <- as.character(newcolData$V1) 
      savecols <- as.vector(savecols) 
      
      newcountData <- countData %>% dplyr::select(one_of(savecols)) 
      
      dds <- DESeqDataSetFromMatrix(countData = newcountData,
                                    colData = newcolData,
                                    design = ~ sex )

      
      dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
      print(dds)
      print(dim(dds))
      dds <- DESeq(dds, parallel = TRUE) # Differential expression analysis
      
      vsd <- as.data.frame(assay(vst(dds, blind=FALSE)))
      
      myfilename = paste0("../results/DEseq2/sex/", i, "_vsd.csv")
      write.csv(vsd, myfilename)

      print(head(vsd))

      # save differential gene expression results
      sex <- createDEGdfsavesex("male", "female", i)
    }

     
    ## then for tissue

    dds <- DESeqDataSetFromMatrix(countData = countData,
                                    colData = colData,
                                    design = ~ tissue )
    dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
    print(dds)
    print(dim(dds))
    dds <- DESeq(dds, parallel = TRUE) # Differential expression analysis
      
    vsd <- as.data.frame(assay(vst(dds, blind=FALSE)))
    write.csv(vsd, "../results/DEseq2/tissue/tissuevsds.csv")
    print(head(vsd))

    # save differential gene expression results
    HP <- createDEGdfsavestissue("pituitary", "hypothalamus")
    HG <- createDEGdfsavestissue("gonad", "hypothalamus")
    PG <- createDEGdfsavestissue("gonad", "pituitary")
      
    str(HP)

FOR DEGs of all TREATMENT comparisons, char y manip
---------------------------------------------------

    createDEGdftreatment <- function(up, down, mytissue){
      
      res <- results(dds, contrast = c("treatment", up, down), 
                     independentFiltering = T, alpha = 0.1)
      
      DEGs <- data.frame(gene = row.names(res),
                            padj = res$padj, 
                            logpadj = -log10(res$padj),
                            lfc = res$log2FoldChange,
                            sextissue = mytissue)
      DEGs <- na.omit(DEGs)
      DEGs <- DEGs %>%
        dplyr::mutate(direction = ifelse(DEGs$lfc > 0 & DEGs$padj < 0.1, 
                                         yes = up, no = ifelse(DEGs$lfc < 0 & DEGs$padj < 0.1, 
                                                               yes = down, no = "NS"))) %>% 
        dplyr::arrange(desc(lfc)) 
      
      DEGs$direction <- factor(DEGs$direction, levels = c(down, "NS", up)) 
      
      # write DEGsframe of only significant genes
      DEGs <- DEGs %>% dplyr::filter(direction != "NS")
      print(str(DEGs))
      
      partialfilename = paste("_", down, "_", up, sep = "")
      myfilename = paste0("../results/DESeq2/treatment/", mytissue, partialfilename, "_DEGs.csv")
      
      write.csv(DEGs, myfilename, row.names = F)
      # return DEGs frome with all data, included NS genes
      #print(head(DEGs))
    }  


    for(i in levels(colData$sextissue)){
      
      newcolData <- subsetcolData2(colData, i)
      
      # save counts that match colData
      savecols <- as.character(newcolData$V1) 
      savecols <- as.vector(savecols) 
      
      newcountData <- countData %>% dplyr::select(one_of(savecols)) 
      
      dds <- DESeqDataSetFromMatrix(countData = newcountData,
                                    colData = newcolData,
                                    design = ~ treatment )
      dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
      print(dds)
      print(dim(dds))
      dds <- DESeq(dds, parallel = TRUE) # Differential expression analysis
      
      vsd <- as.data.frame(assay(vst(dds, blind=FALSE)))
      
      myfilename = paste0("../results/DEseq2/treatment/", i, "_vsd.csv")
      write.csv(vsd, myfilename)
      
      #return(dds)
      #return(vsd)
      #print(head(vsd))

      # save differential gene expression results
      
      inc.d9.early <- createDEGdftreatment("inc.d9", "early",  i) 
      earyl.hatch <- createDEGdftreatment("early", "hatch",  i) 
      
      inc.d17.prolong <- createDEGdftreatment("inc.d17", "prolong",  i) 
      prolong.hatch <- createDEGdftreatment("prolong", "hatch",  i) 
      
      hatch.extend <- createDEGdftreatment("hatch", "extend",  i) 
      extend.n5 <- createDEGdftreatment("extend", "n5",  i) 
     

    }

    FALSE class: DESeqDataSet 
    FALSE dim: 13801 167 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13801): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(167): L.G118_female_gonad_control
    FALSE   R.G106_female_gonad_control ... y97.x_female_gonad_n9
    FALSE   y98.g54_female_gonad_m.hatch
    FALSE colData names(9): V1 bird ... treatment <- factor(treatment,
    FALSE   levels = alllevels) sextissue
    FALSE [1] 13801   167
    FALSE 'data.frame': 22 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13801 levels "A2ML1","A2ML2",..: 6703 564 882 11998 625 848 6820 11569 648 7787 ...
    FALSE  $ padj     : num  0.0452 0.00581 0.09426 0.0452 0.0452 ...
    FALSE  $ logpadj  : num  1.34 2.24 1.03 1.34 1.34 ...
    FALSE  $ lfc      : num  2.195 2.162 1.979 1.475 0.412 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "early","NS","inc.d9": 3 3 3 3 3 3 3 3 3 1 ...
    FALSE NULL
    FALSE 'data.frame': 20 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13801 levels "A2ML1","A2ML2",..: 5308 6634 6800 7808 145 3271 6572 6526 2801 4481 ...
    FALSE  $ padj     : num  1.61e-05 5.44e-10 3.43e-02 6.13e-02 3.43e-02 ...
    FALSE  $ logpadj  : num  4.79 9.26 1.47 1.21 1.47 ...
    FALSE  $ lfc      : num  18.25 15.9 4.98 3.42 1.93 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","early": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 2210 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13801 levels "A2ML1","A2ML2",..: 11007 6512 2238 7032 8765 11154 11724 8547 4207 2732 ...
    FALSE  $ padj     : num  6.36e-05 1.09e-04 4.39e-02 1.54e-02 7.31e-04 ...
    FALSE  $ logpadj  : num  4.2 3.96 1.36 1.81 3.14 ...
    FALSE  $ lfc      : num  5.86 3.53 3.53 2.94 2.88 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "prolong","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 1410 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13801 levels "A2ML1","A2ML2",..: 5308 6634 518 6936 1327 3534 6667 1698 2864 1951 ...
    FALSE  $ padj     : num  1.79e-11 6.60e-25 2.45e-03 7.15e-14 4.22e-05 ...
    FALSE  $ logpadj  : num  10.75 24.18 2.61 13.15 4.38 ...
    FALSE  $ lfc      : num  23.75 22.68 6.45 6.12 6.01 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","prolong": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 82 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13801 levels "A2ML1","A2ML2",..: 11130 8050 6633 7361 1386 7191 2155 7063 9014 5027 ...
    FALSE  $ padj     : num  0.07833 0.04583 0.07486 0.06883 0.00965 ...
    FALSE  $ logpadj  : num  1.11 1.34 1.13 1.16 2.02 ...
    FALSE  $ lfc      : num  5.56 5.53 5.45 2.62 2.26 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "extend","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 15 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13801 levels "A2ML1","A2ML2",..: 5308 6408 2364 6745 8047 6194 2735 6289 6636 11693 ...
    FALSE  $ padj     : num  3.85e-07 1.67e-05 8.91e-02 1.56e-06 6.36e-03 ...
    FALSE  $ logpadj  : num  6.41 4.78 1.05 5.81 2.2 ...
    FALSE  $ lfc      : num  20.34 9.42 4.22 3.79 3.22 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "n5","NS","extend": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE class: DESeqDataSet 
    FALSE dim: 13631 165 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13631): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(165): L.G118_female_hypothalamus_control.NYNO
    FALSE   R.G106_female_hypothalamus_control ...
    FALSE   y97.x_female_hypothalamus_n9 y98.g54_female_hypothalamus_m.hatch
    FALSE colData names(9): V1 bird ... treatment <- factor(treatment,
    FALSE   levels = alllevels) sextissue
    FALSE [1] 13631   165
    FALSE 'data.frame': 1 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13631 levels "A2ML1","A2ML2",..: 6303
    FALSE  $ padj     : num 3.04e-19
    FALSE  $ logpadj  : num 18.5
    FALSE  $ lfc      : num -18.4
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1
    FALSE  $ direction: Factor w/ 3 levels "early","NS","inc.d9": 1
    FALSE NULL
    FALSE 'data.frame': 4911 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13631 levels "A2ML1","A2ML2",..: 3118 6698 5764 11104 9643 5148 12470 10807 11830 6754 ...
    FALSE  $ padj     : num  0.00133 0.05905 0.07117 0.00974 0.04537 ...
    FALSE  $ logpadj  : num  2.88 1.23 1.15 2.01 1.34 ...
    FALSE  $ lfc      : num  3.72 3.28 3.15 3.13 3.09 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","early": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 2663 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13631 levels "A2ML1","A2ML2",..: 5102 8703 8436 8413 5101 5731 9159 6092 10945 7578 ...
    FALSE  $ padj     : num  0.001306 0.086044 0.00266 0.000464 0.002238 ...
    FALSE  $ logpadj  : num  2.88 1.07 2.58 3.33 2.65 ...
    FALSE  $ lfc      : num  4.22 3.97 3.91 3.58 3.36 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "prolong","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 4979 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13631 levels "A2ML1","A2ML2",..: 3118 2999 3010 10232 7927 6767 4913 5640 6165 495 ...
    FALSE  $ padj     : num  0.000294 0.002066 0.086761 0.00645 0.04897 ...
    FALSE  $ logpadj  : num  3.53 2.68 1.06 2.19 1.31 ...
    FALSE  $ lfc      : num  4.09 3.4 3.25 2.82 2.59 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","prolong": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 5499 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13631 levels "A2ML1","A2ML2",..: 8899 37 8703 8705 9159 6677 4643 10011 6423 12713 ...
    FALSE  $ padj     : num  3.61e-04 4.75e-05 4.24e-02 7.29e-02 8.77e-06 ...
    FALSE  $ logpadj  : num  3.44 4.32 1.37 1.14 5.06 ...
    FALSE  $ lfc      : num  6.06 5.19 4.01 3.98 3.79 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "extend","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13631 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 
    FALSE  $ direction: Factor w/ 3 levels "n5","NS","extend": 
    FALSE NULL
    FALSE class: DESeqDataSet 
    FALSE dim: 13554 165 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13554): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(165): L.G118_female_pituitary_control.NYNO
    FALSE   R.G106_female_pituitary_control ... y97.x_female_pituitary_n9
    FALSE   y98.g54_female_pituitary_m.hatch
    FALSE colData names(9): V1 bird ... treatment <- factor(treatment,
    FALSE   levels = alllevels) sextissue
    FALSE [1] 13554   165
    FALSE 'data.frame': 109 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13554 levels "A2ML1","A2ML2",..: 11177 560 949 4879 11766 10388 2514 1594 6471 13437 ...
    FALSE  $ padj     : num  0.0615 0.0301 0.0996 0.0783 0.0996 ...
    FALSE  $ logpadj  : num  1.21 1.52 1 1.11 1 ...
    FALSE  $ lfc      : num  3.55 1.85 1.81 1.64 1.16 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "early","NS","inc.d9": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 4198 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13554 levels "A2ML1","A2ML2",..: 6662 486 2472 9979 12048 5753 6855 4243 8569 4651 ...
    FALSE  $ padj     : num  0.00911 0.02835 0.09407 0.02783 0.04268 ...
    FALSE  $ logpadj  : num  2.04 1.55 1.03 1.56 1.37 ...
    FALSE  $ lfc      : num  3.73 3.22 3.08 2.73 2.68 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","early": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 1126 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13554 levels "A2ML1","A2ML2",..: 12894 6389 11132 6268 5468 5158 8191 2451 6761 6814 ...
    FALSE  $ padj     : num  0.0147 0.0317 0.0291 0.024 0.0762 ...
    FALSE  $ logpadj  : num  1.83 1.5 1.54 1.62 1.12 ...
    FALSE  $ lfc      : num  2.37 2.31 2.12 2.11 1.78 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "prolong","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 844 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13554 levels "A2ML1","A2ML2",..: 5753 8436 796 8110 9141 5394 7906 677 7858 3802 ...
    FALSE  $ padj     : num  0.0661 0.0487 0.0524 0.0045 0.0723 ...
    FALSE  $ logpadj  : num  1.18 1.31 1.28 2.35 1.14 ...
    FALSE  $ lfc      : num  2.4 1.55 1.48 1.44 1.39 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","prolong": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 1985 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13554 levels "A2ML1","A2ML2",..: 3497 12400 1435 3593 10729 10498 617 2532 5080 12635 ...
    FALSE  $ padj     : num  0.004702 0.014454 0.09531 0.021117 0.000208 ...
    FALSE  $ logpadj  : num  2.33 1.84 1.02 1.68 3.68 ...
    FALSE  $ lfc      : num  3.85 3.41 2.38 2.38 2.34 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "extend","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 363 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13554 levels "A2ML1","A2ML2",..: 7609 11166 1841 12948 12028 201 7893 4014 323 6919 ...
    FALSE  $ padj     : num  3.49e-02 5.04e-06 3.76e-02 2.97e-03 1.07e-04 ...
    FALSE  $ logpadj  : num  1.46 5.3 1.42 2.53 3.97 ...
    FALSE  $ lfc      : num  1.72 1.55 1.42 1.32 1.27 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "n5","NS","extend": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE class: DESeqDataSet 
    FALSE dim: 13825 163 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13825): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(163): L.Blu13_male_gonad_control.NYNO
    FALSE   L.G107_male_gonad_control ... y95.g131.x_male_gonad_inc.d9
    FALSE   y98.o50.x_male_gonad_inc.d3
    FALSE colData names(9): V1 bird ... treatment <- factor(treatment,
    FALSE   levels = alllevels) sextissue
    FALSE [1] 13825   163
    FALSE 'data.frame': 267 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13825 levels "A2ML1","A2ML2",..: 567 886 6724 10287 1625 10581 13784 7535 12798 10444 ...
    FALSE  $ padj     : num  0.000161 0.00104 0.004568 0.004146 0.00104 ...
    FALSE  $ logpadj  : num  3.79 2.98 2.34 2.38 2.98 ...
    FALSE  $ lfc      : num  2.46 2.16 1.93 1.89 1.65 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "early","NS","inc.d9": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13825 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "male_gonads": 
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","early": 
    FALSE NULL
    FALSE 'data.frame': 3031 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13825 levels "A2ML1","A2ML2",..: 6016 13557 6533 6820 2988 6689 8272 6945 270 1243 ...
    FALSE  $ padj     : num  0.000329 0.001628 0.000178 0.062787 0.077114 ...
    FALSE  $ logpadj  : num  3.48 2.79 3.75 1.2 1.11 ...
    FALSE  $ lfc      : num  4.33 4.14 3.68 2.72 2.44 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "prolong","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 1937 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13825 levels "A2ML1","A2ML2",..: 9417 2700 6629 3656 10450 3381 324 4258 2103 12452 ...
    FALSE  $ padj     : num  1.33e-11 4.33e-04 2.60e-06 9.46e-06 2.11e-03 ...
    FALSE  $ logpadj  : num  10.88 3.36 5.58 5.02 2.68 ...
    FALSE  $ lfc      : num  5.6 4.98 4.89 4.57 4.48 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","prolong": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 651 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13825 levels "A2ML1","A2ML2",..: 10287 904 567 6129 5107 8485 10581 10964 4641 13433 ...
    FALSE  $ padj     : num  0.0379 0.0239 0.0798 0.0935 0.0118 ...
    FALSE  $ logpadj  : num  1.42 1.62 1.1 1.03 1.93 ...
    FALSE  $ lfc      : num  1.39 1.35 1.28 1.21 1.09 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "extend","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 96 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13825 levels "A2ML1","A2ML2",..: 13070 503 12142 11039 3577 10335 9520 4448 3288 7851 ...
    FALSE  $ padj     : num  0.0289 0.055 0.0964 0.0823 0.055 ...
    FALSE  $ logpadj  : num  1.54 1.26 1.02 1.08 1.26 ...
    FALSE  $ lfc      : num  1.987 0.383 0.374 0.31 0.302 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "n5","NS","extend": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE class: DESeqDataSet 
    FALSE dim: 13595 162 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13595): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(162): L.Blu13_male_hypothalamus_control.NYNO
    FALSE   L.G107_male_hypothalamus_control ...
    FALSE   y95.g131.x_male_hypothalamus_inc.d9
    FALSE   y98.o50.x_male_hypothalamus_inc.d3
    FALSE colData names(9): V1 bird ... treatment <- factor(treatment,
    FALSE   levels = alllevels) sextissue
    FALSE [1] 13595   162
    FALSE 'data.frame': 1678 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13595 levels "A2ML1","A2ML2",..: 9128 8380 1558 6549 10915 4227 1778 4298 11413 13246 ...
    FALSE  $ padj     : num  7.43e-05 1.30e-04 7.43e-05 7.43e-05 1.02e-03 ...
    FALSE  $ logpadj  : num  4.13 3.89 4.13 4.13 2.99 ...
    FALSE  $ lfc      : num  4.38 3.9 3.11 2.6 2.44 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "early","NS","inc.d9": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 4768 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13595 levels "A2ML1","A2ML2",..: 9814 5090 13212 5089 6596 3908 8274 6472 3249 1809 ...
    FALSE  $ padj     : num  0.03192 0.00914 0.00274 0.00273 0.00178 ...
    FALSE  $ logpadj  : num  1.5 2.04 2.56 2.56 2.75 ...
    FALSE  $ lfc      : num  4.41 3.26 3.14 3.06 2.42 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","early": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 1422 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13595 levels "A2ML1","A2ML2",..: 8670 5090 8403 8095 11752 9128 8380 6159 12597 10915 ...
    FALSE  $ padj     : num  0.0939 0.0603 0.0474 0.0463 0.0951 ...
    FALSE  $ logpadj  : num  1.03 1.22 1.32 1.33 1.02 ...
    FALSE  $ lfc      : num  4.68 2.94 2.68 2.62 2.36 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "prolong","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 820 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13595 levels "A2ML1","A2ML2",..: 7491 8887 762 4612 5014 5785 11563 1239 13572 2203 ...
    FALSE  $ padj     : num  0.0802 0.0332 0.0339 0.076 0.0617 ...
    FALSE  $ logpadj  : num  1.1 1.48 1.47 1.12 1.21 ...
    FALSE  $ lfc      : num  2.55 2.34 1.93 1.49 1.41 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","prolong": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 2275 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13595 levels "A2ML1","A2ML2",..: 12447 8403 9128 7945 1558 8792 269 6549 8380 10915 ...
    FALSE  $ padj     : num  4.00e-04 6.00e-04 7.97e-05 3.87e-04 2.63e-04 ...
    FALSE  $ logpadj  : num  3.4 3.22 4.1 3.41 3.58 ...
    FALSE  $ lfc      : num  4.6 4.29 4.19 2.93 2.89 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "extend","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 1 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13595 levels "A2ML1","A2ML2",..: 12508
    FALSE  $ padj     : num 0.0063
    FALSE  $ logpadj  : num 2.2
    FALSE  $ lfc      : num -2.95
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1
    FALSE  $ direction: Factor w/ 3 levels "n5","NS","extend": 1
    FALSE NULL
    FALSE class: DESeqDataSet 
    FALSE dim: 13541 165 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13541): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(165): L.Blu13_male_pituitary_control.NYNO
    FALSE   L.G107_male_pituitary_control ...
    FALSE   y95.g131.x_male_pituitary_inc.d9 y98.o50.x_male_pituitary_inc.d3
    FALSE colData names(9): V1 bird ... treatment <- factor(treatment,
    FALSE   levels = alllevels) sextissue
    FALSE [1] 13541   165
    FALSE 'data.frame': 1577 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13541 levels "A2ML1","A2ML2",..: 6121 1232 7661 6607 10044 6380 7623 12537 5996 1219 ...
    FALSE  $ padj     : num  0.00221 0.01148 0.00437 0.02575 0.00693 ...
    FALSE  $ logpadj  : num  2.66 1.94 2.36 1.59 2.16 ...
    FALSE  $ lfc      : num  5.88 5.43 4.83 4.38 4.37 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "early","NS","inc.d9": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 4579 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13541 levels "A2ML1","A2ML2",..: 8016 11039 9221 13012 8331 11143 324 7369 7602 4565 ...
    FALSE  $ padj     : num  6.76e-06 3.16e-07 1.06e-02 1.43e-06 6.85e-06 ...
    FALSE  $ logpadj  : num  5.17 6.5 1.98 5.84 5.16 ...
    FALSE  $ lfc      : num  8.55 8.44 7.48 7.25 7.19 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","early": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 408 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13541 levels "A2ML1","A2ML2",..: 11039 8331 5934 6964 8536 4252 615 6632 6653 12314 ...
    FALSE  $ padj     : num  0.0021 0.018885 0.086142 0.000151 0.018885 ...
    FALSE  $ logpadj  : num  2.68 1.72 1.06 3.82 1.72 ...
    FALSE  $ lfc      : num  7.12 5.63 5.24 4.6 4.43 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "prolong","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 1154 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13541 levels "A2ML1","A2ML2",..: 12033 13071 13191 154 6892 8855 3427 2352 4123 9127 ...
    FALSE  $ padj     : num  0.0695 0.0627 0.0626 0.0545 0.0342 ...
    FALSE  $ logpadj  : num  1.16 1.2 1.2 1.26 1.47 ...
    FALSE  $ lfc      : num  2.29 1.33 1.19 1.16 1.13 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","prolong": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 1622 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13541 levels "A2ML1","A2ML2",..: 9205 10921 6063 11827 7893 6783 5273 10720 11453 6549 ...
    FALSE  $ padj     : num  0.084656 0.000111 0.0018 0.009349 0.005133 ...
    FALSE  $ logpadj  : num  1.07 3.95 2.74 2.03 2.29 ...
    FALSE  $ lfc      : num  6.12 5.07 3.46 3.3 3.2 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "extend","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE 'data.frame': 11 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13541 levels "A2ML1","A2ML2",..: 3609 10142 11461 11212 3163 1708 9847 1956 200 3162 ...
    FALSE  $ padj     : num  0.099 0.099 0.099 0.0568 0.099 ...
    FALSE  $ logpadj  : num  1 1 1 1.25 1 ...
    FALSE  $ lfc      : num  1.413 -0.516 -0.583 -0.706 -0.87 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "n5","NS","extend": 3 1 1 1 1 1 1 1 1 1 ...
    FALSE NULL

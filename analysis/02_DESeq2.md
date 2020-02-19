    library(tidyverse)

    ## ── Attaching packages ───────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ──────────────────────────── tidyverse_conflicts() ──
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
    library(caret) 

    ## Loading required package: lattice

    ## 
    ## Attaching package: 'caret'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     lift

    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    source("../R/functions.R")  # load custom functions 
    source("../R/themes.R")  # load custom themes and color palletes

    knitr::opts_chunk$set(fig.path = '../figures/DESeq2/',  message=F, comment=F, warning=F)

Warning: This script can take a very long time to run.
------------------------------------------------------

DESeq2 was not designed to run on 100+ samples. But, I really like it,
so I do it anyways. Some of these commands take like 15 min to run using
6 cores.

Characterization
----------------

    # import "colData" which contains sample information and "countData" which contains read counts
    countData <- read.csv("../results/00_countData_characterization.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)
    head(geneinfo)

    FALSE       Name geneid       entrezid
    FALSE 1    EDNRB 408082 NP_001001127.1
    FALSE 2  CYP26A1 408183 NP_001001129.1
    FALSE 3    CFDP1 374073 NP_001001189.1
    FALSE 4    AvBD7 407777 NP_001001194.1
    FALSE 5     KRT5 407779 NP_001001195.1
    FALSE 6 HSD11B1L 408034 NP_001001201.1

    # craete variable that will be critical for subset later on
    colData <- read.csv("../metadata/00_colData_characterization.csv", header = T, row.names = 1)
    colData$sextissue <- as.factor(paste(colData$sex, colData$tissue, sep = "_"))
    colData$treatment <- factor(colData$treatment, levels = charlevels)
    colData$tissue <- factor(colData$tissue, levels = tissuelevel)
    levels(colData$treatment)

    FALSE [1] "control" "bldg"    "lay"     "inc.d3"  "inc.d9"  "inc.d17" "hatch"  
    FALSE [8] "n5"      "n9"

    levels(colData$sex)

    FALSE [1] "female" "male"

    levels(colData$sextissue)

    FALSE [1] "female_gonad"        "female_hypothalamus" "female_pituitary"   
    FALSE [4] "male_gonad"          "male_hypothalamus"   "male_pituitary"

    levels(colData$tissue)

    FALSE [1] "hypothalamus" "pituitary"    "gonad"

Manipulation
------------

### Data wrangle

### Manipuation: comparisons to respective controls

    createDEGdfsaveManip <- function(up, down, mytissue){
      
      res <- results(dds, contrast = c("treatment", up, down), independentFiltering = T, alpha = 0.1)
      
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
      myfilename = paste0("../results/DESeq2/manip/", mytissue, partialfilename, "_DEGs.csv")
      
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
      
      myfilename = paste0("../results/DEseq2/manip/", i, "_vsd.csv")
      write.csv(vsd, myfilename)
      
      #return(dds)
      #return(vsd)
      #print(head(vsd))

      # save differential gene expression results
      inc.d3.m.inc.d3 <- createDEGdfsaveManip("m.inc.d3", "inc.d3",  i) 
      inc.d9.m.inc.d9 <- createDEGdfsaveManip("m.inc.d9", "inc.d9",  i) 
      inc.d17.m.inc.d17 <- createDEGdfsaveManip("m.inc.d17", "inc.d17",  i) 
      hatch.m.n2 <- createDEGdfsaveManip("m.n2", "hatch",  i) 
      inc.d9.m.inc.d8 <- createDEGdfsaveManip("m.inc.d8", "inc.d9",  i) 
      inc.d17.prolong <- createDEGdfsaveManip("prolong", "inc.d17",  i) 
      hatch.extend <- createDEGdfsaveManip("extend", "hatch",  i) 
    }

Manipulation: comparisons to other manipluations
------------------------------------------------

    # import "colData" which contains sample information and "countData" which contains read counts
    countData <- read.csv("../results/00_counts.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)
    head(geneinfo)


    levelsreplace <- c( "m.inc.d8" , "prolong" , "extend")
    levelsremoval <- c( "m.inc.d3" ,    "m.inc.d9" , "m.inc.d17" , "m.n2")
    controlsremovalreplace <- c( "inc.d3" ,    "inc.d9" , "inc.d17" , "hatch")

    manipulationsandcontrols <- c(controlsremoval, levelsreplace, levelsremoval)

    colData <- read.csv("../metadata/00_samples.csv", header = T, row.names = 1)
    colData$sextissue <- as.factor(paste(colData$sex, colData$tissue, sep = "_"))
    colData$treatment <- factor(colData$treatment, levels = manipulationsandcontrols)
    colData$tissue <- factor(colData$tissue, levels = tissuelevel)
    levels(colData$treatment)
    levels(colData$sex)
    levels(colData$sextissue)
    levels(colData$tissue)
    colData <- colData %>% drop_na()

    createDEGdfsaveManip <- function(up, down, mytissue){
      
      res <- results(dds, contrast = c("treatment", up, down), independentFiltering = T, alpha = 0.1)
      
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
      myfilename = paste0("../results/DESeq2/manip/", mytissue, partialfilename, "_DEGs.csv")
      
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
      
      myfilename = paste0("../results/DEseq2/manip/", i, "_vsd.csv")
      write.csv(vsd, myfilename)
      
      #return(dds)
      #return(vsd)
      #print(head(vsd))

      # save differential gene expression results
      hatch.m.n2 <- createDEGdfsaveManip("m.n2", "hatch",  i) 
      prolong.extend <- createDEGdfsaveManip("extend", "prolong",  i) 
      early.extend <- createDEGdfsaveManip("extend", "m.inc.d8",  i) 
      early.prolong <- createDEGdfsaveManip("prolong", "m.inc.d8",  i) 
      m.inc.d9.m.inc.d3 <- createDEGdfsaveManip("m.inc.d9", "m.inc.d3",  i) 
      m.inc.d17.m.inc.d9 <- createDEGdfsaveManip("m.inc.d17", "m.inc.d9",  i) 
      m.inc.d17.m.inc.d3 <- createDEGdfsaveManip("m.inc.d17", "m.inc.d3",  i) 
      m.n2.m.inc.d17 <- createDEGdfsaveManip("m.n2", "m.inc.d17",  i) 
      m.n2.m.inc.d9 <- createDEGdfsaveManip("m.n2", "m.inc.d9",  i)
      m.n2.m.inc.d3 <- createDEGdfsaveManip("m.n2", "m.inc.d3",  i)
      prolong.m.inc.d17 <- createDEGdfsaveManip("prolong", "m.inc.d17",  i) 
      m.inc.d8.m.inc.d9 <- createDEGdfsaveManip("m.inc.d8", "m.inc.d9",  i) 
      extend.m.n2 <- createDEGdfsaveManip("extend", "m.n2",  i) 
    }

PRL-driven DESeq2
-----------------

    PRLdata <- read_csv("../results/PRLvsd.csv")

    PRLsamples <- colData %>%
      mutate(samples = V1)  %>%
      select(samples) 
    PRLsamples

    FALSE                                          samples
    FALSE 1                L.Blu13_male_gonad_control.NYNO
    FALSE 2         L.Blu13_male_hypothalamus_control.NYNO
    FALSE 3            L.Blu13_male_pituitary_control.NYNO
    FALSE 4                      L.G107_male_gonad_control
    FALSE 5               L.G107_male_hypothalamus_control
    FALSE 6                  L.G107_male_pituitary_control
    FALSE 7                    L.G118_female_gonad_control
    FALSE 8        L.G118_female_hypothalamus_control.NYNO
    FALSE 9           L.G118_female_pituitary_control.NYNO
    FALSE 10                  L.R3_male_gonad_control.NYNO
    FALSE 11                L.R3_male_hypothalamus_control
    FALSE 12              L.R3_male_pituitary_control.NYNO
    FALSE 13                       L.R8_male_gonad_control
    FALSE 14                L.R8_male_hypothalamus_control
    FALSE 15                   L.R8_male_pituitary_control
    FALSE 16                      L.W33_male_gonad_control
    FALSE 17          L.W33_male_hypothalamus_control.NYNO
    FALSE 18                  L.W33_male_pituitary_control
    FALSE 19                  L.W3_male_gonad_control.NYNO
    FALSE 20                L.W3_male_hypothalamus_control
    FALSE 21                   L.W3_male_pituitary_control
    FALSE 22                  L.W4_male_gonad_control.NYNO
    FALSE 23                L.W4_male_hypothalamus_control
    FALSE 24                   L.W4_male_pituitary_control
    FALSE 25                   R.G106_female_gonad_control
    FALSE 26            R.G106_female_hypothalamus_control
    FALSE 27               R.G106_female_pituitary_control
    FALSE 28                    R.R20_female_gonad_control
    FALSE 29        R.R20_female_hypothalamus_control.NYNO
    FALSE 30                R.R20_female_pituitary_control
    FALSE 31                     R.R9_female_gonad_control
    FALSE 32              R.R9_female_hypothalamus_control
    FALSE 33            R.R9_female_pituitary_control.NYNO
    FALSE 34                    R.W44_female_gonad_control
    FALSE 35             R.W44_female_hypothalamus_control
    FALSE 36           R.W44_female_pituitary_control.NYNO
    FALSE 37                 R.Y108.W29_male_gonad_control
    FALSE 38     R.Y108.W29_male_hypothalamus_control.NYNO
    FALSE 39             R.Y108.W29_male_pituitary_control
    FALSE 40             blk.s061.pu.y_female_gonad_inc.d9
    FALSE 41      blk.s061.pu.y_female_hypothalamus_inc.d9
    FALSE 42         blk.s061.pu.y_female_pituitary_inc.d9
    FALSE 43                     blk11.x_female_gonad_bldg
    FALSE 44              blk11.x_female_hypothalamus_bldg
    FALSE 45                 blk11.x_female_pituitary_bldg
    FALSE 46                         blk12.x_male_gonad_n5
    FALSE 47             blk12.x_male_hypothalamus_n5.NYNO
    FALSE 48                     blk12.x_male_pituitary_n5
    FALSE 49                    blk17.x_male_gonad_inc.d17
    FALSE 50             blk17.x_male_hypothalamus_inc.d17
    FALSE 51                blk17.x_male_pituitary_inc.d17
    FALSE 52                    blk21.x_female_gonad_hatch
    FALSE 53             blk21.x_female_hypothalamus_hatch
    FALSE 54                blk21.x_female_pituitary_hatch
    FALSE 55                        blk4.x_female_gonad_n9
    FALSE 56                 blk4.x_female_hypothalamus_n9
    FALSE 57                    blk4.x_female_pituitary_n9
    FALSE 58            blu.o.x.ATLAS_female_gonad_control
    FALSE 59     blu.o.x.ATLAS_female_hypothalamus_control
    FALSE 60        blu.o.x.ATLAS_female_pituitary_control
    FALSE 61              blu103.x_female_gonad_hatch.NYNO
    FALSE 62            blu103.x_female_hypothalamus_hatch
    FALSE 63          blu103.x_female_pituitary_hatch.NYNO
    FALSE 64                blu104.w120.x_male_gonad_hatch
    FALSE 65         blu104.w120.x_male_hypothalamus_hatch
    FALSE 66       blu104.w120.x_male_pituitary_hatch.NYNO
    FALSE 67             blu108.w40.o158_male_gonad_inc.d9
    FALSE 68      blu108.w40.o158_male_hypothalamus_inc.d9
    FALSE 69         blu108.w40.o158_male_pituitary_inc.d9
    FALSE 70               blu111.w113.x_male_gonad_inc.d3
    FALSE 71        blu111.w113.x_male_hypothalamus_inc.d3
    FALSE 72           blu111.w113.x_male_pituitary_inc.d3
    FALSE 73              blu113.w124.x_male_gonad_inc.d17
    FALSE 74       blu113.w124.x_male_hypothalamus_inc.d17
    FALSE 75     blu113.w124.x_male_pituitary_inc.d17.NYNO
    FALSE 76               blu114.r38.w198_male_gonad_bldg
    FALSE 77        blu114.r38.w198_male_hypothalamus_bldg
    FALSE 78           blu114.r38.w198_male_pituitary_bldg
    FALSE 79               blu121.w91.x_male_gonad_inc.d17
    FALSE 80        blu121.w91.x_male_hypothalamus_inc.d17
    FALSE 81           blu121.w91.x_male_pituitary_inc.d17
    FALSE 82              blu124.w180.x_female_gonad_hatch
    FALSE 83       blu124.w180.x_female_hypothalamus_hatch
    FALSE 84          blu124.w180.x_female_pituitary_hatch
    FALSE 85                   blu33.y88.x_male_gonad_bldg
    FALSE 86            blu33.y88.x_male_hypothalamus_bldg
    FALSE 87               blu33.y88.x_male_pituitary_bldg
    FALSE 88                     blu36.w16_female_gonad_n9
    FALSE 89              blu36.w16_female_hypothalamus_n9
    FALSE 90                 blu36.w16_female_pituitary_n9
    FALSE 91                     blu37.r65.x_male_gonad_n5
    FALSE 92              blu37.r65.x_male_hypothalamus_n5
    FALSE 93                 blu37.r65.x_male_pituitary_n5
    FALSE 94                blu38.g135.x_female_gonad_bldg
    FALSE 95         blu38.g135.x_female_hypothalamus_bldg
    FALSE 96            blu38.g135.x_female_pituitary_bldg
    FALSE 97               blu39.o26.x_female_gonad_inc.d3
    FALSE 98   blu39.o26.x_female_hypothalamus_inc.d3.NYNO
    FALSE 99      blu39.o26.x_female_pituitary_inc.d3.NYNO
    FALSE 100                   blu41.y100.x_male_gonad_n5
    FALSE 101       blu41.y100.x_male_hypothalamus_n5.NYNO
    FALSE 102               blu41.y100.x_male_pituitary_n5
    FALSE 103              blu47.y96.x_female_gonad_inc.d9
    FALSE 104       blu47.y96.x_female_hypothalamus_inc.d9
    FALSE 105          blu47.y96.x_female_pituitary_inc.d9
    FALSE 106                    blu55.g51_female_gonad_n5
    FALSE 107             blu55.g51_female_hypothalamus_n5
    FALSE 108                blu55.g51_female_pituitary_n5
    FALSE 109                      blu81.r88_male_gonad_n9
    FALSE 110               blu81.r88_male_hypothalamus_n9
    FALSE 111                  blu81.r88_male_pituitary_n9
    FALSE 112                   d.s008.y.blk_male_gonad_n5
    FALSE 113            d.s008.y.blk_male_hypothalamus_n5
    FALSE 114               d.s008.y.blk_male_pituitary_n5
    FALSE 115                   d.s047.blk.o_male_gonad_n5
    FALSE 116            d.s047.blk.o_male_hypothalamus_n5
    FALSE 117               d.s047.blk.o_male_pituitary_n5
    FALSE 118               g.blk.s004.pk_female_gonad_lay
    FALSE 119        g.blk.s004.pk_female_hypothalamus_lay
    FALSE 120           g.blk.s004.pk_female_pituitary_lay
    FALSE 121                      g.s.blk.d_male_gonad_n9
    FALSE 122               g.s.blk.d_male_hypothalamus_n9
    FALSE 123                  g.s.blk.d_male_pituitary_n9
    FALSE 124                     g.s.blk.y_male_gonad_lay
    FALSE 125              g.s.blk.y_male_hypothalamus_lay
    FALSE 126                 g.s.blk.y_male_pituitary_lay
    FALSE 127                 g.s043.pu.blk_male_gonad_lay
    FALSE 128          g.s043.pu.blk_male_hypothalamus_lay
    FALSE 129             g.s043.pu.blk_male_pituitary_lay
    FALSE 130                g.s078.blk.o_female_gonad_lay
    FALSE 131         g.s078.blk.o_female_hypothalamus_lay
    FALSE 132            g.s078.blk.o_female_pituitary_lay
    FALSE 133               g.x.ATLAS_female_gonad_control
    FALSE 134                   g104.w82.x_male_gonad_bldg
    FALSE 135            g104.w82.x_male_hypothalamus_bldg
    FALSE 136               g104.w82.x_male_pituitary_bldg
    FALSE 137             g114.w83.x_male_gonad_hatch.NYNO
    FALSE 138           g114.w83.x_male_hypothalamus_hatch
    FALSE 139         g114.w83.x_male_pituitary_hatch.NYNO
    FALSE 140                g130.y81.x_male_gonad_inc.d17
    FALSE 141         g130.y81.x_male_hypothalamus_inc.d17
    FALSE 142            g130.y81.x_male_pituitary_inc.d17
    FALSE 143               g141.blu27.x_female_gonad_bldg
    FALSE 144        g141.blu27.x_female_hypothalamus_bldg
    FALSE 145           g141.blu27.x_female_pituitary_bldg
    FALSE 146              g142.r40.x_female_gonad_inc.d17
    FALSE 147       g142.r40.x_female_hypothalamus_inc.d17
    FALSE 148          g142.r40.x_female_pituitary_inc.d17
    FALSE 149              g143.blu32.x_male_gonad_inc.d17
    FALSE 150       g143.blu32.x_male_hypothalamus_inc.d17
    FALSE 151          g143.blu32.x_male_pituitary_inc.d17
    FALSE 152                 g146.blu51_male_gonad_inc.d3
    FALSE 153     g146.blu51_male_hypothalamus_inc.d3.NYNO
    FALSE 154             g146.blu51_male_pituitary_inc.d3
    FALSE 155                 g20.w106.x_male_gonad_inc.d3
    FALSE 156          g20.w106.x_male_hypothalamus_inc.d3
    FALSE 157             g20.w106.x_male_pituitary_inc.d3
    FALSE 158                    g52.blu58_male_gonad_bldg
    FALSE 159             g52.blu58_male_hypothalamus_bldg
    FALSE 160                g52.blu58_male_pituitary_bldg
    FALSE 161                     g53.y84_male_gonad_hatch
    FALSE 162              g53.y84_male_hypothalamus_hatch
    FALSE 163                 g53.y84_male_pituitary_hatch
    FALSE 164                g6.w197.x_female_gonad_inc.d3
    FALSE 165         g6.w197.x_female_hypothalamus_inc.d3
    FALSE 166            g6.w197.x_female_pituitary_inc.d3
    FALSE 167                    g75.x_female_gonad_inc.d9
    FALSE 168             g75.x_female_hypothalamus_inc.d9
    FALSE 169                g75.x_female_pituitary_inc.d9
    FALSE 170               l.s120.y.blk_female_gonad_bldg
    FALSE 171        l.s120.y.blk_female_hypothalamus_bldg
    FALSE 172           l.s120.y.blk_female_pituitary_bldg
    FALSE 173                       o.s.w.r_male_gonad_lay
    FALSE 174                o.s.w.r_male_hypothalamus_lay
    FALSE 175                   o.s.w.r_male_pituitary_lay
    FALSE 176                  o152.o120.w42_male_gonad_n5
    FALSE 177           o152.o120.w42_male_hypothalamus_n5
    FALSE 178              o152.o120.w42_male_pituitary_n5
    FALSE 179               o156.w80.x_female_gonad_inc.d3
    FALSE 180        o156.w80.x_female_hypothalamus_inc.d3
    FALSE 181           o156.w80.x_female_pituitary_inc.d3
    FALSE 182         o165.w122.x_female_gonad_inc.d3.NYNO
    FALSE 183       o165.w122.x_female_hypothalamus_inc.d3
    FALSE 184     o165.w122.x_female_pituitary_inc.d3.NYNO
    FALSE 185               o172.w115.x_female_gonad_hatch
    FALSE 186        o172.w115.x_female_hypothalamus_hatch
    FALSE 187      o172.w115.x_female_pituitary_hatch.NYNO
    FALSE 188              o173.w179.x_female_gonad_inc.d3
    FALSE 189       o173.w179.x_female_hypothalamus_inc.d3
    FALSE 190          o173.w179.x_female_pituitary_inc.d3
    FALSE 191               o35.r51.x_female_gonad_inc.d17
    FALSE 192        o35.r51.x_female_hypothalamus_inc.d17
    FALSE 193           o35.r51.x_female_pituitary_inc.d17
    FALSE 194                o38.blu29.x_female_gonad_bldg
    FALSE 195         o38.blu29.x_female_hypothalamus_bldg
    FALSE 196            o38.blu29.x_female_pituitary_bldg
    FALSE 197                   o39.y77.x_male_gonad_hatch
    FALSE 198            o39.y77.x_male_hypothalamus_hatch
    FALSE 199               o39.y77.x_male_pituitary_hatch
    FALSE 200                 o44.blu26.x_male_gonad_hatch
    FALSE 201          o44.blu26.x_male_hypothalamus_hatch
    FALSE 202             o44.blu26.x_male_pituitary_hatch
    FALSE 203                 o48.r197.x_male_gonad_inc.d3
    FALSE 204          o48.r197.x_male_hypothalamus_inc.d3
    FALSE 205             o48.r197.x_male_pituitary_inc.d3
    FALSE 206                      o49.x_male_gonad_inc.d9
    FALSE 207               o49.x_male_hypothalamus_inc.d9
    FALSE 208                  o49.x_male_pituitary_inc.d9
    FALSE 209               o52.blu53_female_gonad_inc.d17
    FALSE 210        o52.blu53_female_hypothalamus_inc.d17
    FALSE 211           o52.blu53_female_pituitary_inc.d17
    FALSE 212                    o57.g59_male_gonad_inc.d9
    FALSE 213             o57.g59_male_hypothalamus_inc.d9
    FALSE 214                o57.g59_male_pituitary_inc.d9
    FALSE 215                    o73.x_female_gonad_inc.d9
    FALSE 216             o73.x_female_hypothalamus_inc.d9
    FALSE 217                o73.x_female_pituitary_inc.d9
    FALSE 218                 pk.s238.blk.w_male_gonad_lay
    FALSE 219          pk.s238.blk.w_male_hypothalamus_lay
    FALSE 220             pk.s238.blk.w_male_pituitary_lay
    FALSE 221                   pk.w.s141.o_male_gonad_lay
    FALSE 222            pk.w.s141.o_male_hypothalamus_lay
    FALSE 223               pk.w.s141.o_male_pituitary_lay
    FALSE 224        r.r.x.ATLAS.R2XR_female_gonad_control
    FALSE 225 r.r.x.ATLAS.R2XR_female_hypothalamus_control
    FALSE 226    r.r.x.ATLAS.R2XR_female_pituitary_control
    FALSE 227             r.r.x.ATLAS_female_gonad_control
    FALSE 228      r.r.x.ATLAS_female_hypothalamus_control
    FALSE 229         r.r.x.ATLAS_female_pituitary_control
    FALSE 230                 r.s005.pk.blk_male_gonad_lay
    FALSE 231          r.s005.pk.blk_male_hypothalamus_lay
    FALSE 232             r.s005.pk.blk_male_pituitary_lay
    FALSE 233                 r.s056.g.o_female_gonad_bldg
    FALSE 234          r.s056.g.o_female_hypothalamus_bldg
    FALSE 235             r.s056.g.o_female_pituitary_bldg
    FALSE 236                   r.s059.d.o_male_gonad_bldg
    FALSE 237            r.s059.d.o_male_hypothalamus_bldg
    FALSE 238               r.s059.d.o_male_pituitary_bldg
    FALSE 239                 r.s116.blk.pu_male_gonad_lay
    FALSE 240          r.s116.blk.pu_male_hypothalamus_lay
    FALSE 241             r.s116.blk.pu_male_pituitary_lay
    FALSE 242                   r.s171.l.w_female_gonad_n9
    FALSE 243            r.s171.l.w_female_hypothalamus_n9
    FALSE 244               r.s171.l.w_female_pituitary_n9
    FALSE 245                   r.y.s007.blk_male_gonad_n9
    FALSE 246            r.y.s007.blk_male_hypothalamus_n9
    FALSE 247               r.y.s007.blk_male_pituitary_n9
    FALSE 248                r176.blu54_male_gonad_inc.d17
    FALSE 249         r176.blu54_male_hypothalamus_inc.d17
    FALSE 250            r176.blu54_male_pituitary_inc.d17
    FALSE 251                  r183.o22_female_gonad_hatch
    FALSE 252           r183.o22_female_hypothalamus_hatch
    FALSE 253              r183.o22_female_pituitary_hatch
    FALSE 254                    r190.o43.x_male_gonad_lay
    FALSE 255             r190.o43.x_male_hypothalamus_lay
    FALSE 256                r190.o43.x_male_pituitary_lay
    FALSE 257                         r195.x_male_gonad_n9
    FALSE 258                  r195.x_male_hypothalamus_n9
    FALSE 259                     r195.x_male_pituitary_n9
    FALSE 260          r27.w111.blu125_female_gonad_inc.d3
    FALSE 261   r27.w111.blu125_female_hypothalamus_inc.d3
    FALSE 262      r27.w111.blu125_female_pituitary_inc.d3
    FALSE 263             r30.w112.r46_female_gonad_inc.d9
    FALSE 264      r30.w112.r46_female_hypothalamus_inc.d9
    FALSE 265         r30.w112.r46_female_pituitary_inc.d9
    FALSE 266               r36.w184.x_female_gonad_inc.d9
    FALSE 267        r36.w184.x_female_hypothalamus_inc.d9
    FALSE 268           r36.w184.x_female_pituitary_inc.d9
    FALSE 269                     r37.w100.x_male_gonad_n9
    FALSE 270              r37.w100.x_male_hypothalamus_n9
    FALSE 271                 r37.w100.x_male_pituitary_n9
    FALSE 272                   r41.w99.x_male_gonad_hatch
    FALSE 273            r41.w99.x_male_hypothalamus_hatch
    FALSE 274               r41.w99.x_male_pituitary_hatch
    FALSE 275                      r45.X_male_gonad_inc.d9
    FALSE 276                  r45.X_male_pituitary_inc.d9
    FALSE 277               r45.x_male_hypothalamus_inc.d9
    FALSE 278              r49.w189.x_female_gonad_inc.d17
    FALSE 279       r49.w189.x_female_hypothalamus_inc.d17
    FALSE 280          r49.w189.x_female_pituitary_inc.d17
    FALSE 281               r6.x_female_gonad_control.NYNO
    FALSE 282        r6.x_female_hypothalamus_control.NYNO
    FALSE 283                r6.x_female_pituitary_control
    FALSE 284                   r72.y83.x_male_gonad_hatch
    FALSE 285            r72.y83.x_male_hypothalamus_hatch
    FALSE 286               r72.y83.x_male_pituitary_hatch
    FALSE 287               r73.g127.x_female_gonad_inc.d3
    FALSE 288        r73.g127.x_female_hypothalamus_inc.d3
    FALSE 289           r73.g127.x_female_pituitary_inc.d3
    FALSE 290                    r83.g45_female_gonad_bldg
    FALSE 291             r83.g45_female_hypothalamus_bldg
    FALSE 292                r83.g45_female_pituitary_bldg
    FALSE 293                    r95.blu99_female_gonad_n9
    FALSE 294             r95.blu99_female_hypothalamus_n9
    FALSE 295                r95.blu99_female_pituitary_n9
    FALSE 296                      s.o.pk_female_gonad_lay
    FALSE 297               s.o.pk_female_hypothalamus_lay
    FALSE 298                  s.o.pk_female_pituitary_lay
    FALSE 299                s.pu148.blk.r_male_gonad_bldg
    FALSE 300         s.pu148.blk.r_male_hypothalamus_bldg
    FALSE 301            s.pu148.blk.r_male_pituitary_bldg
    FALSE 302               s.x.ATLAS_female_gonad_control
    FALSE 303        s.x.ATLAS_female_hypothalamus_control
    FALSE 304           s.x.ATLAS_female_pituitary_control
    FALSE 305               s063.d.blk.l_female_gonad_bldg
    FALSE 306        s063.d.blk.l_female_hypothalamus_bldg
    FALSE 307           s063.d.blk.l_female_pituitary_bldg
    FALSE 308                   s065.l.d.o_male_gonad_bldg
    FALSE 309            s065.l.d.o_male_hypothalamus_bldg
    FALSE 310               s065.l.d.o_male_pituitary_bldg
    FALSE 311                   s066.l.d.r_male_gonad_bldg
    FALSE 312            s066.l.d.r_male_hypothalamus_bldg
    FALSE 313               s066.l.d.r_male_pituitary_bldg
    FALSE 314               s092.blk.r.o_female_gonad_bldg
    FALSE 315        s092.blk.r.o_female_hypothalamus_bldg
    FALSE 316           s092.blk.r.o_female_pituitary_bldg
    FALSE 317                s095.g.blk.o_female_gonad_lay
    FALSE 318         s095.g.blk.o_female_hypothalamus_lay
    FALSE 319            s095.g.blk.o_female_pituitary_lay
    FALSE 320                  s136.d.w.o_female_gonad_lay
    FALSE 321           s136.d.w.o_female_hypothalamus_lay
    FALSE 322              s136.d.w.o_female_pituitary_lay
    FALSE 323                s142.o.pk.pu_female_gonad_lay
    FALSE 324         s142.o.pk.pu_female_hypothalamus_lay
    FALSE 325            s142.o.pk.pu_female_pituitary_lay
    FALSE 326                  s150.w.g.blk_male_gonad_lay
    FALSE 327           s150.w.g.blk_male_hypothalamus_lay
    FALSE 328              s150.w.g.blk_male_pituitary_lay
    FALSE 329               s176.blk.pu.r_female_gonad_lay
    FALSE 330        s176.blk.pu.r_female_hypothalamus_lay
    FALSE 331           s176.blk.pu.r_female_pituitary_lay
    FALSE 332                     s187.l.o.r_male_gonad_n9
    FALSE 333              s187.l.o.r_male_hypothalamus_n9
    FALSE 334                 s187.l.o.r_male_pituitary_n9
    FALSE 335                 s243.blk.pk.r_male_gonad_lay
    FALSE 336          s243.blk.pk.r_male_hypothalamus_lay
    FALSE 337             s243.blk.pk.r_male_pituitary_lay
    FALSE 338                 w191.r1_female_gonad_control
    FALSE 339          w191.r1_female_hypothalamus_control
    FALSE 340             w191.r1_female_pituitary_control
    FALSE 341                      w34.x_male_gonad_inc.d9
    FALSE 342               w34.x_male_hypothalamus_inc.d9
    FALSE 343                  w34.x_male_pituitary_inc.d9
    FALSE 344           x.blk.blk.ATLAS_male_gonad_control
    FALSE 345    x.blk.blk.ATLAS_male_hypothalamus_control
    FALSE 346       x.blk.blk.ATLAS_male_pituitary_control
    FALSE 347                   x.blk16_male_gonad_n9.NYNO
    FALSE 348            x.blk16_male_hypothalamus_n9.NYNO
    FALSE 349                    x.blk16_male_pituitary_n9
    FALSE 350         x.blu.o.ATLAS_male_pituitary_control
    FALSE 351             x.blu101.w43_female_gonad_inc.d9
    FALSE 352      x.blu101.w43_female_hypothalamus_inc.d9
    FALSE 353         x.blu101.w43_female_pituitary_inc.d9
    FALSE 354            x.blu102.w105_female_gonad_inc.d3
    FALSE 355     x.blu102.w105_female_hypothalamus_inc.d3
    FALSE 356        x.blu102.w105_female_pituitary_inc.d3
    FALSE 357         x.blu106.o153_male_gonad_inc.d9.NYNO
    FALSE 358       x.blu106.o153_male_hypothalamus_inc.d9
    FALSE 359     x.blu106.o153_male_pituitary_inc.d9.NYNO
    FALSE 360                x.blu109.w121_female_gonad_n5
    FALSE 361         x.blu109.w121_female_hypothalamus_n5
    FALSE 362            x.blu109.w121_female_pituitary_n5
    FALSE 363           x.blu116.w107_female_gonad_inc.d17
    FALSE 364    x.blu116.w107_female_hypothalamus_inc.d17
    FALSE 365  x.blu116.w107_female_pituitary_inc.d17.NYNO
    FALSE 366              x.blu117.w89_male_gonad_inc.d17
    FALSE 367       x.blu117.w89_male_hypothalamus_inc.d17
    FALSE 368          x.blu117.w89_male_pituitary_inc.d17
    FALSE 369             x.blu122.r66_female_gonad_inc.d9
    FALSE 370      x.blu122.r66_female_hypothalamus_inc.d9
    FALSE 371         x.blu122.r66_female_pituitary_inc.d9
    FALSE 372                    x.blu23.w14_male_gonad_n9
    FALSE 373             x.blu23.w14_male_hypothalamus_n9
    FALSE 374                x.blu23.w14_male_pituitary_n9
    FALSE 375                        x.blu30_male_gonad_n5
    FALSE 376                 x.blu30_male_hypothalamus_n5
    FALSE 377                    x.blu30_male_pituitary_n5
    FALSE 378                x.blu42.o28_male_gonad_inc.d3
    FALSE 379    x.blu42.o28_male_hypothalamus_inc.d3.NYNO
    FALSE 380            x.blu42.o28_male_pituitary_inc.d3
    FALSE 381                 x.blu43.g132_female_gonad_n9
    FALSE 382          x.blu43.g132_female_hypothalamus_n9
    FALSE 383             x.blu43.g132_female_pituitary_n9
    FALSE 384                  x.blu6.y80_female_gonad_lay
    FALSE 385           x.blu6.y80_female_hypothalamus_lay
    FALSE 386              x.blu6.y80_female_pituitary_lay
    FALSE 387                 x.g.ATLAS_male_gonad_control
    FALSE 388          x.g.ATLAS_male_hypothalamus_control
    FALSE 389             x.g.ATLAS_male_pituitary_control
    FALSE 390             x.g.g.ATLAS_female_gonad_control
    FALSE 391               x.g.g.ATLAS_male_gonad_control
    FALSE 392        x.g.g.ATLAS_male_hypothalamus_control
    FALSE 393           x.g.g.ATLAS_male_pituitary_control
    FALSE 394             x.g.g.g.ATLAS_male_gonad_control
    FALSE 395         x.g.g.g.ATLAS_male_pituitary_control
    FALSE 396                 x.g13.w109_male_gonad_inc.d9
    FALSE 397          x.g13.w109_male_hypothalamus_inc.d9
    FALSE 398             x.g13.w109_male_pituitary_inc.d9
    FALSE 399                x.g14.w199_male_gonad_inc.d17
    FALSE 400         x.g14.w199_male_hypothalamus_inc.d17
    FALSE 401            x.g14.w199_male_pituitary_inc.d17
    FALSE 402               x.g147.blu28_male_gonad_inc.d3
    FALSE 403        x.g147.blu28_male_hypothalamus_inc.d3
    FALSE 404           x.g147.blu28_male_pituitary_inc.d3
    FALSE 405                        x.g37_female_gonad_n5
    FALSE 406                 x.g37_female_hypothalamus_n5
    FALSE 407                    x.g37_female_pituitary_n5
    FALSE 408                     x.g4.w50_female_gonad_n9
    FALSE 409              x.g4.w50_female_hypothalamus_n9
    FALSE 410                 x.g4.w50_female_pituitary_n9
    FALSE 411                        x.g43_female_gonad_n5
    FALSE 412                 x.g43_female_hypothalamus_n5
    FALSE 413                    x.g43_female_pituitary_n5
    FALSE 414                        x.g49_female_gonad_n5
    FALSE 415                 x.g49_female_hypothalamus_n5
    FALSE 416               x.g49_female_pituitary_n5.NYNO
    FALSE 417                       x.g70_male_gonad_hatch
    FALSE 418                x.g70_male_hypothalamus_hatch
    FALSE 419                   x.g70_male_pituitary_hatch
    FALSE 420           x.g9.o166_female_gonad_inc.d9.NYNO
    FALSE 421         x.g9.o166_female_hypothalamus_inc.d9
    FALSE 422       x.g9.o166_female_pituitary_inc.d9.NYNO
    FALSE 423              x.o159.w90_female_gonad_inc.d17
    FALSE 424       x.o159.w90_female_hypothalamus_inc.d17
    FALSE 425          x.o159.w90_female_pituitary_inc.d17
    FALSE 426                 x.o160.w102_male_gonad_hatch
    FALSE 427          x.o160.w102_male_hypothalamus_hatch
    FALSE 428             x.o160.w102_male_pituitary_hatch
    FALSE 429                x.o163.w101_male_gonad_inc.d3
    FALSE 430         x.o163.w101_male_hypothalamus_inc.d3
    FALSE 431       x.o163.w101_male_pituitary_inc.d3.NYNO
    FALSE 432                    x.o164.w123_male_gonad_n5
    FALSE 433             x.o164.w123_male_hypothalamus_n5
    FALSE 434           x.o164.w123_male_pituitary_n5.NYNO
    FALSE 435                   x.o175.g21_female_gonad_n5
    FALSE 436            x.o175.g21_female_hypothalamus_n5
    FALSE 437               x.o175.g21_female_pituitary_n5
    FALSE 438                           x.o2_male_gonad_n9
    FALSE 439                    x.o2_male_hypothalamus_n9
    FALSE 440                       x.o2_male_pituitary_n9
    FALSE 441                   x.o30.g134_male_gonad_bldg
    FALSE 442            x.o30.g134_male_hypothalamus_bldg
    FALSE 443               x.o30.g134_male_pituitary_bldg
    FALSE 444               x.o37.blu50_female_gonad_hatch
    FALSE 445   x.o37.blu50_female_hypothalamus_hatch.NYNO
    FALSE 446           x.o37.blu50_female_pituitary_hatch
    FALSE 447                  x.o47.y82_male_gonad_inc.d9
    FALSE 448           x.o47.y82_male_hypothalamus_inc.d9
    FALSE 449              x.o47.y82_male_pituitary_inc.d9
    FALSE 450                          x.o68_male_gonad_n5
    FALSE 451                   x.o68_male_hypothalamus_n5
    FALSE 452                      x.o68_male_pituitary_n5
    FALSE 453                        x.o70_female_gonad_n5
    FALSE 454            x.o70_female_hypothalamus_n5.NYNO
    FALSE 455                    x.o70_female_pituitary_n5
    FALSE 456                      x.r178_male_gonad_hatch
    FALSE 457               x.r178_male_hypothalamus_hatch
    FALSE 458                  x.r178_male_pituitary_hatch
    FALSE 459                         x.r181_male_gonad_n5
    FALSE 460                  x.r181_male_hypothalamus_n5
    FALSE 461                     x.r181_male_pituitary_n5
    FALSE 462                 x.r29.w96_male_gonad_inc.d17
    FALSE 463          x.r29.w96_male_hypothalamus_inc.d17
    FALSE 464             x.r29.w96_male_pituitary_inc.d17
    FALSE 465               x.r33.w183_female_gonad_inc.d3
    FALSE 466        x.r33.w183_female_hypothalamus_inc.d3
    FALSE 467           x.r33.w183_female_pituitary_inc.d3
    FALSE 468                  x.r39.g10_female_gonad_bldg
    FALSE 469           x.r39.g10_female_hypothalamus_bldg
    FALSE 470              x.r39.g10_female_pituitary_bldg
    FALSE 471                 x.r44.w95_female_gonad_hatch
    FALSE 472          x.r44.w95_female_hypothalamus_hatch
    FALSE 473             x.r44.w95_female_pituitary_hatch
    FALSE 474              x.r48.y139_female_gonad_inc.d17
    FALSE 475       x.r48.y139_female_hypothalamus_inc.d17
    FALSE 476     x.r48.y139_female_pituitary_inc.d17.NYNO
    FALSE 477                    x.r50.w97_female_gonad_n5
    FALSE 478             x.r50.w97_female_hypothalamus_n5
    FALSE 479                x.r50.w97_female_pituitary_n5
    FALSE 480                 x.r64.g140_male_gonad_inc.d3
    FALSE 481          x.r64.g140_male_hypothalamus_inc.d3
    FALSE 482             x.r64.g140_male_pituitary_inc.d3
    FALSE 483                  x.r67.blu35_male_gonad_bldg
    FALSE 484           x.r67.blu35_male_hypothalamus_bldg
    FALSE 485         x.r67.blu35_male_pituitary_bldg.NYNO
    FALSE 486                       x.w178_female_gonad_n9
    FALSE 487                x.w178_female_hypothalamus_n9
    FALSE 488                   x.w178_female_pituitary_n9
    FALSE 489                x.w192.o157_male_gonad_inc.d9
    FALSE 490         x.w192.o157_male_hypothalamus_inc.d9
    FALSE 491            x.w192.o157_male_pituitary_inc.d9
    FALSE 492                       x.w51_female_gonad_lay
    FALSE 493                x.w51_female_hypothalamus_lay
    FALSE 494                   x.w51_female_pituitary_lay
    FALSE 495                         x.w6_female_gonad_n9
    FALSE 496                  x.w6_female_hypothalamus_n9
    FALSE 497                     x.w6_female_pituitary_n9
    FALSE 498               x.y.s.ATLAS_male_gonad_control
    FALSE 499           x.y.s.ATLAS_male_pituitary_control
    FALSE 500                   x.y109_female_gonad_inc.d9
    FALSE 501            x.y109_female_hypothalamus_inc.d9
    FALSE 502               x.y109_female_pituitary_inc.d9
    FALSE 503                x.y132.w76_male_gonad_inc.d17
    FALSE 504         x.y132.w76_male_hypothalamus_inc.d17
    FALSE 505            x.y132.w76_male_pituitary_inc.d17
    FALSE 506                  x.y138.w176_female_gonad_n9
    FALSE 507           x.y138.w176_female_hypothalamus_n9
    FALSE 508              x.y138.w176_female_pituitary_n9
    FALSE 509                x.y141.w116_male_gonad_inc.d9
    FALSE 510         x.y141.w116_male_hypothalamus_inc.d9
    FALSE 511            x.y141.w116_male_pituitary_inc.d9
    FALSE 512                     x.y90_female_gonad_hatch
    FALSE 513              x.y90_female_hypothalamus_hatch
    FALSE 514                 x.y90_female_pituitary_hatch
    FALSE 515               x.y93.g126_female_gonad_inc.d9
    FALSE 516        x.y93.g126_female_hypothalamus_inc.d9
    FALSE 517           x.y93.g126_female_pituitary_inc.d9
    FALSE 518                         x.y9_female_gonad_n9
    FALSE 519                  x.y9_female_hypothalamus_n9
    FALSE 520                     x.y9_female_pituitary_n9
    FALSE 521                  y.s156.o.r_female_gonad_lay
    FALSE 522           y.s156.o.r_female_hypothalamus_lay
    FALSE 523              y.s156.o.r_female_pituitary_lay
    FALSE 524              y126.w92.x_female_gonad_inc.d17
    FALSE 525       y126.w92.x_female_hypothalamus_inc.d17
    FALSE 526          y126.w92.x_female_pituitary_inc.d17
    FALSE 527               y128.g23.x_female_gonad_inc.d9
    FALSE 528           y128.g23.x_female_pituitary_inc.d9
    FALSE 529                         y129.x_male_gonad_n9
    FALSE 530                  y129.x_male_hypothalamus_n9
    FALSE 531                     y129.x_male_pituitary_n9
    FALSE 532                    y13.x_female_gonad_inc.d3
    FALSE 533             y13.x_female_hypothalamus_inc.d3
    FALSE 534                y13.x_female_pituitary_inc.d3
    FALSE 535             y130.o170.x_female_gonad_inc.d17
    FALSE 536      y130.o170.x_female_hypothalamus_inc.d17
    FALSE 537         y130.o170.x_female_pituitary_inc.d17
    FALSE 538                    y131.w185.x_male_gonad_n9
    FALSE 539             y131.w185.x_male_hypothalamus_n9
    FALSE 540                y131.w185.x_male_pituitary_n9
    FALSE 541              y133.w77.r58_male_gonad_inc.d17
    FALSE 542       y133.w77.r58_male_hypothalamus_inc.d17
    FALSE 543          y133.w77.r58_male_pituitary_inc.d17
    FALSE 544           y135.blu107.x_female_gonad_inc.d17
    FALSE 545    y135.blu107.x_female_hypothalamus_inc.d17
    FALSE 546  y135.blu107.x_female_pituitary_inc.d17.NYNO
    FALSE 547                  y136.x_female_gonad_inc.d17
    FALSE 548           y136.x_female_hypothalamus_inc.d17
    FALSE 549              y136.x_female_pituitary_inc.d17
    FALSE 550              y140.w119.x_female_gonad_inc.d9
    FALSE 551       y140.w119.x_female_hypothalamus_inc.d9
    FALSE 552          y140.w119.x_female_pituitary_inc.d9
    FALSE 553                 y149.r52.x_male_gonad_inc.d3
    FALSE 554          y149.r52.x_male_hypothalamus_inc.d3
    FALSE 555             y149.r52.x_male_pituitary_inc.d3
    FALSE 556                     y15.x_female_gonad_hatch
    FALSE 557              y15.x_female_hypothalamus_hatch
    FALSE 558                 y15.x_female_pituitary_hatch
    FALSE 559                       y6.o54_female_gonad_n5
    FALSE 560                y6.o54_female_hypothalamus_n5
    FALSE 561                   y6.o54_female_pituitary_n5
    FALSE 562                    y7.g58_female_gonad_hatch
    FALSE 563             y7.g58_female_hypothalamus_hatch
    FALSE 564                y7.g58_female_pituitary_hatch
    FALSE 565                   y94.g133.x_female_gonad_n5
    FALSE 566       y94.g133.x_female_hypothalamus_n5.NYNO
    FALSE 567               y94.g133.x_female_pituitary_n5
    FALSE 568                 y95.g131.x_male_gonad_inc.d9
    FALSE 569          y95.g131.x_male_hypothalamus_inc.d9
    FALSE 570             y95.g131.x_male_pituitary_inc.d9
    FALSE 571                        y97.x_female_gonad_n9
    FALSE 572                 y97.x_female_hypothalamus_n9
    FALSE 573                    y97.x_female_pituitary_n9
    FALSE 574                  y98.o50.x_male_gonad_inc.d3
    FALSE 575           y98.o50.x_male_hypothalamus_inc.d3
    FALSE 576              y98.o50.x_male_pituitary_inc.d3

    # join with hi lo, and drop controls
    colDataPRL <- full_join(PRLdata, PRLsamples) %>% 
      drop_na() %>%
      select(samples, tissue, sex, treatment, hiloPRL) %>% 
      filter(!treatment %in% c("control", "bldg"))

    # sample 100 only, better for deseq2
    colDataPRL <- colDataPRL[sample(1:nrow(colDataPRL), 100,
       replace=FALSE),]

    colDataPRL <- as.data.frame(colDataPRL)
    row.names(colDataPRL) <- colDataPRL$samples

    savecols <- as.character(colDataPRL$samples) 
    savecols <- as.vector(savecols) 
      
    newcountData <- countData %>% dplyr::select(one_of(savecols)) 

    colDataPRL %>%
      group_by(sex, tissue, hiloPRL) %>%
      summarize(n = n())

    FALSE # A tibble: 4 x 4
    FALSE # Groups:   sex, tissue [2]
    FALSE   sex    tissue    hiloPRL     n
    FALSE   <chr>  <chr>     <chr>   <int>
    FALSE 1 female pituitary hi         25
    FALSE 2 female pituitary lo         26
    FALSE 3 male   pituitary hi         27
    FALSE 4 male   pituitary lo         22

    dds <- DESeqDataSetFromMatrix(countData = newcountData,
                                    colData = colDataPRL,
                                    design = ~sex * hiloPRL )
    dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
    print(dds)

    FALSE class: DESeqDataSet 
    FALSE dim: 13419 100 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13419): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(100): y94.g133.x_female_pituitary_n5
    FALSE   g130.y81.x_male_pituitary_inc.d17 ... x.y9_female_pituitary_n9
    FALSE   g53.y84_male_pituitary_hatch
    FALSE colData names(5): samples tissue sex treatment hiloPRL

    print(dim(dds))

    FALSE [1] 13419   100

    dds <- DESeq(dds, parallel = TRUE) # Differential expression analysis
      
    vsd <- as.data.frame(assay(vst(dds, blind=FALSE)))

    ressex <- results( dds, contrast = c("sex", "male", "female") )
    reshiloPRL <- results( dds, contrast = c("hiloPRL", "hi", "lo") )
    hist( ressex$pvalue, breaks=20, col="grey" )

![](../figures/DESeq2/PRL-1.png)

    hist( reshiloPRL$pvalue, breaks=20, col="grey" )

![](../figures/DESeq2/PRL-2.png)

    restoDEGs <- function(res, up, down){
      DEGs <- data.frame(gene = row.names(res),
                            padj = res$padj, 
                            logpadj = -log10(res$padj),
                            lfc = res$log2FoldChange)
      DEGs <- na.omit(DEGs)
      DEGs <- DEGs %>%
        dplyr::mutate(direction = ifelse(DEGs$lfc > 0 & DEGs$padj < 0.1, 
                                         yes = up, no = ifelse(DEGs$lfc < 0 & DEGs$padj < 0.1, 
                                                               yes = down, no = "NS"))) %>% 
        dplyr::arrange(gene)
      DEGs$direction <- factor(DEGs$direction, levels = c(down, "NS", up)) 
      DEGs <- DEGs %>% dplyr::filter(direction != "NS")
      print(head(DEGs))
      return(DEGs)
    }  

    sexdegs <- restoDEGs(ressex, "male", "female")

    FALSE      gene         padj   logpadj        lfc direction
    FALSE 1   A2ML4 1.801089e-11 10.744465  2.4625893      male
    FALSE 2   AAED1 4.650209e-10  9.332528  0.5339629      male
    FALSE 3    ABAT 4.765668e-02  1.321876  0.1468298      male
    FALSE 4   ABCA1 3.204461e-29 28.494245  0.8795504      male
    FALSE 5 ABCB1LB 1.260771e-02  1.899364 -0.1937905    female
    FALSE 6   ABCC5 9.741251e-02  1.011385  0.1447697      male

    PRLdegs <- restoDEGs(reshiloPRL, "hi", "lo")

    FALSE      gene         padj  logpadj        lfc direction
    FALSE 1   ABCA1 1.137419e-02 1.944080 -0.2542339        lo
    FALSE 2   ABCA2 2.146421e-02 1.668285 -0.3544446        lo
    FALSE 3 ABCB1LB 6.367575e-02 1.196026  0.1528690        hi
    FALSE 4   ABCC3 3.714724e-03 2.430073  0.3251622        hi
    FALSE 5   ABCC8 1.148503e-10 9.939868  1.0214231        hi
    FALSE 6   ABCC9 7.174456e-08 7.144211  0.8173118        hi

    plotbarvolcano <- function(df, mylabels, mycolorname){
      p1 <- ggplot(df, aes(x = lfc, y = logpadj, color = direction)) +
        geom_point() +
        scale_color_manual(values = allcolors, name = mycolorname) +
        theme_B3() +
        theme(legend.position = "top") +
        labs(x = "Log fold change", y = "Log adj. p-value") +
        xlim(-5,5)

      p2 <- ggplot(df, aes(x = direction,  fill = direction)) +
        geom_bar() +
        theme_B3() +
        scale_fill_manual(values = allcolors) +
        theme_B3() +
        theme(legend.position = "none") +
        labs(x = mycolorname, y = "Total DEGs")
      
      p <- plot_grid(p1,p2, labels  = mylabels, label_size = 12)
      plot(p)
    }

    top <- plotbarvolcano(sexdegs, c("a","b"), "sex")

![](../figures/DESeq2/PRL-3.png)

    bottom <- plotbarvolcano(PRLdegs, c("c","d"), "PRL")

![](../figures/DESeq2/PRL-4.png)

    ab <- plot_grid(top,bottom, nrow = 2)

    sexPRLdesg <- full_join(sexdegs, PRLdegs, by = "gene")

    sexPRLdesg$direction.x <- as.character(sexPRLdesg$direction.x)

    sexPRLdesg <- full_join(sexdegs, PRLdegs, by = "gene") %>%
      dplyr::select(gene, direction.x, direction.y) %>%
      tidyr::replace_na(list(direction.x = "no sex effect", direction.y = "no PRL effect"))
    sexPRLdesg

    FALSE                gene direction.x direction.y
    FALSE 1             A2ML4        male        <NA>
    FALSE 2             AAED1        male        <NA>
    FALSE 3              ABAT        male        <NA>
    FALSE 4             ABCA1        male          lo
    FALSE 5           ABCB1LB      female          hi
    FALSE 6             ABCC5        male        <NA>
    FALSE 7             ABCG1        male        <NA>
    FALSE 8             ABCG2      female        <NA>
    FALSE 9             ABCG4        male          lo
    FALSE 10           ABHD12        male          lo
    FALSE 11          ABHD17B        male        <NA>
    FALSE 12          ABHD17C        male          lo
    FALSE 13            ABHD2      female          hi
    FALSE 14            ABHD6      female        <NA>
    FALSE 15           ABLIM3      female        <NA>
    FALSE 16              ABR        male          lo
    FALSE 17            ACAA2      female        <NA>
    FALSE 18            ACADL      female          hi
    FALSE 19            ACAP3        male        <NA>
    FALSE 20            ACER2        male          lo
    FALSE 21            ACER3      female        <NA>
    FALSE 22            ACKR2      female        <NA>
    FALSE 23           ACOT11        male        <NA>
    FALSE 24           ACOT13      female        <NA>
    FALSE 25            ACOT9      female        <NA>
    FALSE 26            ACPL2      female        <NA>
    FALSE 27            ACSF2      female        <NA>
    FALSE 28            ACSL4        male        <NA>
    FALSE 29            ACSL5      female        <NA>
    FALSE 30             ACY1      female        <NA>
    FALSE 31           ADAM22      female          hi
    FALSE 32           ADAM23      female          hi
    FALSE 33         ADAMTS19        male        <NA>
    FALSE 34          ADAMTS3      female        <NA>
    FALSE 35          ADAMTS5      female        <NA>
    FALSE 36          ADAMTS6        male        <NA>
    FALSE 37         ADAMTSL2      female          hi
    FALSE 38         ADAMTSL3      female        <NA>
    FALSE 39           ADGRD1        male          lo
    FALSE 40           ADGRF5      female        <NA>
    FALSE 41           ADGRL4      female        <NA>
    FALSE 42           ADGRV1        male        <NA>
    FALSE 43           ADHFE1      female        <NA>
    FALSE 44           ADPRHL      female        <NA>
    FALSE 45            ADRB2        male        <NA>
    FALSE 46            AFAP1      female        <NA>
    FALSE 47          AFAP1L1      female        <NA>
    FALSE 48            AGBL4      female          hi
    FALSE 49            AGGF1        male        <NA>
    FALSE 50           AGPAT2      female        <NA>
    FALSE 51             AGR2      female          hi
    FALSE 52          AGTPBP1        male        <NA>
    FALSE 53           AGTRAP      female        <NA>
    FALSE 54            AGXT2        male        <NA>
    FALSE 55          AGXT2L2      female        <NA>
    FALSE 56             AHI1      female          hi
    FALSE 57             AIM1        male        <NA>
    FALSE 58            AIMP1      female          hi
    FALSE 59              AK3        male        <NA>
    FALSE 60              AK4      female        <NA>
    FALSE 61              AK6        male        <NA>
    FALSE 62              AK7      female        <NA>
    FALSE 63            AKTIP      female        <NA>
    FALSE 64            ALAS1      female        <NA>
    FALSE 65          ALDH1A1      female        <NA>
    FALSE 66          ALDH1A2      female        <NA>
    FALSE 67          ALDH5A1      female        <NA>
    FALSE 68          ALDH7A1        male          hi
    FALSE 69            AMACR        male        <NA>
    FALSE 70            AMPD3        male        <NA>
    FALSE 71             AMPH      female          lo
    FALSE 72            AMY1A      female          hi
    FALSE 73           ANGPT4      female        <NA>
    FALSE 74          ANGPTL4      female        <NA>
    FALSE 75          ANGPTL7      female        <NA>
    FALSE 76          ANKDD1B        male        <NA>
    FALSE 77           ANKRA2        male        <NA>
    FALSE 78          ANKRD10      female        <NA>
    FALSE 79          ANKRD26      female        <NA>
    FALSE 80          ANKRD29        male        <NA>
    FALSE 81          ANKRD32        male        <NA>
    FALSE 82           ANKS1B      female        <NA>
    FALSE 83             ANO1      female        <NA>
    FALSE 84            ANXA1        male          lo
    FALSE 85            AP3B1        male        <NA>
    FALSE 86            AP3S1        male        <NA>
    FALSE 87            AP4S1      female        <NA>
    FALSE 88            APBA1        male        <NA>
    FALSE 89          APBB1IP        male          lo
    FALSE 90             APEH      female        <NA>
    FALSE 91          APOBEC2      female        <NA>
    FALSE 92             APTX        male        <NA>
    FALSE 93               AR      female        <NA>
    FALSE 94            ARAP3      female          lo
    FALSE 95           ARFRP1        male          hi
    FALSE 96         ARHGAP24      female          lo
    FALSE 97          ARHGEF1      female        <NA>
    FALSE 98         ARHGEF16      female        <NA>
    FALSE 99         ARHGEF39        male        <NA>
    FALSE 100          ARID4B      female          hi
    FALSE 101           ARL15        male        <NA>
    FALSE 102           ARL16      female        <NA>
    FALSE 103           ARL4C        male        <NA>
    FALSE 104         ARL6IP5        male          hi
    FALSE 105           ARMC4        male        <NA>
    FALSE 106          ARRDC1        male        <NA>
    FALSE 107          ARRDC3      female        <NA>
    FALSE 108            ARSB        male        <NA>
    FALSE 109            ARSI        male          lo
    FALSE 110            ARSK        male        <NA>
    FALSE 111           ASB15      female        <NA>
    FALSE 112            ASB2        male          lo
    FALSE 113            ASB7        male        <NA>
    FALSE 114           ASIC4      female        <NA>
    FALSE 115            ASL2      female        <NA>
    FALSE 116            ATF1        male        <NA>
    FALSE 117           ATG12        male        <NA>
    FALSE 118           ATG4A      female          hi
    FALSE 119            ATIC      female        <NA>
    FALSE 120             ATM      female          hi
    FALSE 121           ATOH8      female        <NA>
    FALSE 122          ATP10A      female        <NA>
    FALSE 123          ATP12A        male        <NA>
    FALSE 124          ATP1B3        male        <NA>
    FALSE 125          ATP2A3        male        <NA>
    FALSE 126          ATP2C1      female        <NA>
    FALSE 127          ATP5A1        male        <NA>
    FALSE 128         ATP6V0B      female        <NA>
    FALSE 129        ATP6V1C2      female          hi
    FALSE 130          ATRNL1      female        <NA>
    FALSE 131             AUH        male        <NA>
    FALSE 132             AVD        male        <NA>
    FALSE 133          AVPR1B        male        <NA>
    FALSE 134        B4GALNT4      female        <NA>
    FALSE 135         B4GALT1        male        <NA>
    FALSE 136         B4GALT2      female          lo
    FALSE 137           BAMBI      female        <NA>
    FALSE 138           BBIP1      female        <NA>
    FALSE 139           BCAR3      female        <NA>
    FALSE 140           BCAT1      female          hi
    FALSE 141            BCHE      female        <NA>
    FALSE 142            BCO2        male        <NA>
    FALSE 143           BDP1L        male        <NA>
    FALSE 144           BEAN1        male        <NA>
    FALSE 145         BHLHA15      female        <NA>
    FALSE 146            BHMT        male        <NA>
    FALSE 147            BIVM      female          hi
    FALSE 148            BNC2        male        <NA>
    FALSE 149             BOC      female          lo
    FALSE 150             BOK      female          lo
    FALSE 151            BRAP        male          lo
    FALSE 152           BRIX1        male        <NA>
    FALSE 153             BSG      female        <NA>
    FALSE 154           BTBD6      female        <NA>
    FALSE 155           BTBD7      female        <NA>
    FALSE 156             BTC      female        <NA>
    FALSE 157            BTF3        male        <NA>
    FALSE 158      BX255923.2        male        <NA>
    FALSE 159        C10ORF10      female        <NA>
    FALSE 160     C11H19ORF40      female        <NA>
    FALSE 161        C11ORF52      female        <NA>
    FALSE 162      C12H3ORF67        male        <NA>
    FALSE 163    C14H17ORF103      female        <NA>
    FALSE 164     C15H12ORF49      female        <NA>
    FALSE 165     C17H9ORF114      female        <NA>
    FALSE 166        C17ORF58      female        <NA>
    FALSE 167      C1H12ORF63      female        <NA>
    FALSE 168       C1HXORF36      female        <NA>
    FALSE 169        C1ORF198      female        <NA>
    FALSE 170         C1QTNF5      female        <NA>
    FALSE 171             C1R      female        <NA>
    FALSE 172       C20ORF196      female        <NA>
    FALSE 173     C28H19ORF10      female          hi
    FALSE 174           C2CD2      female        <NA>
    FALSE 175       C2H8ORF46        male        <NA>
    FALSE 176      C5H11ORF58      female        <NA>
    FALSE 177      C5H11ORF96        male        <NA>
    FALSE 178     C5H14ORF159      female        <NA>
    FALSE 179         C5ORF30        male        <NA>
    FALSE 180         C5ORF49      female        <NA>
    FALSE 181     C6H10ORF118      female        <NA>
    FALSE 182      C6H10ORF76        male        <NA>
    FALSE 183      C6H10ORF90      female        <NA>
    FALSE 184              C7        male          lo
    FALSE 185             C8B        male        <NA>
    FALSE 186         C8ORF22      female        <NA>
    FALSE 187          C8ORF4      female          lo
    FALSE 188         C9ORF72        male        <NA>
    FALSE 189         C9ORF85        male        <NA>
    FALSE 190            CA10      female        <NA>
    FALSE 191             CA4      female        <NA>
    FALSE 192             CA9        male        <NA>
    FALSE 193           CAAP1        male        <NA>
    FALSE 194          CACFD1      female          hi
    FALSE 195         CACNA1B        male        <NA>
    FALSE 196         CACNA1D      female        <NA>
    FALSE 197           CADM2        male        <NA>
    FALSE 198           CADM3      female        <NA>
    FALSE 199          CALCRL      female        <NA>
    FALSE 200           CAMK4        male        <NA>
    FALSE 201          CAMTA1      female          hi
    FALSE 202           CAPN1        male          lo
    FALSE 203          CAPN14      female        <NA>
    FALSE 204           CAPN6      female        <NA>
    FALSE 205         CAPRIN2      female        <NA>
    FALSE 206           CAPSL        male        <NA>
    FALSE 207          CARD10      female        <NA>
    FALSE 208         CARHSP1        male        <NA>
    FALSE 209         CARNMT1        male        <NA>
    FALSE 210           CASC4        male        <NA>
    FALSE 211           CASD1      female        <NA>
    FALSE 212            CAST        male        <NA>
    FALSE 213            CBX6      female          hi
    FALSE 214           CCBE1        male        <NA>
    FALSE 215         CCDC125        male        <NA>
    FALSE 216         CCDC171        male        <NA>
    FALSE 217          CCDC60      female          hi
    FALSE 218          CCDC78      female        <NA>
    FALSE 219         CCDC85C      female        <NA>
    FALSE 220             CCK      female        <NA>
    FALSE 221            CCL1        male        <NA>
    FALSE 222           CCL28        male        <NA>
    FALSE 223            CCL5        male        <NA>
    FALSE 224           CCM2L      female        <NA>
    FALSE 225           CCNB2        male        <NA>
    FALSE 226           CCND2        male          hi
    FALSE 227           CCNE1      female        <NA>
    FALSE 228            CCNF      female        <NA>
    FALSE 229            CCNH        male        <NA>
    FALSE 230           CCSAP        male        <NA>
    FALSE 231          CCSER1        male        <NA>
    FALSE 232           CD151      female        <NA>
    FALSE 233           CD274        male          lo
    FALSE 234            CD34      female        <NA>
    FALSE 235            CD38      female        <NA>
    FALSE 236            CD47        male        <NA>
    FALSE 237             CD9      female        <NA>
    FALSE 238            CD99      female        <NA>
    FALSE 239             CDA      female        <NA>
    FALSE 240          CDC14B        male        <NA>
    FALSE 241         CDC37L1        male        <NA>
    FALSE 242        CDC42EP1      female        <NA>
    FALSE 243        CDC42SE2        male        <NA>
    FALSE 244          CDCA7L      female        <NA>
    FALSE 245           CDCP1        male          hi
    FALSE 246           CDH20      female        <NA>
    FALSE 247            CDH4      female          lo
    FALSE 248            CDH5      female        <NA>
    FALSE 249           CDHR5        male        <NA>
    FALSE 250            CDK7        male        <NA>
    FALSE 251          CDKN1A      female          hi
    FALSE 252          CDKN1B      female        <NA>
    FALSE 253          CDKN2C        male        <NA>
    FALSE 254            CDO1        male          lo
    FALSE 255           CEBPA      female        <NA>
    FALSE 256          CECR5L        male        <NA>
    FALSE 257           CELF2      female        <NA>
    FALSE 258           CELF3        male        <NA>
    FALSE 259           CELF4        male        <NA>
    FALSE 260           CELF5        male          lo
    FALSE 261           CENPV        male        <NA>
    FALSE 262          CEP120        male        <NA>
    FALSE 263          CEP162      female        <NA>
    FALSE 264         CEP170B      female        <NA>
    FALSE 265           CEP41      female        <NA>
    FALSE 266           CEP78        male          hi
    FALSE 267           CEP83      female          hi
    FALSE 268           CEP85        male        <NA>
    FALSE 269            CERK        male          hi
    FALSE 270           CERS4      female        <NA>
    FALSE 271           CETN3        male          lo
    FALSE 272            CETP        male        <NA>
    FALSE 273          CFAP44      female        <NA>
    FALSE 274           CFHR2        male          lo
    FALSE 275            CHAT      female        <NA>
    FALSE 276         CHCHD10      female        <NA>
    FALSE 277            CHDH      female        <NA>
    FALSE 278            CHGA      female        <NA>
    FALSE 279            CHKA      female          hi
    FALSE 280          CHRDL1        male          lo
    FALSE 281          CHRNA3        male        <NA>
    FALSE 282          CHRNA6        male        <NA>
    FALSE 283          CHRNB4        male          hi
    FALSE 284          CHST11        male        <NA>
    FALSE 285          CHST15      female        <NA>
    FALSE 286           CHST3      female        <NA>
    FALSE 287           CHST9      female        <NA>
    FALSE 288         CIAPIN1      female        <NA>
    FALSE 289           CIDEC      female        <NA>
    FALSE 290           CIRBP      female        <NA>
    FALSE 291           CISD2      female        <NA>
    FALSE 292            CISH        male        <NA>
    FALSE 293          CITED4      female          hi
    FALSE 294             CKB      female        <NA>
    FALSE 295          CKMT1A      female        <NA>
    FALSE 296            CKS2        male          hi
    FALSE 297          CLDN11      female          hi
    FALSE 298           CLIC4        male        <NA>
    FALSE 299           CLIC6        male          lo
    FALSE 300            CLN5      female        <NA>
    FALSE 301            CLTA        male          hi
    FALSE 302             CLU        male          lo
    FALSE 303            CMAS      female          hi
    FALSE 304           CMYA5        male        <NA>
    FALSE 305          CNKSR3      female          lo
    FALSE 306          CNOT10      female        <NA>
    FALSE 307            CNR1        male        <NA>
    FALSE 308            CNR2      female        <NA>
    FALSE 309            CNST      female          hi
    FALSE 310           CNTFR      female          hi
    FALSE 311           CNTLN        male          hi
    FALSE 312         CNTNAP4        male        <NA>
    FALSE 313            COG1      female        <NA>
    FALSE 314         COL11A1      female        <NA>
    FALSE 315         COL12A1        male        <NA>
    FALSE 316         COL14A1        male        <NA>
    FALSE 317         COL18A1        male        <NA>
    FALSE 318         COL26A1      female        <NA>
    FALSE 319        COL4A3BP        male        <NA>
    FALSE 320          COL4A5      female        <NA>
    FALSE 321          COL6A1        male        <NA>
    FALSE 322          COL6A2        male          lo
    FALSE 323          COL6A3      female        <NA>
    FALSE 324         COMMD10        male        <NA>
    FALSE 325          COMMD9      female        <NA>
    FALSE 326            COMT      female        <NA>
    FALSE 327           COPB2      female          hi
    FALSE 328            COPE      female          hi
    FALSE 329           COPG2        male        <NA>
    FALSE 330          COPS7B      female          hi
    FALSE 331            COQ6      female        <NA>
    FALSE 332            COQ9      female          hi
    FALSE 333          CORO2A        male        <NA>
    FALSE 334           CORO7      female        <NA>
    FALSE 335           COTL1      female        <NA>
    FALSE 336           CPEB2      female          hi
    FALSE 337           CPLX1        male        <NA>
    FALSE 338             CPM        male        <NA>
    FALSE 339           CPNE8      female          lo
    FALSE 340          CPPED1      female        <NA>
    FALSE 341             CPQ      female        <NA>
    FALSE 342            CPT2      female        <NA>
    FALSE 343            CR1L      female        <NA>
    FALSE 344           CREB3        male        <NA>
    FALSE 345            CREM        male          hi
    FALSE 346           CRHR1        male        <NA>
    FALSE 347        CRISPLD2        male        <NA>
    FALSE 348             CRK        male        <NA>
    FALSE 349            CRKL      female        <NA>
    FALSE 350            CROT      female        <NA>
    FALSE 351          CRYBA2      female        <NA>
    FALSE 352          CRYBG3      female        <NA>
    FALSE 353      CSGALNACT1      female        <NA>
    FALSE 354           CSMD3      female        <NA>
    FALSE 355         CSNK1G2        male        <NA>
    FALSE 356           CSRP1      female        <NA>
    FALSE 357            CSTB      female        <NA>
    FALSE 358           CTBP2      female        <NA>
    FALSE 359            CTBS      female          hi
    FALSE 360    CTC-487M23.8        male          hi
    FALSE 361     CTC-554D6.1        male        <NA>
    FALSE 362   CTD-2287O16.3        male        <NA>
    FALSE 363          CTDSP1      female        <NA>
    FALSE 364          CTHRC1        male          lo
    FALSE 365            CTIF        male        <NA>
    FALSE 366         CTNNBL1      female        <NA>
    FALSE 367            CTNS      female        <NA>
    FALSE 368           CTPS2      female        <NA>
    FALSE 369            CTSA      female          lo
    FALSE 370           CTSL2      female          hi
    FALSE 371            CTSO      female        <NA>
    FALSE 372         CTTNBP2      female        <NA>
    FALSE 373            CUTA      female        <NA>
    FALSE 374           CWC27        male        <NA>
    FALSE 375          CYB5R2      female        <NA>
    FALSE 376          CYBRD1      female        <NA>
    FALSE 377            CYGB        male          lo
    FALSE 378         CYP26A1        male        <NA>
    FALSE 379          CYP3A7      female        <NA>
    FALSE 380         CYP4B1L      female        <NA>
    FALSE 381           CYR61        male        <NA>
    FALSE 382           CYTH4      female        <NA>
    FALSE 383           CYYR1      female        <NA>
    FALSE 384      CZH18ORF25        male        <NA>
    FALSE 385       CZH5ORF28        male        <NA>
    FALSE 386       CZH5ORF34        male        <NA>
    FALSE 387       CZH5ORF42        male        <NA>
    FALSE 388       CZH5ORF51        male          hi
    FALSE 389       CZH5ORF63        male        <NA>
    FALSE 390        CZH9ORF3        male        <NA>
    FALSE 391       CZH9ORF40        male        <NA>
    FALSE 392       CZH9ORF64        male        <NA>
    FALSE 393       CZH9ORF84        male        <NA>
    FALSE 394           DAAM1        male          lo
    FALSE 395            DAB2        male          lo
    FALSE 396           DACT1      female        <NA>
    FALSE 397           DACT2        male          lo
    FALSE 398           DAPK1        male        <NA>
    FALSE 399            DAW1        male        <NA>
    FALSE 400            DBN1        male          lo
    FALSE 401          DBNDD2      female        <NA>
    FALSE 402          DCAF10        male          lo
    FALSE 403          DCAF12        male        <NA>
    FALSE 404           DCAF4      female        <NA>
    FALSE 405           DCAF6      female        <NA>
    FALSE 406           DCDC2      female        <NA>
    FALSE 407           DCHS1      female        <NA>
    FALSE 408             DCN        male        <NA>
    FALSE 409            DCP2        male        <NA>
    FALSE 410           DCTN3        male        <NA>
    FALSE 411             DDC      female        <NA>
    FALSE 412           DDOST      female          hi
    FALSE 413           DDX3X        male        <NA>
    FALSE 414            DDX4        male        <NA>
    FALSE 415           DEAF1        male        <NA>
    FALSE 416           DECR1      female        <NA>
    FALSE 417         DENND2D      female          lo
    FALSE 418         DENND4C        male        <NA>
    FALSE 419          DEPTOR      female        <NA>
    FALSE 420           DERL2      female          hi
    FALSE 421          DFNB59      female        <NA>
    FALSE 422            DGKE        male        <NA>
    FALSE 423            DGKQ        male        <NA>
    FALSE 424            DGKZ      female        <NA>
    FALSE 425          DHCR24        male        <NA>
    FALSE 426            DHFR        male        <NA>
    FALSE 427           DHRS3        male        <NA>
    FALSE 428           DHRS4      female        <NA>
    FALSE 429           DHX29        male        <NA>
    FALSE 430           DHX32        male        <NA>
    FALSE 431           DIMT1        male        <NA>
    FALSE 432            DIO2        male          lo
    FALSE 433          DIRAS1        male          hi
    FALSE 434           DIS3L      female        <NA>
    FALSE 435            DLAT      female          lo
    FALSE 436             DLD      female        <NA>
    FALSE 437          DLGAP3        male          lo
    FALSE 438           DMGDH        male        <NA>
    FALSE 439           DNAH9        male          lo
    FALSE 440           DNAI1        male        <NA>
    FALSE 441          DNAJA1        male        <NA>
    FALSE 442          DNAJB5        male        <NA>
    FALSE 443          DNAJC1      female          hi
    FALSE 444         DNAJC21        male        <NA>
    FALSE 445         DNAJC25        male        <NA>
    FALSE 446          DNAJC6        male        <NA>
    FALSE 447           DNAL1      female        <NA>
    FALSE 448           DNAL4      female        <NA>
    FALSE 449        DNASE1L3      female        <NA>
    FALSE 450           DNM1L      female        <NA>
    FALSE 451           DOCK8        male        <NA>
    FALSE 452            DOK3        male        <NA>
    FALSE 453          DPAGT1      female        <NA>
    FALSE 454            DPH6      female          hi
    FALSE 455            DPP6        male          lo
    FALSE 456            DPP7      female          hi
    FALSE 457           DRCC1      female        <NA>
    FALSE 458            DRD1        male        <NA>
    FALSE 459            DRD4        male        <NA>
    FALSE 460             DSE        male        <NA>
    FALSE 461           DTHD1        male        <NA>
    FALSE 462           DTWD2        male        <NA>
    FALSE 463            DTX1      female        <NA>
    FALSE 464          DUSP10        male          hi
    FALSE 465           DUSP8        male        <NA>
    FALSE 466             DYM      female        <NA>
    FALSE 467         DYNC2H1      female        <NA>
    FALSE 468          DYNLT3      female        <NA>
    FALSE 469           ECEL1      female        <NA>
    FALSE 470          ECHDC2      female        <NA>
    FALSE 471          ECHDC3      female        <NA>
    FALSE 472           ECHS1      female        <NA>
    FALSE 473            ECI1      female        <NA>
    FALSE 474            EDN1        male        <NA>
    FALSE 475          EEF1A1      female          hi
    FALSE 476            EEF2      female          hi
    FALSE 477          EEFSEC      female        <NA>
    FALSE 478           EEPD1        male          lo
    FALSE 479          EFCAB1      female        <NA>
    FALSE 480           EFHD2      female        <NA>
    FALSE 481           EFNA5        male        <NA>
    FALSE 482           EGFL6      female        <NA>
    FALSE 483           EGFL7      female        <NA>
    FALSE 484            EHD3      female        <NA>
    FALSE 485           EIF3I      female        <NA>
    FALSE 486           EIF3L      female          hi
    FALSE 487          EIF4G2        male          hi
    FALSE 488            ELF1      female        <NA>
    FALSE 489            ELK3      female        <NA>
    FALSE 490          ELOVL1      female        <NA>
    FALSE 491             EMB        male          hi
    FALSE 492            EMCN      female        <NA>
    FALSE 493           EMID1      female        <NA>
    FALSE 494         EMILIN3        male        <NA>
    FALSE 495            ENAH      female        <NA>
    FALSE 496            ENO1      female          hi
    FALSE 497            ENO2      female          hi
    FALSE 498          ENOPH1      female        <NA>
    FALSE 499           ENPP2        male        <NA>
    FALSE 500           ENS-3      female        <NA>
    FALSE 501          ENTPD4        male          hi
    FALSE 502           EPAS1      female        <NA>
    FALSE 503         EPB41L1        male        <NA>
    FALSE 504        EPB41L4A        male        <NA>
    FALSE 505            EPG5        male        <NA>
    FALSE 506           EPHB1      female        <NA>
    FALSE 507           EPHX4      female          hi
    FALSE 508           ERAP1        male          lo
    FALSE 509         ERBB2IP        male        <NA>
    FALSE 510           ERCC5      female        <NA>
    FALSE 511         ERCC6L2        male        <NA>
    FALSE 512           ERCC8        male        <NA>
    FALSE 513          ERICH5        male          lo
    FALSE 514            ERMN        male        <NA>
    FALSE 515           ERMP1        male        <NA>
    FALSE 516          ERRFI1      female        <NA>
    FALSE 517            ESAM      female          lo
    FALSE 518           ESRP1      female          hi
    FALSE 519            ETFA      female          hi
    FALSE 520           ETFDH      female        <NA>
    FALSE 521            ETS1      female          lo
    FALSE 522             EVL      female        <NA>
    FALSE 523            EXD3        male          lo
    FALSE 524           EXOC5        male        <NA>
    FALSE 525           EXPH5      female        <NA>
    FALSE 526            EXT2      female          hi
    FALSE 527             F10        male          lo
    FALSE 528           F13A1        male          hi
    FALSE 529             F2R        male        <NA>
    FALSE 530           F2RL2        male        <NA>
    FALSE 531              F8      female        <NA>
    FALSE 532            FAAH      female        <NA>
    FALSE 533           FADS2      female        <NA>
    FALSE 534         FAM101A        male        <NA>
    FALSE 535         FAM102B        male        <NA>
    FALSE 536         FAM110C        male          lo
    FALSE 537        FAM114A1      female          hi
    FALSE 538         FAM120B      female        <NA>
    FALSE 539         FAM122A        male        <NA>
    FALSE 540         FAM131A        male        <NA>
    FALSE 541         FAM132A      female        <NA>
    FALSE 542         FAM135B      female          hi
    FALSE 543          FAM13B        male        <NA>
    FALSE 544         FAM161B        male        <NA>
    FALSE 545         FAM168A        male        <NA>
    FALSE 546         FAM169A        male        <NA>
    FALSE 547         FAM172A        male        <NA>
    FALSE 548         FAM174A        male        <NA>
    FALSE 549         FAM174B        male        <NA>
    FALSE 550         FAM181A        male          lo
    FALSE 551         FAM195A      female        <NA>
    FALSE 552         FAM19A5      female        <NA>
    FALSE 553          FAM20A        male        <NA>
    FALSE 554          FAM26F      female          lo
    FALSE 555          FAM35A      female        <NA>
    FALSE 556          FAM43A      female          lo
    FALSE 557          FAM73A      female          hi
    FALSE 558          FAM83B      female          hi
    FALSE 559            FAN1      female          hi
    FALSE 560           FANCC        male        <NA>
    FALSE 561           FARP1      female        <NA>
    FALSE 562            FAXC        male        <NA>
    FALSE 563           FBLN1        male          lo
    FALSE 564            FBN2        male        <NA>
    FALSE 565            FBP1        male        <NA>
    FALSE 566          FBXL17        male        <NA>
    FALSE 567           FBXL3        male        <NA>
    FALSE 568           FBXL4      female        <NA>
    FALSE 569          FBXO32        male        <NA>
    FALSE 570          FBXO33      female        <NA>
    FALSE 571           FBXO4        male        <NA>
    FALSE 572           FBXW4      female        <NA>
    FALSE 573           FBXW8      female        <NA>
    FALSE 574           FCGBP      female        <NA>
    FALSE 575           FCHO2        male        <NA>
    FALSE 576            FCN2      female        <NA>
    FALSE 577            FDPS      female        <NA>
    FALSE 578           FDX1L      female        <NA>
    FALSE 579           FEM1C        male        <NA>
    FALSE 580            FEN1        male        <NA>
    FALSE 581             FER        male          hi
    FALSE 582           FGF10        male        <NA>
    FALSE 583            FGF2        male        <NA>
    FALSE 584            FGF6        male          hi
    FALSE 585             FGG        male        <NA>
    FALSE 586            FGGY      female        <NA>
    FALSE 587              FH      female        <NA>
    FALSE 588            FHIT      female        <NA>
    FALSE 589            FHL2      female        <NA>
    FALSE 590         FILIP1L        male          lo
    FALSE 591            FKTN        male        <NA>
    FALSE 592            FLT1      female        <NA>
    FALSE 593            FLT4      female          lo
    FALSE 594          FLVCR2      female        <NA>
    FALSE 595           FNDC1        male          lo
    FALSE 596          FNDC3B      female        <NA>
    FALSE 597           FOCAD        male          hi
    FALSE 598           FOSL2        male          hi
    FALSE 599           FOXA1      female        <NA>
    FALSE 600           FOXN4      female        <NA>
    FALSE 601           FOXP1      female        <NA>
    FALSE 602           FREM1        male        <NA>
    FALSE 603           FREM3        male        <NA>
    FALSE 604           FRMD1      female        <NA>
    FALSE 605           FRMD3        male        <NA>
    FALSE 606          FRMD4B      female          lo
    FALSE 607          FRMPD1        male        <NA>
    FALSE 608          FRMPD4        male          lo
    FALSE 609            FRZB      female        <NA>
    FALSE 610           FSD1L        male        <NA>
    FALSE 611            FSHB      female        <NA>
    FALSE 612             FST        male        <NA>
    FALSE 613             FTL        male          hi
    FALSE 614          FUNDC2      female        <NA>
    FALSE 615            FUT8      female        <NA>
    FALSE 616             FXN        male        <NA>
    FALSE 617             FYB        male          lo
    FALSE 618             FYN      female        <NA>
    FALSE 619            FZD4      female          lo
    FALSE 620            FZD8      female        <NA>
    FALSE 621            FZR1        male        <NA>
    FALSE 622            GAB2        male          hi
    FALSE 623            GAB3      female        <NA>
    FALSE 624          GABRG3      female        <NA>
    FALSE 625           GABRQ      female        <NA>
    FALSE 626            GAD2        male        <NA>
    FALSE 627         GADD45B        male        <NA>
    FALSE 628           GADL1      female          hi
    FALSE 629             GAK        male        <NA>
    FALSE 630            GALC      female        <NA>
    FALSE 631          GALNT1        male        <NA>
    FALSE 632          GALNT5        male          lo
    FALSE 633            GALT        male        <NA>
    FALSE 634           GAP43      female        <NA>
    FALSE 635           GAPDH      female        <NA>
    FALSE 636            GAS1        male          lo
    FALSE 637            GBA2        male        <NA>
    FALSE 638            GBAS      female        <NA>
    FALSE 639           GBGT1      female        <NA>
    FALSE 640              GC      female        <NA>
    FALSE 641            GCC2      female          hi
    FALSE 642            GCH1        male          hi
    FALSE 643           GCNT1        male        <NA>
    FALSE 644           GCNT4        male        <NA>
    FALSE 645           GCNT7      female        <NA>
    FALSE 646         GDAP1L1        male        <NA>
    FALSE 647           GDAP2      female        <NA>
    FALSE 648           GDF11        male        <NA>
    FALSE 649            GDF7        male          lo
    FALSE 650           GDPD1        male          hi
    FALSE 651            GET4      female        <NA>
    FALSE 652            GFM1      female        <NA>
    FALSE 653            GFM2        male        <NA>
    FALSE 654           GHRHR      female        <NA>
    FALSE 655            GIN1        male        <NA>
    FALSE 656           GIPC2      female        <NA>
    FALSE 657            GJA1      female        <NA>
    FALSE 658            GJA3      female        <NA>
    FALSE 659            GJA4      female          lo
    FALSE 660            GJB5        male        <NA>
    FALSE 661              GK        male          hi
    FALSE 662          GLCCI1      female        <NA>
    FALSE 663            GLDC        male          lo
    FALSE 664           GLP1R        male        <NA>
    FALSE 665         GLTSCR2      female        <NA>
    FALSE 666           GLUD1      female        <NA>
    FALSE 667            GMFB      female        <NA>
    FALSE 668           GMPPB      female          hi
    FALSE 669           GNA11      female          hi
    FALSE 670            GNAQ        male        <NA>
    FALSE 671            GNAZ        male        <NA>
    FALSE 672             GNE        male        <NA>
    FALSE 673           GNG10        male        <NA>
    FALSE 674            GNG2        male          lo
    FALSE 675          GNPTAB      female        <NA>
    FALSE 676          GOLGA1        male        <NA>
    FALSE 677           GOLM1        male        <NA>
    FALSE 678          GOLPH3        male        <NA>
    FALSE 679         GORASP2      female          hi
    FALSE 680            GOT2      female          hi
    FALSE 681           GPAT2      female        <NA>
    FALSE 682         GPBP1L1      female        <NA>
    FALSE 683            GPC5        male        <NA>
    FALSE 684            GPC6      female          hi
    FALSE 685          GPD1L2        male        <NA>
    FALSE 686            GPD2        male          lo
    FALSE 687           GPER1        male        <NA>
    FALSE 688         GPR137B      female        <NA>
    FALSE 689           GPR56      female        <NA>
    FALSE 690           GPR64      female        <NA>
    FALSE 691           GPR85      female        <NA>
    FALSE 692           GPSM1      female        <NA>
    FALSE 693            GPT2        male        <NA>
    FALSE 694            GPX8        male        <NA>
    FALSE 695         GRAMD1B        male          lo
    FALSE 696           GRB10      female          hi
    FALSE 697           GREB1      female        <NA>
    FALSE 698           GRHPR        male        <NA>
    FALSE 699           GRIA3      female        <NA>
    FALSE 700           GRIN1      female        <NA>
    FALSE 701          GRIN3A        male          hi
    FALSE 702           GRIP1        male          hi
    FALSE 703           GRTP1      female        <NA>
    FALSE 704           GSTO1      female        <NA>
    FALSE 705           GSTZ1      female        <NA>
    FALSE 706          GTF2H2        male        <NA>
    FALSE 707           HABP4        male        <NA>
    FALSE 708           HAUS6        male        <NA>
    FALSE 709            HCN1        male        <NA>
    FALSE 710            HCN4        male        <NA>
    FALSE 711           HDAC4      female        <NA>
    FALSE 712           HDHD2        male        <NA>
    FALSE 713           HEBP1      female        <NA>
    FALSE 714            HEG1      female        <NA>
    FALSE 715         HERPUD2      female        <NA>
    FALSE 716          HGSNAT      female        <NA>
    FALSE 717           HIBCH      female        <NA>
    FALSE 718           HINT1        male        <NA>
    FALSE 719           HINT2        male        <NA>
    FALSE 720           HIPK1        male          lo
    FALSE 721      HIST1H2B5L      female        <NA>
    FALSE 722       HIST1H2B7      female        <NA>
    FALSE 723           HMGCR        male          hi
    FALSE 724            HNMT      female        <NA>
    FALSE 725          HNRNPK        male        <NA>
    FALSE 726           HOOK3        male          hi
    FALSE 727          HP1BP3      female        <NA>
    FALSE 728          HS3ST1      female        <NA>
    FALSE 729          HS3ST5      female          hi
    FALSE 730          HS6ST3        male          hi
    FALSE 731         HSD17B4        male          hi
    FALSE 732           HSDL2        male        <NA>
    FALSE 733           HSPB1        male          lo
    FALSE 734           HTR1F      female        <NA>
    FALSE 735           HTR2C        male          lo
    FALSE 736            HTR7      female        <NA>
    FALSE 737           HYAL3      female        <NA>
    FALSE 738             HYI      female        <NA>
    FALSE 739            ICA1      female          hi
    FALSE 740             ID1      female        <NA>
    FALSE 741             ID2      female        <NA>
    FALSE 742             ID3      female        <NA>
    FALSE 743            IDH2      female        <NA>
    FALSE 744            IDNK        male        <NA>
    FALSE 745            IDUA        male        <NA>
    FALSE 746         IER3IP1        male        <NA>
    FALSE 747           IFI30      female        <NA>
    FALSE 748         IFITM10        male          lo
    FALSE 749          IFNAR2      female        <NA>
    FALSE 750          IFT172        male        <NA>
    FALSE 751           IFT46      female        <NA>
    FALSE 752           IFT74        male        <NA>
    FALSE 753           IFT81        male        <NA>
    FALSE 754          IGFBP4        male          lo
    FALSE 755          IGFBP7      female        <NA>
    FALSE 756          IGSF11      female          lo
    FALSE 757           IGSF3        male          lo
    FALSE 758          IKBKAP        male        <NA>
    FALSE 759           IL12A        male          lo
    FALSE 760            IL18        male        <NA>
    FALSE 761        IL1RAPL2        male        <NA>
    FALSE 762          IL20RA        male          lo
    FALSE 763          IL31RA        male        <NA>
    FALSE 764           IL6ST        male        <NA>
    FALSE 765            IL7R        male        <NA>
    FALSE 766            IMMT      female        <NA>
    FALSE 767          INDOL1      female        <NA>
    FALSE 768            ING1      female        <NA>
    FALSE 769           INHBA      female        <NA>
    FALSE 770            INIP        male        <NA>
    FALSE 771           INSM1        male          hi
    FALSE 772           IQCB1      female        <NA>
    FALSE 773          IQGAP2        male        <NA>
    FALSE 774          IQSEC3        male        <NA>
    FALSE 775            IQUB      female        <NA>
    FALSE 776           ISCA1        male          hi
    FALSE 777            ISLR        male          lo
    FALSE 778            ISM1        male        <NA>
    FALSE 779           ISOC1        male        <NA>
    FALSE 780          ISYNA1      female        <NA>
    FALSE 781           ITGA2        male          lo
    FALSE 782           ITGA8      female        <NA>
    FALSE 783           ITGB5        male          lo
    FALSE 784           ITPR1      female        <NA>
    FALSE 785           ITPR2      female          lo
    FALSE 786           ITPR3      female        <NA>
    FALSE 787        ITPRIPL2        male          lo
    FALSE 788           ITSN1      female        <NA>
    FALSE 789        IVNS1ABP        male          hi
    FALSE 790           JADE3        male        <NA>
    FALSE 791         JAKMIP1        male        <NA>
    FALSE 792         JAKMIP3        male        <NA>
    FALSE 793            JAM3        male        <NA>
    FALSE 794            JDP2      female        <NA>
    FALSE 795           JMJD8      female        <NA>
    FALSE 796             JMY        male        <NA>
    FALSE 797           KAT2B        male        <NA>
    FALSE 798         KATNAL2        male        <NA>
    FALSE 799          KCNAB1        male        <NA>
    FALSE 800           KCND2        male          lo
    FALSE 801           KCNE3        male        <NA>
    FALSE 802           KCNG4        male          lo
    FALSE 803           KCNH7      female        <NA>
    FALSE 804           KCNJ2      female        <NA>
    FALSE 805           KCNJ3        male        <NA>
    FALSE 806           KCNJ8      female        <NA>
    FALSE 807          KCNK15      female        <NA>
    FALSE 808           KCNK5      female        <NA>
    FALSE 809           KCNN2        male          lo
    FALSE 810           KCNQ5      female        <NA>
    FALSE 811          KCTD12      female        <NA>
    FALSE 812          KCTD14        male          hi
    FALSE 813          KCTD15        male          lo
    FALSE 814           KCTD2        male        <NA>
    FALSE 815           KCTD7        male        <NA>
    FALSE 816          KDELR3      female          hi
    FALSE 817           KDM4C        male        <NA>
    FALSE 818           KDM6A      female        <NA>
    FALSE 819         KHDRBS2        male        <NA>
    FALSE 820        KIAA0020        male        <NA>
    FALSE 821        KIAA0825        male        <NA>
    FALSE 822        KIAA1161        male        <NA>
    FALSE 823        KIAA1191        male          lo
    FALSE 824        KIAA1217        male          lo
    FALSE 825        KIAA1328        male        <NA>
    FALSE 826        KIAA1614      female        <NA>
    FALSE 827        KIAA2026        male        <NA>
    FALSE 828           KIF24        male        <NA>
    FALSE 829           KIF27        male        <NA>
    FALSE 830           KIF2A        male        <NA>
    FALSE 831             KIN      female        <NA>
    FALSE 832              KL      female        <NA>
    FALSE 833           KLF10      female        <NA>
    FALSE 834           KLF11      female        <NA>
    FALSE 835            KLF4        male          lo
    FALSE 836            KLF9        male        <NA>
    FALSE 837           KLHL4      female        <NA>
    FALSE 838           KLHL7      female        <NA>
    FALSE 839           KRT12      female        <NA>
    FALSE 840          LACTB2      female        <NA>
    FALSE 841           LAMP5      female        <NA>
    FALSE 842          LANCL1      female        <NA>
    FALSE 843          LANCL2        male        <NA>
    FALSE 844            LAP3      female          hi
    FALSE 845           LARS2      female        <NA>
    FALSE 846           LAS1L      female        <NA>
    FALSE 847            LCOR      female        <NA>
    FALSE 848          LGALS1      female        <NA>
    FALSE 849            LGI1      female          lo
    FALSE 850            LGI2      female        <NA>
    FALSE 851            LIPG      female        <NA>
    FALSE 852           LMAN1      female          hi
    FALSE 853           LMAN2      female        <NA>
    FALSE 854          LMBR1L      female        <NA>
    FALSE 855          LMBRD2        male        <NA>
    FALSE 856            LMNA      female        <NA>
    FALSE 857           LMNB1        male          hi
    FALSE 858            LMO2      female        <NA>
    FALSE 859           LNPEP        male        <NA>
    FALSE 860    LOC100857116        male        <NA>
    FALSE 861    LOC100857197        male        <NA>
    FALSE 862    LOC100857280      female        <NA>
    FALSE 863    LOC100857745        male        <NA>
    FALSE 864    LOC100857781      female        <NA>
    FALSE 865    LOC100858115        male        <NA>
    FALSE 866    LOC100858127        male        <NA>
    FALSE 867    LOC100858130        male          hi
    FALSE 868    LOC100858139        male        <NA>
    FALSE 869    LOC100858326      female        <NA>
    FALSE 870    LOC100858579        male        <NA>
    FALSE 871    LOC100858626      female        <NA>
    FALSE 872    LOC100858693      female        <NA>
    FALSE 873    LOC100858799        male        <NA>
    FALSE 874    LOC100858919        male          lo
    FALSE 875    LOC100859020      female        <NA>
    FALSE 876    LOC100859272      female          lo
    FALSE 877    LOC100859284        male        <NA>
    FALSE 878    LOC100859434      female        <NA>
    FALSE 879    LOC100859478      female        <NA>
    FALSE 880    LOC100859840      female        <NA>
    FALSE 881    LOC100859848      female        <NA>
    FALSE 882    LOC100859904        male        <NA>
    FALSE 883    LOC101747901        male        <NA>
    FALSE 884    LOC101748065        male        <NA>
    FALSE 885    LOC101748084        male          lo
    FALSE 886    LOC101748153      female        <NA>
    FALSE 887    LOC101748283        male        <NA>
    FALSE 888    LOC101748344      female        <NA>
    FALSE 889    LOC101748527      female        <NA>
    FALSE 890    LOC101748543      female        <NA>
    FALSE 891    LOC101748614        male        <NA>
    FALSE 892    LOC101749230        male        <NA>
    FALSE 893    LOC101749375      female        <NA>
    FALSE 894    LOC101749492        male        <NA>
    FALSE 895    LOC101750188      female        <NA>
    FALSE 896    LOC101750393        male        <NA>
    FALSE 897    LOC101750586      female        <NA>
    FALSE 898    LOC101750711        male        <NA>
    FALSE 899    LOC101751203      female        <NA>
    FALSE 900    LOC101751827      female        <NA>
    FALSE 901    LOC101752289      female        <NA>
    FALSE 902    LOC107049015        male        <NA>
    FALSE 903    LOC107049275      female        <NA>
    FALSE 904    LOC107049309        male        <NA>
    FALSE 905    LOC107049315        male          lo
    FALSE 906    LOC107049317        male        <NA>
    FALSE 907    LOC107049322        male        <NA>
    FALSE 908    LOC107049324        male        <NA>
    FALSE 909    LOC107049327      female        <NA>
    FALSE 910    LOC107049346        male        <NA>
    FALSE 911    LOC107050049      female        <NA>
    FALSE 912    LOC107050292      female        <NA>
    FALSE 913    LOC107050350        male        <NA>
    FALSE 914    LOC107050373        male          lo
    FALSE 915    LOC107050569        male        <NA>
    FALSE 916    LOC107050608        male        <NA>
    FALSE 917    LOC107050724      female          hi
    FALSE 918    LOC107051161        male          lo
    FALSE 919    LOC107052324        male        <NA>
    FALSE 920    LOC107052389        male        <NA>
    FALSE 921    LOC107052863      female        <NA>
    FALSE 922    LOC107053040        male          hi
    FALSE 923    LOC107053893      female        <NA>
    FALSE 924    LOC107054425      female        <NA>
    FALSE 925    LOC107054798        male          hi
    FALSE 926    LOC107055024        male          lo
    FALSE 927    LOC107055115      female        <NA>
    FALSE 928    LOC107055116      female        <NA>
    FALSE 929    LOC107055438      female        <NA>
    FALSE 930    LOC107055637      female        <NA>
    FALSE 931    LOC107055647      female        <NA>
    FALSE 932    LOC107055731      female        <NA>
    FALSE 933    LOC107056154      female        <NA>
    FALSE 934    LOC107056260      female        <NA>
    FALSE 935    LOC107056348      female        <NA>
    FALSE 936    LOC107056412      female        <NA>
    FALSE 937    LOC107056972        male        <NA>
    FALSE 938    LOC107057057      female        <NA>
    FALSE 939    LOC107057238        male          lo
    FALSE 940    LOC107057567      female        <NA>
    FALSE 941       LOC378902      female          hi
    FALSE 942       LOC395159      female          hi
    FALSE 943       LOC415324      female        <NA>
    FALSE 944       LOC415664      female        <NA>
    FALSE 945       LOC415852      female          lo
    FALSE 946       LOC416696      female        <NA>
    FALSE 947       LOC417253      female          hi
    FALSE 948       LOC417414      female        <NA>
    FALSE 949       LOC418189      female        <NA>
    FALSE 950       LOC418298      female        <NA>
    FALSE 951       LOC419677      female        <NA>
    FALSE 952       LOC420716      female        <NA>
    FALSE 953       LOC421106      female          hi
    FALSE 954       LOC421419      female        <NA>
    FALSE 955       LOC421892      female        <NA>
    FALSE 956       LOC421975      female        <NA>
    FALSE 957       LOC422171      female        <NA>
    FALSE 958       LOC422316      female        <NA>
    FALSE 959       LOC422528        male        <NA>
    FALSE 960       LOC422926      female          hi
    FALSE 961       LOC423247      female        <NA>
    FALSE 962       LOC423752      female          lo
    FALSE 963       LOC424167      female        <NA>
    FALSE 964       LOC424300      female        <NA>
    FALSE 965       LOC427201        male        <NA>
    FALSE 966       LOC427400        male        <NA>
    FALSE 967       LOC427439      female        <NA>
    FALSE 968       LOC427826      female        <NA>
    FALSE 969       LOC428479        male          hi
    FALSE 970       LOC430303      female        <NA>
    FALSE 971       LOC768709        male          lo
    FALSE 972       LOC769000      female        <NA>
    FALSE 973       LOC769052      female        <NA>
    FALSE 974       LOC769421      female        <NA>
    FALSE 975       LOC770345        male          lo
    FALSE 976       LOC770429      female        <NA>
    FALSE 977       LOC770492      female        <NA>
    FALSE 978       LOC771735      female        <NA>
    FALSE 979       LOC772243      female        <NA>
    FALSE 980          LPCAT1      female        <NA>
    FALSE 981             LPL        male        <NA>
    FALSE 982           LPPR1        male          lo
    FALSE 983            LRAT        male          lo
    FALSE 984           LRFN3        male          lo
    FALSE 985           LRP12        male          hi
    FALSE 986            LRP2      female        <NA>
    FALSE 987            LRP8        male        <NA>
    FALSE 988          LRPPRC      female          hi
    FALSE 989           LRRC1      female        <NA>
    FALSE 990          LRRC19        male        <NA>
    FALSE 991          LRRC3B        male        <NA>
    FALSE 992          LRRC3C      female        <NA>
    FALSE 993          LRRC70        male        <NA>
    FALSE 994          LRRIQ1      female          hi
    FALSE 995           LRRN1        male          lo
    FALSE 996           LRRN2        male          lo
    FALSE 997          LRRTM1      female        <NA>
    FALSE 998           LRTM2      female          lo
    FALSE 999           LSAMP      female        <NA>
    FALSE 1000          LSM12        male        <NA>
    FALSE 1001           LSP1        male        <NA>
    FALSE 1002          LTBP2      female        <NA>
    FALSE 1003            LUM        male        <NA>
    FALSE 1004          LYRM7        male        <NA>
    FALSE 1005         LYSMD3        male        <NA>
    FALSE 1006          MAGI3        male        <NA>
    FALSE 1007          MAK16      female          hi
    FALSE 1008         MAMDC2      female        <NA>
    FALSE 1009          MANBA      female          hi
    FALSE 1010         MANSC4      female          lo
    FALSE 1011           MAOA      female        <NA>
    FALSE 1012          MAP1B        male        <NA>
    FALSE 1013         MAP2K4      female        <NA>
    FALSE 1014         MAP3K1        male        <NA>
    FALSE 1015        MAP3K15        male        <NA>
    FALSE 1016         MAP7D3      female          hi
    FALSE 1017         MAPK13      female        <NA>
    FALSE 1018       MAPKAPK2        male        <NA>
    FALSE 1019         MARCH3        male        <NA>
    FALSE 1020       MARCKSL1      female        <NA>
    FALSE 1021       MARVELD2        male        <NA>
    FALSE 1022          MAST4        male        <NA>
    FALSE 1023         MBLAC2        male        <NA>
    FALSE 1024            MBP        male        <NA>
    FALSE 1025         MBTPS1      female        <NA>
    FALSE 1026            MCC        male        <NA>
    FALSE 1027          MCCC2        male        <NA>
    FALSE 1028          MCHR1      female        <NA>
    FALSE 1029            MDK      female        <NA>
    FALSE 1030          MEAF6      female        <NA>
    FALSE 1031          MECOM      female        <NA>
    FALSE 1032          MED10      female        <NA>
    FALSE 1033          MED18        male        <NA>
    FALSE 1034           MED8      female        <NA>
    FALSE 1035          MEF2D        male        <NA>
    FALSE 1036           MELK        male          hi
    FALSE 1037           MEST      female          hi
    FALSE 1038    Metazoa_SRP        male        <NA>
    FALSE 1039         METRNL      female          hi
    FALSE 1040        METTL24        male        <NA>
    FALSE 1041         METTL8      female        <NA>
    FALSE 1042         METTL9      female          hi
    FALSE 1043          MFAP2      female        <NA>
    FALSE 1044           MFI2      female        <NA>
    FALSE 1045           MFN1      female        <NA>
    FALSE 1046          MFSD1      female        <NA>
    FALSE 1047          MFSD4        male          hi
    FALSE 1048          MFSD6        male        <NA>
    FALSE 1049          MFSD7        male          lo
    FALSE 1050         MGAT4C      female          hi
    FALSE 1051          MGAT5        male        <NA>
    FALSE 1052           MIB2      female        <NA>
    FALSE 1053          MIER3      female        <NA>
    FALSE 1054          MKRN1      female        <NA>
    FALSE 1055          MLLT3        male        <NA>
    FALSE 1056          MMGT1      female        <NA>
    FALSE 1057           MMP2        male        <NA>
    FALSE 1058          MMRN2      female          lo
    FALSE 1059           MNX1      female        <NA>
    FALSE 1060          MOCS2        male        <NA>
    FALSE 1061          MORC3      female        <NA>
    FALSE 1062           MPDZ        male        <NA>
    FALSE 1063           MPND      female        <NA>
    FALSE 1064           MPP4      female        <NA>
    FALSE 1065         MPPED1        male        <NA>
    FALSE 1066           MREG      female        <NA>
    FALSE 1067         MRPL17        male        <NA>
    FALSE 1068          MRPL2      female        <NA>
    FALSE 1069         MRPL28      female        <NA>
    FALSE 1070         MRPL32      female        <NA>
    FALSE 1071         MRPL35      female        <NA>
    FALSE 1072         MRPL48      female          hi
    FALSE 1073         MRPL50        male        <NA>
    FALSE 1074         MRPS27        male        <NA>
    FALSE 1075         MRPS36      female        <NA>
    FALSE 1076           MSH3        male          hi
    FALSE 1077         MTHFD2        male        <NA>
    FALSE 1078          MTMR1        male          hi
    FALSE 1079         MTMR12        male        <NA>
    FALSE 1080          MTMR9        male        <NA>
    FALSE 1081          MTUS1      female        <NA>
    FALSE 1082           MTX3        male        <NA>
    FALSE 1083            MYC        male          hi
    FALSE 1084           MYCL      female        <NA>
    FALSE 1085          MYCT1      female          lo
    FALSE 1086          MYLIP        male        <NA>
    FALSE 1087           MYLK        male        <NA>
    FALSE 1088         MYO3AL        male        <NA>
    FALSE 1089          MYO7A        male          lo
    FALSE 1090          MYO7B        male        <NA>
    FALSE 1091          MYRIP      female        <NA>
    FALSE 1092          MYT1L        male        <NA>
    FALSE 1093          N4BP1      female        <NA>
    FALSE 1094         N6AMT2      female          hi
    FALSE 1095          NAA35        male        <NA>
    FALSE 1096       NAALADL2      female          hi
    FALSE 1097          NADK2        male        <NA>
    FALSE 1098        NADSYN1      female        <NA>
    FALSE 1099           NAGA      female        <NA>
    FALSE 1100           NANS        male        <NA>
    FALSE 1101           NARS        male        <NA>
    FALSE 1102           NAV2        male          lo
    FALSE 1103           NCAN      female        <NA>
    FALSE 1104          NCBP1        male        <NA>
    FALSE 1105          NCOA4      female        <NA>
    FALSE 1106           NCS1        male          hi
    FALSE 1107          NDRG1      female        <NA>
    FALSE 1108          NDRG4      female          lo
    FALSE 1109          NDST3        male        <NA>
    FALSE 1110          NDST4        male          hi
    FALSE 1111        NDUFA10      female        <NA>
    FALSE 1112        NDUFAF2        male        <NA>
    FALSE 1113         NDUFS4        male        <NA>
    FALSE 1114           NEBL        male          hi
    FALSE 1115         NECAB2      female          hi
    FALSE 1116         NECAP2      female        <NA>
    FALSE 1117          NELFA      female        <NA>
    FALSE 1118          NETO2        male        <NA>
    FALSE 1119           NFIB        male        <NA>
    FALSE 1120          NFIL3        male        <NA>
    FALSE 1121          NFKB1      female        <NA>
    FALSE 1122         NFKBIA      female          lo
    FALSE 1123           NFX1      female        <NA>
    FALSE 1124           NID1        male        <NA>
    FALSE 1125          NIPA2      female        <NA>
    FALSE 1126      NIPSNAP3A        male          lo
    FALSE 1127            NLN        male        <NA>
    FALSE 1128         NMNAT1      female        <NA>
    FALSE 1129          NMRK1        male        <NA>
    FALSE 1130          NOC3L      female        <NA>
    FALSE 1131            NOG      female        <NA>
    FALSE 1132           NOL6        male        <NA>
    FALSE 1133           NOL8      female          hi
    FALSE 1134           NOS3      female        <NA>
    FALSE 1135           NPC2      female          hi
    FALSE 1136          NPDC1        male        <NA>
    FALSE 1137            NPL      female        <NA>
    FALSE 1138          NR2F2        male        <NA>
    FALSE 1139          NR4A3        male          hi
    FALSE 1140          NRARP      female        <NA>
    FALSE 1141           NREP        male        <NA>
    FALSE 1142           NRG4      female        <NA>
    FALSE 1143          NRSN1      female          hi
    FALSE 1144           NSMF      female        <NA>
    FALSE 1145          NTAN1        male          lo
    FALSE 1146           NTF3      female        <NA>
    FALSE 1147           NTN4      female        <NA>
    FALSE 1148          NTRK2        male        <NA>
    FALSE 1149          NTRK3        male        <NA>
    FALSE 1150          NUBPL      female        <NA>
    FALSE 1151         NUDCD1      female          hi
    FALSE 1152         NUDT12        male        <NA>
    FALSE 1153          NUDT2        male        <NA>
    FALSE 1154         NUP155        male        <NA>
    FALSE 1155          NUP88        male        <NA>
    FALSE 1156          NUP93      female        <NA>
    FALSE 1157           NWD2        male        <NA>
    FALSE 1158          NXNL2        male        <NA>
    FALSE 1159          NXPH2        male        <NA>
    FALSE 1160           OAZ2      female        <NA>
    FALSE 1161           ODF2        male        <NA>
    FALSE 1162          OLFM3        male          lo
    FALSE 1163         OLFML1        male          lo
    FALSE 1164          OPRM1      female        <NA>
    FALSE 1165           ORC3      female        <NA>
    FALSE 1166         OSBPL9      female        <NA>
    FALSE 1167         OSGIN1        male        <NA>
    FALSE 1168           OSMR        male        <NA>
    FALSE 1169          OSTF1        male        <NA>
    FALSE 1170           OTOA      female        <NA>
    FALSE 1171          OXCT1        male        <NA>
    FALSE 1172         P2RY11      female        <NA>
    FALSE 1173          P4HA1        male        <NA>
    FALSE 1174       PAFAH1B1        male        <NA>
    FALSE 1175          PAIP1        male        <NA>
    FALSE 1176           PAK1        male          lo
    FALSE 1177           PAK6        male          hi
    FALSE 1178           PALM        male        <NA>
    FALSE 1179          PAMR1      female        <NA>
    FALSE 1180          PAPD4        male        <NA>
    FALSE 1181          PAPD7      female        <NA>
    FALSE 1182          PAPLN      female        <NA>
    FALSE 1183           PARG      female          hi
    FALSE 1184          PARK2      female        <NA>
    FALSE 1185          PARM1      female        <NA>
    FALSE 1186          PARP8        male        <NA>
    FALSE 1187           PAX6        male          lo
    FALSE 1188           PAX7        male          hi
    FALSE 1189           PCCB      female        <NA>
    FALSE 1190         PCDH10        male        <NA>
    FALSE 1191         PCDH19        male        <NA>
    FALSE 1192          PCF11      female          hi
    FALSE 1193          PCGF3        male        <NA>
    FALSE 1194         PCMTD2      female        <NA>
    FALSE 1195         PCYT1A      female        <NA>
    FALSE 1196         PDCD11      female          hi
    FALSE 1197          PDDC1        male        <NA>
    FALSE 1198          PDE3B        male          hi
    FALSE 1199          PDE4D        male          hi
    FALSE 1200          PDE8A        male          hi
    FALSE 1201          PDE8B        male          hi
    FALSE 1202          PDGFD      female          hi
    FALSE 1203         PDGFRA        male        <NA>
    FALSE 1204          PDIA6      female          hi
    FALSE 1205           PDK4      female          hi
    FALSE 1206         PDLIM3        male        <NA>
    FALSE 1207         PDLIM4        male          lo
    FALSE 1208           PDYN      female        <NA>
    FALSE 1209         PDZRN3      female        <NA>
    FALSE 1210         PECAM1      female        <NA>
    FALSE 1211           PECR      female        <NA>
    FALSE 1212          PELI2      female        <NA>
    FALSE 1213           PELO        male        <NA>
    FALSE 1214          PEX5L      female        <NA>
    FALSE 1215           PFKM        male        <NA>
    FALSE 1216          PGBD5      female        <NA>
    FALSE 1217           PGM3      female          hi
    FALSE 1218           PGM5        male          lo
    FALSE 1219            PGR      female        <NA>
    FALSE 1220         PGRMC1        male        <NA>
    FALSE 1221        PHACTR3        male        <NA>
    FALSE 1222           PHAX        male        <NA>
    FALSE 1223           PHC3      female        <NA>
    FALSE 1224          PHF19      female        <NA>
    FALSE 1225         PHF21B      female        <NA>
    FALSE 1226         PHLDA1        male        <NA>
    FALSE 1227          PIAS2        male        <NA>
    FALSE 1228           PIGG        male        <NA>
    FALSE 1229           PIGN      female        <NA>
    FALSE 1230           PIGO        male        <NA>
    FALSE 1231         PIK3C3        male          hi
    FALSE 1232        PIK3IP1      female        <NA>
    FALSE 1233         PIK3R2        male        <NA>
    FALSE 1234           PIM1        male          lo
    FALSE 1235           PIM3        male          hi
    FALSE 1236        PIP5K1B        male          lo
    FALSE 1237           PISD      female        <NA>
    FALSE 1238         PITPNA        male        <NA>
    FALSE 1239        PITPNM2        male        <NA>
    FALSE 1240           PJA2        male        <NA>
    FALSE 1241           PKIA      female          lo
    FALSE 1242         PKNOX2        male          lo
    FALSE 1243         PLA2R1      female        <NA>
    FALSE 1244           PLAA        male        <NA>
    FALSE 1245        PLAC8L1        male          lo
    FALSE 1246          PLCB1      female          hi
    FALSE 1247          PLCB2      female          hi
    FALSE 1248          PLCB4        male        <NA>
    FALSE 1249          PLCE1      female        <NA>
    FALSE 1250          PLCG1        male          lo
    FALSE 1251         PLCXD3        male        <NA>
    FALSE 1252        PLEKHA8        male        <NA>
    FALSE 1253        PLEKHF2        male        <NA>
    FALSE 1254        PLEKHG5      female          lo
    FALSE 1255        PLEKHO2        male          lo
    FALSE 1256          PLIN2        male        <NA>
    FALSE 1257           PLK2        male        <NA>
    FALSE 1258           PLLP        male        <NA>
    FALSE 1259          PLOD1      female        <NA>
    FALSE 1260          PLVAP      female        <NA>
    FALSE 1261         PM20D1      female        <NA>
    FALSE 1262         PM20D2      female        <NA>
    FALSE 1263           PMM1        male        <NA>
    FALSE 1264          PNLIP      female        <NA>
    FALSE 1265          PODXL      female          lo
    FALSE 1266          POLA1      female          hi
    FALSE 1267           POLK        male        <NA>
    FALSE 1268         POLR1C      female        <NA>
    FALSE 1269         POLR1E        male        <NA>
    FALSE 1270         POLR2H      female        <NA>
    FALSE 1271         POLR3F      female          hi
    FALSE 1272         POLR3G        male        <NA>
    FALSE 1273            POR      female        <NA>
    FALSE 1274         PPAP2A        male        <NA>
    FALSE 1275       PPAPDC1A        male        <NA>
    FALSE 1276          PPARD        male        <NA>
    FALSE 1277           PPIC        male          hi
    FALSE 1278          PPIL2      female        <NA>
    FALSE 1279            PPL      female        <NA>
    FALSE 1280       PPP1R14D      female        <NA>
    FALSE 1281        PPP1R1B        male          lo
    FALSE 1282        PPP1R1C      female        <NA>
    FALSE 1283        PPP1R9A      female        <NA>
    FALSE 1284         PPP2CA        male        <NA>
    FALSE 1285        PPP2R5C      female        <NA>
    FALSE 1286         PPP3CB      female        <NA>
    FALSE 1287          PPTC7      female          lo
    FALSE 1288          PPWD1        male        <NA>
    FALSE 1289          PQLC3      female        <NA>
    FALSE 1290           PRCP      female        <NA>
    FALSE 1291          PRDM1      female        <NA>
    FALSE 1292          PRDM5      female          hi
    FALSE 1293          PREPL      female        <NA>
    FALSE 1294         PRKAA1        male        <NA>
    FALSE 1295          PRKG2      female        <NA>
    FALSE 1296        PRKRIP1      female        <NA>
    FALSE 1297           PRLH      female        <NA>
    FALSE 1298          PRLHR      female        <NA>
    FALSE 1299         PROCA1      female          lo
    FALSE 1300          PROM1      female        <NA>
    FALSE 1301           PSAP      female        <NA>
    FALSE 1302          PSAT1        male        <NA>
    FALSE 1303           PSD3        male        <NA>
    FALSE 1304          PSIP1        male        <NA>
    FALSE 1305          PSMC3      female          hi
    FALSE 1306          PTAR1        male        <NA>
    FALSE 1307          PTCD2        male        <NA>
    FALSE 1308          PTCH1        male        <NA>
    FALSE 1309         PTCHD1      female        <NA>
    FALSE 1310          PTGES      female        <NA>
    FALSE 1311          PTGR1        male        <NA>
    FALSE 1312          PTHLH        male        <NA>
    FALSE 1313           PTK2        male        <NA>
    FALSE 1314           PTMS        male        <NA>
    FALSE 1315        PTPLAD2        male        <NA>
    FALSE 1316          PTPLB        male        <NA>
    FALSE 1317          PTPN5        male          hi
    FALSE 1318          PTPRR      female        <NA>
    FALSE 1319         PTPRZ1      female        <NA>
    FALSE 1320            PTS      female        <NA>
    FALSE 1321        PTTG1IP      female          hi
    FALSE 1322           PUM2      female        <NA>
    FALSE 1323          PYCRL      female        <NA>
    FALSE 1324           QARS      female          hi
    FALSE 1325           QPCT      female          hi
    FALSE 1326          QRFPR        male        <NA>
    FALSE 1327          QRSL1        male          hi
    FALSE 1328         RAB11A      female        <NA>
    FALSE 1329      RAB11FIP2        male        <NA>
    FALSE 1330          RAB20      female        <NA>
    FALSE 1331          RAB36        male        <NA>
    FALSE 1332          RAB3A        male        <NA>
    FALSE 1333          RAB3C        male        <NA>
    FALSE 1334          RAB6A      female        <NA>
    FALSE 1335        RABGAP1      female        <NA>
    FALSE 1336           RAD1        male        <NA>
    FALSE 1337          RAD17        male        <NA>
    FALSE 1338         RAD23B        male          lo
    FALSE 1339          RAI14        male        <NA>
    FALSE 1340           RAI2      female        <NA>
    FALSE 1341        RALGAPB        male        <NA>
    FALSE 1342           RALY        male        <NA>
    FALSE 1343          RAMP2      female        <NA>
    FALSE 1344        RAPGEF1        male        <NA>
    FALSE 1345        RAPGEF3      female          lo
    FALSE 1346           RARS      female          hi
    FALSE 1347          RASA1        male        <NA>
    FALSE 1348          RASA4      female          hi
    FALSE 1349         RASAL1      female          hi
    FALSE 1350          RASD1        male        <NA>
    FALSE 1351          RASD2      female        <NA>
    FALSE 1352          RASEF        male        <NA>
    FALSE 1353        RASGRF1        male          lo
    FALSE 1354        RASGRP1        male        <NA>
    FALSE 1355        RASL10A      female          hi
    FALSE 1356        RASL11A      female          hi
    FALSE 1357        RASL11B      female        <NA>
    FALSE 1358         RASSF2        male        <NA>
    FALSE 1359         RASSF5        male          hi
    FALSE 1360         RBFOX1      female          lo
    FALSE 1361          RBM25      female          hi
    FALSE 1362          RBP4B        male        <NA>
    FALSE 1363          RBPMS        male          lo
    FALSE 1364          RCAN2        male          hi
    FALSE 1365         RCBTB2        male          lo
    FALSE 1366           RCL1        male        <NA>
    FALSE 1367          RCOR1      female        <NA>
    FALSE 1368            RDX      female          lo
    FALSE 1369         RECQL5      female        <NA>
    FALSE 1370          REEP1      female          hi
    FALSE 1371          REEP6      female        <NA>
    FALSE 1372           REG4      female          lo
    FALSE 1373          RELL1      female        <NA>
    FALSE 1374          REPS2        male        <NA>
    FALSE 1375          RERGL      female        <NA>
    FALSE 1376         RETSAT      female        <NA>
    FALSE 1377           RFC1        male        <NA>
    FALSE 1378          RFESD        male        <NA>
    FALSE 1379            RFK        male        <NA>
    FALSE 1380           RFNG        male        <NA>
    FALSE 1381           RFX8      female        <NA>
    FALSE 1382           RGMB        male        <NA>
    FALSE 1383           RGP1        male          lo
    FALSE 1384           RGS1        male        <NA>
    FALSE 1385           RGS2        male        <NA>
    FALSE 1386           RGS5      female          lo
    FALSE 1387           RGS6      female        <NA>
    FALSE 1388         RGS7BP        male        <NA>
    FALSE 1389         RHBDL3        male        <NA>
    FALSE 1390        RHOBTB3        male          hi
    FALSE 1391           RHOJ      female        <NA>
    FALSE 1392           RHOU        male          lo
    FALSE 1393           RIC1        male        <NA>
    FALSE 1394         RICTOR        male        <NA>
    FALSE 1395           RIN3      female          lo
    FALSE 1396          RIOK2        male        <NA>
    FALSE 1397          RIPK3      female          hi
    FALSE 1398           RIT2        male        <NA>
    FALSE 1399           RLN1        male        <NA>
    FALSE 1400           RMI1        male        <NA>
    FALSE 1401         RMND5A      female        <NA>
    FALSE 1402          RNF13      female        <NA>
    FALSE 1403         RNF165        male          lo
    FALSE 1404         RNF170        male        <NA>
    FALSE 1405         RNF180        male        <NA>
    FALSE 1406         RNF19A      female        <NA>
    FALSE 1407         RNF19B        male        <NA>
    FALSE 1408           RNF2        male        <NA>
    FALSE 1409          RNF20        male        <NA>
    FALSE 1410         RNF207      female        <NA>
    FALSE 1411          RNF32      female          hi
    FALSE 1412          RNF38        male        <NA>
    FALSE 1413           RNF4        male        <NA>
    FALSE 1414          RNPEP      female        <NA>
    FALSE 1415          ROBO4      female        <NA>
    FALSE 1416           ROR1        male          lo
    FALSE 1417   RP11-145E5.5        male        <NA>
    FALSE 1418 RP11-195F19.29      female          lo
    FALSE 1419  RP11-296A16.1        male          hi
    FALSE 1420   RP11-505K9.4      female        <NA>
    FALSE 1421  RP11-561B11.2      female        <NA>
    FALSE 1422 RP11-574K11.31        male          hi
    FALSE 1423  RP11-598P20.5        male        <NA>
    FALSE 1424           RPA1        male          hi
    FALSE 1425           RPIA      female        <NA>
    FALSE 1426         RPL10A      female        <NA>
    FALSE 1427          RPL12      female        <NA>
    FALSE 1428          RPL15      female          hi
    FALSE 1429          RPL17      female        <NA>
    FALSE 1430        RPL22L1      female        <NA>
    FALSE 1431         RPL27A      female          hi
    FALSE 1432           RPL3      female        <NA>
    FALSE 1433          RPL30      female          hi
    FALSE 1434          RPL37        male        <NA>
    FALSE 1435          RPL39      female          hi
    FALSE 1436           RPL4      female          hi
    FALSE 1437           RPL5      female          hi
    FALSE 1438           RPL6      female          hi
    FALSE 1439          RPL7A      female        <NA>
    FALSE 1440           RPL8      female        <NA>
    FALSE 1441          RPLP0      female        <NA>
    FALSE 1442          RPLP1      female        <NA>
    FALSE 1443           RPN2      female          hi
    FALSE 1444          RPP38      female        <NA>
    FALSE 1445          RPS12      female        <NA>
    FALSE 1446           RPS2      female          hi
    FALSE 1447          RPS23        male        <NA>
    FALSE 1448          RPS24      female        <NA>
    FALSE 1449          RPS25      female        <NA>
    FALSE 1450         RPS27A      female          hi
    FALSE 1451           RPS3      female        <NA>
    FALSE 1452          RPS3A      female          hi
    FALSE 1453          RPS4X      female        <NA>
    FALSE 1454           RPS6        male        <NA>
    FALSE 1455        RPS6KB2        male        <NA>
    FALSE 1456           RPS7      female          hi
    FALSE 1457           RPS8      female          hi
    FALSE 1458           RPSA      female          hi
    FALSE 1459          RRBP1      female          hi
    FALSE 1460          RSPO3      female        <NA>
    FALSE 1461        RTN4RL1      female        <NA>
    FALSE 1462          RUSC2        male        <NA>
    FALSE 1463          SALL3        male        <NA>
    FALSE 1464          SAP30      female        <NA>
    FALSE 1465          SAR1B      female        <NA>
    FALSE 1466          SARM1      female        <NA>
    FALSE 1467          SATB2        male        <NA>
    FALSE 1468         SCAMP1        male        <NA>
    FALSE 1469         SCAMP2      female        <NA>
    FALSE 1470            SCD        male        <NA>
    FALSE 1471          SCFD2      female          hi
    FALSE 1472          SCML2      female        <NA>
    FALSE 1473          SCN3B        male        <NA>
    FALSE 1474         SCNN1D      female        <NA>
    FALSE 1475           SCP2      female        <NA>
    FALSE 1476            SCT        male        <NA>
    FALSE 1477          SDAD1      female        <NA>
    FALSE 1478           SDC4      female        <NA>
    FALSE 1479           SDHA      female        <NA>
    FALSE 1480           SDPR      female        <NA>
    FALSE 1481          SEBOX      female        <NA>
    FALSE 1482        SEC14L5        male          hi
    FALSE 1483          SEC63      female          hi
    FALSE 1484       SECISBP2        male        <NA>
    FALSE 1485          SEL1L      female          hi
    FALSE 1486       SELENBP1      female        <NA>
    FALSE 1487         SEMA3A      female          hi
    FALSE 1488         SEMA4D        male        <NA>
    FALSE 1489         SEMA6A        male        <NA>
    FALSE 1490         SEPT12      female        <NA>
    FALSE 1491          SEPT6      female        <NA>
    FALSE 1492          SEPT9        male        <NA>
    FALSE 1493         SERHL2      female        <NA>
    FALSE 1494       SERPINC1      female        <NA>
    FALSE 1495         SERTM1        male        <NA>
    FALSE 1496         SETBP1        male        <NA>
    FALSE 1497          SETD9        male        <NA>
    FALSE 1498          SEZ6L      female        <NA>
    FALSE 1499          SFRP2        male        <NA>
    FALSE 1500         SFTPA1        male        <NA>
    FALSE 1501           SGCE      female        <NA>
    FALSE 1502           SGK1      female        <NA>
    FALSE 1503           SGK2      female        <NA>
    FALSE 1504          SGMS2      female        <NA>
    FALSE 1505           SGTB        male        <NA>
    FALSE 1506         SH2D3C      female        <NA>
    FALSE 1507         SH2D4A      female        <NA>
    FALSE 1508       SH3BGRL2        male        <NA>
    FALSE 1509         SH3GL1      female        <NA>
    FALSE 1510         SH3GL2      female        <NA>
    FALSE 1511       SH3PXD2A      female        <NA>
    FALSE 1512            SHB        male        <NA>
    FALSE 1513           SHC3        male          hi
    FALSE 1514            SHE      female          lo
    FALSE 1515         SHISA2      female        <NA>
    FALSE 1516         SHISA9      female        <NA>
    FALSE 1517        SIGMAR1        male        <NA>
    FALSE 1518           SIK1        male        <NA>
    FALSE 1519           SIK2        male        <NA>
    FALSE 1520           SIL1      female        <NA>
    FALSE 1521           SIX6        male        <NA>
    FALSE 1522           SKA2        male        <NA>
    FALSE 1523        SKIV2L2        male        <NA>
    FALSE 1524          SKOR1      female        <NA>
    FALSE 1525           SKP2        male        <NA>
    FALSE 1526        SLC12A2        male        <NA>
    FALSE 1527        SLC13A2      female        <NA>
    FALSE 1528        SLC14A2        male        <NA>
    FALSE 1529        SLC15A1        male        <NA>
    FALSE 1530       SLC16A14      female          hi
    FALSE 1531        SLC16A3        male        <NA>
    FALSE 1532        SLC16A6        male        <NA>
    FALSE 1533        SLC18A3      female        <NA>
    FALSE 1534        SLC19A3      female        <NA>
    FALSE 1535       SLC25A20      female          hi
    FALSE 1536       SLC25A25        male        <NA>
    FALSE 1537       SLC25A29      female          hi
    FALSE 1538       SLC25A30      female        <NA>
    FALSE 1539       SLC25A46        male        <NA>
    FALSE 1540       SLC25A47      female        <NA>
    FALSE 1541       SLC25A51        male        <NA>
    FALSE 1542        SLC27A4      female        <NA>
    FALSE 1543        SLC29A4        male          lo
    FALSE 1544        SLC2A13        male          lo
    FALSE 1545         SLC2A9        male        <NA>
    FALSE 1546        SLC30A2      female        <NA>
    FALSE 1547        SLC30A5        male          lo
    FALSE 1548        SLC31A2      female        <NA>
    FALSE 1549        SLC32A1      female        <NA>
    FALSE 1550        SLC35C1      female        <NA>
    FALSE 1551        SLC35D2        male        <NA>
    FALSE 1552        SLC35F5      female          hi
    FALSE 1553        SLC38A8      female        <NA>
    FALSE 1554        SLC38A9        male        <NA>
    FALSE 1555        SLC39A6      female        <NA>
    FALSE 1556        SLC41A3      female          hi
    FALSE 1557        SLC44A1        male        <NA>
    FALSE 1558        SLC44A3      female        <NA>
    FALSE 1559         SLC5A6      female        <NA>
    FALSE 1560         SLC6A2      female        <NA>
    FALSE 1561         SLC6A4      female        <NA>
    FALSE 1562        SLC7A10      female        <NA>
    FALSE 1563         SLC8A2        male          lo
    FALSE 1564         SLC9A5      female        <NA>
    FALSE 1565         SLC9A8      female        <NA>
    FALSE 1566        SLCO4C1        male        <NA>
    FALSE 1567          SLIRP      female        <NA>
    FALSE 1568        SLITRK5        male        <NA>
    FALSE 1569          SMAD2      female          lo
    FALSE 1570          SMAD7        male        <NA>
    FALSE 1571        SMARCA2        male        <NA>
    FALSE 1572        SMARCA5      female          hi
    FALSE 1573        SMARCD3      female        <NA>
    FALSE 1574           SMC2        male          hi
    FALSE 1575           SMC5        male          hi
    FALSE 1576         SMIM15        male        <NA>
    FALSE 1577            SMN        male        <NA>
    FALSE 1578           SMU1        male        <NA>
    FALSE 1579          SMYD1      female          hi
    FALSE 1580          SMYD3      female        <NA>
    FALSE 1581          SNAI2      female        <NA>
    FALSE 1582         SNAPC3        male        <NA>
    FALSE 1583         SNCAIP        male        <NA>
    FALSE 1584           SNRK        male        <NA>
    FALSE 1585        SNRNP48      female        <NA>
    FALSE 1586          SNTB2      female        <NA>
    FALSE 1587           SNW1      female        <NA>
    FALSE 1588           SNX2        male        <NA>
    FALSE 1589          SNX24        male        <NA>
    FALSE 1590          SNX30        male        <NA>
    FALSE 1591            SON      female        <NA>
    FALSE 1592         SORCS1      female        <NA>
    FALSE 1593           SOUL        male          lo
    FALSE 1594          SOX18      female          lo
    FALSE 1595           SOX7      female        <NA>
    FALSE 1596        SPATA17      female        <NA>
    FALSE 1597          SPEF2        male        <NA>
    FALSE 1598          SPERT        male        <NA>
    FALSE 1599           SPG7      female        <NA>
    FALSE 1600          SPON1        male        <NA>
    FALSE 1601          SPON2        male          lo
    FALSE 1602          SPRY1        male          hi
    FALSE 1603          SPRY2      female        <NA>
    FALSE 1604         SPRYD7      female        <NA>
    FALSE 1605         SPTLC1        male        <NA>
    FALSE 1606        SPTY2D1      female          hi
    FALSE 1607          SREK1        male        <NA>
    FALSE 1608       SREK1IP1        male        <NA>
    FALSE 1609            SRF        male          lo
    FALSE 1610         SRFBP1        male        <NA>
    FALSE 1611         SRGAP3        male          lo
    FALSE 1612          SRRM4        male        <NA>
    FALSE 1613          SRSF4        male        <NA>
    FALSE 1614          SSBP2        male        <NA>
    FALSE 1615          SSTR1      female          hi
    FALSE 1616           ST13      female          hi
    FALSE 1617            ST5      female        <NA>
    FALSE 1618     ST6GALNAC5      female          hi
    FALSE 1619        ST8SIA2        male          lo
    FALSE 1620        ST8SIA3        male        <NA>
    FALSE 1621        ST8SIA4        male        <NA>
    FALSE 1622        ST8SIA5        male        <NA>
    FALSE 1623        STARD13      female        <NA>
    FALSE 1624         STARD4        male        <NA>
    FALSE 1625          STAT3        male          lo
    FALSE 1626          STIM2      female          hi
    FALSE 1627         STK17A        male          lo
    FALSE 1628          STK40        male        <NA>
    FALSE 1629         STRADA        male          lo
    FALSE 1630          STT3B      female          hi
    FALSE 1631         STXBP4      female        <NA>
    FALSE 1632          SULF1      female        <NA>
    FALSE 1633         SULT1B      female        <NA>
    FALSE 1634        SULT6B1        male        <NA>
    FALSE 1635          SUMF1      female        <NA>
    FALSE 1636          SUMF2      female          hi
    FALSE 1637          SUSD1        male        <NA>
    FALSE 1638          SUSD3      female        <NA>
    FALSE 1639        SUV39H2      female        <NA>
    FALSE 1640           SV2B        male          hi
    FALSE 1641          SVEP1        male        <NA>
    FALSE 1642            SYK        male          lo
    FALSE 1643           SYN3        male          hi
    FALSE 1644       SYNDIG1L        male        <NA>
    FALSE 1645          SYNE2      female          hi
    FALSE 1646         SYNGR3        male        <NA>
    FALSE 1647          SYPL2      female        <NA>
    FALSE 1648          SYT14      female          hi
    FALSE 1649          SYT15        male          lo
    FALSE 1650          SYT17        male        <NA>
    FALSE 1651           SYT4        male          hi
    FALSE 1652          TACC2      female        <NA>
    FALSE 1653          TACR1      female        <NA>
    FALSE 1654          TAF1B      female        <NA>
    FALSE 1655          TAGLN        male          lo
    FALSE 1656           TARS        male        <NA>
    FALSE 1657        TBC1D16      female        <NA>
    FALSE 1658        TBC1D30        male          hi
    FALSE 1659           TBCA        male        <NA>
    FALSE 1660           TBL2      female          hi
    FALSE 1661          TBX19        male        <NA>
    FALSE 1662           TCN2      female        <NA>
    FALSE 1663          TCTE3      female        <NA>
    FALSE 1664            TDG        male        <NA>
    FALSE 1665           TDP1      female        <NA>
    FALSE 1666          TDRD7        male        <NA>
    FALSE 1667            TEC        male        <NA>
    FALSE 1668          TECRL      female        <NA>
    FALSE 1669           TERT      female        <NA>
    FALSE 1670          TESK1        male          lo
    FALSE 1671           TET1      female          hi
    FALSE 1672          TEX10      female          hi
    FALSE 1673          TEX11      female        <NA>
    FALSE 1674          TEX12        male        <NA>
    FALSE 1675          TEX30      female        <NA>
    FALSE 1676          TFB2M      female          hi
    FALSE 1677          TGFBI        male          lo
    FALSE 1678          TGIF1      female        <NA>
    FALSE 1679           TGM3      female        <NA>
    FALSE 1680             TH        male          lo
    FALSE 1681          THAP1        male        <NA>
    FALSE 1682           THBD      female        <NA>
    FALSE 1683          THBS4        male          lo
    FALSE 1684           THRB        male          hi
    FALSE 1685           TIE1      female          lo
    FALSE 1686        TINAGL1      female          lo
    FALSE 1687         TIPARP      female        <NA>
    FALSE 1688           TJP2        male        <NA>
    FALSE 1689           TLN2        male          lo
    FALSE 1690        TM4SF18      female        <NA>
    FALSE 1691         TM6SF1      female        <NA>
    FALSE 1692         TMBIM1      female        <NA>
    FALSE 1693          TMCC2        male          lo
    FALSE 1694          TMCC3      female        <NA>
    FALSE 1695          TMCO1      female          hi
    FALSE 1696          TMED7      female          hi
    FALSE 1697        TMEM100        male        <NA>
    FALSE 1698       TMEM106B      female        <NA>
    FALSE 1699        TMEM117        male        <NA>
    FALSE 1700        TMEM123        male        <NA>
    FALSE 1701       TMEM132B        male        <NA>
    FALSE 1702       TMEM161B        male        <NA>
    FALSE 1703       TMEM167A        male          hi
    FALSE 1704        TMEM171        male        <NA>
    FALSE 1705        TMEM175        male        <NA>
    FALSE 1706        TMEM177      female        <NA>
    FALSE 1707          TMEM2        male        <NA>
    FALSE 1708       TMEM200A        male          lo
    FALSE 1709        TMEM204      female        <NA>
    FALSE 1710        TMEM209      female        <NA>
    FALSE 1711        TMEM246        male        <NA>
    FALSE 1712       TMEM255B      female          hi
    FALSE 1713        TMEM258      female          hi
    FALSE 1714        TMEM38B        male        <NA>
    FALSE 1715         TMEM47      female        <NA>
    FALSE 1716        TMEM50B      female        <NA>
    FALSE 1717         TMEM51      female        <NA>
    FALSE 1718         TMEM61        male        <NA>
    FALSE 1719         TMEM62        male        <NA>
    FALSE 1720         TMEM64      female        <NA>
    FALSE 1721         TMEM65        male        <NA>
    FALSE 1722         TMEM69      female        <NA>
    FALSE 1723         TMEM8A      female        <NA>
    FALSE 1724         TMEM8B        male        <NA>
    FALSE 1725          TMLHE      female        <NA>
    FALSE 1726          TMOD1        male        <NA>
    FALSE 1727      TMPRSS11F      female        <NA>
    FALSE 1728        TMPRSS3        male        <NA>
    FALSE 1729          TMTC1      female          lo
    FALSE 1730        TNFRSF9      female        <NA>
    FALSE 1731           TNIK        male          lo
    FALSE 1732           TNK2        male          lo
    FALSE 1733          TNPO1        male          hi
    FALSE 1734         TNRC18        male        <NA>
    FALSE 1735           TOB1        male        <NA>
    FALSE 1736         TOLLIP        male        <NA>
    FALSE 1737         TOMM22      female        <NA>
    FALSE 1738          TOMM5        male        <NA>
    FALSE 1739         TOPORS        male        <NA>
    FALSE 1740           TOX2        male        <NA>
    FALSE 1741         TP53I3      female        <NA>
    FALSE 1742          TPCN1        male        <NA>
    FALSE 1743          TPGS2        male        <NA>
    FALSE 1744           TPI1      female        <NA>
    FALSE 1745           TPM2        male        <NA>
    FALSE 1746           TPMT      female        <NA>
    FALSE 1747          TPST2      female        <NA>
    FALSE 1748          TPTE2      female        <NA>
    FALSE 1749          TRAP1      female        <NA>
    FALSE 1750       TRAPPC13        male        <NA>
    FALSE 1751          TRAT1      female        <NA>
    FALSE 1752          TRIB1        male        <NA>
    FALSE 1753         TRIM23        male        <NA>
    FALSE 1754         TRIM36        male        <NA>
    FALSE 1755         TRIM54      female        <NA>
    FALSE 1756         TRIM66      female        <NA>
    FALSE 1757         TRIM71      female        <NA>
    FALSE 1758          TRIP4      female        <NA>
    FALSE 1759        TRMT10B        male        <NA>
    FALSE 1760        TRMT61A      female        <NA>
    FALSE 1761         TROJAN        male        <NA>
    FALSE 1762          TRPA1        male          hi
    FALSE 1763          TRPC4      female        <NA>
    FALSE 1764          TRPC5      female        <NA>
    FALSE 1765          TRPC6      female        <NA>
    FALSE 1766          TRPM2      female        <NA>
    FALSE 1767          TRPM3        male        <NA>
    FALSE 1768          TRPM6        male        <NA>
    FALSE 1769          TRPV6      female          hi
    FALSE 1770        TSC22D1      female        <NA>
    FALSE 1771          TSEN2      female        <NA>
    FALSE 1772         TSGA10      female        <NA>
    FALSE 1773         TSPAN6      female        <NA>
    FALSE 1774          TSTD2        male        <NA>
    FALSE 1775          TTC19      female        <NA>
    FALSE 1776          TTC33        male        <NA>
    FALSE 1777          TTC37        male        <NA>
    FALSE 1778         TTC39B        male        <NA>
    FALSE 1779          TTC7B      female        <NA>
    FALSE 1780           TTF2        male        <NA>
    FALSE 1781          TULP3      female          hi
    FALSE 1782          TUSC1        male        <NA>
    FALSE 1783           TWF2      female          lo
    FALSE 1784          TXNIP      female        <NA>
    FALSE 1785          TYRP1        male          lo
    FALSE 1786           UACA        male          lo
    FALSE 1787           UBA6      female        <NA>
    FALSE 1788         UBE2R2        male        <NA>
    FALSE 1789          UBE3C      female        <NA>
    FALSE 1790          UBE4B      female        <NA>
    FALSE 1791         UBQLN1        male        <NA>
    FALSE 1792         UBXN2B      female        <NA>
    FALSE 1793           UCK1        male        <NA>
    FALSE 1794           UGCG        male        <NA>
    FALSE 1795           UMOD      female        <NA>
    FALSE 1796           UMPS      female        <NA>
    FALSE 1797         UNC119        male          hi
    FALSE 1798         UNC13C        male          hi
    FALSE 1799           UNKL      female          hi
    FALSE 1800          UPK1B      female        <NA>
    FALSE 1801         UQCRC2      female          hi
    FALSE 1802          UQCRQ      female        <NA>
    FALSE 1803          UROC1      female        <NA>
    FALSE 1804          USH1C      female        <NA>
    FALSE 1805           USO1      female          hi
    FALSE 1806          USP12      female        <NA>
    FALSE 1807          USP15      female        <NA>
    FALSE 1808          USP35        male          hi
    FALSE 1809           USP4      female        <NA>
    FALSE 1810          USP46      female        <NA>
    FALSE 1811          UTP15        male        <NA>
    FALSE 1812           VAPA      female        <NA>
    FALSE 1813          VASH2      female        <NA>
    FALSE 1814           VASN      female        <NA>
    FALSE 1815          VAT1L        male        <NA>
    FALSE 1816            VCP      female          hi
    FALSE 1817          VDAC3        male        <NA>
    FALSE 1818          VEGFA      female        <NA>
    FALSE 1819          VEGFC      female        <NA>
    FALSE 1820          VGLL3      female        <NA>
    FALSE 1821        VIPAS39      female        <NA>
    FALSE 1822          VLDLR        male        <NA>
    FALSE 1823          VOPP1      female        <NA>
    FALSE 1824         VPS37B        male        <NA>
    FALSE 1825          VPS41      female        <NA>
    FALSE 1826           VRK3      female        <NA>
    FALSE 1827            VTN      female        <NA>
    FALSE 1828         VWA5B2        male        <NA>
    FALSE 1829           VWA8      female        <NA>
    FALSE 1830          VWC2L        male        <NA>
    FALSE 1831          WASF1        male          lo
    FALSE 1832         WBP2NL      female        <NA>
    FALSE 1833          WDPCP      female        <NA>
    FALSE 1834           WDR1        male        <NA>
    FALSE 1835          WDR19      female        <NA>
    FALSE 1836           WDR3      female          hi
    FALSE 1837          WDR36        male        <NA>
    FALSE 1838          WDR41        male        <NA>
    FALSE 1839          WDR70        male        <NA>
    FALSE 1840          WDR74      female        <NA>
    FALSE 1841          WDR90        male        <NA>
    FALSE 1842        WFIKKN1      female        <NA>
    FALSE 1843           WFS1      female        <NA>
    FALSE 1844           WIF1        male        <NA>
    FALSE 1845           WNT4        male        <NA>
    FALSE 1846          WNT5A      female        <NA>
    FALSE 1847          WNT9A        male        <NA>
    FALSE 1848           WSB2      female          hi
    FALSE 1849          WSCD2      female          hi
    FALSE 1850           XKRX      female          lo
    FALSE 1851            XPA        male        <NA>
    FALSE 1852           XPOT      female          hi
    FALSE 1853          XRCC4        male        <NA>
    FALSE 1854           XYLB      female        <NA>
    FALSE 1855          YIF1B        male        <NA>
    FALSE 1856         YTHDC2        male        <NA>
    FALSE 1857            YY1      female        <NA>
    FALSE 1858          ZBED1        male        <NA>
    FALSE 1859         ZBTB49      female        <NA>
    FALSE 1860          ZBTB5      female        <NA>
    FALSE 1861         ZBTB7C        male        <NA>
    FALSE 1862       ZC3HAV1L      female          hi
    FALSE 1863         ZCCHC6        male        <NA>
    FALSE 1864         ZCCHC7        male        <NA>
    FALSE 1865         ZCCHC9        male        <NA>
    FALSE 1866        ZDHHC21        male        <NA>
    FALSE 1867         ZDHHC6      female        <NA>
    FALSE 1868         ZFAND1      female        <NA>
    FALSE 1869         ZFAND5      female        <NA>
    FALSE 1870          ZFPM1      female        <NA>
    FALSE 1871            ZFR      female        <NA>
    FALSE 1872        ZFYVE16        male        <NA>
    FALSE 1873        ZFYVE21      female        <NA>
    FALSE 1874          ZMYM2      female        <NA>
    FALSE 1875         ZNF131      female        <NA>
    FALSE 1876         ZNF346      female        <NA>
    FALSE 1877         ZNF367        male          hi
    FALSE 1878         ZNF568        male          lo
    FALSE 1879         ZNF608        male          lo
    FALSE 1880         ZNF703      female        <NA>
    FALSE 1881         ZSWIM6        male        <NA>
    FALSE 1882         ZSWIM8      female        <NA>
    FALSE 1883          ZUFSP      female        <NA>
    FALSE 1884          ABCA2        <NA>          lo
    FALSE 1885          ABCC3        <NA>          hi
    FALSE 1886          ABCC8        <NA>          hi
    FALSE 1887          ABCC9        <NA>          hi
    FALSE 1888          ABCE1        <NA>          hi
    FALSE 1889          ABHD3        <NA>          lo
    FALSE 1890          ABHD8        <NA>          lo
    FALSE 1891           ABRA        <NA>          hi
    FALSE 1892     AC025048.1        <NA>          hi
    FALSE 1893          ACAD9        <NA>          hi
    FALSE 1894           ACCS        <NA>          lo
    FALSE 1895           ACHE        <NA>          lo
    FALSE 1896           ACRC        <NA>          lo
    FALSE 1897         ACSBG2        <NA>          hi
    FALSE 1898          ACSL6        <NA>          hi
    FALSE 1899          ACSS2        <NA>          hi
    FALSE 1900          ACTN4        <NA>          lo
    FALSE 1901         ACTR10        <NA>          hi
    FALSE 1902          ACTR2        <NA>          hi
    FALSE 1903         ACVR1B        <NA>          lo
    FALSE 1904         ACVR2B        <NA>          hi
    FALSE 1905         ACVRL1        <NA>          lo
    FALSE 1906          ACYP2        <NA>          hi
    FALSE 1907         ADAM17        <NA>          lo
    FALSE 1908         ADAM19        <NA>          lo
    FALSE 1909         ADAM33        <NA>          lo
    FALSE 1910          ADAM5        <NA>          lo
    FALSE 1911       ADAMTS15        <NA>          hi
    FALSE 1912        ADAMTS7        <NA>          lo
    FALSE 1913        ADAMTS8        <NA>          hi
    FALSE 1914        ADAMTS9        <NA>          lo
    FALSE 1915       ADAMTSL1        <NA>          hi
    FALSE 1916       ADAMTSL5        <NA>          lo
    FALSE 1917          ADAP1        <NA>          lo
    FALSE 1918          ADAP2        <NA>          lo
    FALSE 1919         ADARB1        <NA>          hi
    FALSE 1920         ADARB2        <NA>          hi
    FALSE 1921          ADCK3        <NA>          lo
    FALSE 1922            ADK        <NA>          hi
    FALSE 1923         ADRA1D        <NA>          lo
    FALSE 1924         ADRA2C        <NA>          lo
    FALSE 1925           ADSL        <NA>          hi
    FALSE 1926           ADSS        <NA>          hi
    FALSE 1927        AFAP1L2        <NA>          lo
    FALSE 1928           AFF2        <NA>          hi
    FALSE 1929          AGBL5        <NA>          lo
    FALSE 1930           AGRN        <NA>          lo
    FALSE 1931         AHNAK2        <NA>          lo
    FALSE 1932           AIDA        <NA>          hi
    FALSE 1933          AIFM1        <NA>          hi
    FALSE 1934          AIMP2        <NA>          hi
    FALSE 1935            AK5        <NA>          lo
    FALSE 1936          AKAP6        <NA>          hi
    FALSE 1937          AKAP9        <NA>          hi
    FALSE 1938        AKIRIN2        <NA>          lo
    FALSE 1939           AKT2        <NA>          lo
    FALSE 1940     AL928654.7        <NA>          lo
    FALSE 1941       ALDH18A1        <NA>          hi
    FALSE 1942        ALDH1L2        <NA>          hi
    FALSE 1943        ALDH4A1        <NA>          lo
    FALSE 1944        ALDH6A1        <NA>          hi
    FALSE 1945           ALG1        <NA>          hi
    FALSE 1946          ALG12        <NA>          hi
    FALSE 1947           ALG2        <NA>          hi
    FALSE 1948           ALG5        <NA>          hi
    FALSE 1949           ALG8        <NA>          hi
    FALSE 1950         ALKBH8        <NA>          hi
    FALSE 1951          ALOX5        <NA>          lo
    FALSE 1952         AMBRA1        <NA>          lo
    FALSE 1953           AMFR        <NA>          lo
    FALSE 1954          AMHR2        <NA>          lo
    FALSE 1955         AMIGO2        <NA>          hi
    FALSE 1956         ANAPC4        <NA>          hi
    FALSE 1957         ANGPT2        <NA>          lo
    FALSE 1958        ANGPTL2        <NA>          lo
    FALSE 1959        ANKRD12        <NA>          hi
    FALSE 1960       ANKRD13B        <NA>          lo
    FALSE 1961        ANKRD42        <NA>          hi
    FALSE 1962        ANKRD47        <NA>          lo
    FALSE 1963         ANKRD6        <NA>          lo
    FALSE 1964         ANKRD9        <NA>          lo
    FALSE 1965           ANLN        <NA>          hi
    FALSE 1966           ANO3        <NA>          hi
    FALSE 1967         ANTXRL        <NA>          lo
    FALSE 1968           AOC1        <NA>          lo
    FALSE 1969    AP000350.10        <NA>          lo
    FALSE 1970          AP1AR        <NA>          hi
    FALSE 1971          AP1S3        <NA>          hi
    FALSE 1972          APAF1        <NA>          hi
    FALSE 1973          APBA3        <NA>          lo
    FALSE 1974           APC2        <NA>          lo
    FALSE 1975           API5        <NA>          hi
    FALSE 1976            APP        <NA>          hi
    FALSE 1977          AQP11        <NA>          hi
    FALSE 1978          ARAP2        <NA>          hi
    FALSE 1979           ARF1        <NA>          hi
    FALSE 1980           ARF4        <NA>          hi
    FALSE 1981        ARFGAP1        <NA>          hi
    FALSE 1982        ARFGAP3        <NA>          hi
    FALSE 1983        ARFGEF3        <NA>          hi
    FALSE 1984         ARFIP2        <NA>          lo
    FALSE 1985         ARGLU1        <NA>          hi
    FALSE 1986        ARHGAP1        <NA>          lo
    FALSE 1987       ARHGAP15        <NA>          lo
    FALSE 1988       ARHGAP19        <NA>          hi
    FALSE 1989       ARHGAP20        <NA>          lo
    FALSE 1990       ARHGAP25        <NA>          lo
    FALSE 1991       ARHGAP26        <NA>          lo
    FALSE 1992       ARHGAP27        <NA>          lo
    FALSE 1993       ARHGAP32        <NA>          lo
    FALSE 1994       ARHGAP44        <NA>          lo
    FALSE 1995        ARHGDIG        <NA>          lo
    FALSE 1996      ARHGEF10L        <NA>          lo
    FALSE 1997        ARHGEF3        <NA>          hi
    FALSE 1998         ARID1A        <NA>          lo
    FALSE 1999         ARID1B        <NA>          lo
    FALSE 2000           ARL1        <NA>          hi
    FALSE 2001         ARL2BP        <NA>          hi
    FALSE 2002          ARL5B        <NA>          hi
    FALSE 2003          ARL8A        <NA>          hi
    FALSE 2004         ARMC10        <NA>          hi
    FALSE 2005           ARNT        <NA>          lo
    FALSE 2006         ARPC1B        <NA>          lo
    FALSE 2007           ARSE        <NA>          lo
    FALSE 2008          ASB10        <NA>          lo
    FALSE 2009          ASB11        <NA>          lo
    FALSE 2010          ASB14        <NA>          hi
    FALSE 2011           ASB9        <NA>          hi
    FALSE 2012          ASH1L        <NA>          lo
    FALSE 2013          ASNA1        <NA>          lo
    FALSE 2014           ASPH        <NA>          hi
    FALSE 2015           ASPM        <NA>          hi
    FALSE 2016           ASS1        <NA>          lo
    FALSE 2017          ASTN1        <NA>          lo
    FALSE 2018          ASTN2        <NA>          lo
    FALSE 2019          ATCAY        <NA>          lo
    FALSE 2020          ATG2B        <NA>          hi
    FALSE 2021           ATG7        <NA>          hi
    FALSE 2022          ATHL1        <NA>          lo
    FALSE 2023        ATP13A3        <NA>          hi
    FALSE 2024         ATP1A3        <NA>          lo
    FALSE 2025         ATP1B1        <NA>          lo
    FALSE 2026         ATP2A1        <NA>          lo
    FALSE 2027          ATP4B        <NA>          hi
    FALSE 2028         ATP5C1        <NA>          hi
    FALSE 2029         ATP5G3        <NA>          hi
    FALSE 2030          ATP5I        <NA>          hi
    FALSE 2031          ATP5J        <NA>          hi
    FALSE 2032        ATP6AP2        <NA>          hi
    FALSE 2033       ATP6V0D1        <NA>          hi
    FALSE 2034       ATP6V0E1        <NA>          hi
    FALSE 2035       ATP6V0E2        <NA>          lo
    FALSE 2036        ATP6V1A        <NA>          hi
    FALSE 2037       ATP6V1E1        <NA>          hi
    FALSE 2038         ATP8A1        <NA>          hi
    FALSE 2039         ATP8B2        <NA>          lo
    FALSE 2040         ATP8B3        <NA>          lo
    FALSE 2041            ATR        <NA>          hi
    FALSE 2042           ATRX        <NA>          hi
    FALSE 2043          AURKA        <NA>          hi
    FALSE 2044          AXIN2        <NA>          lo
    FALSE 2045          AZIN1        <NA>          hi
    FALSE 2046       B3GALNT2        <NA>          hi
    FALSE 2047        B3GALT1        <NA>          lo
    FALSE 2048        B3GALT2        <NA>          hi
    FALSE 2049        B3GALT6        <NA>          lo
    FALSE 2050        B4GALT5        <NA>          lo
    FALSE 2051           BAG2        <NA>          lo
    FALSE 2052         BAHCC1        <NA>          lo
    FALSE 2053         BAIAP2        <NA>          lo
    FALSE 2054          BARD1        <NA>          hi
    FALSE 2055          BASP1        <NA>          hi
    FALSE 2056          BAZ1A        <NA>          hi
    FALSE 2057           BBS9        <NA>          hi
    FALSE 2058            BBX        <NA>          hi
    FALSE 2059         BCKDHA        <NA>          lo
    FALSE 2060         BCKDHB        <NA>          hi
    FALSE 2061         BCL11B        <NA>          lo
    FALSE 2062           BCL2        <NA>          hi
    FALSE 2063         BCORL1        <NA>          lo
    FALSE 2064          BECN1        <NA>          hi
    FALSE 2065           BET1        <NA>          hi
    FALSE 2066            BF1        <NA>          lo
    FALSE 2067        BHLHE41        <NA>          lo
    FALSE 2068          BICC1        <NA>          hi
    FALSE 2069          BICD1        <NA>          hi
    FALSE 2070          BIRC5        <NA>          hi
    FALSE 2071           BLB1        <NA>          lo
    FALSE 2072           BLMH        <NA>          hi
    FALSE 2073          BLVRA        <NA>          lo
    FALSE 2074          BLZF1        <NA>          hi
    FALSE 2075           BMI1        <NA>          lo
    FALSE 2076          BMP2K        <NA>          lo
    FALSE 2077         BMPR1B        <NA>          lo
    FALSE 2078          BNIP1        <NA>          hi
    FALSE 2079         BNIP3L        <NA>          hi
    FALSE 2080         BOD1L1        <NA>          hi
    FALSE 2081          BOLA3        <NA>          hi
    FALSE 2082           BORA        <NA>          hi
    FALSE 2083         BPIFB2        <NA>          lo
    FALSE 2084          BRCA1        <NA>          hi
    FALSE 2085         BRI3BP        <NA>          lo
    FALSE 2086          BRIP1        <NA>          hi
    FALSE 2087          BRPF3        <NA>          lo
    FALSE 2088          BRWD3        <NA>          hi
    FALSE 2089          BTBD1        <NA>          hi
    FALSE 2090         BTBD10        <NA>          hi
    FALSE 2091         BTBD17        <NA>          lo
    FALSE 2092          BTBD2        <NA>          lo
    FALSE 2093          BTBD3        <NA>          lo
    FALSE 2094          BTBD9        <NA>          hi
    FALSE 2095         BTF3L4        <NA>          hi
    FALSE 2096           BTG4        <NA>          hi
    FALSE 2097            BTK        <NA>          lo
    FALSE 2098           BUB1        <NA>          hi
    FALSE 2099          BUB1B        <NA>          hi
    FALSE 2100           BUB3        <NA>          hi
    FALSE 2101         BZRAP1        <NA>          lo
    FALSE 2102           BZW1        <NA>          hi
    FALSE 2103    C10H15ORF40        <NA>          hi
    FALSE 2104       C11ORF49        <NA>          lo
    FALSE 2105       C12ORF29        <NA>          hi
    FALSE 2106    C14H16ORF52        <NA>          hi
    FALSE 2107    C14H16ORF88        <NA>          hi
    FALSE 2108        C14ORF2        <NA>          hi
    FALSE 2109        C14ORF4        <NA>          lo
    FALSE 2110       C16ORF45        <NA>          hi
    FALSE 2111       C18ORF21        <NA>          hi
    FALSE 2112       C18ORF32        <NA>          hi
    FALSE 2113     C1H12ORF23        <NA>          hi
    FALSE 2114      C1H12ORF4        <NA>          hi
    FALSE 2115      C1H12ORF5        <NA>          hi
    FALSE 2116     C1H12ORF66        <NA>          hi
    FALSE 2117     C1H21ORF33        <NA>          lo
    FALSE 2118           C1QA        <NA>          lo
    FALSE 2119           C1QB        <NA>          lo
    FALSE 2120           C1QC        <NA>          lo
    FALSE 2121        C1QTNF1        <NA>          lo
    FALSE 2122        C1QTNF2        <NA>          lo
    FALSE 2123   C20H20ORF112        <NA>          lo
    FALSE 2124    C21H1ORF174        <NA>          lo
    FALSE 2125      C2H6ORF62        <NA>          hi
    FALSE 2126      C2H7ORF63        <NA>          hi
    FALSE 2127      C2H9ORF30        <NA>          hi
    FALSE 2128        C2ORF47        <NA>          hi
    FALSE 2129        C2ORF50        <NA>          lo
    FALSE 2130             C3        <NA>          lo
    FALSE 2131      C3H2ORF43        <NA>          lo
    FALSE 2132     C3H6ORF154        <NA>          lo
    FALSE 2133        C3ORF58        <NA>          hi
    FALSE 2134        C3ORF70        <NA>          lo
    FALSE 2135          C4BPA        <NA>          hi
    FALSE 2136      C4H4ORF21        <NA>          hi
    FALSE 2137      C4HXORF57        <NA>          hi
    FALSE 2138    C5H14ORF166        <NA>          hi
    FALSE 2139     C6H10ORF71        <NA>          hi
    FALSE 2140        C7ORF50        <NA>          hi
    FALSE 2141            C8G        <NA>          lo
    FALSE 2142     C8H1ORF173        <NA>          lo
    FALSE 2143      C8H1ORF21        <NA>          hi
    FALSE 2144      C8H1ORF27        <NA>          hi
    FALSE 2145      C9H3ORF55        <NA>          hi
    FALSE 2146            CA6        <NA>          lo
    FALSE 2147         CAB39L        <NA>          hi
    FALSE 2148        CABLES2        <NA>          hi
    FALSE 2149          CABP1        <NA>          hi
    FALSE 2150        CACNA1G        <NA>          lo
    FALSE 2151        CACNA1H        <NA>          lo
    FALSE 2152       CACNA2D1        <NA>          hi
    FALSE 2153       CALCOCO2        <NA>          lo
    FALSE 2154          CALII        <NA>          lo
    FALSE 2155          CALN1        <NA>          hi
    FALSE 2156        CAMK2N1        <NA>          lo
    FALSE 2157         CAMKK1        <NA>          lo
    FALSE 2158         CAMKK2        <NA>          lo
    FALSE 2159         CAMKMT        <NA>          hi
    FALSE 2160          CAMKV        <NA>          hi
    FALSE 2161           CANX        <NA>          hi
    FALSE 2162           CAP2        <NA>          lo
    FALSE 2163          CAPN3        <NA>          lo
    FALSE 2164          CARKD        <NA>          lo
    FALSE 2165           CARS        <NA>          hi
    FALSE 2166          CASC3        <NA>          lo
    FALSE 2167          CASC5        <NA>          hi
    FALSE 2168           CASK        <NA>          hi
    FALSE 2169        CASKIN1        <NA>          lo
    FALSE 2170          CASP1        <NA>          lo
    FALSE 2171       CASP8AP2        <NA>          hi
    FALSE 2172        CBFA2T2        <NA>          lo
    FALSE 2173          CBLL1        <NA>          lo
    FALSE 2174            CBS        <NA>          hi
    FALSE 2175           CBX4        <NA>          lo
    FALSE 2176          CCAR1        <NA>          hi
    FALSE 2177         CCDC15        <NA>          lo
    FALSE 2178        CCDC167        <NA>          hi
    FALSE 2179        CCDC177        <NA>          lo
    FALSE 2180          CCDC3        <NA>          lo
    FALSE 2181         CCDC58        <NA>          hi
    FALSE 2182        CCDC71L        <NA>          lo
    FALSE 2183        CCDC90B        <NA>          hi
    FALSE 2184         CCDC91        <NA>          hi
    FALSE 2185         CCDC97        <NA>          lo
    FALSE 2186          CCNA1        <NA>          lo
    FALSE 2187          CCNA2        <NA>          hi
    FALSE 2188          CCNB3        <NA>          hi
    FALSE 2189           CCNC        <NA>          hi
    FALSE 2190          CCNE2        <NA>          hi
    FALSE 2191          CCNT2        <NA>          hi
    FALSE 2192           CCR5        <NA>          lo
    FALSE 2193         CCSER2        <NA>          hi
    FALSE 2194           CCT2        <NA>          hi
    FALSE 2195           CCT8        <NA>          hi
    FALSE 2196          CD164        <NA>          hi
    FALSE 2197            CD2        <NA>          lo
    FALSE 2198          CD200        <NA>          lo
    FALSE 2199           CD24        <NA>          hi
    FALSE 2200         CD300A        <NA>          lo
    FALSE 2201           CD74        <NA>          lo
    FALSE 2202           CD86        <NA>          lo
    FALSE 2203          CDC20        <NA>          hi
    FALSE 2204         CDC25A        <NA>          hi
    FALSE 2205          CDC45        <NA>          hi
    FALSE 2206           CDC6        <NA>          hi
    FALSE 2207          CDCA3        <NA>          hi
    FALSE 2208          CDCA7        <NA>          hi
    FALSE 2209           CDH1        <NA>          lo
    FALSE 2210          CDH18        <NA>          hi
    FALSE 2211           CDH3        <NA>          lo
    FALSE 2212          CDIP1        <NA>          lo
    FALSE 2213           CDK1        <NA>          hi
    FALSE 2214          CDK18        <NA>          lo
    FALSE 2215          CDK19        <NA>          lo
    FALSE 2216           CDK5        <NA>          lo
    FALSE 2217       CDK5RAP2        <NA>          hi
    FALSE 2218       CDK5RAP3        <NA>          hi
    FALSE 2219           CDK6        <NA>          hi
    FALSE 2220           CDK8        <NA>          hi
    FALSE 2221          CDKN3        <NA>          hi
    FALSE 2222           CDR2        <NA>          hi
    FALSE 2223           CDS2        <NA>          lo
    FALSE 2224           CDT1        <NA>          hi
    FALSE 2225          CECR1        <NA>          lo
    FALSE 2226          CEMIP        <NA>          lo
    FALSE 2227          CENPA        <NA>          hi
    FALSE 2228          CENPC        <NA>          hi
    FALSE 2229          CENPE        <NA>          hi
    FALSE 2230          CENPF        <NA>          hi
    FALSE 2231          CENPI        <NA>          hi
    FALSE 2232          CENPK        <NA>          hi
    FALSE 2233          CENPM        <NA>          hi
    FALSE 2234          CENPN        <NA>          hi
    FALSE 2235          CENPO        <NA>          hi
    FALSE 2236         CEP104        <NA>          lo
    FALSE 2237         CEP128        <NA>          hi
    FALSE 2238         CEP131        <NA>          lo
    FALSE 2239         CEP152        <NA>          hi
    FALSE 2240         CEP170        <NA>          hi
    FALSE 2241         CEP192        <NA>          hi
    FALSE 2242         CEP290        <NA>          hi
    FALSE 2243          CEP55        <NA>          hi
    FALSE 2244         CEP85L        <NA>          lo
    FALSE 2245          CEP95        <NA>          hi
    FALSE 2246        CERS2L2        <NA>          lo
    FALSE 2247          CERS6        <NA>          hi
    FALSE 2248          CETN1        <NA>          hi
    FALSE 2249         CGRRF1        <NA>          hi
    FALSE 2250         CHCHD3        <NA>          hi
    FALSE 2251         CHCHD4        <NA>          hi
    FALSE 2252         CHCHD7        <NA>          hi
    FALSE 2253          CHD1L        <NA>          hi
    FALSE 2254           CHD4        <NA>          lo
    FALSE 2255           CHD5        <NA>          lo
    FALSE 2256          CHERP        <NA>          lo
    FALSE 2257           CHL1        <NA>          lo
    FALSE 2258           CHN1        <NA>          lo
    FALSE 2259           CHPF        <NA>          hi
    FALSE 2260         CHRNA5        <NA>          hi
    FALSE 2261          CHST1        <NA>          lo
    FALSE 2262         CHST10        <NA>          hi
    FALSE 2263          CHST8        <NA>          lo
    FALSE 2264          CHSY1        <NA>          hi
    FALSE 2265          CIART        <NA>          lo
    FALSE 2266            CIT        <NA>          hi
    FALSE 2267         CITED2        <NA>          lo
    FALSE 2268          CKAP2        <NA>          hi
    FALSE 2269          CKAP5        <NA>          hi
    FALSE 2270          CKS1B        <NA>          hi
    FALSE 2271         CLASP2        <NA>          hi
    FALSE 2272          CLDN1        <NA>          lo
    FALSE 2273         CLEC3A        <NA>          hi
    FALSE 2274           CLGN        <NA>          hi
    FALSE 2275          CLIC2        <NA>          lo
    FALSE 2276          CLIP4        <NA>          hi
    FALSE 2277         CLNS1A        <NA>          hi
    FALSE 2278        CLPTM1L        <NA>          hi
    FALSE 2279          CLSPN        <NA>          hi
    FALSE 2280         CLTCL1        <NA>          hi
    FALSE 2281          CLVS2        <NA>          hi
    FALSE 2282           CMIP        <NA>          lo
    FALSE 2283          CMPK2        <NA>          lo
    FALSE 2284          CMSS1        <NA>          hi
    FALSE 2285          CMTM3        <NA>          lo
    FALSE 2286          CMTM7        <NA>          lo
    FALSE 2287          CNIH1        <NA>          hi
    FALSE 2288         CNKSR1        <NA>          lo
    FALSE 2289           CNN2        <NA>          lo
    FALSE 2290          CNNM2        <NA>          lo
    FALSE 2291         CNPPD1        <NA>          hi
    FALSE 2292         CNRIP1        <NA>          hi
    FALSE 2293        COL19A1        <NA>          hi
    FALSE 2294        COL20A1        <NA>          hi
    FALSE 2295         COL4A3        <NA>          hi
    FALSE 2296           COPA        <NA>          hi
    FALSE 2297          COPS5        <NA>          hi
    FALSE 2298          COPS8        <NA>          hi
    FALSE 2299         COQ10B        <NA>          lo
    FALSE 2300           COQ2        <NA>          hi
    FALSE 2301           COQ5        <NA>          hi
    FALSE 2302          COX19        <NA>          lo
    FALSE 2303         COX4I1        <NA>          hi
    FALSE 2304         COX7A2        <NA>          hi
    FALSE 2305             CP        <NA>          lo
    FALSE 2306            CPE        <NA>          hi
    FALSE 2307          CPLX4        <NA>          hi
    FALSE 2308         CRABP1        <NA>          lo
    FALSE 2309          CREB1        <NA>          hi
    FALSE 2310        CREB3L1        <NA>          hi
    FALSE 2311        CREB3L2        <NA>          hi
    FALSE 2312         CRELD2        <NA>          hi
    FALSE 2313          CROCC        <NA>          lo
    FALSE 2314          CRTAP        <NA>          hi
    FALSE 2315           CRYM        <NA>          lo
    FALSE 2316          CSDE1        <NA>          hi
    FALSE 2317           CSF1        <NA>          lo
    FALSE 2318         CSF2RB        <NA>          lo
    FALSE 2319          CSF3R        <NA>          lo
    FALSE 2320     CSGALNACT2        <NA>          hi
    FALSE 2321          CSTF3        <NA>          hi
    FALSE 2322          CTBP1        <NA>          hi
    FALSE 2323  CTD-2410N18.5        <NA>          hi
    FALSE 2324  CTD-2561J22.3        <NA>          lo
    FALSE 2325         CTNNA3        <NA>          hi
    FALSE 2326       CTNNBIP1        <NA>          hi
    FALSE 2327          CTPS1        <NA>          hi
    FALSE 2328           CTR9        <NA>          hi
    FALSE 2329           CTSH        <NA>          lo
    FALSE 2330           CTSK        <NA>          lo
    FALSE 2331           CTSS        <NA>          lo
    FALSE 2332      CTTNBP2NL        <NA>          lo
    FALSE 2333           CUTC        <NA>          hi
    FALSE 2334           CUX1        <NA>          hi
    FALSE 2335           CUX2        <NA>          lo
    FALSE 2336          CWC22        <NA>          hi
    FALSE 2337         CX3CL1        <NA>          lo
    FALSE 2338         CX3CR1        <NA>          lo
    FALSE 2339          CXADR        <NA>          hi
    FALSE 2340         CXCL14        <NA>          lo
    FALSE 2341          CYB5A        <NA>          lo
    FALSE 2342           CYBB        <NA>          lo
    FALSE 2343        CYSLTR1        <NA>          lo
    FALSE 2344          DACH1        <NA>          hi
    FALSE 2345          DACH2        <NA>          hi
    FALSE 2346           DAG1        <NA>          lo
    FALSE 2347          DAPK3        <NA>          lo
    FALSE 2348         DAZAP2        <NA>          lo
    FALSE 2349          DBF4B        <NA>          lo
    FALSE 2350          DCAF8        <NA>          lo
    FALSE 2351            DCK        <NA>          hi
    FALSE 2352            DCX        <NA>          hi
    FALSE 2353           DCXR        <NA>          hi
    FALSE 2354         DDRGK1        <NA>          hi
    FALSE 2355           DDX1        <NA>          hi
    FALSE 2356          DDX24        <NA>          hi
    FALSE 2357          DDX27        <NA>          lo
    FALSE 2358          DDX46        <NA>          hi
    FALSE 2359          DDX59        <NA>          hi
    FALSE 2360        DENND1B        <NA>          hi
    FALSE 2361         DEPDC1        <NA>          hi
    FALSE 2362        DEPDC1B        <NA>          hi
    FALSE 2363         DEPDC6        <NA>          hi
    FALSE 2364           DERA        <NA>          hi
    FALSE 2365          DERL1        <NA>          hi
    FALSE 2366          DESI1        <NA>          hi
    FALSE 2367           DFFB        <NA>          lo
    FALSE 2368         DFNB31        <NA>          hi
    FALSE 2369           DGKH        <NA>          hi
    FALSE 2370           DGKI        <NA>          hi
    FALSE 2371          DGUOK        <NA>          hi
    FALSE 2372          DHX36        <NA>          hi
    FALSE 2373          DHX58        <NA>          lo
    FALSE 2374         DIAPH1        <NA>          hi
    FALSE 2375         DIAPH3        <NA>          hi
    FALSE 2376           DIO3        <NA>          lo
    FALSE 2377           DIS3        <NA>          hi
    FALSE 2378           DKK2        <NA>          hi
    FALSE 2379         DLGAP5        <NA>          hi
    FALSE 2380           DLK1        <NA>          lo
    FALSE 2381           DLL4        <NA>          lo
    FALSE 2382            DMA        <NA>          lo
    FALSE 2383         DMRTA2        <NA>          lo
    FALSE 2384           DMTN        <NA>          lo
    FALSE 2385           DNA2        <NA>          hi
    FALSE 2386         DNAJA3        <NA>          lo
    FALSE 2387        DNAJB11        <NA>          hi
    FALSE 2388        DNAJB12        <NA>          lo
    FALSE 2389         DNAJB9        <NA>          hi
    FALSE 2390        DNAJC12        <NA>          hi
    FALSE 2391         DNAJC2        <NA>          hi
    FALSE 2392          DOC2B        <NA>          lo
    FALSE 2393          DOCK2        <NA>          lo
    FALSE 2394          DOCK5        <NA>          hi
    FALSE 2395          DOCK9        <NA>          lo
    FALSE 2396           DOK5        <NA>          hi
    FALSE 2397           DPF1        <NA>          lo
    FALSE 2398           DPH5        <NA>          hi
    FALSE 2399           DPH7        <NA>          hi
    FALSE 2400           DPM1        <NA>          hi
    FALSE 2401          DPP10        <NA>          hi
    FALSE 2402          DPY30        <NA>          hi
    FALSE 2403           DPYD        <NA>          lo
    FALSE 2404         DPYSL3        <NA>          hi
    FALSE 2405            DR1        <NA>          hi
    FALSE 2406           DRG1        <NA>          hi
    FALSE 2407           DRP2        <NA>          lo
    FALSE 2408           DSC2        <NA>          hi
    FALSE 2409        DSCAML1        <NA>          lo
    FALSE 2410           DTD1        <NA>          hi
    FALSE 2411            DTL        <NA>          hi
    FALSE 2412         DTNBP1        <NA>          hi
    FALSE 2413           DUS2        <NA>          hi
    FALSE 2414         DUSP16        <NA>          lo
    FALSE 2415          DUSP4        <NA>          hi
    FALSE 2416          DUSP7        <NA>          lo
    FALSE 2417            DUT        <NA>          hi
    FALSE 2418       DYNC1LI1        <NA>          hi
    FALSE 2419       DYNC2LI1        <NA>          hi
    FALSE 2420         DYRK1A        <NA>          lo
    FALSE 2421           DYSF        <NA>          lo
    FALSE 2422          DZIP1        <NA>          hi
    FALSE 2423           E2F1        <NA>          hi
    FALSE 2424          EBAG9        <NA>          hi
    FALSE 2425          ECSCR        <NA>          lo
    FALSE 2426           ECT2        <NA>          hi
    FALSE 2427            EDA        <NA>          lo
    FALSE 2428          EDEM1        <NA>          hi
    FALSE 2429          EDEM2        <NA>          hi
    FALSE 2430          EDEM3        <NA>          hi
    FALSE 2431           EDN3        <NA>          lo
    FALSE 2432          EDNRA        <NA>          lo
    FALSE 2433         EEF1B2        <NA>          hi
    FALSE 2434         EFEMP1        <NA>          lo
    FALSE 2435          EFR3A        <NA>          hi
    FALSE 2436          EFR3B        <NA>          hi
    FALSE 2437         EFTUD2        <NA>          lo
    FALSE 2438         EGFLAM        <NA>          lo
    FALSE 2439         EHHADH        <NA>          lo
    FALSE 2440         EIF1AY        <NA>          hi
    FALSE 2441        EIF2AK1        <NA>          hi
    FALSE 2442        EIF2AK3        <NA>          hi
    FALSE 2443        EIF2AK4        <NA>          hi
    FALSE 2444         EIF2S1        <NA>          hi
    FALSE 2445         EIF2S2        <NA>          hi
    FALSE 2446        EIF2S3L        <NA>          hi
    FALSE 2447          EIF3B        <NA>          hi
    FALSE 2448          EIF3D        <NA>          hi
    FALSE 2449          EIF3E        <NA>          hi
    FALSE 2450          EIF3M        <NA>          hi
    FALSE 2451          EIF4H        <NA>          hi
    FALSE 2452           EIF5        <NA>          hi
    FALSE 2453          EIF5B        <NA>          hi
    FALSE 2454         ELAVL1        <NA>          lo
    FALSE 2455         ELAVL2        <NA>          hi
    FALSE 2456         ELAVL4        <NA>          hi
    FALSE 2457           ELF3        <NA>          lo
    FALSE 2458         ELMOD2        <NA>          hi
    FALSE 2459         ELOVL4        <NA>          hi
    FALSE 2460           EMC3        <NA>          hi
    FALSE 2461           EME1        <NA>          hi
    FALSE 2462           EML1        <NA>          lo
    FALSE 2463           EMX2        <NA>          lo
    FALSE 2464         ENDOD1        <NA>          hi
    FALSE 2465          ENOX1        <NA>          lo
    FALSE 2466          ENPP6        <NA>          lo
    FALSE 2467         ENTPD1        <NA>          lo
    FALSE 2468        ENTPD2L        <NA>          lo
    FALSE 2469         ENTPD5        <NA>          hi
    FALSE 2470         ENTPD7        <NA>          hi
    FALSE 2471          EPB41        <NA>          lo
    FALSE 2472        EPB41L5        <NA>          lo
    FALSE 2473          EPHX1        <NA>          lo
    FALSE 2474           EPRS        <NA>          hi
    FALSE 2475           EPS8        <NA>          hi
    FALSE 2476          ERCC6        <NA>          hi
    FALSE 2477         ERGIC1        <NA>          hi
    FALSE 2478            ERH        <NA>          hi
    FALSE 2479           ERI2        <NA>          hi
    FALSE 2480         ERLEC1        <NA>          hi
    FALSE 2481          ERO1L        <NA>          hi
    FALSE 2482         ERO1LB        <NA>          hi
    FALSE 2483          ERP29        <NA>          hi
    FALSE 2484          ERP44        <NA>          hi
    FALSE 2485          ESCO2        <NA>          hi
    FALSE 2486           ESF1        <NA>          hi
    FALSE 2487           ESPN        <NA>          lo
    FALSE 2488           ESR1        <NA>          lo
    FALSE 2489          ESYT1        <NA>          lo
    FALSE 2490           ETF1        <NA>          hi
    FALSE 2491           ETV1        <NA>          hi
    FALSE 2492           EXO1        <NA>          hi
    FALSE 2493        EXOC3L1        <NA>          lo
    FALSE 2494          EXOC4        <NA>          hi
    FALSE 2495         EXOSC7        <NA>          hi
    FALSE 2496          FADS6        <NA>          lo
    FALSE 2497          FAIM2        <NA>          lo
    FALSE 2498        FAM105A        <NA>          lo
    FALSE 2499        FAM107B        <NA>          lo
    FALSE 2500        FAM110B        <NA>          hi
    FALSE 2501         FAM13A        <NA>          lo
    FALSE 2502       FAM160A2        <NA>          lo
    FALSE 2503        FAM175B        <NA>          hi
    FALSE 2504       FAM177A1        <NA>          lo
    FALSE 2505        FAM180B        <NA>          lo
    FALSE 2506        FAM184B        <NA>          hi
    FALSE 2507       FAM189A2        <NA>          lo
    FALSE 2508        FAM18B1        <NA>          hi
    FALSE 2509        FAM193B        <NA>          lo
    FALSE 2510        FAM19A1        <NA>          hi
    FALSE 2511        FAM213A        <NA>          hi
    FALSE 2512        FAM219B        <NA>          lo
    FALSE 2513        FAM222B        <NA>          lo
    FALSE 2514         FAM46A        <NA>          hi
    FALSE 2515         FAM46D        <NA>          hi
    FALSE 2516         FAM65A        <NA>          lo
    FALSE 2517         FAM76A        <NA>          lo
    FALSE 2518         FAM78A        <NA>          lo
    FALSE 2519         FAM83G        <NA>          hi
    FALSE 2520         FAM84A        <NA>          hi
    FALSE 2521         FAM98A        <NA>          hi
    FALSE 2522         FANCD2        <NA>          hi
    FALSE 2523          FANCI        <NA>          hi
    FALSE 2524          FANCL        <NA>          hi
    FALSE 2525           FAR1        <NA>          hi
    FALSE 2526          FARSB        <NA>          hi
    FALSE 2527        FASTKD1        <NA>          hi
    FALSE 2528           FAT2        <NA>          lo
    FALSE 2529           FAT4        <NA>          lo
    FALSE 2530         FBRSL1        <NA>          lo
    FALSE 2531          FBXL7        <NA>          lo
    FALSE 2532         FBXO34        <NA>          lo
    FALSE 2533          FBXO5        <NA>          hi
    FALSE 2534          FBXW2        <NA>          hi
    FALSE 2535         FCER1G        <NA>          lo
    FALSE 2536          FCHO1        <NA>          lo
    FALSE 2537         FCHSD2        <NA>          lo
    FALSE 2538          FFAR4        <NA>          hi
    FALSE 2539           FGD4        <NA>          hi
    FALSE 2540          FGF12        <NA>          lo
    FALSE 2541           FGF7        <NA>          lo
    FALSE 2542           FGL2        <NA>          lo
    FALSE 2543          FHOD3        <NA>          hi
    FALSE 2544           FICD        <NA>          hi
    FALSE 2545         FKBP10        <NA>          lo
    FALSE 2546         FKBP11        <NA>          hi
    FALSE 2547         FKBP1B        <NA>          lo
    FALSE 2548          FKBP4        <NA>          hi
    FALSE 2549           FLI1        <NA>          lo
    FALSE 2550          FLOT2        <NA>          lo
    FALSE 2551            FN1        <NA>          lo
    FALSE 2552          FNBP4        <NA>          lo
    FALSE 2553         FNDC3A        <NA>          hi
    FALSE 2554           FNTB        <NA>          lo
    FALSE 2555          FOXC1        <NA>          lo
    FALSE 2556          FOXM1        <NA>          hi
    FALSE 2557          FOXO3        <NA>          lo
    FALSE 2558          FRAS1        <NA>          lo
    FALSE 2559         FRRS1L        <NA>          hi
    FALSE 2560          FSTL1        <NA>          lo
    FALSE 2561          FTSJ2        <NA>          hi
    FALSE 2562          FUT11        <NA>          lo
    FALSE 2563           FUT9        <NA>          hi
    FALSE 2564          FXYD2        <NA>          hi
    FALSE 2565        FYTTD1L        <NA>          lo
    FALSE 2566           FZD1        <NA>          lo
    FALSE 2567           FZD7        <NA>          lo
    FALSE 2568           FZD9        <NA>          lo
    FALSE 2569           G2E3        <NA>          hi
    FALSE 2570          G3BP2        <NA>          hi
    FALSE 2571           GAB1        <NA>          hi
    FALSE 2572      GABARAPL1        <NA>          hi
    FALSE 2573          GABRE        <NA>          hi
    FALSE 2574         GABRG1        <NA>          lo
    FALSE 2575          GABRP        <NA>          lo
    FALSE 2576        GAL3ST2        <NA>          hi
    FALSE 2577           GALE        <NA>          hi
    FALSE 2578         GALNT7        <NA>          hi
    FALSE 2579           GARS        <NA>          hi
    FALSE 2580         GAS2L3        <NA>          hi
    FALSE 2581          GATA2        <NA>          lo
    FALSE 2582         GATSL3        <NA>          lo
    FALSE 2583           GCLM        <NA>          hi
    FALSE 2584           GCM1        <NA>          lo
    FALSE 2585          GCOM1        <NA>          lo
    FALSE 2586           GCSH        <NA>          hi
    FALSE 2587          GDPD2        <NA>          lo
    FALSE 2588          GFRAL        <NA>          hi
    FALSE 2589           GGA2        <NA>          lo
    FALSE 2590           GGCT        <NA>          hi
    FALSE 2591            GGH        <NA>          hi
    FALSE 2592           GGT5        <NA>          lo
    FALSE 2593             GH        <NA>          lo
    FALSE 2594           GID8        <NA>          hi
    FALSE 2595            GIF        <NA>          lo
    FALSE 2596          GINS1        <NA>          hi
    FALSE 2597          GINS2        <NA>          hi
    FALSE 2598           GLB1        <NA>          hi
    FALSE 2599           GLG1        <NA>          hi
    FALSE 2600         GLIPR2        <NA>          hi
    FALSE 2601          GLIS2        <NA>          lo
    FALSE 2602           GLMN        <NA>          hi
    FALSE 2603          GLRX3        <NA>          hi
    FALSE 2604         GLT8D1        <NA>          hi
    FALSE 2605           GMNN        <NA>          hi
    FALSE 2606           GMPR        <NA>          lo
    FALSE 2607           GMPS        <NA>          hi
    FALSE 2608          GNAO1        <NA>          lo
    FALSE 2609           GNB1        <NA>          hi
    FALSE 2610           GNG4        <NA>          hi
    FALSE 2611           GNL2        <NA>          hi
    FALSE 2612           GNL3        <NA>          hi
    FALSE 2613          GNPTG        <NA>          hi
    FALSE 2614         GOLGA2        <NA>          hi
    FALSE 2615         GOLGA3        <NA>          hi
    FALSE 2616         GOLGA4        <NA>          hi
    FALSE 2617         GOLGA5        <NA>          hi
    FALSE 2618         GOLT1B        <NA>          hi
    FALSE 2619           GOPC        <NA>          hi
    FALSE 2620           GOT1        <NA>          hi
    FALSE 2621          GPD1L        <NA>          lo
    FALSE 2622           GPHN        <NA>          hi
    FALSE 2623            GPI        <NA>          hi
    FALSE 2624          GPNMB        <NA>          hi
    FALSE 2625         GPR123        <NA>          hi
    FALSE 2626         GPR157        <NA>          lo
    FALSE 2627         GPR162        <NA>          lo
    FALSE 2628         GPR180        <NA>          hi
    FALSE 2629          GPR25        <NA>          hi
    FALSE 2630          GPR34        <NA>          lo
    FALSE 2631          GPR75        <NA>          lo
    FALSE 2632         GPRIN3        <NA>          hi
    FALSE 2633           GPX4        <NA>          hi
    FALSE 2634           GPX7        <NA>          hi
    FALSE 2635          GRB14        <NA>          lo
    FALSE 2636           GRB2        <NA>          lo
    FALSE 2637          GRID1        <NA>          hi
    FALSE 2638          GRID2        <NA>          hi
    FALSE 2639        GRID2IP        <NA>          lo
    FALSE 2640          GRIK2        <NA>          hi
    FALSE 2641         GRIN2B        <NA>          lo
    FALSE 2642         GRIN2C        <NA>          lo
    FALSE 2643         GRXCR1        <NA>          hi
    FALSE 2644          GSG1L        <NA>          lo
    FALSE 2645          GSPT1        <NA>          hi
    FALSE 2646         GTF2H4        <NA>          lo
    FALSE 2647       GTF2IRD1        <NA>          lo
    FALSE 2648         GTF3C6        <NA>          hi
    FALSE 2649         GTPBP4        <NA>          hi
    FALSE 2650        GUCY1B4        <NA>          hi
    FALSE 2651         GXYLT2        <NA>          hi
    FALSE 2652           GYS2        <NA>          lo
    FALSE 2653          H2AFJ        <NA>          lo
    FALSE 2654           H6PD        <NA>          lo
    FALSE 2655         HAPLN3        <NA>          lo
    FALSE 2656           HAT1        <NA>          hi
    FALSE 2657          HAUS1        <NA>          hi
    FALSE 2658          HBS1L        <NA>          hi
    FALSE 2659            HCK        <NA>          lo
    FALSE 2660          HCLS1        <NA>          lo
    FALSE 2661          HDLBP        <NA>          hi
    FALSE 2662         HEATR3        <NA>          hi
    FALSE 2663          HELLS        <NA>          hi
    FALSE 2664        HERPUD1        <NA>          hi
    FALSE 2665           HES4        <NA>          lo
    FALSE 2666           HEY1        <NA>          lo
    FALSE 2667           HEYL        <NA>          lo
    FALSE 2668           HHEX        <NA>          lo
    FALSE 2669         HHIPL1        <NA>          lo
    FALSE 2670          HIF3A        <NA>          lo
    FALSE 2671          HIPK3        <NA>          lo
    FALSE 2672          HJURP        <NA>          hi
    FALSE 2673        HLA-DRA        <NA>          lo
    FALSE 2674            HLF        <NA>          lo
    FALSE 2675            HLX        <NA>          lo
    FALSE 2676           HM13        <NA>          hi
    FALSE 2677          HMGB2        <NA>          hi
    FALSE 2678        HMGCLL1        <NA>          lo
    FALSE 2679          HMGN1        <NA>          hi
    FALSE 2680          HMGN5        <NA>          hi
    FALSE 2681           HMMR        <NA>          hi
    FALSE 2682          HNF4A        <NA>          lo
    FALSE 2683      HNRNPA2B1        <NA>          hi
    FALSE 2684        HNRNPAB        <NA>          hi
    FALSE 2685        HNRNPH1        <NA>          hi
    FALSE 2686         HOMER2        <NA>          hi
    FALSE 2687         HOMER3        <NA>          lo
    FALSE 2688          HOOK1        <NA>          hi
    FALSE 2689         HPCAL1        <NA>          hi
    FALSE 2690          HPRT1        <NA>          hi
    FALSE 2691           HRH3        <NA>          hi
    FALSE 2692       HS3ST3A1        <NA>          lo
    FALSE 2693         HS6ST2        <NA>          hi
    FALSE 2694       HSD17B12        <NA>          hi
    FALSE 2695           HSF2        <NA>          lo
    FALSE 2696        HSP90B1        <NA>          hi
    FALSE 2697         HSPA13        <NA>          hi
    FALSE 2698         HSPA14        <NA>          hi
    FALSE 2699          HSPA5        <NA>          hi
    FALSE 2700          HSPG2        <NA>          lo
    FALSE 2701        HTATSF1        <NA>          hi
    FALSE 2702          HTR1B        <NA>          lo
    FALSE 2703          HYOU1        <NA>          hi
    FALSE 2704           IARS        <NA>          hi
    FALSE 2705          IARS2        <NA>          hi
    FALSE 2706           IBTK        <NA>          hi
    FALSE 2707           ICOS        <NA>          lo
    FALSE 2708            IDE        <NA>          hi
    FALSE 2709           IDH1        <NA>          lo
    FALSE 2710          IDH3A        <NA>          hi
    FALSE 2711           IDI1        <NA>          hi
    FALSE 2712          IFI35        <NA>          lo
    FALSE 2713          IFIH1        <NA>          lo
    FALSE 2714          IFIT5        <NA>          lo
    FALSE 2715         IFNGR1        <NA>          hi
    FALSE 2716         IFNGR2        <NA>          lo
    FALSE 2717          IFRD1        <NA>          hi
    FALSE 2718          IFT88        <NA>          hi
    FALSE 2719        IGF2BP2        <NA>          lo
    FALSE 2720          IGSF5        <NA>          hi
    FALSE 2721          IGSF9        <NA>          lo
    FALSE 2722          IKZF1        <NA>          lo
    FALSE 2723        IL12RB2        <NA>          lo
    FALSE 2724        IL13RA1        <NA>          lo
    FALSE 2725         IL18R1        <NA>          lo
    FALSE 2726       IL1RAPL1        <NA>          hi
    FALSE 2727          IL2RG        <NA>          lo
    FALSE 2728           IL34        <NA>          lo
    FALSE 2729            ILK        <NA>          lo
    FALSE 2730          IMPA1        <NA>          hi
    FALSE 2731         IMPACT        <NA>          hi
    FALSE 2732         IMPAD1        <NA>          hi
    FALSE 2733         INCENP        <NA>          hi
    FALSE 2734           INF2        <NA>          lo
    FALSE 2735         INPP5A        <NA>          hi
    FALSE 2736         INPP5D        <NA>          lo
    FALSE 2737          INTS3        <NA>          lo
    FALSE 2738         IPCEF1        <NA>          lo
    FALSE 2739           IPO5        <NA>          hi
    FALSE 2740           IPO7        <NA>          hi
    FALSE 2741           IRF1        <NA>          lo
    FALSE 2742           IRF6        <NA>          lo
    FALSE 2743           IRF7        <NA>          lo
    FALSE 2744           IRF8        <NA>          lo
    FALSE 2745           ISPD        <NA>          hi
    FALSE 2746          ITFG1        <NA>          hi
    FALSE 2747          ITGA3        <NA>          lo
    FALSE 2748          ITGA9        <NA>          lo
    FALSE 2749       ITGB1BP3        <NA>          hi
    FALSE 2750          ITGB2        <NA>          lo
    FALSE 2751          ITM2B        <NA>          hi
    FALSE 2752           JAG2        <NA>          lo
    FALSE 2753           JAK1        <NA>          lo
    FALSE 2754           JAK3        <NA>          lo
    FALSE 2755           KARS        <NA>          hi
    FALSE 2756           KAT5        <NA>          lo
    FALSE 2757          KAT6B        <NA>          lo
    FALSE 2758         KATNA1        <NA>          hi
    FALSE 2759         KATNB1        <NA>          lo
    FALSE 2760            KBP        <NA>          lo
    FALSE 2761         KBTBD3        <NA>          hi
    FALSE 2762          KCNB2        <NA>          hi
    FALSE 2763          KCNH2        <NA>          lo
    FALSE 2764          KCNH8        <NA>          lo
    FALSE 2765         KCNJ16        <NA>          lo
    FALSE 2766          KCNJ5        <NA>          hi
    FALSE 2767          KCNK1        <NA>          hi
    FALSE 2768         KCNK13        <NA>          lo
    FALSE 2769         KCNMA1        <NA>          hi
    FALSE 2770          KCNQ3        <NA>          hi
    FALSE 2771         KCTD21        <NA>          hi
    FALSE 2772         KDELC2        <NA>          hi
    FALSE 2773        KHDRBS3        <NA>          lo
    FALSE 2774       KIAA0040        <NA>          lo
    FALSE 2775       KIAA0195        <NA>          lo
    FALSE 2776       KIAA0556        <NA>          hi
    FALSE 2777       KIAA1107        <NA>          hi
    FALSE 2778       KIAA1524        <NA>          hi
    FALSE 2779      KIAA1549L        <NA>          lo
    FALSE 2780       KIAA1715        <NA>          hi
    FALSE 2781       KIAA1919        <NA>          lo
    FALSE 2782          KIF11        <NA>          hi
    FALSE 2783          KIF15        <NA>          hi
    FALSE 2784         KIF16B        <NA>          hi
    FALSE 2785         KIF18A        <NA>          hi
    FALSE 2786         KIF20A        <NA>          hi
    FALSE 2787          KIF4A        <NA>          hi
    FALSE 2788           KIF7        <NA>          lo
    FALSE 2789        KIRREL3        <NA>          lo
    FALSE 2790           KLF8        <NA>          lo
    FALSE 2791         KLHDC3        <NA>          hi
    FALSE 2792        KLHDC8B        <NA>          lo
    FALSE 2793          KLHL1        <NA>          hi
    FALSE 2794         KLHL11        <NA>          lo
    FALSE 2795         KLHL12        <NA>          lo
    FALSE 2796         KLHL13        <NA>          hi
    FALSE 2797         KLHL21        <NA>          lo
    FALSE 2798          KMT2E        <NA>          lo
    FALSE 2799          KNDC1        <NA>          hi
    FALSE 2800         KNSTRN        <NA>          hi
    FALSE 2801          KNTC1        <NA>          hi
    FALSE 2802          KPNA2        <NA>          hi
    FALSE 2803           KYNU        <NA>          hi
    FALSE 2804        L3MBTL4        <NA>          lo
    FALSE 2805          LAMP2        <NA>          hi
    FALSE 2806          LAMP3        <NA>          lo
    FALSE 2807         LANCL3        <NA>          hi
    FALSE 2808        LAPTM4B        <NA>          hi
    FALSE 2809         LAPTM5        <NA>          lo
    FALSE 2810           LAT2        <NA>          lo
    FALSE 2811            LBH        <NA>          hi
    FALSE 2812           LCAT        <NA>          lo
    FALSE 2813            LCK        <NA>          lo
    FALSE 2814         LCLAT1        <NA>          hi
    FALSE 2815           LCP1        <NA>          lo
    FALSE 2816           LCP2        <NA>          lo
    FALSE 2817           LDHA        <NA>          hi
    FALSE 2818       LEPROTL1        <NA>          hi
    FALSE 2819          LETM1        <NA>          hi
    FALSE 2820           LGMN        <NA>          hi
    FALSE 2821           LGR4        <NA>          lo
    FALSE 2822           LHFP        <NA>          lo
    FALSE 2823         LIMCH1        <NA>          lo
    FALSE 2824          LIN52        <NA>          hi
    FALSE 2825           LIPA        <NA>          lo
    FALSE 2826          LITAF        <NA>          lo
    FALSE 2827          LLGL2        <NA>          lo
    FALSE 2828          LMBR1        <NA>          hi
    FALSE 2829           LMO4        <NA>          lo
    FALSE 2830   LOC100857358        <NA>          lo
    FALSE 2831   LOC100857589        <NA>          lo
    FALSE 2832   LOC100857714        <NA>          lo
    FALSE 2833   LOC100857732        <NA>          lo
    FALSE 2834   LOC100858320        <NA>          lo
    FALSE 2835   LOC100859104        <NA>          hi
    FALSE 2836   LOC100859273        <NA>          lo
    FALSE 2837   LOC100859347        <NA>          lo
    FALSE 2838   LOC100859351        <NA>          lo
    FALSE 2839   LOC100859442        <NA>          hi
    FALSE 2840   LOC100859454        <NA>          lo
    FALSE 2841   LOC100859557        <NA>          hi
    FALSE 2842   LOC100859599        <NA>          lo
    FALSE 2843   LOC100859605        <NA>          lo
    FALSE 2844   LOC100859696        <NA>          lo
    FALSE 2845   LOC101747823        <NA>          lo
    FALSE 2846   LOC101748079        <NA>          lo
    FALSE 2847   LOC101748190        <NA>          lo
    FALSE 2848   LOC101748559        <NA>          hi
    FALSE 2849   LOC101748755        <NA>          lo
    FALSE 2850   LOC101748867        <NA>          lo
    FALSE 2851   LOC101748987        <NA>          lo
    FALSE 2852   LOC101749175        <NA>          hi
    FALSE 2853   LOC101749223        <NA>          lo
    FALSE 2854   LOC101749287        <NA>          hi
    FALSE 2855   LOC101749352        <NA>          lo
    FALSE 2856   LOC101749531        <NA>          hi
    FALSE 2857   LOC101749540        <NA>          hi
    FALSE 2858   LOC101749643        <NA>          lo
    FALSE 2859   LOC101749756        <NA>          lo
    FALSE 2860   LOC101749834        <NA>          hi
    FALSE 2861   LOC101750288        <NA>          hi
    FALSE 2862   LOC101750513        <NA>          lo
    FALSE 2863   LOC101750767        <NA>          hi
    FALSE 2864   LOC101750874        <NA>          lo
    FALSE 2865   LOC101751263        <NA>          lo
    FALSE 2866   LOC101751319        <NA>          lo
    FALSE 2867   LOC101751638        <NA>          hi
    FALSE 2868   LOC101751749        <NA>          hi
    FALSE 2869   LOC101751823        <NA>          hi
    FALSE 2870   LOC101752155        <NA>          lo
    FALSE 2871   LOC107049011        <NA>          lo
    FALSE 2872   LOC107049088        <NA>          hi
    FALSE 2873   LOC107049099        <NA>          hi
    FALSE 2874   LOC107049328        <NA>          lo
    FALSE 2875   LOC107049412        <NA>          lo
    FALSE 2876   LOC107049579        <NA>          lo
    FALSE 2877   LOC107049683        <NA>          lo
    FALSE 2878   LOC107049717        <NA>          lo
    FALSE 2879   LOC107049926        <NA>          lo
    FALSE 2880   LOC107050012        <NA>          lo
    FALSE 2881   LOC107050083        <NA>          lo
    FALSE 2882   LOC107050106        <NA>          lo
    FALSE 2883   LOC107050152        <NA>          lo
    FALSE 2884   LOC107050335        <NA>          lo
    FALSE 2885   LOC107050461        <NA>          lo
    FALSE 2886   LOC107050524        <NA>          lo
    FALSE 2887   LOC107050548        <NA>          lo
    FALSE 2888   LOC107050572        <NA>          lo
    FALSE 2889   LOC107050638        <NA>          lo
    FALSE 2890   LOC107050675        <NA>          lo
    FALSE 2891   LOC107050718        <NA>          lo
    FALSE 2892   LOC107050779        <NA>          lo
    FALSE 2893   LOC107050785        <NA>          lo
    FALSE 2894   LOC107050828        <NA>          lo
    FALSE 2895   LOC107050860        <NA>          lo
    FALSE 2896   LOC107050869        <NA>          lo
    FALSE 2897   LOC107050880        <NA>          lo
    FALSE 2898   LOC107050951        <NA>          lo
    FALSE 2899   LOC107050953        <NA>          hi
    FALSE 2900   LOC107051105        <NA>          lo
    FALSE 2901   LOC107051142        <NA>          lo
    FALSE 2902   LOC107051182        <NA>          hi
    FALSE 2903   LOC107051396        <NA>          lo
    FALSE 2904   LOC107051408        <NA>          lo
    FALSE 2905   LOC107051440        <NA>          lo
    FALSE 2906   LOC107051454        <NA>          lo
    FALSE 2907   LOC107051507        <NA>          lo
    FALSE 2908   LOC107051533        <NA>          lo
    FALSE 2909   LOC107051544        <NA>          lo
    FALSE 2910   LOC107052169        <NA>          lo
    FALSE 2911   LOC107052456        <NA>          lo
    FALSE 2912   LOC107052644        <NA>          lo
    FALSE 2913   LOC107053186        <NA>          lo
    FALSE 2914   LOC107053190        <NA>          lo
    FALSE 2915   LOC107053497        <NA>          hi
    FALSE 2916   LOC107053690        <NA>          lo
    FALSE 2917   LOC107054051        <NA>          hi
    FALSE 2918   LOC107054175        <NA>          lo
    FALSE 2919   LOC107054305        <NA>          hi
    FALSE 2920   LOC107054440        <NA>          lo
    FALSE 2921   LOC107054855        <NA>          lo
    FALSE 2922   LOC107055405        <NA>          lo
    FALSE 2923   LOC107055650        <NA>          hi
    FALSE 2924   LOC107055715        <NA>          hi
    FALSE 2925   LOC107056124        <NA>          hi
    FALSE 2926   LOC107056239        <NA>          lo
    FALSE 2927   LOC107056399        <NA>          lo
    FALSE 2928   LOC107056441        <NA>          lo
    FALSE 2929   LOC107057152        <NA>          lo
    FALSE 2930   LOC107057175        <NA>          hi
    FALSE 2931      LOC395100        <NA>          lo
    FALSE 2932      LOC395647        <NA>          hi
    FALSE 2933      LOC407092        <NA>          hi
    FALSE 2934      LOC414835        <NA>          hi
    FALSE 2935      LOC415641        <NA>          hi
    FALSE 2936      LOC416263        <NA>          lo
    FALSE 2937      LOC416633        <NA>          lo
    FALSE 2938      LOC417013        <NA>          hi
    FALSE 2939      LOC417372        <NA>          hi
    FALSE 2940      LOC418554        <NA>          lo
    FALSE 2941      LOC418667        <NA>          lo
    FALSE 2942      LOC418811        <NA>          hi
    FALSE 2943      LOC419112        <NA>          lo
    FALSE 2944      LOC420419        <NA>          lo
    FALSE 2945      LOC420516        <NA>          lo
    FALSE 2946      LOC420860        <NA>          hi
    FALSE 2947      LOC421441        <NA>          lo
    FALSE 2948      LOC421740        <NA>          lo
    FALSE 2949      LOC421856        <NA>          hi
    FALSE 2950      LOC421935        <NA>          lo
    FALSE 2951      LOC422151        <NA>          hi
    FALSE 2952      LOC422224        <NA>          hi
    FALSE 2953      LOC422249        <NA>          hi
    FALSE 2954      LOC422301        <NA>          lo
    FALSE 2955      LOC422513        <NA>          lo
    FALSE 2956      LOC423793        <NA>          hi
    FALSE 2957      LOC423967        <NA>          hi
    FALSE 2958      LOC424109        <NA>          hi
    FALSE 2959      LOC424201        <NA>          hi
    FALSE 2960      LOC424401        <NA>          hi
    FALSE 2961      LOC424919        <NA>          lo
    FALSE 2962      LOC424998        <NA>          hi
    FALSE 2963      LOC425431        <NA>          lo
    FALSE 2964      LOC426106        <NA>          lo
    FALSE 2965      LOC427665        <NA>          hi
    FALSE 2966      LOC428250        <NA>          lo
    FALSE 2967      LOC429167        <NA>          lo
    FALSE 2968      LOC429518        <NA>          lo
    FALSE 2969      LOC429809        <NA>          lo
    FALSE 2970      LOC429955        <NA>          lo
    FALSE 2971      LOC431656        <NA>          lo
    FALSE 2972      LOC768374        <NA>          lo
    FALSE 2973      LOC769039        <NA>          lo
    FALSE 2974      LOC769174        <NA>          hi
    FALSE 2975      LOC770612        <NA>          lo
    FALSE 2976      LOC771545        <NA>          lo
    FALSE 2977      LOC771758        <NA>          hi
    FALSE 2978      LOC771811        <NA>          hi
    FALSE 2979      LOC776275        <NA>          lo
    FALSE 2980          LOXL2        <NA>          lo
    FALSE 2981          LOXL3        <NA>          lo
    FALSE 2982          LPAR5        <NA>          lo
    FALSE 2983          LPHN2        <NA>          lo
    FALSE 2984          LPHN3        <NA>          lo
    FALSE 2985            LPO        <NA>          lo
    FALSE 2986          LPPR5        <NA>          hi
    FALSE 2987           LPXN        <NA>          lo
    FALSE 2988           LRP3        <NA>          lo
    FALSE 2989           LRP5        <NA>          hi
    FALSE 2990        LRRC10B        <NA>          lo
    FALSE 2991         LRRC47        <NA>          lo
    FALSE 2992         LRRC59        <NA>          hi
    FALSE 2993          LRRC7        <NA>          hi
    FALSE 2994        LRRC75A        <NA>          lo
    FALSE 2995         LRRC8D        <NA>          lo
    FALSE 2996          LRRK1        <NA>          lo
    FALSE 2997         LRRTM3        <NA>          hi
    FALSE 2998         LRRTM4        <NA>          hi
    FALSE 2999          LRWD1        <NA>          lo
    FALSE 3000           LSM1        <NA>          hi
    FALSE 3001          LTA4H        <NA>          hi
    FALSE 3002           LTN1        <NA>          hi
    FALSE 3003           LY86        <NA>          lo
    FALSE 3004           LY96        <NA>          lo
    FALSE 3005         LYPD6B        <NA>          lo
    FALSE 3006         LYPLA1        <NA>          hi
    FALSE 3007        LYPLAL1        <NA>          hi
    FALSE 3008          LYRM4        <NA>          hi
    FALSE 3009         LYSMD2        <NA>          hi
    FALSE 3010         LYSMD4        <NA>          hi
    FALSE 3011        MAB21L3        <NA>          hi
    FALSE 3012         MAD2L1        <NA>          hi
    FALSE 3013       MAD2L1BP        <NA>          hi
    FALSE 3014           MADD        <NA>          lo
    FALSE 3015           MAEA        <NA>          hi
    FALSE 3016           MAFB        <NA>          lo
    FALSE 3017         MAMLD1        <NA>          lo
    FALSE 3018         MAN1A1        <NA>          hi
    FALSE 3019         MAN1B1        <NA>          lo
    FALSE 3020           MANF        <NA>          hi
    FALSE 3021          MAP10        <NA>          lo
    FALSE 3022         MAP2K1        <NA>          hi
    FALSE 3023         MAP2K2        <NA>          lo
    FALSE 3024         MAP2K3        <NA>          hi
    FALSE 3025         MAP3K7        <NA>          hi
    FALSE 3026           MAP4        <NA>          lo
    FALSE 3027         MAP6D1        <NA>          lo
    FALSE 3028         MAP7D2        <NA>          hi
    FALSE 3029      MAPK1IP1L        <NA>          lo
    FALSE 3030          MAPK8        <NA>          hi
    FALSE 3031          MARC2        <NA>          hi
    FALSE 3032         MARCKS        <NA>          lo
    FALSE 3033          MARCO        <NA>          lo
    FALSE 3034          MAST3        <NA>          lo
    FALSE 3035          MASTL        <NA>          hi
    FALSE 3036         MB21D2        <NA>          hi
    FALSE 3037           MBD2        <NA>          hi
    FALSE 3038           MCAM        <NA>          lo
    FALSE 3039          MCFD2        <NA>          hi
    FALSE 3040          MCM10        <NA>          hi
    FALSE 3041           MCM2        <NA>          hi
    FALSE 3042           MCM5        <NA>          hi
    FALSE 3043           MCM6        <NA>          hi
    FALSE 3044           MCM8        <NA>          hi
    FALSE 3045         MCOLN1        <NA>          lo
    FALSE 3046            MCU        <NA>          hi
    FALSE 3047           MDC1        <NA>          lo
    FALSE 3048           MDH1        <NA>          hi
    FALSE 3049           MDM4        <NA>          lo
    FALSE 3050          MED22        <NA>          hi
    FALSE 3051          MED28        <NA>          hi
    FALSE 3052          MED30        <NA>          hi
    FALSE 3053          MEIS1        <NA>          hi
    FALSE 3054          MERTK        <NA>          lo
    FALSE 3055            MET        <NA>          lo
    FALSE 3056         METAP1        <NA>          hi
    FALSE 3057        METAP1D        <NA>          hi
    FALSE 3058       METTL21A        <NA>          hi
    FALSE 3059        METTL22        <NA>          hi
    FALSE 3060         MFSD11        <NA>          hi
    FALSE 3061         MFSD2A        <NA>          hi
    FALSE 3062          MFSD5        <NA>          lo
    FALSE 3063           MIA3        <NA>          hi
    FALSE 3064           MIB1        <NA>          hi
    FALSE 3065        MICALL1        <NA>          lo
    FALSE 3066        MICALL2        <NA>          lo
    FALSE 3067        MID1IP1        <NA>          lo
    FALSE 3068         MIPOL1        <NA>          hi
    FALSE 3069          MIS12        <NA>          hi
    FALSE 3070          MITD1        <NA>          hi
    FALSE 3071           MITF        <NA>          lo
    FALSE 3072          MKI67        <NA>          hi
    FALSE 3073          MKNK2        <NA>          lo
    FALSE 3074           MLEC        <NA>          hi
    FALSE 3075            MLN        <NA>          hi
    FALSE 3076            MLX        <NA>          lo
    FALSE 3077          MMP16        <NA>          hi
    FALSE 3078         MMS22L        <NA>          hi
    FALSE 3079            MN1        <NA>          lo
    FALSE 3080           MND1        <NA>          hi
    FALSE 3081         MOSPD1        <NA>          hi
    FALSE 3082          MOV10        <NA>          lo
    FALSE 3083          MPEG1        <NA>          lo
    FALSE 3084      MPHOSPH10        <NA>          hi
    FALSE 3085          MRAP2        <NA>          lo
    FALSE 3086           MRC2        <NA>          lo
    FALSE 3087         MRPL15        <NA>          hi
    FALSE 3088         MRPL47        <NA>          hi
    FALSE 3089         MRPS22        <NA>          hi
    FALSE 3090          MRPS9        <NA>          hi
    FALSE 3091        MSANTD1        <NA>          hi
    FALSE 3092        MSANTD2        <NA>          lo
    FALSE 3093        MSANTD4        <NA>          hi
    FALSE 3094           MSH2        <NA>          hi
    FALSE 3095           MSL1        <NA>          lo
    FALSE 3096          MSMO1        <NA>          hi
    FALSE 3097          MSTO1        <NA>          lo
    FALSE 3098            MT4        <NA>          lo
    FALSE 3099           MTBP        <NA>          hi
    FALSE 3100           MTF1        <NA>          lo
    FALSE 3101          MTFR1        <NA>          hi
    FALSE 3102           MTG2        <NA>          hi
    FALSE 3103          MTIF2        <NA>          hi
    FALSE 3104          MTMR3        <NA>          lo
    FALSE 3105           MTPN        <NA>          hi
    FALSE 3106            MTR        <NA>          hi
    FALSE 3107           MTRR        <NA>          hi
    FALSE 3108            MX1        <NA>          lo
    FALSE 3109           MXD3        <NA>          hi
    FALSE 3110        MYBBP1A        <NA>          hi
    FALSE 3111          MYBL1        <NA>          hi
    FALSE 3112         MYBPC1        <NA>          lo
    FALSE 3113         MYCBP2        <NA>          hi
    FALSE 3114          MYD88        <NA>          lo
    FALSE 3115          MYEF2        <NA>          hi
    FALSE 3116         MYO15A        <NA>          lo
    FALSE 3117         MYO18A        <NA>          lo
    FALSE 3118          MYO1C        <NA>          lo
    FALSE 3119          MYO1D        <NA>          lo
    FALSE 3120          MYO9B        <NA>          lo
    FALSE 3121           MYOC        <NA>          lo
    FALSE 3122           MYOT        <NA>          hi
    FALSE 3123        N4BP2L1        <NA>          lo
    FALSE 3124          NAA25        <NA>          hi
    FALSE 3125        NAALAD2        <NA>          lo
    FALSE 3126          NABP1        <NA>          hi
    FALSE 3127          NAMPT        <NA>          lo
    FALSE 3128           NAPG        <NA>          hi
    FALSE 3129           NASP        <NA>          hi
    FALSE 3130          NAT10        <NA>          hi
    FALSE 3131          NAT8L        <NA>          lo
    FALSE 3132           NBAS        <NA>          hi
    FALSE 3133         NBEAL2        <NA>          lo
    FALSE 3134           NBR1        <NA>          hi
    FALSE 3135          NCAM2        <NA>          hi
    FALSE 3136         NCAPD3        <NA>          hi
    FALSE 3137          NCAPG        <NA>          hi
    FALSE 3138          NCAPH        <NA>          hi
    FALSE 3139           NCF1        <NA>          lo
    FALSE 3140           NCF2        <NA>          lo
    FALSE 3141         NCKAP1        <NA>          hi
    FALSE 3142            NCL        <NA>          hi
    FALSE 3143           NCLN        <NA>          hi
    FALSE 3144          NCOA3        <NA>          lo
    FALSE 3145          NDC80        <NA>          hi
    FALSE 3146          NDNL2        <NA>          hi
    FALSE 3147         NDUFB5        <NA>          hi
    FALSE 3148         NDUFB9        <NA>          hi
    FALSE 3149         NDUFC2        <NA>          hi
    FALSE 3150         NDUFS1        <NA>          hi
    FALSE 3151         NDUFV2        <NA>          hi
    FALSE 3152         NECAP1        <NA>          hi
    FALSE 3153           NEFH        <NA>          lo
    FALSE 3154           NEK2        <NA>          hi
    FALSE 3155           NEK6        <NA>          hi
    FALSE 3156           NEK9        <NA>          lo
    FALSE 3157          NELL2        <NA>          lo
    FALSE 3158           NEO1        <NA>          lo
    FALSE 3159           NET1        <NA>          lo
    FALSE 3160          NETO1        <NA>          hi
    FALSE 3161           NEU3        <NA>          lo
    FALSE 3162         NEURL1        <NA>          hi
    FALSE 3163          NFASC        <NA>          lo
    FALSE 3164         NFE2L1        <NA>          lo
    FALSE 3165           NFIC        <NA>          lo
    FALSE 3166           NFIX        <NA>          lo
    FALSE 3167          NFKB2        <NA>          lo
    FALSE 3168         NFKBIE        <NA>          lo
    FALSE 3169          NGLY1        <NA>          hi
    FALSE 3170         NHEDC2        <NA>          hi
    FALSE 3171          NHSL2        <NA>          lo
    FALSE 3172           NIP7        <NA>          hi
    FALSE 3173         NIPAL3        <NA>          hi
    FALSE 3174          NISCH        <NA>          hi
    FALSE 3175         NKAIN3        <NA>          hi
    FALSE 3176         NKX2-2        <NA>          lo
    FALSE 3177          NLRC5        <NA>          lo
    FALSE 3178         NMRAL1        <NA>          lo
    FALSE 3179           NMT2        <NA>          hi
    FALSE 3180           NNF1        <NA>          hi
    FALSE 3181            NNT        <NA>          hi
    FALSE 3182          NOL11        <NA>          hi
    FALSE 3183           NOL4        <NA>          lo
    FALSE 3184           NOL9        <NA>          lo
    FALSE 3185           NOM1        <NA>          hi
    FALSE 3186          NOP14        <NA>          hi
    FALSE 3187          NOP16        <NA>          hi
    FALSE 3188          NOP58        <NA>          hi
    FALSE 3189           NOS2        <NA>          lo
    FALSE 3190         NOTCH2        <NA>          lo
    FALSE 3191          NPHP1        <NA>          lo
    FALSE 3192          NPHP4        <NA>          lo
    FALSE 3193           NPM1        <NA>          hi
    FALSE 3194           NPM3        <NA>          hi
    FALSE 3195           NPR2        <NA>          lo
    FALSE 3196          NPRL2        <NA>          lo
    FALSE 3197           NPTN        <NA>          lo
    FALSE 3198          NR0B1        <NA>          lo
    FALSE 3199          NR1D1        <NA>          lo
    FALSE 3200          NR1H3        <NA>          lo
    FALSE 3201          NR5A1        <NA>          lo
    FALSE 3202           NRAS        <NA>          lo
    FALSE 3203          NRBP2        <NA>          lo
    FALSE 3204           NRF1        <NA>          lo
    FALSE 3205           NRG3        <NA>          lo
    FALSE 3206           NRGN        <NA>          lo
    FALSE 3207          NRIP3        <NA>          lo
    FALSE 3208            NRK        <NA>          lo
    FALSE 3209           NRP2        <NA>          lo
    FALSE 3210          NSUN2        <NA>          hi
    FALSE 3211         NT5C1B        <NA>          hi
    FALSE 3212          NT5C2        <NA>          hi
    FALSE 3213         NT5C3A        <NA>          hi
    FALSE 3214         NT5DC3        <NA>          lo
    FALSE 3215           NTN1        <NA>          hi
    FALSE 3216          NTN4L        <NA>          hi
    FALSE 3217          NTNG1        <NA>          lo
    FALSE 3218          NUCB2        <NA>          hi
    FALSE 3219         NUDCD2        <NA>          hi
    FALSE 3220          NUDT6        <NA>          hi
    FALSE 3221           NUF2        <NA>          hi
    FALSE 3222           NUMB        <NA>          lo
    FALSE 3223          NUP58        <NA>          hi
    FALSE 3224           NUS1        <NA>          hi
    FALSE 3225            NVL        <NA>          hi
    FALSE 3226            NXN        <NA>          lo
    FALSE 3227           OASL        <NA>          lo
    FALSE 3228           OAZ1        <NA>          hi
    FALSE 3229            OC3        <NA>          lo
    FALSE 3230         OCIAD1        <NA>          hi
    FALSE 3231          OFCC1        <NA>          hi
    FALSE 3232           OGFR        <NA>          lo
    FALSE 3233            OGN        <NA>          lo
    FALSE 3234           OIT3        <NA>          hi
    FALSE 3235           OLA1        <NA>          hi
    FALSE 3236          OLFM1        <NA>          hi
    FALSE 3237          OLFM2        <NA>          hi
    FALSE 3238        OLFML2A        <NA>          lo
    FALSE 3239            OMP        <NA>          hi
    FALSE 3240           OPN3        <NA>          lo
    FALSE 3241          OPRL1        <NA>          lo
    FALSE 3242           ORC1        <NA>          hi
    FALSE 3243        OSBPL10        <NA>          lo
    FALSE 3244        OSBPL11        <NA>          hi
    FALSE 3245          OSER1        <NA>          hi
    FALSE 3246           OSTC        <NA>          hi
    FALSE 3247          OSTM1        <NA>          lo
    FALSE 3248          OTUD1        <NA>          lo
    FALSE 3249         OTUD6B        <NA>          hi
    FALSE 3250         OTUD7A        <NA>          hi
    FALSE 3251           OVST        <NA>          hi
    FALSE 3252          OVSTL        <NA>          hi
    FALSE 3253           OXR1        <NA>          hi
    FALSE 3254          OXSR1        <NA>          lo
    FALSE 3255          P2RX5        <NA>          lo
    FALSE 3256           P4HB        <NA>          hi
    FALSE 3257        PAK1IP1        <NA>          hi
    FALSE 3258          PALD1        <NA>          lo
    FALSE 3259            PAM        <NA>          hi
    FALSE 3260         PAPSS1        <NA>          hi
    FALSE 3261         PAPSS2        <NA>          lo
    FALSE 3262          PARK7        <NA>          hi
    FALSE 3263          PARVA        <NA>          hi
    FALSE 3264          PARVG        <NA>          lo
    FALSE 3265           PASK        <NA>          hi
    FALSE 3266            PBK        <NA>          hi
    FALSE 3267         PCASP2        <NA>          lo
    FALSE 3268          PCBP2        <NA>          lo
    FALSE 3269          PCBP3        <NA>          lo
    FALSE 3270         PCDH18        <NA>          hi
    FALSE 3271          PCDH7        <NA>          hi
    FALSE 3272          PCDH9        <NA>          hi
    FALSE 3273        PCDHA12        <NA>          hi
    FALSE 3274         PCDHA3        <NA>          hi
    FALSE 3275        PCDHGC3        <NA>          lo
    FALSE 3276          PCGF6        <NA>          hi
    FALSE 3277           PCLO        <NA>          hi
    FALSE 3278           PCM1        <NA>          hi
    FALSE 3279           PCNA        <NA>          hi
    FALSE 3280           PCP4        <NA>          hi
    FALSE 3281          PCSK1        <NA>          hi
    FALSE 3282          PCSK6        <NA>          hi
    FALSE 3283          PCSK7        <NA>          lo
    FALSE 3284         PCYT1B        <NA>          hi
    FALSE 3285          PDCD4        <NA>          hi
    FALSE 3286         PDE10A        <NA>          hi
    FALSE 3287         PDE11A        <NA>          hi
    FALSE 3288          PDE12        <NA>          lo
    FALSE 3289          PDE1C        <NA>          hi
    FALSE 3290          PDE5A        <NA>          lo
    FALSE 3291          PDGFB        <NA>          lo
    FALSE 3292          PDGFC        <NA>          lo
    FALSE 3293         PDGFRB        <NA>          lo
    FALSE 3294         PDGFRL        <NA>          hi
    FALSE 3295          PDIA3        <NA>          hi
    FALSE 3296          PDIA4        <NA>          hi
    FALSE 3297         PDLIM7        <NA>          lo
    FALSE 3298         PDXDC1        <NA>          hi
    FALSE 3299          PEAR1        <NA>          lo
    FALSE 3300           PEMT        <NA>          lo
    FALSE 3301           PER3        <NA>          lo
    FALSE 3302          PEX10        <NA>          hi
    FALSE 3303           PEX2        <NA>          hi
    FALSE 3304           PEX6        <NA>          lo
    FALSE 3305          PGAM1        <NA>          hi
    FALSE 3306          PGAP1        <NA>          hi
    FALSE 3307           PGM2        <NA>          hi
    FALSE 3308           PHC1        <NA>          lo
    FALSE 3309          PHF14        <NA>          hi
    FALSE 3310           PHF2        <NA>          lo
    FALSE 3311          PHGDH        <NA>          hi
    FALSE 3312           PHIP        <NA>          hi
    FALSE 3313         PHLDB1        <NA>          lo
    FALSE 3314         PI4K2B        <NA>          hi
    FALSE 3315          PICK1        <NA>          hi
    FALSE 3316           PIGA        <NA>          hi
    FALSE 3317           PIGB        <NA>          hi
    FALSE 3318           PIGK        <NA>          hi
    FALSE 3319           PIGM        <NA>          hi
    FALSE 3320           PIGR        <NA>          lo
    FALSE 3321        PIK3AP1        <NA>          lo
    FALSE 3322        PIK3C2B        <NA>          lo
    FALSE 3323        PIK3C2G        <NA>          lo
    FALSE 3324         PIK3CG        <NA>          hi
    FALSE 3325         PIK3R6        <NA>          lo
    FALSE 3326          PINK1        <NA>          lo
    FALSE 3327         PINLYP        <NA>          lo
    FALSE 3328        PIP4K2A        <NA>          lo
    FALSE 3329         PITRM1        <NA>          hi
    FALSE 3330          PITX2        <NA>          lo
    FALSE 3331          PITX3        <NA>          lo
    FALSE 3332         PIWIL1        <NA>          lo
    FALSE 3333           PKD2        <NA>          lo
    FALSE 3334         PKDCCb        <NA>          lo
    FALSE 3335         PKDREJ        <NA>          lo
    FALSE 3336        PLA2G4E        <NA>          lo
    FALSE 3337         PLA2G7        <NA>          hi
    FALSE 3338           PLAU        <NA>          lo
    FALSE 3339          PLBD2        <NA>          lo
    FALSE 3340          PLCH1        <NA>          lo
    FALSE 3341           PLD4        <NA>          lo
    FALSE 3342        PLEKHA6        <NA>          lo
    FALSE 3343        PLEKHG1        <NA>          hi
    FALSE 3344        PLEKHG3        <NA>          lo
    FALSE 3345        PLEKHM1        <NA>          lo
    FALSE 3346           PLK1        <NA>          hi
    FALSE 3347           PLP1        <NA>          hi
    FALSE 3348           PLS1        <NA>          hi
    FALSE 3349           PLS3        <NA>          hi
    FALSE 3350         PLXND1        <NA>          lo
    FALSE 3351            PM5        <NA>          hi
    FALSE 3352            PML        <NA>          lo
    FALSE 3353           PMM2        <NA>          hi
    FALSE 3354           PMS2        <NA>          hi
    FALSE 3355         PNPLA8        <NA>          hi
    FALSE 3356          PNRC2        <NA>          hi
    FALSE 3357         POFUT2        <NA>          hi
    FALSE 3358           POGZ        <NA>          lo
    FALSE 3359         POLR3H        <NA>          hi
    FALSE 3360         POM121        <NA>          lo
    FALSE 3361        POMGNT1        <NA>          lo
    FALSE 3362         POU3F2        <NA>          hi
    FALSE 3363         POU3F3        <NA>          hi
    FALSE 3364           PPA1        <NA>          hi
    FALSE 3365         PPAP2B        <NA>          lo
    FALSE 3366           PPIB        <NA>          hi
    FALSE 3367       PPP1R16B        <NA>          lo
    FALSE 3368        PPP1R3B        <NA>          lo
    FALSE 3369          PPP6C        <NA>          lo
    FALSE 3370         PPP6R3        <NA>          hi
    FALSE 3371           PRC1        <NA>          hi
    FALSE 3372          PRDX3        <NA>          hi
    FALSE 3373          PRDX4        <NA>          hi
    FALSE 3374        PRELID1        <NA>          hi
    FALSE 3375       PRICKLE2        <NA>          lo
    FALSE 3376         PRIMA1        <NA>          lo
    FALSE 3377        PRIMPOL        <NA>          hi
    FALSE 3378         PRKAB2        <NA>          lo
    FALSE 3379         PRKAG3        <NA>          lo
    FALSE 3380        PRKAR2A        <NA>          hi
    FALSE 3381          PRKCB        <NA>          lo
    FALSE 3382          PRKCD        <NA>          lo
    FALSE 3383          PRKDC        <NA>          hi
    FALSE 3384            PRL        <NA>          hi
    FALSE 3385          PRMT3        <NA>          hi
    FALSE 3386          PROS1        <NA>          hi
    FALSE 3387          PRPS2        <NA>          hi
    FALSE 3388          PRR15        <NA>          lo
    FALSE 3389           PRR7        <NA>          lo
    FALSE 3390         PRSS12        <NA>          hi
    FALSE 3391         PRSS23        <NA>          hi
    FALSE 3392            PSD        <NA>          lo
    FALSE 3393           PSD2        <NA>          lo
    FALSE 3394          PSMA3        <NA>          hi
    FALSE 3395          PSMB3        <NA>          hi
    FALSE 3396          PSMC1        <NA>          hi
    FALSE 3397         PSMD12        <NA>          hi
    FALSE 3398          PSMD6        <NA>          hi
    FALSE 3399          PSMD7        <NA>          hi
    FALSE 3400          PSMF1        <NA>          lo
    FALSE 3401          PSMG2        <NA>          hi
    FALSE 3402           PSPH        <NA>          hi
    FALSE 3403          PTBP1        <NA>          lo
    FALSE 3404         PTCHD2        <NA>          lo
    FALSE 3405           PTEN        <NA>          lo
    FALSE 3406         PTGER3        <NA>          lo
    FALSE 3407          PTGFR        <NA>          hi
    FALSE 3408          PTGS1        <NA>          lo
    FALSE 3409         PTPN14        <NA>          lo
    FALSE 3410          PTPRC        <NA>          lo
    FALSE 3411          PTPRO        <NA>          lo
    FALSE 3412          PTPRS        <NA>          lo
    FALSE 3413         PTPRVP        <NA>          lo
    FALSE 3414          PTRH2        <NA>          hi
    FALSE 3415          PUS10        <NA>          hi
    FALSE 3416          PXMP4        <NA>          lo
    FALSE 3417            PXN        <NA>          lo
    FALSE 3418          PYCR2        <NA>          hi
    FALSE 3419           PYGL        <NA>          lo
    FALSE 3420           QDPR        <NA>          hi
    FALSE 3421      RAB11FIP5        <NA>          lo
    FALSE 3422          RAB1A        <NA>          hi
    FALSE 3423          RAB26        <NA>          hi
    FALSE 3424          RAB28        <NA>          hi
    FALSE 3425          RAB2A        <NA>          hi
    FALSE 3426          RAB30        <NA>          hi
    FALSE 3427          RAB32        <NA>          lo
    FALSE 3428          RAB3B        <NA>          lo
    FALSE 3429       RAB3GAP2        <NA>          hi
    FALSE 3430        RAB3IL1        <NA>          lo
    FALSE 3431          RAB5C        <NA>          lo
    FALSE 3432       RABGAP1L        <NA>          hi
    FALSE 3433          RABL3        <NA>          hi
    FALSE 3434        RACGAP1        <NA>          hi
    FALSE 3435          RAD21        <NA>          hi
    FALSE 3436         RAD54B        <NA>          hi
    FALSE 3437          RADIL        <NA>          hi
    FALSE 3438         RALGDS        <NA>          lo
    FALSE 3439         RANBP1        <NA>          hi
    FALSE 3440        RANBP17        <NA>          lo
    FALSE 3441         RANBP2        <NA>          hi
    FALSE 3442        RANGAP1        <NA>          hi
    FALSE 3443          RAP1A        <NA>          hi
    FALSE 3444       RAP1GAP2        <NA>          lo
    FALSE 3445        RAPGEF5        <NA>          lo
    FALSE 3446           RARB        <NA>          lo
    FALSE 3447        RARRES1        <NA>          lo
    FALSE 3448          RASA2        <NA>          hi
    FALSE 3449         RASSF1        <NA>          hi
    FALSE 3450         RASSF3        <NA>          hi
    FALSE 3451         RASSF6        <NA>          hi
    FALSE 3452            RB1        <NA>          hi
    FALSE 3453         RB1CC1        <NA>          hi
    FALSE 3454          RBBP8        <NA>          hi
    FALSE 3455          RBM20        <NA>          lo
    FALSE 3456          RBM38        <NA>          lo
    FALSE 3457          RBM41        <NA>          hi
    FALSE 3458          RBM47        <NA>          lo
    FALSE 3459          RBMS2        <NA>          lo
    FALSE 3460         RBPMS2        <NA>          lo
    FALSE 3461          RCAN3        <NA>          lo
    FALSE 3462           RCN1        <NA>          hi
    FALSE 3463           RDM1        <NA>          lo
    FALSE 3464          REEP5        <NA>          hi
    FALSE 3465            REL        <NA>          lo
    FALSE 3466           RER1        <NA>          hi
    FALSE 3467           REST        <NA>          hi
    FALSE 3468            RET        <NA>          hi
    FALSE 3469           REV1        <NA>          lo
    FALSE 3470           RFC3        <NA>          hi
    FALSE 3471           RFC5        <NA>          hi
    FALSE 3472           RGCC        <NA>          hi
    FALSE 3473           RGL1        <NA>          lo
    FALSE 3474          RGS11        <NA>          lo
    FALSE 3475          RGS18        <NA>          lo
    FALSE 3476           RGS3        <NA>          lo
    FALSE 3477           RGS4        <NA>          lo
    FALSE 3478           RGS7        <NA>          hi
    FALSE 3479         RHBDD2        <NA>          lo
    FALSE 3480         RHBDF1        <NA>          lo
    FALSE 3481         RHBDL1        <NA>          lo
    FALSE 3482           RHOA        <NA>          hi
    FALSE 3483          RHOT2        <NA>          lo
    FALSE 3484          RIBC2        <NA>          lo
    FALSE 3485          RINT1        <NA>          hi
    FALSE 3486          RIOK1        <NA>          hi
    FALSE 3487         RNF145        <NA>          lo
    FALSE 3488         RNF149        <NA>          lo
    FALSE 3489         RNF182        <NA>          hi
    FALSE 3490         RNF213        <NA>          lo
    FALSE 3491          RNF34        <NA>          lo
    FALSE 3492           RNH1        <NA>          hi
    FALSE 3493           ROM1        <NA>          lo
    FALSE 3494 RP11-1035H13.3        <NA>          hi
    FALSE 3495 RP11-152F13.10        <NA>          lo
    FALSE 3496 RP11-155D18.12        <NA>          lo
    FALSE 3497   RP11-156P1.2        <NA>          hi
    FALSE 3498   RP11-166N6.3        <NA>          lo
    FALSE 3499  RP11-196G11.1        <NA>          hi
    FALSE 3500  RP11-407N17.3        <NA>          hi
    FALSE 3501           RPF2        <NA>          hi
    FALSE 3502           RPGR        <NA>          hi
    FALSE 3503          RPH3A        <NA>          lo
    FALSE 3504          RPL35        <NA>          hi
    FALSE 3505           RPL7        <NA>          hi
    FALSE 3506           RPN1        <NA>          hi
    FALSE 3507          RPP25        <NA>          lo
    FALSE 3508          RPP30        <NA>          hi
    FALSE 3509          RPP40        <NA>          hi
    FALSE 3510         RPRD1B        <NA>          lo
    FALSE 3511          RPS21        <NA>          hi
    FALSE 3512        RPS6KA6        <NA>          hi
    FALSE 3513          RRAGB        <NA>          lo
    FALSE 3514          RRAS2        <NA>          lo
    FALSE 3515           RRM2        <NA>          hi
    FALSE 3516          RRP15        <NA>          hi
    FALSE 3517          RRP1B        <NA>          hi
    FALSE 3518          RSAD2        <NA>          lo
    FALSE 3519         RSL1D1        <NA>          hi
    FALSE 3520           RSU1        <NA>          lo
    FALSE 3521           RTCA        <NA>          hi
    FALSE 3522           RTN1        <NA>          lo
    FALSE 3523           RTN4        <NA>          hi
    FALSE 3524        RTN4IP1        <NA>          hi
    FALSE 3525          RTN4R        <NA>          lo
    FALSE 3526           RTTN        <NA>          hi
    FALSE 3527          RUFY1        <NA>          hi
    FALSE 3528        RUNX1T1        <NA>          hi
    FALSE 3529          RUSC1        <NA>          lo
    FALSE 3530        S100A10        <NA>          lo
    FALSE 3531          S1PR2        <NA>          lo
    FALSE 3532         SAMD11        <NA>          lo
    FALSE 3533         SAMM50        <NA>          hi
    FALSE 3534         SAMSN1        <NA>          lo
    FALSE 3535          SARNP        <NA>          lo
    FALSE 3536          SASH3        <NA>          lo
    FALSE 3537           SAV1        <NA>          lo
    FALSE 3538           SBK1        <NA>          lo
    FALSE 3539         SCAMP3        <NA>          lo
    FALSE 3540           SCAP        <NA>          lo
    FALSE 3541         SCARA3        <NA>          lo
    FALSE 3542         SCARB1        <NA>          lo
    FALSE 3543          SCFD1        <NA>          hi
    FALSE 3544           SCG3        <NA>          hi
    FALSE 3545         SCHIP1        <NA>          lo
    FALSE 3546           SCIN        <NA>          lo
    FALSE 3547          SCN2B        <NA>          lo
    FALSE 3548          SCN5A        <NA>          hi
    FALSE 3549          SCRN3        <NA>          hi
    FALSE 3550          SCYL2        <NA>          hi
    FALSE 3551         SDF2L1        <NA>          hi
    FALSE 3552           SDF4        <NA>          hi
    FALSE 3553         SEC11C        <NA>          hi
    FALSE 3554          SEC13        <NA>          hi
    FALSE 3555         SEC22C        <NA>          hi
    FALSE 3556         SEC23B        <NA>          hi
    FALSE 3557        SEC23IP        <NA>          hi
    FALSE 3558         SEC24D        <NA>          hi
    FALSE 3559        SEC61A2        <NA>          hi
    FALSE 3560         SEC61B        <NA>          hi
    FALSE 3561         SEC61G        <NA>          hi
    FALSE 3562           SELK        <NA>          hi
    FALSE 3563         SELPLG        <NA>          lo
    FALSE 3564         SEMA3B        <NA>          lo
    FALSE 3565         SEMA3F        <NA>          lo
    FALSE 3566         SEMA3G        <NA>          lo
    FALSE 3567         SEMA4B        <NA>          lo
    FALSE 3568         SEMA4C        <NA>          lo
    FALSE 3569         SEMA4G        <NA>          lo
    FALSE 3570         SEMA6B        <NA>          lo
    FALSE 3571          SENP1        <NA>          lo
    FALSE 3572          SEP15        <NA>          hi
    FALSE 3573          SEPP1        <NA>          hi
    FALSE 3574         SEPP1L        <NA>          hi
    FALSE 3575        SERINC1        <NA>          hi
    FALSE 3576        SERINC2        <NA>          lo
    FALSE 3577          SERP1        <NA>          hi
    FALSE 3578      SERPINB10        <NA>          lo
    FALSE 3579       SERPINB6        <NA>          hi
    FALSE 3580       SERPING1        <NA>          lo
    FALSE 3581          SETD3        <NA>          hi
    FALSE 3582          SETD8        <NA>          hi
    FALSE 3583         SETDB1        <NA>          lo
    FALSE 3584           SEZ6        <NA>          lo
    FALSE 3585          SF3B4        <NA>          lo
    FALSE 3586          SF3B6        <NA>          lo
    FALSE 3587           SFPQ        <NA>          hi
    FALSE 3588          SFRP5        <NA>          hi
    FALSE 3589          SFXN4        <NA>          hi
    FALSE 3590           SGCB        <NA>          lo
    FALSE 3591           SGCZ        <NA>          hi
    FALSE 3592         SGK223        <NA>          lo
    FALSE 3593          SGMS1        <NA>          hi
    FALSE 3594          SGOL1        <NA>          hi
    FALSE 3595          SGSM3        <NA>          hi
    FALSE 3596         SH3BP1        <NA>          lo
    FALSE 3597         SH3BP2        <NA>          lo
    FALSE 3598         SH3GL3        <NA>          lo
    FALSE 3599       SH3PXD2B        <NA>          lo
    FALSE 3600         SH3TC1        <NA>          hi
    FALSE 3601           SHC2        <NA>          hi
    FALSE 3602         SHCBP1        <NA>          hi
    FALSE 3603          SHOX2        <NA>          hi
    FALSE 3604          SHPRH        <NA>          hi
    FALSE 3605        SHROOM4        <NA>          lo
    FALSE 3606          SIRT5        <NA>          hi
    FALSE 3607           SIX3        <NA>          hi
    FALSE 3608           SKA1        <NA>          hi
    FALSE 3609           SKA3        <NA>          hi
    FALSE 3610         SKIDA1        <NA>          lo
    FALSE 3611           SKIL        <NA>          hi
    FALSE 3612            SLA        <NA>          lo
    FALSE 3613         SLAMF1        <NA>          lo
    FALSE 3614           SLBP        <NA>          hi
    FALSE 3615        SLC10A4        <NA>          lo
    FALSE 3616        SLC12A9        <NA>          lo
    FALSE 3617        SLC16A1        <NA>          hi
    FALSE 3618        SLC16A2        <NA>          lo
    FALSE 3619        SLC16A9        <NA>          lo
    FALSE 3620        SLC17A5        <NA>          hi
    FALSE 3621        SLC17A9        <NA>          hi
    FALSE 3622         SLC1A2        <NA>          hi
    FALSE 3623         SLC1A4        <NA>          hi
    FALSE 3624       SLC22A15        <NA>          hi
    FALSE 3625        SLC23A2        <NA>          hi
    FALSE 3626        SLC25A1        <NA>          lo
    FALSE 3627       SLC25A12        <NA>          hi
    FALSE 3628       SLC25A21        <NA>          hi
    FALSE 3629       SLC25A22        <NA>          hi
    FALSE 3630        SLC25A3        <NA>          hi
    FALSE 3631       SLC25A34        <NA>          hi
    FALSE 3632       SLC25A37        <NA>          lo
    FALSE 3633       SLC25A44        <NA>          lo
    FALSE 3634         SLC2A1        <NA>          lo
    FALSE 3635         SLC2A5        <NA>          lo
    FALSE 3636        SLC30A7        <NA>          hi
    FALSE 3637        SLC33A1        <NA>          hi
    FALSE 3638        SLC35B1        <NA>          hi
    FALSE 3639        SLC35B3        <NA>          hi
    FALSE 3640        SLC35D1        <NA>          lo
    FALSE 3641        SLC35F2        <NA>          hi
    FALSE 3642        SLC36A4        <NA>          hi
    FALSE 3643        SLC37A3        <NA>          hi
    FALSE 3644        SLC37A4        <NA>          lo
    FALSE 3645        SLC38A2        <NA>          hi
    FALSE 3646        SLC38A3        <NA>          lo
    FALSE 3647       SLC39A13        <NA>          hi
    FALSE 3648        SLC40A1        <NA>          lo
    FALSE 3649        SLC41A2        <NA>          hi
    FALSE 3650        SLC44A5        <NA>          lo
    FALSE 3651         SLC4A2        <NA>          lo
    FALSE 3652        SLC5A10        <NA>          hi
    FALSE 3653         SLC6A1        <NA>          lo
    FALSE 3654        SLC6A14        <NA>          lo
    FALSE 3655        SLC6A17        <NA>          lo
    FALSE 3656         SLC7A3        <NA>          hi
    FALSE 3657         SLC7A5        <NA>          hi
    FALSE 3658         SLC7A6        <NA>          hi
    FALSE 3659         SLC9A1        <NA>          lo
    FALSE 3660         SLC9A2        <NA>          lo
    FALSE 3661       SLC9A3R1        <NA>          lo
    FALSE 3662       SLC9A3R2        <NA>          lo
    FALSE 3663         SLC9A6        <NA>          hi
    FALSE 3664        SLCO4A1        <NA>          hi
    FALSE 3665        SLITRK6        <NA>          hi
    FALSE 3666            SLK        <NA>          hi
    FALSE 3667          SLMO2        <NA>          hi
    FALSE 3668         SLX4IP        <NA>          hi
    FALSE 3669       SMARCAL1        <NA>          hi
    FALSE 3670        SMARCD1        <NA>          lo
    FALSE 3671        SMARCD2        <NA>          lo
    FALSE 3672          SMCO4        <NA>          hi
    FALSE 3673           SMG8        <NA>          hi
    FALSE 3674         SMIM18        <NA>          hi
    FALSE 3675          SMIM5        <NA>          lo
    FALSE 3676         SMNDC1        <NA>          hi
    FALSE 3677            SMO        <NA>          lo
    FALSE 3678          SMPD3        <NA>          hi
    FALSE 3679           SMTN        <NA>          lo
    FALSE 3680         SMTNL1        <NA>          hi
    FALSE 3681         SNRPA1        <NA>          hi
    FALSE 3682          SNRPC        <NA>          lo
    FALSE 3683          SNRPF        <NA>          hi
    FALSE 3684          SNTB1        <NA>          lo
    FALSE 3685          SNX14        <NA>          hi
    FALSE 3686          SNX20        <NA>          lo
    FALSE 3687          SNX27        <NA>          lo
    FALSE 3688          SNX33        <NA>          lo
    FALSE 3689           SNX6        <NA>          hi
    FALSE 3690           SOD1        <NA>          hi
    FALSE 3691          SOGA3        <NA>          hi
    FALSE 3692          SORT1        <NA>          lo
    FALSE 3693        SOSTDC1        <NA>          hi
    FALSE 3694           SOX2        <NA>          lo
    FALSE 3695           SOX4        <NA>          hi
    FALSE 3696            SP1        <NA>          lo
    FALSE 3697          SPAM1        <NA>          hi
    FALSE 3698         SPATA5        <NA>          hi
    FALSE 3699          SPC24        <NA>          hi
    FALSE 3700          SPC25        <NA>          hi
    FALSE 3701          SPCS1        <NA>          hi
    FALSE 3702          SPCS2        <NA>          hi
    FALSE 3703          SPCS3        <NA>          hi
    FALSE 3704          SPHK1        <NA>          lo
    FALSE 3705           SPI1        <NA>          lo
    FALSE 3706          SPNS3        <NA>          lo
    FALSE 3707         SPOCK3        <NA>          hi
    FALSE 3708          SPOPL        <NA>          hi
    FALSE 3709           SPP1        <NA>          lo
    FALSE 3710          SPPL3        <NA>          lo
    FALSE 3711         SPRED1        <NA>          hi
    FALSE 3712           SPTB        <NA>          lo
    FALSE 3713           SQLE        <NA>          hi
    FALSE 3714         SRCIN1        <NA>          lo
    FALSE 3715           SRGN        <NA>          lo
    FALSE 3716            SRM        <NA>          hi
    FALSE 3717          SRP14        <NA>          hi
    FALSE 3718          SRP54        <NA>          hi
    FALSE 3719          SRP72        <NA>          hi
    FALSE 3720           SRPR        <NA>          hi
    FALSE 3721          SRPRB        <NA>          hi
    FALSE 3722          SRSF3        <NA>          lo
    FALSE 3723         SS18L2        <NA>          hi
    FALSE 3724            SS2        <NA>          hi
    FALSE 3725            SSB        <NA>          hi
    FALSE 3726           SSPN        <NA>          lo
    FALSE 3727           SSPO        <NA>          hi
    FALSE 3728           SSR1        <NA>          hi
    FALSE 3729           SSR2        <NA>          hi
    FALSE 3730           SSR3        <NA>          hi
    FALSE 3731           SSR4        <NA>          hi
    FALSE 3732          SSTR4        <NA>          lo
    FALSE 3733          SSTR5        <NA>          lo
    FALSE 3734           ST18        <NA>          hi
    FALSE 3735     ST6GALNAC6        <NA>          lo
    FALSE 3736          STAB1        <NA>          lo
    FALSE 3737           STAC        <NA>          hi
    FALSE 3738       STAMBPL1        <NA>          hi
    FALSE 3739          STAT2        <NA>          lo
    FALSE 3740         STAT5B        <NA>          lo
    FALSE 3741           STC1        <NA>          hi
    FALSE 3742         STEAP1        <NA>          hi
    FALSE 3743         STEAP2        <NA>          hi
    FALSE 3744          STK10        <NA>          lo
    FALSE 3745          STK16        <NA>          hi
    FALSE 3746           STK3        <NA>          hi
    FALSE 3747           STK4        <NA>          hi
    FALSE 3748          STMN1        <NA>          hi
    FALSE 3749          STMN3        <NA>          lo
    FALSE 3750          STRAP        <NA>          hi
    FALSE 3751          STX12        <NA>          hi
    FALSE 3752          STX17        <NA>          lo
    FALSE 3753         STXBP5        <NA>          hi
    FALSE 3754           SUB1        <NA>          hi
    FALSE 3755         SUB1L1        <NA>          hi
    FALSE 3756           SUCO        <NA>          hi
    FALSE 3757          SUGP1        <NA>          lo
    FALSE 3758         SUPT5H        <NA>          lo
    FALSE 3759         SUPT6H        <NA>          lo
    FALSE 3760          SURF4        <NA>          hi
    FALSE 3761          SYDE1        <NA>          lo
    FALSE 3762         SYNGR1        <NA>          lo
    FALSE 3763          SYNPO        <NA>          lo
    FALSE 3764          SYNPR        <NA>          hi
    FALSE 3765           SYT1        <NA>          lo
    FALSE 3766          SYT13        <NA>          hi
    FALSE 3767          TACC1        <NA>          lo
    FALSE 3768          TACR3        <NA>          hi
    FALSE 3769         TALDO1        <NA>          hi
    FALSE 3770         TAPBPL        <NA>          lo
    FALSE 3771           TARP        <NA>          lo
    FALSE 3772        TAX1BP1        <NA>          hi
    FALSE 3773        TBC1D14        <NA>          lo
    FALSE 3774        TBC1D15        <NA>          hi
    FALSE 3775       TBC1D22B        <NA>          lo
    FALSE 3776        TBC1D32        <NA>          hi
    FALSE 3777         TBC1D9        <NA>          hi
    FALSE 3778           TBCE        <NA>          hi
    FALSE 3779          TBX20        <NA>          lo
    FALSE 3780           TBX3        <NA>          lo
    FALSE 3781           TBX5        <NA>          lo
    FALSE 3782         TBXAS1        <NA>          hi
    FALSE 3783        TCERG1L        <NA>          hi
    FALSE 3784           TCP1        <NA>          hi
    FALSE 3785          TEAD1        <NA>          hi
    FALSE 3786            TEK        <NA>          lo
    FALSE 3787          TEKT1        <NA>          lo
    FALSE 3788          TENM2        <NA>          hi
    FALSE 3789           TESC        <NA>          hi
    FALSE 3790          TESK2        <NA>          lo
    FALSE 3791             TF        <NA>          lo
    FALSE 3792         TFAP2B        <NA>          hi
    FALSE 3793          TFCP2        <NA>          lo
    FALSE 3794          TFDP1        <NA>          hi
    FALSE 3795           TFE3        <NA>          lo
    FALSE 3796          TFPI2        <NA>          hi
    FALSE 3797           TGDS        <NA>          hi
    FALSE 3798          TGFB1        <NA>          lo
    FALSE 3799          TGFB3        <NA>          lo
    FALSE 3800         TGFBR2        <NA>          lo
    FALSE 3801           TGM2        <NA>          lo
    FALSE 3802           TGS1        <NA>          hi
    FALSE 3803          THADA        <NA>          hi
    FALSE 3804          THAP5        <NA>          hi
    FALSE 3805          THOC1        <NA>          hi
    FALSE 3806          THOC2        <NA>          hi
    FALSE 3807         THSD7B        <NA>          lo
    FALSE 3808           THY1        <NA>          lo
    FALSE 3809         TIMM44        <NA>          hi
    FALSE 3810           TJP3        <NA>          lo
    FALSE 3811            TK1        <NA>          hi
    FALSE 3812           TLL1        <NA>          lo
    FALSE 3813          TLR2A        <NA>          lo
    FALSE 3814           TLR7        <NA>          lo
    FALSE 3815          TM2D3        <NA>          hi
    FALSE 3816         TM9SF3        <NA>          hi
    FALSE 3817           TMC5        <NA>          hi
    FALSE 3818         TMED10        <NA>          hi
    FALSE 3819         TMEFF2        <NA>          lo
    FALSE 3820        TMEM104        <NA>          lo
    FALSE 3821       TMEM120A        <NA>          lo
    FALSE 3822       TMEM120B        <NA>          hi
    FALSE 3823        TMEM135        <NA>          lo
    FALSE 3824        TMEM154        <NA>          lo
    FALSE 3825        TMEM159        <NA>          hi
    FALSE 3826        TMEM181        <NA>          hi
    FALSE 3827       TMEM185A        <NA>          hi
    FALSE 3828        TMEM186        <NA>          hi
    FALSE 3829        TMEM196        <NA>          hi
    FALSE 3830        TMEM240        <NA>          hi
    FALSE 3831        TMEM241        <NA>          hi
    FALSE 3832        TMEM243        <NA>          lo
    FALSE 3833        TMEM39A        <NA>          hi
    FALSE 3834        TMEM39B        <NA>          lo
    FALSE 3835         TMEM40        <NA>          lo
    FALSE 3836        TMEM50A        <NA>          hi
    FALSE 3837        TMEM52B        <NA>          hi
    FALSE 3838         TMEM57        <NA>          hi
    FALSE 3839         TMEM59        <NA>          hi
    FALSE 3840           TMF1        <NA>          hi
    FALSE 3841          TMOD3        <NA>          hi
    FALSE 3842           TMPO        <NA>          hi
    FALSE 3843       TMPRSS13        <NA>          hi
    FALSE 3844        TMPRSS9        <NA>          lo
    FALSE 3845         TMSB10        <NA>          hi
    FALSE 3846          TMTC3        <NA>          hi
    FALSE 3847           TMX3        <NA>          hi
    FALSE 3848             TN        <NA>          lo
    FALSE 3849        TNFAIP3        <NA>          lo
    FALSE 3850        TNFAIP6        <NA>          lo
    FALSE 3851      TNFAIP8L1        <NA>          lo
    FALSE 3852       TNFRSF25        <NA>          lo
    FALSE 3853        TNFSF11        <NA>          lo
    FALSE 3854          TNKS2        <NA>          hi
    FALSE 3855          TNNI1        <NA>          lo
    FALSE 3856           TOB2        <NA>          lo
    FALSE 3857         TOMM20        <NA>          hi
    FALSE 3858          TOP2A        <NA>          hi
    FALSE 3859        TP53I11        <NA>          lo
    FALSE 3860       TP53INP1        <NA>          lo
    FALSE 3861           TP73        <NA>          lo
    FALSE 3862           TPP2        <NA>          hi
    FALSE 3863         TPRG1L        <NA>          lo
    FALSE 3864          TPRKB        <NA>          hi
    FALSE 3865           TPT1        <NA>          hi
    FALSE 3866          TRA2B        <NA>          hi
    FALSE 3867          TRABD        <NA>          hi
    FALSE 3868          TRAF1        <NA>          lo
    FALSE 3869       TRAF3IP2        <NA>          lo
    FALSE 3870       TRAF3IP3        <NA>          lo
    FALSE 3871          TRAF4        <NA>          lo
    FALSE 3872          TRAF7        <NA>          lo
    FALSE 3873         TRAFD1        <NA>          lo
    FALSE 3874          TRAK2        <NA>          lo
    FALSE 3875          TRAM1        <NA>          hi
    FALSE 3876        TREM-B2        <NA>          lo
    FALSE 3877         TRIM24        <NA>          lo
    FALSE 3878         TRIM35        <NA>          hi
    FALSE 3879       TRIM39.1        <NA>          lo
    FALSE 3880         TRIOBP        <NA>          lo
    FALSE 3881         TRIP11        <NA>          hi
    FALSE 3882         TRMT11        <NA>          hi
    FALSE 3883         TRMT1L        <NA>          hi
    FALSE 3884          TRPC1        <NA>          hi
    FALSE 3885          TRPV2        <NA>          lo
    FALSE 3886          TRPV4        <NA>          lo
    FALSE 3887          TSHZ1        <NA>          hi
    FALSE 3888            TSN        <NA>          hi
    FALSE 3889        TSPAN13        <NA>          hi
    FALSE 3890         TSPAN8        <NA>          lo
    FALSE 3891         TSPAN9        <NA>          lo
    FALSE 3892           TSPO        <NA>          lo
    FALSE 3893           TSR3        <NA>          hi
    FALSE 3894          TSTA3        <NA>          hi
    FALSE 3895          TTC13        <NA>          hi
    FALSE 3896          TTC17        <NA>          hi
    FALSE 3897         TTC21B        <NA>          hi
    FALSE 3898          TTC27        <NA>          hi
    FALSE 3899           TTC3        <NA>          hi
    FALSE 3900           TTC9        <NA>          lo
    FALSE 3901           TTF1        <NA>          hi
    FALSE 3902            TTK        <NA>          hi
    FALSE 3903            TTL        <NA>          hi
    FALSE 3904         TTLL11        <NA>          hi
    FALSE 3905          TTLL4        <NA>          lo
    FALSE 3906            TUB        <NA>          hi
    FALSE 3907          TUBD1        <NA>          hi
    FALSE 3908          TUBG1        <NA>          hi
    FALSE 3909        TUBGCP5        <NA>          hi
    FALSE 3910          TUSC2        <NA>          lo
    FALSE 3911          TUSC3        <NA>          hi
    FALSE 3912        TXNDC11        <NA>          hi
    FALSE 3913         TXNDC5        <NA>          hi
    FALSE 3914          TXNL1        <NA>          hi
    FALSE 3915           TYMS        <NA>          hi
    FALSE 3916           TYW5        <NA>          hi
    FALSE 3917          U2AF1        <NA>          hi
    FALSE 3918         UAP1L1        <NA>          lo
    FALSE 3919           UBA5        <NA>          hi
    FALSE 3920        UBASH3B        <NA>          lo
    FALSE 3921          UBE2A        <NA>          hi
    FALSE 3922          UBE2C        <NA>          hi
    FALSE 3923         UBE2E1        <NA>          lo
    FALSE 3924         UBE2V1        <NA>          hi
    FALSE 3925         UBLCP1        <NA>          hi
    FALSE 3926           UBR3        <NA>          hi
    FALSE 3927           UBR7        <NA>          hi
    FALSE 3928         UBXN2A        <NA>          hi
    FALSE 3929          UBXN6        <NA>          lo
    FALSE 3930          UCHL5        <NA>          hi
    FALSE 3931           UCK2        <NA>          hi
    FALSE 3932           UFC1        <NA>          hi
    FALSE 3933           UFL1        <NA>          hi
    FALSE 3934           UFM1        <NA>          hi
    FALSE 3935          UGGT1        <NA>          hi
    FALSE 3936          UGGT2        <NA>          hi
    FALSE 3937          UGT8L        <NA>          lo
    FALSE 3938      UHRF1BP1L        <NA>          hi
    FALSE 3939           ULK3        <NA>          lo
    FALSE 3940         UMODL1        <NA>          hi
    FALSE 3941          UNC5C        <NA>          hi
    FALSE 3942          UNC79        <NA>          hi
    FALSE 3943        UNC93B1        <NA>          lo
    FALSE 3944           URAH        <NA>          lo
    FALSE 3945           URI1        <NA>          lo
    FALSE 3946           USP1        <NA>          hi
    FALSE 3947          USP13        <NA>          hi
    FALSE 3948          USP16        <NA>          hi
    FALSE 3949           USP2        <NA>          lo
    FALSE 3950          USP32        <NA>          lo
    FALSE 3951          USP47        <NA>          hi
    FALSE 3952            UST        <NA>          lo
    FALSE 3953          UTP23        <NA>          hi
    FALSE 3954         VANGL2        <NA>          lo
    FALSE 3955           VAT1        <NA>          lo
    FALSE 3956          VDAC1        <NA>          hi
    FALSE 3957          VEPH1        <NA>          hi
    FALSE 3958          VGLL4        <NA>          lo
    FALSE 3959            VIM        <NA>          lo
    FALSE 3960           VIMP        <NA>          hi
    FALSE 3961          VIPR1        <NA>          hi
    FALSE 3962          VMA21        <NA>          hi
    FALSE 3963           VNN1        <NA>          lo
    FALSE 3964         VPS13A        <NA>          hi
    FALSE 3965          VPS45        <NA>          lo
    FALSE 3966          VPS50        <NA>          hi
    FALSE 3967          VPS53        <NA>          lo
    FALSE 3968          VPS72        <NA>          lo
    FALSE 3969           VRK2        <NA>          lo
    FALSE 3970          VSIG4        <NA>          lo
    FALSE 3971          WBP11        <NA>          lo
    FALSE 3972          WBP1L        <NA>          lo
    FALSE 3973           WBP2        <NA>          lo
    FALSE 3974        WBSCR17        <NA>          lo
    FALSE 3975          WDFY4        <NA>          lo
    FALSE 3976          WDHD1        <NA>          hi
    FALSE 3977          WDR17        <NA>          hi
    FALSE 3978          WDR35        <NA>          hi
    FALSE 3979           WDR4        <NA>          hi
    FALSE 3980          WDR43        <NA>          hi
    FALSE 3981          WDR75        <NA>          hi
    FALSE 3982          WDR82        <NA>          lo
    FALSE 3983          WDR89        <NA>          hi
    FALSE 3984          WDTC1        <NA>          lo
    FALSE 3985          WFDC1        <NA>          lo
    FALSE 3986          WFDC2        <NA>          lo
    FALSE 3987          WIPF3        <NA>          lo
    FALSE 3988          WIPI1        <NA>          hi
    FALSE 3989           WNK1        <NA>          lo
    FALSE 3990          WNT16        <NA>          lo
    FALSE 3991           WNT3        <NA>          lo
    FALSE 3992          WSCD1        <NA>          hi
    FALSE 3993           WWP2        <NA>          lo
    FALSE 3994          WWTR1        <NA>          lo
    FALSE 3995           XAF1        <NA>          lo
    FALSE 3996           XBP1        <NA>          hi
    FALSE 3997           XKR8        <NA>          lo
    FALSE 3998           XPO4        <NA>          hi
    FALSE 3999           XPO6        <NA>          lo
    FALSE 4000          XRCC5        <NA>          hi
    FALSE 4001          XRCC6        <NA>          hi
    FALSE 4002           XRN1        <NA>          hi
    FALSE 4003           YAF2        <NA>          hi
    FALSE 4004           YBX1        <NA>          hi
    FALSE 4005          YIPF1        <NA>          hi
    FALSE 4006          YIPF4        <NA>          hi
    FALSE 4007          YIPF5        <NA>          hi
    FALSE 4008         YME1L1        <NA>          hi
    FALSE 4009          YWHAH        <NA>          hi
    FALSE 4010          YWHAQ        <NA>          hi
    FALSE 4011          ZBED4        <NA>          hi
    FALSE 4012         ZBTB17        <NA>          lo
    FALSE 4013         ZBTB47        <NA>          lo
    FALSE 4014         ZBTB48        <NA>          lo
    FALSE 4015        ZC3H11A        <NA>          lo
    FALSE 4016        ZC3H12A        <NA>          lo
    FALSE 4017          ZC3H6        <NA>          lo
    FALSE 4018        ZC3HAV1        <NA>          lo
    FALSE 4019        ZCCHC14        <NA>          lo
    FALSE 4020          ZCRB1        <NA>          hi
    FALSE 4021         ZDHHC3        <NA>          hi
    FALSE 4022           ZEB2        <NA>          lo
    FALSE 4023        ZFP36L1        <NA>          lo
    FALSE 4024         ZNF148        <NA>          hi
    FALSE 4025         ZNF277        <NA>          hi
    FALSE 4026         ZNF335        <NA>          lo
    FALSE 4027         ZNF362        <NA>          lo
    FALSE 4028        ZNF385B        <NA>          lo
    FALSE 4029         ZNF503        <NA>          lo
    FALSE 4030         ZNF541        <NA>          lo
    FALSE 4031         ZNF644        <NA>          hi
    FALSE 4032         ZNF652        <NA>          lo
    FALSE 4033         ZNF653        <NA>          lo
    FALSE 4034         ZNF710        <NA>          lo
    FALSE 4035         ZNF740        <NA>          lo
    FALSE 4036        ZNF804B        <NA>          hi
    FALSE 4037          ZNFX1        <NA>          lo
    FALSE 4038         ZWILCH        <NA>          hi
    FALSE 4039           ZXDC        <NA>          lo
    FALSE 4040            ZYX        <NA>          lo

    write.csv(sexPRLdesg, "../results/DEG-PRLsex.csv")

    sexPRLdesg2 <- sexPRLdesg %>%
      mutate(higherin = paste(direction.x, direction.y, sep = "_")) %>%
      mutate(direction = paste(higherin, "PRL", sep = "")) %>%
      select(gene, direction) %>%
      group_by(direction) %>%
      summarize(DEGs = str_c(gene , collapse = " "))
    sexPRLdesg2

    FALSE # A tibble: 8 x 2
    FALSE   direction    DEGs                                                        
    FALSE   <chr>        <chr>                                                       
    FALSE 1 female_hiPRL ABCB1LB ABHD2 ACADL ADAM22 ADAM23 ADAMTSL2 AGBL4 AGR2 AHI1 …
    FALSE 2 female_loPRL AMPH ARAP3 ARHGAP24 B4GALT2 BOC BOK C8ORF4 CDH4 CNKSR3 CPNE…
    FALSE 3 female_NAPRL ABCG2 ABHD6 ABLIM3 ACAA2 ACER3 ACKR2 ACOT13 ACOT9 ACPL2 ACS…
    FALSE 4 male_hiPRL   ALDH7A1 ARFRP1 ARL6IP5 CCND2 CDCP1 CEP78 CERK CHRNB4 CKS2 C…
    FALSE 5 male_loPRL   ABCA1 ABCG4 ABHD12 ABHD17C ABR ACER2 ADGRD1 ANXA1 APBB1IP A…
    FALSE 6 male_NAPRL   A2ML4 AAED1 ABAT ABCC5 ABCG1 ABHD17B ACAP3 ACOT11 ACSL4 ADA…
    FALSE 7 NA_hiPRL     ABCC3 ABCC8 ABCC9 ABCE1 ABRA AC025048.1 ACAD9 ACSBG2 ACSL6 …
    FALSE 8 NA_loPRL     ABCA2 ABHD3 ABHD8 ACCS ACHE ACRC ACTN4 ACVR1B ACVRL1 ADAM17…

    write.csv(sexPRLdesg2, "../results/PRLsexinteraction.csv")

    c <- ggplot(sexPRLdesg, aes(x = direction.x, fill = direction.y)) +
      geom_bar(position = "dodge") +
      scale_fill_manual(values = allcolors, na.value= "#bdbdbd",
                        labels= c("lo", "hi", "NS"),
                        name = "PRL") +
      scale_x_discrete(breaks = c("female", "male", NA),
                       labels = c("female", "male", "NS"),
                       name = "Sex")  +
      labs(y = "Total DEGs") +
      theme_B3() +
      theme(legend.position = "top")



    plot_grid(ab, c, nrow = 1, rel_widths = c(2,1.2), labels = c(" ", "d"), label_size = 12)

![](../figures/DESeq2/PRL-5.png)

    head(sexPRLdesg)

    FALSE      gene direction.x direction.y
    FALSE 1   A2ML4        male        <NA>
    FALSE 2   AAED1        male        <NA>
    FALSE 3    ABAT        male        <NA>
    FALSE 4   ABCA1        male          lo
    FALSE 5 ABCB1LB      female          hi
    FALSE 6   ABCC5        male        <NA>

    PRLhi <- sexPRLdesg %>%
      filter(direction.y == "hi")
    PRLhi

    FALSE                gene direction.x direction.y
    FALSE 1           ABCB1LB      female          hi
    FALSE 2             ABHD2      female          hi
    FALSE 3             ACADL      female          hi
    FALSE 4            ADAM22      female          hi
    FALSE 5            ADAM23      female          hi
    FALSE 6          ADAMTSL2      female          hi
    FALSE 7             AGBL4      female          hi
    FALSE 8              AGR2      female          hi
    FALSE 9              AHI1      female          hi
    FALSE 10            AIMP1      female          hi
    FALSE 11          ALDH7A1        male          hi
    FALSE 12            AMY1A      female          hi
    FALSE 13           ARFRP1        male          hi
    FALSE 14           ARID4B      female          hi
    FALSE 15          ARL6IP5        male          hi
    FALSE 16            ATG4A      female          hi
    FALSE 17              ATM      female          hi
    FALSE 18         ATP6V1C2      female          hi
    FALSE 19            BCAT1      female          hi
    FALSE 20             BIVM      female          hi
    FALSE 21      C28H19ORF10      female          hi
    FALSE 22           CACFD1      female          hi
    FALSE 23           CAMTA1      female          hi
    FALSE 24             CBX6      female          hi
    FALSE 25           CCDC60      female          hi
    FALSE 26            CCND2        male          hi
    FALSE 27            CDCP1        male          hi
    FALSE 28           CDKN1A      female          hi
    FALSE 29            CEP78        male          hi
    FALSE 30            CEP83      female          hi
    FALSE 31             CERK        male          hi
    FALSE 32             CHKA      female          hi
    FALSE 33           CHRNB4        male          hi
    FALSE 34           CITED4      female          hi
    FALSE 35             CKS2        male          hi
    FALSE 36           CLDN11      female          hi
    FALSE 37             CLTA        male          hi
    FALSE 38             CMAS      female          hi
    FALSE 39             CNST      female          hi
    FALSE 40            CNTFR      female          hi
    FALSE 41            CNTLN        male          hi
    FALSE 42            COPB2      female          hi
    FALSE 43             COPE      female          hi
    FALSE 44           COPS7B      female          hi
    FALSE 45             COQ9      female          hi
    FALSE 46            CPEB2      female          hi
    FALSE 47             CREM        male          hi
    FALSE 48             CTBS      female          hi
    FALSE 49     CTC-487M23.8        male          hi
    FALSE 50            CTSL2      female          hi
    FALSE 51        CZH5ORF51        male          hi
    FALSE 52            DDOST      female          hi
    FALSE 53            DERL2      female          hi
    FALSE 54           DIRAS1        male          hi
    FALSE 55           DNAJC1      female          hi
    FALSE 56             DPH6      female          hi
    FALSE 57             DPP7      female          hi
    FALSE 58           DUSP10        male          hi
    FALSE 59           EEF1A1      female          hi
    FALSE 60             EEF2      female          hi
    FALSE 61            EIF3L      female          hi
    FALSE 62           EIF4G2        male          hi
    FALSE 63              EMB        male          hi
    FALSE 64             ENO1      female          hi
    FALSE 65             ENO2      female          hi
    FALSE 66           ENTPD4        male          hi
    FALSE 67            EPHX4      female          hi
    FALSE 68            ESRP1      female          hi
    FALSE 69             ETFA      female          hi
    FALSE 70             EXT2      female          hi
    FALSE 71            F13A1        male          hi
    FALSE 72         FAM114A1      female          hi
    FALSE 73          FAM135B      female          hi
    FALSE 74           FAM73A      female          hi
    FALSE 75           FAM83B      female          hi
    FALSE 76             FAN1      female          hi
    FALSE 77              FER        male          hi
    FALSE 78             FGF6        male          hi
    FALSE 79            FOCAD        male          hi
    FALSE 80            FOSL2        male          hi
    FALSE 81              FTL        male          hi
    FALSE 82             GAB2        male          hi
    FALSE 83            GADL1      female          hi
    FALSE 84             GCC2      female          hi
    FALSE 85             GCH1        male          hi
    FALSE 86            GDPD1        male          hi
    FALSE 87               GK        male          hi
    FALSE 88            GMPPB      female          hi
    FALSE 89            GNA11      female          hi
    FALSE 90          GORASP2      female          hi
    FALSE 91             GOT2      female          hi
    FALSE 92             GPC6      female          hi
    FALSE 93            GRB10      female          hi
    FALSE 94           GRIN3A        male          hi
    FALSE 95            GRIP1        male          hi
    FALSE 96            HMGCR        male          hi
    FALSE 97            HOOK3        male          hi
    FALSE 98           HS3ST5      female          hi
    FALSE 99           HS6ST3        male          hi
    FALSE 100         HSD17B4        male          hi
    FALSE 101            ICA1      female          hi
    FALSE 102           INSM1        male          hi
    FALSE 103           ISCA1        male          hi
    FALSE 104        IVNS1ABP        male          hi
    FALSE 105          KCTD14        male          hi
    FALSE 106          KDELR3      female          hi
    FALSE 107            LAP3      female          hi
    FALSE 108           LMAN1      female          hi
    FALSE 109           LMNB1        male          hi
    FALSE 110    LOC100858130        male          hi
    FALSE 111    LOC107050724      female          hi
    FALSE 112    LOC107053040        male          hi
    FALSE 113    LOC107054798        male          hi
    FALSE 114       LOC378902      female          hi
    FALSE 115       LOC395159      female          hi
    FALSE 116       LOC417253      female          hi
    FALSE 117       LOC421106      female          hi
    FALSE 118       LOC422926      female          hi
    FALSE 119       LOC428479        male          hi
    FALSE 120           LRP12        male          hi
    FALSE 121          LRPPRC      female          hi
    FALSE 122          LRRIQ1      female          hi
    FALSE 123           MAK16      female          hi
    FALSE 124           MANBA      female          hi
    FALSE 125          MAP7D3      female          hi
    FALSE 126            MELK        male          hi
    FALSE 127            MEST      female          hi
    FALSE 128          METRNL      female          hi
    FALSE 129          METTL9      female          hi
    FALSE 130           MFSD4        male          hi
    FALSE 131          MGAT4C      female          hi
    FALSE 132          MRPL48      female          hi
    FALSE 133            MSH3        male          hi
    FALSE 134           MTMR1        male          hi
    FALSE 135             MYC        male          hi
    FALSE 136          N6AMT2      female          hi
    FALSE 137        NAALADL2      female          hi
    FALSE 138            NCS1        male          hi
    FALSE 139           NDST4        male          hi
    FALSE 140            NEBL        male          hi
    FALSE 141          NECAB2      female          hi
    FALSE 142            NOL8      female          hi
    FALSE 143            NPC2      female          hi
    FALSE 144           NR4A3        male          hi
    FALSE 145           NRSN1      female          hi
    FALSE 146          NUDCD1      female          hi
    FALSE 147            PAK6        male          hi
    FALSE 148            PARG      female          hi
    FALSE 149            PAX7        male          hi
    FALSE 150           PCF11      female          hi
    FALSE 151          PDCD11      female          hi
    FALSE 152           PDE3B        male          hi
    FALSE 153           PDE4D        male          hi
    FALSE 154           PDE8A        male          hi
    FALSE 155           PDE8B        male          hi
    FALSE 156           PDGFD      female          hi
    FALSE 157           PDIA6      female          hi
    FALSE 158            PDK4      female          hi
    FALSE 159            PGM3      female          hi
    FALSE 160          PIK3C3        male          hi
    FALSE 161            PIM3        male          hi
    FALSE 162           PLCB1      female          hi
    FALSE 163           PLCB2      female          hi
    FALSE 164           POLA1      female          hi
    FALSE 165          POLR3F      female          hi
    FALSE 166            PPIC        male          hi
    FALSE 167           PRDM5      female          hi
    FALSE 168           PSMC3      female          hi
    FALSE 169           PTPN5        male          hi
    FALSE 170         PTTG1IP      female          hi
    FALSE 171            QARS      female          hi
    FALSE 172            QPCT      female          hi
    FALSE 173           QRSL1        male          hi
    FALSE 174            RARS      female          hi
    FALSE 175           RASA4      female          hi
    FALSE 176          RASAL1      female          hi
    FALSE 177         RASL10A      female          hi
    FALSE 178         RASL11A      female          hi
    FALSE 179          RASSF5        male          hi
    FALSE 180           RBM25      female          hi
    FALSE 181           RCAN2        male          hi
    FALSE 182           REEP1      female          hi
    FALSE 183         RHOBTB3        male          hi
    FALSE 184           RIPK3      female          hi
    FALSE 185           RNF32      female          hi
    FALSE 186   RP11-296A16.1        male          hi
    FALSE 187  RP11-574K11.31        male          hi
    FALSE 188            RPA1        male          hi
    FALSE 189           RPL15      female          hi
    FALSE 190          RPL27A      female          hi
    FALSE 191           RPL30      female          hi
    FALSE 192           RPL39      female          hi
    FALSE 193            RPL4      female          hi
    FALSE 194            RPL5      female          hi
    FALSE 195            RPL6      female          hi
    FALSE 196            RPN2      female          hi
    FALSE 197            RPS2      female          hi
    FALSE 198          RPS27A      female          hi
    FALSE 199           RPS3A      female          hi
    FALSE 200            RPS7      female          hi
    FALSE 201            RPS8      female          hi
    FALSE 202            RPSA      female          hi
    FALSE 203           RRBP1      female          hi
    FALSE 204           SCFD2      female          hi
    FALSE 205         SEC14L5        male          hi
    FALSE 206           SEC63      female          hi
    FALSE 207           SEL1L      female          hi
    FALSE 208          SEMA3A      female          hi
    FALSE 209            SHC3        male          hi
    FALSE 210        SLC16A14      female          hi
    FALSE 211        SLC25A20      female          hi
    FALSE 212        SLC25A29      female          hi
    FALSE 213         SLC35F5      female          hi
    FALSE 214         SLC41A3      female          hi
    FALSE 215         SMARCA5      female          hi
    FALSE 216            SMC2        male          hi
    FALSE 217            SMC5        male          hi
    FALSE 218           SMYD1      female          hi
    FALSE 219           SPRY1        male          hi
    FALSE 220         SPTY2D1      female          hi
    FALSE 221           SSTR1      female          hi
    FALSE 222            ST13      female          hi
    FALSE 223      ST6GALNAC5      female          hi
    FALSE 224           STIM2      female          hi
    FALSE 225           STT3B      female          hi
    FALSE 226           SUMF2      female          hi
    FALSE 227            SV2B        male          hi
    FALSE 228            SYN3        male          hi
    FALSE 229           SYNE2      female          hi
    FALSE 230           SYT14      female          hi
    FALSE 231            SYT4        male          hi
    FALSE 232         TBC1D30        male          hi
    FALSE 233            TBL2      female          hi
    FALSE 234            TET1      female          hi
    FALSE 235           TEX10      female          hi
    FALSE 236           TFB2M      female          hi
    FALSE 237            THRB        male          hi
    FALSE 238           TMCO1      female          hi
    FALSE 239           TMED7      female          hi
    FALSE 240        TMEM167A        male          hi
    FALSE 241        TMEM255B      female          hi
    FALSE 242         TMEM258      female          hi
    FALSE 243           TNPO1        male          hi
    FALSE 244           TRPA1        male          hi
    FALSE 245           TRPV6      female          hi
    FALSE 246           TULP3      female          hi
    FALSE 247          UNC119        male          hi
    FALSE 248          UNC13C        male          hi
    FALSE 249            UNKL      female          hi
    FALSE 250          UQCRC2      female          hi
    FALSE 251            USO1      female          hi
    FALSE 252           USP35        male          hi
    FALSE 253             VCP      female          hi
    FALSE 254            WDR3      female          hi
    FALSE 255            WSB2      female          hi
    FALSE 256           WSCD2      female          hi
    FALSE 257            XPOT      female          hi
    FALSE 258        ZC3HAV1L      female          hi
    FALSE 259          ZNF367        male          hi
    FALSE 260           ABCC3        <NA>          hi
    FALSE 261           ABCC8        <NA>          hi
    FALSE 262           ABCC9        <NA>          hi
    FALSE 263           ABCE1        <NA>          hi
    FALSE 264            ABRA        <NA>          hi
    FALSE 265      AC025048.1        <NA>          hi
    FALSE 266           ACAD9        <NA>          hi
    FALSE 267          ACSBG2        <NA>          hi
    FALSE 268           ACSL6        <NA>          hi
    FALSE 269           ACSS2        <NA>          hi
    FALSE 270          ACTR10        <NA>          hi
    FALSE 271           ACTR2        <NA>          hi
    FALSE 272          ACVR2B        <NA>          hi
    FALSE 273           ACYP2        <NA>          hi
    FALSE 274        ADAMTS15        <NA>          hi
    FALSE 275         ADAMTS8        <NA>          hi
    FALSE 276        ADAMTSL1        <NA>          hi
    FALSE 277          ADARB1        <NA>          hi
    FALSE 278          ADARB2        <NA>          hi
    FALSE 279             ADK        <NA>          hi
    FALSE 280            ADSL        <NA>          hi
    FALSE 281            ADSS        <NA>          hi
    FALSE 282            AFF2        <NA>          hi
    FALSE 283            AIDA        <NA>          hi
    FALSE 284           AIFM1        <NA>          hi
    FALSE 285           AIMP2        <NA>          hi
    FALSE 286           AKAP6        <NA>          hi
    FALSE 287           AKAP9        <NA>          hi
    FALSE 288        ALDH18A1        <NA>          hi
    FALSE 289         ALDH1L2        <NA>          hi
    FALSE 290         ALDH6A1        <NA>          hi
    FALSE 291            ALG1        <NA>          hi
    FALSE 292           ALG12        <NA>          hi
    FALSE 293            ALG2        <NA>          hi
    FALSE 294            ALG5        <NA>          hi
    FALSE 295            ALG8        <NA>          hi
    FALSE 296          ALKBH8        <NA>          hi
    FALSE 297          AMIGO2        <NA>          hi
    FALSE 298          ANAPC4        <NA>          hi
    FALSE 299         ANKRD12        <NA>          hi
    FALSE 300         ANKRD42        <NA>          hi
    FALSE 301            ANLN        <NA>          hi
    FALSE 302            ANO3        <NA>          hi
    FALSE 303           AP1AR        <NA>          hi
    FALSE 304           AP1S3        <NA>          hi
    FALSE 305           APAF1        <NA>          hi
    FALSE 306            API5        <NA>          hi
    FALSE 307             APP        <NA>          hi
    FALSE 308           AQP11        <NA>          hi
    FALSE 309           ARAP2        <NA>          hi
    FALSE 310            ARF1        <NA>          hi
    FALSE 311            ARF4        <NA>          hi
    FALSE 312         ARFGAP1        <NA>          hi
    FALSE 313         ARFGAP3        <NA>          hi
    FALSE 314         ARFGEF3        <NA>          hi
    FALSE 315          ARGLU1        <NA>          hi
    FALSE 316        ARHGAP19        <NA>          hi
    FALSE 317         ARHGEF3        <NA>          hi
    FALSE 318            ARL1        <NA>          hi
    FALSE 319          ARL2BP        <NA>          hi
    FALSE 320           ARL5B        <NA>          hi
    FALSE 321           ARL8A        <NA>          hi
    FALSE 322          ARMC10        <NA>          hi
    FALSE 323           ASB14        <NA>          hi
    FALSE 324            ASB9        <NA>          hi
    FALSE 325            ASPH        <NA>          hi
    FALSE 326            ASPM        <NA>          hi
    FALSE 327           ATG2B        <NA>          hi
    FALSE 328            ATG7        <NA>          hi
    FALSE 329         ATP13A3        <NA>          hi
    FALSE 330           ATP4B        <NA>          hi
    FALSE 331          ATP5C1        <NA>          hi
    FALSE 332          ATP5G3        <NA>          hi
    FALSE 333           ATP5I        <NA>          hi
    FALSE 334           ATP5J        <NA>          hi
    FALSE 335         ATP6AP2        <NA>          hi
    FALSE 336        ATP6V0D1        <NA>          hi
    FALSE 337        ATP6V0E1        <NA>          hi
    FALSE 338         ATP6V1A        <NA>          hi
    FALSE 339        ATP6V1E1        <NA>          hi
    FALSE 340          ATP8A1        <NA>          hi
    FALSE 341             ATR        <NA>          hi
    FALSE 342            ATRX        <NA>          hi
    FALSE 343           AURKA        <NA>          hi
    FALSE 344           AZIN1        <NA>          hi
    FALSE 345        B3GALNT2        <NA>          hi
    FALSE 346         B3GALT2        <NA>          hi
    FALSE 347           BARD1        <NA>          hi
    FALSE 348           BASP1        <NA>          hi
    FALSE 349           BAZ1A        <NA>          hi
    FALSE 350            BBS9        <NA>          hi
    FALSE 351             BBX        <NA>          hi
    FALSE 352          BCKDHB        <NA>          hi
    FALSE 353            BCL2        <NA>          hi
    FALSE 354           BECN1        <NA>          hi
    FALSE 355            BET1        <NA>          hi
    FALSE 356           BICC1        <NA>          hi
    FALSE 357           BICD1        <NA>          hi
    FALSE 358           BIRC5        <NA>          hi
    FALSE 359            BLMH        <NA>          hi
    FALSE 360           BLZF1        <NA>          hi
    FALSE 361           BNIP1        <NA>          hi
    FALSE 362          BNIP3L        <NA>          hi
    FALSE 363          BOD1L1        <NA>          hi
    FALSE 364           BOLA3        <NA>          hi
    FALSE 365            BORA        <NA>          hi
    FALSE 366           BRCA1        <NA>          hi
    FALSE 367           BRIP1        <NA>          hi
    FALSE 368           BRWD3        <NA>          hi
    FALSE 369           BTBD1        <NA>          hi
    FALSE 370          BTBD10        <NA>          hi
    FALSE 371           BTBD9        <NA>          hi
    FALSE 372          BTF3L4        <NA>          hi
    FALSE 373            BTG4        <NA>          hi
    FALSE 374            BUB1        <NA>          hi
    FALSE 375           BUB1B        <NA>          hi
    FALSE 376            BUB3        <NA>          hi
    FALSE 377            BZW1        <NA>          hi
    FALSE 378     C10H15ORF40        <NA>          hi
    FALSE 379        C12ORF29        <NA>          hi
    FALSE 380     C14H16ORF52        <NA>          hi
    FALSE 381     C14H16ORF88        <NA>          hi
    FALSE 382         C14ORF2        <NA>          hi
    FALSE 383        C16ORF45        <NA>          hi
    FALSE 384        C18ORF21        <NA>          hi
    FALSE 385        C18ORF32        <NA>          hi
    FALSE 386      C1H12ORF23        <NA>          hi
    FALSE 387       C1H12ORF4        <NA>          hi
    FALSE 388       C1H12ORF5        <NA>          hi
    FALSE 389      C1H12ORF66        <NA>          hi
    FALSE 390       C2H6ORF62        <NA>          hi
    FALSE 391       C2H7ORF63        <NA>          hi
    FALSE 392       C2H9ORF30        <NA>          hi
    FALSE 393         C2ORF47        <NA>          hi
    FALSE 394         C3ORF58        <NA>          hi
    FALSE 395           C4BPA        <NA>          hi
    FALSE 396       C4H4ORF21        <NA>          hi
    FALSE 397       C4HXORF57        <NA>          hi
    FALSE 398     C5H14ORF166        <NA>          hi
    FALSE 399      C6H10ORF71        <NA>          hi
    FALSE 400         C7ORF50        <NA>          hi
    FALSE 401       C8H1ORF21        <NA>          hi
    FALSE 402       C8H1ORF27        <NA>          hi
    FALSE 403       C9H3ORF55        <NA>          hi
    FALSE 404          CAB39L        <NA>          hi
    FALSE 405         CABLES2        <NA>          hi
    FALSE 406           CABP1        <NA>          hi
    FALSE 407        CACNA2D1        <NA>          hi
    FALSE 408           CALN1        <NA>          hi
    FALSE 409          CAMKMT        <NA>          hi
    FALSE 410           CAMKV        <NA>          hi
    FALSE 411            CANX        <NA>          hi
    FALSE 412            CARS        <NA>          hi
    FALSE 413           CASC5        <NA>          hi
    FALSE 414            CASK        <NA>          hi
    FALSE 415        CASP8AP2        <NA>          hi
    FALSE 416             CBS        <NA>          hi
    FALSE 417           CCAR1        <NA>          hi
    FALSE 418         CCDC167        <NA>          hi
    FALSE 419          CCDC58        <NA>          hi
    FALSE 420         CCDC90B        <NA>          hi
    FALSE 421          CCDC91        <NA>          hi
    FALSE 422           CCNA2        <NA>          hi
    FALSE 423           CCNB3        <NA>          hi
    FALSE 424            CCNC        <NA>          hi
    FALSE 425           CCNE2        <NA>          hi
    FALSE 426           CCNT2        <NA>          hi
    FALSE 427          CCSER2        <NA>          hi
    FALSE 428            CCT2        <NA>          hi
    FALSE 429            CCT8        <NA>          hi
    FALSE 430           CD164        <NA>          hi
    FALSE 431            CD24        <NA>          hi
    FALSE 432           CDC20        <NA>          hi
    FALSE 433          CDC25A        <NA>          hi
    FALSE 434           CDC45        <NA>          hi
    FALSE 435            CDC6        <NA>          hi
    FALSE 436           CDCA3        <NA>          hi
    FALSE 437           CDCA7        <NA>          hi
    FALSE 438           CDH18        <NA>          hi
    FALSE 439            CDK1        <NA>          hi
    FALSE 440        CDK5RAP2        <NA>          hi
    FALSE 441        CDK5RAP3        <NA>          hi
    FALSE 442            CDK6        <NA>          hi
    FALSE 443            CDK8        <NA>          hi
    FALSE 444           CDKN3        <NA>          hi
    FALSE 445            CDR2        <NA>          hi
    FALSE 446            CDT1        <NA>          hi
    FALSE 447           CENPA        <NA>          hi
    FALSE 448           CENPC        <NA>          hi
    FALSE 449           CENPE        <NA>          hi
    FALSE 450           CENPF        <NA>          hi
    FALSE 451           CENPI        <NA>          hi
    FALSE 452           CENPK        <NA>          hi
    FALSE 453           CENPM        <NA>          hi
    FALSE 454           CENPN        <NA>          hi
    FALSE 455           CENPO        <NA>          hi
    FALSE 456          CEP128        <NA>          hi
    FALSE 457          CEP152        <NA>          hi
    FALSE 458          CEP170        <NA>          hi
    FALSE 459          CEP192        <NA>          hi
    FALSE 460          CEP290        <NA>          hi
    FALSE 461           CEP55        <NA>          hi
    FALSE 462           CEP95        <NA>          hi
    FALSE 463           CERS6        <NA>          hi
    FALSE 464           CETN1        <NA>          hi
    FALSE 465          CGRRF1        <NA>          hi
    FALSE 466          CHCHD3        <NA>          hi
    FALSE 467          CHCHD4        <NA>          hi
    FALSE 468          CHCHD7        <NA>          hi
    FALSE 469           CHD1L        <NA>          hi
    FALSE 470            CHPF        <NA>          hi
    FALSE 471          CHRNA5        <NA>          hi
    FALSE 472          CHST10        <NA>          hi
    FALSE 473           CHSY1        <NA>          hi
    FALSE 474             CIT        <NA>          hi
    FALSE 475           CKAP2        <NA>          hi
    FALSE 476           CKAP5        <NA>          hi
    FALSE 477           CKS1B        <NA>          hi
    FALSE 478          CLASP2        <NA>          hi
    FALSE 479          CLEC3A        <NA>          hi
    FALSE 480            CLGN        <NA>          hi
    FALSE 481           CLIP4        <NA>          hi
    FALSE 482          CLNS1A        <NA>          hi
    FALSE 483         CLPTM1L        <NA>          hi
    FALSE 484           CLSPN        <NA>          hi
    FALSE 485          CLTCL1        <NA>          hi
    FALSE 486           CLVS2        <NA>          hi
    FALSE 487           CMSS1        <NA>          hi
    FALSE 488           CNIH1        <NA>          hi
    FALSE 489          CNPPD1        <NA>          hi
    FALSE 490          CNRIP1        <NA>          hi
    FALSE 491         COL19A1        <NA>          hi
    FALSE 492         COL20A1        <NA>          hi
    FALSE 493          COL4A3        <NA>          hi
    FALSE 494            COPA        <NA>          hi
    FALSE 495           COPS5        <NA>          hi
    FALSE 496           COPS8        <NA>          hi
    FALSE 497            COQ2        <NA>          hi
    FALSE 498            COQ5        <NA>          hi
    FALSE 499          COX4I1        <NA>          hi
    FALSE 500          COX7A2        <NA>          hi
    FALSE 501             CPE        <NA>          hi
    FALSE 502           CPLX4        <NA>          hi
    FALSE 503           CREB1        <NA>          hi
    FALSE 504         CREB3L1        <NA>          hi
    FALSE 505         CREB3L2        <NA>          hi
    FALSE 506          CRELD2        <NA>          hi
    FALSE 507           CRTAP        <NA>          hi
    FALSE 508           CSDE1        <NA>          hi
    FALSE 509      CSGALNACT2        <NA>          hi
    FALSE 510           CSTF3        <NA>          hi
    FALSE 511           CTBP1        <NA>          hi
    FALSE 512   CTD-2410N18.5        <NA>          hi
    FALSE 513          CTNNA3        <NA>          hi
    FALSE 514        CTNNBIP1        <NA>          hi
    FALSE 515           CTPS1        <NA>          hi
    FALSE 516            CTR9        <NA>          hi
    FALSE 517            CUTC        <NA>          hi
    FALSE 518            CUX1        <NA>          hi
    FALSE 519           CWC22        <NA>          hi
    FALSE 520           CXADR        <NA>          hi
    FALSE 521           DACH1        <NA>          hi
    FALSE 522           DACH2        <NA>          hi
    FALSE 523             DCK        <NA>          hi
    FALSE 524             DCX        <NA>          hi
    FALSE 525            DCXR        <NA>          hi
    FALSE 526          DDRGK1        <NA>          hi
    FALSE 527            DDX1        <NA>          hi
    FALSE 528           DDX24        <NA>          hi
    FALSE 529           DDX46        <NA>          hi
    FALSE 530           DDX59        <NA>          hi
    FALSE 531         DENND1B        <NA>          hi
    FALSE 532          DEPDC1        <NA>          hi
    FALSE 533         DEPDC1B        <NA>          hi
    FALSE 534          DEPDC6        <NA>          hi
    FALSE 535            DERA        <NA>          hi
    FALSE 536           DERL1        <NA>          hi
    FALSE 537           DESI1        <NA>          hi
    FALSE 538          DFNB31        <NA>          hi
    FALSE 539            DGKH        <NA>          hi
    FALSE 540            DGKI        <NA>          hi
    FALSE 541           DGUOK        <NA>          hi
    FALSE 542           DHX36        <NA>          hi
    FALSE 543          DIAPH1        <NA>          hi
    FALSE 544          DIAPH3        <NA>          hi
    FALSE 545            DIS3        <NA>          hi
    FALSE 546            DKK2        <NA>          hi
    FALSE 547          DLGAP5        <NA>          hi
    FALSE 548            DNA2        <NA>          hi
    FALSE 549         DNAJB11        <NA>          hi
    FALSE 550          DNAJB9        <NA>          hi
    FALSE 551         DNAJC12        <NA>          hi
    FALSE 552          DNAJC2        <NA>          hi
    FALSE 553           DOCK5        <NA>          hi
    FALSE 554            DOK5        <NA>          hi
    FALSE 555            DPH5        <NA>          hi
    FALSE 556            DPH7        <NA>          hi
    FALSE 557            DPM1        <NA>          hi
    FALSE 558           DPP10        <NA>          hi
    FALSE 559           DPY30        <NA>          hi
    FALSE 560          DPYSL3        <NA>          hi
    FALSE 561             DR1        <NA>          hi
    FALSE 562            DRG1        <NA>          hi
    FALSE 563            DSC2        <NA>          hi
    FALSE 564            DTD1        <NA>          hi
    FALSE 565             DTL        <NA>          hi
    FALSE 566          DTNBP1        <NA>          hi
    FALSE 567            DUS2        <NA>          hi
    FALSE 568           DUSP4        <NA>          hi
    FALSE 569             DUT        <NA>          hi
    FALSE 570        DYNC1LI1        <NA>          hi
    FALSE 571        DYNC2LI1        <NA>          hi
    FALSE 572           DZIP1        <NA>          hi
    FALSE 573            E2F1        <NA>          hi
    FALSE 574           EBAG9        <NA>          hi
    FALSE 575            ECT2        <NA>          hi
    FALSE 576           EDEM1        <NA>          hi
    FALSE 577           EDEM2        <NA>          hi
    FALSE 578           EDEM3        <NA>          hi
    FALSE 579          EEF1B2        <NA>          hi
    FALSE 580           EFR3A        <NA>          hi
    FALSE 581           EFR3B        <NA>          hi
    FALSE 582          EIF1AY        <NA>          hi
    FALSE 583         EIF2AK1        <NA>          hi
    FALSE 584         EIF2AK3        <NA>          hi
    FALSE 585         EIF2AK4        <NA>          hi
    FALSE 586          EIF2S1        <NA>          hi
    FALSE 587          EIF2S2        <NA>          hi
    FALSE 588         EIF2S3L        <NA>          hi
    FALSE 589           EIF3B        <NA>          hi
    FALSE 590           EIF3D        <NA>          hi
    FALSE 591           EIF3E        <NA>          hi
    FALSE 592           EIF3M        <NA>          hi
    FALSE 593           EIF4H        <NA>          hi
    FALSE 594            EIF5        <NA>          hi
    FALSE 595           EIF5B        <NA>          hi
    FALSE 596          ELAVL2        <NA>          hi
    FALSE 597          ELAVL4        <NA>          hi
    FALSE 598          ELMOD2        <NA>          hi
    FALSE 599          ELOVL4        <NA>          hi
    FALSE 600            EMC3        <NA>          hi
    FALSE 601            EME1        <NA>          hi
    FALSE 602          ENDOD1        <NA>          hi
    FALSE 603          ENTPD5        <NA>          hi
    FALSE 604          ENTPD7        <NA>          hi
    FALSE 605            EPRS        <NA>          hi
    FALSE 606            EPS8        <NA>          hi
    FALSE 607           ERCC6        <NA>          hi
    FALSE 608          ERGIC1        <NA>          hi
    FALSE 609             ERH        <NA>          hi
    FALSE 610            ERI2        <NA>          hi
    FALSE 611          ERLEC1        <NA>          hi
    FALSE 612           ERO1L        <NA>          hi
    FALSE 613          ERO1LB        <NA>          hi
    FALSE 614           ERP29        <NA>          hi
    FALSE 615           ERP44        <NA>          hi
    FALSE 616           ESCO2        <NA>          hi
    FALSE 617            ESF1        <NA>          hi
    FALSE 618            ETF1        <NA>          hi
    FALSE 619            ETV1        <NA>          hi
    FALSE 620            EXO1        <NA>          hi
    FALSE 621           EXOC4        <NA>          hi
    FALSE 622          EXOSC7        <NA>          hi
    FALSE 623         FAM110B        <NA>          hi
    FALSE 624         FAM175B        <NA>          hi
    FALSE 625         FAM184B        <NA>          hi
    FALSE 626         FAM18B1        <NA>          hi
    FALSE 627         FAM19A1        <NA>          hi
    FALSE 628         FAM213A        <NA>          hi
    FALSE 629          FAM46A        <NA>          hi
    FALSE 630          FAM46D        <NA>          hi
    FALSE 631          FAM83G        <NA>          hi
    FALSE 632          FAM84A        <NA>          hi
    FALSE 633          FAM98A        <NA>          hi
    FALSE 634          FANCD2        <NA>          hi
    FALSE 635           FANCI        <NA>          hi
    FALSE 636           FANCL        <NA>          hi
    FALSE 637            FAR1        <NA>          hi
    FALSE 638           FARSB        <NA>          hi
    FALSE 639         FASTKD1        <NA>          hi
    FALSE 640           FBXO5        <NA>          hi
    FALSE 641           FBXW2        <NA>          hi
    FALSE 642           FFAR4        <NA>          hi
    FALSE 643            FGD4        <NA>          hi
    FALSE 644           FHOD3        <NA>          hi
    FALSE 645            FICD        <NA>          hi
    FALSE 646          FKBP11        <NA>          hi
    FALSE 647           FKBP4        <NA>          hi
    FALSE 648          FNDC3A        <NA>          hi
    FALSE 649           FOXM1        <NA>          hi
    FALSE 650          FRRS1L        <NA>          hi
    FALSE 651           FTSJ2        <NA>          hi
    FALSE 652            FUT9        <NA>          hi
    FALSE 653           FXYD2        <NA>          hi
    FALSE 654            G2E3        <NA>          hi
    FALSE 655           G3BP2        <NA>          hi
    FALSE 656            GAB1        <NA>          hi
    FALSE 657       GABARAPL1        <NA>          hi
    FALSE 658           GABRE        <NA>          hi
    FALSE 659         GAL3ST2        <NA>          hi
    FALSE 660            GALE        <NA>          hi
    FALSE 661          GALNT7        <NA>          hi
    FALSE 662            GARS        <NA>          hi
    FALSE 663          GAS2L3        <NA>          hi
    FALSE 664            GCLM        <NA>          hi
    FALSE 665            GCSH        <NA>          hi
    FALSE 666           GFRAL        <NA>          hi
    FALSE 667            GGCT        <NA>          hi
    FALSE 668             GGH        <NA>          hi
    FALSE 669            GID8        <NA>          hi
    FALSE 670           GINS1        <NA>          hi
    FALSE 671           GINS2        <NA>          hi
    FALSE 672            GLB1        <NA>          hi
    FALSE 673            GLG1        <NA>          hi
    FALSE 674          GLIPR2        <NA>          hi
    FALSE 675            GLMN        <NA>          hi
    FALSE 676           GLRX3        <NA>          hi
    FALSE 677          GLT8D1        <NA>          hi
    FALSE 678            GMNN        <NA>          hi
    FALSE 679            GMPS        <NA>          hi
    FALSE 680            GNB1        <NA>          hi
    FALSE 681            GNG4        <NA>          hi
    FALSE 682            GNL2        <NA>          hi
    FALSE 683            GNL3        <NA>          hi
    FALSE 684           GNPTG        <NA>          hi
    FALSE 685          GOLGA2        <NA>          hi
    FALSE 686          GOLGA3        <NA>          hi
    FALSE 687          GOLGA4        <NA>          hi
    FALSE 688          GOLGA5        <NA>          hi
    FALSE 689          GOLT1B        <NA>          hi
    FALSE 690            GOPC        <NA>          hi
    FALSE 691            GOT1        <NA>          hi
    FALSE 692            GPHN        <NA>          hi
    FALSE 693             GPI        <NA>          hi
    FALSE 694           GPNMB        <NA>          hi
    FALSE 695          GPR123        <NA>          hi
    FALSE 696          GPR180        <NA>          hi
    FALSE 697           GPR25        <NA>          hi
    FALSE 698          GPRIN3        <NA>          hi
    FALSE 699            GPX4        <NA>          hi
    FALSE 700            GPX7        <NA>          hi
    FALSE 701           GRID1        <NA>          hi
    FALSE 702           GRID2        <NA>          hi
    FALSE 703           GRIK2        <NA>          hi
    FALSE 704          GRXCR1        <NA>          hi
    FALSE 705           GSPT1        <NA>          hi
    FALSE 706          GTF3C6        <NA>          hi
    FALSE 707          GTPBP4        <NA>          hi
    FALSE 708         GUCY1B4        <NA>          hi
    FALSE 709          GXYLT2        <NA>          hi
    FALSE 710            HAT1        <NA>          hi
    FALSE 711           HAUS1        <NA>          hi
    FALSE 712           HBS1L        <NA>          hi
    FALSE 713           HDLBP        <NA>          hi
    FALSE 714          HEATR3        <NA>          hi
    FALSE 715           HELLS        <NA>          hi
    FALSE 716         HERPUD1        <NA>          hi
    FALSE 717           HJURP        <NA>          hi
    FALSE 718            HM13        <NA>          hi
    FALSE 719           HMGB2        <NA>          hi
    FALSE 720           HMGN1        <NA>          hi
    FALSE 721           HMGN5        <NA>          hi
    FALSE 722            HMMR        <NA>          hi
    FALSE 723       HNRNPA2B1        <NA>          hi
    FALSE 724         HNRNPAB        <NA>          hi
    FALSE 725         HNRNPH1        <NA>          hi
    FALSE 726          HOMER2        <NA>          hi
    FALSE 727           HOOK1        <NA>          hi
    FALSE 728          HPCAL1        <NA>          hi
    FALSE 729           HPRT1        <NA>          hi
    FALSE 730            HRH3        <NA>          hi
    FALSE 731          HS6ST2        <NA>          hi
    FALSE 732        HSD17B12        <NA>          hi
    FALSE 733         HSP90B1        <NA>          hi
    FALSE 734          HSPA13        <NA>          hi
    FALSE 735          HSPA14        <NA>          hi
    FALSE 736           HSPA5        <NA>          hi
    FALSE 737         HTATSF1        <NA>          hi
    FALSE 738           HYOU1        <NA>          hi
    FALSE 739            IARS        <NA>          hi
    FALSE 740           IARS2        <NA>          hi
    FALSE 741            IBTK        <NA>          hi
    FALSE 742             IDE        <NA>          hi
    FALSE 743           IDH3A        <NA>          hi
    FALSE 744            IDI1        <NA>          hi
    FALSE 745          IFNGR1        <NA>          hi
    FALSE 746           IFRD1        <NA>          hi
    FALSE 747           IFT88        <NA>          hi
    FALSE 748           IGSF5        <NA>          hi
    FALSE 749        IL1RAPL1        <NA>          hi
    FALSE 750           IMPA1        <NA>          hi
    FALSE 751          IMPACT        <NA>          hi
    FALSE 752          IMPAD1        <NA>          hi
    FALSE 753          INCENP        <NA>          hi
    FALSE 754          INPP5A        <NA>          hi
    FALSE 755            IPO5        <NA>          hi
    FALSE 756            IPO7        <NA>          hi
    FALSE 757            ISPD        <NA>          hi
    FALSE 758           ITFG1        <NA>          hi
    FALSE 759        ITGB1BP3        <NA>          hi
    FALSE 760           ITM2B        <NA>          hi
    FALSE 761            KARS        <NA>          hi
    FALSE 762          KATNA1        <NA>          hi
    FALSE 763          KBTBD3        <NA>          hi
    FALSE 764           KCNB2        <NA>          hi
    FALSE 765           KCNJ5        <NA>          hi
    FALSE 766           KCNK1        <NA>          hi
    FALSE 767          KCNMA1        <NA>          hi
    FALSE 768           KCNQ3        <NA>          hi
    FALSE 769          KCTD21        <NA>          hi
    FALSE 770          KDELC2        <NA>          hi
    FALSE 771        KIAA0556        <NA>          hi
    FALSE 772        KIAA1107        <NA>          hi
    FALSE 773        KIAA1524        <NA>          hi
    FALSE 774        KIAA1715        <NA>          hi
    FALSE 775           KIF11        <NA>          hi
    FALSE 776           KIF15        <NA>          hi
    FALSE 777          KIF16B        <NA>          hi
    FALSE 778          KIF18A        <NA>          hi
    FALSE 779          KIF20A        <NA>          hi
    FALSE 780           KIF4A        <NA>          hi
    FALSE 781          KLHDC3        <NA>          hi
    FALSE 782           KLHL1        <NA>          hi
    FALSE 783          KLHL13        <NA>          hi
    FALSE 784           KNDC1        <NA>          hi
    FALSE 785          KNSTRN        <NA>          hi
    FALSE 786           KNTC1        <NA>          hi
    FALSE 787           KPNA2        <NA>          hi
    FALSE 788            KYNU        <NA>          hi
    FALSE 789           LAMP2        <NA>          hi
    FALSE 790          LANCL3        <NA>          hi
    FALSE 791         LAPTM4B        <NA>          hi
    FALSE 792             LBH        <NA>          hi
    FALSE 793          LCLAT1        <NA>          hi
    FALSE 794            LDHA        <NA>          hi
    FALSE 795        LEPROTL1        <NA>          hi
    FALSE 796           LETM1        <NA>          hi
    FALSE 797            LGMN        <NA>          hi
    FALSE 798           LIN52        <NA>          hi
    FALSE 799           LMBR1        <NA>          hi
    FALSE 800    LOC100859104        <NA>          hi
    FALSE 801    LOC100859442        <NA>          hi
    FALSE 802    LOC100859557        <NA>          hi
    FALSE 803    LOC101748559        <NA>          hi
    FALSE 804    LOC101749175        <NA>          hi
    FALSE 805    LOC101749287        <NA>          hi
    FALSE 806    LOC101749531        <NA>          hi
    FALSE 807    LOC101749540        <NA>          hi
    FALSE 808    LOC101749834        <NA>          hi
    FALSE 809    LOC101750288        <NA>          hi
    FALSE 810    LOC101750767        <NA>          hi
    FALSE 811    LOC101751638        <NA>          hi
    FALSE 812    LOC101751749        <NA>          hi
    FALSE 813    LOC101751823        <NA>          hi
    FALSE 814    LOC107049088        <NA>          hi
    FALSE 815    LOC107049099        <NA>          hi
    FALSE 816    LOC107050953        <NA>          hi
    FALSE 817    LOC107051182        <NA>          hi
    FALSE 818    LOC107053497        <NA>          hi
    FALSE 819    LOC107054051        <NA>          hi
    FALSE 820    LOC107054305        <NA>          hi
    FALSE 821    LOC107055650        <NA>          hi
    FALSE 822    LOC107055715        <NA>          hi
    FALSE 823    LOC107056124        <NA>          hi
    FALSE 824    LOC107057175        <NA>          hi
    FALSE 825       LOC395647        <NA>          hi
    FALSE 826       LOC407092        <NA>          hi
    FALSE 827       LOC414835        <NA>          hi
    FALSE 828       LOC415641        <NA>          hi
    FALSE 829       LOC417013        <NA>          hi
    FALSE 830       LOC417372        <NA>          hi
    FALSE 831       LOC418811        <NA>          hi
    FALSE 832       LOC420860        <NA>          hi
    FALSE 833       LOC421856        <NA>          hi
    FALSE 834       LOC422151        <NA>          hi
    FALSE 835       LOC422224        <NA>          hi
    FALSE 836       LOC422249        <NA>          hi
    FALSE 837       LOC423793        <NA>          hi
    FALSE 838       LOC423967        <NA>          hi
    FALSE 839       LOC424109        <NA>          hi
    FALSE 840       LOC424201        <NA>          hi
    FALSE 841       LOC424401        <NA>          hi
    FALSE 842       LOC424998        <NA>          hi
    FALSE 843       LOC427665        <NA>          hi
    FALSE 844       LOC769174        <NA>          hi
    FALSE 845       LOC771758        <NA>          hi
    FALSE 846       LOC771811        <NA>          hi
    FALSE 847           LPPR5        <NA>          hi
    FALSE 848            LRP5        <NA>          hi
    FALSE 849          LRRC59        <NA>          hi
    FALSE 850           LRRC7        <NA>          hi
    FALSE 851          LRRTM3        <NA>          hi
    FALSE 852          LRRTM4        <NA>          hi
    FALSE 853            LSM1        <NA>          hi
    FALSE 854           LTA4H        <NA>          hi
    FALSE 855            LTN1        <NA>          hi
    FALSE 856          LYPLA1        <NA>          hi
    FALSE 857         LYPLAL1        <NA>          hi
    FALSE 858           LYRM4        <NA>          hi
    FALSE 859          LYSMD2        <NA>          hi
    FALSE 860          LYSMD4        <NA>          hi
    FALSE 861         MAB21L3        <NA>          hi
    FALSE 862          MAD2L1        <NA>          hi
    FALSE 863        MAD2L1BP        <NA>          hi
    FALSE 864            MAEA        <NA>          hi
    FALSE 865          MAN1A1        <NA>          hi
    FALSE 866            MANF        <NA>          hi
    FALSE 867          MAP2K1        <NA>          hi
    FALSE 868          MAP2K3        <NA>          hi
    FALSE 869          MAP3K7        <NA>          hi
    FALSE 870          MAP7D2        <NA>          hi
    FALSE 871           MAPK8        <NA>          hi
    FALSE 872           MARC2        <NA>          hi
    FALSE 873           MASTL        <NA>          hi
    FALSE 874          MB21D2        <NA>          hi
    FALSE 875            MBD2        <NA>          hi
    FALSE 876           MCFD2        <NA>          hi
    FALSE 877           MCM10        <NA>          hi
    FALSE 878            MCM2        <NA>          hi
    FALSE 879            MCM5        <NA>          hi
    FALSE 880            MCM6        <NA>          hi
    FALSE 881            MCM8        <NA>          hi
    FALSE 882             MCU        <NA>          hi
    FALSE 883            MDH1        <NA>          hi
    FALSE 884           MED22        <NA>          hi
    FALSE 885           MED28        <NA>          hi
    FALSE 886           MED30        <NA>          hi
    FALSE 887           MEIS1        <NA>          hi
    FALSE 888          METAP1        <NA>          hi
    FALSE 889         METAP1D        <NA>          hi
    FALSE 890        METTL21A        <NA>          hi
    FALSE 891         METTL22        <NA>          hi
    FALSE 892          MFSD11        <NA>          hi
    FALSE 893          MFSD2A        <NA>          hi
    FALSE 894            MIA3        <NA>          hi
    FALSE 895            MIB1        <NA>          hi
    FALSE 896          MIPOL1        <NA>          hi
    FALSE 897           MIS12        <NA>          hi
    FALSE 898           MITD1        <NA>          hi
    FALSE 899           MKI67        <NA>          hi
    FALSE 900            MLEC        <NA>          hi
    FALSE 901             MLN        <NA>          hi
    FALSE 902           MMP16        <NA>          hi
    FALSE 903          MMS22L        <NA>          hi
    FALSE 904            MND1        <NA>          hi
    FALSE 905          MOSPD1        <NA>          hi
    FALSE 906       MPHOSPH10        <NA>          hi
    FALSE 907          MRPL15        <NA>          hi
    FALSE 908          MRPL47        <NA>          hi
    FALSE 909          MRPS22        <NA>          hi
    FALSE 910           MRPS9        <NA>          hi
    FALSE 911         MSANTD1        <NA>          hi
    FALSE 912         MSANTD4        <NA>          hi
    FALSE 913            MSH2        <NA>          hi
    FALSE 914           MSMO1        <NA>          hi
    FALSE 915            MTBP        <NA>          hi
    FALSE 916           MTFR1        <NA>          hi
    FALSE 917            MTG2        <NA>          hi
    FALSE 918           MTIF2        <NA>          hi
    FALSE 919            MTPN        <NA>          hi
    FALSE 920             MTR        <NA>          hi
    FALSE 921            MTRR        <NA>          hi
    FALSE 922            MXD3        <NA>          hi
    FALSE 923         MYBBP1A        <NA>          hi
    FALSE 924           MYBL1        <NA>          hi
    FALSE 925          MYCBP2        <NA>          hi
    FALSE 926           MYEF2        <NA>          hi
    FALSE 927            MYOT        <NA>          hi
    FALSE 928           NAA25        <NA>          hi
    FALSE 929           NABP1        <NA>          hi
    FALSE 930            NAPG        <NA>          hi
    FALSE 931            NASP        <NA>          hi
    FALSE 932           NAT10        <NA>          hi
    FALSE 933            NBAS        <NA>          hi
    FALSE 934            NBR1        <NA>          hi
    FALSE 935           NCAM2        <NA>          hi
    FALSE 936          NCAPD3        <NA>          hi
    FALSE 937           NCAPG        <NA>          hi
    FALSE 938           NCAPH        <NA>          hi
    FALSE 939          NCKAP1        <NA>          hi
    FALSE 940             NCL        <NA>          hi
    FALSE 941            NCLN        <NA>          hi
    FALSE 942           NDC80        <NA>          hi
    FALSE 943           NDNL2        <NA>          hi
    FALSE 944          NDUFB5        <NA>          hi
    FALSE 945          NDUFB9        <NA>          hi
    FALSE 946          NDUFC2        <NA>          hi
    FALSE 947          NDUFS1        <NA>          hi
    FALSE 948          NDUFV2        <NA>          hi
    FALSE 949          NECAP1        <NA>          hi
    FALSE 950            NEK2        <NA>          hi
    FALSE 951            NEK6        <NA>          hi
    FALSE 952           NETO1        <NA>          hi
    FALSE 953          NEURL1        <NA>          hi
    FALSE 954           NGLY1        <NA>          hi
    FALSE 955          NHEDC2        <NA>          hi
    FALSE 956            NIP7        <NA>          hi
    FALSE 957          NIPAL3        <NA>          hi
    FALSE 958           NISCH        <NA>          hi
    FALSE 959          NKAIN3        <NA>          hi
    FALSE 960            NMT2        <NA>          hi
    FALSE 961            NNF1        <NA>          hi
    FALSE 962             NNT        <NA>          hi
    FALSE 963           NOL11        <NA>          hi
    FALSE 964            NOM1        <NA>          hi
    FALSE 965           NOP14        <NA>          hi
    FALSE 966           NOP16        <NA>          hi
    FALSE 967           NOP58        <NA>          hi
    FALSE 968            NPM1        <NA>          hi
    FALSE 969            NPM3        <NA>          hi
    FALSE 970           NSUN2        <NA>          hi
    FALSE 971          NT5C1B        <NA>          hi
    FALSE 972           NT5C2        <NA>          hi
    FALSE 973          NT5C3A        <NA>          hi
    FALSE 974            NTN1        <NA>          hi
    FALSE 975           NTN4L        <NA>          hi
    FALSE 976           NUCB2        <NA>          hi
    FALSE 977          NUDCD2        <NA>          hi
    FALSE 978           NUDT6        <NA>          hi
    FALSE 979            NUF2        <NA>          hi
    FALSE 980           NUP58        <NA>          hi
    FALSE 981            NUS1        <NA>          hi
    FALSE 982             NVL        <NA>          hi
    FALSE 983            OAZ1        <NA>          hi
    FALSE 984          OCIAD1        <NA>          hi
    FALSE 985           OFCC1        <NA>          hi
    FALSE 986            OIT3        <NA>          hi
    FALSE 987            OLA1        <NA>          hi
    FALSE 988           OLFM1        <NA>          hi
    FALSE 989           OLFM2        <NA>          hi
    FALSE 990             OMP        <NA>          hi
    FALSE 991            ORC1        <NA>          hi
    FALSE 992         OSBPL11        <NA>          hi
    FALSE 993           OSER1        <NA>          hi
    FALSE 994            OSTC        <NA>          hi
    FALSE 995          OTUD6B        <NA>          hi
    FALSE 996          OTUD7A        <NA>          hi
    FALSE 997            OVST        <NA>          hi
    FALSE 998           OVSTL        <NA>          hi
    FALSE 999            OXR1        <NA>          hi
    FALSE 1000           P4HB        <NA>          hi
    FALSE 1001        PAK1IP1        <NA>          hi
    FALSE 1002            PAM        <NA>          hi
    FALSE 1003         PAPSS1        <NA>          hi
    FALSE 1004          PARK7        <NA>          hi
    FALSE 1005          PARVA        <NA>          hi
    FALSE 1006           PASK        <NA>          hi
    FALSE 1007            PBK        <NA>          hi
    FALSE 1008         PCDH18        <NA>          hi
    FALSE 1009          PCDH7        <NA>          hi
    FALSE 1010          PCDH9        <NA>          hi
    FALSE 1011        PCDHA12        <NA>          hi
    FALSE 1012         PCDHA3        <NA>          hi
    FALSE 1013          PCGF6        <NA>          hi
    FALSE 1014           PCLO        <NA>          hi
    FALSE 1015           PCM1        <NA>          hi
    FALSE 1016           PCNA        <NA>          hi
    FALSE 1017           PCP4        <NA>          hi
    FALSE 1018          PCSK1        <NA>          hi
    FALSE 1019          PCSK6        <NA>          hi
    FALSE 1020         PCYT1B        <NA>          hi
    FALSE 1021          PDCD4        <NA>          hi
    FALSE 1022         PDE10A        <NA>          hi
    FALSE 1023         PDE11A        <NA>          hi
    FALSE 1024          PDE1C        <NA>          hi
    FALSE 1025         PDGFRL        <NA>          hi
    FALSE 1026          PDIA3        <NA>          hi
    FALSE 1027          PDIA4        <NA>          hi
    FALSE 1028         PDXDC1        <NA>          hi
    FALSE 1029          PEX10        <NA>          hi
    FALSE 1030           PEX2        <NA>          hi
    FALSE 1031          PGAM1        <NA>          hi
    FALSE 1032          PGAP1        <NA>          hi
    FALSE 1033           PGM2        <NA>          hi
    FALSE 1034          PHF14        <NA>          hi
    FALSE 1035          PHGDH        <NA>          hi
    FALSE 1036           PHIP        <NA>          hi
    FALSE 1037         PI4K2B        <NA>          hi
    FALSE 1038          PICK1        <NA>          hi
    FALSE 1039           PIGA        <NA>          hi
    FALSE 1040           PIGB        <NA>          hi
    FALSE 1041           PIGK        <NA>          hi
    FALSE 1042           PIGM        <NA>          hi
    FALSE 1043         PIK3CG        <NA>          hi
    FALSE 1044         PITRM1        <NA>          hi
    FALSE 1045         PLA2G7        <NA>          hi
    FALSE 1046        PLEKHG1        <NA>          hi
    FALSE 1047           PLK1        <NA>          hi
    FALSE 1048           PLP1        <NA>          hi
    FALSE 1049           PLS1        <NA>          hi
    FALSE 1050           PLS3        <NA>          hi
    FALSE 1051            PM5        <NA>          hi
    FALSE 1052           PMM2        <NA>          hi
    FALSE 1053           PMS2        <NA>          hi
    FALSE 1054         PNPLA8        <NA>          hi
    FALSE 1055          PNRC2        <NA>          hi
    FALSE 1056         POFUT2        <NA>          hi
    FALSE 1057         POLR3H        <NA>          hi
    FALSE 1058         POU3F2        <NA>          hi
    FALSE 1059         POU3F3        <NA>          hi
    FALSE 1060           PPA1        <NA>          hi
    FALSE 1061           PPIB        <NA>          hi
    FALSE 1062         PPP6R3        <NA>          hi
    FALSE 1063           PRC1        <NA>          hi
    FALSE 1064          PRDX3        <NA>          hi
    FALSE 1065          PRDX4        <NA>          hi
    FALSE 1066        PRELID1        <NA>          hi
    FALSE 1067        PRIMPOL        <NA>          hi
    FALSE 1068        PRKAR2A        <NA>          hi
    FALSE 1069          PRKDC        <NA>          hi
    FALSE 1070            PRL        <NA>          hi
    FALSE 1071          PRMT3        <NA>          hi
    FALSE 1072          PROS1        <NA>          hi
    FALSE 1073          PRPS2        <NA>          hi
    FALSE 1074         PRSS12        <NA>          hi
    FALSE 1075         PRSS23        <NA>          hi
    FALSE 1076          PSMA3        <NA>          hi
    FALSE 1077          PSMB3        <NA>          hi
    FALSE 1078          PSMC1        <NA>          hi
    FALSE 1079         PSMD12        <NA>          hi
    FALSE 1080          PSMD6        <NA>          hi
    FALSE 1081          PSMD7        <NA>          hi
    FALSE 1082          PSMG2        <NA>          hi
    FALSE 1083           PSPH        <NA>          hi
    FALSE 1084          PTGFR        <NA>          hi
    FALSE 1085          PTRH2        <NA>          hi
    FALSE 1086          PUS10        <NA>          hi
    FALSE 1087          PYCR2        <NA>          hi
    FALSE 1088           QDPR        <NA>          hi
    FALSE 1089          RAB1A        <NA>          hi
    FALSE 1090          RAB26        <NA>          hi
    FALSE 1091          RAB28        <NA>          hi
    FALSE 1092          RAB2A        <NA>          hi
    FALSE 1093          RAB30        <NA>          hi
    FALSE 1094       RAB3GAP2        <NA>          hi
    FALSE 1095       RABGAP1L        <NA>          hi
    FALSE 1096          RABL3        <NA>          hi
    FALSE 1097        RACGAP1        <NA>          hi
    FALSE 1098          RAD21        <NA>          hi
    FALSE 1099         RAD54B        <NA>          hi
    FALSE 1100          RADIL        <NA>          hi
    FALSE 1101         RANBP1        <NA>          hi
    FALSE 1102         RANBP2        <NA>          hi
    FALSE 1103        RANGAP1        <NA>          hi
    FALSE 1104          RAP1A        <NA>          hi
    FALSE 1105          RASA2        <NA>          hi
    FALSE 1106         RASSF1        <NA>          hi
    FALSE 1107         RASSF3        <NA>          hi
    FALSE 1108         RASSF6        <NA>          hi
    FALSE 1109            RB1        <NA>          hi
    FALSE 1110         RB1CC1        <NA>          hi
    FALSE 1111          RBBP8        <NA>          hi
    FALSE 1112          RBM41        <NA>          hi
    FALSE 1113           RCN1        <NA>          hi
    FALSE 1114          REEP5        <NA>          hi
    FALSE 1115           RER1        <NA>          hi
    FALSE 1116           REST        <NA>          hi
    FALSE 1117            RET        <NA>          hi
    FALSE 1118           RFC3        <NA>          hi
    FALSE 1119           RFC5        <NA>          hi
    FALSE 1120           RGCC        <NA>          hi
    FALSE 1121           RGS7        <NA>          hi
    FALSE 1122           RHOA        <NA>          hi
    FALSE 1123          RINT1        <NA>          hi
    FALSE 1124          RIOK1        <NA>          hi
    FALSE 1125         RNF182        <NA>          hi
    FALSE 1126           RNH1        <NA>          hi
    FALSE 1127 RP11-1035H13.3        <NA>          hi
    FALSE 1128   RP11-156P1.2        <NA>          hi
    FALSE 1129  RP11-196G11.1        <NA>          hi
    FALSE 1130  RP11-407N17.3        <NA>          hi
    FALSE 1131           RPF2        <NA>          hi
    FALSE 1132           RPGR        <NA>          hi
    FALSE 1133          RPL35        <NA>          hi
    FALSE 1134           RPL7        <NA>          hi
    FALSE 1135           RPN1        <NA>          hi
    FALSE 1136          RPP30        <NA>          hi
    FALSE 1137          RPP40        <NA>          hi
    FALSE 1138          RPS21        <NA>          hi
    FALSE 1139        RPS6KA6        <NA>          hi
    FALSE 1140           RRM2        <NA>          hi
    FALSE 1141          RRP15        <NA>          hi
    FALSE 1142          RRP1B        <NA>          hi
    FALSE 1143         RSL1D1        <NA>          hi
    FALSE 1144           RTCA        <NA>          hi
    FALSE 1145           RTN4        <NA>          hi
    FALSE 1146        RTN4IP1        <NA>          hi
    FALSE 1147           RTTN        <NA>          hi
    FALSE 1148          RUFY1        <NA>          hi
    FALSE 1149        RUNX1T1        <NA>          hi
    FALSE 1150         SAMM50        <NA>          hi
    FALSE 1151          SCFD1        <NA>          hi
    FALSE 1152           SCG3        <NA>          hi
    FALSE 1153          SCN5A        <NA>          hi
    FALSE 1154          SCRN3        <NA>          hi
    FALSE 1155          SCYL2        <NA>          hi
    FALSE 1156         SDF2L1        <NA>          hi
    FALSE 1157           SDF4        <NA>          hi
    FALSE 1158         SEC11C        <NA>          hi
    FALSE 1159          SEC13        <NA>          hi
    FALSE 1160         SEC22C        <NA>          hi
    FALSE 1161         SEC23B        <NA>          hi
    FALSE 1162        SEC23IP        <NA>          hi
    FALSE 1163         SEC24D        <NA>          hi
    FALSE 1164        SEC61A2        <NA>          hi
    FALSE 1165         SEC61B        <NA>          hi
    FALSE 1166         SEC61G        <NA>          hi
    FALSE 1167           SELK        <NA>          hi
    FALSE 1168          SEP15        <NA>          hi
    FALSE 1169          SEPP1        <NA>          hi
    FALSE 1170         SEPP1L        <NA>          hi
    FALSE 1171        SERINC1        <NA>          hi
    FALSE 1172          SERP1        <NA>          hi
    FALSE 1173       SERPINB6        <NA>          hi
    FALSE 1174          SETD3        <NA>          hi
    FALSE 1175          SETD8        <NA>          hi
    FALSE 1176           SFPQ        <NA>          hi
    FALSE 1177          SFRP5        <NA>          hi
    FALSE 1178          SFXN4        <NA>          hi
    FALSE 1179           SGCZ        <NA>          hi
    FALSE 1180          SGMS1        <NA>          hi
    FALSE 1181          SGOL1        <NA>          hi
    FALSE 1182          SGSM3        <NA>          hi
    FALSE 1183         SH3TC1        <NA>          hi
    FALSE 1184           SHC2        <NA>          hi
    FALSE 1185         SHCBP1        <NA>          hi
    FALSE 1186          SHOX2        <NA>          hi
    FALSE 1187          SHPRH        <NA>          hi
    FALSE 1188          SIRT5        <NA>          hi
    FALSE 1189           SIX3        <NA>          hi
    FALSE 1190           SKA1        <NA>          hi
    FALSE 1191           SKA3        <NA>          hi
    FALSE 1192           SKIL        <NA>          hi
    FALSE 1193           SLBP        <NA>          hi
    FALSE 1194        SLC16A1        <NA>          hi
    FALSE 1195        SLC17A5        <NA>          hi
    FALSE 1196        SLC17A9        <NA>          hi
    FALSE 1197         SLC1A2        <NA>          hi
    FALSE 1198         SLC1A4        <NA>          hi
    FALSE 1199       SLC22A15        <NA>          hi
    FALSE 1200        SLC23A2        <NA>          hi
    FALSE 1201       SLC25A12        <NA>          hi
    FALSE 1202       SLC25A21        <NA>          hi
    FALSE 1203       SLC25A22        <NA>          hi
    FALSE 1204        SLC25A3        <NA>          hi
    FALSE 1205       SLC25A34        <NA>          hi
    FALSE 1206        SLC30A7        <NA>          hi
    FALSE 1207        SLC33A1        <NA>          hi
    FALSE 1208        SLC35B1        <NA>          hi
    FALSE 1209        SLC35B3        <NA>          hi
    FALSE 1210        SLC35F2        <NA>          hi
    FALSE 1211        SLC36A4        <NA>          hi
    FALSE 1212        SLC37A3        <NA>          hi
    FALSE 1213        SLC38A2        <NA>          hi
    FALSE 1214       SLC39A13        <NA>          hi
    FALSE 1215        SLC41A2        <NA>          hi
    FALSE 1216        SLC5A10        <NA>          hi
    FALSE 1217         SLC7A3        <NA>          hi
    FALSE 1218         SLC7A5        <NA>          hi
    FALSE 1219         SLC7A6        <NA>          hi
    FALSE 1220         SLC9A6        <NA>          hi
    FALSE 1221        SLCO4A1        <NA>          hi
    FALSE 1222        SLITRK6        <NA>          hi
    FALSE 1223            SLK        <NA>          hi
    FALSE 1224          SLMO2        <NA>          hi
    FALSE 1225         SLX4IP        <NA>          hi
    FALSE 1226       SMARCAL1        <NA>          hi
    FALSE 1227          SMCO4        <NA>          hi
    FALSE 1228           SMG8        <NA>          hi
    FALSE 1229         SMIM18        <NA>          hi
    FALSE 1230         SMNDC1        <NA>          hi
    FALSE 1231          SMPD3        <NA>          hi
    FALSE 1232         SMTNL1        <NA>          hi
    FALSE 1233         SNRPA1        <NA>          hi
    FALSE 1234          SNRPF        <NA>          hi
    FALSE 1235          SNX14        <NA>          hi
    FALSE 1236           SNX6        <NA>          hi
    FALSE 1237           SOD1        <NA>          hi
    FALSE 1238          SOGA3        <NA>          hi
    FALSE 1239        SOSTDC1        <NA>          hi
    FALSE 1240           SOX4        <NA>          hi
    FALSE 1241          SPAM1        <NA>          hi
    FALSE 1242         SPATA5        <NA>          hi
    FALSE 1243          SPC24        <NA>          hi
    FALSE 1244          SPC25        <NA>          hi
    FALSE 1245          SPCS1        <NA>          hi
    FALSE 1246          SPCS2        <NA>          hi
    FALSE 1247          SPCS3        <NA>          hi
    FALSE 1248         SPOCK3        <NA>          hi
    FALSE 1249          SPOPL        <NA>          hi
    FALSE 1250         SPRED1        <NA>          hi
    FALSE 1251           SQLE        <NA>          hi
    FALSE 1252            SRM        <NA>          hi
    FALSE 1253          SRP14        <NA>          hi
    FALSE 1254          SRP54        <NA>          hi
    FALSE 1255          SRP72        <NA>          hi
    FALSE 1256           SRPR        <NA>          hi
    FALSE 1257          SRPRB        <NA>          hi
    FALSE 1258         SS18L2        <NA>          hi
    FALSE 1259            SS2        <NA>          hi
    FALSE 1260            SSB        <NA>          hi
    FALSE 1261           SSPO        <NA>          hi
    FALSE 1262           SSR1        <NA>          hi
    FALSE 1263           SSR2        <NA>          hi
    FALSE 1264           SSR3        <NA>          hi
    FALSE 1265           SSR4        <NA>          hi
    FALSE 1266           ST18        <NA>          hi
    FALSE 1267           STAC        <NA>          hi
    FALSE 1268       STAMBPL1        <NA>          hi
    FALSE 1269           STC1        <NA>          hi
    FALSE 1270         STEAP1        <NA>          hi
    FALSE 1271         STEAP2        <NA>          hi
    FALSE 1272          STK16        <NA>          hi
    FALSE 1273           STK3        <NA>          hi
    FALSE 1274           STK4        <NA>          hi
    FALSE 1275          STMN1        <NA>          hi
    FALSE 1276          STRAP        <NA>          hi
    FALSE 1277          STX12        <NA>          hi
    FALSE 1278         STXBP5        <NA>          hi
    FALSE 1279           SUB1        <NA>          hi
    FALSE 1280         SUB1L1        <NA>          hi
    FALSE 1281           SUCO        <NA>          hi
    FALSE 1282          SURF4        <NA>          hi
    FALSE 1283          SYNPR        <NA>          hi
    FALSE 1284          SYT13        <NA>          hi
    FALSE 1285          TACR3        <NA>          hi
    FALSE 1286         TALDO1        <NA>          hi
    FALSE 1287        TAX1BP1        <NA>          hi
    FALSE 1288        TBC1D15        <NA>          hi
    FALSE 1289        TBC1D32        <NA>          hi
    FALSE 1290         TBC1D9        <NA>          hi
    FALSE 1291           TBCE        <NA>          hi
    FALSE 1292         TBXAS1        <NA>          hi
    FALSE 1293        TCERG1L        <NA>          hi
    FALSE 1294           TCP1        <NA>          hi
    FALSE 1295          TEAD1        <NA>          hi
    FALSE 1296          TENM2        <NA>          hi
    FALSE 1297           TESC        <NA>          hi
    FALSE 1298         TFAP2B        <NA>          hi
    FALSE 1299          TFDP1        <NA>          hi
    FALSE 1300          TFPI2        <NA>          hi
    FALSE 1301           TGDS        <NA>          hi
    FALSE 1302           TGS1        <NA>          hi
    FALSE 1303          THADA        <NA>          hi
    FALSE 1304          THAP5        <NA>          hi
    FALSE 1305          THOC1        <NA>          hi
    FALSE 1306          THOC2        <NA>          hi
    FALSE 1307         TIMM44        <NA>          hi
    FALSE 1308            TK1        <NA>          hi
    FALSE 1309          TM2D3        <NA>          hi
    FALSE 1310         TM9SF3        <NA>          hi
    FALSE 1311           TMC5        <NA>          hi
    FALSE 1312         TMED10        <NA>          hi
    FALSE 1313       TMEM120B        <NA>          hi
    FALSE 1314        TMEM159        <NA>          hi
    FALSE 1315        TMEM181        <NA>          hi
    FALSE 1316       TMEM185A        <NA>          hi
    FALSE 1317        TMEM186        <NA>          hi
    FALSE 1318        TMEM196        <NA>          hi
    FALSE 1319        TMEM240        <NA>          hi
    FALSE 1320        TMEM241        <NA>          hi
    FALSE 1321        TMEM39A        <NA>          hi
    FALSE 1322        TMEM50A        <NA>          hi
    FALSE 1323        TMEM52B        <NA>          hi
    FALSE 1324         TMEM57        <NA>          hi
    FALSE 1325         TMEM59        <NA>          hi
    FALSE 1326           TMF1        <NA>          hi
    FALSE 1327          TMOD3        <NA>          hi
    FALSE 1328           TMPO        <NA>          hi
    FALSE 1329       TMPRSS13        <NA>          hi
    FALSE 1330         TMSB10        <NA>          hi
    FALSE 1331          TMTC3        <NA>          hi
    FALSE 1332           TMX3        <NA>          hi
    FALSE 1333          TNKS2        <NA>          hi
    FALSE 1334         TOMM20        <NA>          hi
    FALSE 1335          TOP2A        <NA>          hi
    FALSE 1336           TPP2        <NA>          hi
    FALSE 1337          TPRKB        <NA>          hi
    FALSE 1338           TPT1        <NA>          hi
    FALSE 1339          TRA2B        <NA>          hi
    FALSE 1340          TRABD        <NA>          hi
    FALSE 1341          TRAM1        <NA>          hi
    FALSE 1342         TRIM35        <NA>          hi
    FALSE 1343         TRIP11        <NA>          hi
    FALSE 1344         TRMT11        <NA>          hi
    FALSE 1345         TRMT1L        <NA>          hi
    FALSE 1346          TRPC1        <NA>          hi
    FALSE 1347          TSHZ1        <NA>          hi
    FALSE 1348            TSN        <NA>          hi
    FALSE 1349        TSPAN13        <NA>          hi
    FALSE 1350           TSR3        <NA>          hi
    FALSE 1351          TSTA3        <NA>          hi
    FALSE 1352          TTC13        <NA>          hi
    FALSE 1353          TTC17        <NA>          hi
    FALSE 1354         TTC21B        <NA>          hi
    FALSE 1355          TTC27        <NA>          hi
    FALSE 1356           TTC3        <NA>          hi
    FALSE 1357           TTF1        <NA>          hi
    FALSE 1358            TTK        <NA>          hi
    FALSE 1359            TTL        <NA>          hi
    FALSE 1360         TTLL11        <NA>          hi
    FALSE 1361            TUB        <NA>          hi
    FALSE 1362          TUBD1        <NA>          hi
    FALSE 1363          TUBG1        <NA>          hi
    FALSE 1364        TUBGCP5        <NA>          hi
    FALSE 1365          TUSC3        <NA>          hi
    FALSE 1366        TXNDC11        <NA>          hi
    FALSE 1367         TXNDC5        <NA>          hi
    FALSE 1368          TXNL1        <NA>          hi
    FALSE 1369           TYMS        <NA>          hi
    FALSE 1370           TYW5        <NA>          hi
    FALSE 1371          U2AF1        <NA>          hi
    FALSE 1372           UBA5        <NA>          hi
    FALSE 1373          UBE2A        <NA>          hi
    FALSE 1374          UBE2C        <NA>          hi
    FALSE 1375         UBE2V1        <NA>          hi
    FALSE 1376         UBLCP1        <NA>          hi
    FALSE 1377           UBR3        <NA>          hi
    FALSE 1378           UBR7        <NA>          hi
    FALSE 1379         UBXN2A        <NA>          hi
    FALSE 1380          UCHL5        <NA>          hi
    FALSE 1381           UCK2        <NA>          hi
    FALSE 1382           UFC1        <NA>          hi
    FALSE 1383           UFL1        <NA>          hi
    FALSE 1384           UFM1        <NA>          hi
    FALSE 1385          UGGT1        <NA>          hi
    FALSE 1386          UGGT2        <NA>          hi
    FALSE 1387      UHRF1BP1L        <NA>          hi
    FALSE 1388         UMODL1        <NA>          hi
    FALSE 1389          UNC5C        <NA>          hi
    FALSE 1390          UNC79        <NA>          hi
    FALSE 1391           USP1        <NA>          hi
    FALSE 1392          USP13        <NA>          hi
    FALSE 1393          USP16        <NA>          hi
    FALSE 1394          USP47        <NA>          hi
    FALSE 1395          UTP23        <NA>          hi
    FALSE 1396          VDAC1        <NA>          hi
    FALSE 1397          VEPH1        <NA>          hi
    FALSE 1398           VIMP        <NA>          hi
    FALSE 1399          VIPR1        <NA>          hi
    FALSE 1400          VMA21        <NA>          hi
    FALSE 1401         VPS13A        <NA>          hi
    FALSE 1402          VPS50        <NA>          hi
    FALSE 1403          WDHD1        <NA>          hi
    FALSE 1404          WDR17        <NA>          hi
    FALSE 1405          WDR35        <NA>          hi
    FALSE 1406           WDR4        <NA>          hi
    FALSE 1407          WDR43        <NA>          hi
    FALSE 1408          WDR75        <NA>          hi
    FALSE 1409          WDR89        <NA>          hi
    FALSE 1410          WIPI1        <NA>          hi
    FALSE 1411          WSCD1        <NA>          hi
    FALSE 1412           XBP1        <NA>          hi
    FALSE 1413           XPO4        <NA>          hi
    FALSE 1414          XRCC5        <NA>          hi
    FALSE 1415          XRCC6        <NA>          hi
    FALSE 1416           XRN1        <NA>          hi
    FALSE 1417           YAF2        <NA>          hi
    FALSE 1418           YBX1        <NA>          hi
    FALSE 1419          YIPF1        <NA>          hi
    FALSE 1420          YIPF4        <NA>          hi
    FALSE 1421          YIPF5        <NA>          hi
    FALSE 1422         YME1L1        <NA>          hi
    FALSE 1423          YWHAH        <NA>          hi
    FALSE 1424          YWHAQ        <NA>          hi
    FALSE 1425          ZBED4        <NA>          hi
    FALSE 1426          ZCRB1        <NA>          hi
    FALSE 1427         ZDHHC3        <NA>          hi
    FALSE 1428         ZNF148        <NA>          hi
    FALSE 1429         ZNF277        <NA>          hi
    FALSE 1430         ZNF644        <NA>          hi
    FALSE 1431        ZNF804B        <NA>          hi
    FALSE 1432         ZWILCH        <NA>          hi

    PRLwgcna <- read_csv("../results/08_PRL_associated.csv") %>%
      dplyr::rename("gene" =  "x") %>%
      mutate(wgcna = "yes")
    PRLwgcna

    FALSE # A tibble: 98 x 2
    FALSE    gene     wgcna
    FALSE    <chr>    <chr>
    FALSE  1 ANKLE1   yes  
    FALSE  2 ARHGAP19 yes  
    FALSE  3 ASPM     yes  
    FALSE  4 AURKA    yes  
    FALSE  5 BIRC5    yes  
    FALSE  6 BRCA1    yes  
    FALSE  7 BUB1     yes  
    FALSE  8 CASC5    yes  
    FALSE  9 CCNA2    yes  
    FALSE 10 CCNB3    yes  
    FALSE # … with 88 more rows

    PRLwgcnaDEGs <- left_join(PRLwgcna, sexPRLdesg) %>%
      select(-direction.x) %>%
      dplyr::rename("PRL" =  "direction.y") 
    PRLwgcnaDEGs %>%
      group_by(wgcna, PRL) %>%
      summarize(n = n())

    FALSE # A tibble: 2 x 3
    FALSE # Groups:   wgcna [1]
    FALSE   wgcna PRL       n
    FALSE   <chr> <fct> <int>
    FALSE 1 yes   hi       91
    FALSE 2 yes   <NA>      7

    head(PRLwgcnaDEGs)

    FALSE # A tibble: 6 x 3
    FALSE   gene     wgcna PRL  
    FALSE   <chr>    <chr> <fct>
    FALSE 1 ANKLE1   yes   <NA> 
    FALSE 2 ARHGAP19 yes   hi   
    FALSE 3 ASPM     yes   hi   
    FALSE 4 AURKA    yes   hi   
    FALSE 5 BIRC5    yes   hi   
    FALSE 6 BRCA1    yes   hi

    # 90 of the 98 genes in the wgcna module are also DEGs, with increase expression in hi PRL

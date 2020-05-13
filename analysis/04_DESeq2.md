    library(tidyverse)

    ## ── Attaching packages ───────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

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

    source("../R/themes.R")
    source("../R/functions.R")

### Hypothesis related data wrangle

-   PRL WGCNA
-   hi v lo PRL
-   extral environment

<!-- -->

    countData <- read.csv("../results/00_counts.csv", header = T, row.names = 1)

    ## sample info
    colData <- read_csv("../metadata/00_colData.csv") %>%
      mutate(tissue = factor(tissue, levels = tissuelevels),
             treatment = factor(treatment, levels = alllevels)) %>%
      mutate(sextissue = as.factor(paste(sex,tissue, sep = "_"))) %>%
      column_to_rownames(var = "X1") %>%
      mutate("samples" = V1) 

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_character(),
    ##   V1 = col_character(),
    ##   bird = col_character(),
    ##   sex = col_character(),
    ##   tissue = col_character(),
    ##   treatment = col_character(),
    ##   group = col_character(),
    ##   study = col_character()
    ## )

    head(colData)

    ##                                       V1    bird  sex       tissue
    ## 1        L.Blu13_male_gonad_control.NYNO L.Blu13 male       gonads
    ## 2 L.Blu13_male_hypothalamus_control.NYNO L.Blu13 male hypothalamus
    ## 3    L.Blu13_male_pituitary_control.NYNO L.Blu13 male    pituitary
    ## 4              L.G107_male_gonad_control  L.G107 male       gonads
    ## 5       L.G107_male_hypothalamus_control  L.G107 male hypothalamus
    ## 6          L.G107_male_pituitary_control  L.G107 male    pituitary
    ##   treatment                     group           study         sextissue
    ## 1   control       male.gonads.control charcterization       male_gonads
    ## 2   control male.hypothalamus.control charcterization male_hypothalamus
    ## 3   control    male.pituitary.control charcterization    male_pituitary
    ## 4   control       male.gonads.control charcterization       male_gonads
    ## 5   control male.hypothalamus.control charcterization male_hypothalamus
    ## 6   control    male.pituitary.control charcterization    male_pituitary
    ##                                  samples
    ## 1        L.Blu13_male_gonad_control.NYNO
    ## 2 L.Blu13_male_hypothalamus_control.NYNO
    ## 3    L.Blu13_male_pituitary_control.NYNO
    ## 4              L.G107_male_gonad_control
    ## 5       L.G107_male_hypothalamus_control
    ## 6          L.G107_male_pituitary_control

    ## calcualte hi lo PRL, add to col Data
    PRLvsd <- read_csv("../results/03_candidatevsd.csv") %>%
      filter(gene == "PRL", tissue == "pituitary") 

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_double(),
    ##   sex = col_character(),
    ##   tissue = col_character(),
    ##   treatment = col_character(),
    ##   gene = col_character(),
    ##   samples = col_character(),
    ##   counts = col_double()
    ## )

    PRLvsd %>%
      group_by(sex) %>%
      summarize(avecounts = mean(counts))

    ## # A tibble: 2 x 2
    ##   sex    avecounts
    ##   <chr>      <dbl>
    ## 1 female      18.9
    ## 2 male        18.6

    PRLvsd <- PRLvsd %>%
      mutate(hiloPRL = ifelse(counts >= 18.5, "hi", "lo")) %>%
      mutate(bird = sapply(strsplit(samples, '\\_'), "[", 1)) %>%
      select(bird, hiloPRL) 
    colData <- left_join(colData, PRLvsd, by = "bird")   
    unique(colData$hiloPRL)

    ## [1] "lo" "hi" NA

    ## add external hypothesis

    colData <- colData %>%
      mutate(external = fct_collapse(treatment,
                                     eggs = c("lay", "inc.d3", "inc.d9", "inc.d17", "prolong"),
                                     chicks = c("hatch", "n5", "n9", "extend" , "m.inc.d8" ),
                                     nest = c("bldg", "m.inc.d3",  "m.inc.d9",  "m.inc.d17",  "m.n2"),
                                     control = c("control"))) %>%
      drop_na()

    ## Warning: Unknown levels in `f`: m.inc.d8

Hypothesis
----------

    createDEGdfhypothesis <- function(variable, up, down, mytissue){
      
      res <- results(dds, contrast = c(variable, up, down), 
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
      myfilename = paste0("../results/DESeq2/hypothesis/", mytissue, partialfilename, "_DEGs.csv")
      
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
                                    design = ~ hiloPRL + external )
      dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
      print(dds)
      print(dim(dds))
      dds <- DESeq(dds, parallel = TRUE) # Differential expression analysis
      
      vsd <- as.data.frame(assay(varianceStabilizingTransformation(dds, blind=FALSE)))
      
      myfilename = paste0("../results/DEseq2/hypothesis/", i, "_vsd.csv")
      write.csv(vsd, myfilename)
      
      createDEGdfhypothesis("hiloPRL", "hi", "lo", i)
      createDEGdfhypothesis("external", "chicks", "eggs", i)
      createDEGdfhypothesis("external", "chicks", "nest", i)
      createDEGdfhypothesis("external", "eggs", "nest", i)

    }

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## class: DESeqDataSet 
    ## dim: 13801 166 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(13801): A2ML1 A2ML2 ... ZYX ZZZ3
    ## rowData names(0):
    ## colnames(166): L.G118_female_gonad_control
    ##   R.G106_female_gonad_control ... y97.x_female_gonad_n9
    ##   y98.g54_female_gonad_m.hatch
    ## colData names(11): V1 bird ... hiloPRL external
    ## [1] 13801   166

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 338 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    ## 'data.frame':    3985 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13801 levels "A2ML1","A2ML2",..: 12307 9389 1370 7234 3094 1345 965 9834 6782 11242 ...
    ##  $ padj     : num  0.02778 0.00735 0.03899 0.02876 0.00624 ...
    ##  $ logpadj  : num  1.56 2.13 1.41 1.54 2.21 ...
    ##  $ lfc      : num  1.97 1.88 1.76 1.59 1.58 ...
    ##  $ sextissue: Factor w/ 1 level "female_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "lo","NS","hi": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL
    ## 'data.frame':    635 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13801 levels "A2ML1","A2ML2",..: 6196 2541 2364 8830 6745 10145 11686 1889 6871 1528 ...
    ##  $ padj     : num  5.06e-07 2.18e-02 5.37e-05 2.03e-02 7.18e-12 ...
    ##  $ logpadj  : num  6.3 1.66 4.27 1.69 11.14 ...
    ##  $ lfc      : num  3.97 3.71 3.38 3.07 3.06 ...
    ##  $ sextissue: Factor w/ 1 level "female_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "eggs","NS","chicks": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL
    ## 'data.frame':    687 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13801 levels "A2ML1","A2ML2",..: 2541 6408 6745 2364 8715 2238 4154 2741 8739 3539 ...
    ##  $ padj     : num  5.31e-04 1.16e-02 6.91e-11 1.07e-03 1.89e-02 ...
    ##  $ logpadj  : num  3.27 1.93 10.16 2.97 1.72 ...
    ##  $ lfc      : num  5.04 3.08 3 2.88 2.79 ...
    ##  $ sextissue: Factor w/ 1 level "female_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "nest","NS","chicks": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL
    ## 'data.frame':    1942 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13801 levels "A2ML1","A2ML2",..: 925 8833 6634 4375 8045 3094 8739 2741 4198 4207 ...
    ##  $ padj     : num  8.05e-02 4.37e-10 6.64e-02 1.34e-06 2.46e-03 ...
    ##  $ logpadj  : num  1.09 9.36 1.18 5.87 2.61 ...
    ##  $ lfc      : num  2.95 2.81 2.73 2.58 2.52 ...
    ##  $ sextissue: Factor w/ 1 level "female_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "nest","NS","eggs": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## class: DESeqDataSet 
    ## dim: 13631 165 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(13631): A2ML1 A2ML2 ... ZYX ZZZ3
    ## rowData names(0):
    ## colnames(165): L.G118_female_hypothalamus_control.NYNO
    ##   R.G106_female_hypothalamus_control ...
    ##   y97.x_female_hypothalamus_n9 y98.g54_female_hypothalamus_m.hatch
    ## colData names(11): V1 bird ... hiloPRL external
    ## [1] 13631   165

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 20 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    ## 'data.frame':    7 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13631 levels "A2ML1","A2ML2",..: 334 10440 5612 13482 6910 1053 1004
    ##  $ padj     : num  0.0921 0.0791 0.0791 0.0921 0.0791 ...
    ##  $ logpadj  : num  1.04 1.1 1.1 1.04 1.1 ...
    ##  $ lfc      : num  0.2458 0.2328 0.1687 0.1585 -0.0846 ...
    ##  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1 1 1 1 1 1 1
    ##  $ direction: Factor w/ 3 levels "lo","NS","hi": 3 3 3 3 1 1 1
    ## NULL
    ## 'data.frame':    9 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13631 levels "A2ML1","A2ML2",..: 37 12122 9284 3985 12786 441 12846 5334 9711
    ##  $ padj     : num  0.02164 0.02704 0.02704 0.00124 0.00391 ...
    ##  $ logpadj  : num  1.66 1.57 1.57 2.91 2.41 ...
    ##  $ lfc      : num  2.976 2 1.538 0.434 0.416 ...
    ##  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1 1 1 1 1 1 1 1 1
    ##  $ direction: Factor w/ 3 levels "eggs","NS","chicks": 3 3 3 3 3 3 3 3 3
    ## NULL
    ## 'data.frame':    1680 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13631 levels "A2ML1","A2ML2",..: 13570 3267 8413 10946 5266 10443 2618 401 12435 10567 ...
    ##  $ padj     : num  0.0123 0.0653 0.0466 0.0497 0.0252 ...
    ##  $ logpadj  : num  1.91 1.19 1.33 1.3 1.6 ...
    ##  $ lfc      : num  2.221 1.411 1.25 0.834 0.77 ...
    ##  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "nest","NS","chicks": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL
    ## 'data.frame':    1030 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13631 levels "A2ML1","A2ML2",..: 2950 6179 3940 4268 3895 10523 4616 1129 1719 8713 ...
    ##  $ padj     : num  0.071486 0.007518 0.000056 0.017826 0.078943 ...
    ##  $ logpadj  : num  1.15 2.12 4.25 1.75 1.1 ...
    ##  $ lfc      : num  0.98 0.949 0.861 0.796 0.722 ...
    ##  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "nest","NS","eggs": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## class: DESeqDataSet 
    ## dim: 13554 165 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(13554): A2ML1 A2ML2 ... ZYX ZZZ3
    ## rowData names(0):
    ## colnames(165): L.G118_female_pituitary_control.NYNO
    ##   R.G106_female_pituitary_control ... y97.x_female_pituitary_n9
    ##   y98.g54_female_pituitary_m.hatch
    ## colData names(11): V1 bird ... hiloPRL external
    ## [1] 13554   165

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 85 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    ## 'data.frame':    5962 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13554 levels "A2ML1","A2ML2",..: 4028 5632 6111 9427 12401 3920 10852 2177 1919 1776 ...
    ##  $ padj     : num  7.52e-16 4.80e-32 4.90e-12 9.64e-06 6.07e-06 ...
    ##  $ logpadj  : num  15.12 31.32 11.31 5.02 5.22 ...
    ##  $ lfc      : num  3.25 3.22 3.14 3.08 3.02 ...
    ##  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "lo","NS","hi": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL
    ## 'data.frame':    418 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13554 levels "A2ML1","A2ML2",..: 9230 10639 13084 9871 4690 9857 6528 6406 12726 11028 ...
    ##  $ padj     : num  0.023647 0.002279 0.000204 0.00204 0.012958 ...
    ##  $ logpadj  : num  1.63 2.64 3.69 2.69 1.89 ...
    ##  $ lfc      : num  2.05 1.85 1.59 1.56 1.52 ...
    ##  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "eggs","NS","chicks": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL
    ## 'data.frame':    2127 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13554 levels "A2ML1","A2ML2",..: 9769 8420 9517 2969 4690 7609 6125 3645 10958 479 ...
    ##  $ padj     : num  9.22e-02 5.52e-02 4.39e-02 1.49e-08 2.11e-04 ...
    ##  $ logpadj  : num  1.04 1.26 1.36 7.83 3.67 ...
    ##  $ lfc      : num  2.12 2.06 1.88 1.83 1.76 ...
    ##  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "nest","NS","chicks": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL
    ## 'data.frame':    3465 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13554 levels "A2ML1","A2ML2",..: 9493 7911 12401 9517 8611 9979 6125 6539 3645 2484 ...
    ##  $ padj     : num  2.13e-10 1.31e-06 1.32e-03 1.83e-04 2.04e-05 ...
    ##  $ logpadj  : num  9.67 5.88 2.88 3.74 4.69 ...
    ##  $ lfc      : num  3.41 2.91 2.58 2.51 2.2 ...
    ##  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "nest","NS","eggs": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## class: DESeqDataSet 
    ## dim: 13825 163 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(13825): A2ML1 A2ML2 ... ZYX ZZZ3
    ## rowData names(0):
    ## colnames(163): L.Blu13_male_gonad_control.NYNO
    ##   L.G107_male_gonad_control ... y95.g131.x_male_gonad_inc.d9
    ##   y98.o50.x_male_gonad_inc.d3
    ## colData names(11): V1 bird ... hiloPRL external
    ## [1] 13825   163

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 391 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    ## 'data.frame':    1439 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13825 levels "A2ML1","A2ML2",..: 2700 1683 3994 1586 12902 4973 11166 11899 3201 12452 ...
    ##  $ padj     : num  0.0193 0.0224 0.0295 0.0176 0.0273 ...
    ##  $ logpadj  : num  1.71 1.65 1.53 1.75 1.56 ...
    ##  $ lfc      : num  1.72 1.56 1.48 1.39 1.38 ...
    ##  $ sextissue: Factor w/ 1 level "male_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "lo","NS","hi": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL
    ## 'data.frame':    156 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13825 levels "A2ML1","A2ML2",..: 12978 11108 9427 6722 1315 3284 12404 2574 7066 1921 ...
    ##  $ padj     : num  0.00407 0.01368 0.06396 0.09328 0.09858 ...
    ##  $ logpadj  : num  2.39 1.86 1.19 1.03 1.01 ...
    ##  $ lfc      : num  1.342 0.855 0.779 0.694 0.574 ...
    ##  $ sextissue: Factor w/ 1 level "male_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "eggs","NS","chicks": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL
    ## 'data.frame':    0 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13825 levels "A2ML1","A2ML2",..: 
    ##  $ padj     : num 
    ##  $ logpadj  : num 
    ##  $ lfc      : num 
    ##  $ sextissue: Factor w/ 1 level "male_gonads": 
    ##  $ direction: Factor w/ 3 levels "nest","NS","chicks": 
    ## NULL
    ## 'data.frame':    40 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13825 levels "A2ML1","A2ML2",..: 959 2150 10740 6533 10118 2368 1950 4972 6068 454 ...
    ##  $ padj     : num  0.0845 0.0437 0.0384 0.0845 0.0557 ...
    ##  $ logpadj  : num  1.07 1.36 1.42 1.07 1.25 ...
    ##  $ lfc      : num  2 1.87 1.53 1.48 1.09 ...
    ##  $ sextissue: Factor w/ 1 level "male_gonads": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "nest","NS","eggs": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## class: DESeqDataSet 
    ## dim: 13595 161 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(13595): A2ML1 A2ML2 ... ZYX ZZZ3
    ## rowData names(0):
    ## colnames(161): L.Blu13_male_hypothalamus_control.NYNO
    ##   L.G107_male_hypothalamus_control ...
    ##   y95.g131.x_male_hypothalamus_inc.d9
    ##   y98.o50.x_male_hypothalamus_inc.d3
    ## colData names(11): V1 bird ... hiloPRL external
    ## [1] 13595   161

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 21 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    ## 'data.frame':    0 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13595 levels "A2ML1","A2ML2",..: 
    ##  $ padj     : num 
    ##  $ logpadj  : num 
    ##  $ lfc      : num 
    ##  $ sextissue: Factor w/ 1 level "male_hypothalamus": 
    ##  $ direction: Factor w/ 3 levels "lo","NS","hi": 
    ## NULL
    ## 'data.frame':    2 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13595 levels "A2ML1","A2ML2",..: 12508 3977
    ##  $ padj     : num  0.0689 0.0139
    ##  $ logpadj  : num  1.16 1.86
    ##  $ lfc      : num  1.342 0.458
    ##  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1 1
    ##  $ direction: Factor w/ 3 levels "eggs","NS","chicks": 3 3
    ## NULL
    ## 'data.frame':    13 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13595 levels "A2ML1","A2ML2",..: 9128 12508 7945 1300 242 5199 10834 1051 10410 560 ...
    ##  $ padj     : num  0.083136 0.000971 0.083136 0.083136 0.091794 ...
    ##  $ logpadj  : num  1.08 3.01 1.08 1.08 1.04 ...
    ##  $ lfc      : num  1.804 1.637 1.441 0.632 0.521 ...
    ##  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "nest","NS","chicks": 3 3 3 3 3 3 1 1 1 1 ...
    ## NULL
    ## 'data.frame':    4 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13595 levels "A2ML1","A2ML2",..: 8404 10661 3977 9526
    ##  $ padj     : num  0.040647 0.040647 0.000702 0.011918
    ##  $ logpadj  : num  1.39 1.39 3.15 1.92
    ##  $ lfc      : num  0.593 0.358 -0.429 -1.384
    ##  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1 1 1 1
    ##  $ direction: Factor w/ 3 levels "nest","NS","eggs": 3 3 1 1
    ## NULL

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## class: DESeqDataSet 
    ## dim: 13541 165 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(13541): A2ML1 A2ML2 ... ZYX ZZZ3
    ## rowData names(0):
    ## colnames(165): L.Blu13_male_pituitary_control.NYNO
    ##   L.G107_male_pituitary_control ...
    ##   y95.g131.x_male_pituitary_inc.d9 y98.o50.x_male_pituitary_inc.d3
    ## colData names(11): V1 bird ... hiloPRL external
    ## [1] 13541   165

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 6 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 6 workers

    ## -- replacing outliers and refitting for 200 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    ## 'data.frame':    4722 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13541 levels "A2ML1","A2ML2",..: 5617 6097 4022 10838 3504 2171 10406 1880 1884 1912 ...
    ##  $ padj     : num  7.56e-61 7.15e-21 1.13e-23 1.70e-34 5.05e-33 ...
    ##  $ logpadj  : num  60.1 20.1 22.9 33.8 32.3 ...
    ##  $ lfc      : num  3.79 3.68 3.67 3.52 3.47 ...
    ##  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "lo","NS","hi": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL
    ## 'data.frame':    90 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13541 levels "A2ML1","A2ML2",..: 7882 12711 200 1789 476 4309 13304 6737 11749 1956 ...
    ##  $ padj     : num  0.0163 0.0185 0.0234 0.0589 0.0332 ...
    ##  $ logpadj  : num  1.79 1.73 1.63 1.23 1.48 ...
    ##  $ lfc      : num  1.285 1 0.776 0.75 0.704 ...
    ##  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "eggs","NS","chicks": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL
    ## 'data.frame':    672 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13541 levels "A2ML1","A2ML2",..: 6111 1308 9481 5091 1161 3641 6527 12934 3577 7601 ...
    ##  $ padj     : num  0.001878 0.082398 0.000393 0.000332 0.000121 ...
    ##  $ logpadj  : num  2.73 1.08 3.41 3.48 3.92 ...
    ##  $ lfc      : num  1.5 1.33 1.27 1.18 1.11 ...
    ##  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "nest","NS","chicks": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL
    ## 'data.frame':    1180 obs. of  6 variables:
    ##  $ gene     : Factor w/ 13541 levels "A2ML1","A2ML2",..: 9297 6964 4022 1308 6111 6378 2962 8188 4869 9177 ...
    ##  $ padj     : num  0.000983 0.014127 0.066324 0.09881 0.049459 ...
    ##  $ logpadj  : num  3.01 1.85 1.18 1.01 1.31 ...
    ##  $ lfc      : num  1.445 1.145 1.051 1.021 0.968 ...
    ##  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction: Factor w/ 3 levels "nest","NS","eggs": 3 3 3 3 3 3 3 3 3 3 ...
    ## NULL

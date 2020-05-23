    library(tidyverse)

    ## ── Attaching packages ─────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

### Hypothesis-related data wrangle

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
                                     chicks = c("hatch", "n5", "n9", "extend" , "early" ),
                                     nest = c("bldg", "m.inc.d3",  "m.inc.d9",  "m.inc.d17",  "m.n2"),
                                     control = c("control"))) %>%
      drop_na()

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
      #print(str(DEGs))
      
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
      
      #createDEGdfhypothesis("hiloPRL", "hi", "lo", i)
      #createDEGdfhypothesis("external", "chicks", "eggs", i)
      #createDEGdfhypothesis("external", "chicks", "nest", i)
      #createDEGdfhypothesis("external", "eggs", "nest", i)

    }

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

    write.csv(colData, "../metadata/04_colData.csv", row.names = F)

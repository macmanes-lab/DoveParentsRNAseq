    library(tidyverse)

    ## ── Attaching packages ─────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(readr)
    library(corrr)

    source("../R/themes.R")
    source("../R/functions.R")

Wrangle data for hyotheses testing
----------------------------------

    # col data with internal (hiloPRL) and external  hypothesese
    colData <- read_csv("../metadata/04_colData.csv") %>%
      column_to_rownames(var = "V1")

    ## Parsed with column specification:
    ## cols(
    ##   V1 = col_character(),
    ##   bird = col_character(),
    ##   sex = col_character(),
    ##   tissue = col_character(),
    ##   treatment = col_character(),
    ##   group = col_character(),
    ##   study = col_character(),
    ##   sextissue = col_character(),
    ##   samples = col_character(),
    ##   hiloPRL = col_character(),
    ##   external = col_character()
    ## )

    tail(colData)

    ##                                          bird    sex       tissue
    ## y98.g54_female_gonad_m.hatch          y98.g54 female       gonads
    ## y98.g54_female_hypothalamus_m.hatch   y98.g54 female hypothalamus
    ## y98.g54_female_pituitary_m.hatch      y98.g54 female    pituitary
    ## y98.o50.x_male_gonad_inc.d3         y98.o50.x   male       gonads
    ## y98.o50.x_male_hypothalamus_inc.d3  y98.o50.x   male hypothalamus
    ## y98.o50.x_male_pituitary_inc.d3     y98.o50.x   male    pituitary
    ##                                     treatment                    group
    ## y98.g54_female_gonad_m.hatch             m.n2       female.gonads.m.n2
    ## y98.g54_female_hypothalamus_m.hatch      m.n2 female.hypothalamus.m.n2
    ## y98.g54_female_pituitary_m.hatch         m.n2    female.pituitary.m.n2
    ## y98.o50.x_male_gonad_inc.d3            inc.d3       male.gonads.inc.d3
    ## y98.o50.x_male_hypothalamus_inc.d3     inc.d3 male.hypothalamus.inc.d3
    ## y98.o50.x_male_pituitary_inc.d3        inc.d3    male.pituitary.inc.d3
    ##                                               study           sextissue
    ## y98.g54_female_gonad_m.hatch           manipulation       female_gonads
    ## y98.g54_female_hypothalamus_m.hatch    manipulation female_hypothalamus
    ## y98.g54_female_pituitary_m.hatch       manipulation    female_pituitary
    ## y98.o50.x_male_gonad_inc.d3         charcterization         male_gonads
    ## y98.o50.x_male_hypothalamus_inc.d3  charcterization   male_hypothalamus
    ## y98.o50.x_male_pituitary_inc.d3     charcterization      male_pituitary
    ##                                                                 samples
    ## y98.g54_female_gonad_m.hatch               y98.g54_female_gonad_m.hatch
    ## y98.g54_female_hypothalamus_m.hatch y98.g54_female_hypothalamus_m.hatch
    ## y98.g54_female_pituitary_m.hatch       y98.g54_female_pituitary_m.hatch
    ## y98.o50.x_male_gonad_inc.d3                 y98.o50.x_male_gonad_inc.d3
    ## y98.o50.x_male_hypothalamus_inc.d3   y98.o50.x_male_hypothalamus_inc.d3
    ## y98.o50.x_male_pituitary_inc.d3         y98.o50.x_male_pituitary_inc.d3
    ##                                     hiloPRL external
    ## y98.g54_female_gonad_m.hatch             hi     nest
    ## y98.g54_female_hypothalamus_m.hatch      hi     nest
    ## y98.g54_female_pituitary_m.hatch         hi     nest
    ## y98.o50.x_male_gonad_inc.d3              lo     eggs
    ## y98.o50.x_male_hypothalamus_inc.d3       lo     eggs
    ## y98.o50.x_male_pituitary_inc.d3          lo     eggs

### variance stabilized gene expression (vsd)

    vsd_path <- "../results/DEseq2/hypothesis/"   # path to the data
    vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
    vsd_pathfiles <- paste0(vsd_path, vsd_files)
    vsd_files

    ## [1] "female_gonad_vsd.csv"        "female_gonads_vsd.csv"      
    ## [3] "female_hypothalamus_vsd.csv" "female_pituitary_vsd.csv"   
    ## [5] "male_gonad_vsd.csv"          "male_gonads_vsd.csv"        
    ## [7] "male_hypothalamus_vsd.csv"   "male_pituitary_vsd.csv"

    ## before pivoting, check names of df with
    ## head(names(allvsd)) and tail(names(allvsd))

    allvsd <- vsd_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name")  %>% 
      dplyr::rename("gene" = "X1") %>% 
      pivot_longer(cols = blk.s061.pu.y_female_gonad_inc.d9:y98.o50.x_male_pituitary_inc.d3, 
                   names_to = "samples", values_to = "counts") 
    allvsd %>% select(-file_name) %>% head()

    ## # A tibble: 6 x 3
    ##   gene  samples                           counts
    ##   <chr> <chr>                              <dbl>
    ## 1 A2ML1 blk.s061.pu.y_female_gonad_inc.d9   6.66
    ## 2 A2ML1 blk21.x_female_gonad_hatch          7.97
    ## 3 A2ML1 blk4.x_female_gonad_n9              7.37
    ## 4 A2ML1 blu103.x_female_gonad_hatch.NYNO    6.53
    ## 5 A2ML1 blu124.w180.x_female_gonad_hatch    7.21
    ## 6 A2ML1 blu36.w16_female_gonad_n9           7.34

    head(allvsd)

    ## # A tibble: 6 x 4
    ##   file_name                            gene  samples                 counts
    ##   <chr>                                <chr> <chr>                    <dbl>
    ## 1 ../results/DEseq2/hypothesis/female… A2ML1 blk.s061.pu.y_female_g…   6.66
    ## 2 ../results/DEseq2/hypothesis/female… A2ML1 blk21.x_female_gonad_h…   7.97
    ## 3 ../results/DEseq2/hypothesis/female… A2ML1 blk4.x_female_gonad_n9    7.37
    ## 4 ../results/DEseq2/hypothesis/female… A2ML1 blu103.x_female_gonad_…   6.53
    ## 5 ../results/DEseq2/hypothesis/female… A2ML1 blu124.w180.x_female_g…   7.21
    ## 6 ../results/DEseq2/hypothesis/female… A2ML1 blu36.w16_female_gonad…   7.34

All DEGs (but currently only char DEGs)
---------------------------------------

    DEG_path <- "../results/DEseq2/hypothesis/"   # path to the data
    DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
    DEG_pathfiles <- paste0(DEG_path, DEG_files)
    #DEG_files

    allDEG <- DEG_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name") %>% 
      mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/hypothesis/'), "[", 2))  %>% 
      mutate(DEG = sapply(strsplit(as.character(DEG),'_DEGs.csv'), "[", 1))  %>% 
      mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
      mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%
      select(DEG, sex, tissue, direction, everything())   %>%
      mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%
      mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%
      mutate(comparison = paste(down,up, sep = "_")) %>%
      mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
      mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%
      dplyr::select(sex,tissue,comparison, direction, gene, lfc, padj, logpadj)  %>%
      mutate(tissue = factor(tissue)) %>% 
      drop_na()

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_character(),
    ##   logpadj = col_character(),
    ##   lfc = col_character(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_character(),
    ##   logpadj = col_character(),
    ##   lfc = col_character(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_character(),
    ##   logpadj = col_character(),
    ##   lfc = col_character(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_character(),
    ##   logpadj = col_character(),
    ##   lfc = col_character(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    head(allDEG)

    ## # A tibble: 6 x 8
    ##   sex    tissue comparison  direction gene           lfc     padj logpadj
    ##   <chr>  <fct>  <chr>       <chr>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  eggs_chicks chicks    OVALX         6.17 7.05e- 3    2.15
    ## 2 female gonad  eggs_chicks chicks    LOC101749216  5.19 2.45e- 6    5.61
    ## 3 female gonad  eggs_chicks chicks    BPIFB2        4.95 5.23e- 3    2.28
    ## 4 female gonad  eggs_chicks chicks    OvoDA1        4.90 2.69e- 4    3.57
    ## 5 female gonad  eggs_chicks chicks    CDC20B        4.07 2.99e-10    9.52
    ## 6 female gonad  eggs_chicks chicks    WFIKKN2       4.01 1.90e- 2    1.72

gene lists
----------

    # candidate genes from GO and literature
    parentalcaregenes <- read_csv("../metadata/03_parentalcaregenes.csv") %>% pull(gene)

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_double(),
    ##   geneid = col_double(),
    ##   NCBI = col_character(),
    ##   literature = col_character(),
    ##   GO = col_character(),
    ##   gene = col_character()
    ## )

    parentalcaregenes

    ##  [1] "ADRA2A" "AVP"    "AVPR1A" "BRINP1" "COMT"   "CREBRF" "CRH"   
    ##  [8] "CRHBP"  "CRHR1"  "CRHR2"  "DBH"    "DRD1"   "DRD4"   "ESR1"  
    ## [15] "ESR2"   "FOS"    "GNAQ"   "HTR2C"  "KALRN"  "MBD2"   "MEST"  
    ## [22] "NPAS3"  "NPAS3"  "NR3C1"  "OPRK1"  "OPRM1"  "OXT"    "PGR"   
    ## [29] "PRL"    "PRLR"   "PTEN"   "SLC6A4" "ZFX"

    ## genes WGCNA prl module
    WGCNAgenes <- read_csv("../results/05_PRLmodule.csv") %>% pull(x)

    ## Parsed with column specification:
    ## cols(
    ##   x = col_character()
    ## )

    WGCNAgenes

    ##  [1] "ACOT12"       "AGR2"         "ARAP2"        "C2H8ORF46"   
    ##  [5] "CALN1"        "CD24"         "CDKN3"        "CDT1"        
    ##  [9] "CHRNB4"       "COL20A1"      "CREB3L1"      "CRELD2"      
    ## [13] "CREM"         "DUSP10"       "F13A1"        "FAM19A1"     
    ## [17] "FAM83G"       "FGF6"         "FKBP11"       "FOSL2"       
    ## [21] "FXYD2"        "GABRE"        "GRXCR1"       "LAMC2"       
    ## [25] "LAPTM4B"      "LBH"          "LOC101750767" "LOC107053040"
    ## [29] "LOC107054798" "LOC422224"    "LOC428479"    "LPPR5"       
    ## [33] "MFSD2A"       "MFSD4"        "MLN"          "MYC"         
    ## [37] "NPM3"         "NR4A3"        "NTN1"         "OMP"         
    ## [41] "PAX7"         "PDE11A"       "PDE8B"        "POU3F3"      
    ## [45] "PRL"          "PROCR"        "PYCR2"        "RASSF5"      
    ## [49] "SDF2L1"       "SLC16A1"      "SLC5A10"      "SLC7A5"      
    ## [53] "SOST"         "STC1"         "TFAP2B"       "TFPI2"       
    ## [57] "USP35"        "VEPH1"

    unique(allDEG$comparison)

    ## [1] "eggs_chicks" "lo_hi"       "nest_chicks" "nest_eggs"

    ## hilogenes

    hilogenes <- allDEG %>%
      filter(comparison == "lo_hi")  %>%
      drop_na() %>%
      arrange(padj) %>% 
      top_n(20) %>% pull(gene)

    ## Selecting by logpadj

    eggchickgenes <- allDEG %>%
      filter(comparison != "lo_hi")  %>%
      drop_na() %>%
      arrange(padj)  %>% 
      top_n(20) %>% pull(gene)

    ## Selecting by logpadj

### vsd for all and candidate, hypothesis, and data-driven genes

    candidategenes <- c(parentalcaregenes, WGCNAgenes, hilogenes, eggchickgenes)

    getcandidatevsd2 <- function(whichgenes, whichtissue, whichsex){
      candidates  <- allvsd %>%
        filter(gene %in% whichgenes) %>%
        dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
        dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/hypothesis/'), "[", 2)) %>%
        dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                      tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                      treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
        dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
        dplyr::select(sex, tissue, treatment, gene, samples, counts) %>%
        filter(tissue == whichtissue, sex %in% whichsex)  %>%
        drop_na()
      #candidates$treatment <- factor(candidates$treatment, levels = alllevels)
      return(candidates)
    }

    hypvsd <- getcandidatevsd2(candidategenes, "hypothalamus", sexlevels)
    pitvsd <- getcandidatevsd2(candidategenes, "pituitary", sexlevels)
    gonvsd <- getcandidatevsd2(candidategenes, "gonad", sexlevels)
    candidatevsd <- rbind(hypvsd, pitvsd)
    candidatevsd <- rbind(candidatevsd, gonvsd)

    head(candidatevsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue      treatment gene  samples                         counts
    ##   <chr>  <chr>       <chr>     <chr> <chr>                            <dbl>
    ## 1 female hypothalam… control   ABCA4 L.G118_female_hypothalamus_con…   5.91
    ## 2 female hypothalam… control   ABCA4 R.G106_female_hypothalamus_con…   5.88
    ## 3 female hypothalam… control   ABCA4 R.R20_female_hypothalamus_cont…   5.68
    ## 4 female hypothalam… control   ABCA4 R.R9_female_hypothalamus_contr…   5.66
    ## 5 female hypothalam… control   ABCA4 R.W44_female_hypothalamus_cont…   5.67
    ## 6 female hypothalam… prolong   ABCA4 blk.s031.pu.d_female_hypothala…   5.74

    write.csv(candidatevsd, "../results/06_candidatevsd.csv", row.names = F)
    write.csv(allDEG, "../results/06_allDEG.csv", row.names = F)

manipulation box plots for candidate genes
------------------------------------------

    library(tidyverse)

    ## ── Attaching packages ─────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(ggtext)
    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    library(ggpubr)

    ## Loading required package: magrittr

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     set_names

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

    ## 
    ## Attaching package: 'ggpubr'

    ## The following object is masked from 'package:cowplot':
    ## 
    ##     get_legend

    library(kableExtra)

    ## 
    ## Attaching package: 'kableExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

    source("../R/themes.R")

    knitr::opts_chunk$set(echo = TRUE, message = F, fig.path = "../figures/")

Candidate Genes (from literature and relevant GO terms)
-------------------------------------------------------

    source("../R/wrangledata.R")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

variance stabilized gene expression (vsd)
-----------------------------------------

    vsd_path <- "../results/DEseq2/manip/"   # path to the data
    vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
    vsd_pathfiles <- paste0(vsd_path, vsd_files)
    vsd_files

    ## [1] "female_gonad_vsd.csv"        "female_hypothalamus_vsd.csv"
    ## [3] "female_pituitary_vsd.csv"    "male_gonad_vsd.csv"         
    ## [5] "male_hypothalamus_vsd.csv"   "male_pituitary_vsd.csv"

    allvsd <- vsd_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name")  %>% 
      dplyr::rename("gene" = "X1") %>% 
      pivot_longer(cols = blk.s031.pu.d_female_gonad_prolong:y98.o50.x_male_pituitary_inc.d3, 
                   names_to = "samples", values_to = "counts") 

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    allvsd <- as.tibble(allvsd)

    ## Warning: `as.tibble()` is deprecated, use `as_tibble()` (but mind the new semantics).
    ## This warning is displayed once per session.

    getcandidatevsd <- function(whichgenes, whichtissue, whichsex){
      candidates  <- allvsd %>%
        filter(gene %in% whichgenes) %>%
        dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
        dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/manip/'), "[", 2)) %>%
        dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                    tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                    treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
        #dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1))# %>%
       # dplyr::recode(treatment, "m.hatch" = "m.n2", .default = levels(treatment)) 
        
         dplyr::mutate(treatment = recode(treatment, "m.hatch" = "m.n2"))%>%
        
        dplyr::select(sex, tissue, treatment, gene, samples, counts) %>%
        filter(tissue == whichtissue, sex %in% whichsex)  %>%
        drop_na()
      candidates$treatment <- factor(candidates$treatment, levels = alllevels)
      return(candidates)
    }

    hypvsd <- getcandidatevsd(candidategenes, "hypothalamus", sexlevels)
    pitvsd <- getcandidatevsd(candidategenes, "pituitary", sexlevels)
    gonvsd <- getcandidatevsd(candidategenes, "gonad", sexlevels)

    candidatevsd <- rbind(hypvsd, pitvsd)
    candidatevsd <- rbind(candidatevsd, gonvsd)

    head(candidatevsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue      treatment gene  samples                         counts
    ##   <chr>  <chr>       <fct>     <chr> <chr>                            <dbl>
    ## 1 female hypothalam… prolong   ADRA… blk.s031.pu.d_female_hypothala…   9.20
    ## 2 female hypothalam… m.n2      ADRA… blk.s032.g.w_female_hypothalam…   9.07
    ## 3 female hypothalam… m.inc.d3  ADRA… blk.s049.y.g_female_hypothalam…   9.21
    ## 4 female hypothalam… m.inc.d3  ADRA… blk.s060.pu.w_female_hypothala…   9.13
    ## 5 female hypothalam… inc.d9    ADRA… blk.s061.pu.y_female_hypothala…   9.19
    ## 6 female hypothalam… m.inc.d8  ADRA… blk.y.l.s109_female_hypothalam…   9.31

    levels(candidatevsd$treatment)

    ##  [1] "control"   "bldg"      "lay"       "inc.d3"    "m.inc.d3" 
    ##  [6] "inc.d9"    "m.inc.d8"  "m.inc.d9"  "inc.d17"   "m.inc.d17"
    ## [11] "prolong"   "hatch"     "m.n2"      "extend"    "n5"       
    ## [16] "n9"

DEGs
----

    DEG_path <- "../results/DEseq2/manip/"   # path to the data
    DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
    DEG_pathfiles <- paste0(DEG_path, DEG_files)
    #DEG_files

    allDEG <- DEG_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name") %>% 
      mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/manip/'), "[", 2))  %>% 
      mutate(DEG = sapply(strsplit(as.character(DEG),'_diffexp.csv'), "[", 1))  %>% 
      mutate(tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1)) %>%
      mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%
      mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%
      mutate(comparison = paste(down,up, sep = "_")) %>%
      mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
      mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%
    dplyr::select(sex,tissue,comparison, direction, gene, lfc, padj, logpadj) 
    head(allDEG)

    ## # A tibble: 6 x 8
    ##   sex    tissue comparison   direction gene           lfc     padj logpadj
    ##   <chr>  <chr>  <chr>        <chr>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  hatch_extend extend    LOC107053414 19.2  5.48e-13   12.3 
    ## 2 female gonad  hatch_extend extend    COL10A1       5.68 6.14e- 4    3.21
    ## 3 female gonad  hatch_extend extend    LOC107049904  5.30 6.20e- 2    1.21
    ## 4 female gonad  hatch_extend extend    SULT1C3       5.22 1.91e- 3    2.72
    ## 5 female gonad  hatch_extend extend    KRT20         4.87 1.46e- 2    1.84
    ## 6 female gonad  hatch_extend extend    LOC107057630  4.14 8.24e- 2    1.08

    allDEG$tissue <- factor(allDEG$tissue , levels = tissuelevel)
    allDEG$comparison <- factor(allDEG$comparison , levels = comparisonlevels)
    allDEG$direction <- factor(allDEG$direction, levels = alllevels)

    head(allDEG)

    ## # A tibble: 6 x 8
    ##   sex    tissue comparison   direction gene           lfc     padj logpadj
    ##   <chr>  <fct>  <fct>        <fct>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  hatch_extend extend    LOC107053414 19.2  5.48e-13   12.3 
    ## 2 female gonad  hatch_extend extend    COL10A1       5.68 6.14e- 4    3.21
    ## 3 female gonad  hatch_extend extend    LOC107049904  5.30 6.20e- 2    1.21
    ## 4 female gonad  hatch_extend extend    SULT1C3       5.22 1.91e- 3    2.72
    ## 5 female gonad  hatch_extend extend    KRT20         4.87 1.46e- 2    1.84
    ## 6 female gonad  hatch_extend extend    LOC107057630  4.14 8.24e- 2    1.08

    candidateDEGS <- allDEG %>%
      filter(gene %in% candidategenes) %>%
      mutate(posneg = ifelse(lfc >= 0, "+", "-"),
             sex = recode(sex, "female" = "F", "male" = "M" ),
             tissue = recode(tissue, 
                             "hypothalamus" = "H",
                             "pituitary" = "P", "gonad" = "G")) %>%
      mutate(res = paste(sex, tissue, posneg, sep = "")) %>%
      select(gene, res, comparison)  %>%
      group_by(gene,  comparison) %>%
      summarize(res = str_c(res, collapse = " ")) %>%
      pivot_wider(names_from = comparison, values_from = res) %>%
      select(gene, contains("m.")) 
    candidateDEGS

    ## # A tibble: 32 x 16
    ## # Groups:   gene [32]
    ##    gene  m.inc.d17_prolo… inc.d17_m.inc.d… inc.d9_m.inc.d8 hatch_m.n2
    ##    <chr> <chr>            <chr>            <chr>           <chr>     
    ##  1 ADRA… MG+              <NA>             <NA>            <NA>      
    ##  2 AVP   <NA>             FH- MH-          MP+             <NA>      
    ##  3 AVPR… <NA>             <NA>             <NA>            FH+ MP+   
    ##  4 BRIN… <NA>             FG- FH+          <NA>            MH+       
    ##  5 COMT  <NA>             FH- MH-          MH-             FH- MH-   
    ##  6 CREB… <NA>             FP+ MP+          MP+             FH+ FP+ M…
    ##  7 CRH   <NA>             <NA>             <NA>            <NA>      
    ##  8 CRHBP <NA>             FH+ MH+          MH+             FH+ MH+   
    ##  9 CRHR1 <NA>             <NA>             MH+             FH+ FP+ M…
    ## 10 CRHR2 <NA>             FH+ MH+          MH+             FH+ MH+   
    ## # … with 22 more rows, and 11 more variables: m.inc.d3_m.inc.d17 <chr>,
    ## #   m.inc.d3_m.n2 <chr>, inc.d3_m.inc.d3 <chr>, m.inc.d3_m.inc.d9 <chr>,
    ## #   m.inc.d8_extend <chr>, m.inc.d8_prolong <chr>,
    ## #   m.inc.d9_m.inc.d8 <chr>, m.inc.d9_m.n2 <chr>, inc.d9_m.inc.d9 <chr>,
    ## #   m.n2_extend <chr>, m.inc.d9_m.inc.d17 <chr>

    table2 <- left_join(candidateDEGS, GOgenesWide) %>%
      select(gene, m.inc.d17_prolong:m.inc.d9_m.inc.d17, parentalcare, parentalbehavior, NCBI) %>%
      mutate(parentalcare = if_else(is.na(parentalcare), " ", "X"),
             parentalbehavior = if_else(is.na(parentalbehavior), " ", "X")) %>%
      rename("Literature" = "parentalcare", "GO" =  "parentalbehavior")
    table2$numDEGs <- rowSums(is.na(table2)) # count NAs to know how many are NS
    table2 <- table2 %>% 
      mutate(sig = ifelse(numDEGs == 6, "NS", "DEG")) %>% 
      arrange(sig, gene)  %>%  select(-sig, -numDEGs)
    table2[is.na(table2)] <- " " # replace NA with blank space so it's pretty
    table2

    ## # A tibble: 33 x 19
    ## # Groups:   gene [32]
    ##    gene  m.inc.d17_prolo… inc.d17_m.inc.d… inc.d9_m.inc.d8 hatch_m.n2
    ##    <chr> <chr>            <chr>            <chr>           <chr>     
    ##  1 ADRA… "MG+"            " "              " "             " "       
    ##  2 AVP   " "              "FH- MH-"        "MP+"           " "       
    ##  3 AVPR… " "              " "              " "             "FH+ MP+" 
    ##  4 BRIN… " "              "FG- FH+"        " "             "MH+"     
    ##  5 COMT  " "              "FH- MH-"        "MH-"           "FH- MH-" 
    ##  6 CREB… " "              "FP+ MP+"        "MP+"           "FH+ FP+ …
    ##  7 CRH   " "              " "              " "             " "       
    ##  8 CRHBP " "              "FH+ MH+"        "MH+"           "FH+ MH+" 
    ##  9 CRHR1 " "              " "              "MH+"           "FH+ FP+ …
    ## 10 CRHR2 " "              "FH+ MH+"        "MH+"           "FH+ MH+" 
    ## # … with 23 more rows, and 14 more variables: m.inc.d3_m.inc.d17 <chr>,
    ## #   m.inc.d3_m.n2 <chr>, inc.d3_m.inc.d3 <chr>, m.inc.d3_m.inc.d9 <chr>,
    ## #   m.inc.d8_extend <chr>, m.inc.d8_prolong <chr>,
    ## #   m.inc.d9_m.inc.d8 <chr>, m.inc.d9_m.n2 <chr>, inc.d9_m.inc.d9 <chr>,
    ## #   m.n2_extend <chr>, m.inc.d9_m.inc.d17 <chr>, Literature <chr>,
    ## #   GO <chr>, NCBI <chr>

    write.csv(table2, "../results/table2.csv")

    # for victoria
    prlmanip <- candidatevsd %>%
      filter(gene == "PRL", 
             treatment %in% c("extend", "m.inc.d8"))
    write.csv(prlmanip, "../results/prlmanip.csv")

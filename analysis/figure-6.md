manipulation box plots for candidate genes
------------------------------------------

    library(tidyverse)

    ## ── Attaching packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

    head(allvsd)

    ## # A tibble: 6 x 4
    ##   file_name                         gene  samples                    counts
    ##   <chr>                             <chr> <chr>                       <dbl>
    ## 1 ../results/DEseq2/manip/female_g… A2ML1 blk.s031.pu.d_female_gona…   7.60
    ## 2 ../results/DEseq2/manip/female_g… A2ML1 blk.s032.g.w_female_gonad…   7.49
    ## 3 ../results/DEseq2/manip/female_g… A2ML1 blk.s049.y.g_female_gonad…   7.81
    ## 4 ../results/DEseq2/manip/female_g… A2ML1 blk.s060.pu.w_female_gona…   7.16
    ## 5 ../results/DEseq2/manip/female_g… A2ML1 blk.s061.pu.y_female_gona…   7.19
    ## 6 ../results/DEseq2/manip/female_g… A2ML1 blk.y.l.s109_female_gonad…   7.50

    # for plotting gene expression over time
    getcandidatevsd <- function(whichgenes, whichtissue, whichsex){
      candidates  <- allvsd %>%
        filter(gene %in% whichgenes) %>%
        dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
        dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/manip/'), "[", 2)) %>%
        dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                    tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                    treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
        dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
         dplyr::mutate(treatment = recode(treatment, "m.hatch" = "m.n2",
                                                      "extend.hatch" = "m.n2",
                                          "inc.prolong" = "prolong")) %>%
        
        dplyr::mutate(day = fct_collapse(treatment,
                                         "1" = c("inc.d3", "m.inc.d3" ),
                                         "2" = c("inc.d9", "m.inc.d9", "m.inc.d8"  ),
                                         "2" = c("inc.d17", "m.inc.d17", "prolong"  ),
                                         "4" = c("hatch", "m.n2", "extend"  )
                                         )) %>%
        dplyr::mutate(day = as.numeric(day)) %>%
        dplyr::mutate(manip = fct_collapse(treatment,
                                           "reference" = charlevels,
                                           "removal" = c("m.inc.d3", "m.inc.d9", "m.inc.d17","m.n2"),
                                           "replacement" = c("m.inc.d8", "prolong", "extend"))) %>%
        dplyr::mutate(manip = factor(manip, levels = c("reference", "removal", "replacement"))) %>%
        dplyr::select(sex, tissue, treatment, day,  manip, gene, samples, counts) %>%
        filter(tissue == whichtissue, sex %in% whichsex)  %>%
        drop_na()
      candidates$treatment <- factor(candidates$treatment, levels = alllevels)
      return(candidates)
    }

    hypvsd <- getcandidatevsd(candidategenes, "hypothalamus", sexlevels)

    ## Warning: Unknown levels in `f`: control, bldg, lay, n5, n9

    pitvsd <- getcandidatevsd(candidategenes, "pituitary", sexlevels)

    ## Warning: Unknown levels in `f`: control, bldg, lay, n5, n9

    gonvsd <- getcandidatevsd(candidategenes, "gonad", sexlevels)

    ## Warning: Unknown levels in `f`: control, bldg, lay, n5, n9

    candidatevsd <- rbind(hypvsd, pitvsd)
    candidatevsd <- rbind(candidatevsd, gonvsd)
    head(candidatevsd)

    ## # A tibble: 6 x 8
    ##   sex    tissue    treatment   day manip   gene  samples             counts
    ##   <chr>  <chr>     <fct>     <dbl> <fct>   <chr> <chr>                <dbl>
    ## 1 female hypothal… prolong       2 replac… ADRA… blk.s031.pu.d_fema…   9.20
    ## 2 female hypothal… m.n2          1 removal ADRA… blk.s032.g.w_femal…   9.07
    ## 3 female hypothal… m.inc.d3      3 removal ADRA… blk.s049.y.g_femal…   9.21
    ## 4 female hypothal… m.inc.d3      3 removal ADRA… blk.s060.pu.w_fema…   9.13
    ## 5 female hypothal… inc.d9        2 refere… ADRA… blk.s061.pu.y_fema…   9.19
    ## 6 female hypothal… m.inc.d8      2 replac… ADRA… blk.y.l.s109_femal…   9.31

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

    table2 %>% select(gene,inc.d3_m.inc.d3, inc.d9_m.inc.d9, inc.d17_m.inc.d17, hatch_m.n2)

    ## # A tibble: 33 x 5
    ## # Groups:   gene [32]
    ##    gene   inc.d3_m.inc.d3 inc.d9_m.inc.d9 inc.d17_m.inc.d17 hatch_m.n2   
    ##    <chr>  <chr>           <chr>           <chr>             <chr>        
    ##  1 ADRA2A " "             " "             " "               " "          
    ##  2 AVP    " "             " "             "FH- MH-"         " "          
    ##  3 AVPR1A " "             " "             " "               "FH+ MP+"    
    ##  4 BRINP1 "FH+"           " "             "FG- FH+"         "MH+"        
    ##  5 COMT   "FH-"           " "             "FH- MH-"         "FH- MH-"    
    ##  6 CREBRF " "             "MG+"           "FP+ MP+"         "FH+ FP+ MP+"
    ##  7 CRH    "FH+"           " "             " "               " "          
    ##  8 CRHBP  " "             " "             "FH+ MH+"         "FH+ MH+"    
    ##  9 CRHR1  " "             " "             " "               "FH+ FP+ MH+"
    ## 10 CRHR2  "FH+"           " "             "FH+ MH+"         "FH+ MH+"    
    ## # … with 23 more rows

    plotcandidatemanip <- function(whichtissue, whichsex, whichgenes){
      
      p <- candidatevsd %>%
      filter(tissue == whichtissue, sex == whichsex) %>%
      filter(gene %in% whichgenes) %>%
      filter(treatment %in% c(levelsremoval, controlsremoval)) %>%
      mutate(treatment = factor(treatment, levels = alllevels3)) %>%
      ggplot(aes(y =  counts, x = treatment, fill = treatment, color = sex)) +
      geom_boxplot() +
      facet_wrap(~gene, scales = "free_y", nrow = 2) +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(face = "italic")) +
      scale_color_manual(values = allcolors) +
      scale_fill_manual(values = allcolors) 
      
      
      return(p)
    }

    p1 <- plotcandidatemanip("hypothalamus", "female", c("AVP", "BRINP1", "COMT", "CREBRF", "CRH", "CRHBP", "CRHR1", "CRHR2", "OPRK1")) +
      labs(subtitle = "hypothalamus") +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank()) 
    p2 <- plotcandidatemanip("pituitary", "female", c("CREBRF", "CRHR1", "DRD1", "DRD4", "GNAQ")) +
      labs(subtitle = "pituitary")
    p3 <- plotcandidatemanip("gonad", "female", c("BRINP1", "OPRK1"))  +
      labs(subtitle = "gonads")


    p23 <- plot_grid(p2,p3, nrow = 1, rel_widths = c(3,1))

    plot_grid(p1,p23, nrow = 2, rel_heights = c(1,1.2))

![](../figures/fig6-1.png)

    ## hypothesis

    ## not working
    meanse <- candidatevsd %>%
      dplyr::group_by(sex, tissue, treatment, day, manip, gene) %>%
      dplyr::summarise(m = mean(counts), 
                       se = sd(counts)/sqrt(length(counts))) %>%
      dplyr::mutate(m = round(m,0))
    head(meanse)

    meanse %>%
      filter(gene == "PRL") %>%
      ggplot(aes(x = day, y = m, color = manip)) +
      geom_errorbar(aes(ymin=m-se, ymax=m+se), width=.1) +
        geom_point(size = 1) +
      facet_wrap(sex~tissue, scales = "free")

    write.csv(table2, "../results/table2.csv")

    # for victoria
    prlmanip <- candidatevsd %>%
      filter(gene == "PRL", 
             treatment %in% c("extend", "m.inc.d8"))
    write.csv(prlmanip, "../results/prlmanip.csv")

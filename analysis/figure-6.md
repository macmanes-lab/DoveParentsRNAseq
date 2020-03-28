manipulation box plots for candidate genes
------------------------------------------

    library(tidyverse)

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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


    candidategenesslim <- c("OXT", "AVP", "GNRH1", "GNRHR", 
                        "AR",  "CYP19A1", 
                         "AVPR1A", "AVPR1B", "AVPR2","VIP",
                      "DRD1", "DRD2", 
                      "PRL", "PRLR",  
                        "ESR1","ESR2", "LBH",  
                       
                        "FOS", "JUN", "EGR1", "BDNF"
                          ) 

    candidategenes <- candidategenesslim

variance stabilized gene expression (vsd)
-----------------------------------------

    geneids <- read_csv("../metadata/00_geneinfo.csv")

    ## Warning: Missing column names filled in: 'X1' [1]

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

    hypvsd <- getcandidatevsd(candidategenesslim, "hypothalamus", sexlevels)
    pitvsd <- getcandidatevsd(candidategenesslim, "pituitary", sexlevels)
    gonvsd <- getcandidatevsd(candidategenesslim, "gonad", sexlevels)
    head(hypvsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue     treatment gene  samples                          counts
    ##   <chr>  <chr>      <fct>     <chr> <chr>                             <dbl>
    ## 1 female hypothala… prolong   AR    blk.s031.pu.d_female_hypothalam…   7.60
    ## 2 female hypothala… m.n2      AR    blk.s032.g.w_female_hypothalamu…   7.65
    ## 3 female hypothala… m.inc.d3  AR    blk.s049.y.g_female_hypothalamu…   7.80
    ## 4 female hypothala… m.inc.d3  AR    blk.s060.pu.w_female_hypothalam…   7.69
    ## 5 female hypothala… inc.d9    AR    blk.s061.pu.y_female_hypothalam…   7.64
    ## 6 female hypothala… m.inc.d8  AR    blk.y.l.s109_female_hypothalamu…   7.71

    head(pitvsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue   treatment gene  samples                            counts
    ##   <chr>  <chr>    <fct>     <chr> <chr>                               <dbl>
    ## 1 female pituita… prolong   AR    blk.s031.pu.d_female_pituitary_pr…  10.1 
    ## 2 female pituita… m.n2      AR    blk.s032.g.w_female_pituitary_m.h…   9.48
    ## 3 female pituita… m.inc.d3  AR    blk.s049.y.g_female_pituitary_m.i…  10.5 
    ## 4 female pituita… <NA>      AR    blk.s060.pu.w_female_pituitary_m.…  10.5 
    ## 5 female pituita… inc.d9    AR    blk.s061.pu.y_female_pituitary_in…  10.2 
    ## 6 female pituita… m.inc.d8  AR    blk.y.l.s109_female_pituitary_m.i…   9.72

    head(gonvsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue treatment gene  samples                             counts
    ##   <chr>  <chr>  <fct>     <chr> <chr>                                <dbl>
    ## 1 female gonad  prolong   AR    blk.s031.pu.d_female_gonad_prolong    9.63
    ## 2 female gonad  m.n2      AR    blk.s032.g.w_female_gonad_m.hatch     9.45
    ## 3 female gonad  m.inc.d3  AR    blk.s049.y.g_female_gonad_m.inc.d3    8.89
    ## 4 female gonad  m.inc.d3  AR    blk.s060.pu.w_female_gonad_m.inc.d3   9.32
    ## 5 female gonad  inc.d9    AR    blk.s061.pu.y_female_gonad_inc.d9     9.52
    ## 6 female gonad  m.inc.d8  AR    blk.y.l.s109_female_gonad_m.inc.d8    9.80

    candidatevsd <- rbind(hypvsd, pitvsd)
    candidatevsd <- rbind(candidatevsd, gonvsd)

    head(candidatevsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue     treatment gene  samples                          counts
    ##   <chr>  <chr>      <fct>     <chr> <chr>                             <dbl>
    ## 1 female hypothala… prolong   AR    blk.s031.pu.d_female_hypothalam…   7.60
    ## 2 female hypothala… m.n2      AR    blk.s032.g.w_female_hypothalamu…   7.65
    ## 3 female hypothala… m.inc.d3  AR    blk.s049.y.g_female_hypothalamu…   7.80
    ## 4 female hypothala… m.inc.d3  AR    blk.s060.pu.w_female_hypothalam…   7.69
    ## 5 female hypothala… inc.d9    AR    blk.s061.pu.y_female_hypothalam…   7.64
    ## 6 female hypothala… m.inc.d8  AR    blk.y.l.s109_female_hypothalamu…   7.71

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

    table1 <- allDEG %>%
      mutate(updown = ifelse(lfc > 0, "+", "-"))  %>%
      mutate(geneupdown = paste(gene, updown, sep = "")) %>%
      filter(gene %in% candidategenes) %>%
      arrange(geneupdown) %>%
      group_by(sex, tissue, comparison) %>%
      summarize(genes = str_c(geneupdown, collapse = " ")) %>%
      pivot_wider(names_from = comparison, values_from = genes ) %>%
      select(sex, tissue, contains("m."))  %>%
      arrange( tissue, sex)

    (table1$inc.d3_m.inc.d3)

    ## [1] "BDNF+ DRD1+ EGR1+ ESR2- FOS+ PRLR+"
    ## [2] NA                                  
    ## [3] "DRD1+"                             
    ## [4] NA                                  
    ## [5] NA                                  
    ## [6] NA

    (table1$inc.d9_m.inc.d9)

    ## [1] NA     NA     "VIP+" NA     NA     NA

    (table1$inc.d17_m.inc.d17)

    ## [1] "AVP- BDNF+ DRD1+ EGR1+ ESR2- FOS+ OXT- PRLR+"    
    ## [2] "AR- AVP- BDNF+ DRD1+ EGR1+ ESR2- LBH- PRL- PRLR+"
    ## [3] "AR+ BDNF+ LBH- PRL- PRLR-"                       
    ## [4] "DRD1+ ESR2- LBH-"                                
    ## [5] "JUN-"                                            
    ## [6] "LBH-"

    (table1$hatch_m.n2)

    ## [1] "AVPR1A+ CYP19A1+ DRD1+ EGR1+ GNRH1+ JUN- PRL- PRLR+"
    ## [2] "DRD1+ EGR1+ GNRH1+ PRL- PRLR+"                      
    ## [3] "AR+ AVPR1B+ AVPR2+ EGR1+ GNRHR+ JUN+ LBH- PRL-"     
    ## [4] "AVPR1A+ LBH- PRL- VIP+"                             
    ## [5] NA                                                   
    ## [6] NA

Figs
----

    makeboxplotsmanip <- function(df, whichgene, mysubtitle, whichsex, whichtimepoint){
      p <- df %>%
        filter(treatment %in% alllevels,
               gene %in% whichgene,
               sex %in% whichsex) %>%
        filter(treatment %in% whichtimepoint) %>%
        ggplot(aes(x = treatment, y = counts, fill = treatment, color = sex)) +
        geom_boxplot(outlier.shape = NA) + 
        geom_jitter(size = 0.5, width = 0.1) +
        facet_wrap(~gene, scales = "free_y", nrow = 1) +
        theme_B3() +
        scale_fill_manual(values = allcolors) +
        scale_color_manual(values = allcolors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none",
              plot.caption = element_text(face = "italic")) +
        labs(y = "gene expression" , x = "Parental Stages and Manipulations", subtitle = mysubtitle) +
        theme(strip.text = element_text(face = "italic"))
      return(p)
    }

    ## hyp 
    a2 <- makeboxplotsmanip(hypvsd, c("DRD1","EGR1","PRLR"), "Female hypothalamus", "female",
                      c("inc.d3", "m.inc.d3", "inc.d17", "m.inc.d17",  "hatch", "m.n2"))  +
      theme(axis.title.x = element_blank())
    b2 <- makeboxplotsmanip(hypvsd, c("BDNF","ESR2", "FOS"), " ", "female",
                      c("inc.d3", "m.inc.d3", "inc.d17", "m.inc.d17"))  +
      theme(axis.title.x = element_blank())

    c <- makeboxplotsmanip(hypvsd, c("AVP","OXT"), " ", "female", 
                           c("inc.d17", "m.inc.d17"))  +
      theme(axis.title = element_blank())
    d <- makeboxplotsmanip(hypvsd, c("AVPR1A","CYP19A1","GNRH1"), " ", "female",
                       c("hatch", "m.n2"))  +
      theme(axis.title = element_blank())
    e <- makeboxplotsmanip(hypvsd, c("JUN"), " ", "female", 
                     c("hatch", "m.n2"))  +
      theme(axis.title = element_blank()) + labs(subtitle = " ")



    f1 <- makeboxplotsmanip(hypvsd, c("DRD1", "EGR1", "PRLR"), "Male hypothalamus", "male",
                      c( "inc.d17", "m.inc.d17", "hatch", "m.n2"))  +
      theme(axis.title.x = element_blank())

    f2 <- makeboxplotsmanip(hypvsd, c( "PRL"), " ", "male",
                      c( "inc.d17", "m.inc.d17", "hatch", "m.n2"))  +
      theme(axis.title = element_blank())


    h <- makeboxplotsmanip(hypvsd, c("AR", "AVP", "BDNF", "ESR1", "LBH"), " ", "male",
                      c( "inc.d17", "m.inc.d17"))  + labs(caption = " ")

    i <- makeboxplotsmanip(hypvsd, c("GNRH1"), " ", "male",
                      c("hatch", "m.n2"))  +
      theme(axis.title = element_blank()) +
      labs(x = " ",
           caption = "Candidate genes that were not differentially expressed following removal: AVPR1B, AVPR2, ESR1, GNRHR, VIP")

    ab2 <- plot_grid(a2,c, rel_widths = c(18,5), labels = "auto", label_size = 8)
    cde <- plot_grid(b2,d,e, rel_widths = c(12,6,2), nrow =  1, align = "h",  
                     labels = c("c", "d", " "), label_size = 8)
    f12 <- plot_grid(f1,f2, rel_widths = c(3,1), labels = c("e",  " "), label_size = 8)
    hi <- plot_grid(h,i, ncol = 2, rel_widths = c(10,2),labels = c("f",  "g"), label_size = 8, align = "h")

    plot_grid(ab2, cde, f12,hi, nrow = 4, rel_heights = c(1,1,1,1.25))

![](../figures/fig6-1.png)

    write.csv(candidatevsd, "../../musicalgenes/data/candidatecounts.csv")

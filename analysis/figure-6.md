Figure 4
========

    library(tidyverse)

    ## ── Attaching packages ──────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    source("../R/themes.R")

    knitr::opts_chunk$set(fig.path = '../figures/',message=F, warning=FALSE)

    geneids <- read_csv("../metadata/00_geneinfo.csv")

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
      dplyr::rename("gene" = "X1")  %>% 
      pivot_longer(cols = blk.s031.pu.d_female_gonad_prolong:y98.o50.x_male_pituitary_inc.d3, 
                   names_to = "samples", values_to = "counts") 

    hypreplace <- read_csv("../results/DEseq2/manip/female_hypothalamus_hatch_extend_DEGs.csv") %>%
      arrange(desc(logpadj)) %>% head(., n = 1) %>% pull(gene)
    hypreplace

    ## [1] "CAPN2"

    hypremove <- read_csv("../results/DEseq2/manip/female_hypothalamus_m.n2_hatch_DEGs.csv") %>%
      arrange(desc(logpadj)) %>% head(., n = 1) %>% pull(gene)
    hypremove

    ## [1] "CTNNAL1"

    pitreplace <- read_csv("../results/DEseq2/manip/female_pituitary_hatch_extend_DEGs.csv") %>%
      arrange(desc(logpadj)) %>% head(., n = 1) %>% pull(gene)
    pitreplace

    ## [1] "YIPF1"

    pitremove <- read_csv("../results/DEseq2/manip/female_pituitary_m.n2_hatch_DEGs.csv") %>%
      arrange(desc(logpadj)) %>% head(., n = 1) %>% pull(gene)
    pitremove

    ## [1] "ANKH"

    gonreplace <- read_csv("../results/DEseq2/manip/female_gonad_inc.d17_prolong_DEGs.csv") %>%
      arrange(desc(logpadj)) %>% head(., n = 1) %>% pull(gene)
    gonreplace

    ## [1] "TEKT3"

    gonremove <- read_csv("../results/DEseq2/manip/female_gonad_inc.d17_m.inc.d17_DEGs.csv") %>%
      arrange(desc(logpadj)) %>% head(., n = 1) %>% pull(gene)
    gonremove

    ## [1] "LOC100858655"

    plottopgenes <- function(whichgenes, whichtissue, whichsex, mysubtitle, whichtstages){
      candidates  <- allvsd %>%
        filter(gene %in% whichgenes) %>%
        dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
        dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/manip/'), "[", 2)) %>%
        dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                    tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                    treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
        dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
        dplyr::select(sex, tissue, treatment, gene, samples, counts) %>%
        filter(tissue == whichtissue, sex %in% whichsex) 
      
      candidates$treatment <- factor(candidates$treatment, levels = alllevels2)
      
      p <- candidates %>%
        filter(treatment %in% whichtstages) %>%
        ggplot(aes(x = treatment, y = counts, fill = treatment, color = sex)) +
        geom_boxplot() + facet_wrap(~gene, scales = "free_y") + 
        theme_B3() +
        scale_fill_manual(values = allcolors) +
        scale_color_manual(values = allcolors) +
        theme(axis.text.x = element_blank(),
              legend.position = "none",
              axis.title.x = element_blank(),
              strip.text = element_text(face = "italic")) +
        labs(y = mysubtitle) 
      return(p)
      
    }


    a <- plottopgenes(hypreplace, "hypothalamus", c("female"), "hypothalamus", c("hatch", "extend")) 
    b <- plottopgenes(pitreplace, "pituitary", c("female"), "pituitary", c("hatch", "extend")) 
    c <- plottopgenes(gonreplace, "gonad", c("female"), "gonads", c("inc.d17", "prolong"))

    fig6a <- plot_grid(a,b,c, ncol = 1) 

    d <- plottopgenes(hypreplace, "hypothalamus", c("female", "male"), NULL, alllevels2) 
    e <- plottopgenes(pitreplace, "pituitary", c("female", "male"), NULL, alllevels2) 
    f <- plottopgenes(gonreplace, "gonad", c("female", "male"), NULL, alllevels2) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    fig6b <- plot_grid(d,e,f, ncol = 1)

    fig6 <- plot_grid(fig6a, fig6b, rel_widths = c(1,2))
    fig6

![](../figures/fig6-1.png)

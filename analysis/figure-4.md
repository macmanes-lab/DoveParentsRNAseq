Figure 4
========

    library(tidyverse)

    ## ── Attaching packages ─────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

    vsd_path <- "../results/DEseq2/"   # path to the data
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
      pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
                   names_to = "samples", values_to = "counts") 

    hypgenes <- read_csv("../results/DEseq2/female_hypothalamus_hatch_n5_DEGs.csv") %>%
      arrange(desc(logpadj)) %>% head(., n = 2) %>% pull(gene)
    hypgenes

    ## [1] "CAPN2"        "LOC107054855"

    pitgenes <- read_csv("../results/DEseq2/female_pituitary_inc.d9_inc.d17_DEGs.csv") %>%
      arrange(desc(logpadj)) %>% head(., n = 2) %>% pull(gene)
    pitgenes

    ## [1] "LBH"  "CDK1"

    gongenes <- read_csv("../results/DEseq2/female_gonad_lay_inc.d3_DEGs.csv") %>%
      arrange(desc(logpadj)) %>% head(., n = 2) %>% pull(gene)
    gongenes

    ## [1] "OVSTL"  "BPIFB2"

    candidategenes <- c("PRL", "PRLR" ,  "AVPR1B" , "MYC" ,  
                    "MAP2KA", "MSH" ,   "PALB2" , "BRCA1", 
                    "OXT", "ESR1", "ESR2", "CISH", 
                    "POU2F1", "RAP1GAP", "CEBPB","ARHGAP32"
                    )

    returnvsds <- function(whichgenes){
      
      df  <- allvsd %>%
        filter(gene %in% whichgenes) %>%
        dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
        dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/'), "[", 2)) %>%
        dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                    tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                    treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
        dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
        dplyr::select(sex, tissue, treatment, gene, samples, counts) 
      
      df$treatment <- factor(df$treatment, levels = alllevels)
      print(head(df))
      return(df)
    }

    candidatevsd <- returnvsds(candidategenes)

    ## # A tibble: 6 x 6
    ##   sex    tissue treatment gene     samples                           counts
    ##   <chr>  <chr>  <fct>     <chr>    <chr>                              <dbl>
    ## 1 female gonad  control   ARHGAP32 L.G118_female_gonad_control         8.03
    ## 2 female gonad  control   ARHGAP32 R.G106_female_gonad_control         8.36
    ## 3 female gonad  control   ARHGAP32 R.R20_female_gonad_control          8.21
    ## 4 female gonad  control   ARHGAP32 R.R9_female_gonad_control           7.74
    ## 5 female gonad  control   ARHGAP32 R.W44_female_gonad_control          8.68
    ## 6 female gonad  inc.d9    ARHGAP32 blk.s061.pu.y_female_gonad_inc.d9   8.81

    plottopgenes <- function(whichgenes, whichtissue, whichsex, mysubtitle, whichtstages){
      candidates  <- allvsd %>%
        filter(gene %in% whichgenes) %>%
        dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
        dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/'), "[", 2)) %>%
        dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                    tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                    treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
        dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
        dplyr::select(sex, tissue, treatment, gene, samples, counts) %>%
        filter(tissue == whichtissue, sex %in% whichsex) 
      
      candidates$treatment <- factor(candidates$treatment, levels = alllevels)
      
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


    a <- plottopgenes(hypgenes, "hypothalamus", c("female"), "hypothalamus", c("hatch", "n5")) 
    b <- plottopgenes(pitgenes, "pituitary", c("female"), "pituitary", c("inc.d9", "inc.d17")) 
    c <- plottopgenes(gongenes, "gonad", c("female"), "gonads", c("lay", "inc.d3"))

    fig4a <- plot_grid(a,b,c, ncol = 1) 

    d <- plottopgenes(hypgenes, "hypothalamus", c("female", "male"), NULL, charlevels) 
    e <- plottopgenes(pitgenes, "pituitary", c("female", "male"), NULL, charlevels) 
    f <- plottopgenes(gongenes, "gonad", c("female", "male"), NULL, charlevels) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    fig4b <- plot_grid(d,e,f, ncol = 1)

    fig4 <- plot_grid(fig4a, fig4b, rel_widths = c(1,2))
    fig4

![](../figures/fig4-1.png)

    write.csv(candidatevsd, "../../musicalgenes/data/candidatecounts.csv")

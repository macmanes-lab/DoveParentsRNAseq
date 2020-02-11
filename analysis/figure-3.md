Plots with Prolactin
====================

    library(tidyverse)

    ## ── Attaching packages ───────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    library(ggimage)

    ## 
    ## Attaching package: 'ggimage'

    ## The following object is masked from 'package:cowplot':
    ## 
    ##     theme_nothing

    library(apaTables)
    library(factoextra)

    ## Welcome! Related Books: `Practical Guide To Cluster Analysis in R` at https://goo.gl/13EFCZ

    source("../R/themes.R") 
    source("../R/functions.R")
    source("../R/icons.R")

    ## Warning: Column `icons` joining factor and character vector, coercing into
    ## character vector

    source("../R/wrangledata.R")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )

    ## See spec(...) for full column specifications.

    knitr::opts_chunk$set(fig.path = '../figures/',message=F, warning=FALSE)

Data: PCA, hormones, PRL expression
-----------------------------------

    # pca
    charpca <- subsetmakepca(tissuelevels, charlevels, sexlevels)   
    charfviz <- makefvizdf(tissuelevels, charlevels, sexlevels)

    # hormones
    prolactin <- read_csv("../results/07_hormones.csv") %>%
        filter(study == "characterization", hormone %in% c("prolactin"))  %>% 
        droplevels() 
    prolactin$treatment <- factor(prolactin$treatment, levels = alllevels)

    # PRL vsd

    vsd_path <- "../results/DEseq2/"   # path to the data
    vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
    vsd_pathfiles <- paste0(vsd_path, vsd_files)
    vsd_files

    ## [1] "female_gonad_vsd.csv"        "female_hypothalamus_vsd.csv"
    ## [3] "female_pituitary_vsd.csv"    "male_gonad_vsd.csv"         
    ## [5] "male_hypothalamus_vsd.csv"   "male_pituitary_vsd.csv"

    PRLvsd <- vsd_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name")  %>% 
      dplyr::rename("gene" = "X1") %>% 
      pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
                   names_to = "samples", values_to = "counts") %>%
      filter(gene == "PRL")  %>%
      dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
      dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/'), "[", 2)) %>%
      dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                    tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                    treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
      dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
      dplyr::select(sex, tissue, treatment, gene, samples, counts)
    PRLvsd$treatment <- factor(PRLvsd$treatment, levels = alllevels)

    PRLhyp <- PRLvsd %>% filter(tissue == "hypothalamus")
    PRLpit <- PRLvsd %>% filter(tissue == "pituitary")
    PRLgon <- PRLvsd %>% filter(tissue == "gonad")

    a <- plotcolorfulpcs(charpca,charpca$treatment, allcolors) + labs(subtitle = " ") +
      theme(legend.position = c(0.6,0.1), 
            legend.direction = "horizontal", 
            legend.key.size = unit(0.5, 'lines')) + 
      guides(color = FALSE) +
      labs(subtitle = " ")   

    b <- plotprolactin(PRLhyp, PRLhyp$counts, "PRL", "hypothalamus") + 
      theme(legend.position = c(0.7,0.9), axis.text.x = element_blank(), 
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "italic"),
            legend.title = element_blank()) + guides(fill = F)

    c <- plotfriz(charfviz) + labs(subtitle = "  ")

    d <- plotprolactin(PRLpit, PRLpit$counts, "PRL", "pituitary") + 
      theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank(),
            axis.title.y = element_text(face = "italic"))

    e <- plotprolactin(prolactin, prolactin$plasma_conc, "prolactin (ng/mL)", "blood") +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) 

    f <- plotprolactin(PRLgon, PRLgon$counts, "PRL", "gonads") + 
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y = element_text(face = "italic"))


    plot_grid(a,b,c,d,e,f, labels = "auto", ncol = 2, rel_heights = c(1,1,1.2))

![](../figures/fig3-1.png)
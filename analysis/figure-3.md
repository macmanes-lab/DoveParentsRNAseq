Figure 3
========

    library(tidyverse)

    ## ── Attaching packages ────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────────── tidyverse_conflicts() ──
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
    charpca <- subsetmakepca("pituitary", charlevels, sexlevels)    
    charfviz <- makefvizdf("pituitary", charlevels, sexlevels)

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
      dplyr::select(sex, tissue, treatment, gene, samples, counts) %>%
      drop_na()
    PRLvsd$treatment <- factor(PRLvsd$treatment, levels = alllevels) 

    PRLhyp <- PRLvsd %>% filter(tissue == "hypothalamus") 
    PRLpit <- PRLvsd %>% filter(tissue == "pituitary")
    PRLgon <- PRLvsd %>% filter(tissue == "gonad")

    a <- plotcolorfulpcs(charpca,charpca$treatment, allcolors) + labs(subtitle = " ") +
      theme(legend.position = "none", 
            legend.direction = "horizontal", 
            legend.key.size = unit(0.5, 'lines')) + 
      guides(color = FALSE) +
      labs(subtitle = "pituitary")   

    b <- plotfriz(charfviz) + labs(subtitle = "pituitary") +
      #ylim(-250000,500000) + xlim(-250000,500000) + 
      theme(axis.text = element_blank())

    c <- plotprolactin(PRLpit, PRLpit$counts, "PRL", "pituitary") + 
      theme(legend.position = "none", 
             axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y = element_text(face = "italic"))

    d <- plotprolactin(prolactin, prolactin$plasma_conc, "prolactin (ng/mL)", "blood") +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) 

    abcd <- plot_grid(a,b,c,d, labels = c("b", "c", "d", "e"), ncol = 2, label_size = 8)

    expdesign <- png::readPNG("../figures/images/parentalstages_fig1a.png")
    expdesign <- ggdraw() +  draw_image(expdesign, scale = 1)

    abcde <- plot_grid(expdesign, abcd, nrow = 2, labels = c("a", "b"), label_size = 8, rel_heights = c(0.3,1))
    abcde

![](../figures/fig3-1.png)

    PRLvsd2 <- PRLvsd %>%
      filter(tissue == "pituitary")
    head(PRLvsd2)

    ## # A tibble: 6 x 6
    ##   sex    tissue    treatment gene  samples                           counts
    ##   <chr>  <chr>     <fct>     <chr> <chr>                              <dbl>
    ## 1 female pituitary control   PRL   L.G118_female_pituitary_control.…   17.9
    ## 2 female pituitary control   PRL   R.G106_female_pituitary_control     17.0
    ## 3 female pituitary control   PRL   R.R20_female_pituitary_control      18.6
    ## 4 female pituitary control   PRL   R.R9_female_pituitary_control.NY…   16.8
    ## 5 female pituitary control   PRL   R.W44_female_pituitary_control.N…   18.6
    ## 6 female pituitary inc.d9    PRL   blk.s061.pu.y_female_pituitary_i…   17.6

    PRLvsd2 %>%
      group_by(sex) %>%
      summarize(median = median(counts))

    ## # A tibble: 2 x 2
    ##   sex    median
    ##   <chr>   <dbl>
    ## 1 female   18.1
    ## 2 male     17.9

    a <- ggplot(PRLvsd2, aes(x = sex, y = counts, fill = sex))  +
      geom_boxplot() +
      labs(x = "Sex", y = "Prolactin expression")+
      theme_B3() + theme(legend.position = "none") +
      scale_fill_manual(values = sexcolors) +
      geom_hline(yintercept=18, linetype="dashed", color = "black") 


    PRLvsd3 <- PRLvsd2 %>%
      mutate(hiloPRL = ifelse(counts >= 18, "hi", "lo"))  %>%
      drop_na()
    PRLvsd3$hiloPRL <- factor(PRLvsd3$hiloPRL, levels = c("lo", "hi"))


    b <- ggplot(PRLvsd3, aes(x = sex, y = counts, fill = sex, color = hiloPRL))  +
      geom_boxplot() +
      labs(x = "Prolactin, binned into `lo` and `hi`", y = "Prolactin expression")+
      theme_B3() + theme(legend.position = "none") +
      scale_fill_manual(values = allcolors) +
      scale_color_manual(values = allcolors)

    plot_grid(a,b, nrow = 2, rel_widths = c(1,1))

![](../figures/determinePRLhiglo-1.png)

    PRLvsd3 %>%
      group_by(sex, tissue, hiloPRL) %>%
      summarize(n = n())

    ## # A tibble: 4 x 4
    ## # Groups:   sex, tissue [2]
    ##   sex    tissue    hiloPRL     n
    ##   <chr>  <chr>     <fct>   <int>
    ## 1 female pituitary lo         47
    ## 2 female pituitary hi         49
    ## 3 male   pituitary lo         51
    ## 4 male   pituitary hi         46

    write.csv(PRLvsd3,"../results/PRLvsd.csv", row.names = F)

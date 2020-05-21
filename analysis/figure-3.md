Figure 4: All things prolactin
==============================

    library(tidyverse)

    ## ── Attaching packages ────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

    library(corrr)
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

    ## The following object is masked from 'package:ggimage':
    ## 
    ##     theme_transparent

    ## The following object is masked from 'package:cowplot':
    ## 
    ##     get_legend

    library(ggrepel)
    library(ggsignif)
    library(forcats)

    source("../R/themes.R") 
    source("../R/functions.R")

    knitr::opts_chunk$set(fig.path = '../figures/',message=F, warning=FALSE)

import limma counts and sample info
-----------------------------------

    countData <- read_csv("../results/01_limma.csv") %>%
      column_to_rownames(var = "X1")
    countData <- as.data.frame(t(countData))
    head(countData[1:3])

    ##                                            A2ML1     A2ML2      A2ML3
    ## L.Blu13_male_gonad_control.NYNO         42.65677 4.4397552  211.96191
    ## L.Blu13_male_hypothalamus_control.NYNO 201.26331 4.8633048 6919.87335
    ## L.Blu13_male_pituitary_control.NYNO    161.22614 0.3158851  212.10845
    ## L.G107_male_gonad_control               43.22441 2.0404458  203.74048
    ## L.G107_male_hypothalamus_control       382.33084 4.9190817 9531.52550
    ## L.G107_male_pituitary_control           85.34910 0.3577761   69.02124

    colData <- read_csv("../metadata/04_colData.csv") %>%
      mutate(treatment = factor(treatment, levels = alllevels),
             tissue = factor(tissue, levels = tissuelevels),
             temprowname = V1) %>% 
      column_to_rownames(var = "temprowname") 
    head(colData)

    ##                                                                            V1
    ## L.Blu13_male_gonad_control.NYNO               L.Blu13_male_gonad_control.NYNO
    ## L.Blu13_male_hypothalamus_control.NYNO L.Blu13_male_hypothalamus_control.NYNO
    ## L.Blu13_male_pituitary_control.NYNO       L.Blu13_male_pituitary_control.NYNO
    ## L.G107_male_gonad_control                           L.G107_male_gonad_control
    ## L.G107_male_hypothalamus_control             L.G107_male_hypothalamus_control
    ## L.G107_male_pituitary_control                   L.G107_male_pituitary_control
    ##                                           bird  sex       tissue treatment
    ## L.Blu13_male_gonad_control.NYNO        L.Blu13 male       gonads   control
    ## L.Blu13_male_hypothalamus_control.NYNO L.Blu13 male hypothalamus   control
    ## L.Blu13_male_pituitary_control.NYNO    L.Blu13 male    pituitary   control
    ## L.G107_male_gonad_control               L.G107 male       gonads   control
    ## L.G107_male_hypothalamus_control        L.G107 male hypothalamus   control
    ## L.G107_male_pituitary_control           L.G107 male    pituitary   control
    ##                                                            group
    ## L.Blu13_male_gonad_control.NYNO              male.gonads.control
    ## L.Blu13_male_hypothalamus_control.NYNO male.hypothalamus.control
    ## L.Blu13_male_pituitary_control.NYNO       male.pituitary.control
    ## L.G107_male_gonad_control                    male.gonads.control
    ## L.G107_male_hypothalamus_control       male.hypothalamus.control
    ## L.G107_male_pituitary_control             male.pituitary.control
    ##                                                  study         sextissue
    ## L.Blu13_male_gonad_control.NYNO        charcterization       male_gonads
    ## L.Blu13_male_hypothalamus_control.NYNO charcterization male_hypothalamus
    ## L.Blu13_male_pituitary_control.NYNO    charcterization    male_pituitary
    ## L.G107_male_gonad_control              charcterization       male_gonads
    ## L.G107_male_hypothalamus_control       charcterization male_hypothalamus
    ## L.G107_male_pituitary_control          charcterization    male_pituitary
    ##                                                                       samples
    ## L.Blu13_male_gonad_control.NYNO               L.Blu13_male_gonad_control.NYNO
    ## L.Blu13_male_hypothalamus_control.NYNO L.Blu13_male_hypothalamus_control.NYNO
    ## L.Blu13_male_pituitary_control.NYNO       L.Blu13_male_pituitary_control.NYNO
    ## L.G107_male_gonad_control                           L.G107_male_gonad_control
    ## L.G107_male_hypothalamus_control             L.G107_male_hypothalamus_control
    ## L.G107_male_pituitary_control                   L.G107_male_pituitary_control
    ##                                        hiloPRL external
    ## L.Blu13_male_gonad_control.NYNO             lo  control
    ## L.Blu13_male_hypothalamus_control.NYNO      lo  control
    ## L.Blu13_male_pituitary_control.NYNO         lo  control
    ## L.G107_male_gonad_control                   lo  control
    ## L.G107_male_hypothalamus_control            lo  control
    ## L.G107_male_pituitary_control               lo  control

    # check ready for analysis
    # row.names(countData) == row.names(colData)
    head(row.names(countData) == row.names(colData))

    ## [1] TRUE TRUE TRUE TRUE TRUE TRUE

PCA data
--------

    # pca relies on count data from limma

    pca1f <- subsetmakepca("pituitary", charlevelsnocontrol, "female")  
    pca2f <- makefvizdf("pituitary", charlevelsnocontrol, "female") 

    pca1m <- subsetmakepca("pituitary", charlevelsnocontrol, "male")    
    pca2m <- makefvizdf("pituitary", charlevelsnocontrol, "male")

    pca1 <- subsetmakepca("pituitary", charlevelsnocontrol, sexlevels)  
    pca2 <- makefvizdf("pituitary", charlevelsnocontrol, sexlevels) 

candidate genes, hypotheses genes, and data-driven genes
--------------------------------------------------------

    # candidate counts for boxplots
    candidatevsd <- read_csv("../results/06_candidatevsd.csv") %>%
      mutate(treatment = factor(treatment)) %>%
      mutate(tissue = factor(tissue, levels = tissuelevel),
            treatment = factor(treatment, levels = charlevelsnocontrol)) %>%
      drop_na()

    ## prlpit
    PRLpit <- candidatevsd %>% 
      filter(gene == "PRL", tissue == "pituitary")

    # list of pvalue and DEGs
    allDEG <- read_csv("../results/06_allDEG.csv")  %>%
      mutate(direction = factor(direction, levels = hypothesislevels))

    # candidate genes from GO and literature
    parentalcaregenes <- read_csv("../metadata/03_parentalcaregenes.csv") %>% 
      distinct(literature)  %>% drop_na() %>% pull(literature)

    ## genes WGCNA prl module
    WGCNAgenes <- read_csv("../results/05_PRLmodule.csv") %>% pull(x)

    # DEGs infor for PRL and external
    hilodf <- allDEG %>%
      filter(comparison == "lo_hi")  %>%
      arrange(padj) 

    eggchickdf <- allDEG %>%
      filter(!comparison %in% c("lo_hi", "nest_eggs", "nest_chicks"))  %>%
      arrange(padj)  
      
    eggchickgenes <- eggchickdf %>%  top_n(20) %>% pull(gene)
    hilogenes <- hilodf %>%  top_n(20) %>% pull(gene)

    ## for correlations
    pithilo <- makecorrdf("female", "pituitary", hilogenes)

    ## [1] "The top two correlations in this tissue and sex are:"
    ## # A tibble: 2 x 5
    ##   tissue    sex    gene1 gene2  corr
    ##   <chr>     <chr>  <chr> <chr> <dbl>
    ## 1 pituitary female CKAP2 KPNA2 0.975
    ## 2 pituitary female CDK1  KPNA2 0.971

    pitWGCNA <- makecorrdf("female", "pituitary", WGCNAgenes[1:20])

    ## [1] "The top two correlations in this tissue and sex are:"
    ## # A tibble: 2 x 5
    ##   tissue    sex    gene1   gene2   corr
    ##   <chr>     <chr>  <chr>   <chr>  <dbl>
    ## 1 pituitary female CRELD2  FKBP11 0.912
    ## 2 pituitary female CREB3L1 FKBP11 0.910

    ## for scatter correlation

    candidatevsdwide <- candidatevsd  %>%
        pivot_wider(names_from = gene, values_from = counts) 
    pitwide <- subsetcandidatevsdwide(sexlevels, "pituitary")

    a1 <- plotpc12(pca1, pca2, pca1$treatment, allcolors, "Pituitary gene expression", NULL)   
    a2 <- plotprolactin(PRLpit, PRLpit$counts, "PRL", NULL) + 
      theme(axis.text.x = element_text(angle = 25))  +
      facet_wrap(~sex)
    a <- plot_grid(a1,a2, labels = c("A", "B"), label_size = 8, rel_widths = c(1,1))

    b1 <- png::readPNG("../figures/images/fig_fig3b.png")
    b1 <- ggdraw() +  draw_image(b1, scale = 1)
    b2 <- makenewbargraph("pituitary", "female","eggs_chicks", 0, 4000) + labs(subtitle = "female")  + labs(title = " ") +
      scale_x_discrete(labels = "eggs vs. chicks")
    b3 <- makenewbargraph("pituitary", "male", "eggs_chicks", 0, 4000)  + labs(subtitle = "male")  + labs(title = " ") +
      theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
            axis.text.y = element_blank(), axis.line.y = element_blank()) +
      scale_x_discrete(labels = "eggs vs. chicks")
    b4 <- makenewbargraph("pituitary", "male", "lo_hi", 0, 4000)  +
      scale_x_discrete(labels = "lo vs. hi PRL")
    b5 <- makenewbargraph("pituitary", "female", "lo_hi", 0, 4000)   + 
      theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
            axis.text.y = element_blank(), axis.line.y = element_blank()) +
      scale_x_discrete(labels = "lo vs. hi PRL")
    b25 <- plot_grid(b2,b3,b4,b5, nrow = 2, rel_widths = c(1.2,1), rel_heights = c(1.2,1))

    b6 <-  plottopDEGs(eggchickdf, "female",  "LFC eggs vs. chicks", "Top 10 female DEGs") + labs(title = " ")
    b7 <-  plottopDEGs(eggchickdf, "male",  "LFC eggs vs. chicks", "Top 10 male DEGs") + labs(title = " ")
    b8 <-  plottopDEGs(hilodf, "female",  "LFC lo vs. hi PRL", NULL)
    b9 <-  plottopDEGs(hilodf, "male",  "LFC lo vs. hi PRL", NULL)
    b15 <- plot_grid(b1,b25,nrow = 1, rel_widths = c(1.2,1))
    b69 <- plot_grid(b6,b7,b8,b9, align = "v" , rel_heights = c(1.2,1))
    b <- plot_grid(b15,b69, labels = c("C"), label_size = 8, rel_widths = c(1.2,1))

    c1 <- plotcorrplot(pithilo, NULL) + theme(legend.position = "none") + labs(subtitle = "Top 20 lo vs. hi PRL DEGs") 
    c2 <- plotcorrplot(pitWGCNA, NULL) + theme(legend.position = "none") + labs(subtitle = "20 genes from WGCNA module with PRL") + theme(legend.position = "right")
    c3 <- candidateboxplot("pituitary", c("KPNA2"), sexlevels)
    c4 <- candidateboxplot("pituitary", c("CDK1"), sexlevels) + 
      theme(axis.text.x = element_text(angle = 45, hjust =1), strip.text = element_blank())  + labs(x = " ")
    c34 <- plot_grid(c3,c4, nrow = 2, rel_heights = c(1,1.2))

    c <- plot_grid(c1,c2,c34, nrow = 1, labels = c("D", "", "E"), label_size = 8, rel_widths = c(1.1,1.4,1))

    fig3 <- plot_grid(a, b, c, nrow = 3, rel_heights = c(1.2,2,1.8))
    fig3

![](../figures/fig3-1.png)

    pdf(file="../figures/fig3-1.pdf", width=7.25, height=7.25)
    plot(fig3)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    write.csv(PRLpit,"../results/PRLvsd.csv", row.names = F)

Figure 3: All things prolactin
==============================

    library(tidyverse)

    ## ── Attaching packages ──────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────── tidyverse_conflicts() ──
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
    source("../R/genelists.R")

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

    ## Parsed with column specification:
    ## cols(
    ##   x = col_character()
    ## )

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
      filter(treatment %in% charlevels) %>%
      mutate(tissue = factor(tissue, levels = tissuelevels),
             treatment = factor(treatment, levels = charlevels),
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
    pca1f <- subsetmakepca("pituitary", charlevels, "female")   
    pca2f <- makefvizdf("pituitary", charlevels, "female")  

    pca1m <- subsetmakepca("pituitary", charlevels, "male") 
    pca2m <- makefvizdf("pituitary", charlevels, "male")

    pca1 <- subsetmakepca("pituitary", charlevels, sexlevels)   
    pca2 <- makefvizdf("pituitary", charlevels, sexlevels)  

candidate genes, hypotheses genes, and data-driven genes
--------------------------------------------------------

    # all genes lists in "./R/genelists.R"

    # candidate counts for boxplots
    candidatevsd <- read_csv("../results/06_candidatevsd.csv") %>%
      filter(treatment %in% charlevels) %>%
      mutate(treatment = factor(treatment)) %>%
      mutate(tissue = factor(tissue, levels = tissuelevel),
            treatment = factor(treatment, levels = charlevels)) %>%
      drop_na()
    head(candidatevsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue      treatment gene  samples                         counts
    ##   <chr>  <fct>       <fct>     <chr> <chr>                            <dbl>
    ## 1 female hypothalam… control   ABCA4 L.G118_female_hypothalamus_con…   6.49
    ## 2 female hypothalam… control   ABCA4 R.G106_female_hypothalamus_con…   6.47
    ## 3 female hypothalam… control   ABCA4 R.R20_female_hypothalamus_cont…   6.31
    ## 4 female hypothalam… control   ABCA4 R.R9_female_hypothalamus_contr…   6.30
    ## 5 female hypothalam… control   ABCA4 R.W44_female_hypothalamus_cont…   6.31
    ## 6 female hypothalam… inc.d9    ABCA4 blk.s061.pu.y_female_hypothala…   6.49

    ## prlpit
    PRLpit <- candidatevsd %>% 
      filter(gene == "PRL", tissue == "pituitary")

    # list of pvalue and DEGs
    allDEG <- read_csv("../results/06_allDEG.csv")  %>%
      mutate(direction = factor(direction, levels = hypothesislevels))

    # DEGs infor for PRL and external
    hilodf <- allDEG %>%
      filter(comparison == "lo_hi")  %>%
      arrange(padj) 
    hilogenes <- hilodf %>%  top_n(20) %>% pull(gene)

    eggchickdf <- allDEG %>%
      filter(!comparison %in% c("lo_hi", "nest_eggs", "nest_chicks"))  %>%
      arrange(padj)  
    eggchickgenes <- eggchickdf %>%  top_n(20) %>% pull(gene)


    ## data-driven PRL genes
    pithilo <- makecorrdf("female", "pituitary", hilogenes)
    pitWGCNA <- makecorrdf("female", "pituitary", WGCNAgenes[32:52])


    ## https://www.gynecologiconcology-online.net/article/S0090-8258(19)30069-1/fulltext
    cancergenes <- makecorrdf("female", "pituitary", c(shaidgenes,suszynskaagenes, "PRL"))


    ## for scatter correlation (not currently usesed)
    candidatevsdwide <- candidatevsd  %>%
        pivot_wider(names_from = gene, values_from = counts) 
    pitwide <- subsetcandidatevsdwide(sexlevels, "pituitary")

    ## number of genes total vs plotted

    eggchickdf %>% filter(tissue == "pituitary") %>% nrow()

    ## [1] 508

    hilodf %>% filter(tissue == "pituitary") %>% nrow()

    ## [1] 10684

    length(WGCNAgenes)

    ## [1] 58

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

    c1 <- plotcorrplot(pithilo, NULL) + theme(legend.position = "none") + labs(subtitle = "Top 20 DEGs, low vs. high PRL") 
    c2 <- plotcorrplot(pitWGCNA, NULL) + theme(legend.position = "none") + labs(subtitle = "20 genes in WGCNA module with PRL") 
    c3 <- plotcorrplot(cancergenes, NULL) + theme(legend.position = "right") + labs(subtitle = "Breast, ovary, & piuitary cancer genes") 
    c <- plot_grid(c1,c2,c3, nrow = 1, labels = c("D"), label_size = 8, rel_widths = c(1,1,1.3))

    fig3 <- plot_grid(a, b, c, nrow = 3, rel_heights = c(1.2,2,1.8))
    fig3

![](../figures/fig3-1.png)

    write.csv(PRLpit,"../results/PRLvsd.csv", row.names = F)

    pdf(file="../figures/fig3-1.pdf", width=7.25, height=7.25)
    plot(fig3)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    ## PRL and PRLR in three tissues 

    candidatevsd <- candidatevsd %>%
      filter(treatment %in%  c("early", "extend", charlevels))

    p1 <- candidateboxplot("hypothalamus", c("PRL"), sexlevels) + labs(subtitle = "hypothalamus", title = "PRL")
    p2 <- candidateboxplot("pituitary", c("PRL"), sexlevels) + labs(subtitle = "pituitary")
    p3 <- candidateboxplot("gonad", c("PRL"), sexlevels) + labs(subtitle = "gonads") + theme(axis.text.x = element_text(angle = 45))

    p4 <- candidateboxplot("hypothalamus", c("PRLR"), sexlevels) + labs(subtitle = " ", title = " ") 
    p5 <- candidateboxplot("pituitary", c("PRLR"), sexlevels) + labs(subtitle = " ")
    p6 <- candidateboxplot("gonad", c("PRLR"), sexlevels) + theme(axis.text.x = element_text(angle = 45))+ labs(subtitle = " ")

    plot_grid(p1,p4,p2,p5,p3,p6, ncol = 2, rel_heights = c(1.2,1,1.2))

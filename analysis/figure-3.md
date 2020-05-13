Figure 4: All things prolactin
==============================

    library(tidyverse)

    ## ── Attaching packages ────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

    pca1f <- subsetmakepca("pituitary", alllevels, "female")    
    pca2f <- makefvizdf("pituitary", alllevels, "female")   

    pca1m <- subsetmakepca("pituitary", alllevels, "male")  
    pca2m <- makefvizdf("pituitary", alllevels, "male")

    pca1 <- subsetmakepca("pituitary", alllevels, sexlevels)    
    pca2 <- makefvizdf("pituitary", alllevels, sexlevels)   

candidate genes, hypotheses genes, and data-driven genes
--------------------------------------------------------

    # candidate counts for boxplots
    candidatevsd <- read_csv("../results/06_candidatevsd.csv") %>%
      mutate(treatment = factor(treatment)) %>%
      mutate(tissue = factor(tissue, levels = tissuelevel),
            treatment = factor(treatment, levels = alllevels)) %>%
      drop_na()

    ## prlpit
    PRLpit <- candidatevsd %>% 
      filter(gene == "PRL", tissue == "pituitary")

    # list of pvalue and DEGs
    allDEG <- read_csv("../results/06_allDEG.csv") 
    allDEG2 <- allDEG

    # candidate genes from GO and literature
    parentalcaregenes <- read_csv("../metadata/03_parentalcaregenes.csv") %>% pull(gene)

    ## genes WGCNA prl module
    WGCNAgenes <- read_csv("../results/05_PRLmodule.csv") %>% pull(x)

    # DEGs infor for PRL and external
    hilodf <- allDEG %>%
      filter(comparison == "lo_hi")  %>%
      arrange(padj) 

    eggchickdf <- allDEG %>%
      filter(comparison != "lo_hi")  %>%
      arrange(padj)  
      
    eggchickgenes <- eggchickdf %>%  top_n(20) %>% pull(gene)
    hilogenes <- hilodf %>%  top_n(20) %>% pull(gene)

    ## for correlations
    pithilo <- makecorrdf("female", "pituitary", hilogenes)

    ## # A tibble: 6 x 19
    ##   rowname   CD24   SCG3  OLFM1  ARAP2 LOC428479   STC1 LAPTM4B CENPI KNTC1
    ##   <chr>    <dbl>  <dbl>  <dbl>  <dbl>     <dbl>  <dbl>   <dbl> <dbl> <dbl>
    ## 1 CD24    NA      0.827  0.749  0.755     0.767  0.748   0.813 0.480 0.518
    ## 2 SCG3     0.827 NA      0.790  0.808     0.790  0.786   0.770 0.531 0.530
    ## 3 OLFM1    0.749  0.790 NA      0.753     0.781  0.811   0.846 0.520 0.589
    ## 4 ARAP2    0.755  0.808  0.753 NA         0.809  0.787   0.852 0.603 0.652
    ## 5 LOC428…  0.767  0.790  0.781  0.809    NA      0.887   0.875 0.617 0.649
    ## 6 STC1     0.748  0.786  0.811  0.787     0.887 NA       0.909 0.669 0.726
    ## # … with 9 more variables: UBE2C <dbl>, CDK1 <dbl>, LOC423793 <dbl>,
    ## #   BUB1 <dbl>, KPNA2 <dbl>, RACGAP1 <dbl>, CCNB3 <dbl>, CKAP2 <dbl>,
    ## #   RRM2 <dbl>

    piteggchick <- makecorrdf("female", "pituitary", eggchickgenes)

    ## # A tibble: 6 x 20
    ##   rowname TMEM8A    VTN  SARM1  DDX10  NPDC1  GREB1 SPTY2D1 TRIOBP NETO2
    ##   <chr>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>  <dbl> <dbl>
    ## 1 TMEM8A  NA      0.742  0.718  0.705  0.780  0.759   0.660  0.635 0.594
    ## 2 VTN      0.742 NA      0.831  0.660  0.596  0.729   0.696  0.628 0.579
    ## 3 SARM1    0.718  0.831 NA      0.580  0.569  0.780   0.617  0.595 0.699
    ## 4 DDX10    0.705  0.660  0.580 NA      0.679  0.544   0.722  0.579 0.563
    ## 5 NPDC1    0.780  0.596  0.569  0.679 NA      0.536   0.519  0.710 0.576
    ## 6 GREB1    0.759  0.729  0.780  0.544  0.536 NA       0.714  0.418 0.708
    ## # … with 10 more variables: MXRA7 <dbl>, PNPLA3 <dbl>, TH <dbl>,
    ## #   ABCA4 <dbl>, TMEM201 <dbl>, LOC107056420 <dbl>, ANGPTL4 <dbl>,
    ## #   FADS6 <dbl>, PCNA <dbl>, H2AFY2 <dbl>

    pitWGCNA <- makecorrdf("female", "pituitary", WGCNAgenes[1:20])

    ## # A tibble: 6 x 21
    ##   rowname  ACOT12 C2H8ORF46 CHRNB4  F13A1 FAM83G    CD24 FOSL2  AGR2
    ##   <chr>     <dbl>     <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl> <dbl>
    ## 1 ACOT12  NA          0.143  0.179  0.325  0.211  0.0955 0.265 0.195
    ## 2 C2H8OR…  0.143     NA      0.358  0.187  0.410  0.383  0.264 0.371
    ## 3 CHRNB4   0.179      0.358 NA      0.353  0.433  0.652  0.347 0.545
    ## 4 F13A1    0.325      0.187  0.353 NA      0.461  0.325  0.464 0.535
    ## 5 FAM83G   0.211      0.410  0.433  0.461 NA      0.534  0.541 0.468
    ## 6 CD24     0.0955     0.383  0.652  0.325  0.534 NA      0.380 0.460
    ## # … with 12 more variables: COL20A1 <dbl>, DUSP10 <dbl>, CDT1 <dbl>,
    ## #   CALN1 <dbl>, CDKN3 <dbl>, CREM <dbl>, ARAP2 <dbl>, FKBP11 <dbl>,
    ## #   FAM19A1 <dbl>, CRELD2 <dbl>, FGF6 <dbl>, CREB3L1 <dbl>

    parentalcare <- makecorrdf("female", "pituitary", parentalcaregenes[1:20])

    ## # A tibble: 6 x 20
    ##   rowname  AVPR1A   CRHR1      FOS    MBD2     ESR1   CRHR2    DRD1
    ##   <chr>     <dbl>   <dbl>    <dbl>   <dbl>    <dbl>   <dbl>   <dbl>
    ## 1 AVPR1A  NA       0.736   0.488    0.110   0.326    0.0841 0.0720 
    ## 2 CRHR1    0.736  NA       0.369    0.179   0.0909   0.200  0.109  
    ## 3 FOS      0.488   0.369  NA       -0.0201 -0.00459  0.147  0.0473 
    ## 4 MBD2     0.110   0.179  -0.0201  NA      -0.0590  -0.0434 0.00966
    ## 5 ESR1     0.326   0.0909 -0.00459 -0.0590 NA        0.0945 0.0972 
    ## 6 CRHR2    0.0841  0.200   0.147   -0.0434  0.0945  NA      0.0969 
    ## # … with 12 more variables: ESR2 <dbl>, ADRA2A <dbl>, DBH <dbl>,
    ## #   KALRN <dbl>, GNAQ <dbl>, DRD4 <dbl>, HTR2C <dbl>, BRINP1 <dbl>,
    ## #   CREBRF <dbl>, AVP <dbl>, COMT <dbl>, CRHBP <dbl>

    a1 <- plotpc12(pca1, pca2, pca1$treatment, allcolors, "Pituitary gene expression", " ")   
    a2 <- plotprolactin(PRLpit, PRLpit$counts, "PRL", " ") + 
      theme(axis.text.x = element_text(angle = 25))  +
      theme(legend.position = c(0.8,0.2)) + 
      guides(fill = FALSE) 
    a <- plot_grid(a1,a2, labels = c("A"), label_size = 8, rel_widths = c(1,1.1))

    b1 <- png::readPNG("../figures/images/fig_fig3b.png")
    b1 <- ggdraw() +  draw_image(b1, scale = 1)
    b2 <- makenewbargraph("pituitary", "female","eggs_chicks", 0, 400) + labs(subtitle = "female")  + labs(title = " ")  
    b3 <- makenewbargraph("pituitary", "male", "eggs_chicks", 0, 400)  + labs(subtitle = "male")  + labs(title = " ") +
      theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
            axis.text.y = element_blank(), axis.line.y = element_blank())
    b4 <- makenewbargraph("pituitary", "male", "lo_hi", 0, 4000) 
    b5 <- makenewbargraph("pituitary", "female", "lo_hi", 0, 4000)   + 
      theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
            axis.text.y = element_blank(), axis.line.y = element_blank())
    b25 <- plot_grid(b2,b3,b4,b5, nrow = 2, rel_widths = c(1.2,1), rel_heights = c(1.2,1))
    b6 <-  plottopDEGs(eggchickdf, "female", "#3A80B9", "LFC eggs vs. chicks", "Top 10 female DEGs") + labs(title = " ")
    b7 <-  plottopDEGs(hilodf, "male", "#3A80B9", "LFC eggs vs. chicks", "Top 10 male DEGs") + labs(title = " ")
    b8 <-  plottopDEGs(hilodf, "female", "#19757A", "LFC lo vs. hi PRL", NULL)
    b9 <-  plottopDEGs(hilodf, "male", "#19757A", "LFC lo vs. hi PRL", NULL)
    b15 <- plot_grid(b1,b25,nrow = 1, rel_widths = c(1.2,1))
    b69 <- plot_grid(b6,b7,b8,b9, align = "v" , rel_heights = c(1.2,1))
    b <- plot_grid(b15,b69, labels = c("B"), label_size = 8, rel_widths = c(1.2,1))

    c1 <- plotcorrplot(piteggchick, NULL) + theme(legend.position = "none") + labs(subtitle = "Top 20 eggs vs. chicks DEGs")
    c2 <- plotcorrplot(pithilo, NULL) + theme(legend.position = "none") + labs(subtitle = "Top 20 lo vs. PRL DEGs") 
    d <- plotcorrplot(pitWGCNA, NULL) + theme(legend.position = "none") + labs(subtitle = "20 genes from WGCNA module with PRL") + theme(legend.position = "right")
    cd <- plot_grid(c1,c2,d, nrow = 1, align = "h", labels = c("C", "", "D"), label_size = 8, rel_widths = c(1,1,1.3))


    fig3 <- plot_grid(a, b, cd, nrow = 3, rel_heights = c(1.2,2,1.8))
    fig3

![](../figures/fig3-1.png)

    pdf(file="../figures/fig3-1.pdf", width=7.25, height=7.25)
    plot(fig3)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    write.csv(PRLpit,"../results/PRLvsd.csv", row.names = F)

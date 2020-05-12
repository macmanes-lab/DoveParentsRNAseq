Figure 4: All things prolactin
==============================

    library(tidyverse)

    ## ── Attaching packages ──────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

    source("../R/themes.R") 
    source("../R/functions.R")

    source("../R/datawrangling.R")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_double(),
    ##   gene = col_character(),
    ##   geneid = col_double(),
    ##   NCBI = col_character()
    ## )

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_character(),
    ##   logpadj = col_character(),
    ##   lfc = col_character(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_character(),
    ##   logpadj = col_character(),
    ##   lfc = col_character(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_character(),
    ##   logpadj = col_character(),
    ##   lfc = col_character(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
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

    colData <- read_csv("../metadata/00_colData.csv") %>%
      mutate(treatment = factor(treatment, levels = alllevels),
             tissue = factor(tissue, levels = tissuelevels)) %>% 
      column_to_rownames(var = "X1") 
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
    ##                                                  study
    ## L.Blu13_male_gonad_control.NYNO        charcterization
    ## L.Blu13_male_hypothalamus_control.NYNO charcterization
    ## L.Blu13_male_pituitary_control.NYNO    charcterization
    ## L.G107_male_gonad_control              charcterization
    ## L.G107_male_hypothalamus_control       charcterization
    ## L.G107_male_pituitary_control          charcterization

    # check ready for analysis
    # row.names(countData) == row.names(colData)
    head(row.names(countData) == row.names(colData))

    ## [1] TRUE TRUE TRUE TRUE TRUE TRUE

PCA data
--------

    # pca
    pca1f <- subsetmakepca("pituitary", charlevels, "female")   
    pca2f <- makefvizdf("pituitary", charlevels, "female")  

    pca1m <- subsetmakepca("pituitary", charlevels, "male") 
    pca2m <- makefvizdf("pituitary", charlevels, "male")

    pca1 <- subsetmakepca("pituitary", charlevels, sexlevels)   
    pca2 <- makefvizdf("pituitary", charlevels, sexlevels)  

Hi Lo PRL and WGCNA PRL
=======================

    # `candidatevsd` loaded with `source("../R/wrangledata.R")`

    PRLpit <- candidatevsd %>% filter(tissue == "pituitary", gene == "PRL") %>%
      mutate(hiloPRL = ifelse(counts >= 18, "hi", "lo"))  %>%
      drop_na()
    PRLpit$hiloPRL <- factor(PRLpit$hiloPRL, levels = c("lo", "hi"))

    PRLpit %>%
      group_by(sex) %>%
      summarize(median = median(counts))

    ## # A tibble: 2 x 2
    ##   sex    median
    ##   <chr>   <dbl>
    ## 1 female   18.9
    ## 2 male     18.5

    PRLpit %>%
      group_by(sex, tissue, hiloPRL) %>%
      summarize(n = n())

    ## # A tibble: 4 x 4
    ## # Groups:   sex, tissue [2]
    ##   sex    tissue    hiloPRL     n
    ##   <chr>  <chr>     <fct>   <int>
    ## 1 female pituitary lo         54
    ## 2 female pituitary hi        111
    ## 3 male   pituitary lo         70
    ## 4 male   pituitary hi         95

    PRLpitF <- PRLpit %>% filter( sex == "female" )
    PRLpitM <- PRLpit %>% filter( sex == "male")

### WGCNA prolactin module

    ## genes WGCNA prl module

    WGCNAgenes <- read_csv("../results/PRLmodule.csv") %>% pull(x)

    WGCNAvsd <- getcandidatevsd(WGCNAgenes, "pituitary", sexlevels) 

Internal versus external hyotheses
----------------------------------

    #deg wrangling moved to datawrangling.R

    hypothesisDEGs <- c(hypothesisDEGs, WGCNAgenes) %>% unique()

    candidatevsd <- getcandidatevsd(hypothesisDEGs, "pituitary", sexlevels) 

    #candidateboxplot("pituitary", c("PNOC"), "female") + labs(x = NULL )  + theme( strip.text = element_blank())
    #candidateboxplot("pituitary", c("KPNA2"), "male") + labs(x = NULL )  + theme( strip.text = element_blank())


    pithilo <- makecorrdf("female", "pituitary", hilogenes)

    ## # A tibble: 6 x 17
    ##   rowname  OLFM1  STT3B  ARAP2    MLN LOC101749834    PRL  FGF6 LAPTM4B
    ##   <chr>    <dbl>  <dbl>  <dbl>  <dbl>        <dbl>  <dbl> <dbl>   <dbl>
    ## 1 OLFM1   NA      0.832  0.737  0.741        0.709  0.858 0.676   0.847
    ## 2 STT3B    0.832 NA      0.769  0.766        0.741  0.862 0.759   0.865
    ## 3 ARAP2    0.737  0.769 NA      0.721        0.653  0.815 0.681   0.838
    ## 4 MLN      0.741  0.766  0.721 NA            0.737  0.759 0.778   0.839
    ## 5 LOC101…  0.709  0.741  0.653  0.737       NA      0.659 0.755   0.782
    ## 6 PRL      0.858  0.862  0.815  0.759        0.659 NA     0.683   0.899
    ## # … with 8 more variables: CENPI <dbl>, CDK1 <dbl>, KPNA2 <dbl>,
    ## #   BUB1 <dbl>, LOC423793 <dbl>, CKAP2 <dbl>, RRM2 <dbl>, RACGAP1 <dbl>

    piteggchick <- makecorrdf("female", "pituitary", eggchickgenes)

    ## # A tibble: 6 x 20
    ##   rowname  SF3A1 KIAA1522 ZC3H18  WIPF3 LOC101748741  UBAC1 SEBOX   VTN
    ##   <chr>    <dbl>    <dbl>  <dbl>  <dbl>        <dbl>  <dbl> <dbl> <dbl>
    ## 1 SF3A1   NA        0.686  0.764  0.508        0.709  0.572 0.261 0.182
    ## 2 KIAA15…  0.686   NA      0.668  0.555        0.541  0.498 0.522 0.420
    ## 3 ZC3H18   0.764    0.668 NA      0.567        0.497  0.405 0.412 0.285
    ## 4 WIPF3    0.508    0.555  0.567 NA            0.476  0.380 0.557 0.450
    ## 5 LOC101…  0.709    0.541  0.497  0.476       NA      0.519 0.224 0.209
    ## 6 UBAC1    0.572    0.498  0.405  0.380        0.519 NA     0.226 0.257
    ## # … with 11 more variables: CTNND2 <dbl>, FBXO21 <dbl>, NKRF <dbl>,
    ## #   NEMF <dbl>, UIMC1 <dbl>, LOC422319 <dbl>, TULP1 <dbl>, ZBTB16 <dbl>,
    ## #   INSM1 <dbl>, FKBP5 <dbl>, MYBPC3 <dbl>

    pitWGCNA <- makecorrdf("female", "pituitary", WGCNAgenes[1:20])

    ## # A tibble: 6 x 21
    ##   rowname  ACOT12 C2H8ORF46 CHRNB4  F13A1 FAM83G    CD24  AGR2 FOSL2 DUSP10
    ##   <chr>     <dbl>     <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl> <dbl>  <dbl>
    ## 1 ACOT12  NA          0.111  0.156  0.341  0.218  0.0815 0.212 0.238  0.271
    ## 2 C2H8OR…  0.111     NA      0.370  0.103  0.427  0.405  0.363 0.266  0.276
    ## 3 CHRNB4   0.156      0.370 NA      0.250  0.422  0.638  0.530 0.351  0.420
    ## 4 F13A1    0.341      0.103  0.250 NA      0.413  0.286  0.487 0.409  0.434
    ## 5 FAM83G   0.218      0.427  0.422  0.413 NA      0.522  0.490 0.528  0.521
    ## 6 CD24     0.0815     0.405  0.638  0.286  0.522 NA      0.440 0.393  0.524
    ## # … with 11 more variables: COL20A1 <dbl>, CDT1 <dbl>, CREM <dbl>,
    ## #   CDKN3 <dbl>, CALN1 <dbl>, ARAP2 <dbl>, FAM19A1 <dbl>, FKBP11 <dbl>,
    ## #   CRELD2 <dbl>, FGF6 <dbl>, CREB3L1 <dbl>

    a1 <- plotpc12(pca1, pca2, pca1$treatment, allcolors, "Pituitary gene expression", " ")   
    a2 <- plotprolactin(PRLpit, PRLpit$counts, "PRL", " ") + theme(axis.text.x = element_text())  +
      theme(legend.position = c(0.8,0.2)) + guides(fill = FALSE)

    a <- plot_grid(a1,a2, labels = c("A"), label_size = 8)


    b1 <- png::readPNG("../figures/images/fig_fig3b.png")
    b1 <- ggdraw() +  draw_image(b1, scale = 1)


    b2 <- makenewbargraph("pituitary", "female","eggs vs. chicks", 0, 2200) + labs(subtitle = "female")  + labs(title = " ")  
    b3 <- makenewbargraph("pituitary", "male", "eggs vs. chicks", 0, 2200)  + labs(subtitle = "male")  + labs(title = " ") +
      theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
            axis.text.y = element_blank(), axis.line.y = element_blank())

    b4 <- makenewbargraph("pituitary", "male", "lo vs. hi PRL   ", 0, 2200) 
    b5 <- makenewbargraph("pituitary", "female", "lo vs. hi PRL   ", 0, 2200)   + 
      theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), 
            axis.text.y = element_blank(), axis.line.y = element_blank())

    b25 <- plot_grid(b2,b3,b4,b5, nrow = 2, rel_widths = c(1.2,1), rel_heights = c(1.2,1))

    b6 <-  plottopDEGs(envDEGs, "female", "#3A80B9", "LFC eggs vs. chicks", "Top 10 female DEGs") + labs(title = " ")
    b7 <-  plottopDEGs(PRLDEGs, "male", "#3A80B9", "LFC eggs vs. chicks", "Top 10 male DEGs") + labs(title = " ")
    b8 <-  plottopDEGs(PRLDEGs, "female", "#19757A", "LFC lo vs. hi PRL", NULL)
    b9 <-  plottopDEGs(PRLDEGs, "male", "#19757A", "LFC lo vs. hi PRL", NULL)


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

Figure 4: All things prolactin
==============================

    library(tidyverse)

    ## ── Attaching packages ─────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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
    source("../R/wrangledata.R")

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
    ##   X1 = col_double(),
    ##   Name = col_character(),
    ##   geneid = col_double(),
    ##   entrezid = col_character()
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
    ## 1 female   18.1
    ## 2 male     17.9

    PRLpit %>%
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
    ##   rowname  OLFM1  STT3B    PRL  ARAP2 LOC101749834    MLN LAPTM4B  FGF6
    ##   <chr>    <dbl>  <dbl>  <dbl>  <dbl>        <dbl>  <dbl>   <dbl> <dbl>
    ## 1 OLFM1   NA      0.831  0.839  0.722        0.752  0.782   0.827 0.692
    ## 2 STT3B    0.831 NA      0.904  0.809        0.752  0.807   0.863 0.749
    ## 3 PRL      0.839  0.904 NA      0.817        0.756  0.839   0.928 0.720
    ## 4 ARAP2    0.722  0.809  0.817 NA            0.710  0.741   0.840 0.724
    ## 5 LOC101…  0.752  0.752  0.756  0.710       NA      0.768   0.798 0.756
    ## 6 MLN      0.782  0.807  0.839  0.741        0.768 NA       0.868 0.794
    ## # … with 8 more variables: CENPI <dbl>, CDK1 <dbl>, BUB1 <dbl>,
    ## #   LOC423793 <dbl>, KPNA2 <dbl>, RACGAP1 <dbl>, CKAP2 <dbl>, RRM2 <dbl>

    piteggchick <- makecorrdf("female", "pituitary", eggchickgenes)

    ## # A tibble: 6 x 20
    ##   rowname  SF3A1 LOC101748741 FBXO21  UBAC1   NEMF KIAA1522  NKRF ZC3H18
    ##   <chr>    <dbl>        <dbl>  <dbl>  <dbl>  <dbl>    <dbl> <dbl>  <dbl>
    ## 1 SF3A1   NA            0.711  0.556  0.610  0.684    0.725 0.634  0.780
    ## 2 LOC101…  0.711       NA      0.575  0.589  0.482    0.543 0.524  0.504
    ## 3 FBXO21   0.556        0.575 NA      0.605  0.467    0.380 0.439  0.313
    ## 4 UBAC1    0.610        0.589  0.605 NA      0.373    0.497 0.430  0.373
    ## 5 NEMF     0.684        0.482  0.467  0.373 NA        0.512 0.424  0.547
    ## 6 KIAA15…  0.725        0.543  0.380  0.497  0.512   NA     0.426  0.708
    ## # … with 11 more variables: UIMC1 <dbl>, ZBTB16 <dbl>, CTNND2 <dbl>,
    ## #   TULP1 <dbl>, WIPF3 <dbl>, INSM1 <dbl>, LOC422319 <dbl>, VTN <dbl>,
    ## #   SEBOX <dbl>, FKBP5 <dbl>, MYBPC3 <dbl>

    pitWGCNA <- makecorrdf("female", "pituitary", WGCNAgenes[1:20])

    ## # A tibble: 6 x 21
    ##   rowname ACOT12 C2H8ORF46 CHRNB4  F13A1 COL20A1 FAM83G  CD24 FOSL2  AGR2
    ##   <chr>    <dbl>     <dbl>  <dbl>  <dbl>   <dbl>  <dbl> <dbl> <dbl> <dbl>
    ## 1 ACOT12  NA         0.242  0.248  0.299   0.241  0.201 0.133 0.251 0.263
    ## 2 C2H8OR…  0.242    NA      0.371  0.158   0.388  0.408 0.360 0.370 0.408
    ## 3 CHRNB4   0.248     0.371 NA      0.336   0.325  0.413 0.618 0.265 0.596
    ## 4 F13A1    0.299     0.158  0.336 NA       0.361  0.433 0.366 0.484 0.528
    ## 5 COL20A1  0.241     0.388  0.325  0.361  NA      0.542 0.652 0.413 0.461
    ## 6 FAM83G   0.201     0.408  0.413  0.433   0.542 NA     0.513 0.588 0.480
    ## # … with 11 more variables: CALN1 <dbl>, DUSP10 <dbl>, CDT1 <dbl>,
    ## #   CREM <dbl>, ARAP2 <dbl>, CDKN3 <dbl>, FAM19A1 <dbl>, FKBP11 <dbl>,
    ## #   FGF6 <dbl>, CRELD2 <dbl>, CREB3L1 <dbl>

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

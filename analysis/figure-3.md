Figure 4: All things prolactin
==============================

    library(tidyverse)

    ## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

    DEG_path <- "../results/DEseq2/hypothesis/"   # path to the data
    DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
    DEG_pathfiles <- paste0(DEG_path, DEG_files)
    #DEG_files

    allDEG2 <- DEG_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name") %>% 
      mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/hypothesis/'), "[", 2))  %>% 
      mutate(DEG = sapply(strsplit(as.character(DEG),'_diffexp.csv'), "[", 1))  %>% 
      mutate(tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1)) %>%
      mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%
      mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%
      mutate(comparison = paste(down,up, sep = "_")) %>%
      mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
      mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%
    dplyr::select(sex,tissue,comparison, direction, gene, lfc, padj, logpadj) 
    head(allDEG2)

    ## # A tibble: 6 x 8
    ##   sex    tissue comparison  direction gene           lfc     padj logpadj
    ##   <chr>  <chr>  <chr>       <chr>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  eggs_chicks chicks    OVALX         6.17 7.05e- 3    2.15
    ## 2 female gonad  eggs_chicks chicks    LOC101749216  5.19 2.45e- 6    5.61
    ## 3 female gonad  eggs_chicks chicks    BPIFB2        4.95 5.23e- 3    2.28
    ## 4 female gonad  eggs_chicks chicks    OvoDA1        4.90 2.69e- 4    3.57
    ## 5 female gonad  eggs_chicks chicks    CDC20B        4.07 2.99e-10    9.52
    ## 6 female gonad  eggs_chicks chicks    WFIKKN2       4.01 1.90e- 2    1.72

    allDEG2$tissue <- factor(allDEG2$tissue, levels = tissuelevel)

    allDEG2$comparison <- factor(allDEG2$comparison, levels = c("eggs_chicks", "lo_hi"))
    allDEG2 <- allDEG2 %>% mutate(comparison = fct_recode(comparison, "lo vs. hi PRL   " = "lo_hi",
                                                          "eggs vs. chicks" = "eggs_chicks"))
    allDEG2$direction <- factor(allDEG2$direction, levels = c("eggs", "chicks", "lo", "hi"))


    PRLDEGs <- allDEG2 %>%
      filter(tissue == "pituitary", comparison == "lo vs. hi PRL   ",
             direction == "hi") %>%
      arrange(desc(logpadj))
    PRLDEGs

    ## # A tibble: 3,022 x 8
    ##    sex    tissue   comparison      direction gene      lfc     padj logpadj
    ##    <chr>  <fct>    <fct>           <fct>     <chr>   <dbl>    <dbl>   <dbl>
    ##  1 male   pituita… "lo vs. hi PRL… hi        PRL      2.66 5.24e-31    30.3
    ##  2 female pituita… "lo vs. hi PRL… hi        PRL      2.55 1.81e-28    27.7
    ##  3 female pituita… "lo vs. hi PRL… hi        LAPTM4B  1.44 1.42e-20    19.8
    ##  4 male   pituita… "lo vs. hi PRL… hi        KPNA2    4.44 1.73e-20    19.8
    ##  5 male   pituita… "lo vs. hi PRL… hi        RACGAP1  3.57 1.73e-20    19.8
    ##  6 female pituita… "lo vs. hi PRL… hi        ARAP2    2.10 1.47e-19    18.8
    ##  7 male   pituita… "lo vs. hi PRL… hi        CKAP2    3.91 1.68e-19    18.8
    ##  8 female pituita… "lo vs. hi PRL… hi        LOC423…  3.60 2.93e-18    17.5
    ##  9 female pituita… "lo vs. hi PRL… hi        STT3B    1.05 1.08e-17    17.0
    ## 10 female pituita… "lo vs. hi PRL… hi        KPNA2    4.39 1.09e-17    17.0
    ## # … with 3,012 more rows

    envDEGs <- allDEG2 %>%
      filter(tissue == "pituitary", comparison == "eggs vs. chicks",
             direction == "chicks") %>%
      arrange(desc(logpadj))
    envDEGs

    ## # A tibble: 820 x 8
    ##    sex    tissue   comparison     direction gene        lfc    padj logpadj
    ##    <chr>  <fct>    <fct>          <fct>     <chr>     <dbl>   <dbl>   <dbl>
    ##  1 male   pituita… eggs vs. chic… chicks    FKBP5     0.676 3.59e-4    3.45
    ##  2 male   pituita… eggs vs. chic… chicks    NEMF      0.241 3.59e-4    3.45
    ##  3 female pituita… eggs vs. chic… chicks    VTN       1.98  9.75e-4    3.01
    ##  4 female pituita… eggs vs. chic… chicks    LOC10174… 0.268 9.75e-4    3.01
    ##  5 female pituita… eggs vs. chic… chicks    NKRF      0.187 1.63e-3    2.79
    ##  6 female pituita… eggs vs. chic… chicks    ZBTB16    0.777 1.94e-3    2.71
    ##  7 male   pituita… eggs vs. chic… chicks    TULP1     1.71  2.02e-3    2.69
    ##  8 female pituita… eggs vs. chic… chicks    CTNND2    0.261 2.48e-3    2.60
    ##  9 female pituita… eggs vs. chic… chicks    SEBOX     2.54  2.82e-3    2.55
    ## 10 male   pituita… eggs vs. chic… chicks    INSM1     0.551 3.15e-3    2.50
    ## # … with 810 more rows

    hilogenes <- PRLDEGs %>% head(20) %>% distinct(gene) %>% pull(gene)
    eggchickgenes <- envDEGs %>% head(20) %>% distinct(gene) %>% pull(gene)

    hypothesisDEGs <- c(hilogenes, eggchickgenes, WGCNAgenes) %>% unique()

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
    a2 <- plotprolactin(PRLpit, PRLpit$counts, "PRL", " ") + theme(axis.text.x = element_text())  

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

    # most sig DEGs

    c1 <-  plottopDEGs(envDEGs, "female", "#3A80B9", "LFC eggs vs. chicks", "Top 10 female DEGs") + labs(title = " ")
    c2 <-  plottopDEGs(PRLDEGs, "male", "#3A80B9", "LFC eggs vs. chicks", "Top 10 male DEGs") + labs(title = " ")
    c3 <-  plottopDEGs(PRLDEGs, "female", "#19757A", "LFC lo vs. hi PRL", NULL)
    c4 <-  plottopDEGs(PRLDEGs, "male", "#19757A", "LFC lo vs. hi PRL", NULL)


    b <- plot_grid(b1,b25,nrow = 1, rel_widths = c(1.2,1))

    c <- plot_grid(c1,c2,c3,c4, align = "v" , rel_heights = c(1.2,1))

    bc <- plot_grid(b,c, labels = c("B","C"), label_size = 8, rel_widths = c(1.2,1))

    d1 <-plotcorrplot(piteggchick, NULL) + theme(legend.position = "none") + labs(subtitle = "Top 20 eggs vs. chicks DEGs")
    d2 <- plotcorrplot(pithilo, NULL) + theme(legend.position = "none") + labs(subtitle = "Top 20 lo vs. PRL DEGs")
    d3 <-plotcorrplot(pitWGCNA, NULL) + theme(legend.position = "none") + labs(subtitle = "20 genes from WGCNA module with PRL")

    d <- plot_grid(d1,d2,d3, nrow = 1, align = "h", labels = c("D"), label_size = 8)


    fig3 <- plot_grid(a, bc, d, nrow = 3, rel_heights = c(1.2,2,1.8))
    fig3

![](../figures/fig3-1.png)

    #pdf(file="../figures/fig3.pdf", width=7.25, height=7.25)
    #plot(fig3)
    #dev.off()

    write.csv(PRLpit,"../results/PRLvsd.csv", row.names = F)

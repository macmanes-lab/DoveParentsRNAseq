Figure 4: All things prolactin
==============================

    library(tidyverse)

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

subset PRL vsd
==============

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
      arrange(desc(comparison))
    PRLDEGs

    ## # A tibble: 3,022 x 8
    ##    sex    tissue   comparison     direction gene       lfc     padj logpadj
    ##    <chr>  <fct>    <fct>          <fct>     <chr>    <dbl>    <dbl>   <dbl>
    ##  1 female pituita… "lo vs. hi PR… hi        KPNA2     4.39 1.09e-17   17.0 
    ##  2 female pituita… "lo vs. hi PR… hi        FOXM1     4.19 3.60e- 8    7.44
    ##  3 female pituita… "lo vs. hi PR… hi        FGF6      4.13 2.41e-16   15.6 
    ##  4 female pituita… "lo vs. hi PR… hi        PRC1      4.08 1.13e- 3    2.95
    ##  5 female pituita… "lo vs. hi PR… hi        LOC1017…  4.01 2.31e-17   16.6 
    ##  6 female pituita… "lo vs. hi PR… hi        SHCBP1    3.83 9.89e-13   12.0 
    ##  7 female pituita… "lo vs. hi PR… hi        RRM2      3.77 6.24e-15   14.2 
    ##  8 female pituita… "lo vs. hi PR… hi        CKAP2     3.77 4.76e-17   16.3 
    ##  9 female pituita… "lo vs. hi PR… hi        CDK1      3.68 6.24e-15   14.2 
    ## 10 female pituita… "lo vs. hi PR… hi        CCNB3     3.61 1.65e-12   11.8 
    ## # … with 3,012 more rows

    envDEGs <- allDEG2 %>%
      filter(tissue == "pituitary", comparison == "eggs vs. chicks",
             direction == "chicks") %>%
      arrange(desc(comparison))
    envDEGs

    ## # A tibble: 820 x 8
    ##    sex    tissue    comparison      direction gene     lfc     padj logpadj
    ##    <chr>  <fct>     <fct>           <fct>     <chr>  <dbl>    <dbl>   <dbl>
    ##  1 female pituitary eggs vs. chicks chicks    PNOC    3.88 0.0308      1.51
    ##  2 female pituitary eggs vs. chicks chicks    NCAN    3.07 0.0682      1.17
    ##  3 female pituitary eggs vs. chicks chicks    NPY     3.03 0.0197      1.70
    ##  4 female pituitary eggs vs. chicks chicks    CCK     2.82 0.0274      1.56
    ##  5 female pituitary eggs vs. chicks chicks    NKX2-1  2.79 0.0906      1.04
    ##  6 female pituitary eggs vs. chicks chicks    AGRP    2.57 0.0750      1.13
    ##  7 female pituitary eggs vs. chicks chicks    SEBOX   2.54 0.00282     2.55
    ##  8 female pituitary eggs vs. chicks chicks    MC4R    2.44 0.0344      1.46
    ##  9 female pituitary eggs vs. chicks chicks    HAND2   2.00 0.0310      1.51
    ## 10 female pituitary eggs vs. chicks chicks    VTN     1.98 0.000975    3.01
    ## # … with 810 more rows

    ## genes WGCNA prl module

    PRLgenes <- read_csv("../results/PRLmodule.csv") %>% pull(x)

    candidatevsd <- getcandidatevsd(PRLgenes, "pituitary", sexlevels) 

    df <- candidatevsd %>%
        pivot_wider(names_from = gene, values_from = counts) %>%
        select(-sex,-tissue, -treatment, -samples) %>%
        correlate() %>%
      focus(PRL)  %>%
      arrange(rowname)
    df

    ## # A tibble: 57 x 2
    ##    rowname     PRL
    ##    <chr>     <dbl>
    ##  1 ACOT12    0.137
    ##  2 AGR2      0.633
    ##  3 ARAP2     0.648
    ##  4 C2H8ORF46 0.253
    ##  5 CALN1     0.744
    ##  6 CD24      0.880
    ##  7 CDKN3     0.544
    ##  8 CDT1      0.499
    ##  9 CHRNB4    0.591
    ## 10 COL20A1   0.662
    ## # … with 47 more rows

    ## top correlations
    c <- df %>%
      arrange(desc(PRL)) %>% 
      head(10)  %>%
      mutate(PRLrounded = round(PRL, 2))%>% 
      ggplot(aes(x = reorder(rowname, PRL), y = PRL)) +
      geom_bar(stat = "identity") +
      theme_B3() +
      theme(axis.text.y = element_text(face = "italic")) +
      labs(x = " ", y = "Correlation with PRL", subtitle = "Top 10 co-regulated genes") +
      coord_flip() +
      scale_y_continuous(expand = c(0, 0))

    a1 <- plotpc12(pca1f, pca2f, pca1f$treatment, allcolors, "Female pituitary gene expression", " ")   
    a2 <- plotprolactin(PRLpitF, PRLpitF$counts, "PRL", " ") + theme(axis.text.x = element_text())  

    a3 <- plotpc12(pca1m, pca2m, pca1m$treatment, allcolors, "Male pituitary gene expression", NULL)
    a4 <- plotprolactin(PRLpitM, PRLpitM$counts, "PRL", " ") + theme(axis.text.x = element_text())

    a <- plot_grid(a1,a2,a3,a4, labels = c("A"), label_size = 8, rel_heights = c(1.1,1))


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


    fig3 <- plot_grid(a, bc, nrow = 2, rel_heights = c(1,1))
    fig3

![](../figures/fig3-1.png)

    #pdf(file="../figures/fig3.pdf", width=7.25, height=7.25)
    #plot(fig3)
    #dev.off()

    write.csv(PRLpit,"../results/PRLvsd.csv", row.names = F)

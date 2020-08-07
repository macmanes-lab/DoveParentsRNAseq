Fig 2
=====

    library(tidyverse)

    ## ── Attaching packages ────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(cowplot)

    ## 
    ## ********************************************************

    ## Note: As of version 1.0.0, cowplot does not change the

    ##   default ggplot2 theme anymore. To recover the previous

    ##   behavior, execute:
    ##   theme_set(theme_cowplot())

    ## ********************************************************

    library(ggsignif)

    source("../R/themes.R")
    source("../R/functions.R")

    knitr::opts_chunk$set(echo = TRUE, fig.path = '../figures/')

Treatment specific DEGs
-----------------------

    allDEG <- read_csv("../results/03_allDEG.csv") %>%
      mutate(tissue = factor(tissue, levels = tissuelevel),
             direction = factor(direction, levels = alllevels),
             label = gsub("_","\nvs. ", comparison)) 

    ## Parsed with column specification:
    ## cols(
    ##   sex = col_character(),
    ##   tissue = col_character(),
    ##   comparison = col_character(),
    ##   direction = col_character(),
    ##   gene = col_character(),
    ##   lfc = col_double(),
    ##   padj = col_double(),
    ##   logpadj = col_double()
    ## )

    head(allDEG)

    ## # A tibble: 6 x 9
    ##   sex    tissue comparison  direction gene       lfc     padj logpadj label     
    ##   <chr>  <fct>  <chr>       <fct>     <chr>    <dbl>    <dbl>   <dbl> <chr>     
    ## 1 female gonad  bldg_extend extend    CDK3     19.5  2.69e-20   19.6  "bldg\nvs…
    ## 2 female gonad  bldg_extend extend    CRISP2    6.52 3.48e- 3    2.46 "bldg\nvs…
    ## 3 female gonad  bldg_extend extend    KRT20     5.47 3.67e- 4    3.44 "bldg\nvs…
    ## 4 female gonad  bldg_extend extend    CLDN34    5.01 5.14e- 3    2.29 "bldg\nvs…
    ## 5 female gonad  bldg_extend extend    LOC1070…  4.89 7.51e- 2    1.12 "bldg\nvs…
    ## 6 female gonad  bldg_extend extend    OMD       3.53 5.38e- 4    3.27 "bldg\nvs…

    # for suppl figures
    DEGcontrol <- allDEG %>% 
      filter(grepl("control", comparison),
             !grepl("m.|early|extend|prolong", comparison))  %>%
      mutate(comparison = factor(comparison, levels = comparisonlevelscontrol))


    DEGbldg <- allDEG %>% 
      filter(grepl("bldg", comparison),
             !grepl("m.|early|extend|prolong", comparison))  %>%
      mutate(comparison = factor(comparison, levels = comparisonlevelsbldg))  %>%
      drop_na()

    DEGchar <- allDEG %>% 
      filter(comparison %in% comparisonlevelschar,
             !grepl("control|bldg", comparison)) %>%
      mutate(comparison = factor(comparison, levels = comparisonlevelschar))

    b <- makebargraph(DEGcontrol, "hypothalamus","DEGs w/ + LFC", 0, 4800, comparisonlabelscontrol) +
      labs(subtitle = "Control versus all other reproductive and parental stages")

    c <- makebargraph(DEGbldg, "hypothalamus", "DEGs w/ + LFC", 0, 480, comparisonlabelsbldg) +
      labs(subtitle = "Nest-building versus all other parental stages")
    d <- makebargraph(DEGchar, "hypothalamus","DEGs w/ + LFC", 0, 1700, comparisonlabelscharnobldg) +
      labs(subtitle = "Comparison of sequential parental stages")

    e <- makebargraph(DEGbldg, "pituitary", "DEGs w/ + LFC", 0, 2200, comparisonlabelsbldg)  

    f <- makebargraph(DEGchar, "pituitary","DEGs w/ + LFC", 0, 1500, comparisonlabelscharnobldg) 

    g <- makebargraph(DEGbldg, "gonad", "DEGs w/ + LFC", 0, 2800, comparisonlabelsbldg) 
    h <- makebargraph(DEGchar, "gonad","DEGs w/ + LFC", 0, 1200, comparisonlabelscharnobldg) 
    cd <- plot_grid(c,d,rel_widths = c(1,1))
    ef <- plot_grid(e,f,rel_widths = c(1,1),
                    labels = c("B Pitutary"), label_size = 8 ,hjust = 0)
    gh <- plot_grid(g,h,rel_widths = c(1,1), 
                    labels = c("C Gonads"), label_size = 8, hjust = 0)


    a <- plot.volcano("hypothalamus", sexlevels,  "control_bldg") + 
      facet_wrap(~sex) + labs(subtitle = " ")

    ## Warning in sex == whichsex: longer object length is not a multiple of shorter
    ## object length

    ab <- plot_grid(a,b,rel_widths = c(1,2.5), 
                    labels = c("A Hypothalamus"), label_size = 8, hjust = 0)

    fig <- plot_grid(ab,cd,ef, gh, ncol = 1)
    fig

![](../figures/fig2-new-1.png)

Save files
----------

    pdf(file="../figures/fig2-new-1.pdf", width=7, height=7)
    plot(fig)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    png("../figures/fig2-new-1.png", width = 7, height = 7, 
        units = 'in', res = 300)
    plot(fig) # Make plot
    dev.off()

    ## quartz_off_screen 
    ##                 2

    sessionInfo()

    ## R version 3.6.0 (2019-04-26)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS  10.15.4
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggsignif_0.5.0     cowplot_1.0.0.9000 forcats_0.4.0      stringr_1.4.0     
    ##  [5] dplyr_0.8.3        purrr_0.3.3        readr_1.3.1        tidyr_1.0.0       
    ##  [9] tibble_2.1.3       ggplot2_3.3.0.9000 tidyverse_1.3.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.5 xfun_0.15        haven_2.2.0      lattice_0.20-38 
    ##  [5] colorspace_1.4-1 vctrs_0.2.2      generics_0.0.2   htmltools_0.3.6 
    ##  [9] yaml_2.2.1       utf8_1.1.4       rlang_0.4.4      pillar_1.4.3    
    ## [13] glue_1.3.1       withr_2.1.2      DBI_1.1.0        dbplyr_1.4.2    
    ## [17] modelr_0.1.5     readxl_1.3.1     lifecycle_0.1.0  munsell_0.5.0   
    ## [21] gtable_0.3.0     cellranger_1.1.0 rvest_0.3.5      evaluate_0.14   
    ## [25] labeling_0.3     knitr_1.29       fansi_0.4.1      broom_0.5.2     
    ## [29] Rcpp_1.0.3       scales_1.1.0     backports_1.1.5  jsonlite_1.6.1  
    ## [33] farver_2.0.3     fs_1.3.1         hms_0.5.3        digest_0.6.24   
    ## [37] stringi_1.4.6    grid_3.6.0       cli_2.0.1        tools_3.6.0     
    ## [41] magrittr_1.5     crayon_1.3.4     pkgconfig_2.0.3  ellipsis_0.3.0  
    ## [45] xml2_1.2.2       reprex_0.3.0     lubridate_1.7.4  assertthat_0.2.1
    ## [49] rmarkdown_1.15   httr_1.4.1       rstudioapi_0.11  R6_2.4.1        
    ## [53] nlme_3.1-140     compiler_3.6.0

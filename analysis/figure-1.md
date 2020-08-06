Experimental design, tSNE analysis, and bar chars
=================================================

import limma counts and sample info
-----------------------------------

    countData <- read_csv("../results/01_limma.csv") %>%
      column_to_rownames(var = "X1")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )

    ## See spec(...) for full column specifications.

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

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_character(),
    ##   V1 = col_character(),
    ##   bird = col_character(),
    ##   sex = col_character(),
    ##   tissue = col_character(),
    ##   treatment = col_character(),
    ##   group = col_character(),
    ##   study = col_character()
    ## )

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

tsne
----

    # prep for tsne useing count data from limma and the custom `subsetmaketsne` function
    hyptsne <- subsetmaketsne("hypothalamus", alllevels, sexlevels)
    pittsne <- subsetmaketsne("pituitary", alllevels, sexlevels)
    gontsne <- subsetmaketsne("gonads", alllevels, sexlevels)

Figure 1
--------

    a <- png::readPNG("../figures/images/fig_fig1a.png")
    a <- ggdraw() +  draw_image(a, scale = 1)

    e1 <- plottsneelipsev3(hyptsne, hyptsne$treatment, allcolors) + 
      labs(subtitle = "hypothalamus") 
    e2 <- plottsneelipsev3(pittsne, pittsne$treatment, allcolors ) + 
      labs(subtitle = "pituitary") 
      theme(strip.text = element_blank())

    ## List of 1
    ##  $ strip.text: list()
    ##   ..- attr(*, "class")= chr [1:2] "element_blank" "element"
    ##  - attr(*, "class")= chr [1:2] "theme" "gg"
    ##  - attr(*, "complete")= logi FALSE
    ##  - attr(*, "validate")= logi TRUE

    e3 <- plottsneelipsev3(gontsne, gontsne$treatment, allcolors ) + 
      labs(subtitle = "gonads") 

    e <- plot_grid(e1,e2,e3, ncol = 3, labels = c("D"), label_size = 8)

    fig1 <- plot_grid(a, e, nrow = 2, rel_heights = c(2.4,1.6))
    fig1

![](../figures/fig1-1.png)

Save files
----------

    pdf(file="../figures/fig1-1.pdf", width=7, height=4.5)
    plot(fig1)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    png("../figures/fig1-1.png", width = 7, height = 4.5, 
        units = 'in', res = 300)
    plot(fig1) # Make plot
    dev.off()

    ## quartz_off_screen 
    ##                 2

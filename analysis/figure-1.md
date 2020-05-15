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

    # uses count data from limma

    # prep for tsne
    chartsne <- subsetmaketsne(tissuelevels, charlevels, sexlevels)
    hyptsne <- subsetmaketsne("hypothalamus", charlevels, sexlevels)
    pittsne <- subsetmaketsne("pituitary", charlevels, sexlevels)
    gontsne <- subsetmaketsne("gonads", charlevels, sexlevels)

    ftsne <-  subsetmaketsne(tissuelevels, charlevels, "female")
    mtsne <-  subsetmaketsne(tissuelevels, charlevels, "male")

Treatment specific DEGs
-----------------------

    allDEG <- read_csv("../results/03_allDEGs.csv") %>%
      mutate(comparison = factor(comparison, levels = comparisonlevelschar))

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_double(),
    ##   sex = col_character(),
    ##   tissue = col_character(),
    ##   comparison = col_character(),
    ##   direction = col_character(),
    ##   gene = col_character(),
    ##   lfc = col_double(),
    ##   padj = col_double(),
    ##   logpadj = col_double()
    ## )

    str(allDEG)

    ## Classes 'spec_tbl_df', 'tbl_df', 'tbl' and 'data.frame': 52792 obs. of  9 variables:
    ##  $ X1        : num  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ sex       : chr  "female" "female" "female" "female" ...
    ##  $ tissue    : chr  "gonad" "gonad" "gonad" "gonad" ...
    ##  $ comparison: Factor w/ 7 levels "bldg_lay","lay_inc.d3",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ direction : chr  "lay" "lay" "lay" "lay" ...
    ##  $ gene      : chr  "CDK3" "LOC107053414" "MUC" "OVSTL" ...
    ##  $ lfc       : num  18.63 9.72 5.83 5.45 5.04 ...
    ##  $ padj      : num  5.14e-18 5.07e-04 2.00e-03 8.53e-12 8.43e-04 ...
    ##  $ logpadj   : num  17.29 3.29 2.7 11.07 3.07 ...

    unique(allDEG$comparison)

    ## [1] bldg_lay       <NA>           hatch_n5       inc.d17_hatch 
    ## [5] inc.d3_inc.d9  inc.d9_inc.d17 lay_inc.d3     n5_n9         
    ## 7 Levels: bldg_lay lay_inc.d3 inc.d3_inc.d9 ... n5_n9

Sex and Tissue -related DEGs
----------------------------

    ## other DEG stuff
    hypsex <- read_csv("../results/DEseq2/sex/hypothalamus_female_male_DEGs.csv")

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    pitsex <- read_csv("../results/DEseq2/sex/pituitary_female_male_DEGs.csv")

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    gonsex <- read_csv("../results/DEseq2/sex/gonad_female_male_DEGs.csv")

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   tissue = col_character(),
    ##   direction = col_character()
    ## )

    hyppit <- read_csv("../results/DEseq2/tissue/hypothalamus_pituitary_DEGs.csv") 

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   direction = col_character()
    ## )

    hypgon <- read_csv("../results/DEseq2/tissue/hypothalamus_gonad_DEGs.csv")

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   direction = col_character()
    ## )

    pitgon <- read_csv("../results/DEseq2/tissue/pituitary_gonad_DEGs.csv")

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   direction = col_character()
    ## )

    hyppit$direction <- factor(hyppit$direction, levels = c("hypothalamus", "pituitary"))
    hypgon$direction <- factor(hypgon$direction, levels = c("hypothalamus",  "gonad"))
    pitgon$direction <- factor(pitgon$direction, levels = c( "pituitary", "gonad"))

make figure
-----------

    a <- png::readPNG("../figures/images/fig_fig1a.png")
    a <- ggdraw() +  draw_image(a, scale = 1)

    b1 <- plottsneelipse(chartsne, chartsne$tissue, allcolors)  + labs(y = "tSNE 2 ", subtitle = " ")  
    b2 <- sexbarplots(hyppit, 0, 8100) + labs(subtitle = " ", y = "DEGs w/ + LFC") +
      theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
      scale_x_discrete(labels=c("hypothalamus" = "hyp", "pituitary" = "pit" ))

    ## [1] 12419

    b3 <- sexbarplots(hypgon, 0, 8100) + labs(subtitle = " " ) +
      theme(axis.text.y = element_blank(), axis.line.y = element_blank(), 
            axis.ticks = element_blank(), axis.title.y = element_blank()) +
        scale_x_discrete(labels=c("hypothalamus" = "hyp",  "gonad" = "gon"))

    ## [1] 13122

    b4 <- sexbarplots(pitgon, 0, 8100) + labs(subtitle = " " ) +
      theme(axis.text.y = element_blank(), axis.line.y = element_blank(),  
            axis.ticks = element_blank(), axis.title.y = element_blank()) +
        scale_x_discrete(labels=c(  "pituitary" = "pit", "gonad" = "gon"))

    ## [1] 12999

    c1 <- plottsneelipse(chartsne, chartsne$sex, allcolors)   + labs(y = "tSNE 2 ", subtitle = " ")    
    c2 <- sexbarplots(hypsex, 0, 8100) + labs(y = "DEGs w/ + LFC", subtitle = "hypothalamus" ) +
      theme(axis.text.y = element_blank(), axis.ticks = element_blank())

    ## [1] 2206

    c3 <- sexbarplots(pitsex, 0, 8100) + labs(subtitle = "pituitary" ) + 
      theme(axis.text.y = element_blank(), axis.line.y = element_blank(),  
            axis.ticks = element_blank(), axis.title.y = element_blank())

    ## [1] 3649

    c4 <- sexbarplots(gonsex, 0, 8100) + labs(subtitle = "gonads" ) +
      theme(axis.text.y = element_blank(), axis.line.y = element_blank(),  
            axis.ticks = element_blank(), axis.title.y = element_blank())

    ## [1] 12972

    bc <- plot_grid(b1,b2,b3,b4, c1,c2,c3,c4, nrow = 1, rel_widths = c(1.5,1.1,0.9,0.9,1.5,1.1,0.9,0.9),
                    labels = c("B", "", "", "", "C"), label_size = 8)

    d1 <- plottsneelipsev2(hyptsne, hyptsne$treatment, allcolors) + 
      labs(x = NULL, subtitle = "hypothalamus")  + 
      facet_wrap(~sex, scales = "free")
    d3 <- plottsneelipsev2(pittsne, pittsne$treatment, allcolors ) + 
      labs(x = NULL, subtitle = "pituitary") + 
      facet_wrap(~sex, scales = "free") +
      theme(strip.text = element_blank())
    d5 <- plottsneelipsev2(gontsne, gontsne$treatment, allcolors ) + 
       facet_wrap(~sex, scales = "free") +
      theme(strip.text = element_blank()) +
      labs(x = "tSNE1 \n \n DEGs = differentially expressed genes \n + LFC = positive log fold change", 
           subtitle = "gonads")
    d2 <- makebargraph("hypothalamus","DEGs w/ + LFC", 0, 1800) + 
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())  
    d4 <- makebargraph("pituitary","DEGs w/ + LFC", 0, 1800)  +  
      theme(axis.text.x = element_blank(), 
            axis.title.x = element_blank(), 
            strip.text.x = element_blank())   
    d6 <- makebargraph("gonad","DEGs w/ + LFC", 0, 1800) +  
      theme(strip.text.x = element_blank()) +
      scale_x_discrete(labels = comparisonlabelschar)+
      labs(subtitle = " ")
    d <- plot_grid(d1,d2,d3,d4,d5,d6, ncol = 2, rel_heights = c(1,0.9,1.3), rel_widths = c(1,2),
                      labels = c("D"), label_size = 8)

    fig1 <- plot_grid(a,bc, d, nrow = 3, rel_heights = c(0.7,0.7,2),
                      labels = c("A"), label_size = 8)
    fig1

![](../figures/fig1-1.png)

    pdf(file="../figures/fig1-1.pdf", width=7.25, height=7.25)
    plot(fig1)
    dev.off()

    ## quartz_off_screen 
    ##                 2

supple table 1 of all but control-bldg DEGs
===========================================

    suppletable1 <- allDEG %>%
      #filter(comparison != "control_bldg") %>%
      group_by(sex, tissue, comparison) %>%
      arrange( tissue, sex, direction, gene)

    ## Warning: Factor `comparison` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`

    ## Warning: Factor `comparison` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`

    head(suppletable1)

    ## Warning: Factor `comparison` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`

    ## Warning: Factor `comparison` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`

    ## # A tibble: 6 x 9
    ## # Groups:   sex, tissue, comparison [1]
    ##      X1 sex    tissue comparison direction gene     lfc      padj logpadj
    ##   <dbl> <chr>  <chr>  <fct>      <chr>     <chr>  <dbl>     <dbl>   <dbl>
    ## 1  2816 female gonad  <NA>       bldg      A2ML3  0.679 0.0784       1.11
    ## 2  1761 female gonad  <NA>       bldg      A4GALT 1.19  0.00315      2.50
    ## 3  3674 female gonad  <NA>       bldg      AAAS   0.422 0.0740       1.13
    ## 4  3035 female gonad  <NA>       bldg      AAR2   0.606 0.0000263    4.58
    ## 5  2191 female gonad  <NA>       bldg      AARSD1 0.941 0.000315     3.50
    ## 6  3041 female gonad  <NA>       bldg      AASS   0.605 0.00520      2.28

    write_csv(suppletable1, "../results/suppletable1.csv")

    # save file for musical genes https://raynamharris.shinyapps.io/musicalgenes/
    #write.csv(allDEG, "../../musicalgenes/data/allDEG.csv")

    # tsne files too big to save
    #write.csv(chartsne, "../../musicalgenes/data/tsne.csv")

Experimental design, tSNE analysis, and bar chars
=================================================

tsne
----

    # prep for tsne
    chartsne <- subsetmaketsne(tissuelevels, charlevels, sexlevels)
    hyptsne <- subsetmaketsne("hypothalamus", charlevels, sexlevels)
    pittsne <- subsetmaketsne("pituitary", charlevels, sexlevels)
    gontsne <- subsetmaketsne("gonads", charlevels, sexlevels)

    ftsne <-  subsetmaketsne(tissuelevels, charlevels, "female")
    mtsne <-  subsetmaketsne(tissuelevels, charlevels, "male")

DEGs
----

    DEG_path <- "../results/DEseq2/"   # path to the data
    DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
    DEG_pathfiles <- paste0(DEG_path, DEG_files)
    #DEG_files

    allDEG <- DEG_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name") %>% 
      mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/'), "[", 2))  %>% 
      mutate(DEG = sapply(strsplit(as.character(DEG),'_diffexp.csv'), "[", 1))  %>% 
      mutate(tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1)) %>%
      mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%
      mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%
      mutate(comparison = paste(down,up, sep = "_")) %>%
      mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
      mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%
    dplyr::select(sex,tissue,comparison, direction, gene, lfc, padj, logpadj) 

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

    head(allDEG)

    ## # A tibble: 6 x 8
    ##   sex    tissue comparison direction gene           lfc     padj logpadj
    ##   <chr>  <chr>  <chr>      <chr>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  bldg_lay   lay       LOC107053414  9.72 4.76e- 4    3.32
    ## 2 female gonad  bldg_lay   lay       MUC           5.82 2.72e- 3    2.57
    ## 3 female gonad  bldg_lay   lay       OVSTL         5.45 9.58e-10    9.02
    ## 4 female gonad  bldg_lay   lay       AOC1          4.41 2.75e- 3    2.56
    ## 5 female gonad  bldg_lay   lay       ETNPPL        4.25 4.76e- 5    4.32
    ## 6 female gonad  bldg_lay   lay       GKN2          3.99 1.63e- 2    1.79

    allDEG$tissue <- factor(allDEG$tissue , levels = tissuelevel)
    allDEG$comparison <- factor(allDEG$comparison , levels = comparisonlevels)
    allDEG$direction <- factor(allDEG$direction, levels = charlevels)



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

make figure
-----------

    a <- png::readPNG("../figures/images/fig_fig1a.png")
    a <- ggdraw() +  draw_image(a, scale = 1)

    b <- plottsneelipse(chartsne, chartsne$tissue, allcolors)   + labs(y = "tSNE 2 ", subtitle = " ")  
    c <- plottsneelipse(chartsne, chartsne$sex, allcolors)   + labs(y = "tSNE 2 ", subtitle = " ")    

    d1 <- sexbarplots(hypsex, 0, 7000) + labs(y = "DEGS", subtitle = "hypothalamus" )
    d2 <- sexbarplots(pitsex, 0, 7000) + labs(subtitle = "pituitary" ) + 
      theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks = element_blank())
    d3 <- sexbarplots(gonsex, 0, 7000) + labs(subtitle = "gonads" ) +
      theme(axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks = element_blank())


    e <- plottsneelipsev2(hyptsne, hyptsne$treatment, allcolors) + 
      labs(x = NULL, subtitle = "hypothalamus")  + 
      facet_wrap(~sex, scales = "free")

    g <- plottsneelipsev2(pittsne, pittsne$treatment, allcolors ) + 
      labs(x = NULL, subtitle = "pituitary") + 
      facet_wrap(~sex, scales = "free") +
      theme(strip.text = element_blank())

    i <- plottsneelipsev2(gontsne, gontsne$treatment, allcolors ) + 
       facet_wrap(~sex, scales = "free") +
      theme(strip.text = element_blank()) +
      labs(x = "tSNE1 \n \n \n \n", subtitle = "gonads")


    # hyp
    f <- makebargraph("hypothalamus","DEGs", 0, 4300) + 
      theme(axis.text.x = element_blank(), 
            axis.title.x = element_blank())  +
      labs(subtitle = "hypothalamus")
    # pit
    h <- makebargraph("pituitary","DEGs", 0, 4300)  +  
      theme(axis.text.x = element_blank(), 
            axis.title.x = element_blank(), 
            strip.text.x = element_blank())  +
      labs(subtitle = "pituitary")
    # gon
    j <- makebargraph("gonad","DEGs", 0, 4300) +  
      theme(strip.text.x = element_blank()) +
      scale_x_discrete(labels = comparisonlabels)+
      labs(subtitle = "gonads")

    bcd <- plot_grid(b,c,d1,d2,d3, nrow = 1, labels = c("B", "C", "D"), 
                       label_size = 8, rel_widths = c(1.2,1.2,1.2,1,1))

    etoj <- plot_grid(e,f,g,h,i,j, ncol = 2, rel_heights = c(1,1,1.4), rel_widths = c(2.4,3.2),
                      labels = c("E", "F", "G", "H", "I", "J"), label_size = 8)

    ## Warning in MASS::cov.trob(data[, vars]): Probable convergence failure

    fig1 <- plot_grid(a, bcd, etoj, nrow = 3, rel_heights = c(0.8,0.8,2),
                      labels = c("A"), label_size = 8)
    fig1

![](../figures/fig1-1.png)

supple table 1 of all but control-bldg DEGs
===========================================

    suppletable1 <- allDEG %>%
      #filter(comparison != "control_bldg") %>%
      group_by(sex, tissue, comparison) %>%
      arrange( tissue, sex, direction, gene)
    head(suppletable1)

    ## # A tibble: 6 x 8
    ## # Groups:   sex, tissue, comparison [1]
    ##   sex    tissue      comparison   direction gene     lfc       padj logpadj
    ##   <chr>  <fct>       <fct>        <fct>     <chr>  <dbl>      <dbl>   <dbl>
    ## 1 female hypothalam… control_bldg control   A2ML1 -0.617    6.25e-2    1.20
    ## 2 female hypothalam… control_bldg control   A2ML3 -0.502    2.73e-6    5.56
    ## 3 female hypothalam… control_bldg control   AAED1 -0.568    5.86e-3    2.23
    ## 4 female hypothalam… control_bldg control   AATK  -0.312    6.48e-2    1.19
    ## 5 female hypothalam… control_bldg control   ABAT  -0.440    4.10e-8    7.39
    ## 6 female hypothalam… control_bldg control   ABCA2 -0.409    4.45e-2    1.35

    write_csv(suppletable1, "../results/suppletable1.csv")

    # save file for musical genes https://raynamharris.shinyapps.io/musicalgenes/
    #write.csv(allDEG, "../../musicalgenes/data/allDEG.csv")

    # tsne files too big to save
    #write.csv(chartsne, "../../musicalgenes/data/tsne.csv")

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

    makebargraph <- function(whichtissue, myylab, lowlim, higherlim){
      p <- allDEG %>%
        filter(tissue == whichtissue,
               comparison != "control_bldg") %>%
      ggplot(aes(x = comparison,  fill = direction)) +
        geom_bar(position = "dodge") +
        facet_grid(tissue~sex) +
        theme_B3() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")  +
        guides(fill = guide_legend(nrow = 1)) +
        labs(x = "Sequential parental care stage comparisons", 
             y = myylab,
             subtitle = " ") +
      scale_fill_manual(values = allcolors,
                           name = " ",
                           drop = FALSE) +
      scale_color_manual(values = allcolors) +
      geom_text(stat='count', aes(label=..count..), vjust =-0.5, 
                position = position_dodge(width = 1),
                size = 2, color = "black")  +
      ylim(lowlim, higherlim)
      return(p)
    }

Internal versus external
------------------------

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
    head(allDEG)

    ## # A tibble: 6 x 8
    ##   sex    tissue comparison direction gene           lfc     padj logpadj
    ##   <chr>  <fct>  <fct>      <fct>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  bldg_lay   lay       LOC107053414  9.72 4.76e- 4    3.32
    ## 2 female gonad  bldg_lay   lay       MUC           5.82 2.72e- 3    2.57
    ## 3 female gonad  bldg_lay   lay       OVSTL         5.45 9.58e-10    9.02
    ## 4 female gonad  bldg_lay   lay       AOC1          4.41 2.75e- 3    2.56
    ## 5 female gonad  bldg_lay   lay       ETNPPL        4.25 4.76e- 5    4.32
    ## 6 female gonad  bldg_lay   lay       GKN2          3.99 1.63e- 2    1.79

    allDEG2$tissue <- factor(allDEG2$tissue, levels = tissuelevel)

    allDEG2$comparison <- factor(allDEG2$comparison, levels = c("eggs_chicks", "lo_hi"))
    allDEG2 <- allDEG2 %>% mutate(comparison = fct_recode(comparison, "lo vs. hi PRL   " = "lo_hi",
                                                          "eggs vs. chicks" = "eggs_chicks"))
    allDEG2$direction <- factor(allDEG2$direction, levels = c("eggs", "chicks", "lo", "hi"))

DEGS
----

    # hyp
    p1 <- makebargraph("hypothalamus","DEGs", 0, 1250) + theme(axis.text.x = element_blank(), 
                                                                   axis.title.x = element_blank())

    # pit
    p2 <- makebargraph("pituitary","DEGs", 0, 1250)  +  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
                                                                  strip.text.x = element_blank())

    # gon
    p3 <- makebargraph("gonad","DEGs", 0, 1250) +  theme(strip.text.x = element_blank())

    bcd <- plot_grid(p1,p2,p3, nrow = 3, rel_heights = c(1.2,1,1.5), labels = c("b", "c", "d"), label_size = 8)


    expdesign <- png::readPNG("../figures/images/DoveParentsRNAseq_DEGS.png")
    expdesign <- ggdraw() +  draw_image(expdesign, scale = 1)

    plot_grid(expdesign, bcd, labels = c("a", " "), label_size = 8, nrow = 2, rel_heights = c(0.2,1))

![](../figures/fig2-1.png)

    makenewbargraph <- function(whichtissue, whichsex,  whichcomparison, lowlim, higherlim){
      p <- allDEG2 %>%
        filter(tissue == whichtissue,
               comparison == whichcomparison,
               sex == whichsex) %>%
        ggplot(aes(x = comparison,  fill = direction)) +
        geom_bar(position = "dodge", drop = FALSE) +
        theme_B3() +
        theme(legend.position = "none")  +
        guides(fill = guide_legend(nrow = 1)) +
        labs( y = whichtissue) +
      geom_text(stat='count', aes(label=..count..), vjust =-0.5, 
               position = position_dodge(width = 1),
               size = 2, color = "black")  + 
          ylim(lowlim, higherlim) +
      scale_fill_manual(values = allcolors, name = "higher in")   + 
        theme(axis.text.x = element_blank()) 
      return(p)
    }


    b11 <- makenewbargraph("hypothalamus", "female","eggs vs. chicks", 0, 2500) +  labs(subtitle = "females", x = NULL, y = "Hyp. DEGS")

    ## Warning: Ignoring unknown parameters: drop

    b21 <- makenewbargraph("pituitary", "female","eggs vs. chicks", 0, 2500)  + labs(x = NULL, y = "Pit. DEGS")

    ## Warning: Ignoring unknown parameters: drop

    b31 <- makenewbargraph("gonad", "female", "eggs vs. chicks", 0, 2500) + labs(x = "eggs vs. chicks", y = "Gon. DEGS")   

    ## Warning: Ignoring unknown parameters: drop

    b112131 <- plot_grid(b11,b21,b31, nrow = 3, rel_heights = c(1.1,1,1.1))


    b12 <- makenewbargraph("hypothalamus", "male",  "eggs vs. chicks", 0, 2500) + labs(subtitle = "males", x = NULL)+ 
      theme(axis.title.y = element_blank(), axis.text.y = element_blank())

    ## Warning: Ignoring unknown parameters: drop

    b22 <- makenewbargraph("pituitary", "male", "eggs vs. chicks", 0, 2500) + labs(x = NULL)+ 
      theme(axis.title.y = element_blank(), axis.text.y = element_blank())

    ## Warning: Ignoring unknown parameters: drop

    b32 <- makenewbargraph("gonad", "male", "eggs vs. chicks", 0, 2500) + labs(x = "eggs vs. chicks")   + 
      theme(axis.title.y = element_blank(), axis.text.y = element_blank())

    ## Warning: Ignoring unknown parameters: drop

    b122232 <- plot_grid(b12, b22,b32, nrow = 3, rel_heights = c(1.1,1,1.1)) 


    c11 <- makenewbargraph("hypothalamus", "female", "lo vs. hi PRL   ", 0, 2500) +  labs(subtitle = "females", x = NULL) + 
      theme(axis.title.y = element_blank())

    ## Warning: Ignoring unknown parameters: drop

    c21 <- makenewbargraph("pituitary", "female", "lo vs. hi PRL   ", 0, 2500)  + labs(x = NULL)+ 
      theme(axis.title.y = element_blank())

    ## Warning: Ignoring unknown parameters: drop

    c31 <- makenewbargraph("gonad", "female", "lo vs. hi PRL   ", 0, 2500) + labs(x = "lo vs. hi PRL")  + 
      theme(axis.title.y = element_blank()) 

    ## Warning: Ignoring unknown parameters: drop

    c112131 <- plot_grid(c11,c21,c31, nrow = 3, rel_heights = c(1.1,1,1.1))


    c12 <- makenewbargraph("hypothalamus", "male",  "lo vs. hi PRL   ", 0, 2500) + labs(subtitle = "males", x = NULL)+ 
      theme(axis.title.y = element_blank(), axis.text.y = element_blank())

    ## Warning: Ignoring unknown parameters: drop

    c22 <- makenewbargraph("pituitary", "male", "lo vs. hi PRL   ", 0, 2500) + labs(x = NULL)+ 
      theme(axis.title.y = element_blank(), axis.text.y = element_blank())

    ## Warning: Ignoring unknown parameters: drop

    c32 <- makenewbargraph("gonad", "male", "lo vs. hi PRL   ", 0, 2500) + labs(x = "lo vs. hi PRL")   + 
      theme(axis.title.y = element_blank(), axis.text.y = element_blank())

    ## Warning: Ignoring unknown parameters: drop

    c122232 <- plot_grid(c12, c22,c32, nrow = 3, rel_heights = c(1.1,1,1.1)) 


    hypothesisbars <-  plot_grid(b112131, b122232, c112131, c122232, nrow = 1, rel_widths = c(1.3, 1, 1.1, 1),
              labels = c("b", " ", "c", " "), label_size = 8)

    expdesign2 <- png::readPNG("../figures/images/DoveParentsRNAseq_hypothesis.png")
    expdesign2 <- ggdraw() +  draw_image(expdesign2, scale = 1)


    plot_grid(expdesign2, hypothesisbars, rel_heights = c(0.5,1), nrow = 2, labels = c("a", " "), label_size = 8)

![](../figures/fig2b-1.png)

total degs
----------

    allDEG %>%
      group_by(sex, tissue, comparison) %>%
      summarize(totalDEGs = n()) %>%
      arrange(tissue, comparison)

    ## # A tibble: 37 x 4
    ## # Groups:   sex, tissue [6]
    ##    sex    tissue       comparison     totalDEGs
    ##    <chr>  <fct>        <fct>              <int>
    ##  1 female hypothalamus control_bldg        5683
    ##  2 male   hypothalamus control_bldg        6683
    ##  3 female hypothalamus bldg_lay               1
    ##  4 male   hypothalamus bldg_lay               1
    ##  5 female hypothalamus inc.d3_inc.d9          1
    ##  6 female hypothalamus inc.d9_inc.d17         5
    ##  7 male   hypothalamus inc.d9_inc.d17        87
    ##  8 female hypothalamus inc.d17_hatch          3
    ##  9 male   hypothalamus inc.d17_hatch          6
    ## 10 female hypothalamus hatch_n5            1927
    ## # … with 27 more rows

candidate genes
---------------

    geneids <- read_csv("../metadata/00_geneinfo.csv")

    ## Warning: Missing column names filled in: 'X1' [1]

    candidategenes <- c("OXT", "AVP", "GNRH1", "GNRHR", "CGNRH-R",
                        "AR", "POMC", "AGRP",
                           "CRH", "AVPR1A", "AVPR1B", "AVPR2",
                           "CYP19A1", "DRD1", "DRD2", "PRL", "PRLR", "SOX9", 
                        "ESR1","ESR2", "LBH", "CDK1", "BRCA1",
                        "PTEN", "CREBBP", "FOS", "JUN", "EGR1",
                         "BDNF", "GRM2",
                        "KCNJ5", "CISH", "PTGER3", "CEBPD", "ZBTB16") 

    table1 <- allDEG %>%
      filter(gene %in% candidategenes,
             comparison != "group_bldg") %>%
        arrange(gene) %>%

      group_by(sex, tissue, comparison) %>%
      summarize(genes = str_c(gene, collapse = " ")) %>%
      pivot_wider(names_from = comparison, values_from = genes ) %>%
      select(sex, tissue, bldg_lay, lay_inc.d3, inc.d3_inc.d9,
            inc.d9_inc.d17,inc.d17_hatch, hatch_n5, n5_n9)  %>%
      arrange( tissue, sex)
    table1

    ## # A tibble: 6 x 9
    ## # Groups:   sex, tissue [6]
    ##   sex   tissue bldg_lay lay_inc.d3 inc.d3_inc.d9 inc.d9_inc.d17
    ##   <chr> <fct>  <chr>    <chr>      <chr>         <chr>         
    ## 1 fema… hypot… <NA>     <NA>       <NA>          <NA>          
    ## 2 male  hypot… <NA>     <NA>       <NA>          AR            
    ## 3 fema… pitui… ESR1 GN… ESR1 GNRH… <NA>          BRCA1 CDK1 KC…
    ## 4 male  pitui… <NA>     <NA>       <NA>          BRCA1 CDK1 CI…
    ## 5 fema… gonad  <NA>     AGRP AVPR… AVPR1A        SOX9          
    ## 6 male  gonad  <NA>     SOX9       <NA>          <NA>          
    ## # … with 3 more variables: inc.d17_hatch <chr>, hatch_n5 <chr>,
    ## #   n5_n9 <chr>

    write_csv(table1, "../results/table1.csv")

    suppletable1 <- allDEG %>%
      filter(comparison != "control_bldg") %>%
      group_by(sex, tissue, comparison) %>%
      arrange( tissue, sex, direction, gene)
    head(suppletable1)

    ## # A tibble: 6 x 8
    ## # Groups:   sex, tissue, comparison [3]
    ##   sex    tissue    comparison   direction gene         lfc     padj logpadj
    ##   <chr>  <fct>     <fct>        <fct>     <chr>      <dbl>    <dbl>   <dbl>
    ## 1 female hypothal… bldg_lay     bldg      HEMGN     -1.37  1.93e- 2    1.72
    ## 2 female hypothal… inc.d3_inc.… inc.d3    LOC1070… -17.2   9.56e-16   15.0 
    ## 3 female hypothal… inc.d9_inc.… inc.d9    CFAP44    -0.708 4.54e- 2    1.34
    ## 4 female hypothal… inc.d9_inc.… inc.d9    GMNN      -0.512 4.54e- 2    1.34
    ## 5 female hypothal… inc.d9_inc.… inc.d17   IGLL1      4.20  2.45e- 2    1.61
    ## 6 female hypothal… inc.d9_inc.… inc.d17   LOC1070…  17.8   9.74e-19   18.0

    suppletable1 %>%
      group_by(tissue) %>%
      summarize(totalDEGs = n())

    ## # A tibble: 3 x 2
    ##   tissue       totalDEGs
    ##   <fct>            <int>
    ## 1 hypothalamus      2032
    ## 2 pituitary         4440
    ## 3 gonad             3770

    write_csv(suppletable1, "../results/suppletable1.csv")

    allDEG2 %>%
      filter(gene %in% candidategenes,
             comparison != "group_bldg") %>%
      arrange(gene) %>%
      group_by(sex, tissue, comparison) %>%
      summarize(genes = str_c(gene, collapse = " ")) %>%
      pivot_wider(names_from = comparison, values_from = genes ) %>%
      arrange( tissue, sex)

    ## # A tibble: 4 x 4
    ## # Groups:   sex, tissue [6]
    ##   sex    tissue    `eggs vs. chicks`           `lo vs. hi PRL   `          
    ##   <chr>  <fct>     <chr>                       <chr>                       
    ## 1 male   hypothal… <NA>                        AR LBH                      
    ## 2 female pituitary AGRP AVPR1A BRCA1 CDK1 CEB… AGRP AR BRCA1 CDK1 ESR1 KCN…
    ## 3 male   pituitary CEBPD GRM2 ZBTB16           BRCA1 CDK1 KCNJ5 LBH PRL PR…
    ## 4 female gonad     AGRP AVPR1A                 <NA>

    suppletable2 <- allDEG2 %>%
      filter(comparison != "control_bldg") %>%
      group_by(sex, tissue, comparison) %>%
      arrange( tissue, sex, direction, gene)
    head(suppletable2)

    ## # A tibble: 6 x 8
    ## # Groups:   sex, tissue, comparison [2]
    ##   sex    tissue      comparison     direction gene      lfc    padj logpadj
    ##   <chr>  <fct>       <fct>          <fct>     <chr>   <dbl>   <dbl>   <dbl>
    ## 1 female hypothalam… eggs vs. chic… chicks    ALS2CL  0.456 0.00594    2.23
    ## 2 female hypothalam… eggs vs. chic… chicks    FKBP5   0.521 0.00668    2.18
    ## 3 male   hypothalam… eggs vs. chic… eggs      ABHD1… -0.292 0.0968     1.01
    ## 4 male   hypothalam… eggs vs. chic… eggs      ACAA2  -0.197 0.0768     1.11
    ## 5 male   hypothalam… eggs vs. chic… eggs      ACY1   -0.318 0.0852     1.07
    ## 6 male   hypothalam… eggs vs. chic… eggs      AES    -0.201 0.0961     1.02

    write_csv(suppletable2, "../results/suppletable1.csv")

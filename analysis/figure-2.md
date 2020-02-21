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
                size = 2.5, color = "black")  +
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

    makenewbargraph <- function(whichtissue, myxlab, whichcomparison, lowlim, higherlim){
      p <- allDEG2 %>%
        filter(tissue == whichtissue,
               comparison == whichcomparison) %>%
      ggplot(aes(x = comparison,  fill = direction)) +
        geom_bar(position = "dodge") +
        facet_grid(tissue~sex) +
        theme_B3() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")  +
        guides(fill = guide_legend(nrow = 1)) +
        labs(x = myxlab, 
             y = "Total DEGs") +
      geom_text(stat='count', aes(label=..count..), vjust =-0.5, 
               position = position_dodge(width = 1),
               size = 2.5, color = "black")  + 
          ylim(lowlim, higherlim)
      return(p)
    }

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

    bcd <- plot_grid(p1,p2,p3, nrow = 3, rel_heights = c(1.2,1,1.5), labels = c("b", "c", "d"), label_size = 12)


    expdesign <- png::readPNG("../figures/images/DoveParentsRNAseq_DEGS.png")
    expdesign <- ggdraw() +  draw_image(expdesign, scale = 1)

    plot_grid(expdesign, bcd, labels = c("a", " "), label_size = 12, nrow = 2, rel_heights = c(0.2,1))

![](../figures/fig2-1.png)

    b <- makenewbargraph(tissuelevel, "External stimuli", "eggs vs. chicks", 0, 1250)  + theme(strip.text.y = element_blank()) +
      scale_fill_manual(values = c("#80cdc1", "#018571"))
    c <- makenewbargraph(tissuelevel, "Pituitary PRL expression   ", "lo vs. hi PRL   ", 0, 1250) + labs(y = NULL) + theme(axis.text.y = element_blank()) +
      scale_fill_manual(values = c("#dfc27d", "#a6611a"))

     
    plot_grid(b, c, rel_widths = c(1,1), nrow = 1, labels = "auto", label_size = 12)

![](../figures/fig2-2.png)

total degs
----------

    allDEG %>%
      group_by(tissue, comparison) %>%
      summarize(totalDEGs = n()) %>%
      arrange(tissue, comparison)

    ## # A tibble: 23 x 3
    ## # Groups:   tissue [3]
    ##    tissue       comparison     totalDEGs
    ##    <fct>        <fct>              <int>
    ##  1 hypothalamus control_bldg       12366
    ##  2 hypothalamus bldg_lay               2
    ##  3 hypothalamus inc.d3_inc.d9          1
    ##  4 hypothalamus inc.d9_inc.d17        92
    ##  5 hypothalamus inc.d17_hatch          9
    ##  6 hypothalamus hatch_n5            1927
    ##  7 hypothalamus n5_n9                  1
    ##  8 pituitary    control_bldg       12789
    ##  9 pituitary    bldg_lay             279
    ## 10 pituitary    lay_inc.d3           358
    ## # … with 13 more rows

candidate genes
---------------

    candidategenes <- c("OXT", "AVP", "GNRH1",  "AR", "POMC", "AGRP",
                           "CRH", "AVPR1A", "AVPR1B", "AVPR2",
                           "CYP19A1", "DRD1", "DRD2", "PRL", "PRLR", "SOX9") 


    table1 <- allDEG %>%
      filter(gene %in% candidategenes) %>%
      group_by(sex, tissue, comparison) %>%
      summarize(genes = str_c(gene, collapse = " ")) %>%
      pivot_wider(names_from = comparison, values_from = genes ) %>%
      select(sex, tissue, control_bldg, lay_inc.d3, inc.d3_inc.d9,
            inc.d9_inc.d17, hatch_n5)  %>%
      arrange( tissue, sex)
    table1

    ## # A tibble: 6 x 7
    ## # Groups:   sex, tissue [6]
    ##   sex   tissue control_bldg lay_inc.d3 inc.d3_inc.d9 inc.d9_inc.d17
    ##   <chr> <fct>  <chr>        <chr>      <chr>         <chr>         
    ## 1 fema… hypot… DRD1 AR PRL… <NA>       <NA>          <NA>          
    ## 2 male  hypot… CRH AR AVPR… <NA>       <NA>          AR            
    ## 3 fema… pitui… AVPR2 AR DR… <NA>       <NA>          PRL           
    ## 4 male  pitui… AR AVPR2 AV… <NA>       <NA>          PRL           
    ## 5 fema… gonad  SOX9 AVPR1A… AVPR1A PR… AVPR1A        SOX9          
    ## 6 male  gonad  AR SOX9 PRL… SOX9       <NA>          <NA>          
    ## # … with 1 more variable: hatch_n5 <chr>

    write_csv(table1, "../results/table1.csv")

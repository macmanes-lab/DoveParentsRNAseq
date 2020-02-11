DEGs
----

    DEG_path <- "../results/DEseq2/manip/"   # path to the data
    DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
    DEG_pathfiles <- paste0(DEG_path, DEG_files)
    DEG_files

    ##  [1] "female_gonad_hatch_extend_DEGs.csv"            
    ##  [2] "female_gonad_inc.d17_m.inc.d17_DEGs.csv"       
    ##  [3] "female_gonad_inc.d17_prolong_DEGs.csv"         
    ##  [4] "female_gonad_inc.d3_m.inc.d3_DEGs.csv"         
    ##  [5] "female_gonad_inc.d9_m.inc.d8_DEGs.csv"         
    ##  [6] "female_gonad_inc.d9_m.inc.d9_DEGs.csv"         
    ##  [7] "female_gonad_m.n2_hatch_DEGs.csv"              
    ##  [8] "female_hypothalamus_hatch_extend_DEGs.csv"     
    ##  [9] "female_hypothalamus_inc.d17_m.inc.d17_DEGs.csv"
    ## [10] "female_hypothalamus_inc.d17_prolong_DEGs.csv"  
    ## [11] "female_hypothalamus_inc.d3_m.inc.d3_DEGs.csv"  
    ## [12] "female_hypothalamus_inc.d9_m.inc.d8_DEGs.csv"  
    ## [13] "female_hypothalamus_inc.d9_m.inc.d9_DEGs.csv"  
    ## [14] "female_hypothalamus_m.n2_hatch_DEGs.csv"       
    ## [15] "female_pituitary_hatch_extend_DEGs.csv"        
    ## [16] "female_pituitary_inc.d17_m.inc.d17_DEGs.csv"   
    ## [17] "female_pituitary_inc.d17_prolong_DEGs.csv"     
    ## [18] "female_pituitary_inc.d3_m.inc.d3_DEGs.csv"     
    ## [19] "female_pituitary_inc.d9_m.inc.d8_DEGs.csv"     
    ## [20] "female_pituitary_inc.d9_m.inc.d9_DEGs.csv"     
    ## [21] "female_pituitary_m.n2_hatch_DEGs.csv"          
    ## [22] "male_gonad_hatch_extend_DEGs.csv"              
    ## [23] "male_gonad_inc.d17_m.inc.d17_DEGs.csv"         
    ## [24] "male_gonad_inc.d17_prolong_DEGs.csv"           
    ## [25] "male_gonad_inc.d3_m.inc.d3_DEGs.csv"           
    ## [26] "male_gonad_inc.d9_m.inc.d8_DEGs.csv"           
    ## [27] "male_gonad_inc.d9_m.inc.d9_DEGs.csv"           
    ## [28] "male_gonad_m.n2_hatch_DEGs.csv"                
    ## [29] "male_hypothalamus_hatch_extend_DEGs.csv"       
    ## [30] "male_hypothalamus_inc.d17_m.inc.d17_DEGs.csv"  
    ## [31] "male_hypothalamus_inc.d17_prolong_DEGs.csv"    
    ## [32] "male_hypothalamus_inc.d3_m.inc.d3_DEGs.csv"    
    ## [33] "male_hypothalamus_inc.d9_m.inc.d8_DEGs.csv"    
    ## [34] "male_hypothalamus_inc.d9_m.inc.d9_DEGs.csv"    
    ## [35] "male_hypothalamus_m.n2_hatch_DEGs.csv"         
    ## [36] "male_pituitary_hatch_extend_DEGs.csv"          
    ## [37] "male_pituitary_inc.d17_m.inc.d17_DEGs.csv"     
    ## [38] "male_pituitary_inc.d17_prolong_DEGs.csv"       
    ## [39] "male_pituitary_inc.d3_m.inc.d3_DEGs.csv"       
    ## [40] "male_pituitary_inc.d9_m.inc.d8_DEGs.csv"       
    ## [41] "male_pituitary_inc.d9_m.inc.d9_DEGs.csv"       
    ## [42] "male_pituitary_m.n2_hatch_DEGs.csv"

    allDEG <- DEG_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name") %>% 
      mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/manip/'), "[", 2))  %>% 
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
    ##   sex    tissue comparison   direction gene           lfc     padj logpadj
    ##   <chr>  <chr>  <chr>        <chr>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  hatch_extend extend    LOC107053414 19.2  5.48e-13   12.3 
    ## 2 female gonad  hatch_extend extend    COL10A1       5.68 6.14e- 4    3.21
    ## 3 female gonad  hatch_extend extend    LOC107049904  5.30 6.20e- 2    1.21
    ## 4 female gonad  hatch_extend extend    SULT1C3       5.22 1.91e- 3    2.72
    ## 5 female gonad  hatch_extend extend    KRT20         4.87 1.46e- 2    1.84
    ## 6 female gonad  hatch_extend extend    LOC107057630  4.14 8.24e- 2    1.08

    comparisonlevels <- c( "inc.d9_m.inc.d8", "inc.d17_prolong", "hatch_extend" , 
                           "inc.d3_m.inc.d3" ,  "inc.d9_m.inc.d9", 
                            "inc.d17_m.inc.d17", "m.n2_hatch" 
                            )

    levelsreplace <- c( "m.inc.d8" , "prolong" , "extend")
    levelsremoval <- c( "m.inc.d3" ,    "m.inc.d9" , "m.inc.d17" , "m.n2")
    controlsremovalreplace <- c( "inc.d3" ,    "inc.d9" , "inc.d17" , "hatch")

    manipulationsandcontrols <- c(controlsremoval, levelsreplace, levelsremoval)



    allDEG$tissue <- factor(allDEG$tissue , levels = tissuelevel)
    allDEG$comparison <- factor(allDEG$comparison , levels = comparisonlevels)
    allDEG$direction <- factor(allDEG$direction, levels = manipulationsandcontrols)


    allDEG %>%
      group_by(sex, tissue, comparison) %>%
      summarize(n = n())

    ## # A tibble: 41 x 4
    ## # Groups:   sex, tissue [6]
    ##    sex    tissue       comparison            n
    ##    <chr>  <fct>        <fct>             <int>
    ##  1 female hypothalamus inc.d9_m.inc.d8      40
    ##  2 female hypothalamus inc.d17_prolong    3665
    ##  3 female hypothalamus hatch_extend       6330
    ##  4 female hypothalamus inc.d3_m.inc.d3    1734
    ##  5 female hypothalamus inc.d9_m.inc.d9      90
    ##  6 female hypothalamus inc.d17_m.inc.d17  4030
    ##  7 female hypothalamus m.n2_hatch         6041
    ##  8 female pituitary    inc.d9_m.inc.d8     254
    ##  9 female pituitary    inc.d17_prolong    1419
    ## 10 female pituitary    hatch_extend       2325
    ## # … with 31 more rows

    expdesign <- png::readPNG("../figures/images/fig5a_fig5a.png")
    expdesign <- ggdraw() +  draw_image(expdesign, scale = 1)

    replacementplots <- allDEG %>%
      filter(comparison %in% c("inc.d9_m.inc.d8", "inc.d17_prolong", "hatch_extend" ))  %>%
    ggplot(aes(x = comparison,  fill = direction)) +
      geom_bar(position = "dodge") +
      facet_grid(tissue~sex) +
      theme_B3() +
      theme(#axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.x = element_blank(),
            legend.position = "bottom",
            strip.text.y =  element_blank())  +
      guides(fill = guide_legend(nrow = 1)) +
      labs(x = "Offspring replacement", y = "Total differentially expressed genes") +
      scale_fill_manual(values = allcolors,
                           name = "increased in", drop = FALSE) +
      scale_color_manual(values = sexcolors) +
      geom_text(stat='count', aes(label=..count..), vjust =-0.5, 
                position = position_dodge(width = 1),
                size = 2.5, color = "black") +
      guides(fill =  guide_legend(nrow =2)) +
      ylim(0,3500)

    mylegend <- get_legend(replacementplots)

    removalplots <- allDEG %>%
      filter(comparison %in% c("inc.d3_m.inc.d3" ,  "inc.d9_m.inc.d9", 
                            "inc.d17_m.inc.d17", "m.n2_hatch" ))  %>%
    ggplot(aes(x = comparison,  fill = direction)) +
      geom_bar(position = "dodge") +
      facet_grid(tissue~sex) +
      theme_B3() +
      theme(#axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.x = element_blank(),
             axis.text.y = element_blank(),
            legend.position = "none",
            strip.placement = "right")  +
      guides(fill = guide_legend(nrow = 1)) +
      labs(x = "Offspring removal", y = NULL) +
      scale_fill_manual(values = allcolors,
                           name = "increased in", drop = FALSE) +
      scale_color_manual(values = sexcolors) +
      geom_text(stat='count', aes(label=..count..), vjust =-0.5, 
                position = position_dodge(width = 1),
                size = 2.5, color = "black") +
      ylim(0,3500)

    barplots <- plot_grid( replacementplots + theme(legend.position = "none"),  removalplots,  ncol = 2, rel_widths = c(1,1.2))

    fig5 <- plot_grid(expdesign, barplots, mylegend, ncol =1, rel_heights = c(0.25,1,0.1))
    fig5

![](../figures/fig5-1.png)

DEGs
----

    DEG_path <- "../results/DEseq2/"   # path to the data
    DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
    DEG_pathfiles <- paste0(DEG_path, DEG_files)
    DEG_files

    ##  [1] "female_gonad_bldg_lay_DEGs.csv"             
    ##  [2] "female_gonad_control_bldg_DEGs.csv"         
    ##  [3] "female_gonad_hatch_n5_DEGs.csv"             
    ##  [4] "female_gonad_inc.d17_hatch_DEGs.csv"        
    ##  [5] "female_gonad_inc.d3_inc.d9_DEGs.csv"        
    ##  [6] "female_gonad_inc.d9_inc.d17_DEGs.csv"       
    ##  [7] "female_gonad_lay_inc.d3_DEGs.csv"           
    ##  [8] "female_gonad_n5_n9_DEGs.csv"                
    ##  [9] "female_hypothalamus_bldg_lay_DEGs.csv"      
    ## [10] "female_hypothalamus_control_bldg_DEGs.csv"  
    ## [11] "female_hypothalamus_hatch_n5_DEGs.csv"      
    ## [12] "female_hypothalamus_inc.d17_hatch_DEGs.csv" 
    ## [13] "female_hypothalamus_inc.d3_inc.d9_DEGs.csv" 
    ## [14] "female_hypothalamus_inc.d9_inc.d17_DEGs.csv"
    ## [15] "female_hypothalamus_lay_inc.d3_DEGs.csv"    
    ## [16] "female_hypothalamus_n5_n9_DEGs.csv"         
    ## [17] "female_pituitary_bldg_lay_DEGs.csv"         
    ## [18] "female_pituitary_control_bldg_DEGs.csv"     
    ## [19] "female_pituitary_hatch_n5_DEGs.csv"         
    ## [20] "female_pituitary_inc.d17_hatch_DEGs.csv"    
    ## [21] "female_pituitary_inc.d3_inc.d9_DEGs.csv"    
    ## [22] "female_pituitary_inc.d9_inc.d17_DEGs.csv"   
    ## [23] "female_pituitary_lay_inc.d3_DEGs.csv"       
    ## [24] "female_pituitary_n5_n9_DEGs.csv"            
    ## [25] "male_gonad_bldg_lay_DEGs.csv"               
    ## [26] "male_gonad_control_bldg_DEGs.csv"           
    ## [27] "male_gonad_hatch_n5_DEGs.csv"               
    ## [28] "male_gonad_inc.d17_hatch_DEGs.csv"          
    ## [29] "male_gonad_inc.d3_inc.d9_DEGs.csv"          
    ## [30] "male_gonad_inc.d9_inc.d17_DEGs.csv"         
    ## [31] "male_gonad_lay_inc.d3_DEGs.csv"             
    ## [32] "male_gonad_n5_n9_DEGs.csv"                  
    ## [33] "male_hypothalamus_bldg_lay_DEGs.csv"        
    ## [34] "male_hypothalamus_control_bldg_DEGs.csv"    
    ## [35] "male_hypothalamus_hatch_n5_DEGs.csv"        
    ## [36] "male_hypothalamus_inc.d17_hatch_DEGs.csv"   
    ## [37] "male_hypothalamus_inc.d3_inc.d9_DEGs.csv"   
    ## [38] "male_hypothalamus_inc.d9_inc.d17_DEGs.csv"  
    ## [39] "male_hypothalamus_lay_inc.d3_DEGs.csv"      
    ## [40] "male_hypothalamus_n5_n9_DEGs.csv"           
    ## [41] "male_pituitary_bldg_lay_DEGs.csv"           
    ## [42] "male_pituitary_control_bldg_DEGs.csv"       
    ## [43] "male_pituitary_hatch_n5_DEGs.csv"           
    ## [44] "male_pituitary_inc.d17_hatch_DEGs.csv"      
    ## [45] "male_pituitary_inc.d3_inc.d9_DEGs.csv"      
    ## [46] "male_pituitary_inc.d9_inc.d17_DEGs.csv"     
    ## [47] "male_pituitary_lay_inc.d3_DEGs.csv"         
    ## [48] "male_pituitary_n5_n9_DEGs.csv"

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


    allDEG %>%
      group_by(sex, tissue, comparison) %>%
      summarize(n = n())

    ## # A tibble: 37 x 4
    ## # Groups:   sex, tissue [6]
    ##    sex    tissue       comparison         n
    ##    <chr>  <fct>        <fct>          <int>
    ##  1 female hypothalamus control_bldg    5683
    ##  2 female hypothalamus bldg_lay           1
    ##  3 female hypothalamus inc.d3_inc.d9      1
    ##  4 female hypothalamus inc.d9_inc.d17     5
    ##  5 female hypothalamus inc.d17_hatch      3
    ##  6 female hypothalamus hatch_n5        1927
    ##  7 female pituitary    control_bldg    6129
    ##  8 female pituitary    bldg_lay         279
    ##  9 female pituitary    lay_inc.d3       353
    ## 10 female pituitary    inc.d9_inc.d17  2004
    ## # … with 27 more rows

    p1 <- allDEG %>%
      filter(comparison == "control_bldg") %>%
    ggplot(aes(x = comparison,  fill = direction)) +
      geom_bar(position = "dodge") +
      facet_grid(tissue~sex) +
      theme_B3() +
      theme(#axis.text.x = element_text(angle = 45, hjust = 1),
             axis.text.x = element_blank(),
            legend.position = "none",
            strip.text.y = element_blank())  +
      guides(fill = guide_legend(nrow = 1)) +
      labs(x = " ", y = "Total Differentially Expressed Genes") +
      scale_fill_manual(values = colorscharnew,
                           name = " ",
                           drop = FALSE) +
      scale_color_manual(values = sexcolors)

    p2 <- allDEG %>%
      filter(comparison != "control_bldg") %>%
    ggplot(aes(x = comparison,  fill = direction)) +
      geom_bar(position = "dodge") +
      facet_grid(tissue~sex) +
      theme_B3() +
      theme(#axis.text.x = element_text(angle = 45, hjust = 1),
             axis.text.x = element_blank(),
            legend.position = "none",
            strip.placement = "right")  +
      guides(fill = guide_legend(nrow = 1)) +
      labs(x = "Sequential Parental Care Stages  ", y = NULL) +
      scale_fill_manual(values = colorscharnew,
                           name = " ", drop = FALSE) +
      scale_color_manual(values = sexcolors) +
      geom_text(stat='count', aes(label=..count..), vjust =-0.5, 
                position = position_dodge(width = 1),
                size = 2.5, color = "black") 
      
    forlegend <- allDEG %>%
    ggplot(aes(x = comparison,  fill = direction)) +
      geom_bar(position = "dodge") +
       scale_fill_manual(values = colorscharnew,
                           name = "increased in", drop = FALSE) +
      scale_color_manual(values = sexcolors) +
      guides(color = F) +
      guides(fill= guide_legend(nrow =1)) +
      theme_B3() +
      theme(legend.position = "bottom")

    mylegend <- get_legend(forlegend)

    barplots <- plot_grid(p1,p2, nrow=1, rel_widths = c(0.9,3))

    fig2 <- plot_grid(barplots, mylegend, nrow = 2, rel_heights = c(1,0.1))
    fig2

![](../figures/fig2-1.png)

    write_csv(allDEG, "../results/DEGs-char.csv")

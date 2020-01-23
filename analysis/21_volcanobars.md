DEGs
----

    DEG_path <- "../results/"   # path to the data
    DEG_files <- dir(DEG_path, pattern = "03_DEGs.*") # get file names
    DEG_pathfiles <- paste0(DEG_path, DEG_files)
    DEG_files

    ##  [1] "03_DEGs.gonads.bldg.lay.csv"            
    ##  [2] "03_DEGs.gonads.control.bldg.csv"        
    ##  [3] "03_DEGs.gonads.hatch.n5.csv"            
    ##  [4] "03_DEGs.gonads.inc.d17.hatch.csv"       
    ##  [5] "03_DEGs.gonads.inc.d3.inc.d9.csv"       
    ##  [6] "03_DEGs.gonads.inc.d9.inc.d17.csv"      
    ##  [7] "03_DEGs.gonads.lay.inc.d3.csv"          
    ##  [8] "03_DEGs.gonads.n5.n9.csv"               
    ##  [9] "03_DEGs.hypothalamus.bldg.lay.csv"      
    ## [10] "03_DEGs.hypothalamus.control.bldg.csv"  
    ## [11] "03_DEGs.hypothalamus.hatch.n5.csv"      
    ## [12] "03_DEGs.hypothalamus.inc.d17.hatch.csv" 
    ## [13] "03_DEGs.hypothalamus.inc.d3.inc.d9.csv" 
    ## [14] "03_DEGs.hypothalamus.inc.d9.inc.d17.csv"
    ## [15] "03_DEGs.hypothalamus.lay.inc.d3.csv"    
    ## [16] "03_DEGs.hypothalamus.n5.n9.csv"         
    ## [17] "03_DEGs.pituitary.bldg.lay.csv"         
    ## [18] "03_DEGs.pituitary.control.bldg.csv"     
    ## [19] "03_DEGs.pituitary.hatch.n5.csv"         
    ## [20] "03_DEGs.pituitary.inc.d17.hatch.csv"    
    ## [21] "03_DEGs.pituitary.inc.d3.inc.d9.csv"    
    ## [22] "03_DEGs.pituitary.inc.d9.inc.d17.csv"   
    ## [23] "03_DEGs.pituitary.lay.inc.d3.csv"       
    ## [24] "03_DEGs.pituitary.m.inc.d3.m.inc.d8.csv"
    ## [25] "03_DEGs.pituitary.n5.n9.csv"

    allDEG <- DEG_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name") %>% 
      mutate(DEG = sapply(strsplit(as.character(file_name),'../results/03_DEGs.'), "[", 2))  %>% 
      mutate(DEG = sapply(strsplit(as.character(DEG),'.csv'), "[", 1))  %>% 
      mutate(tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1)) %>%
       mutate(DEG = gsub("m.inc", "m_inc", DEG)) %>%
       mutate(DEG = gsub("inc.d", "inc_d", DEG))  %>%
      mutate(down = sapply(strsplit(as.character(DEG),'\\.'), "[", 2)) %>%
      mutate(up = sapply(strsplit(as.character(DEG),'\\.'), "[", 3)) %>%
      mutate(comparison = paste(down,up, sep = ".")) %>%
      select(tissue, comparison, direction, down, up, gene, lfc, padj, logpadj) 
    head(allDEG)

    ## # A tibble: 6 x 9
    ##   tissue comparison direction down  up    gene         lfc     padj logpadj
    ##   <chr>  <chr>      <chr>     <chr> <chr> <chr>      <dbl>    <dbl>   <dbl>
    ## 1 gonads bldg.lay   lay       bldg  lay   LOC107053…  9.81 4.65e- 9    8.33
    ## 2 gonads bldg.lay   lay       bldg  lay   MUC         5.81 3.36e- 4    3.47
    ## 3 gonads bldg.lay   lay       bldg  lay   OVSTL       5.42 2.22e-15   14.7 
    ## 4 gonads bldg.lay   lay       bldg  lay   AOC1        4.34 7.38e- 7    6.13
    ## 5 gonads bldg.lay   lay       bldg  lay   ETNPPL      4.24 1.97e- 5    4.71
    ## 6 gonads bldg.lay   lay       bldg  lay   GKN2        4.01 3.87e- 5    4.41

    allDEG$tissue <- factor(allDEG$tissue , levels = tissuelevels)
    allDEG$comparison <- factor(allDEG$comparison , levels = comparisonlevels)
    allDEG$direction <- factor(allDEG$direction, levels = charlevels)

    allDEG <- allDEG %>% drop_na()

    ggplot(allDEG, aes(x = lfc, y = logpadj)) +
      geom_point(aes(color = direction), alpha = 0.75) +
      facet_grid(tissue ~ comparison)  +
      scale_color_manual(values = colorscharnew,
                           name = " ",
                           drop = FALSE) +
      theme_B3() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")  +
      guides(color = guide_legend(nrow = 1))

![](../figures/favegenes/volcanobars,%20DEGs-1.png)

    ggplot(allDEG, aes(x = lfc)) +
      geom_histogram(aes( fill = direction)) +
      facet_grid(tissue ~ comparison)  +
      scale_color_manual(values = colorscharnew,
                           name = " ",
                           drop = FALSE) +
      theme_B3() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")  +
      guides(fill = guide_legend(nrow = 1)) 

![](../figures/favegenes/volcanobars,%20DEGs-2.png)

    ggplot(allDEG, aes(x = comparison,  fill = direction)) +
      geom_bar(position = "dodge") +
      facet_wrap(~tissue, nrow = 3) +
      theme_B3() +
      theme(axis.text.x = element_blank(),
            legend.position = "bottom",
            strip.placement = "right")  +
      guides(fill = guide_legend(nrow = 1)) +
      labs(x = "differntial expression between sequential stages", y = "total DEGs") +
      scale_fill_manual(values = colorscharnew,
                           name = " ",
                           drop = FALSE) +
      geom_text(stat='count', aes(label=..count..), vjust=-0.5, 
                position = position_dodge(width = 1),
                size = 3) +
      ylim(0,5000)

![](../figures/favegenes/volcanobars,%20DEGs-3.png)

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

    overlappingDEGs <- allDEG %>%
      filter(tissue == "hypothalamus",
             comparison != "control_bldg") %>%
      mutate(group = paste(sex, tissue,comparison, sep = ".")) %>%
      group_by(gene) %>%
      summarize(n = n()) %>%
      arrange(desc(n))
    overlappingDEGs 

    ## # A tibble: 2,008 x 2
    ##    gene          n
    ##    <chr>     <int>
    ##  1 C1QB          2
    ##  2 C4H4ORF50     2
    ##  3 CAPN2         2
    ##  4 CCDC60        2
    ##  5 COL1A1        2
    ##  6 CSF3R         2
    ##  7 DYTN          2
    ##  8 IGJ           2
    ##  9 IGLL1         2
    ## 10 KCNC4         2
    ## # … with 1,998 more rows

    HYPcandidategenes <- c("OXT", "AVP", "GNRH1",  "AR", "POMC", "AGRP",
                           "CRH", "AVPR1A", "AVPR1B", "AVPR2",
                           "CYP19A1", "DRD1", "DRD2", "PRL", "PRLR") 

    suppletable1 <- allDEG %>%
      filter(gene %in% HYPcandidategenes) %>%
      group_by(sex, tissue, comparison) %>%
      summarize(genes = str_c(gene, collapse = " ")) %>%
      pivot_wider(names_from = comparison, values_from = genes ) %>%
      select(sex, tissue, control_bldg, lay_inc.d3, inc.d3_inc.d9,
             inc.d9_inc.d17, hatch_n5)  %>%
      arrange(sex, tissue)
    suppletable1

    ## # A tibble: 6 x 7
    ## # Groups:   sex, tissue [6]
    ##   sex   tissue control_bldg lay_inc.d3 inc.d3_inc.d9 inc.d9_inc.d17
    ##   <chr> <fct>  <chr>        <chr>      <chr>         <chr>         
    ## 1 fema… hypot… DRD1 AR PRL… <NA>       <NA>          <NA>          
    ## 2 fema… pitui… AVPR2 AR DR… <NA>       <NA>          PRL           
    ## 3 fema… gonad  AVPR1A AR C… AVPR1A PR… AVPR1A        <NA>          
    ## 4 male  hypot… CRH AR AVPR… <NA>       <NA>          AR            
    ## 5 male  pitui… AR AVPR2 AV… <NA>       <NA>          PRL           
    ## 6 male  gonad  AR PRLR PRL… <NA>       <NA>          <NA>          
    ## # … with 1 more variable: hatch_n5 <chr>

    write.csv(suppletable1, "../results/suppltable-1.csv", row.names = F)

variance stabilized gene expression (vsd)
-----------------------------------------

    geneids <- read_csv("../metadata/00_geneinfo.csv")

    ## Warning: Missing column names filled in: 'X1' [1]

    vsd_path <- "../results/DEseq2/"   # path to the data
    vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
    vsd_pathfiles <- paste0(vsd_path, vsd_files)
    vsd_files

    ## [1] "female_gonad_vsd.csv"        "female_hypothalamus_vsd.csv"
    ## [3] "female_pituitary_vsd.csv"    "male_gonad_vsd.csv"         
    ## [5] "male_hypothalamus_vsd.csv"   "male_pituitary_vsd.csv"

    allvsd <- vsd_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name")  %>% 
      dplyr::rename("gene" = "X1") %>% 
      pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
                   names_to = "samples", values_to = "counts") 

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    p1 <- allDEG %>%
      filter(comparison == "control_bldg") %>%
      filter(tissue == "hypothalamus") %>%
    ggplot(aes(x = comparison,  fill = direction)) +
      geom_bar(position = "dodge") +
      facet_grid(tissue~sex) +
      theme_B3() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
             #axis.text.x = element_blank(),
            legend.position = "none",
            strip.text.y = element_blank())  +
      guides(fill = guide_legend(nrow = 1)) +
      labs(x = " ", y = "Hypothalamic DEGs") +
      scale_fill_manual(values = colorscharnew,
                           name = " ",
                           drop = FALSE) +
      scale_color_manual(values = sexcolors) +
      geom_text(stat='count', aes(label=..count..), vjust =-0.5, 
                position = position_dodge(width = 1),
                size = 2.5, color = "black")  +
      ylim(0,4200)

    p2 <- allDEG %>%
      filter(comparison != "control_bldg") %>%
        filter(tissue == "hypothalamus") %>%
    ggplot(aes(x = comparison,  fill = direction)) +
      geom_bar(position = "dodge") +
      facet_wrap(~sex) +
      theme_B3() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
             #axis.text.x = element_blank(),
            legend.position = "none",
            strip.placement = "right")  +
      guides(fill = guide_legend(nrow = 1)) +
      labs(x = "Sequential Parental Care Stages  ", y = NULL) +
      scale_fill_manual(values = colorscharnew,
                           name = " ", drop = FALSE) +
      scale_color_manual(values = sexcolors) +
      geom_text(stat='count', aes(label=..count..), vjust =-0.5, 
                position = position_dodge(width = 1),
                size = 2.5, color = "black")  +
      ylim(0,1150)
      
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


    plottopgenes <- function(whichgenes, whichtissue, whichsex, mysubtitle, whichtstages){
      candidates  <- allvsd %>%
        filter(gene %in% whichgenes) %>%
        dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
        dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/'), "[", 2)) %>%
        dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                    tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                    treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
        dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
        dplyr::select(sex, tissue, treatment, gene, samples, counts) %>%
        filter(tissue == whichtissue, sex %in% whichsex)  %>%
        drop_na()
      
      candidates$treatment <- factor(candidates$treatment, levels = alllevels)
      
      p <- candidates %>%
        filter(treatment %in% whichtstages) %>%
        ggplot(aes(x = treatment, y = counts, fill = treatment, color = sex)) +
        geom_boxplot() + 
        facet_wrap(~gene, scales = "free_y") + 
        theme_B3() +
        scale_fill_manual(values = allcolors) +
        scale_color_manual(values = allcolors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none",
              strip.text = element_text(face = "italic")) +
        labs(y = mysubtitle, x = "Parental stage") 
      return(p)
      
    }

    hypDEGs <- c("OXT", "AVP", "GNRH1",  "AR")
    topgenes <- plottopgenes(hypDEGs, "hypothalamus", c("female", "male"), "Hypothalamic expression", charlevels) 

    plot_grid(barplots, topgenes, nrow = 2, rel_heights = c(0.4,0.6))

![](../figures/fig2-1.png)

    kable(suppletable1)

<table>
<thead>
<tr>
<th style="text-align:left;">
sex
</th>
<th style="text-align:left;">
tissue
</th>
<th style="text-align:left;">
control\_bldg
</th>
<th style="text-align:left;">
lay\_inc.d3
</th>
<th style="text-align:left;">
inc.d3\_inc.d9
</th>
<th style="text-align:left;">
inc.d9\_inc.d17
</th>
<th style="text-align:left;">
hatch\_n5
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
hypothalamus
</td>
<td style="text-align:left;">
DRD1 AR PRLR PRL CYP19A1 POMC AGRP
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
DRD1 CYP19A1 POMC
</td>
</tr>
<tr>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
pituitary
</td>
<td style="text-align:left;">
AVPR2 AR DRD1 CYP19A1 PRLR PRL OXT AVP
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
PRL
</td>
<td style="text-align:left;">
AVPR2 PRL
</td>
</tr>
<tr>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
gonad
</td>
<td style="text-align:left;">
AVPR1A AR CYP19A1 POMC
</td>
<td style="text-align:left;">
AVPR1A PRLR AGRP
</td>
<td style="text-align:left;">
AVPR1A
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
male
</td>
<td style="text-align:left;">
hypothalamus
</td>
<td style="text-align:left;">
CRH AR AVPR2 CYP19A1 PRL GNRH1 AGRP POMC
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
AR
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
male
</td>
<td style="text-align:left;">
pituitary
</td>
<td style="text-align:left;">
AR AVPR2 AVPR1B CYP19A1 PRL POMC PRLR OXT AVP
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
PRL
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
male
</td>
<td style="text-align:left;">
gonad
</td>
<td style="text-align:left;">
AR PRLR PRL POMC
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
</tbody>
</table>

    #write_csv(allDEG, "../results/DEGs-char.csv")

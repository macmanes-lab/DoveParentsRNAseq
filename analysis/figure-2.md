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

    candidategenes <- c("OXT", "AVP", "GNRH1",  "AR", "POMC", "AGRP",
                           "CRH", "AVPR1A", "AVPR1B", "AVPR2",
                           "CYP19A1", "DRD1", "DRD2", "PRL", "PRLR",
                        "SOX9") 


    nocontrol <- allDEG %>%
      filter(comparison != "control_bldg")

    malegon <- nocontrol %>%
      filter(sex == "male", tissue == "gonad")

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

    getcandidatevsd <- function(whichgenes, whichtissue, whichsex){
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
      return(candidates)
    }

    hypvsd <- getcandidatevsd(candidategenes, "hypothalamus", sexlevels)
    pitvsd <- getcandidatevsd(candidategenes, "pituitary", sexlevels)
    gonvsd <- getcandidatevsd(candidategenes, "gonad", sexlevels)
    head(hypvsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue      treatment gene  samples                         counts
    ##   <chr>  <chr>       <fct>     <chr> <chr>                            <dbl>
    ## 1 female hypothalam… control   AGRP  L.G118_female_hypothalamus_con…   6.02
    ## 2 female hypothalam… control   AGRP  R.G106_female_hypothalamus_con…   5.77
    ## 3 female hypothalam… control   AGRP  R.R20_female_hypothalamus_cont…   5.46
    ## 4 female hypothalam… control   AGRP  R.R9_female_hypothalamus_contr…   5.68
    ## 5 female hypothalam… control   AGRP  R.W44_female_hypothalamus_cont…   5.81
    ## 6 female hypothalam… inc.d9    AGRP  blk.s061.pu.y_female_hypothala…   5.46

    head(pitvsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue    treatment gene  samples                           counts
    ##   <chr>  <chr>     <fct>     <chr> <chr>                              <dbl>
    ## 1 female pituitary control   AGRP  L.G118_female_pituitary_control.…   5.63
    ## 2 female pituitary control   AGRP  R.G106_female_pituitary_control     5.34
    ## 3 female pituitary control   AGRP  R.R20_female_pituitary_control      5.34
    ## 4 female pituitary control   AGRP  R.R9_female_pituitary_control.NY…   5.73
    ## 5 female pituitary control   AGRP  R.W44_female_pituitary_control.N…   5.34
    ## 6 female pituitary inc.d9    AGRP  blk.s061.pu.y_female_pituitary_i…   5.34

    head(gonvsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue treatment gene  samples                           counts
    ##   <chr>  <chr>  <fct>     <chr> <chr>                              <dbl>
    ## 1 female gonad  control   AGRP  L.G118_female_gonad_control        13.2 
    ## 2 female gonad  control   AGRP  R.G106_female_gonad_control         8.60
    ## 3 female gonad  control   AGRP  R.R20_female_gonad_control         13.6 
    ## 4 female gonad  control   AGRP  R.R9_female_gonad_control           9.82
    ## 5 female gonad  control   AGRP  R.W44_female_gonad_control         13.2 
    ## 6 female gonad  inc.d9    AGRP  blk.s061.pu.y_female_gonad_inc.d9   5.85

    makebargraph <- function(whichtissue, myylab, lowlim){
      p <- allDEG %>%
        filter(tissue == whichtissue,
               comparison != "control_bldg") %>%
      ggplot(aes(x = comparison,  fill = direction)) +
        geom_bar(position = "dodge") +
        facet_grid(tissue~sex) +
        theme_B3() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none",
              strip.text.y = element_blank())  +
        guides(fill = guide_legend(nrow = 1)) +
        labs(x = "Sequential parental care stage comparisons", 
             y = myylab,
             subtitle = " ") +
      scale_fill_manual(values = allcolors,
                           name = " ",
                           drop = FALSE) +
      scale_color_manual(values = allcolors)  #+
     # geom_text(stat='count', aes(label=..count..), vjust =-0.5, 
     #           position = position_dodge(width = 1),
     #           size = 2.5, color = "black")  +
     # ylim(lowlim, higherlim)
      return(p)
    }


    makeboxplots <- function(df, whichgene, myylab, whichsex){
      p <- df %>%
        filter(treatment %in% charlevels,
               gene %in% whichgene,
               sex %in% whichsex) %>%
        filter(treatment != "control") %>%
        ggplot(aes(x = treatment, y = counts, fill = treatment, color = sex)) +
        geom_boxplot() + 
        #geom_point() +
        facet_wrap(~sex) +
        theme_B3() +
        scale_fill_manual(values = allcolors) +
        scale_color_manual(values = allcolors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none") +
        labs(y = myylab , x = "Parental stage", subtitle = "") +
        theme(axis.title.y = element_markdown())
      return(p)
    }

HPG
---

    allDEG %>%
      filter(gene %in% candidategenes, 
             tissue == "hypothalamus") %>%
      group_by(sex, tissue, comparison) %>%
      summarize(genes = str_c(gene, collapse = " ")) %>%
      pivot_wider(names_from = comparison, values_from = genes ) %>%
      select(sex, tissue, control_bldg,  
            inc.d9_inc.d17, hatch_n5)  %>%
      arrange(sex, tissue)

    ## # A tibble: 2 x 5
    ## # Groups:   sex, tissue [6]
    ##   sex    tissue     control_bldg                inc.d9_inc.d17 hatch_n5    
    ##   <chr>  <fct>      <chr>                       <chr>          <chr>       
    ## 1 female hypothala… DRD1 AR PRLR PRL CYP19A1 P… <NA>           DRD1 CYP19A…
    ## 2 male   hypothala… CRH AR AVPR2 CYP19A1 PRL G… AR             <NA>

    write.csv(table1, "../results/table-1.csv", row.names = F)


    # hyp
    p1 <- makebargraph("hypothalamus","Hypothalamic DEGs") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    p2 <- makeboxplots(hypvsd, "AR","*AR* expression", "male") +
                  geom_signif(comparisons=list(c("inc.d9", "inc.d17")), annotations= "*", 
                  y_position = 8, tip_length = 0.005, vjust = 0.5, color= "black" , textsize = 3 ) +
                  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

    p3 <- makeboxplots(hypvsd, "DRD1","*DRD1* expression", "female") +
                  geom_signif(comparisons=list(c("hatch", "n5")), annotations= "*", 
                  y_position = 10, tip_length = 0.005, vjust = 0.5, color= "black" , textsize = 3 ) +
                  #geom_signif(comparisons=list(c("control", "bldg")), annotations= "*", 
                  #y_position = 10, tip_length = 0.005, vjust = 0.5, color= "black" , textsize = 3 ) + 
                  theme(axis.text.x = element_blank(), axis.title.x = element_blank())

    p23 <- plot_grid(p3,p2, labels = c("b", "c"), label_size = 12)

    # pit

    p4 <- makebargraph("pituitary","Pituitary DEGs")  +  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), strip.text = element_blank())
    p5 <- makeboxplots(pitvsd, "PRL","*PRL* expression", "male") +  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), strip.text = element_blank())+
                  #geom_signif(comparisons=list(c("control", "bldg")), annotations= "*", 
                  #y_position = 21.5, tip_length = 0.005, vjust = 0.5, color= "black" , textsize = 3 ) +
                  geom_signif(comparisons=list(c("inc.d9", "inc.d17")), annotations= "*", 
                  y_position = 21.5, tip_length = 0.005, vjust = 0.5, color= "black" , textsize = 3 ) 
    p6 <- makeboxplots(pitvsd, "PRL","*PRL* expression", "female")  +  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), strip.text = element_blank()) +
                  geom_signif(comparisons=list(c("inc.d9", "inc.d17")), annotations= "*", 
                  y_position = 21, tip_length = 0.005, vjust = 0.5, color= "black" , textsize = 3 ) +
                  #geom_signif(comparisons=list(c("control", "bldg")), annotations= "*", 
                  #y_position = 21, tip_length = 0.005, vjust = 0.5, color= "black" , textsize = 3 ) +
                  geom_signif(comparisons=list(c("hatch", "n5")), annotations= "*", 
                  y_position = 21, tip_length = 0.005, vjust = 0.5, color= "black" , textsize = 3 )

    p56 <- plot_grid(p6,p5, labels = c("e", "f"), label_size = 12)

    # gon

    p7 <- makebargraph("gonad","Gonadal DEGs") +  theme(strip.text = element_blank())
    p8 <- makeboxplots(gonvsd, "SOX9","*SOX9* expression", "male")  +  theme(strip.text = element_blank()) +
                  geom_signif(comparisons=list(c("lay", "inc.d3")), annotations= "*", 
                  y_position = 7.5, tip_length = 0.005, vjust = 0.5, color= "black" , textsize = 3 )
    p9 <- makeboxplots(gonvsd, "AVPR1A","*AVPR1A* expression", "female") +  theme(strip.text = element_blank()) +
                  #geom_signif(comparisons=list(c("control", "bldg")), annotations= "*", 
                  #y_position = 9.5, tip_length = 0.005, vjust = 0.5, color= "black" , textsize = 3 ) +
                  geom_signif(comparisons=list(c("lay", "inc.d3")), annotations= "*", 
                  y_position = 9.5, tip_length = 0.005, vjust = 0.5, color= "black" , textsize = 3 ) +
                  geom_signif(comparisons=list(c("inc.d3", "inc.d9")), annotations= "*", 
                  y_position = 9.5, tip_length = 0.005, vjust = 0.5, color= "black" , textsize = 3 )

    p89 <- plot_grid(p9,p8, labels = c("h", "i"), label_size = 12)


    a <- plot_grid(p1, p23, nrow = 1, rel_heights = c(0.5,0.5), labels = c("a", " "), label_size = 12)
    b <- plot_grid(p4, p56, nrow = 1, rel_heights = c(0.5,0.5), labels = c("d", " "), label_size = 12)
    c <- plot_grid(p7, p89, nrow = 1, rel_heights = c(0.5,0.5), labels = c("g", " "), label_size = 12)
    plot_grid(a,b,c, nrow = 3, rel_heights = c(1,0.9,1.3))

![](../figures/fig2-1.png)

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

Figure 4. What is the effect of removal?
========================================

import limma counts and sample info
-----------------------------------

    countData <- read_csv("../results/01_limma.csv") %>%
      column_to_rownames(var = "X1")
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

    # prep for tsne
    hyptsne <- subsetmaketsne("hypothalamus", removallevels, sexlevels)
    pittsne <- subsetmaketsne("pituitary", removallevels, sexlevels)
    gontsne <- subsetmaketsne("gonads", removallevels, sexlevels)

    addgroupings <- function(df){
      
      df <- df %>% mutate(earlylate = fct_collapse(treatment, 
                                      early = c("early", "m.inc.d3", "m.inc.d9", "inc.d3", 
                                             "inc.d9",  "bldg"),
                                      late = c("inc.d17",   "m.inc.d17", "prolong" ,  "hatch" , 
                                             "m.n2", "extend", "n5"),
                                      reference = "control"),
                    extint = fct_collapse(treatment, 
                                      nest = c("m.inc.d3", "m.inc.d9", "m.n2", "m.inc.d17", "bldg"),
                                      eggs = c("inc.d3", "inc.d9", "inc.d17", "prolong"),
                                      chicks = c("early", "hatch", "extend", "n5"),
                                      reference = "control")) 
      return(df)
    }
      

    pittsne <- addgroupings(pittsne)

    t4 <- plottsneelipsev2(pittsne, pittsne$treatment, allcolors) +  labs(title = "pituitary ", subtitle = "parental stage and manipulations") +
      facet_wrap(~sex) 
    t5 <- plottsneelipsev2(pittsne, pittsne$extint, colorhypothesis) +  labs(title = "  ", subtitle = "external enviornment")   +
      facet_wrap(~sex) + theme(axis.title.y = element_blank())
    t6 <- plottsneelipsev2(pittsne, pittsne$earlylate, colorhypothesis) +  labs(title = "  ", subtitle = "internal clock")   +
      facet_wrap(~sex) + theme(axis.title.y = element_blank())
    t456 <- plot_grid(t4, t5, t6, nrow = 1, rel_widths = c(1.1,1,1))

variance stabilized gene expression (vsd)
-----------------------------------------

    vsd_path <- "../results/DEseq2/treatment/"   # path to the data
    vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
    vsd_pathfiles <- paste0(vsd_path, vsd_files)
    vsd_files

    ## [1] "female_gonad_vsd.csv"        "female_gonads_vsd.csv"      
    ## [3] "female_hypothalamus_vsd.csv" "female_pituitary_vsd.csv"   
    ## [5] "male_gonad_vsd.csv"          "male_gonads_vsd.csv"        
    ## [7] "male_hypothalamus_vsd.csv"   "male_pituitary_vsd.csv"

    allvsd <- vsd_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name")  %>% 
      dplyr::rename("gene" = "X1") 
    ## before pivoting, check names of df with
    ## head(names(allvsd)) and tail(names(allvsd))
    allvsd <- allvsd %>% 
      pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
                   names_to = "samples", values_to = "counts") 
    allvsd <- as_tibble(allvsd)
    head(allvsd)

    ## # A tibble: 6 x 4
    ##   file_name                           gene  samples                  counts
    ##   <chr>                               <chr> <chr>                     <dbl>
    ## 1 ../results/DEseq2/treatment/female… A2ML1 L.G118_female_gonad_con…  14.7 
    ## 2 ../results/DEseq2/treatment/female… A2ML1 R.G106_female_gonad_con…   7.79
    ## 3 ../results/DEseq2/treatment/female… A2ML1 R.R20_female_gonad_cont…  14.9 
    ## 4 ../results/DEseq2/treatment/female… A2ML1 R.R9_female_gonad_contr…  11.7 
    ## 5 ../results/DEseq2/treatment/female… A2ML1 R.W44_female_gonad_cont…  13.2 
    ## 6 ../results/DEseq2/treatment/female… A2ML1 blk.s031.pu.d_female_go…   7.42

    tail(allvsd)

    ## # A tibble: 6 x 4
    ##   file_name                            gene  samples                 counts
    ##   <chr>                                <chr> <chr>                    <dbl>
    ## 1 ../results/DEseq2/treatment/male_pi… ZZZ3  y18.x_male_pituitary_m…  10.0 
    ## 2 ../results/DEseq2/treatment/male_pi… ZZZ3  y4.x_male_pituitary_m.…   9.88
    ## 3 ../results/DEseq2/treatment/male_pi… ZZZ3  y55.x_male_pituitary_m…   9.91
    ## 4 ../results/DEseq2/treatment/male_pi… ZZZ3  y63.x_male_pituitary_m…   9.92
    ## 5 ../results/DEseq2/treatment/male_pi… ZZZ3  y95.g131.x_male_pituit…   9.79
    ## 6 ../results/DEseq2/treatment/male_pi… ZZZ3  y98.o50.x_male_pituita…   9.89

    # for plotting gene expression over time
    getcandidatevsdmanip <- function(whichgenes, whichtissue){
      
      candidates  <- allvsd %>%
        filter(gene %in% whichgenes) %>%
        dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1),
                      sextissue = sapply(strsplit(sextissue, '../results/DEseq2/treatment/'), "[", 2),
                      sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                    tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                    treatment = sapply(strsplit(samples, '\\_'), "[", 4),
                    treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1),
                    treatment = recode(treatment, "m.hatch" = "m.n2", "extend.hatch" = "m.n2",
                                                       "inc.prolong" = "prolong", "m.inc.d8" = "early"),
                    treatment =  factor(treatment, levels = alllevels),
                    earlylate = fct_collapse(treatment, 
                                      early = c("early", "m.inc.d3", "m.inc.d9", "inc.d3", 
                                             "inc.d9",  "bldg"),
                                      late = c("inc.d17",   "m.inc.d17", "prolong" ,  "hatch" , 
                                             "m.n2", "extend", "n5"),
                                      reference = "control"),
                    extint = fct_collapse(treatment, 
                                      nest = c("m.inc.d3", "m.inc.d9", "m.n2", "m.inc.d17", "bldg"),
                                      eggs = c("inc.d3", "inc.d9", "inc.d17", "prolong"),
                                      chicks = c("early", "hatch", "extend", "n5"),
                                      reference = "control"))  %>%
        dplyr::select(gene, counts, sex, tissue, treatment, earlylate, extint, samples)  
      return(candidates)
    }

    candidatevsd <- getcandidatevsdmanip(candidategenes)
    hypvsd <- candidatevsd %>% filter(tissue %in% "hypothalamus")  
    pitvsd <- candidatevsd %>% filter(tissue %in% "pituitary") 
    gonvsd <- candidatevsd %>% filter(tissue %in% "gonad") 

    # not all genes appear to be expressed in the gonads... these are
    filter(gonvsd, treatment == "control") %>% head()

    ## # A tibble: 6 x 8
    ##   gene   counts sex    tissue treatment earlylate extint  samples          
    ##   <chr>   <dbl> <chr>  <chr>  <fct>     <fct>     <fct>   <chr>            
    ## 1 ACOT12   4.69 female gonad  control   reference refere… L.G118_female_go…
    ## 2 ACOT12   5.28 female gonad  control   reference refere… R.G106_female_go…
    ## 3 ACOT12   5.39 female gonad  control   reference refere… R.R20_female_gon…
    ## 4 ACOT12   5.66 female gonad  control   reference refere… R.R9_female_gona…
    ## 5 ACOT12   5.20 female gonad  control   reference refere… R.W44_female_gon…
    ## 6 ACOT12   5.65 female gonad  control   reference refere… blu.o.x.ATLAS_fe…

    plotcandidateremoval <- function(df, whichgene, whichtissue){
      
      p <- df %>%
        dplyr::filter(gene  == whichgene, 
                      treatment %in% removallevels)  %>%
        ggplot(aes(y =  counts, x = treatment, fill = treatment, color = sex)) +
        geom_boxplot(outlier.shape = NA, lwd=0.5) +
        geom_jitter(size = 0.25, aes(color = sex)) +
        theme_B3() +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              plot.subtitle = element_text(face = "italic")) +
        scale_color_manual(values = allcolors) +
        scale_fill_manual(values = allcolors)  +
        labs(title = paste(whichtissue, sep = " "), subtitle = whichgene, y = "gene expression") +
        geom_signif(comparisons = list(c( "inc.d3", "m.inc.d3"),
                                       c( "inc.d9", "m.inc.d9"),
                                       c( "inc.d17", "m.inc.d17"),
                                       c( "hatch", "m.n2")),
                    map_signif_level=TRUE,
                    textsize = 1.5, family = 'Helvetica',
                    vjust = 1.5, size = 0.5) 
      return(p)
      
    }

    a <- plotcandidateremoval(hypvsd, "HTR2C", "hypothalamus") + theme(axis.title.x = element_blank())
    b <- plotcandidateremoval(pitvsd, "PRL", "pituitary")+ theme(axis.title.x = element_blank(), strip.text = element_blank())
    c <- plotcandidateremoval(gonvsd, "ESR2", "gonads") + theme(strip.text = element_blank())
    boxplots1 <- plot_grid(a,b,c, nrow = 3)


    plotcandidatetiming <- function(df, whichgene, whichtissue){
      
      p <- df %>%
        dplyr::filter(gene  == whichgene, 
                      treatment %in% timinglevels,
                      !treatment %in%  c("control", "bldg", "inc.d3"))  %>%
        ggplot(aes(y =  counts, x = treatment, fill = treatment, color = sex)) +
        geom_boxplot(outlier.shape = NA, lwd=0.5) +
        geom_jitter(size = 0.25, aes(color = sex)) +
        theme_B3() +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              plot.subtitle = element_text(face = "italic")) +
        scale_color_manual(values = allcolors) +
        scale_fill_manual(values = allcolors)  +
        labs(title = paste(whichtissue, sep = " "), subtitle = whichgene, y = "gene expression") +
        geom_signif(comparisons = list(c("inc.d9", "early"), 
                                      c( "inc.d17", "prolong"), c( "prolong", "hatch"),
                                      c( "hatch", "extend"), c( "n5", "extend")),
                    map_signif_level=TRUE,
                    textsize = 1.5, family = 'Helvetica',
                    vjust = 1.5, size = 0.5) 
      return(p)
      
    }

    e <- plotcandidatetiming(hypvsd, "HTR2C", " ") + theme(axis.title.x = element_blank()) + labs(y = NULL, subtitle = " ", title = " ") 
    f <- plotcandidatetiming(pitvsd, "PRL", "pituitary")+ theme(axis.title.x = element_blank(), strip.text = element_blank()) + labs(y = NULL, subtitle = " ", title = " ") 
    g <- plotcandidatetiming(gonvsd, "ESR2", "gonads") + theme(strip.text = element_blank()) + labs(y = NULL, subtitle = " ", title = " ") 
    boxplots2 <- plot_grid(e,f,g, nrow = 3)
    boxplots2

![](../figures/boxplots-1.png)

DEGs
----

    DEG_path <- "../results/DEseq2/treatment/"   # path to the data
    DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
    DEG_pathfiles <- paste0(DEG_path, DEG_files)
    #DEG_files

    allDEG <- DEG_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name") %>% 
      mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/treatment/'), "[", 2))  %>% 
      mutate(DEG = sapply(strsplit(as.character(DEG),'_diffexp.csv'), "[", 1))  %>% 
      mutate(tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1)) %>%
      mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%
      mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%
      mutate(comparison = paste(down,up, sep = "_")) %>%
      mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
      mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%
      dplyr::select(sex,tissue,comparison, direction, gene, lfc, padj, logpadj) # %>%
     # mutate(tissue = factor(tissue, levels = tissuelevel),
      #       comparison = factor(comparison , levels = comparisonlevelsremoval),
      #       direction = factor(direction, levels = alllevels))
    head(allDEG)

    ## # A tibble: 6 x 8
    ##   sex    tissue comparison  direction gene           lfc     padj logpadj
    ##   <chr>  <chr>  <chr>       <chr>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  bldg_extend extend    CDK3         19.5  2.69e-20   19.6 
    ## 2 female gonad  bldg_extend extend    CRISP2        6.52 3.48e- 3    2.46
    ## 3 female gonad  bldg_extend extend    KRT20         5.47 3.67e- 4    3.44
    ## 4 female gonad  bldg_extend extend    CLDN34        5.01 5.14e- 3    2.29
    ## 5 female gonad  bldg_extend extend    LOC107049005  4.89 7.51e- 2    1.12
    ## 6 female gonad  bldg_extend extend    OMD           3.53 5.38e- 4    3.27

    candidateDEGS <- allDEG %>%
      filter(gene %in% parentalcaregenes) %>%
      mutate(posneg = ifelse(lfc >= 0, "+", "-"),
             sex = recode(sex, "female" = "F", "male" = "M" ),
             tissue = recode(tissue, 
                             "hypothalamus" = "H",
                             "pituitary" = "P", 
                             "gonad" = "G", "gonads" = "G")) %>%
      mutate(res = paste(sex, tissue, posneg, sep = "")) %>%
      select(gene, res, comparison)  %>%
      group_by(gene,  comparison) %>%
      summarize(res = str_c(res, collapse = " ")) %>%
      pivot_wider(names_from = comparison, values_from = res) %>%
      select(gene, inc.d3_m.inc.d3, inc.d9_m.inc.d9, inc.d17_m.inc.d17, hatch_m.n2,
             inc.d17_prolong, hatch_prolong,extend_hatch, n5_extend, 
             early_inc.d9 ,hatch_early)
    head(candidateDEGS[1:5])

    ## # A tibble: 6 x 5
    ## # Groups:   gene [6]
    ##   gene   inc.d3_m.inc.d3 inc.d9_m.inc.d9 inc.d17_m.inc.d17 hatch_m.n2
    ##   <chr>  <chr>           <chr>           <chr>             <chr>     
    ## 1 ADRA2A <NA>            <NA>            <NA>              <NA>      
    ## 2 AVP    <NA>            <NA>            FH-               <NA>      
    ## 3 AVPR1A <NA>            <NA>            <NA>              FH+ MP+   
    ## 4 BRINP1 <NA>            <NA>            FG- FH+           <NA>      
    ## 5 COMT   FH-             <NA>            FH-               FH- MH-   
    ## 6 CREBRF <NA>            <NA>            FP+ MP+           FP+ MP+

make figure
-----------

    design <- png::readPNG("../figures/images/fig_fig4a.png")
    design <- ggdraw() +  draw_image(design, scale = 1)


    boxplots <- plot_grid(boxplots1, boxplots2, nrow = 1, rel_widths = c(1.5,1),
                          labels = "AUTO", label_size = 8)
    fig4 <- plot_grid(boxplots, design, t456, ncol = 1, rel_heights = c(2.75,1.3,1),
                          labels = c(" ", " " , "C"), label_size = 8)
    fig4

![](../figures/fig4-1.png)

suppl fig 1
-----------

    comparisonlevelsmanip <- c("inc.d3_m.inc.d3", "inc.d9_m.inc.d9", "inc.d17_m.inc.d17", "hatch_m.n2",
             "early_inc.d9" ,"hatch_early", "inc.d17_prolong", "hatch_prolong","extend_hatch", "n5_extend")

    # for fig 1
    DEGmanip <- allDEG %>% 
      filter(comparison %in% comparisonlevelsmanip) %>%
      mutate(comparison = factor(comparison, levels = comparisonlevelsmanip))

    DEGmanip %>% 
      group_by(comparison, sex, tissue) %>%
      summarize(totalDEGs = n())  %>%
      arrange(desc(totalDEGs))

    ## # A tibble: 59 x 4
    ## # Groups:   comparison, sex [20]
    ##    comparison        sex    tissue       totalDEGs
    ##    <fct>             <chr>  <chr>            <int>
    ##  1 extend_hatch      female hypothalamus      5499
    ##  2 hatch_m.n2        female hypothalamus      5154
    ##  3 hatch_prolong     female hypothalamus      4979
    ##  4 hatch_early       female hypothalamus      4911
    ##  5 hatch_early       male   hypothalamus      4768
    ##  6 hatch_early       male   pituitary         4579
    ##  7 hatch_early       female pituitary         4198
    ##  8 inc.d17_m.inc.d17 female pituitary         3933
    ##  9 hatch_m.n2        female pituitary         3289
    ## 10 inc.d17_m.inc.d17 female hypothalamus      3149
    ## # … with 49 more rows

    s1a <- makebargraphsuppl(DEGmanip, "hypothalamus","DEGs w/ + LFC", 0, 3200) + 
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank()) + 
      labs(subtitle = "hypothalamus") + geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey") 

    s1b <- makebargraphsuppl(DEGmanip, "pituitary","DEGs w/ + LFC", 0, 3200) + 
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(), strip.text = element_blank()) + 
      labs(subtitle = "pituitary")  + geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey") 


    s1c <- makebargraphsuppl(DEGmanip, "gonad","DEGs w/ + LFC", 0, 3200) + 
      theme(strip.text = element_blank()) + 
      labs(subtitle = "gonad", x = "Comparison of parental stage and manipulation groups") + geom_vline(xintercept = 4.5, linetype = "dashed", color = "grey") + 
      scale_x_discrete(breaks= comparisonlevelsmanip,
                          labels= comparisonlevelsmanip,
                       drop=FALSE)

    contrasts <- png::readPNG("../figures/images/fig_fig4sup1.png")
    contrasts <- ggdraw() +  draw_image(contrasts, scale = 1)


    supplfig4 <- plot_grid(contrasts, s1a, s1b, s1c, ncol = 1, rel_heights = c(0.5, 1,1,1.5))
    supplfig4

![](../figures/supplfig-4-1.png)

write files
-----------

    pdf(file="../figures/fig4-1.pdf", width=7.25, height=7.25)
    plot(fig4)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/supplfig-4-1.pdf", width=5, height=5)
    plot(supplfig4)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    write.csv(candidateDEGS, "../results/table2.csv", row.names = F)

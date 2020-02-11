    geneids <- read_csv("../metadata/00_geneinfo.csv")

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

    hypgenes <- read_csv("../results/DEseq2/female_hypothalamus_hatch_n5_DEGs.csv") %>%
      arrange(desc(logpadj)) %>% head(., n = 2) %>% pull(gene)
    hypgenes

    ## [1] "CAPN2"        "LOC107054855"

    pitgenes <- read_csv("../results/DEseq2/female_pituitary_inc.d9_inc.d17_DEGs.csv") %>%
      arrange(desc(logpadj)) %>% head(., n = 2) %>% pull(gene)
    pitgenes

    ## [1] "LBH"  "CDK1"

    gongenes <- read_csv("../results/DEseq2/female_gonad_lay_inc.d3_DEGs.csv") %>%
      arrange(desc(logpadj)) %>% head(., n = 2) %>% pull(gene)
    gongenes

    ## [1] "OVSTL"  "BPIFB2"

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
        filter(tissue == whichtissue, sex %in% whichsex) 
      
      candidates$treatment <- factor(candidates$treatment, levels = alllevels)
      
      p <- candidates %>%
        filter(treatment %in% whichtstages) %>%
        ggplot(aes(x = treatment, y = counts, fill = treatment, color = sex)) +
        geom_boxplot() + facet_wrap(~gene, scales = "free_y") + 
        theme_B3() +
        scale_fill_manual(values = allcolors) +
        scale_color_manual(values = allcolors) +
        theme(axis.text.x = element_blank(),
              legend.position = "none",
              axis.title.x = element_blank(),
              strip.text = element_text(face = "italic")) +
        labs(y = mysubtitle) 
      return(p)
      
    }


    a <- plottopgenes(hypgenes, "hypothalamus", c("female"), "hypothalamus", c("hatch", "n5")) 
    b <- plottopgenes(pitgenes, "pituitary", c("female"), "pituitary", c("inc.d9", "inc.d17")) 
    c <- plottopgenes(gongenes, "gonad", c("female"), "gonads", c("lay", "inc.d3"))

    fig4a <- plot_grid(a,b,c, ncol = 1) 

    d <- plottopgenes(hypgenes, "hypothalamus", c("female", "male"), NULL, charlevels) 
    e <- plottopgenes(pitgenes, "pituitary", c("female", "male"), NULL, charlevels) 
    f <- plottopgenes(gongenes, "gonad", c("female", "male"), NULL, charlevels) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    fig4b <- plot_grid(d,e,f, ncol = 1)

    fig4 <- plot_grid(fig4a, fig4b, rel_widths = c(1,2))
    fig4

![](../figures/fig4-1.png)

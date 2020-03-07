Figure 3
========

    library(tidyverse)

    ## ── Attaching packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    library(ggimage)

    ## 
    ## Attaching package: 'ggimage'

    ## The following object is masked from 'package:cowplot':
    ## 
    ##     theme_nothing

    library(apaTables)
    library(factoextra)

    ## Welcome! Related Books: `Practical Guide To Cluster Analysis in R` at https://goo.gl/13EFCZ

    source("../R/themes.R") 
    source("../R/functions.R")
    source("../R/wrangledata.R")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )

    ## See spec(...) for full column specifications.

    knitr::opts_chunk$set(fig.path = '../figures/',message=F, warning=FALSE)

Data: PCA, hormones, PRL expression
-----------------------------------

    # pca
    charpca <- subsetmakepca("pituitary", charlevels, sexlevels)    
    charfviz <- makefvizdf("pituitary", charlevels, sexlevels)

    # hormones
    prolactin <- read_csv("../results/07_hormones.csv") %>%
        filter(study == "characterization", hormone %in% c("prolactin"))  %>% 
        droplevels() 
    prolactin$treatment <- factor(prolactin$treatment, levels = alllevels)

    # PRL vsd

    vsd_path <- "../results/DEseq2/"   # path to the data
    vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
    vsd_pathfiles <- paste0(vsd_path, vsd_files)
    vsd_files

    ## [1] "female_gonad_vsd.csv"        "female_hypothalamus_vsd.csv"
    ## [3] "female_pituitary_vsd.csv"    "male_gonad_vsd.csv"         
    ## [5] "male_hypothalamus_vsd.csv"   "male_pituitary_vsd.csv"

    PRLvsd <- vsd_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name")  %>% 
      dplyr::rename("gene" = "X1") %>% 
      pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
                   names_to = "samples", values_to = "counts") %>%
      filter(gene == "PRL")  %>%
      dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
      dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/'), "[", 2)) %>%
      dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                    tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                    treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
      dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
      dplyr::select(sex, tissue, treatment, gene, samples, counts) %>%
      drop_na()
    PRLvsd$treatment <- factor(PRLvsd$treatment, levels = alllevels) 

    PRLhyp <- PRLvsd %>% filter(tissue == "hypothalamus") 
    PRLpit <- PRLvsd %>% filter(tissue == "pituitary")
    PRLgon <- PRLvsd %>% filter(tissue == "gonad")

    PRLvsd2 <- PRLvsd %>%
      filter(tissue == "pituitary")
    head(PRLvsd2)

    ## # A tibble: 6 x 6
    ##   sex    tissue    treatment gene  samples                           counts
    ##   <chr>  <chr>     <fct>     <chr> <chr>                              <dbl>
    ## 1 female pituitary control   PRL   L.G118_female_pituitary_control.…   17.9
    ## 2 female pituitary control   PRL   R.G106_female_pituitary_control     17.0
    ## 3 female pituitary control   PRL   R.R20_female_pituitary_control      18.6
    ## 4 female pituitary control   PRL   R.R9_female_pituitary_control.NY…   16.8
    ## 5 female pituitary control   PRL   R.W44_female_pituitary_control.N…   18.6
    ## 6 female pituitary inc.d9    PRL   blk.s061.pu.y_female_pituitary_i…   17.6

    PRLvsd2 %>%
      group_by(sex) %>%
      summarize(median = median(counts))

    ## # A tibble: 2 x 2
    ##   sex    median
    ##   <chr>   <dbl>
    ## 1 female   18.1
    ## 2 male     17.9

    PRLvsd3 <- PRLvsd2 %>%
      mutate(hiloPRL = ifelse(counts >= 18, "hi", "lo"))  %>%
      drop_na()
    PRLvsd3$hiloPRL <- factor(PRLvsd3$hiloPRL, levels = c("lo", "hi"))


    PRLvsd3 %>%
      group_by(sex, tissue, hiloPRL) %>%
      summarize(n = n())

    ## # A tibble: 4 x 4
    ## # Groups:   sex, tissue [2]
    ##   sex    tissue    hiloPRL     n
    ##   <chr>  <chr>     <fct>   <int>
    ## 1 female pituitary lo         47
    ## 2 female pituitary hi         49
    ## 3 male   pituitary lo         51
    ## 4 male   pituitary hi         46

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
    head(allDEG2)

    ## # A tibble: 6 x 8
    ##   sex    tissue comparison  direction gene           lfc     padj logpadj
    ##   <chr>  <chr>  <chr>       <chr>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  eggs_chicks chicks    OVALX         6.17 7.05e- 3    2.15
    ## 2 female gonad  eggs_chicks chicks    LOC101749216  5.19 2.45e- 6    5.61
    ## 3 female gonad  eggs_chicks chicks    BPIFB2        4.95 5.23e- 3    2.28
    ## 4 female gonad  eggs_chicks chicks    OvoDA1        4.90 2.69e- 4    3.57
    ## 5 female gonad  eggs_chicks chicks    CDC20B        4.07 2.99e-10    9.52
    ## 6 female gonad  eggs_chicks chicks    WFIKKN2       4.01 1.90e- 2    1.72

    allDEG2$tissue <- factor(allDEG2$tissue, levels = tissuelevel)

    allDEG2$comparison <- factor(allDEG2$comparison, levels = c("eggs_chicks", "lo_hi"))
    allDEG2 <- allDEG2 %>% mutate(comparison = fct_recode(comparison, "lo vs. hi PRL   " = "lo_hi",
                                                          "eggs vs. chicks" = "eggs_chicks"))
    allDEG2$direction <- factor(allDEG2$direction, levels = c("eggs", "chicks", "lo", "hi"))

    expdesign <- png::readPNG("../figures/images/fig_fig3a.png")
    expdesign <- ggdraw() +  draw_image(expdesign, scale = 1)


    a <- plotcolorfulpcs(charpca,charpca$treatment, allcolors) + labs(subtitle = " ") +
      theme(legend.position = "none", 
            legend.direction = "horizontal", 
            legend.key.size = unit(0.5, 'lines')) + 
      guides(color = FALSE) +
      labs(subtitle = "pituitary")   

    b <- plotfriz(charfviz) + labs(subtitle = "pituitary") +
      #ylim(-250000,500000) + xlim(-250000,500000) + 
      theme(axis.text = element_blank())

    c <- plotprolactin(PRLpit, PRLpit$counts, "PRL", "pituitary") + 
      theme(legend.position = "none", 
             axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.y = element_text(face = "italic"))

    d <- plotprolactin(prolactin, prolactin$plasma_conc, "prolactin (ng/mL)", "blood") +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) 

    bcde <- plot_grid(a,b,c,d, labels = c("b", "c", "d", "e"), ncol = 2, label_size = 8)



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


    b11 <- makenewbargraph("hypothalamus", "female","eggs vs. chicks", 0, 2500) +  labs(subtitle = "females", x = NULL, y = " ") + 
      theme( axis.text = element_blank())
    b21 <- makenewbargraph("pituitary", "female","eggs vs. chicks", 0, 2500)  + labs(x = NULL, y = " ") + 
      theme(axis.text = element_blank())
    b31 <- makenewbargraph("gonad", "female", "eggs vs. chicks", 0, 2500) + labs(x = "eggs vs. chicks", y = " ")    + 
      theme(axis.text = element_blank())
    b112131 <- plot_grid(b11,b21,b31, nrow = 3, rel_heights = c(1.1,1,1.1))


    b12 <- makenewbargraph("hypothalamus", "male",  "eggs vs. chicks", 0, 2500) + labs(subtitle = "males", x = NULL)+ 
      theme(axis.title.y = element_blank(), axis.text = element_blank())
    b22 <- makenewbargraph("pituitary", "male", "eggs vs. chicks", 0, 2500) + labs(x = NULL)+ 
      theme(axis.title.y = element_blank(), axis.text = element_blank())
    b32 <- makenewbargraph("gonad", "male", "eggs vs. chicks", 0, 2500) + labs(x = "eggs vs. chicks")   + 
      theme(axis.title.y = element_blank(), axis.text = element_blank())
    b122232 <- plot_grid(b12, b22,b32, nrow = 3, rel_heights = c(1.1,1,1.1)) 




    c11 <- makenewbargraph("hypothalamus", "female", "lo vs. hi PRL   ", 0, 2500) +  labs(subtitle = "females", x = NULL, y = "hypothalamus DEGs") + 
      theme( axis.text.x = element_blank()) 
    c21 <- makenewbargraph("pituitary", "female", "lo vs. hi PRL   ", 0, 2500)  + labs(x = NULL, y = "pituitary DEGs")+ 
      theme( axis.text.x = element_blank())
    c31 <- makenewbargraph("gonad", "female", "lo vs. hi PRL   ", 0, 2500) + labs(x = "lo vs. hi PRL", y = "gonad DEGs")  + 
      theme(axis.text.x = element_blank()) 
    c112131 <- plot_grid(c11,c21,c31, nrow = 3, rel_heights = c(1.1,1,1.1))


    c12 <- makenewbargraph("hypothalamus", "male",  "lo vs. hi PRL   ", 0, 2500) + labs(subtitle = "males", x = NULL)+ 
      theme(axis.title.y = element_blank(), axis.text = element_blank())
    c22 <- makenewbargraph("pituitary", "male", "lo vs. hi PRL   ", 0, 2500) + labs(x = NULL)+ 
      theme(axis.title.y = element_blank(), axis.text = element_blank())
    c32 <- makenewbargraph("gonad", "male", "lo vs. hi PRL   ", 0, 2500) + labs(x = "lo vs. hi PRL")   + 
      theme(axis.title.y = element_blank(), axis.text = element_blank())
    c122232 <- plot_grid(c12, c22,c32, nrow = 3, rel_heights = c(1.1,1,1.1)) 


    hypothesisbars <-  plot_grid(c112131, c122232, b112131, b122232, nrow = 1, rel_widths = c(1.3, 1, 1.1, 1),
              labels = c("f", " ", "g", " "), label_size = 8)

    dataplots <- plot_grid(bcde, hypothesisbars)


    fig3 <- plot_grid(expdesign, dataplots, rel_heights = c(0.25,1), nrow = 2, labels = c(" ", "e"), label_size = 8)
    fig3

![](../figures/fig3-1.png)

    #write.csv(PRLvsd3,"../results/PRLvsd.csv", row.names = F)

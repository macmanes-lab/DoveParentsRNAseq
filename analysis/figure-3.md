all things box plots
====================

candidate genes
---------------

    library(tidyverse)

    ## ── Attaching packages ─────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(ggtext)
    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    library(ggpubr)

    ## Loading required package: magrittr

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     set_names

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

    ## 
    ## Attaching package: 'ggpubr'

    ## The following object is masked from 'package:cowplot':
    ## 
    ##     get_legend

    source("../R/themes.R")

    knitr::opts_chunk$set(echo = TRUE, message = F, fig.path = "../figures/")

    candidategenes <- c("OXT", "AVP", "GNRH1", "GNRHR", "CGNRH-R",
                        "AR", "POMC", "AGRP", 
                           "CRH", "AVPR1A", "AVPR1B", "AVPR2","VIP",
                           "CYP19A1", "DRD1", "DRD2", "PRL", "PRLR", "SOX9", 
                        "ESR1","ESR2", "LBH", "CDK1", "BRCA1",
                        "PTEN", "CREBBP", "FOS", "JUN", "EGR1",
                         "BDNF", "GRM2","GRIA1",
                        "KCNJ5", "CISH", "PTGER3", "CEBPD", "ZBTB16", 
                        "DIO3", "DIO2", "DIO1",
                        
                         "VIPR1", "VIPR2", "NPY", "LH", "FSH", "FSHB",
                        "JAK2", "HCRT", "TRH", "TSHB",
                        "MC3R", "MC4R", "MC5R",  "NR3C1", "NR3C2",
                        "STAT5B","NPVF") 

    candidategenesslim <- c("OXT", "AVP", "GNRH1", "GNRHR", 
                        "AR",  "CYP19A1", 
                         "AVPR1A", "AVPR1B", "AVPR2","VIP",
                      "DRD1", "DRD2", 
                      "PRL", "PRLR",  
                        "ESR1","ESR2", "LBH",  
                       
                        "FOS", "JUN", "EGR1", "BDNF"
                          ) 



    # all candidate genes: AR, AVP, AVPR1A, AVPR1B, AVPR2, BDNF, CYP19A1, DRD1, EGR1, ESR1, ESR2, FOS, GNRH1, GNRHR, JUN, LBH, OXT, PRL, PRLR, VIP


    # DEG candidate genes: AR, AVPR1A, AVPR1B, AVPR2, BDNF, CYP19A1, DRD1, ESR1, FOS, GNRHR, LBH, PRL, PRLR, 

    # Candidate genes that were not differentially expressed: AVP, AVPR1B, ESR2,  GNRH1, JUN, OXT, VIP

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

    hypvsd <- getcandidatevsd(candidategenesslim, "hypothalamus", sexlevels)
    pitvsd <- getcandidatevsd(candidategenesslim, "pituitary", sexlevels)
    gonvsd <- getcandidatevsd(candidategenesslim, "gonad", sexlevels)
    head(hypvsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue      treatment gene  samples                         counts
    ##   <chr>  <chr>       <fct>     <chr> <chr>                            <dbl>
    ## 1 female hypothalam… control   AR    L.G118_female_hypothalamus_con…   7.28
    ## 2 female hypothalam… control   AR    R.G106_female_hypothalamus_con…   7.32
    ## 3 female hypothalam… control   AR    R.R20_female_hypothalamus_cont…   7.74
    ## 4 female hypothalam… control   AR    R.R9_female_hypothalamus_contr…   7.31
    ## 5 female hypothalam… control   AR    R.W44_female_hypothalamus_cont…   7.02
    ## 6 female hypothalam… inc.d9    AR    blk.s061.pu.y_female_hypothala…   7.46

    head(pitvsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue    treatment gene  samples                           counts
    ##   <chr>  <chr>     <fct>     <chr> <chr>                              <dbl>
    ## 1 female pituitary control   AR    L.G118_female_pituitary_control.…   9.61
    ## 2 female pituitary control   AR    R.G106_female_pituitary_control     9.41
    ## 3 female pituitary control   AR    R.R20_female_pituitary_control      8.57
    ## 4 female pituitary control   AR    R.R9_female_pituitary_control.NY…   8.85
    ## 5 female pituitary control   AR    R.W44_female_pituitary_control.N…   9.42
    ## 6 female pituitary inc.d9    AR    blk.s061.pu.y_female_pituitary_i…   9.75

    head(gonvsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue treatment gene  samples                           counts
    ##   <chr>  <chr>  <fct>     <chr> <chr>                              <dbl>
    ## 1 female gonad  control   AR    L.G118_female_gonad_control         8.76
    ## 2 female gonad  control   AR    R.G106_female_gonad_control         9.27
    ## 3 female gonad  control   AR    R.R20_female_gonad_control          8.17
    ## 4 female gonad  control   AR    R.R9_female_gonad_control           8.29
    ## 5 female gonad  control   AR    R.W44_female_gonad_control          8.83
    ## 6 female gonad  inc.d9    AR    blk.s061.pu.y_female_gonad_inc.d9   9.06

    candidatevsd <- rbind(hypvsd, pitvsd)
    candidatevsd <- rbind(candidatevsd, gonvsd)

    unique(candidatevsd$gene)

    ##  [1] "AR"      "AVP"     "AVPR1A"  "AVPR1B"  "AVPR2"   "BDNF"    "CYP19A1"
    ##  [8] "DRD1"    "EGR1"    "ESR1"    "ESR2"    "FOS"     "GNRH1"   "GNRHR"  
    ## [15] "JUN"     "LBH"     "OXT"     "PRL"     "PRLR"    "VIP"

Figs
----

    makeboxplots <- function(df, whichgene, mysubtitle, whichsex, whichtimepoint){
      p <- df %>%
        filter(treatment %in% charlevels,
               gene %in% whichgene,
               sex %in% whichsex) %>%
        filter(treatment %in% whichtimepoint) %>%
        ggplot(aes(x = treatment, y = counts, fill = treatment, color = sex)) +
        geom_boxplot(outlier.shape = NA) + 
        geom_jitter(size = 0.5, width = 0.1) +
        facet_wrap(~gene, scales = "free_y", nrow = 1) +
        theme_B3() +
        scale_fill_manual(values = allcolors) +
        scale_color_manual(values = allcolors) +
        theme(#axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none") +
        labs(y = "gene expression" , x = "Sequential stages with candidate DEGs", subtitle = mysubtitle) +
        theme(strip.text = element_text(face = "italic"))
      return(p)
    }

    ## hyp 
    a <- makeboxplots(hypvsd, c("BDNF","CYP19A1", "DRD1", "EGR1"), 
                      "Female hypothalamus", "female", c("hatch", "n5"))  +
      theme(axis.title.x = element_blank())

    b <- makeboxplots(hypvsd, c("AR"),"Male hypothalamus", "male", 
                      c( "inc.d9", "inc.d17"))  +
      theme(axis.title = element_blank())

    ab <- plot_grid(a,b, rel_widths = c(4,1),  labels = "auto", label_size = 8)

    c <- makeboxplots(pitvsd, c("ESR1","GNRHR"), "Female pituitary", "female", c("bldg", "lay", "inc.d3"))   +
      theme(axis.title.x = element_blank())
    d <- makeboxplots(pitvsd, c("LBH", "PRL" ), " ", "female", c("inc.d9", "inc.d17", "hatch", "n5"))  +
      theme(axis.title  = element_blank())
    e <- makeboxplots(pitvsd, c( "AVPR2"), " ", "female", c( "hatch", "n5"))  +
      theme(axis.title  = element_blank())

    cde <- plot_grid(c,e,d, rel_widths = c(6,2.25,7), nrow = 1,  labels = c("c"), label_size = 8)


    g <- makeboxplots(gonvsd, c( "AVPR1A"), "Female gonads", "female", c( "lay", "inc.d3", "inc.d9"))  +
      labs(caption = " ", x = " ")
    h <- makeboxplots(gonvsd, c( "EGR1", "FOS", "PRLR"), " ", "female", c( "lay", "inc.d3"))  +
      theme(axis.title.y = element_blank()) + labs(caption = " ")
    f <- makeboxplots(pitvsd, c("LBH", "PRL"), "Male pituitary", "male", 
                      c("inc.d9", "inc.d17")) + 
      theme(plot.caption = element_text(face = "italic")) + 
      labs(x = " ",
           caption = "Candidate genes that were not differentially expressed from bldg to n5: AVP, AVPR1B, ESR2,  GNRH1, JUN, OXT, VIP.")

    ghf <- plot_grid(g,h,f, nrow = 1, rel_widths = c(3,6,4),  labels = c("d", " ", "e"), label_size = 8)

    plot_grid(ab, cde, ghf, nrow = 3, rel_heights = c(1,1,1.25))

![](../figures/fig3-1.png)

    write.csv(candidatevsd, "../../musicalgenes/data/candidatecounts.csv")

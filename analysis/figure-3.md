Correlations in gene expression
===============================

    library(tidyverse)

    ## ── Attaching packages ─────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    library(corrr)
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
    source("../R/functions.R")
    source("../R/wrangledata.R")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_double(),
    ##   Name = col_character(),
    ##   geneid = col_double(),
    ##   entrezid = col_character()
    ## )

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    knitr::opts_chunk$set(echo = TRUE, message = F, fig.path = "../figures/")

data wrangle
------------

    # load `candidatevsd` with `source("../R/wrangledata.R")`
    head(candidatevsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue      treatment gene  samples                         counts
    ##   <chr>  <chr>       <fct>     <chr> <chr>                            <dbl>
    ## 1 female hypothalam… control   ADRA… L.G118_female_hypothalamus_con…   8.87
    ## 2 female hypothalam… control   ADRA… R.G106_female_hypothalamus_con…   8.73
    ## 3 female hypothalam… control   ADRA… R.R20_female_hypothalamus_cont…   9.11
    ## 4 female hypothalam… control   ADRA… R.R9_female_hypothalamus_contr…   8.63
    ## 5 female hypothalam… control   ADRA… R.W44_female_hypothalamus_cont…   9.17
    ## 6 female hypothalam… inc.d9    ADRA… blk.s061.pu.y_female_hypothala…   9.04

    curleychampagnegenes <- curleychampagnegenes %>% pull(gene)

    # widen for correlations 

    candidatevsdwide <- candidatevsd  %>%
        pivot_wider(names_from = gene, values_from = counts) 
    head(candidatevsdwide)

    ## # A tibble: 6 x 125
    ##   sex   tissue treatment samples ADRA2A   AHR ANXA1    AR  AREG ATP2C2
    ##   <chr> <chr>  <fct>     <chr>    <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
    ## 1 fema… hypot… control   L.G118…   8.87  9.61  6.97  7.28  6.67   10.7
    ## 2 fema… hypot… control   R.G106…   8.73  9.43  6.49  7.32  7.14   10.8
    ## 3 fema… hypot… control   R.R20_…   9.11  9.26  8.43  7.74  7.13   10.6
    ## 4 fema… hypot… control   R.R9_f…   8.63  9.53  6.69  7.31  7.02   10.9
    ## 5 fema… hypot… control   R.W44_…   9.17  9.56  7.99  7.02  6.98   10.9
    ## 6 fema… hypot… inc.d9    blk.s0…   9.04  9.73  6.12  7.46  6.94   10.7
    ## # … with 115 more variables: AVP <dbl>, AVPR1A <dbl>, BMP4 <dbl>,
    ## #   BMP7 <dbl>, BRCA2 <dbl>, BRINP1 <dbl>, BTRC <dbl>, CCND1 <dbl>,
    ## #   CD44 <dbl>, CDKN1B <dbl>, CEBPB <dbl>, CHUK <dbl>, COMT <dbl>,
    ## #   CREBRF <dbl>, CRH <dbl>, CRHBP <dbl>, CRHR1 <dbl>, CRHR2 <dbl>,
    ## #   CSF1 <dbl>, CTNNB1 <dbl>, CYP19A1 <dbl>, CYP7B1 <dbl>, DBH <dbl>,
    ## #   DEAF1 <dbl>, DRD1 <dbl>, DRD4 <dbl>, EAF2 <dbl>, EPHA2 <dbl>,
    ## #   ERBB4 <dbl>, ESR1 <dbl>, ESR2 <dbl>, ETV5 <dbl>, FEM1B <dbl>,
    ## #   FGF10 <dbl>, FGF2 <dbl>, FGFR2 <dbl>, FKBP4 <dbl>, FOS <dbl>,
    ## #   FOXA1 <dbl>, FRS2 <dbl>, GLI2 <dbl>, GNAQ <dbl>, HIF1A <dbl>,
    ## #   HTR2C <dbl>, ID2 <dbl>, ID4 <dbl>, IGF1 <dbl>, IGF1R <dbl>,
    ## #   JAK2 <dbl>, KALRN <dbl>, LATS1 <dbl>, LBH <dbl>, LRP5 <dbl>,
    ## #   LRP6 <dbl>, MAPK1 <dbl>, MBD2 <dbl>, MED1 <dbl>, MEST <dbl>,
    ## #   MMP2 <dbl>, MST1 <dbl>, MSX1 <dbl>, MSX2 <dbl>, NCOA3 <dbl>,
    ## #   `NKX3-1` <dbl>, NOG <dbl>, NPAS3 <dbl>, NR3C1 <dbl>, NTN1 <dbl>,
    ## #   OPRK1 <dbl>, OPRM1 <dbl>, ORAI1 <dbl>, OXT <dbl>, PGR <dbl>,
    ## #   PHB2 <dbl>, PLAG1 <dbl>, PML <dbl>, PRL <dbl>, PRLR <dbl>, PSAP <dbl>,
    ## #   PTCH1 <dbl>, PTEN <dbl>, PTHLH <dbl>, PYGO2 <dbl>, RLN1 <dbl>,
    ## #   ROBO1 <dbl>, RREB1 <dbl>, RTN4 <dbl>, RXRA <dbl>, SCRIB <dbl>,
    ## #   SERPINB5 <dbl>, SERPINF1 <dbl>, SFRP1 <dbl>, SHH <dbl>, SLC12A2 <dbl>,
    ## #   SLC6A4 <dbl>, SLIT2 <dbl>, SMO <dbl>, SOSTDC1 <dbl>, SOX9 <dbl>,
    ## #   SRC <dbl>, …

    makecorrdf <- function(whichsex, whichtissue, whichgenes){
     
      corrrdf <- candidatevsd %>%
        filter(sex == whichsex, tissue == whichtissue,
               gene %in% whichgenes) %>%
        pivot_wider(names_from = gene, values_from = counts) %>%
        select(-sex,-tissue, -treatment, -samples) %>%
        correlate() %>%
        rearrange()
      print(head(corrrdf))
      return(corrrdf)
    }


    c1 <- makecorrdf("female", "hypothalamus", curleychampagnegenes)  

    ## # A tibble: 6 x 23
    ##   rowname  HTR2C   DRD1  CRHBP  CRHR2   MEST ADRA2A AVPR1A CRHR1   CRH
    ##   <chr>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>
    ## 1 HTR2C   NA      0.851  0.795  0.835  0.492  0.411  0.400 0.551 0.346
    ## 2 DRD1     0.851 NA      0.778  0.730  0.521  0.401  0.439 0.353 0.356
    ## 3 CRHBP    0.795  0.778 NA      0.756  0.583  0.470  0.464 0.362 0.266
    ## 4 CRHR2    0.835  0.730  0.756 NA      0.589  0.344  0.401 0.539 0.499
    ## 5 MEST     0.492  0.521  0.583  0.589 NA      0.376  0.361 0.356 0.335
    ## 6 ADRA2A   0.411  0.401  0.470  0.344  0.376 NA      0.338 0.257 0.183
    ## # … with 13 more variables: OPRM1 <dbl>, ESR1 <dbl>, PGR <dbl>, AVP <dbl>,
    ## #   FOS <dbl>, OXT <dbl>, NR3C1 <dbl>, ESR2 <dbl>, DRD4 <dbl>,
    ## #   SLC6A4 <dbl>, PRLR <dbl>, PRL <dbl>, COMT <dbl>

    c2 <- makecorrdf("female", "pituitary", curleychampagnegenes)  

    ## # A tibble: 6 x 22
    ##   rowname   NR3C1   ADRA2A   AVPR1A    ESR1    ESR2   CRHR1   DRD1  SLC6A4
    ##   <chr>     <dbl>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>  <dbl>   <dbl>
    ## 1 NR3C1   NA       0.398    0.110    0.159   0.340   0.0968 0.324   0.294 
    ## 2 ADRA2A   0.398  NA        0.00466  0.0909  0.153  -0.163  0.204   0.227 
    ## 3 AVPR1A   0.110   0.00466 NA        0.355  -0.0372  0.735  0.0625 -0.0823
    ## 4 ESR1     0.159   0.0909   0.355   NA      -0.0485  0.0993 0.0927  0.0187
    ## 5 ESR2     0.340   0.153   -0.0372  -0.0485 NA       0.0363 0.0739  0.0248
    ## 6 CRHR1    0.0968 -0.163    0.735    0.0993  0.0363 NA      0.0734  0.0147
    ## # … with 13 more variables: MEST <dbl>, FOS <dbl>, PRL <dbl>, CRHR2 <dbl>,
    ## #   OPRM1 <dbl>, PGR <dbl>, DRD4 <dbl>, HTR2C <dbl>, PRLR <dbl>,
    ## #   COMT <dbl>, CRHBP <dbl>, AVP <dbl>, OXT <dbl>

    c3 <- makecorrdf("female", "gonad", curleychampagnegenes)  

    ## # A tibble: 6 x 23
    ##   rowname    PGR    ESR1  OPRM1   COMT   MEST     AVP   PRLR    PRL  ADRA2A
    ##   <chr>    <dbl>   <dbl>  <dbl>  <dbl>  <dbl>   <dbl>  <dbl>  <dbl>   <dbl>
    ## 1 PGR     NA      0.866   0.369  0.307  0.323  0.163  0.455  0.220   0.535 
    ## 2 ESR1     0.866 NA       0.269  0.213  0.210  0.0618 0.429  0.187   0.587 
    ## 3 OPRM1    0.369  0.269  NA      0.312  0.330  0.379  0.0995 0.204  -0.0102
    ## 4 COMT     0.307  0.213   0.312 NA      0.216  0.308  0.304  0.0823  0.0887
    ## 5 MEST     0.323  0.210   0.330  0.216 NA      0.279  0.0254 0.169  -0.173 
    ## 6 AVP      0.163  0.0618  0.379  0.308  0.279 NA      0.158  0.219  -0.152 
    ## # … with 13 more variables: NR3C1 <dbl>, OXT <dbl>, CRH <dbl>,
    ## #   AVPR1A <dbl>, SLC6A4 <dbl>, DRD4 <dbl>, FOS <dbl>, CRHR2 <dbl>,
    ## #   HTR2C <dbl>, DRD1 <dbl>, CRHR1 <dbl>, CRHBP <dbl>, ESR2 <dbl>

    c4 <- makecorrdf("male", "hypothalamus", curleychampagnegenes)  

    ## # A tibble: 6 x 23
    ##   rowname  HTR2C   DRD1  CRHR2  CRHBP    CRH  CRHR1  DRD4  MEST  ESR1   FOS
    ##   <chr>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 HTR2C   NA      0.831  0.777  0.759  0.578  0.513 0.479 0.351 0.236 0.106
    ## 2 DRD1     0.831 NA      0.724  0.668  0.475  0.461 0.575 0.351 0.140 0.162
    ## 3 CRHR2    0.777  0.724 NA      0.742  0.336  0.441 0.352 0.329 0.171 0.206
    ## 4 CRHBP    0.759  0.668  0.742 NA      0.425  0.376 0.385 0.426 0.264 0.143
    ## 5 CRH      0.578  0.475  0.336  0.425 NA      0.601 0.198 0.175 0.488 0.195
    ## 6 CRHR1    0.513  0.461  0.441  0.376  0.601 NA     0.346 0.187 0.270 0.366
    ## # … with 12 more variables: NR3C1 <dbl>, OPRM1 <dbl>, ADRA2A <dbl>,
    ## #   AVPR1A <dbl>, PGR <dbl>, PRLR <dbl>, SLC6A4 <dbl>, AVP <dbl>,
    ## #   ESR2 <dbl>, OXT <dbl>, PRL <dbl>, COMT <dbl>

    c5 <- makecorrdf("male", "pituitary", curleychampagnegenes)  

    ## # A tibble: 6 x 22
    ##   rowname  ADRA2A     OXT    PRLR    COMT    OPRM1      AVP    DRD4    ESR2
    ##   <chr>     <dbl>   <dbl>   <dbl>   <dbl>    <dbl>    <dbl>   <dbl>   <dbl>
    ## 1 ADRA2A  NA       0.172  -0.0533 -0.0310  0.411    0.0213  -0.198   0.0827
    ## 2 OXT      0.172  NA       0.100   0.154   0.0562   0.618    0.0378 -0.127 
    ## 3 PRLR    -0.0533  0.100  NA       0.431   0.0204   0.0125   0.516   0.0330
    ## 4 COMT    -0.0310  0.154   0.431  NA       0.0848   0.0128   0.0659  0.147 
    ## 5 OPRM1    0.411   0.0562  0.0204  0.0848 NA        0.00625 -0.281   0.0835
    ## 6 AVP      0.0213  0.618   0.0125  0.0128  0.00625 NA        0.0832 -0.258 
    ## # … with 13 more variables: ESR1 <dbl>, CRHR2 <dbl>, SLC6A4 <dbl>,
    ## #   PGR <dbl>, CRHBP <dbl>, MEST <dbl>, DRD1 <dbl>, PRL <dbl>,
    ## #   NR3C1 <dbl>, HTR2C <dbl>, FOS <dbl>, CRHR1 <dbl>, AVPR1A <dbl>

    c6 <- makecorrdf("male", "gonad", curleychampagnegenes)  

    ## # A tibble: 6 x 23
    ##   rowname     ESR1    FOS  ADRA2A    DRD1   HTR2C   SLC6A4    DRD4   MEST
    ##   <chr>      <dbl>  <dbl>   <dbl>   <dbl>   <dbl>    <dbl>   <dbl>  <dbl>
    ## 1 ESR1    NA        0.427  0.264   0.257   0.135  -0.00566  0.127  0.0499
    ## 2 FOS      0.427   NA      0.377   0.161   0.216   0.193    0.0816 0.0519
    ## 3 ADRA2A   0.264    0.377 NA       0.0863 -0.119   0.0421   0.0870 0.140 
    ## 4 DRD1     0.257    0.161  0.0863 NA       0.0688 -0.0200   0.103  0.0404
    ## 5 HTR2C    0.135    0.216 -0.119   0.0688 NA       0.0991  -0.140  0.154 
    ## 6 SLC6A4  -0.00566  0.193  0.0421 -0.0200  0.0991 NA        0.0757 0.117 
    ## # … with 14 more variables: OXT <dbl>, PGR <dbl>, ESR2 <dbl>, COMT <dbl>,
    ## #   NR3C1 <dbl>, CRH <dbl>, AVPR1A <dbl>, CRHBP <dbl>, PRL <dbl>,
    ## #   AVP <dbl>, OPRM1 <dbl>, PRLR <dbl>, CRHR1 <dbl>, CRHR2 <dbl>

    subsetcandidatevsdwide <- function(whichsex, whichtissue){
      df <- candidatevsdwide %>%
        filter(sex == whichsex, tissue == whichtissue,
               treatment != "control") 
      return(df)
    }

    FH <- subsetcandidatevsdwide("female", "hypothalamus")
    FP <- subsetcandidatevsdwide("female", "pituitary")
    FG <- subsetcandidatevsdwide("female", "gonad")
    MH <- subsetcandidatevsdwide("male", "hypothalamus")
    MP <- subsetcandidatevsdwide("male", "pituitary")
    MG <- subsetcandidatevsdwide("male", "gonad") 

### figure

    plotcorrplot <- function(df, subtitle){  
      p <- df %>%
        rplot2() +
        theme_B3() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              axis.text = element_text(face = "italic"))

      return(p)
    }

    p1 <- plotcorrplot(c1, "FH") + theme(legend.position = "none", axis.text.x = element_blank()) + labs(subtitle = "Hypothalamus")
    p2 <- plotcorrplot(c2, "FP") + theme(legend.position = "none", axis.text.x = element_blank()) + 
      labs(subtitle = "Pituitary") 
    p3 <- plotcorrplot(c3, "FG") + theme(legend.position = "none") + labs(subtitle = "Gonad")  

    p4 <- plotcorrplot(c4, "MH") + theme(legend.position = "none", axis.text.y = element_blank(), axis.text.x = element_blank()) + labs(subtitle = " ")
    p5 <- plotcorrplot(c5, "MP") + theme(legend.position = "none", axis.text.y = element_blank(), axis.text.x = element_blank())  +
      labs(subtitle = " ") 
    p6 <- plotcorrplot(c6, " ") + theme(legend.position = "none", axis.text.y = element_blank()) + labs(subtitle = " ") 

    a <- plot_grid(p1 + labs(title = "Female"),
                   p4 + labs(title = "Male"),
                   p2, p5, p3, p6, nrow = 3, rel_widths = c(1.2,1),
                   rel_heights = c(1,1,1.2))

    scattercorrelations <- function(df, gene1, gene2, mylinecolor, myxlab, myylab){
      p <- ggplot(df, aes(x = gene1, y = gene2)) +
        labs(subtitle = " ", x = myxlab, y = myylab) + 
        scale_color_manual(values = allcolors) +
        geom_smooth(method = "lm",  color = mylinecolor) +
        geom_point(aes(color = treatment))  +
        theme_B3() +  
        theme(axis.title = element_text(face = "italic"), 
                legend.position = "none") +
        stat_cor( size = 2)
      return(p)
    }


    # "female" = "#969696", "male" = "#525252"

    p7 <- scattercorrelations(FH, FH$DRD1, FH$HTR2C, "#969696", "DRD1", "HTR2C") + labs(subtitle = "Hypothalamus")
    p8 <- scattercorrelations(FP, FP$AVPR1A, FP$CRHR1, "#969696", "AVPR1A", "CRHR1")  + labs(subtitle = "Pituitary")
    p9 <- scattercorrelations(FG, FG$ESR2, FG$CRHBP, "#969696", "ESR2", "CRHBP")   + labs(subtitle = "Gonad")

    p10 <- scattercorrelations(MH, MH$DRD1, MH$HTR2C, "#525252", "DRD1", "HTR2C") + theme(axis.title.y = element_blank())
    p11 <- scattercorrelations(MP, MP$AVPR1A, MP$CRHR1, "#525252", "AVPR1A", "CRHR1") + theme(axis.title.y = element_blank())
    p12 <- scattercorrelations(MG, MG$ESR2, MG$CRHBP, "#525252", "ESR2", "CRHBP")  + theme(axis.title.y = element_blank())

    b <- plot_grid(p7 + labs(title = "Female"), 
                   p10 + labs(title = "Male"), 
                   p8, p11, p9, p12, nrow = 3, rel_heights = c(1,1,1.2))

    ab <- plot_grid(a,b, rel_widths = c(1.5,1), labels = c("A", "B"), label_size = 8)


    p13 <- plotcorrplot(c1, "FH") + theme(legend.position = "bottom")
    p14 <- scattercorrelations(FH, FH$DRD1, FH$HTR2C, "#969696", "DRD1", "HTR2C") + 
      theme(legend.position = "bottom",
            legend.title = element_blank())

    legend1 <- p13 %>% get_legend()
    legend2 <- p14 %>% get_legend()

    legends <- plot_grid(legend1, legend2, rel_widths = c(1.5,1))

    plot_grid(ab, legends, nrow = 2, rel_heights = c(1,0.1))

![](../figures/fig3-1.png)

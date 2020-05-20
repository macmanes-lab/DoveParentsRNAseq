Candidate gene analysis
=======================

    library(tidyverse)
    library(ggtext)
    library(cowplot)
    library(ggpubr)
    library(knitr)
    library(kableExtra)
    library(corrr)
    library(ggsignif)
    library(magick)
    library(scales)
    library(ggimage)


    source("../R/themes.R")
    source("../R/functions.R")

    knitr::opts_chunk$set(echo = TRUE, message = F, fig.path = "../figures/")

Candidate genes
---------------

    parentalcaregenes <- read_csv("../metadata/03_parentalcaregenes.csv") %>% select(-X1)

    ## Warning: Missing column names filled in: 'X1' [1]

    head(parentalcaregenes)

    ## # A tibble: 6 x 5
    ##   geneid NCBI           literature GO     gene  
    ##    <dbl> <chr>          <chr>      <chr>  <chr> 
    ## 1 428980 XP_004942333.2 ADRA2A     <NA>   ADRA2A
    ## 2 396101 NP_990516.1    AVP        AVP    AVP   
    ## 3 771773 NP_001103908.1 AVPR1A     AVPR1A AVPR1A
    ## 4 395098 NP_989780.1    <NA>       BRINP1 BRINP1
    ## 5 416783 XP_001233014.1 COMT       <NA>   COMT  
    ## 6 416206 XP_001231574.1 <NA>       CREBRF CREBRF

    curleychampagnegenes <- parentalcaregenes %>% distinct(literature) %>% drop_na() %>% pull(literature)
    GOgenes <- parentalcaregenes %>% distinct(GO)  %>% drop_na() %>% pull(GO)
    candidategenes <- parentalcaregenes %>% pull(gene) 
    candidategenes

    ##  [1] "ADRA2A" "AVP"    "AVPR1A" "BRINP1" "COMT"   "CREBRF" "CRH"   
    ##  [8] "CRHBP"  "CRHR1"  "CRHR2"  "DBH"    "DRD1"   "DRD4"   "ESR1"  
    ## [15] "ESR2"   "FOS"    "GNAQ"   "HTR2C"  "KALRN"  "MBD2"   "MEST"  
    ## [22] "NPAS3"  "NPAS3"  "NR3C1"  "OPRK1"  "OPRM1"  "OXT"    "PGR"   
    ## [29] "PRL"    "PRLR"   "PTEN"   "SLC6A4" "ZFX"

Candidate DEGs
--------------

    # summary DEG results from DESeq2
    candidateDEGS <- read_csv("../results/suppltable1.csv") %>%
      filter(gene %in% candidategenes) %>%
      mutate(posneg = ifelse(lfc >= 0, "+", "-"),
             sex = recode(sex, "female" = "F", "male" = "M" ),
             tissue = recode(tissue, 
                             "hypothalamus" = "H",
                             "pituitary" = "P", "gonad" = "G")) %>%
      mutate(res = paste(sex, tissue, posneg, sep = "")) %>%
      select(gene, res, comparison)  %>%
      group_by(gene,  comparison) %>%
      summarize(res = str_c(res, collapse = " ")) %>%
      pivot_wider(names_from = comparison, values_from = res) %>%
      select(gene, bldg_lay, lay_inc.d3, inc.d3_inc.d9, inc.d9_inc.d17, hatch_n5, n5_n9)
    candidateDEGS

    ## # A tibble: 23 x 7
    ## # Groups:   gene [23]
    ##    gene   bldg_lay lay_inc.d3 inc.d3_inc.d9 inc.d9_inc.d17 hatch_n5 n5_n9
    ##    <chr>  <chr>    <chr>      <chr>         <chr>          <chr>    <chr>
    ##  1 ADRA2A <NA>     <NA>       <NA>          MH+            <NA>     <NA> 
    ##  2 AVP    <NA>     <NA>       <NA>          MH+            <NA>     FG+  
    ##  3 AVPR1A <NA>     FG+        FG-           <NA>           <NA>     <NA> 
    ##  4 BRINP1 <NA>     <NA>       <NA>          FG+            <NA>     <NA> 
    ##  5 COMT   <NA>     <NA>       <NA>          <NA>           FH-      <NA> 
    ##  6 CREBRF FG+      FG-        <NA>          <NA>           FP+      <NA> 
    ##  7 CRHBP  <NA>     <NA>       <NA>          <NA>           FH+      <NA> 
    ##  8 CRHR2  <NA>     <NA>       <NA>          <NA>           FH+      <NA> 
    ##  9 DRD1   <NA>     <NA>       <NA>          <NA>           FH+      <NA> 
    ## 10 DRD4   <NA>     FP-        <NA>          <NA>           FH+      <NA> 
    ## # … with 13 more rows

    ## table 1 summary candidate genes
    table1 <- left_join(candidateDEGS, parentalcaregenes) %>%
      select(gene, bldg_lay:n5_n9, literature, GO, NCBI) %>%
      mutate(literature = if_else(is.na(literature), " ", "X"),
             GO = if_else(is.na(GO), " ", "X"))
    table1[is.na(table1)] <- " " # replace NA with blank space so it's pretty
    head(table1)

    ## # A tibble: 6 x 10
    ## # Groups:   gene [6]
    ##   gene  bldg_lay lay_inc.d3 inc.d3_inc.d9 inc.d9_inc.d17 hatch_n5 n5_n9
    ##   <chr> <chr>    <chr>      <chr>         <chr>          <chr>    <chr>
    ## 1 ADRA… " "      " "        " "           "MH+"          " "      " "  
    ## 2 AVP   " "      " "        " "           "MH+"          " "      "FG+"
    ## 3 AVPR… " "      "FG+"      "FG-"         " "            " "      " "  
    ## 4 BRIN… " "      " "        " "           "FG+"          " "      " "  
    ## 5 COMT  " "      " "        " "           " "            "FH-"    " "  
    ## 6 CREB… "FG+"    "FG-"      " "           " "            "FP+"    " "  
    ## # … with 3 more variables: literature <chr>, GO <chr>, NCBI <chr>

Candidate VSDs
--------------

    # load `candidatevsd` with `source("../R/wrangledata.R")`
    candidatevsd <- read_csv("../results/03_candidatevsd.csv") %>% 
      select(-X1) %>%
      filter(treatment %in% charlevels)

    ## Warning: Missing column names filled in: 'X1' [1]

    tail(candidatevsd)

    ## # A tibble: 6 x 6
    ##   sex   tissue treatment gene  samples                         counts
    ##   <chr> <chr>  <chr>     <chr> <chr>                            <dbl>
    ## 1 male  gonad  n9        ZFX   y129.x_male_gonad_n9              10.3
    ## 2 male  gonad  n9        ZFX   y131.w185.x_male_gonad_n9         10.4
    ## 3 male  gonad  inc.d17   ZFX   y133.w77.r58_male_gonad_inc.d17   10.3
    ## 4 male  gonad  inc.d3    ZFX   y149.r52.x_male_gonad_inc.d3      10.5
    ## 5 male  gonad  inc.d9    ZFX   y95.g131.x_male_gonad_inc.d9      10.5
    ## 6 male  gonad  inc.d3    ZFX   y98.o50.x_male_gonad_inc.d3       10.2

    candidatevsdwide <- candidatevsd  %>%
        pivot_wider(names_from = gene, values_from = counts) 
    FH <- subsetcandidatevsdwide("female", "hypothalamus")
    FP <- subsetcandidatevsdwide("female", "pituitary")
    FG <- subsetcandidatevsdwide("female", "gonad")
    MH <- subsetcandidatevsdwide("male", "hypothalamus")
    MP <- subsetcandidatevsdwide("male", "pituitary")
    MG <- subsetcandidatevsdwide("male", "gonad") 

Candidate Correlations - SUppl fig 1
------------------------------------

    hyp1 <- makecorrdf("female", "hypothalamus", curleychampagnegenes)  

    ## # A tibble: 6 x 23
    ##   rowname ADRA2A    AVP AVPR1A   COMT    CRH  CRHBP   CRHR1  CRHR2   DRD1
    ##   <chr>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>  <dbl>  <dbl>
    ## 1 ADRA2A  NA      0.413  0.334 -0.452  0.184  0.464  0.257   0.339  0.396
    ## 2 AVP      0.413 NA      0.300 -0.110  0.353  0.144  0.0892  0.151  0.158
    ## 3 AVPR1A   0.334  0.300 NA     -0.356  0.292  0.460  0.251   0.398  0.436
    ## 4 COMT    -0.452 -0.110 -0.356 NA     -0.336 -0.732 -0.373  -0.640 -0.855
    ## 5 CRH      0.184  0.353  0.292 -0.336 NA      0.262  0.387   0.498  0.352
    ## 6 CRHBP    0.464  0.144  0.460 -0.732  0.262 NA      0.358   0.755  0.776
    ## # … with 13 more variables: DRD4 <dbl>, ESR1 <dbl>, ESR2 <dbl>, FOS <dbl>,
    ## #   HTR2C <dbl>, MEST <dbl>, NR3C1 <dbl>, OPRM1 <dbl>, OXT <dbl>,
    ## #   PGR <dbl>, PRL <dbl>, PRLR <dbl>, SLC6A4 <dbl>

    pit1 <- makecorrdf("female", "pituitary", curleychampagnegenes)  

    ## # A tibble: 6 x 22
    ##   rowname   ADRA2A    AVP   AVPR1A   COMT  CRHBP  CRHR1   CRHR2    DRD1
    ##   <chr>      <dbl>  <dbl>    <dbl>  <dbl>  <dbl>  <dbl>   <dbl>   <dbl>
    ## 1 ADRA2A  NA       -0.267  0.00324 -0.110 -0.154 -0.163 0.00357  0.203 
    ## 2 AVP     -0.267   NA     -0.202    0.484  0.543 -0.165 0.0195  -0.109 
    ## 3 AVPR1A   0.00324 -0.202 NA       -0.120 -0.200  0.739 0.0330   0.0636
    ## 4 COMT    -0.110    0.484 -0.120   NA      0.401 -0.241 0.108   -0.0232
    ## 5 CRHBP   -0.154    0.543 -0.200    0.401 NA     -0.177 0.101   -0.0707
    ## 6 CRHR1   -0.163   -0.165  0.739   -0.241 -0.177 NA     0.0809   0.0738
    ## # … with 13 more variables: DRD4 <dbl>, ESR1 <dbl>, ESR2 <dbl>, FOS <dbl>,
    ## #   HTR2C <dbl>, MEST <dbl>, NR3C1 <dbl>, OPRM1 <dbl>, OXT <dbl>,
    ## #   PGR <dbl>, PRL <dbl>, PRLR <dbl>, SLC6A4 <dbl>

    gon1 <- makecorrdf("female", "gonad", curleychampagnegenes)  

    ## # A tibble: 6 x 23
    ##   rowname  ADRA2A     AVP  AVPR1A    COMT     CRH   CRHBP   CRHR1    CRHR2
    ##   <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>
    ## 1 ADRA2A  NA      -0.152  -0.102   0.0897 -0.0290 -0.172  -0.316   2.28e-5
    ## 2 AVP     -0.152  NA      -0.111   0.308  -0.0685 -0.373  -0.287  -3.43e-1
    ## 3 AVPR1A  -0.102  -0.111  NA       0.0102  0.116   0.0483 -0.0593  1.36e-1
    ## 4 COMT     0.0897  0.308   0.0102 NA       0.0224 -0.369  -0.331  -2.00e-1
    ## 5 CRH     -0.0290 -0.0685  0.116   0.0224 NA      -0.0264 -0.0286 -1.32e-1
    ## 6 CRHBP   -0.172  -0.373   0.0483 -0.369  -0.0264 NA       0.755   7.28e-1
    ## # … with 14 more variables: DRD1 <dbl>, DRD4 <dbl>, ESR1 <dbl>,
    ## #   ESR2 <dbl>, FOS <dbl>, HTR2C <dbl>, MEST <dbl>, NR3C1 <dbl>,
    ## #   OPRM1 <dbl>, OXT <dbl>, PGR <dbl>, PRL <dbl>, PRLR <dbl>, SLC6A4 <dbl>

    hyp2 <- makecorrdf("male", "hypothalamus", curleychampagnegenes)  

    ## # A tibble: 6 x 23
    ##   rowname ADRA2A     AVP AVPR1A   COMT    CRH   CRHBP   CRHR1   CRHR2
    ##   <chr>    <dbl>   <dbl>  <dbl>  <dbl>  <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 ADRA2A  NA      0.369   0.328 -0.313  0.306  0.169   0.261   0.133 
    ## 2 AVP      0.369 NA       0.344 -0.214  0.252 -0.0639  0.145  -0.112 
    ## 3 AVPR1A   0.328  0.344  NA     -0.219  0.211  0.271   0.0491  0.0158
    ## 4 COMT    -0.313 -0.214  -0.219 NA     -0.663 -0.407  -0.639  -0.466 
    ## 5 CRH      0.306  0.252   0.211 -0.663 NA      0.430   0.600   0.342 
    ## 6 CRHBP    0.169 -0.0639  0.271 -0.407  0.430 NA       0.381   0.743 
    ## # … with 14 more variables: DRD1 <dbl>, DRD4 <dbl>, ESR1 <dbl>,
    ## #   ESR2 <dbl>, FOS <dbl>, HTR2C <dbl>, MEST <dbl>, NR3C1 <dbl>,
    ## #   OPRM1 <dbl>, OXT <dbl>, PGR <dbl>, PRL <dbl>, PRLR <dbl>, SLC6A4 <dbl>

    pit2 <- makecorrdf("male", "pituitary", curleychampagnegenes)  

    ## # A tibble: 6 x 23
    ##   rowname   ADRA2A     AVP  AVPR1A    COMT      CRH   CRHBP    CRHR1
    ##   <chr>      <dbl>   <dbl>   <dbl>   <dbl>    <dbl>   <dbl>    <dbl>
    ## 1 ADRA2A  NA        0.0226 -0.252  -0.0329  0.00379 -0.0130 -0.345  
    ## 2 AVP      0.0226  NA      -0.100   0.0156  0.173    0.0232 -0.0955 
    ## 3 AVPR1A  -0.252   -0.100  NA      -0.120  -0.0780   0.104   0.765  
    ## 4 COMT    -0.0329   0.0156 -0.120  NA      -0.0466   0.0426 -0.137  
    ## 5 CRH      0.00379  0.173  -0.0780 -0.0466 NA        0.289  -0.00907
    ## 6 CRHBP   -0.0130   0.0232  0.104   0.0426  0.289   NA       0.0998 
    ## # … with 15 more variables: CRHR2 <dbl>, DRD1 <dbl>, DRD4 <dbl>,
    ## #   ESR1 <dbl>, ESR2 <dbl>, FOS <dbl>, HTR2C <dbl>, MEST <dbl>,
    ## #   NR3C1 <dbl>, OPRM1 <dbl>, OXT <dbl>, PGR <dbl>, PRL <dbl>, PRLR <dbl>,
    ## #   SLC6A4 <dbl>

    gon2 <- makecorrdf("male", "gonad", curleychampagnegenes) 

    ## # A tibble: 6 x 23
    ##   rowname  ADRA2A     AVP  AVPR1A    COMT     CRH   CRHBP   CRHR1   CRHR2
    ##   <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 ADRA2A  NA      -0.165  -0.0284  0.104  -0.0520 -0.182  -0.0646 -0.116 
    ## 2 AVP     -0.165  NA      -0.0463  0.121   0.108   0.105   0.119   0.162 
    ## 3 AVPR1A  -0.0284 -0.0463 NA      -0.104  -0.0361  0.0318  0.0645  0.184 
    ## 4 COMT     0.104   0.121  -0.104  NA       0.116   0.0501  0.193   0.0447
    ## 5 CRH     -0.0520  0.108  -0.0361  0.116  NA       0.175   0.0575  0.368 
    ## 6 CRHBP   -0.182   0.105   0.0318  0.0501  0.175  NA       0.0388  0.337 
    ## # … with 14 more variables: DRD1 <dbl>, DRD4 <dbl>, ESR1 <dbl>,
    ## #   ESR2 <dbl>, FOS <dbl>, HTR2C <dbl>, MEST <dbl>, NR3C1 <dbl>,
    ## #   OPRM1 <dbl>, OXT <dbl>, PGR <dbl>, PRL <dbl>, PRLR <dbl>, SLC6A4 <dbl>

    hyp3 <- makecorrdf("female", "hypothalamus", GOgenes)  

    ## # A tibble: 6 x 17
    ##   rowname     AVP  AVPR1A  BRINP1  CREBRF     DBH    DRD1    GNAQ   KALRN
    ##   <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 AVP     NA       0.300   0.0751 -0.138  -0.0363  0.158  -0.0370 -0.0529
    ## 2 AVPR1A   0.300  NA       0.193   0.0303  0.0341  0.436   0.263   0.220 
    ## 3 BRINP1   0.0751  0.193  NA       0.0123 -0.0441  0.144  -0.235   0.0666
    ## 4 CREBRF  -0.138   0.0303  0.0123 NA      -0.123   0.0291 -0.0655  0.531 
    ## 5 DBH     -0.0363  0.0341 -0.0441 -0.123  NA      -0.0363  0.0105 -0.0238
    ## 6 DRD1     0.158   0.436   0.144   0.0291 -0.0363 NA       0.570   0.548 
    ## # … with 8 more variables: MBD2 <dbl>, NPAS3 <dbl>, NR3C1 <dbl>,
    ## #   OPRK1 <dbl>, OXT <dbl>, PRL <dbl>, PTEN <dbl>, ZFX <dbl>

    pit3 <- makecorrdf("female", "pituitary", GOgenes)  

    ## # A tibble: 6 x 17
    ##   rowname     AVP  AVPR1A  BRINP1 CREBRF     DBH    DRD1   GNAQ   KALRN
    ##   <chr>     <dbl>   <dbl>   <dbl>  <dbl>   <dbl>   <dbl>  <dbl>   <dbl>
    ## 1 AVP     NA      -0.202   0.271   0.136 -0.0172 -0.109  0.0581  0.0494
    ## 2 AVPR1A  -0.202  NA       0.158  -0.218  0.0422  0.0636 0.0611  0.0223
    ## 3 BRINP1   0.271   0.158  NA       0.102  0.0636  0.292  0.109  -0.0817
    ## 4 CREBRF   0.136  -0.218   0.102  NA     -0.109   0.221  0.0461  0.0633
    ## 5 DBH     -0.0172  0.0422  0.0636 -0.109 NA       0.0598 0.0465  0.159 
    ## 6 DRD1    -0.109   0.0636  0.292   0.221  0.0598 NA      0.169  -0.104 
    ## # … with 8 more variables: MBD2 <dbl>, NPAS3 <dbl>, NR3C1 <dbl>,
    ## #   OPRK1 <dbl>, OXT <dbl>, PRL <dbl>, PTEN <dbl>, ZFX <dbl>

    gon3 <- makecorrdf("female", "gonad", GOgenes)  

    ## # A tibble: 6 x 17
    ##   rowname    AVP   AVPR1A   BRINP1   CREBRF     DBH    DRD1    GNAQ
    ##   <chr>    <dbl>    <dbl>    <dbl>    <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 AVP     NA     -0.111   -0.341   -0.141    0.107  -0.415   0.216 
    ## 2 AVPR1A  -0.111 NA       -0.00646  0.0406  -0.0923  0.109  -0.0345
    ## 3 BRINP1  -0.341 -0.00646 NA       -0.00823 -0.0136  0.673  -0.665 
    ## 4 CREBRF  -0.141  0.0406  -0.00823 NA        0.123   0.203   0.267 
    ## 5 DBH      0.107 -0.0923  -0.0136   0.123   NA      -0.0631  0.139 
    ## 6 DRD1    -0.415  0.109    0.673    0.203   -0.0631 NA      -0.574 
    ## # … with 9 more variables: KALRN <dbl>, MBD2 <dbl>, NPAS3 <dbl>,
    ## #   NR3C1 <dbl>, OPRK1 <dbl>, OXT <dbl>, PRL <dbl>, PTEN <dbl>, ZFX <dbl>

    hyp4 <- makecorrdf("male", "hypothalamus", GOgenes)  

    ## # A tibble: 6 x 17
    ##   rowname      AVP  AVPR1A  BRINP1  CREBRF     DBH     DRD1     GNAQ  KALRN
    ##   <chr>      <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>    <dbl>  <dbl>
    ## 1 AVP     NA        0.344   0.165   0.191  -0.0205  0.00251 -0.103   0.0633
    ## 2 AVPR1A   0.344   NA       0.0754  0.161   0.121   0.230    0.220   0.118 
    ## 3 BRINP1   0.165    0.0754 NA       0.157   0.0512  0.486    0.116   0.211 
    ## 4 CREBRF   0.191    0.161   0.157  NA      -0.0682  0.349    0.00434 0.659 
    ## 5 DBH     -0.0205   0.121   0.0512 -0.0682 NA       0.204    0.165   0.163 
    ## 6 DRD1     0.00251  0.230   0.486   0.349   0.204  NA        0.538   0.652 
    ## # … with 8 more variables: MBD2 <dbl>, NPAS3 <dbl>, NR3C1 <dbl>,
    ## #   OPRK1 <dbl>, OXT <dbl>, PRL <dbl>, PTEN <dbl>, ZFX <dbl>

    pit4 <- makecorrdf("male", "pituitary", GOgenes)  

    ## # A tibble: 6 x 17
    ##   rowname     AVP  AVPR1A  BRINP1  CREBRF     DBH    DRD1    GNAQ    KALRN
    ##   <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>
    ## 1 AVP     NA      -0.100   0.267  -0.0207  0.330  -0.0301 -0.231   0.0183 
    ## 2 AVPR1A  -0.100  NA       0.120  -0.0383 -0.216   0.0697 -0.0175 -0.0286 
    ## 3 BRINP1   0.267   0.120  NA       0.170   0.237   0.0928 -0.159  -0.261  
    ## 4 CREBRF  -0.0207 -0.0383  0.170  NA       0.181   0.115  -0.106   0.204  
    ## 5 DBH      0.330  -0.216   0.237   0.181  NA      -0.0880 -0.221  -0.00964
    ## 6 DRD1    -0.0301  0.0697  0.0928  0.115  -0.0880 NA       0.0782  0.213  
    ## # … with 8 more variables: MBD2 <dbl>, NPAS3 <dbl>, NR3C1 <dbl>,
    ## #   OPRK1 <dbl>, OXT <dbl>, PRL <dbl>, PTEN <dbl>, ZFX <dbl>

    gon4 <- makecorrdf("male", "gonad", GOgenes)  

    ## # A tibble: 6 x 17
    ##   rowname     AVP  AVPR1A  BRINP1  CREBRF     DBH    DRD1     GNAQ   KALRN
    ##   <chr>     <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl>   <dbl>
    ## 1 AVP     NA      -0.0463 -0.203   0.252  -0.0800 -0.0761  0.280   -0.0528
    ## 2 AVPR1A  -0.0463 NA       0.0153  0.125   0.353  -0.0576  0.127    0.153 
    ## 3 BRINP1  -0.203   0.0153 NA      -0.0125  0.188   0.128  -0.119   -0.0148
    ## 4 CREBRF   0.252   0.125  -0.0125 NA       0.126   0.151   0.255    0.107 
    ## 5 DBH     -0.0800  0.353   0.188   0.126  NA       0.216   0.140    0.137 
    ## 6 DRD1    -0.0761 -0.0576  0.128   0.151   0.216  NA      -0.00687  0.118 
    ## # … with 8 more variables: MBD2 <dbl>, NPAS3 <dbl>, NR3C1 <dbl>,
    ## #   OPRK1 <dbl>, OXT <dbl>, PRL <dbl>, PTEN <dbl>, ZFX <dbl>

    b1 <- plotcorrplot(hyp1, "females") + labs(y = "Hypothalamus", title = "Parental Care Literature")  + 
      scale_x_discrete(position = "top")   + theme(axis.text.x = element_text(vjust = -0.25))
    b2 <- plotcorrplot(pit1, NULL)  + labs(y = "Pituitary") + 
      theme(axis.text.x = element_blank())
    b3 <- plotcorrplot(gon1, NULL)  + labs(y = "Gonad") + theme(axis.text.x = element_text(vjust = 0.25))
     
    b4 <- plotcorrplot(hyp2, "males") + labs( title =  " ") + theme( axis.text.y = element_blank(), axis.text.x = element_text(vjust = -0.25)) +
      scale_x_discrete(position = "top") 
    b5 <- plotcorrplot(pit2, NULL)    + 
      theme(axis.text.x = element_blank(), axis.text.y = element_blank())
    b6 <- plotcorrplot(gon2, NULL)  + theme( axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.25))


    b7 <- plotcorrplot(hyp3, "females") + labs( title = "Parental Care GO") + 
       scale_x_discrete(position = "top") + theme(axis.text.x = element_text(vjust = -0.25))
    b8 <- plotcorrplot(pit3, NULL)   + 
      theme(axis.text.x = element_blank())
    b9 <- plotcorrplot(gon3, NULL)  + theme(axis.text.x = element_text(vjust = 0.25))

    b10 <- plotcorrplot(hyp4, "males") + labs( title =  " ") + 
      theme(axis.text.y = element_blank(),
            axis.text.x = element_text(vjust = -0.25)) + scale_x_discrete(position = "top") 
    b11 <- plotcorrplot(pit4, NULL)    + 
      theme(axis.text.x = element_blank() , axis.text.y = element_blank())
    b12 <- plotcorrplot(gon4, NULL)   + theme( axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.25))

    forlegend <- plotcorrplot(gon4, NULL) + theme(legend.position = "bottom") 
    mylegend <- get_legend(forlegend)

    allcorplots <- plot_grid(b1,b4,b7,b10,
              b2,b5,b8,b11,
              b3,b6,b9,b12, rel_heights = c(1.4,1,1.2),
              labels = c("A", " ", "B"), label_size = 8, rel_widths = c(1.4,1.2,1.2,1))

    supplfig2 <- plot_grid(allcorplots, mylegend, ncol = 1, rel_heights = c(1,0.1))
    supplfig2

![](../figures/supplfig-2-1.png)

    pdf(file="../figures/supplfig-2.pdf", width=7.25, height=7.25)
    plot(supplfig2)
    dev.off()

    ## quartz_off_screen 
    ##                 2

Figure
------

    a <- png::readPNG("../figures/images/fig_fig2a.png")
    a <- ggdraw() +  draw_image(a, scale = 1)




    c1 <- scattercorrelations(FH, FH$DRD4, FH$HTR2C, "#969696", "DRD4", "HTR2C") + labs(subtitle = "female", title =  "Top correlation")  
    c3 <- scattercorrelations(FP, FP$AVPR1A, FP$CRHR1, "#969696", "AVPR1A", "CRHR1")  + labs(subtitle = NULL) 
    c5 <- scattercorrelations(FG, FG$ESR1, FG$PGR, "#969696", "ESR1", "PGR")   + labs(subtitle = NULL) 

    c2 <- scattercorrelations(MH, MH$DRD4, MH$HTR2C, "#525252", "DRD4", "HTR2C")  + labs(subtitle = "male", title = "", y = NULL)  
    c4 <- scattercorrelations(MP, MP$AVPR1A, MP$CRHR1, "#525252", "AVPR1A", "CRHR1")  + labs(subtitle = NULL, y = NULL)  
    c6 <- scattercorrelations(MG, MG$ESR1, MG$PGR, "#525252", "ESR1", "PGR")  + labs(subtitle = NULL, y = NULL)  


    d1 <- candidateboxplot("hypothalamus", c("DRD4"), sexlevels) + labs(x = NULL, title = "Temporal expression pattern" )  
    d2 <- candidateboxplot("hypothalamus", c("HTR2C"), sexlevels) + labs(x = NULL ) + theme(strip.text = element_blank())
    d3 <- candidateboxplot("pituitary", c("AVPR1A"), sexlevels) + labs(x = NULL )  + theme( strip.text = element_blank())
    d4 <- candidateboxplot("pituitary", c("CRHR1"), sexlevels) + labs(x = NULL  )  + theme(strip.text = element_blank())
    d5 <- candidateboxplot("gonad", c("ESR1"), sexlevels) + labs(x = NULL) + theme(strip.text = element_blank())
    d6 <- candidateboxplot("gonad", c("PGR"), sexlevels) + labs(x = NULL ) + theme(strip.text = element_blank())

    #b <- plot_grid(b1,b2,b3, ncol = 1, rel_heights = c(1.1,1,1))
    c <- plot_grid(c1,c2,c3,c4,c5,c6, ncol = 2, rel_heights = c(1.1,1,1), rel_widths = c(1.1,1))
    d <- plot_grid(d1,d2,d3,d4,d5,d6, ncol = 1, rel_heights = c(1.2,1, 1,1, 1,1))

    bcd <- plot_grid(c,d, ncol = 2, rel_widths = c(1,1), labels = c("A", "B"), label_size = 8)

    fig2 <- plot_grid(bcd, a, ncol = 1, rel_heights = c(1,0.1))
    fig2

![](../figures/fig2-1.png)

    pdf(file="../figures/fig2-1.pdf", width=7.25, height=7.25)
    plot(fig2)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    #write.csv(candidatevsd, "../../musicalgenes/data/candidatecounts.csv")
    #write.csv(candidatevsd, "../results/candidatecounts.csv")
    write.csv(table1, "../results/table1.csv")

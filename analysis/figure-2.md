Candidate gene analysis
=======================

    library(tidyverse)
    library(ggtext)
    library(cowplot)
    library(ggpubr)
    library(knitr)

    ## Warning: package 'knitr' was built under R version 3.6.2

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

    candidategenes <- parentalcaregenes %>% pull(gene) 
    candidategenes

    ##  [1] "ADRA2A"  "AVP"     "AVPR1A"  "BRINP1"  "CGNRH-R" "COMT"    "CREBRF" 
    ##  [8] "CRH"     "CRHBP"   "CRHR1"   "CRHR2"   "DBH"     "DRD1"    "DRD4"   
    ## [15] "ESR1"    "ESR2"    "FOS"     "FSHB"    "FSHR"    "GAL"     "GNAQ"   
    ## [22] "GNRH1"   "GNRHR"   "HTR2C"   "KALRN"   "MBD2"    "MEST"    "NPAS3"  
    ## [29] "NPAS3"   "NR3C1"   "OPRK1"   "OPRM1"   "OXT"     "PGR"     "PRL"    
    ## [36] "PRLR"    "PTEN"    "SLC6A4"  "TH"      "THRB"    "VIP"     "ZFX"

    GOgenes <- parentalcaregenes %>% distinct(GO)  %>% drop_na() %>% pull(GO)
    GOgenes

    ##  [1] "AVP"    "AVPR1A" "BRINP1" "CREBRF" "DBH"    "DRD1"   "GNAQ"   "KALRN" 
    ##  [9] "MBD2"   "NPAS3"  "NR3C1"  "OPRK1"  "OXT"    "PRL"    "PTEN"   "ZFX"

    literaturegenes <- parentalcaregenes %>% 
      distinct(literature) %>% 
      drop_na() %>% 
      pull(literature)
    literaturegenes

    ##  [1] "ADRA2A"  "AVP"     "AVPR1A"  "CGNRH-R" "COMT"    "CRH"     "CRHBP"  
    ##  [8] "CRHR1"   "CRHR2"   "DRD1"    "DRD4"    "ESR1"    "ESR2"    "FOS"    
    ## [15] "FSHB"    "FSHR"    "GAL"     "GNRH1"   "GNRHR"   "HTR2C"   "MEST"   
    ## [22] "NR3C1"   "OPRM1"   "OXT"     "PGR"     "PRL"     "PRLR"    "SLC6A4" 
    ## [29] "TH"      "THRB"    "VIP"

    litNotGO <-  literaturegenes[!literaturegenes %in% GOgenes]
    litNotGOknown <- litNotGO[!litNotGO %in% c("ADRA2A", "MEST", "SLC6A4")]
    litNotGOknown

    ##  [1] "CGNRH-R" "COMT"    "CRH"     "CRHBP"   "CRHR1"   "CRHR2"   "DRD4"   
    ##  [8] "ESR1"    "ESR2"    "FOS"     "FSHB"    "FSHR"    "GAL"     "GNRH1"  
    ## [15] "GNRHR"   "HTR2C"   "OPRM1"   "PGR"     "PRLR"    "TH"      "THRB"   
    ## [22] "VIP"

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
      select(gene, lay_inc.d3, inc.d3_inc.d9, inc.d9_inc.d17, hatch_n5, n5_n9)
    candidateDEGS

    ## # A tibble: 28 x 6
    ## # Groups:   gene [28]
    ##    gene   lay_inc.d3 inc.d3_inc.d9 inc.d9_inc.d17 hatch_n5 n5_n9
    ##    <chr>  <chr>      <chr>         <chr>          <chr>    <chr>
    ##  1 ADRA2A <NA>       <NA>          MH+            <NA>     <NA> 
    ##  2 AVP    <NA>       <NA>          MH+            <NA>     FG+  
    ##  3 AVPR1A FG+        FG-           <NA>           <NA>     <NA> 
    ##  4 BRINP1 <NA>       <NA>          FG+            <NA>     <NA> 
    ##  5 COMT   <NA>       <NA>          <NA>           FH-      <NA> 
    ##  6 CREBRF FG-        <NA>          <NA>           FP+      <NA> 
    ##  7 CRHBP  <NA>       <NA>          <NA>           FH+      <NA> 
    ##  8 CRHR2  <NA>       <NA>          <NA>           FH+      <NA> 
    ##  9 DRD1   <NA>       <NA>          <NA>           FH+      <NA> 
    ## 10 DRD4   FP-        <NA>          <NA>           FH+      <NA> 
    ## # … with 18 more rows

    ## table 1 summary candidate genes
    table1 <- left_join(candidateDEGS, parentalcaregenes) %>%
      select(gene, lay_inc.d3:n5_n9, literature, GO, NCBI) %>%
      mutate(literature = if_else(is.na(literature), " ", "X"),
             GO = if_else(is.na(GO), " ", "X"))
    (table1)

    ## # A tibble: 28 x 9
    ## # Groups:   gene [28]
    ##    gene  lay_inc.d3 inc.d3_inc.d9 inc.d9_inc.d17 hatch_n5 n5_n9 literature GO   
    ##    <chr> <chr>      <chr>         <chr>          <chr>    <chr> <chr>      <chr>
    ##  1 ADRA… <NA>       <NA>          MH+            <NA>     <NA>  "X"        " "  
    ##  2 AVP   <NA>       <NA>          MH+            <NA>     FG+   "X"        "X"  
    ##  3 AVPR… FG+        FG-           <NA>           <NA>     <NA>  "X"        "X"  
    ##  4 BRIN… <NA>       <NA>          FG+            <NA>     <NA>  " "        "X"  
    ##  5 COMT  <NA>       <NA>          <NA>           FH-      <NA>  "X"        " "  
    ##  6 CREB… FG-        <NA>          <NA>           FP+      <NA>  " "        "X"  
    ##  7 CRHBP <NA>       <NA>          <NA>           FH+      <NA>  "X"        " "  
    ##  8 CRHR2 <NA>       <NA>          <NA>           FH+      <NA>  "X"        " "  
    ##  9 DRD1  <NA>       <NA>          <NA>           FH+      <NA>  "X"        "X"  
    ## 10 DRD4  FP-        <NA>          <NA>           FH+      <NA>  "X"        " "  
    ## # … with 18 more rows, and 1 more variable: NCBI <chr>

Candidate VSDs
--------------

    # load `candidatevsd` with `source("../R/wrangledata.R")`
    candidatevsd <- read_csv("../results/03_candidatevsd.csv") %>% 
      select(-X1) %>%
      filter(treatment %in% charlevels) %>%
      mutate(treatment = factor(treatment, levels = charlevels)) %>%
      drop_na()

    ## Warning: Missing column names filled in: 'X1' [1]

    head(candidatevsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue       treatment gene  samples                             counts
    ##   <chr>  <chr>        <fct>     <chr> <chr>                                <dbl>
    ## 1 female hypothalamus control   ADRA… L.G118_female_hypothalamus_control…   8.95
    ## 2 female hypothalamus control   ADRA… R.G106_female_hypothalamus_control    8.81
    ## 3 female hypothalamus control   ADRA… R.R20_female_hypothalamus_control.…   9.18
    ## 4 female hypothalamus control   ADRA… R.R9_female_hypothalamus_control      8.72
    ## 5 female hypothalamus control   ADRA… R.W44_female_hypothalamus_control     9.23
    ## 6 female hypothalamus inc.d9    ADRA… blk.s061.pu.y_female_hypothalamus_…   9.11

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

    hyp1 <- makecorrdf("female", "hypothalamus", litNotGOknown)  
    pit1 <- makecorrdf("female", "pituitary", litNotGOknown)  
    gon1 <- makecorrdf("female", "gonad", litNotGOknown)  

    hyp2 <- makecorrdf("male", "hypothalamus", litNotGOknown)  
    pit2 <- makecorrdf("male", "pituitary", litNotGOknown)  
    gon2 <- makecorrdf("male", "gonad", litNotGOknown) 


    hyp3 <- makecorrdf("female", "hypothalamus", GOgenes)  
    pit3 <- makecorrdf("female", "pituitary", GOgenes)  
    gon3 <- makecorrdf("female", "gonad", GOgenes)  

    hyp4 <- makecorrdf("male", "hypothalamus", GOgenes)  
    pit4 <- makecorrdf("male", "pituitary", GOgenes)  
    gon4 <- makecorrdf("male", "gonad", GOgenes)  

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

    forlegend <- plotcorrplot(gon4, NULL) + theme(legend.position = "bottom") + labs(color = "Correlation")
    mylegend <- get_legend(forlegend)

    allcorplots <- plot_grid(b1,b4,b7,b10,
              b2,b5,b8,b11,
              b3,b6,b9,b12, rel_heights = c(1.4,1,1.2),
              labels = c("A", " ", "B"), label_size = 8, rel_widths = c(1.4,1.2,1.2,1))

    supplfig2 <- plot_grid(allcorplots, mylegend, ncol = 1, rel_heights = c(1,0.1))
    supplfig2

![](../figures/supplfig-2-1.png)

Figure
------

    c1 <- scattercorrelations(FH, FH$DRD1, "DRD1", FH$HTR2C, "HTR2C", "#969696" ) +  labs(title = " ", subtitle = " " )  
    c3 <- scattercorrelations(FP, FP$AVPR1A, "AVPR1A", FP$CRHR1,  "CRHR1", "#969696")  + labs(title = " ", subtitle = " ") 
    c5 <- scattercorrelations(FG, FG$PGR, "PGR",  FG$ESR1,  "ESR1", "#969696")   + labs(title = " ", subtitle = " ") 

    c2 <- scattercorrelations(MH, MH$DRD1, "DRD1", MH$HTR2C, "HTR2C", "#525252")  + labs(title = " ", subtitle = " " )  
    c4 <- scattercorrelations(MP, MP$AVPR1A,  "AVPR1A",  MP$CRHR1, "CRHR1", "#525252")  + labs(title = " ", subtitle = " " ) 
    c6 <- scattercorrelations(MG, MG$PGR, "PGR", MG$ESR1, "ESR1", "#525252")   + labs(title = " ", subtitle = " " ) 

    d1 <- candidateboxplot("hypothalamus", c("DRD1"), "female") + labs(x = NULL, title =  "Hypothlamic expression", subtitle = "females" )  
    d2 <- candidateboxplot("hypothalamus", c("HTR2C"), "female") + labs(x = NULL )  + theme(axis.text.x = element_text(angle = 45, vjust = 1))
    d3 <- candidateboxplot("pituitary", c("AVPR1A"), "female") + labs(x = NULL , title = "Pituitary expression", subtitle = "females") 
    d4 <- candidateboxplot("pituitary", c("CRHR1"), "female") + labs(x = NULL  )  + theme(axis.text.x = element_text(angle = 45, vjust = 1))
    d5 <- candidateboxplot("gonad", c("PGR"), "female") + labs(x = NULL, title = "Gonadal expression", subtitle = "females")
    d6 <- candidateboxplot("gonad", c("ESR1"), "female") + labs(x = NULL ) + theme(axis.text.x = element_text(angle = 45, vjust = 1))

    e1 <- candidateboxplot("hypothalamus", c("DRD1"), "male") + labs(x = NULL, title = " " , subtitle = "males" )  
    e2 <- candidateboxplot("hypothalamus", c("HTR2C"), "male") + labs(x = NULL, title = " " )  + theme(axis.text.x = element_text(angle = 45, vjust = 1))
    e3 <- candidateboxplot("pituitary", c("AVPR1A"), "male") + labs(x = NULL , subtitle = "males", title = " " ) 
    e4 <- candidateboxplot("pituitary", c("CRHR1"), "male") + labs(x = NULL, title = " "  )  + theme(axis.text.x = element_text(angle = 45, vjust = 1))
    e5 <- candidateboxplot("gonad", c("PGR"), "male") + labs(x = NULL, subtitle = "males", title = " " ) 
    e6 <- candidateboxplot("gonad", c("ESR1"), "male") + labs(x = NULL , title = " ")   + theme(axis.text.x = element_text(angle = 45, vjust = 1)) 

    #b <- plot_grid(b1,b2,b3, ncol = 1, rel_heights = c(1.1,1,1))
    c123 <- plot_grid(c1,c3,c5, ncol = 1, rel_heights = c(1.1,1,1))
    c456 <- plot_grid(c2,c4,c6, ncol = 1, rel_heights = c(1.1,1,1))
    d <- plot_grid(d1,d2,d3,d4,d5,d6, ncol = 1, rel_heights = c(1.2,1, 1,1, 1,1), labels = c("A", " ",  "B", " ", "C"), label_size = 8)
    e <- plot_grid(e1,e2,e3,e4,e5,e6, ncol = 1, rel_heights = c(1.2,1, 1,1, 1,1))

    fig2 <- plot_grid(d, c123,e,c456,nrow  = 1)
    fig2

![](../figures/fig2-1.png)

    a <- candidateboxplot("hypothalamus", c("COMT"), "female") + 
      labs(title =  "Sequential", 
           subtitle = "Female hypothalamus" ) 

    b <- externalboxplots("hypothalamus", c("COMT"), "female") + 
      labs(title =  "External", 
           subtitle = " " ) 

    c <- candidateboxplot("pituitary", c("PRL"), "female") + 
      labs(subtitle = "Female pituitary" ) 

    d <- externalboxplots("pituitary", c("PRL"), "female") + 
       labs(subtitle = " " ) 

    e <- candidateboxplot("gonad", c("AVPR1A"), "female") + 
      labs( subtitle = "Female gonads" ) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    f <- externalboxplots("pituitary", c("AVPR1A"), "female") + 
       labs( subtitle = " " ) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    females <- plot_grid(a,b,c,d,e,f, ncol = 2, rel_widths = c(2,1))

    a2 <- candidateboxplot("hypothalamus", c("COMT"), "male") + 
      labs(title =  "Sequential", 
           subtitle = "Male hypothalamus" ) 

    b2 <- externalboxplots("hypothalamus", c("COMT"), "male") + 
      labs(title =  "External", 
           subtitle = "  " ) 

    c2 <- candidateboxplot("pituitary", c("PRL"), "male") + 
      labs(subtitle = "Male pituitary" ) 

    d2 <- externalboxplots("pituitary", c("PRL"), "male") + 
       labs(subtitle = " " ) 

    e2 <- candidateboxplot("gonad", c("AVPR1A"), "male") + 
      labs( subtitle = "Male gonads" ) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    f2 <- externalboxplots("pituitary", c("AVPR1A"), "male") + 
       labs( subtitle = " " )  + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    males <- plot_grid(a2,b2,c2,d2,e2,f2, ncol = 2, rel_widths = c(2,1))

    ## Warning in wilcox.test.default(c(5.74459037236268, 5.69190747435115,
    ## 5.57627376474667, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(5.71639281686726, 5.45090433098895,
    ## 5.72401881911303, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(5.64556372505326, 5.05272358504361,
    ## 5.53492367908641, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(5.05272358504361, 5.72204002358928,
    ## 5.58427891014697, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(5.63670921041675, 5.57616171696985,
    ## 5.05272358504361, : cannot compute exact p-value with ties

    ## Warning in wilcox.test.default(c(5.54670900060151, 5.68124741136788,
    ## 5.90118219899707, : cannot compute exact p-value with ties

    newfig2 <- plot_grid(females, males, labels = c("A", "B"),
                         label_size  = 8)
    newfig2

![](../figures/fig2-2.png)

    #write.csv(candidatevsd, "../../musicalgenes/data/candidatecounts.csv")
    #write.csv(candidatevsd, "../results/candidatecounts.csv")
    write.csv(table1, "../results/table1-v2.csv")


    pdf(file="../figures/fig2-1.pdf", width=7.25, height=7.25)
    plot(fig2)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/supplfig-2.pdf", width=7.25, height=7.25)
    plot(supplfig2)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/fig2-2.pdf", width=7.25, height=7.25)
    plot(newfig2)
    dev.off()

    ## quartz_off_screen 
    ##                 2

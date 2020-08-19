Table 1. Candidate genes and differentially expressed genes
===========================================================

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

Differentially expressed genes (DEGs)
-------------------------------------

    # summary DEG results from DESeq2
    allDEGS <- read_csv("../results/suppltable1.csv") %>%
       mutate(posneg = ifelse(lfc >= 0, "+", "-"),
             sex = recode(sex, "female" = "F", "male" = "M" ),
             tissue = recode(tissue, 
                             "hypothalamus" = "H",
                             "pituitary" = "P", "gonad" = "G"),
             comparison = recode(comparison, 
                             "bldg_lay" = "Bldg2Lay",
                             "lay_inc.d3" = "Lay2Inc3", 
                             "inc.d3_inc.d9" = "Inc3toInc9",
                             "inc.d9_inc.d17" = "Inc9toInc17",
                             "inc.d17_hatch" = "In17toHatch",
                             "hatch_n5" = "Hatch2N5",
                             "n5_n9" = "N5toN9"),
             group = paste(sex, tissue, sep = "")) %>%
      mutate(res = paste("(", posneg, ")", sep = "")) %>%
      mutate(compres = paste(comparison, res, sep = "")) %>%
      select(gene, group, compres) 
    allDEGS

    ## # A tibble: 14,795 x 3
    ##    gene    group compres    
    ##    <chr>   <chr> <chr>      
    ##  1 HBG2    FH    Bldg2Lay(-)
    ##  2 HEMGN   FH    Bldg2Lay(-)
    ##  3 BLB1    FH    Lay2Inc3(-)
    ##  4 C1QA    FH    Lay2Inc3(-)
    ##  5 CFD     FH    Lay2Inc3(-)
    ##  6 CLIC2   FH    Lay2Inc3(-)
    ##  7 HLA-DRA FH    Lay2Inc3(-)
    ##  8 IRF1    FH    Lay2Inc3(-)
    ##  9 SCIN    FH    Lay2Inc3(-)
    ## 10 SLC35B2 FH    Lay2Inc3(-)
    ## # … with 14,785 more rows

    DEGsbytissue <- function(whichtissueshorthand, whichsex){
      
      df <- allDEGS %>%
      filter(grepl(whichtissueshorthand, group),
             grepl(whichsex,group)) %>%
      group_by(group,gene) %>%
      summarize(results = str_c(compres, collapse = "; ")) %>%
      mutate(n = (str_count(results, pattern = ";"))+1) %>%
      arrange(desc(n))  %>%
      filter(n>1)
      #group_by(comparisons) %>%
      #summarize(genes = str_c(gene, collapse = "; "))  

      print(head(df))
      return(df)
      
    }

    hmDEGs <- DEGsbytissue("H", "M")

    ## # A tibble: 6 x 4
    ## # Groups:   group [1]
    ##   group gene   results                                 n
    ##   <chr> <chr>  <chr>                               <dbl>
    ## 1 MH    TPH2   Bldg2Lay(-); Hatch2N5(+); N5toN9(-)     3
    ## 2 MH    ALPK2  Lay2Inc3(+); Inc9toInc17(-)             2
    ## 3 MH    ATP2A1 Hatch2N5(-); In17toHatch(+)             2
    ## 4 MH    CERKL  Inc9toInc17(-); In17toHatch(+)          2
    ## 5 MH    COX14  Bldg2Lay(+); Lay2Inc3(-)                2
    ## 6 MH    CSF3R  In17toHatch(-); Inc9toInc17(+)          2

    pmDEGs <- DEGsbytissue("P", "M")

    ## # A tibble: 6 x 4
    ## # Groups:   group [1]
    ##   group gene        results                                         n
    ##   <chr> <chr>       <chr>                                       <dbl>
    ## 1 MP    C11H19ORF40 Lay2Inc3(-); Inc3toInc9(+); Inc9toInc17(-)      3
    ## 2 MP    CREM        Inc9toInc17(+); Hatch2N5(-); N5toN9(-)          3
    ## 3 MP    ERLEC1      Inc9toInc17(+); Hatch2N5(-); In17toHatch(+)     3
    ## 4 MP    FTL         Inc9toInc17(+); Hatch2N5(-); In17toHatch(+)     3
    ## 5 MP    OLFM1       Inc9toInc17(+); In17toHatch(+); N5toN9(-)       3
    ## 6 MP    ADARB1      Inc9toInc17(+); Hatch2N5(-)                     2

    gmDEGs <- DEGsbytissue("G", "M")

    ## # A tibble: 6 x 4
    ## # Groups:   group [1]
    ##   group gene        results                      n
    ##   <chr> <chr>       <chr>                    <dbl>
    ## 1 MG    ANXA5       Bldg2Lay(-); Lay2Inc3(+)     2
    ## 2 MG    C28H19ORF10 Bldg2Lay(-); Lay2Inc3(+)     2
    ## 3 MG    CHD9        Bldg2Lay(+); Lay2Inc3(-)     2
    ## 4 MG    EIF3H       Bldg2Lay(-); Lay2Inc3(+)     2
    ## 5 MG    EIF3I       Bldg2Lay(-); Lay2Inc3(+)     2
    ## 6 MG    GLYR1       Bldg2Lay(-); Lay2Inc3(+)     2

    hfDEGs <- DEGsbytissue("H", "F")

    ## # A tibble: 6 x 4
    ## # Groups:   group [1]
    ##   group gene    results                                         n
    ##   <chr> <chr>   <chr>                                       <dbl>
    ## 1 FH    IGJ     In17toHatch(-); Inc9toInc17(+); Hatch2N5(+)     3
    ## 2 FH    ABCB5   Hatch2N5(-); In17toHatch(+)                     2
    ## 3 FH    ABHD16A Inc9toInc17(+); Hatch2N5(-)                     2
    ## 4 FH    ADAM9   Inc9toInc17(+); Hatch2N5(-)                     2
    ## 5 FH    ADGRV1  Inc9toInc17(-); Hatch2N5(+)                     2
    ## 6 FH    AGGF1   Inc9toInc17(-); Hatch2N5(+)                     2

    pfDEGs <- DEGsbytissue("P", "F")

    ## # A tibble: 6 x 4
    ## # Groups:   group [1]
    ##   group gene        results                                                   n
    ##   <chr> <chr>       <chr>                                                 <dbl>
    ## 1 FP    ADARB1      Bldg2Lay(-); Lay2Inc3(+); Inc9toInc17(+); Hatch2N5(-)     4
    ## 2 FP    ALG2        Bldg2Lay(-); Lay2Inc3(+); Inc9toInc17(+); Hatch2N5(-)     4
    ## 3 FP    ARF1        Bldg2Lay(-); Lay2Inc3(+); Inc9toInc17(+); Hatch2N5(-)     4
    ## 4 FP    ARF4        Bldg2Lay(-); Lay2Inc3(+); Inc9toInc17(+); Hatch2N5(-)     4
    ## 5 FP    C28H19ORF10 Bldg2Lay(-); Lay2Inc3(+); Inc9toInc17(+); Hatch2N5(-)     4
    ## 6 FP    CRELD2      Bldg2Lay(-); Lay2Inc3(+); Inc9toInc17(+); Hatch2N5(-)     4

    gfDEGs <- DEGsbytissue("G", "F")

    ## # A tibble: 6 x 4
    ## # Groups:   group [1]
    ##   group gene        results                                                    n
    ##   <chr> <chr>       <chr>                                                  <dbl>
    ## 1 FG    LOC1070534… Bldg2Lay(+); Lay2Inc3(-); In17toHatch(-); Inc9toInc17…     5
    ## 2 FG    COL10A1     Bldg2Lay(+); Lay2Inc3(-); Inc3toInc9(-); N5toN9(+)         4
    ## 3 FG    ANP32A      Bldg2Lay(-); Lay2Inc3(+); Inc9toInc17(-)                   3
    ## 4 FG    ANXA5       Bldg2Lay(-); Lay2Inc3(+); Inc9toInc17(-)                   3
    ## 5 FG    CA8         Bldg2Lay(-); Lay2Inc3(-); N5toN9(+)                        3
    ## 6 FG    CALB1       Bldg2Lay(+); Lay2Inc3(-); N5toN9(+)                        3

    maleDEGs <- rbind(hmDEGs, pmDEGs) %>%
      rbind(., gmDEGs)
    maleDEGs

    ## # A tibble: 168 x 4
    ## # Groups:   group [3]
    ##    group gene         results                                 n
    ##    <chr> <chr>        <chr>                               <dbl>
    ##  1 MH    TPH2         Bldg2Lay(-); Hatch2N5(+); N5toN9(-)     3
    ##  2 MH    ALPK2        Lay2Inc3(+); Inc9toInc17(-)             2
    ##  3 MH    ATP2A1       Hatch2N5(-); In17toHatch(+)             2
    ##  4 MH    CERKL        Inc9toInc17(-); In17toHatch(+)          2
    ##  5 MH    COX14        Bldg2Lay(+); Lay2Inc3(-)                2
    ##  6 MH    CSF3R        In17toHatch(-); Inc9toInc17(+)          2
    ##  7 MH    HAX1         Bldg2Lay(+); Lay2Inc3(-)                2
    ##  8 MH    LOC100858707 In17toHatch(-); Inc9toInc17(+)          2
    ##  9 MH    LOC101748402 Hatch2N5(-); In17toHatch(+)             2
    ## 10 MH    NEB          Hatch2N5(-); In17toHatch(+)             2
    ## # … with 158 more rows

    femaleDEGs <- rbind(hfDEGs, pfDEGs) %>%
      rbind(., gfDEGs)
    femaleDEGs

    ## # A tibble: 1,924 x 4
    ## # Groups:   group [3]
    ##    group gene    results                                         n
    ##    <chr> <chr>   <chr>                                       <dbl>
    ##  1 FH    IGJ     In17toHatch(-); Inc9toInc17(+); Hatch2N5(+)     3
    ##  2 FH    ABCB5   Hatch2N5(-); In17toHatch(+)                     2
    ##  3 FH    ABHD16A Inc9toInc17(+); Hatch2N5(-)                     2
    ##  4 FH    ADAM9   Inc9toInc17(+); Hatch2N5(-)                     2
    ##  5 FH    ADGRV1  Inc9toInc17(-); Hatch2N5(+)                     2
    ##  6 FH    AGGF1   Inc9toInc17(-); Hatch2N5(+)                     2
    ##  7 FH    ALDH1A2 Inc9toInc17(-); Hatch2N5(+)                     2
    ##  8 FH    ANKRD24 Inc9toInc17(+); Hatch2N5(-)                     2
    ##  9 FH    APBB1   Inc9toInc17(+); Hatch2N5(-)                     2
    ## 10 FH    ARHGEF3 Inc9toInc17(+); Hatch2N5(-)                     2
    ## # … with 1,914 more rows

    summarizedDEGs <- rbind(femaleDEGs,maleDEGs) %>%
      group_by(group, results) %>%
      summarize(genes = str_c(gene, collapse = "; ")) 
    summarizedDEGs 

    ## # A tibble: 109 x 3
    ## # Groups:   group [6]
    ##    group results                    genes                                       
    ##    <chr> <chr>                      <chr>                                       
    ##  1 FG    Bldg2Lay(-); In17toHatch(… LOC101749216                                
    ##  2 FG    Bldg2Lay(-); Inc3toInc9(+) PI15                                        
    ##  3 FG    Bldg2Lay(-); Inc9toInc17(… DMA                                         
    ##  4 FG    Bldg2Lay(-); Inc9toInc17(… LOC107049309; PRPS2                         
    ##  5 FG    Bldg2Lay(-); Lay2Inc3(-)   IL13RA2; SLC31A1                            
    ##  6 FG    Bldg2Lay(-); Lay2Inc3(-);… CA8; NRG2                                   
    ##  7 FG    Bldg2Lay(-); Lay2Inc3(+)   ADH5; ALDH7A1; ALKBH1; ANAPC15; ANOS1; AP00…
    ##  8 FG    Bldg2Lay(-); Lay2Inc3(+);… ANP32A; ANXA5; RGS7                         
    ##  9 FG    Bldg2Lay(-); Lay2Inc3(+);… FDX1L                                       
    ## 10 FG    Bldg2Lay(-); N5toN9(-)     AQP8; CHGB; CHRNA3                          
    ## # … with 99 more rows

    write.csv(summarizedDEGs,"../results/summarizedDEGsv2.csv")

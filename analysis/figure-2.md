Candidate gene analysis
=======================

    library(tidyverse)

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

    library(knitr)
    library(kableExtra)

    ## 
    ## Attaching package: 'kableExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

    source("../R/themes.R")

    knitr::opts_chunk$set(echo = TRUE, message = F, fig.path = "../figures/")

Candidate Genes (from literature and relevant GO terms)
-------------------------------------------------------

    # all genes in this parental care dataset 

    geneids <- read_csv("../metadata/00_geneinfo.csv") %>%
      select(-X1) %>%
      rename("NCBI" = "entrezid",
             "gene" = "Name")

    ## Warning: Missing column names filled in: 'X1' [1]

    # From Ch 17 Evolution of parental care
    curleychampagnegenes <- c("AVPR1A", "ADRA2A",    
                             "COMT", "CRH", "CRHBP", "CRHR1", "CRHR2", 
                              "DRD1", "DRD4","ESR1", "ESR2", #ERa ERbeta
                              "FOS", "HTR2C", "MEST", "NR3C1", #GR
                              "OPRM1", "PGR", "PRL", "PRLR",  "SLC6A4") #5HTT
    curleychampagnegenes <- as.data.frame(curleychampagnegenes) %>%
      mutate(GO = "parentalcare", gene = curleychampagnegenes)  %>%
      select(GO, gene)
    curleychampagnegenes

    ##              GO   gene
    ## 1  parentalcare AVPR1A
    ## 2  parentalcare ADRA2A
    ## 3  parentalcare   COMT
    ## 4  parentalcare    CRH
    ## 5  parentalcare  CRHBP
    ## 6  parentalcare  CRHR1
    ## 7  parentalcare  CRHR2
    ## 8  parentalcare   DRD1
    ## 9  parentalcare   DRD4
    ## 10 parentalcare   ESR1
    ## 11 parentalcare   ESR2
    ## 12 parentalcare    FOS
    ## 13 parentalcare  HTR2C
    ## 14 parentalcare   MEST
    ## 15 parentalcare  NR3C1
    ## 16 parentalcare  OPRM1
    ## 17 parentalcare    PGR
    ## 18 parentalcare    PRL
    ## 19 parentalcare   PRLR
    ## 20 parentalcare SLC6A4

    # Candidate GO terms 
    GO_path <- "../metadata/goterms/"   # path to the data
    GO_files <- dir(GO_path, pattern = "*.txt") # get file names
    GO_pathfiles <- paste0(GO_path, GO_files)

    GOgenesLong <- GO_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_table(.x, col_types = cols(), col_names = FALSE), .id = "file_name") %>% 
      mutate(GO = sapply(strsplit(as.character(file_name),'../metadata/goterms/'), "[", 2)) %>% 
      mutate(GO = sapply(strsplit(as.character(GO),'.txt'), "[", 1)) %>% 
      mutate(gene = sapply(strsplit(as.character(X1), "[\\\\]|[^[:print:]]" ), "[", 2)) %>% 
      select(GO, gene)  %>%
      filter(gene != "Symbol") %>%
      distinct(GO,gene)  %>%
      mutate(gene = toupper(gene)) %>%
      rbind(., curleychampagnegenes) %>%
      arrange(gene) %>%
      left_join(., geneids, by = "gene") %>%
      drop_na() 
    GOgenesLong

    ## # A tibble: 37 x 4
    ##    GO               gene   geneid NCBI          
    ##    <chr>            <chr>   <dbl> <chr>         
    ##  1 parentalcare     ADRA2A 428980 XP_004942333.2
    ##  2 parentalbehavior AVP    396101 NP_990516.1   
    ##  3 parentalbehavior AVPR1A 771773 NP_001103908.1
    ##  4 parentalcare     AVPR1A 771773 NP_001103908.1
    ##  5 parentalbehavior BRINP1 395098 NP_989780.1   
    ##  6 parentalcare     COMT   416783 XP_001233014.1
    ##  7 parentalbehavior CREBRF 416206 XP_001231574.1
    ##  8 parentalcare     CRH    404297 NP_001116503.1
    ##  9 parentalcare     CRHBP  427214 XP_003643006.2
    ## 10 parentalcare     CRHR1  374218 NP_989652.1   
    ## # … with 27 more rows

    GOgenesWide <- GOgenesLong %>% 
      pivot_wider(
        names_from = GO,
        values_from = gene) %>%
      mutate(numGOs = 4 - rowSums(is.na(.))) %>%
      arrange(desc(numGOs)) %>%
      select(-numGOs) %>%
      left_join(geneids, by = c("geneid","NCBI"))
    GOgenesWide

    ## # A tibble: 33 x 5
    ##    geneid NCBI           parentalcare parentalbehavior gene  
    ##     <dbl> <chr>          <chr>        <chr>            <chr> 
    ##  1 771773 NP_001103908.1 AVPR1A       AVPR1A           AVPR1A
    ##  2 427633 NP_001138320.1 DRD1         DRD1             DRD1  
    ##  3 416343 XP_015149519.1 NR3C1        NR3C1            NR3C1 
    ##  4 396453 NP_990797.2    PRL          PRL              PRL   
    ##  5 428980 XP_004942333.2 ADRA2A       <NA>             ADRA2A
    ##  6 396101 NP_990516.1    <NA>         AVP              AVP   
    ##  7 395098 NP_989780.1    <NA>         BRINP1           BRINP1
    ##  8 416783 XP_001233014.1 COMT         <NA>             COMT  
    ##  9 416206 XP_001231574.1 <NA>         CREBRF           CREBRF
    ## 10 404297 NP_001116503.1 CRH          <NA>             CRH   
    ## # … with 23 more rows

    # genes I can't find in dataset 
    # NOS1 or OXTR or FOX1B or HTR5A

Candidate DEGs
--------------

    candidategenes <- GOgenesLong %>% distinct(gene) %>% pull(gene)

    # summary DEG results from DESeq2
    candidateDEGS <- read_csv("../results/suppletable1.csv") %>%
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

    ## # A tibble: 32 x 7
    ## # Groups:   gene [32]
    ##    gene   bldg_lay lay_inc.d3 inc.d3_inc.d9 inc.d9_inc.d17 hatch_n5 n5_n9
    ##    <chr>  <chr>    <chr>      <chr>         <chr>          <chr>    <chr>
    ##  1 ADRA2A <NA>     <NA>       <NA>          <NA>           <NA>     <NA> 
    ##  2 AVP    <NA>     <NA>       <NA>          <NA>           <NA>     <NA> 
    ##  3 AVPR1A <NA>     FG+        FG-           <NA>           <NA>     <NA> 
    ##  4 BRINP1 <NA>     <NA>       <NA>          <NA>           <NA>     <NA> 
    ##  5 COMT   <NA>     <NA>       <NA>          <NA>           <NA>     <NA> 
    ##  6 CREBRF FG+      FG-        <NA>          <NA>           <NA>     <NA> 
    ##  7 CRH    <NA>     <NA>       <NA>          <NA>           <NA>     <NA> 
    ##  8 CRHBP  <NA>     <NA>       <NA>          <NA>           FH+      <NA> 
    ##  9 CRHR1  <NA>     <NA>       <NA>          <NA>           <NA>     <NA> 
    ## 10 CRHR2  <NA>     <NA>       <NA>          <NA>           FH+      <NA> 
    ## # … with 22 more rows

    table1 <- left_join(candidateDEGS, GOgenesWide) %>%
      select(gene, bldg_lay:n5_n9, parentalcare, parentalbehavior, NCBI) %>%
      mutate(parentalcare = if_else(is.na(parentalcare), " ", "X"),
             parentalbehavior = if_else(is.na(parentalbehavior), " ", "X")) %>%
      rename("Literature" = "parentalcare", "GO" =  "parentalbehavior")
    table1$numDEGs <- rowSums(is.na(table1)) # count NAs to know how many are NS
    table1 <- table1 %>% 
      mutate(sig = ifelse(numDEGs == 6, "NS", "DEG")) %>% 
      arrange(sig, gene)  %>%  select(-sig, -numDEGs)
    table1[is.na(table1)] <- " " # replace NA with blank space so it's pretty
    #table1

    kable(table1)

<table>
<thead>
<tr>
<th style="text-align:left;">
gene
</th>
<th style="text-align:left;">
bldg\_lay
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
<th style="text-align:left;">
n5\_n9
</th>
<th style="text-align:left;">
Literature
</th>
<th style="text-align:left;">
GO
</th>
<th style="text-align:left;">
NCBI
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
AVPR1A
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FG+
</td>
<td style="text-align:left;">
FG-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
NP\_001103908.1
</td>
</tr>
<tr>
<td style="text-align:left;">
CREBRF
</td>
<td style="text-align:left;">
FG+
</td>
<td style="text-align:left;">
FG-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
XP\_001231574.1
</td>
</tr>
<tr>
<td style="text-align:left;">
CRHBP
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FH+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
XP\_003643006.2
</td>
</tr>
<tr>
<td style="text-align:left;">
CRHR2
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FH+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
NP\_989785.1
</td>
</tr>
<tr>
<td style="text-align:left;">
DRD1
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FH+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
NP\_001138320.1
</td>
</tr>
<tr>
<td style="text-align:left;">
ESR1
</td>
<td style="text-align:left;">
FP+
</td>
<td style="text-align:left;">
FP-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
XP\_015139536.1
</td>
</tr>
<tr>
<td style="text-align:left;">
FOS
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FG+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
NP\_990839.1
</td>
</tr>
<tr>
<td style="text-align:left;">
GNAQ
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FG-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FH+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
NP\_001026598.1
</td>
</tr>
<tr>
<td style="text-align:left;">
HTR2C
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FH+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
XP\_004940707.1
</td>
</tr>
<tr>
<td style="text-align:left;">
MEST
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FP+
</td>
<td style="text-align:left;">
FP-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
XP\_015142671.1
</td>
</tr>
<tr>
<td style="text-align:left;">
NR3C1
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FG-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
XP\_015149519.1
</td>
</tr>
<tr>
<td style="text-align:left;">
OPRK1
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FH+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
XP\_426087.2
</td>
</tr>
<tr>
<td style="text-align:left;">
OPRM1
</td>
<td style="text-align:left;">
FG-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FG+
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
XP\_003641008.2
</td>
</tr>
<tr>
<td style="text-align:left;">
PGR
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FH+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
NP\_990593.1
</td>
</tr>
<tr>
<td style="text-align:left;">
PRL
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FP+ MP+
</td>
<td style="text-align:left;">
FP-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
NP\_990797.2
</td>
</tr>
<tr>
<td style="text-align:left;">
PRLR
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FG-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
XP\_015132722.1
</td>
</tr>
<tr>
<td style="text-align:left;">
PTEN
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FP-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
XP\_015134187.1
</td>
</tr>
<tr>
<td style="text-align:left;">
ADRA2A
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
XP\_004942333.2
</td>
</tr>
<tr>
<td style="text-align:left;">
AVP
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
NP\_990516.1
</td>
</tr>
<tr>
<td style="text-align:left;">
BRINP1
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
NP\_989780.1
</td>
</tr>
<tr>
<td style="text-align:left;">
COMT
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
XP\_001233014.1
</td>
</tr>
<tr>
<td style="text-align:left;">
CRH
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
NP\_001116503.1
</td>
</tr>
<tr>
<td style="text-align:left;">
CRHR1
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
NP\_989652.1
</td>
</tr>
<tr>
<td style="text-align:left;">
DBH
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
XP\_415429.5
</td>
</tr>
<tr>
<td style="text-align:left;">
DRD4
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
NP\_001136321.1
</td>
</tr>
<tr>
<td style="text-align:left;">
ESR2
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
NP\_990125.1
</td>
</tr>
<tr>
<td style="text-align:left;">
KALRN
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
XP\_015145468.1
</td>
</tr>
<tr>
<td style="text-align:left;">
MBD2
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
NP\_001012403.1
</td>
</tr>
<tr>
<td style="text-align:left;">
NPAS3
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
XP\_015143131.1
</td>
</tr>
<tr>
<td style="text-align:left;">
NPAS3
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
XP\_015143132.1
</td>
</tr>
<tr>
<td style="text-align:left;">
OXT
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
XP\_004936337.1
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC6A4
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
XP\_015151186.1
</td>
</tr>
<tr>
<td style="text-align:left;">
ZFX
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
X
</td>
<td style="text-align:left;">
XP\_015127980.1
</td>
</tr>
</tbody>
</table>

variance stabilized gene expression (vsd)
-----------------------------------------

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

    candidategenes <- GOgenesLong %>% distinct(gene) %>% pull(gene)

    hypvsd <- getcandidatevsd(candidategenes, "hypothalamus", sexlevels)
    pitvsd <- getcandidatevsd(candidategenes, "pituitary", sexlevels)
    gonvsd <- getcandidatevsd(candidategenes, "gonad", sexlevels)
    candidatevsd <- rbind(hypvsd, pitvsd)
    candidatevsd <- rbind(candidatevsd, gonvsd)
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

    # significant candidate genes
    hypDEGs <- c("CRHBP", "CRHR2", "DRD1", "GNAQ", 
                 "HTR2C", "OPRK1", "PGR")
    pitDEGs <- c("ESR1", "MEST", "PRL", "PTEN")
    gonadDEGs <- c("AVPR1A", "CREBRF", "FOS", "GNAQ", 
                   "NR3C1", "OPRM1", "PRLR")

    candidateboxplot <- function(whichtissue, whichgenes, whichsex){
      
      p <- candidatevsd %>%
        filter(tissue %in% whichtissue,
              gene %in% whichgenes,
              sex %in% whichsex) %>%
        mutate(treatment = factor(treatment, levels = charlevels)) %>%
        mutate(treatmentNum = as.numeric(treatment)) %>%
        ggplot(aes(x = treatmentNum, y = counts)) +
        geom_boxplot(aes(fill = treatment, color = sex), outlier.shape = NA) +
        #geom_smooth(aes(color = sex)) +
        geom_jitter(size = 0.25, aes(color = sex)) +
        facet_wrap(~gene, scales = "free_y", nrow = 1) +
        scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9),
                           labels = charlevels) +
        scale_fill_manual(values = allcolors) +
        scale_color_manual(values = allcolors) +
        theme_B3() + 
        theme(legend.position = "none",
              strip.text = element_text(face = "italic"),
              axis.text.x = element_blank()) +
        labs(subtitle = whichsex,
             y = "gene expression",
             x = NULL) 
      return(p)
    }

    a1 <- candidateboxplot("hypothalamus", hypDEGs, "female") + labs(x = NULL, title = "hypothalamus") 
    a2 <- candidateboxplot("hypothalamus", hypDEGs, "male") + labs(x = NULL) +
      theme(strip.text = element_blank())

    b1 <- candidateboxplot("pituitary", pitDEGs, "female") + labs(x = NULL, title = "pituitary")
    b2 <- candidateboxplot("pituitary", pitDEGs, "male") + labs(x = NULL) +
      theme(strip.text = element_blank())

    c1 <- candidateboxplot("gonad", gonadDEGs, "female") + labs(x = NULL, title = "gonad")
    c2 <- candidateboxplot("gonad", gonadDEGs, "male") + 
      theme(strip.text = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))


    plot_grid(a1,a2, b1,b2, c1,c2, nrow = 6, labels = c("A", " ","B", " ", "C", " "), label_size = 8,
              rel_heights = c(1, 0.75, 1, 0.75, 1, 1))

![](../figures/fig2-1.png)

    #write.csv(candidatevsd, "../../musicalgenes/data/candidatecounts.csv")
    #write.csv(candidatevsd, "../results/candidatecounts.csv")
    write.csv(table1, "../results/table1.csv")

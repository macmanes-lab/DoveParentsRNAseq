Candidate gene analysis
=======================

    library(tidyverse)
    library(ggtext)
    library(cowplot)
    library(ggpubr)
    library(knitr)
    library(kableExtra)
    library(corrr)

    source("../R/themes.R")
    source("../R/functions.R")
    source("../R/wrangledata.R")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    knitr::opts_chunk$set(echo = TRUE, message = F, fig.path = "../figures/")

Candidate DEGs
--------------

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

    ## # A tibble: 114 x 7
    ## # Groups:   gene [114]
    ##    gene   bldg_lay lay_inc.d3 inc.d3_inc.d9 inc.d9_inc.d17 hatch_n5 n5_n9
    ##    <chr>  <chr>    <chr>      <chr>         <chr>          <chr>    <chr>
    ##  1 ADRA2A <NA>     <NA>       <NA>          <NA>           <NA>     <NA> 
    ##  2 AHR    <NA>     <NA>       <NA>          <NA>           <NA>     <NA> 
    ##  3 ANXA1  <NA>     <NA>       <NA>          FG-            <NA>     <NA> 
    ##  4 AR     <NA>     <NA>       <NA>          MH+            <NA>     <NA> 
    ##  5 AREG   <NA>     <NA>       <NA>          <NA>           <NA>     FP+  
    ##  6 ATP2C2 <NA>     <NA>       <NA>          <NA>           <NA>     FG+  
    ##  7 AVP    <NA>     <NA>       <NA>          <NA>           <NA>     <NA> 
    ##  8 AVPR1A <NA>     FG+        FG-           <NA>           <NA>     <NA> 
    ##  9 BMP4   <NA>     <NA>       <NA>          <NA>           <NA>     <NA> 
    ## 10 BRCA2  <NA>     <NA>       <NA>          <NA>           <NA>     <NA> 
    ## # … with 104 more rows

    ## table 1 summary candidate genes
    table1 <- left_join(candidateDEGS, parentalcaregenes) %>%
      select(gene, bldg_lay:n5_n9, parentalcare, parentalbehavior, NCBI) %>%
      mutate(parentalcare = if_else(is.na(parentalcare), " ", "X"),
             parentalbehavior = if_else(is.na(parentalbehavior), " ", "X")) %>%
      rename("Literature" = "parentalcare", "GO" =  "parentalbehavior")
    table1$numDEGs <- rowSums(is.na(table1)) # count NAs to know how many are NS
    table1 <- table1 %>% 
      mutate(sig = ifelse(numDEGs == 6, "NS", "DEG")) %>% 
      arrange(sig, gene)  %>%  select(-sig, -numDEGs)
    table1[is.na(table1)] <- " " # replace NA with blank space so it's pretty
    table1

    ## # A tibble: 115 x 10
    ## # Groups:   gene [114]
    ##    gene  bldg_lay lay_inc.d3 inc.d3_inc.d9 inc.d9_inc.d17 hatch_n5 n5_n9
    ##    <chr> <chr>    <chr>      <chr>         <chr>          <chr>    <chr>
    ##  1 AHR   " "      " "        " "           " "            " "      " "  
    ##  2 AVPR… " "      "FG+"      "FG-"         " "            " "      " "  
    ##  3 BMP4  " "      " "        " "           " "            " "      " "  
    ##  4 BRCA2 " "      " "        " "           " "            " "      " "  
    ##  5 BTRC  " "      " "        " "           " "            " "      " "  
    ##  6 CCND1 " "      " "        " "           " "            " "      " "  
    ##  7 CEBPB " "      " "        " "           " "            " "      " "  
    ##  8 CREB… "FG+"    "FG-"      " "           " "            " "      " "  
    ##  9 CRHBP " "      " "        " "           " "            "FH+"    " "  
    ## 10 CRHR2 " "      " "        " "           " "            "FH+"    " "  
    ## # … with 105 more rows, and 3 more variables: Literature <chr>, GO <chr>,
    ## #   NCBI <chr>

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
AHR
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
</td>
<td style="text-align:left;">
</td>
</tr>
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
BMP4
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
BRCA2
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
BTRC
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
CCND1
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
CEBPB
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
</td>
<td style="text-align:left;">
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
CTNNB1
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
CYP7B1
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
</td>
<td style="text-align:left;">
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
EAF2
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
ERBB4
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
</td>
<td style="text-align:left;">
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
ETV5
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
FGF10
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
FGF2
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
FKBP4
</td>
<td style="text-align:left;">
FG-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FP+
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
FOXA1
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
</td>
<td style="text-align:left;">
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
HOXA5
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
</td>
<td style="text-align:left;">
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
ID2
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
IGF1
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
FH+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
JAK2
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
LBH
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
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
LRP6
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
</td>
<td style="text-align:left;">
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
MST1
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
MSX2
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
NCOA3
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
</td>
<td style="text-align:left;">
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
ORAI1
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FG-
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
PHB2
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
PLAG1
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
PML
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
</td>
<td style="text-align:left;">
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
PSAP
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
FH-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
PTHLH
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
ROBO1
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
SCRIB
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
SERPINB5
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
FH+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
SERPINF1
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
SFRP1
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
SHH
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
SLC12A2
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
FH-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
SOSTDC1
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
FH+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
SOX9
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
MG-
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
</td>
</tr>
<tr>
<td style="text-align:left;">
STK11
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
SULF1
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
TBX3
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
TFAP2C
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
TGFB1
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
TNC
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
TNFSF11
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
UBE3A
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
WDR77
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
FP- MP-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
WNT5A
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
</td>
<td style="text-align:left;">
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
ANXA1
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
</td>
</tr>
<tr>
<td style="text-align:left;">
AR
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
MH+
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
</tr>
<tr>
<td style="text-align:left;">
AREG
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
FP+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
ATP2C2
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
FG+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
X
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
CD44
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
</tr>
<tr>
<td style="text-align:left;">
CDKN1B
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
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
CSF1
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
FH-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
CYP19A1
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
</td>
<td style="text-align:left;">
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
EPHA2
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
FEM1B
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
MP-
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
</tr>
<tr>
<td style="text-align:left;">
FGFR2
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
FH-
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
FRS2
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
</td>
</tr>
<tr>
<td style="text-align:left;">
GLI2
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
</tr>
<tr>
<td style="text-align:left;">
HIF1A
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FG- MG-
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
</tr>
<tr>
<td style="text-align:left;">
ID4
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FP- FG-
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
</tr>
<tr>
<td style="text-align:left;">
IGF1R
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
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
LATS1
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
LRP5
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
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
MAPK1
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
</td>
<td style="text-align:left;">
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
MED1
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
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
MMP2
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
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
NOG
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
NTN1
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
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
X
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
PTCH1
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
FP+
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
PYGO2
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
FP- MP-
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
</tr>
<tr>
<td style="text-align:left;">
RLN1
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
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
RREB1
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
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
RTN4
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
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
RXRA
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
SLIT2
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
</td>
<td style="text-align:left;">
</td>
</tr>
<tr>
<td style="text-align:left;">
SMO
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
</tr>
<tr>
<td style="text-align:left;">
STAT6
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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
</tr>
<tr>
<td style="text-align:left;">
WNT4
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
</td>
<td style="text-align:left;">
</td>
<td style="text-align:left;">
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

Candidate Correlations
----------------------

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

    hyp1 <- makecorrdf("female", "hypothalamus", curleychampagnegenes)  

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

    pit1 <- makecorrdf("female", "pituitary", curleychampagnegenes)  

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

    gon1 <- makecorrdf("female", "gonad", curleychampagnegenes)  

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

    hyp2 <- makecorrdf("male", "hypothalamus", curleychampagnegenes)  

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

    pit2 <- makecorrdf("male", "pituitary", curleychampagnegenes)  

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

    gon2 <- makecorrdf("male", "gonad", curleychampagnegenes)  

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

    hyp3 <- makecorrdf(sexlevels, "hypothalamus", curleychampagnegenes)  

    ## Warning in sex == whichsex: longer object length is not a multiple of
    ## shorter object length

    ## # A tibble: 6 x 23
    ##   rowname  HTR2C   DRD1  CRHR2  CRHBP  CRHR1   DRD4   CRH  MEST   ESR1
    ##   <chr>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl>  <dbl>
    ## 1 HTR2C   NA      0.864  0.851  0.731  0.526  0.370 0.460 0.498 0.385 
    ## 2 DRD1     0.864 NA      0.775  0.683  0.472  0.533 0.349 0.238 0.0413
    ## 3 CRHR2    0.851  0.775 NA      0.658  0.357  0.407 0.171 0.548 0.217 
    ## 4 CRHBP    0.731  0.683  0.658 NA      0.403  0.166 0.350 0.370 0.228 
    ## 5 CRHR1    0.526  0.472  0.357  0.403 NA      0.309 0.717 0.259 0.202 
    ## 6 DRD4     0.370  0.533  0.407  0.166  0.309 NA     0.139 0.542 0.0486
    ## # … with 13 more variables: AVPR1A <dbl>, FOS <dbl>, OPRM1 <dbl>,
    ## #   SLC6A4 <dbl>, NR3C1 <dbl>, PRLR <dbl>, ADRA2A <dbl>, PGR <dbl>,
    ## #   ESR2 <dbl>, AVP <dbl>, OXT <dbl>, PRL <dbl>, COMT <dbl>

    pit3 <- makecorrdf(sexlevels, "pituitary", curleychampagnegenes)  

    ## Warning in sex == whichsex: longer object length is not a multiple of
    ## shorter object length

    ## # A tibble: 6 x 22
    ##   rowname    AVP    PGR    OXT  CRHBP   PRLR   DRD4    COMT   CRHR2    MEST
    ##   <chr>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 AVP     NA      0.204  0.872  0.730  0.550  0.395  0.234  0.402    0.150 
    ## 2 PGR      0.204 NA      0.277  0.269  0.214  0.328 -0.0943 0.00797 -0.135 
    ## 3 OXT      0.872  0.277 NA      0.701  0.488  0.355  0.284  0.457    0.0521
    ## 4 CRHBP    0.730  0.269  0.701 NA      0.413  0.305  0.279  0.444   -0.0792
    ## 5 PRLR     0.550  0.214  0.488  0.413 NA      0.566  0.219  0.251    0.137 
    ## 6 DRD4     0.395  0.328  0.355  0.305  0.566 NA     -0.0451 0.0684  -0.206 
    ## # … with 12 more variables: FOS <dbl>, ADRA2A <dbl>, OPRM1 <dbl>,
    ## #   ESR1 <dbl>, PRL <dbl>, ESR2 <dbl>, HTR2C <dbl>, NR3C1 <dbl>,
    ## #   CRHR1 <dbl>, SLC6A4 <dbl>, DRD1 <dbl>, AVPR1A <dbl>

    gon3 <- makecorrdf(sexlevels, "gonad", curleychampagnegenes)  

    ## Warning in sex == whichsex: longer object length is not a multiple of
    ## shorter object length

    ## # A tibble: 6 x 23
    ##   rowname    FOS ADRA2A   DRD4  HTR2C   ESR2  CRHR1 CRHBP  ESR1   PGR NR3C1
    ##   <chr>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl>
    ## 1 FOS     NA      0.931  0.885  0.886  0.893  0.873 0.888 0.836 0.777 0.725
    ## 2 ADRA2A   0.931 NA      0.875  0.847  0.851  0.828 0.832 0.912 0.860 0.703
    ## 3 DRD4     0.885  0.875 NA      0.849  0.835  0.865 0.829 0.819 0.767 0.689
    ## 4 HTR2C    0.886  0.847  0.849 NA      0.905  0.919 0.919 0.723 0.669 0.640
    ## 5 ESR2     0.893  0.851  0.835  0.905 NA      0.907 0.958 0.701 0.628 0.624
    ## 6 CRHR1    0.873  0.828  0.865  0.919  0.907 NA     0.923 0.677 0.639 0.629
    ## # … with 12 more variables: AVPR1A <dbl>, COMT <dbl>, MEST <dbl>,
    ## #   SLC6A4 <dbl>, DRD1 <dbl>, AVP <dbl>, OXT <dbl>, CRH <dbl>, PRL <dbl>,
    ## #   CRHR2 <dbl>, PRLR <dbl>, OPRM1 <dbl>

    candidatevsdwide <- candidatevsd  %>%
        pivot_wider(names_from = gene, values_from = counts) 
    FH <- subsetcandidatevsdwide("female", "hypothalamus")
    FP <- subsetcandidatevsdwide("female", "pituitary")
    FG <- subsetcandidatevsdwide("female", "gonad")
    MH <- subsetcandidatevsdwide("male", "hypothalamus")
    MP <- subsetcandidatevsdwide("male", "pituitary")
    MG <- subsetcandidatevsdwide("male", "gonad") 

### Candidate gene functions

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
        facet_wrap(~sex, scales = "free_y", nrow = 2) +
        scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9),
                           labels = charlevels) +
        scale_fill_manual(values = allcolors) +
        scale_color_manual(values = allcolors) +
        theme_B3() + 
        theme(legend.position = "none",
              axis.title.y = element_text(face = "italic"),
              axis.text.x = element_blank()) +
        labs(y = whichgenes,
             x = NULL) 
      return(p)
    }


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

    ## plot sig candidate genes
    a <- candidateboxplot("hypothalamus", "DRD1", sexlevels) + labs(x = NULL, title = "Hypothalamus") 

    g <- candidateboxplot("pituitary", "PRL", sexlevels) + labs(x = NULL, title = "Pituitary")

    m <- candidateboxplot("gonad", "CREBRF", sexlevels) + labs(x = NULL, title = "Gonad") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))


    c <- plotcorrplot(hyp3, "females and males") + theme(legend.position = "none")  + labs(title = "")
    i <- plotcorrplot(pit3, "females and males") + theme(legend.position = "none") + labs(title = "")
    o <- plotcorrplot(gon3, "females and males") + theme(legend.position = "bottom")  + labs(title = "")


    e <- scattercorrelations(FH, FH$DRD1, FH$HTR2C, "#969696", "DRD1", "HTR2C") + labs(subtitle = "female")+ labs(title = "")
    k <- scattercorrelations(FP, FP$AVPR1A, FP$CRHR1, "#969696", "AVPR1A", "CRHR1")  + labs(subtitle = "female")+ labs(title = "")
    q <- scattercorrelations(FG, FG$ESR2, FG$CRHBP, "#969696", "ESR2", "CRHBP")   + labs(subtitle = "female")+ labs(title = "")

    f <- scattercorrelations(MH, MH$DRD1, MH$HTR2C, "#525252", "DRD1", "HTR2C") + theme(axis.title.y = element_blank()) + labs(subtitle = "male") + labs(title = "")
    l <- scattercorrelations(MP, MP$AVPR1A, MP$CRHR1, "#525252", "AVPR1A", "CRHR1") + theme(axis.title.y = element_blank()) + labs(subtitle = "male") + labs(title = "")
    r <- scattercorrelations(MG, MG$ESR2, MG$CRHBP, "#525252", "ESR2", "CRHBP")  + theme(axis.title.y = element_blank()) + labs(subtitle = "male") + labs(title = "")

    plot_grid(a,c,e,f,g,i,k,l,m,o,q,r, nrow = 3, rel_widths = c(1,1.2,0.8,0.75), rel_heights = c(1,1,1.2),
              labels = "AUTO", label_size =  8)

![](../figures/fig2-1.png)

    #write.csv(candidatevsd, "../../musicalgenes/data/candidatecounts.csv")
    #write.csv(candidatevsd, "../results/candidatecounts.csv")
    write.csv(table1, "../results/table1.csv")
    write.csv(candidategenes, "../results/candidategenes.csv")

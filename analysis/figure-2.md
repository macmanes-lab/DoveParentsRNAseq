Candidate gene analysis
=======================

    library(tidyverse)

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

    ## # A tibble: 17 x 7
    ## # Groups:   gene [17]
    ##    gene   bldg_lay lay_inc.d3 inc.d3_inc.d9 inc.d9_inc.d17 hatch_n5 n5_n9
    ##    <chr>  <chr>    <chr>      <chr>         <chr>          <chr>    <chr>
    ##  1 AVPR1A <NA>     FG+        FG-           <NA>           <NA>     <NA> 
    ##  2 CREBRF FG+      FG-        <NA>          <NA>           <NA>     <NA> 
    ##  3 CRHBP  <NA>     <NA>       <NA>          <NA>           FH+      <NA> 
    ##  4 CRHR2  <NA>     <NA>       <NA>          <NA>           FH+      <NA> 
    ##  5 DRD1   <NA>     <NA>       <NA>          <NA>           FH+      <NA> 
    ##  6 ESR1   FP+      FP-        <NA>          <NA>           <NA>     <NA> 
    ##  7 FOS    <NA>     FG+        <NA>          <NA>           <NA>     <NA> 
    ##  8 GNAQ   <NA>     FG-        <NA>          <NA>           FH+      <NA> 
    ##  9 HTR2C  <NA>     <NA>       <NA>          <NA>           FH+      <NA> 
    ## 10 MEST   <NA>     <NA>       <NA>          FP+            FP-      <NA> 
    ## 11 NR3C1  <NA>     FG-        <NA>          <NA>           <NA>     <NA> 
    ## 12 OPRK1  <NA>     <NA>       <NA>          <NA>           FH+      <NA> 
    ## 13 OPRM1  FG-      <NA>       <NA>          <NA>           <NA>     FG+  
    ## 14 PGR    <NA>     <NA>       <NA>          <NA>           FH+      <NA> 
    ## 15 PRL    <NA>     <NA>       <NA>          FP+ MP+        FP-      <NA> 
    ## 16 PRLR   <NA>     FG-        <NA>          <NA>           <NA>     <NA> 
    ## 17 PTEN   <NA>     FP-        <NA>          <NA>           <NA>     <NA>

    table1 <- left_join(candidateDEGS, GOgenesWide) %>%
      select(gene, bldg_lay:n5_n9, parentalcare, parentalbehavior, NCBI, ) %>%
      mutate(parentalcare = if_else(is.na(parentalcare), " ", "X"),
             parentalbehavior = if_else(is.na(parentalbehavior), " ", "X")) %>%
      rename("Literature" = "parentalcare", "GO" =  "parentalbehavior")
    table1[is.na(table1)] <- " "
    table1

    ## # A tibble: 17 x 10
    ## # Groups:   gene [17]
    ##    gene  bldg_lay lay_inc.d3 inc.d3_inc.d9 inc.d9_inc.d17 hatch_n5 n5_n9
    ##    <chr> <chr>    <chr>      <chr>         <chr>          <chr>    <chr>
    ##  1 AVPR… " "      "FG+"      "FG-"         " "            " "      " "  
    ##  2 CREB… "FG+"    "FG-"      " "           " "            " "      " "  
    ##  3 CRHBP " "      " "        " "           " "            "FH+"    " "  
    ##  4 CRHR2 " "      " "        " "           " "            "FH+"    " "  
    ##  5 DRD1  " "      " "        " "           " "            "FH+"    " "  
    ##  6 ESR1  "FP+"    "FP-"      " "           " "            " "      " "  
    ##  7 FOS   " "      "FG+"      " "           " "            " "      " "  
    ##  8 GNAQ  " "      "FG-"      " "           " "            "FH+"    " "  
    ##  9 HTR2C " "      " "        " "           " "            "FH+"    " "  
    ## 10 MEST  " "      " "        " "           "FP+"          "FP-"    " "  
    ## 11 NR3C1 " "      "FG-"      " "           " "            " "      " "  
    ## 12 OPRK1 " "      " "        " "           " "            "FH+"    " "  
    ## 13 OPRM1 "FG-"    " "        " "           " "            " "      "FG+"
    ## 14 PGR   " "      " "        " "           " "            "FH+"    " "  
    ## 15 PRL   " "      " "        " "           "FP+ MP+"      "FP-"    " "  
    ## 16 PRLR  " "      "FG-"      " "           " "            " "      " "  
    ## 17 PTEN  " "      "FP-"      " "           " "            " "      " "  
    ## # … with 3 more variables: Literature <chr>, GO <chr>, NCBI <chr>

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
</tbody>
</table>

variance stabilized gene expression (vsd)
-----------------------------------------

    vsd_path <- "../results/DEseq2/"   # path to the data
    vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
    vsd_pathfiles <- paste0(vsd_path, vsd_files)
    vsd_files

    allvsd <- vsd_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name")  %>% 
      dplyr::rename("gene" = "X1") %>% 
      pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
                   names_to = "samples", values_to = "counts") 

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

    head(allvsd)
    candidategenes <- GOgenesLong %>% distinct(gene) %>% pull(gene)

    hypvsd <- getcandidatevsd(candidategenes, "hypothalamus", sexlevels)
    pitvsd <- getcandidatevsd(candidategenes, "pituitary", sexlevels)
    gonvsd <- getcandidatevsd(candidategenes, "gonad", sexlevels)
    candidatevsd <- rbind(hypvsd, pitvsd)
    candidatevsd <- rbind(candidatevsd, gonvsd)
    head(candidatevsd)

    newboxplot <- function(df, mygenes, whichsex, whichstage){
      
      p <- df %>%
        filter(gene %in% mygenes, 
               sex %in% whichsex,
               treatment %in% whichstage) %>%
        droplevels() %>%
        mutate(tissue = factor(tissue, levels = tissuelevel)) %>%
        ggplot(aes(x = treatment, y = counts, fill = treatment, color = sex)) +
        geom_boxplot() +
        geom_jitter(size = 0.5, width = 0.1) +
        facet_grid(tissue~gene, scales = "free_y") +
        theme_B3() +
        scale_color_manual(values = allcolors) +
        scale_fill_manual(values = allcolors) +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text.x = element_text(face = "italic"),
              axis.title.x = element_blank())  +
        labs(x = "Parental stage",
             y = "variance stabilized expression")
      return(p)
    }

    a1 <- newboxplot(hypvsd, c("PGR",  "CYP19A1", "LHCGR" ), "female", c( "hatch", "n5")) + 
      theme(strip.text.y = element_blank()) + labs(subtitle = "Females")
    a2 <- newboxplot(hypvsd, c("AR"), "male", c( "hatch", "n5")) + 
      theme(axis.title.y = element_blank()) + labs(subtitle = "Males")

    b1 <- newboxplot(gonvsd, c("HSD3B2"), "female", c( "bldg", "lay"))  +
      theme(strip.text.y = element_blank(),
            axis.title.y = element_blank()) + labs(subtitle = "Females")
    b2 <- newboxplot(gonvsd, c("VIPR1",  "PRLR" ), "female", c( "lay", "inc.d3"))  + 
      theme(axis.title.y = element_blank()) + labs(subtitle = " ")

    ab <- plot_grid(a1,a2,b1,b2, rel_widths = c(6,2.5,2,4), labels = c("a", " ", "b", " "), label_size = 8, nrow = 1)

    c1 <- newboxplot(pitvsd, c( "ESR1", "GNRHR"), "female", c("bldg", "lay",  "inc.d3")) 
    c2 <- newboxplot(pitvsd, c( "PRL"), "female", c("inc.d9", "inc.d17", "hatch", "n5")) 
    c3 <- newboxplot(pitvsd, c("VIPR1",  "PRL" ), "male", c( "inc.d9", "inc.d17")) 


    c <- plot_grid(c1 + theme(strip.text.y = element_blank()) + labs(subtitle = "Females"),
              c2 + theme(strip.text.y = element_blank(),
                        axis.title.y = element_blank()) + labs(subtitle = " ") ,
              c3 + theme(axis.title.y = element_blank()) + labs(subtitle = "Males"), 
              nrow = 1, rel_widths = c(6,4,4),
              labels = c("c"), label_size = 8)

    plot_grid(ab,c, nrow = 2)

    # correlations
    library(corrr)


    hormones <- read_csv("../results/07_hormoneswide.csv") %>%
      select(-study, -treatment, -sex)
    head(hormones)

    head(candidatevsd)

    vsdhormones <- candidatevsd %>%
      pivot_wider(names_from = gene, values_from = counts) %>%
      mutate(bird_id = sapply(strsplit(samples,"\\_"), "[", 1)) %>%
      left_join(.,hormones, by =  "bird_id")
    vsdhormones


    vsdhormoneshyp <- vsdhormones %>% filter(tissue == "hypothalamus" )
    vsdhormonespit <- vsdhormones %>% filter(tissue == "pituitary" )
    vsdhormonesgon <- vsdhormones %>% filter(tissue == "gonad" )

    plotcorrs <- function(whichtissue, whichsex){
      
      x <- vsdhormones  %>%
        filter(tissue == whichtissue, 
               sex == whichsex) %>%
        select(-sex, - tissue, - treatment, -samples, -bird_id, -moltbin) %>%
        correlate()
      print(fashion(x))
      
      p <-rplot(x, colors = c("#0571b0", "#92c5de", "#f7f7f7","#f7f7f7","#f7f7f7", "#f4a582", "#ca0020")) +
        theme_B3() +
        theme(legend.position = "bottom",
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text = element_text(face = "italic")) +
        labs(subtitle = whichtissue)
      return(p)
    }

    a <- plotcorrs("hypothalamus", "female") + theme(legend.position = "none")
    b <- plotcorrs("pituitary", "female") + theme(axis.text.y = element_blank()) + theme(legend.position = "none")
    c <- plotcorrs("gonad", "female") + theme(axis.text.y = element_blank()) + theme(legend.position = "none")

    d <- plotcorrs("hypothalamus", "male") 
    e <- plotcorrs("pituitary", "male") + theme(axis.text.y = element_blank())
    f <- plotcorrs("gonad", "male") + theme(axis.text.y = element_blank())

    plot_grid(a,b,c, d,e,f, nrow = 2, rel_widths = c(1.2,1,1), rel_heights = c(1,1.2),
              labels = c("females", " ", " ",
                         "males", " ", ""),
              label_size = 8)

    plotspecificcorrs <- function(df, myx, myy, myxlab, myylab){
      p <- df %>%
        mutate(tissue = factor(tissue, levels = tissuelevel)) %>%
        ggplot(aes(x = myx, y = myy)) +
        geom_point(aes(color = treatment)) +
        geom_smooth(method = "lm", aes(color = sex)) +
        facet_wrap(~sex, scales = "free", nrow = 2) +
        scale_color_manual(values = allcolors) +
        labs(x = myxlab, y = myylab) +
        theme_B3() +
        theme(legend.position = "none",
              axis.title = element_text(face = "italic"))
      return(p)
    }


    a <- plotspecificcorrs(vsdhormoneshyp, vsdhormoneshyp$FSHB, vsdhormoneshyp$PRL, "FSHB", "PRL") + labs(subtitle = "Hypothalamus")
    b <- plotspecificcorrs(vsdhormonespit, vsdhormonespit$prolactin, vsdhormonespit$PRL, "circulating prolactin", "PRL") +
      theme(axis.title.x = element_text(face = "plain")) + labs(subtitle = "Pituitary")
    c <- plotspecificcorrs(vsdhormonesgon, vsdhormonesgon$ESR1, vsdhormonesgon$PGR, "ESR1", "PGR")  + labs(subtitle = "Gonads")

    plot_grid(a,b,c,nrow = 1, labels = "auto", label_size = 8)

    #write.csv(candidatevsd, "../../musicalgenes/data/candidatecounts.csv")
    #write.csv(candidatevsd, "../results/candidatecounts.csv")
    write.csv(table1, "../results/table1.csv")

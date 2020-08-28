    library(tidyverse)

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(readr)
    library(VennDiagram)

    ## Loading required package: grid

    ## Loading required package: futile.logger

    source("../R/themes.R")

    knitr::opts_chunk$set(echo = TRUE, message = F, fig.path = '../figures/')

### candidate gene anlayses

    # all genes in this parental care dataset 
    geneids <- read_csv("../metadata/00_geneinfo.csv") %>% select(-X1) 

    ## Warning: Missing column names filled in: 'X1' [1]

    # genes from Ch 17 Evolution of parental care
    curleychampagnegenes <- c("AVP", "AVPR1A", "ADRA2A",    
                              "COMT", "CRH", "CRHBP", "CRHR1", "CRHR2", 
                              "DRD1", "DRD4","ESR1", "ESR2", #ERa ERbeta
                              "FOS", "HTR2C", "MEST", "NR3C1", #GR
                              "OPRM1", "OXT" , "PGR", "PRL", "PRLR",  "SLC6A4") #5HTT

    datadrivengenes <- c("COX1", "COX2", "CYTB", "ATP2B4", "LOC107055658", 
                         "CHGA", "TF", "OVALY", "ANXA5")

    literaturegenes <- c(curleychampagnegenes, datadrivengenes,
                         "VIP", "GNIH", "GAL", "TH", "THRB",
                         "GNRHR", "GNRH1", "CGNRH-R", "FSHB", "FSHR")
    literaturegenes

    ##  [1] "AVP"          "AVPR1A"       "ADRA2A"       "COMT"         "CRH"         
    ##  [6] "CRHBP"        "CRHR1"        "CRHR2"        "DRD1"         "DRD4"        
    ## [11] "ESR1"         "ESR2"         "FOS"          "HTR2C"        "MEST"        
    ## [16] "NR3C1"        "OPRM1"        "OXT"          "PGR"          "PRL"         
    ## [21] "PRLR"         "SLC6A4"       "COX1"         "COX2"         "CYTB"        
    ## [26] "ATP2B4"       "LOC107055658" "CHGA"         "TF"           "OVALY"       
    ## [31] "ANXA5"        "VIP"          "GNIH"         "GAL"          "TH"          
    ## [36] "THRB"         "GNRHR"        "GNRH1"        "CGNRH-R"      "FSHB"        
    ## [41] "FSHR"

    literaturegenes <- as.data.frame(literaturegenes) %>%
      mutate(GO = "literature", gene = literaturegenes)  %>%
      select(GO, gene)

    # candidate gens from parental care GO terms 
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
      rbind(., literaturegenes) %>%
      arrange(gene) %>%
      left_join(., geneids, by = "gene") %>%
      arrange(gene) %>%
      drop_na() 
    head(GOgenesLong)

    ## # A tibble: 6 x 4
    ##   GO               gene   geneid NCBI          
    ##   <chr>            <chr>   <dbl> <chr>         
    ## 1 literature       ADRA2A 428980 XP_004942333.2
    ## 2 literature       ANXA5  428767 NP_001026709.1
    ## 3 literature       ATP2B4 419934 XP_015154447.1
    ## 4 literature       ATP2B4 419934 XP_015154449.1
    ## 5 parentalbehavior AVP    396101 NP_990516.1   
    ## 6 literature       AVP    396101 NP_990516.1

    parentalcaregenes <- GOgenesLong %>% 
      pivot_wider(
        names_from = GO,
        values_from = gene,
        values_fill = list(gene = "")) %>%
      left_join(geneids, by = c("geneid","NCBI")) %>%
      dplyr::rename("GO"= "parentalbehavior" )
    head(parentalcaregenes)

    ## # A tibble: 6 x 5
    ##   geneid NCBI           literature GO       gene  
    ##    <dbl> <chr>          <chr>      <chr>    <chr> 
    ## 1 428980 XP_004942333.2 ADRA2A     ""       ADRA2A
    ## 2 428767 NP_001026709.1 ANXA5      ""       ANXA5 
    ## 3 419934 XP_015154447.1 ATP2B4     ""       ATP2B4
    ## 4 419934 XP_015154449.1 ATP2B4     ""       ATP2B4
    ## 5 396101 NP_990516.1    AVP        "AVP"    AVP   
    ## 6 771773 NP_001103908.1 AVPR1A     "AVPR1A" AVPR1A

### variance stabilized gene expression (vsd)

    vsd_path <- "../results/DEseq2/treatment/"   # path to the data
    vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
    vsd_pathfiles <- paste0(vsd_path, vsd_files)
    vsd_files

    ## [1] "female_gonads_vsd.csv"       "female_hypothalamus_vsd.csv"
    ## [3] "female_pituitary_vsd.csv"    "male_gonads_vsd.csv"        
    ## [5] "male_hypothalamus_vsd.csv"   "male_pituitary_vsd.csv"

    ## before pivoting, check names of df with
    ## head(names(allvsd)) and tail(names(allvsd))

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

    allvsd %>% select(-file_name) %>% head()

    ## # A tibble: 6 x 3
    ##   gene  samples                            counts
    ##   <chr> <chr>                               <dbl>
    ## 1 A2ML1 L.G118_female_gonad_control         14.7 
    ## 2 A2ML1 R.G106_female_gonad_control          7.79
    ## 3 A2ML1 R.R20_female_gonad_control          14.9 
    ## 4 A2ML1 R.R9_female_gonad_control           11.7 
    ## 5 A2ML1 R.W44_female_gonad_control          13.2 
    ## 6 A2ML1 blk.s031.pu.d_female_gonad_prolong   7.42

    tail(names(allvsd))

    ## [1] "file_name" "gene"      "samples"   "counts"

### vsd for all and candidate genes

    candidategenes <- GOgenesLong %>% distinct(gene) %>% pull(gene)

    getcandidatevsd <- function(whichgenes, whichtissue, whichsex){
      candidates  <- allvsd %>%
        filter(gene %in% whichgenes) %>%
        dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
        dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/treatment/'), "[", 2)) %>%
        dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                      tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                      treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
        dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
        dplyr::select(sex, tissue, treatment, gene, samples, counts) %>%
        filter(tissue == whichtissue, sex %in% whichsex)  %>%
        drop_na()
      #candidates$treatment <- factor(candidates$treatment, levels = alllevels)
      return(candidates)
    }


    hypvsd <- getcandidatevsd(candidategenes, "hypothalamus", sexlevels)
    pitvsd <- getcandidatevsd(candidategenes, "pituitary", sexlevels)
    gonvsd <- getcandidatevsd(candidategenes, "gonad", sexlevels)
    candidatevsd <- rbind(hypvsd, pitvsd)
    candidatevsd <- rbind(candidatevsd, gonvsd)

    head(candidatevsd)

    ## # A tibble: 6 x 6
    ##   sex    tissue      treatment gene   samples                             counts
    ##   <chr>  <chr>       <chr>     <chr>  <chr>                                <dbl>
    ## 1 female hypothalam… control   ADRA2A L.G118_female_hypothalamus_control…   8.95
    ## 2 female hypothalam… control   ADRA2A R.G106_female_hypothalamus_control    8.81
    ## 3 female hypothalam… control   ADRA2A R.R20_female_hypothalamus_control.…   9.18
    ## 4 female hypothalam… control   ADRA2A R.R9_female_hypothalamus_control      8.72
    ## 5 female hypothalam… control   ADRA2A R.W44_female_hypothalamus_control     9.23
    ## 6 female hypothalam… prolong   ADRA2A blk.s031.pu.d_female_hypothalamus_…   9.13

All DEGs
--------

    DEG_path <- "../results/DEseq2/treatment/"   # path to the data
    DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
    DEG_pathfiles <- paste0(DEG_path, DEG_files)
    #DEG_files

    allDEG <- DEG_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name") %>% 
      mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/treatment/'), "[", 2))  %>% 
      mutate(DEG = sapply(strsplit(as.character(DEG),'_diffexp.csv'), "[", 1))  %>% 
      mutate(tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1)) %>%
      mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%
      mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%
      mutate(comparison = paste(down,up, sep = "_")) %>%
      mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
      mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%
    dplyr::select(sex,tissue,comparison, direction, gene, lfc, padj, logpadj) 
    head(allDEG)

    unique(allDEG$comparison)

    tempDEGs <- DEG_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name") %>% 
      mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/treatment/'), "[", 2))  %>% 
      mutate(DEG = sapply(strsplit(as.character(DEG),'_diffexp.csv'), "[", 1))  %>% 
      mutate(tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1)) %>%
      mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%
      mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%
      mutate(comparison = paste(down,up, sep = "_")) %>%
      mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
      mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2))
    head(tempDEGs)


    # sequential in 2 or more comparisons
    wideDEGscandidates <- tempDEGs %>%
      filter(gene %in% candidategenes) %>%
      mutate(posneg = ifelse(lfc >= 0, "+", "-"),
             sex = recode(sex, "female" = "F", "male" = "M" ),
             tissue = recode(tissue, 
                             "hypothalamus" = "H",
                             "pituitary" = "P", "gonad" = "G"),
             group = paste(sex, tissue, sep = "")) %>%
      mutate(res = paste("(", posneg, ")", sep = "")) %>%
      mutate(compres = paste(group, res, sep = "")) %>%
      group_by(gene, comparison) %>%
      summarize(results = str_c(compres, collapse = "; ")) %>%
      pivot_wider(names_from = comparison, values_from = results) 
    wideDEGscandidates

    unique(wideDEGscandidates$compres)

    ## Warning: Unknown or uninitialised column: 'compres'.

### manipuluation

    head(allDEG)

    ## # A tibble: 6 x 8
    ##   sex    tissue       comparison direction gene        lfc   padj logpadj
    ##   <chr>  <chr>        <chr>      <chr>     <chr>     <dbl>  <dbl>   <dbl>
    ## 1 female hypothalamus bldg_early early     PPP2R5C    2.03 0.0592    1.23
    ## 2 female hypothalamus bldg_early early     LOC769052  1.94 0.0413    1.38
    ## 3 female hypothalamus bldg_early early     PNO1       1.87 0.0429    1.37
    ## 4 female hypothalamus bldg_early early     UNKL       1.80 0.0426    1.37
    ## 5 female hypothalamus bldg_early early     SPAM1      1.63 0.0365    1.44
    ## 6 female hypothalamus bldg_early early     SLCO4C1    1.62 0.0883    1.05

    filteredDEGs <- allDEG %>%
      filter(lfc > 0.14 | lfc < -0.14 )
    dim(allDEG)

    ## [1] 327163      8

    dim(filteredDEGs)

    ## [1] 327163      8

    makemanipulationDEGtables <- function(whichcomparisons){
      df <- filteredDEGs %>%
        filter(comparison %in% whichcomparisons) %>%
         mutate(posneg = ifelse(lfc >= 0, "+", "-"),
             sex = recode(sex, "female" = "F", "male" = "M" ),
             tissue = recode(tissue, 
                             "hypothalamus" = "H",
                             "pituitary" = "P", "gonad" = "G"),
             group = paste(sex, tissue, sep = "")) %>%
      mutate(res = paste("(", posneg, ")", sep = "")) %>%
      mutate(compres = paste(comparison, res, sep = "")) %>%
      select(gene, group, compres) %>%
      group_by(group,gene) %>%
      summarize(results = str_c(compres, collapse = "; ")) %>%
      mutate(n = (str_count(results, pattern = ";"))+1) %>%
      arrange(desc(n))  %>%
      filter(n>1) %>%
      group_by(group,results) %>%
      summarize(genes = str_c(gene, collapse = "; "))  %>%
      #arrange(group, results) %>%
      mutate(n = (str_count(genes, pattern = ";"))+1) %>%
      arrange(desc(n))
      print(head(df))
    return(df) 
    }

    removalcomps <- c("inc.d3_m.inc.d3", "inc.d9_m.inc.d9", "inc.d17_m.inc.d17", "hatch_m.n2")
    earlys <- c("inc.d9_early", "inc.d17_early", "hatch_early", "n5_early")
    prolongs <- c("inc.d9_prolong","inc.d17_prolong", "hatch_prolong",  "n5_prolong")
    extends <- c("inc.d9_extend","inc.d17_extend", "hatch_extend",  "n5_extend")

    removalDEGs <- makemanipulationDEGtables(removalcomps)

    ## # A tibble: 0 x 4
    ## # Groups:   group [0]
    ## # … with 4 variables: group <chr>, results <chr>, genes <chr>, n <dbl>

    earlyDEGs <- makemanipulationDEGtables(earlys)

    ## # A tibble: 6 x 4
    ## # Groups:   group [3]
    ##   group results                        genes                                   n
    ##   <chr> <chr>                          <chr>                               <dbl>
    ## 1 FP    hatch_early(-); inc.d17_early… ABCC8; ABCC9; ABCE1; ABRACL; AC005…   669
    ## 2 MH    hatch_early(+); inc.d17_early… ABCC8; ABCC9; ABCG1; ACAP2; ACAP3;…   504
    ## 3 MH    hatch_early(-); inc.d17_early… ABHD17B; ACE2; ACOT12; ACYP2; ADIP…   353
    ## 4 FP    hatch_early(+); inc.d17_early… ACACA; ADCY9; AFF4; AGRN; AHDC1; A…   251
    ## 5 MP    hatch_early(+); inc.d17_early… A2ML3; AAK1; ABAT; ACCN2; ACSBG1; …   180
    ## 6 MP    hatch_early(-); inc.d17_early… AC005943.2; ADARB1; AGR2; AIMP1; A…   179

    prolongDEGs <- makemanipulationDEGtables(prolongs)

    ## # A tibble: 6 x 4
    ## # Groups:   group [2]
    ##   group  results                         genes                                 n
    ##   <chr>  <chr>                           <chr>                             <dbl>
    ## 1 Mgona… hatch_prolong(+); inc.d17_prol… ABLIM2; ACOX2; ADD3; AGBL2; AGPA…   363
    ## 2 FH     hatch_prolong(+); inc.d17_prol… ABCG1; ACVR2A; ADGRV1; ADRA2C; A…   149
    ## 3 Mgona… inc.d17_prolong(+); inc.d9_pro… ABCA3; ABCA5; ABCD3; ABCG4; ABHD…   130
    ## 4 FH     hatch_prolong(-); inc.d17_prol… ACYP2; ANLN; APOBEC2; APOH; ARHG…   107
    ## 5 Mgona… hatch_prolong(+); inc.d17_prol… AC113404.1; ACSF2; AEBP1; AHDC1;…   106
    ## 6 Mgona… hatch_prolong(+); inc.d9_prolo… ATRNL1; B3GNT4; C8H1ORF168; CCDC…    35

    extendDEGs <- makemanipulationDEGtables(extends)

    ## # A tibble: 6 x 4
    ## # Groups:   group [4]
    ##   group results                genes                                           n
    ##   <chr> <chr>                  <chr>                                       <dbl>
    ## 1 FP    hatch_extend(-); inc.… AC005943.2; AIMP2; ATOX1; AURKA; BOLA3; CE…    62
    ## 2 FH    hatch_extend(+); inc.… ADAM23; ADGRV1; ANKRD10; ASB14; C5ORF30; C…    44
    ## 3 MH    hatch_extend(-); inc.… CALB2; CBLN2; FZD10; GBX2; HAPLN4; HCRT; I…    27
    ## 4 MP    hatch_extend(-); inc.… AC005943.2; APITD1; BIRC5; CKS2; COX4I1; C…    27
    ## 5 FH    hatch_extend(-); inc.… ANXA6; C18ORF32; CALB2; CBLN2; GBX2; GPR17…    17
    ## 6 FP    inc.d17_extend(-); in… AEBP1; BTC; C1QTNF4; CALML4; CD34; CD99; C…    17

    manipDEGs <- rbind(removalDEGs, earlyDEGs) %>%
      rbind(., prolongDEGs) %>%
      rbind(., extendDEGs) %>%
      arrange(desc(n))
    head(manipDEGs)

    ## # A tibble: 6 x 4
    ## # Groups:   group [4]
    ##   group  results                         genes                                 n
    ##   <chr>  <chr>                           <chr>                             <dbl>
    ## 1 FP     hatch_early(-); inc.d17_early(… ABCC8; ABCC9; ABCE1; ABRACL; AC0…   669
    ## 2 MH     hatch_early(+); inc.d17_early(… ABCC8; ABCC9; ABCG1; ACAP2; ACAP…   504
    ## 3 Mgona… hatch_prolong(+); inc.d17_prol… ABLIM2; ACOX2; ADD3; AGBL2; AGPA…   363
    ## 4 MH     hatch_early(-); inc.d17_early(… ABHD17B; ACE2; ACOT12; ACYP2; AD…   353
    ## 5 FP     hatch_early(+); inc.d17_early(… ACACA; ADCY9; AFF4; AGRN; AHDC1;…   251
    ## 6 MP     hatch_early(+); inc.d17_early(… A2ML3; AAK1; ABAT; ACCN2; ACSBG1…   180

save files
----------

    write.csv(allDEG, "../results/03_allDEG.csv", row.names = F)
    write.csv(candidatevsd, "../results/03_candidatevsd.csv")
    #write.csv(wideDEGsSequential, "../results/table1.csv", row.names = F)
    write.csv(manipDEGs, "../results/manipDEGs.csv")
    write.csv(wideDEGscandidates, "../results/03_wideDEGscandidates.csv", row.names = F)

### removal overlap

    filteredDEGs <- allDEG %>%
      filter(lfc > 0.14 | lfc < -0.14 ) %>%
      filter(tissue == "hypothalamus")

    removalcomps <- c("inc.d3_m.inc.d3", "inc.d9_m.inc.d9", "inc.d17_m.inc.d17", "hatch_m.n2")

    hatch_m.n2DEGs <- filteredDEGs %>%
      filter(comparison == "hatch_m.n2") %>%
      drop_na() %>%
      pull(gene)

    inc.d17_m.inc.d17DEGs <- filteredDEGs %>%
      filter(comparison == "inc.d17_m.inc.d17") %>%
      drop_na() %>%
      pull(gene)

    inc.d3_m.inc.d3DEGs <- filteredDEGs %>%
      filter(comparison == "inc.d3_m.inc.d3") %>%
      drop_na() %>%
      pull(gene)

    inc.d9_m.inc.d9DEGs <- filteredDEGs %>%
      filter(comparison == "inc.d9_m.inc.d9") %>%
      drop_na() %>%
      pull(gene)


    intersect(hatch_m.n2DEGs, inc.d17_m.inc.d17DEGs) %>%
      intersect(.,inc.d3_m.inc.d3DEGs) %>%
      intersect(.,inc.d9_m.inc.d9DEGs)

    ## character(0)

    venn.diagram(
      x = list(inc.d3_m.inc.d3DEGs,hatch_m.n2DEGs,
                inc.d9_m.inc.d9DEGs, inc.d17_m.inc.d17DEGs),
      category.names = c("Inc 3" , "Hatch" , "Inc 9", " Inc 17"),
      filename = '../figures/venn-hyp.png',
      output=FALSE,
      print.mode = "raw",
      imagetype = "png",
      col=c("#CDCDCD", '#262625', '#959595','#959595'),
      fill = c("#CDCDCD", "#262625", "#959595", "#959595"),
      main = "Hypothalamus - Offspring Removal",
      sub = "Differentially expressed genes relative to internal controls"
      )

    ## [1] 1

chicks / nesting care
---------------------

    filteredDEGs <- allDEG %>%
      filter(lfc > 1.1 | lfc < -1.1 ) %>%
      filter(tissue == "hypothalamus")

    chicks <- c("bldg_hatch", "bldg_n5", "bldg_n9",
                "bldg_extend", "bldg_early")

    bldg_hatchDEGs <- filteredDEGs %>%
      filter(comparison == "bldg_hatch") %>%
      drop_na() %>%
      pull(gene)

    bldg_n5DEGs <- filteredDEGs %>%
      filter(comparison == "bldg_n5") %>%
      drop_na() %>%
      pull(gene)

    bldg_n9DEGs <- filteredDEGs %>%
      filter(comparison == "bldg_n9") %>%
      drop_na() %>%
      pull(gene)

    bldg_extendDEGs <- filteredDEGs %>%
      filter(comparison == "bldg_extend") %>%
      drop_na() %>%
      pull(gene)

    bldg_earlyDEGs <- filteredDEGs %>%
      filter(comparison == "bldg_early") %>%
      drop_na() %>%
      pull(gene)


    # chicks
    venn.diagram(
      x = list(bldg_n5DEGs, bldg_hatchDEGs, 
               bldg_earlyDEGs,bldg_extendDEGs,  bldg_n9DEGs),
      category.names = c( "N5" ,"Hatch" , "Early", " Extend", "N9"),
      filename = '../figures/venn-chicks-pit.png',
      output=FALSE,
      print.mode = "raw",
      imagetype = "png",
      height = 1500,
        width = 1750,
      resolution = 500,
      cex = 0.75,
      cat.cex = 0.75,
      col=c('#3182bd', "#6baed6",  '#cbc9e2', '#6a51a3', '#08519c'),
      fill = c("#3182bd", "#6baed6", "#cbc9e2", '#6a51a3', "#08519c" )
      )

    ## [1] 1

eggs / incubation
-----------------

    filteredDEGs <- allDEG %>%
      filter(lfc > 1.1 | lfc < -1.1 ) %>%
      filter(tissue == "hypothalamus") 

    eggs <- c("bldg_lay", 
             "bldg_inc.d3", "bldg_inc.d9", "bldg_inc.d17",
             "bldg_prolong")

    bldg_layDEGs <- filteredDEGs %>%
      filter(comparison == "bldg_lay") %>%
      drop_na() %>%
      pull(gene)

    bldg_inc.d3DEGs <- filteredDEGs %>%
      filter(comparison == "bldg_inc.d3") %>%
      drop_na() %>%
      pull(gene)

    bldg_inc.d9DEGs <- filteredDEGs %>%
      filter(comparison == "bldg_inc.d9") %>%
      drop_na() %>%
      pull(gene)

    bldg_inc.d17DEGs <- filteredDEGs %>%
      filter(comparison == "bldg_inc.d17") %>%
      drop_na() %>%
      pull(gene)

    bldg_prolongDEGs <- filteredDEGs %>%
      filter(comparison == "bldg_prolong") %>%
      drop_na() %>%
      pull(gene)

    # chicks
    venn.diagram(
      x = list(bldg_inc.d3DEGs, bldg_layDEGs,  bldg_prolongDEGs,
               bldg_inc.d17DEGs, bldg_inc.d9DEGs),
      category.names = c("Inc3" , "Lay" , "Prolong", " Inc17", "Inc9"),
      filename = '../figures/venn-eggs-pit.png',
      output=FALSE,
      print.mode = "raw",
      imagetype = "png",
      height = 1500,
        width = 1750,
      resolution = 500,
      cex = 0.75,
      cat.cex = 0.75,
      col=c( '#78c679', "#fed98e",  '#9e9ac8','#006837', '#31a354'),
      fill = c("#78c679", "#fed98e",  "#9e9ac8",'#006837', "#31a354")
      )

    ## [1] 1

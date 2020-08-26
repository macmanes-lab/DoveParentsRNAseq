    library(tidyverse)

    ## ── Attaching packages ───────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ──────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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

    ## [1] "female_gonad_vsd.csv"        "female_gonads_vsd.csv"      
    ## [3] "female_hypothalamus_vsd.csv" "female_pituitary_vsd.csv"   
    ## [5] "male_gonad_vsd.csv"          "male_gonads_vsd.csv"        
    ## [7] "male_hypothalamus_vsd.csv"   "male_pituitary_vsd.csv"

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

    wideDEGs <- tempDEGs %>%
      mutate(lfc = round(lfc, 2),
             padj =  round(padj,2)) %>%
      mutate(sextissue = paste(sex, tissue, sep = "_"),
             lfcpadj = paste(lfc, padj, sep = ", "))  %>%
      dplyr::select(sextissue,comparison, gene, lfcpadj) %>%
      pivot_wider(names_from = comparison, values_from = lfcpadj) 

    wideDEGscandidates <- wideDEGs %>%
      filter(gene %in% GOgenesLong$gene)
    wideDEGscandidates

    wideDEGsBldg <- wideDEGscandidates %>%
      select(sextissue,gene, contains("bldg"))
    wideDEGsBldg

    wideDEGsControl <- wideDEGscandidates %>%
      select(sextissue,gene, contains("control"))
    wideDEGsControl


    # sequential in 2 or more comparisons
    wideDEGsSequential <- wideDEGscandidates %>%
      select(sextissue,gene, "control_bldg", "bldg_lay", "lay_inc.d3", "inc.d3_inc.d9",
             "inc.d9_inc.d17", "inc.d17_hatch", "hatch_n5", "n5_n9") %>%
      mutate(totalNA = rowSums(is.na(.))) %>%
      filter(totalNA < 8) %>%
      arrange(totalNA, gene, sextissue) %>%
      #arrange( gene, sextissue) %>%
      select(-totalNA) %>% 
      replace(., is.na(.), "")
    wideDEGsSequential

### manipuluation

    head(allDEG)

    ## # A tibble: 6 x 8
    ##   sex    tissue comparison  direction gene           lfc     padj logpadj
    ##   <chr>  <chr>  <chr>       <chr>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  bldg_extend extend    CDK3         19.5  2.69e-20   19.6 
    ## 2 female gonad  bldg_extend extend    CRISP2        6.52 3.48e- 3    2.46
    ## 3 female gonad  bldg_extend extend    KRT20         5.47 3.67e- 4    3.44
    ## 4 female gonad  bldg_extend extend    CLDN34        5.01 5.14e- 3    2.29
    ## 5 female gonad  bldg_extend extend    LOC107049005  4.89 7.51e- 2    1.12
    ## 6 female gonad  bldg_extend extend    OMD           3.53 5.38e- 4    3.27

    filteredDEGs <- allDEG %>%
      filter(lfc > 1.1 | lfc < -1.1 )
    dim(allDEG)

    ## [1] 810213      8

    dim(filteredDEGs)

    ## [1] 123444      8

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
    early <- c("early_inc.d9", "hatch_early")
    late <- c("hatch_prolong", "inc.d17_prolong")


    removalDEGs <- makemanipulationDEGtables(removalcomps)

    ## # A tibble: 6 x 4
    ## # Groups:   group [4]
    ##   group results                  genes                                         n
    ##   <chr> <chr>                    <chr>                                     <dbl>
    ## 1 FP    hatch_m.n2(-); inc.d17_… ADAMTSL2; AGR2; B3GNTL2; C12ORF57; CD164…    59
    ## 2 FP    hatch_m.n2(+); inc.d17_… ADAMTS2; ANGPTL7; ANKRD34B; AQP3; ASIC4;…    46
    ## 3 FH    hatch_m.n2(+); inc.d17_… BET3L; BRS3; CCDC60; CPA6; CYP2J2L6; CYP…    43
    ## 4 MP    inc.d17_m.inc.d17(-); i… ANLN; APOD; APOH; AQP4; CDH19; CLDN11; C…    30
    ## 5 FH    hatch_m.n2(-); inc.d17_… ADPRHL1; C7ORF72; CALB2; CBLN2; CRYGS; F…    25
    ## 6 MH    hatch_m.n2(+); inc.d17_… AKAP5; BET3L; BRS3; CCDC60; ESRP2; FGF3;…    19

    hatchDEGs <- makemanipulationDEGtables(early)

    ## # A tibble: 6 x 4
    ## # Groups:   group [3]
    ##   group results               genes                                            n
    ##   <chr> <chr>                 <chr>                                        <dbl>
    ## 1 MP    early_inc.d9(-); hat… A2ML3; ACCN2; ACSBG1; ADAMTS19; ADCYAP1; AD…   270
    ## 2 MH    early_inc.d9(+); hat… ALDH1A1; CALB2; CALCA; CBLN2; CCK; CRYBB2; …    30
    ## 3 MP    early_inc.d9(+); hat… ANXA5; ATP2B4; C10H15ORF60; KYNU; LOC101747…    13
    ## 4 MH    early_inc.d9(-); hat… IGJ; IGLL1; LOC420300; PPP1R1B; SATB2; SH3R…     6
    ## 5 FP    early_inc.d9(+); hat… ANXA5; AVD; SLC9A3                               3
    ## 6 FP    early_inc.d9(-); hat… NINJ2; PDE6H                                     2

    prolongDEGs <- makemanipulationDEGtables(late)

    ## # A tibble: 6 x 4
    ## # Groups:   group [4]
    ##   group results                 genes                                          n
    ##   <chr> <chr>                   <chr>                                      <dbl>
    ## 1 MG    hatch_prolong(+); inc.… ABLIM2; ACOX2; ACSF2; ADAM19; ADAM33; ADA…   431
    ## 2 FG    hatch_prolong(+); inc.… ABCA4; ABCB11; ACOT12; ACR; ACRBP; ADAM20…   372
    ## 3 FH    hatch_prolong(-); inc.… ACP5; ADIRF; APOBEC2; AvBD5; BAIAP2L2; BL…    54
    ## 4 FH    hatch_prolong(+); inc.… ABCG8; ANKK1; ARSI; BRS3; CCDC60; CNTNAP4…    53
    ## 5 MH    hatch_prolong(-); inc.… APOBEC2; CALCA; FXYD2; GBX2; HCRT; LOC101…    13
    ## 6 FG    hatch_prolong(-); inc.… C11ORF34; C23H1ORF63; CDCA2; DPP9; HSP25;…    11

    manipDEGs <- rbind(removalDEGs, hatchDEGs) %>%
      rbind(., prolongDEGs) %>%
      arrange(desc(n))
    head(manipDEGs)

    ## # A tibble: 6 x 4
    ## # Groups:   group [5]
    ##   group results                 genes                                          n
    ##   <chr> <chr>                   <chr>                                      <dbl>
    ## 1 MG    hatch_prolong(+); inc.… ABLIM2; ACOX2; ACSF2; ADAM19; ADAM33; ADA…   431
    ## 2 FG    hatch_prolong(+); inc.… ABCA4; ABCB11; ACOT12; ACR; ACRBP; ADAM20…   372
    ## 3 MP    early_inc.d9(-); hatch… A2ML3; ACCN2; ACSBG1; ADAMTS19; ADCYAP1; …   270
    ## 4 FP    hatch_m.n2(-); inc.d17… ADAMTSL2; AGR2; B3GNTL2; C12ORF57; CD164L…    59
    ## 5 FH    hatch_prolong(-); inc.… ACP5; ADIRF; APOBEC2; AvBD5; BAIAP2L2; BL…    54
    ## 6 FH    hatch_prolong(+); inc.… ABCG8; ANKK1; ARSI; BRS3; CCDC60; CNTNAP4…    53

save files
----------

    write.csv(allDEG, "../results/03_allDEG.csv", row.names = F)
    write.csv(candidatevsd, "../results/03_candidatevsd.csv")
    write.csv(wideDEGsSequential, "../results/table1.csv", row.names = F)

    write.csv(removalDEGs, "../results/removalDEGs.csv")
    write.csv(hatchDEGs, "../results/hatchDEGs.csv")
    write.csv(prolongDEGs, "../results/prolongDEGs.csv")

    write.csv(manipDEGs, "../results/manipDEGs.csv")

### overlapp

    removalcomps <- c("inc.d3_m.inc.d3", "inc.d9_m.inc.d9", "inc.d17_m.inc.d17", "hatch_m.n2")
    early <- c("early_inc.d9", "hatch_early")
    late <- c("hatch_prolong", "inc.d17_prolong")


    filteredDEGs %>%
      filter(comparison %in% removalcomps)  %>%
      select(sex:gene) %>%
      pivot_wider(names_from = comparison, values_from = direction)

    ## # A tibble: 1,847 x 7
    ##    sex   tissue gene  hatch_m.n2 inc.d17_m.inc.d… inc.d3_m.inc.d3
    ##    <chr> <chr>  <chr> <chr>      <chr>            <chr>          
    ##  1 fema… gonad  IQCG  m.n2       <NA>             <NA>           
    ##  2 fema… gonad  LOC1… m.n2       <NA>             m.inc.d3       
    ##  3 fema… gonad  PI15  m.n2       m.inc.d17        m.inc.d3       
    ##  4 fema… gonad  REG4  m.n2       <NA>             <NA>           
    ##  5 fema… gonad  MMP13 m.n2       <NA>             <NA>           
    ##  6 fema… gonad  SLC2… m.n2       <NA>             <NA>           
    ##  7 fema… gonad  GMNC  m.n2       m.inc.d17        <NA>           
    ##  8 fema… gonad  AGT   m.n2       <NA>             <NA>           
    ##  9 fema… gonad  JPH2  m.n2       <NA>             <NA>           
    ## 10 fema… gonad  PHGDH m.n2       <NA>             <NA>           
    ## # … with 1,837 more rows, and 1 more variable: inc.d9_m.inc.d9 <chr>

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
      


     
    # Generate 3 sets of 200 words
    # Chart

    venn.diagram(
      x = list(hatch_m.n2DEGs, inc.d17_m.inc.d17DEGs, 
               inc.d3_m.inc.d3DEGs, inc.d9_m.inc.d9DEGs),
      category.names = c("Inc 3" , "Hatch" , "Inc 9", " Inc 17"),
      filename = '../figures/venn.png',
      output=FALSE,
      col=c("#CDCDCD", '#262625', '#959595','#959595'),
      fill = c("#CDCDCD", "#262625", "#959595", "#959595"),
      main = "No. of genes that respond to offspring removal",
      sub = "9 shared genes: \nPHGDH, NRG2, LHX9, SSTR3, VTN, SARM1, ALK, LPAR1, and APOH" )

    ## [1] 1

    intersect(hatch_m.n2DEGs, inc.d17_m.inc.d17DEGs) %>%
      intersect(.,inc.d3_m.inc.d3DEGs) %>%
      intersect(.,inc.d9_m.inc.d9DEGs)

    ## [1] "PHGDH" "NRG2"  "LHX9"  "SSTR3" "VTN"   "SARM1" "ALK"   "LPAR1" "APOH"

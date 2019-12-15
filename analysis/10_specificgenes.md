analyses downstream of calculating tissue specific variance in gene expression
==============================================================================

    vsd.hyp <- readvsd("../results/04_vsd_hyp.csv")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )

    ## See spec(...) for full column specifications.

    colData.hyp <- readcolData("../results/04_colData_hyp.csv")

    ## Parsed with column specification:
    ## cols(
    ##   V1 = col_character(),
    ##   sex = col_character(),
    ##   treatment = col_character(),
    ##   sextissue = col_character(),
    ##   hypothesis = col_character()
    ## )

    vsd.pit <- readvsd("../results/04_vsd_pit.csv")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    colData.pit <- readcolData("../results/04_colData_pit.csv")

    ## Parsed with column specification:
    ## cols(
    ##   V1 = col_character(),
    ##   sex = col_character(),
    ##   treatment = col_character(),
    ##   sextissue = col_character(),
    ##   hypothesis = col_character()
    ## )

    vsd.gon <- readvsd("../results/04_vsd_gon.csv")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    colData.gon <- readcolData("../results/04_colData_gon.csv")

    ## Parsed with column specification:
    ## cols(
    ##   V1 = col_character(),
    ##   sex = col_character(),
    ##   treatment = col_character(),
    ##   sextissue = col_character(),
    ##   hypothesis = col_character()
    ## )

    colData.hyp <- colData.hyp %>% 
      mutate(tissue = sapply(strsplit(sample,'\\_'), "[", 3)) %>%
      mutate(bird = sapply(strsplit(sample,'\\_'), "[", 1)) %>%
      select(sample, bird, sex, tissue, treatment)

    colData.pit <- colData.pit %>%  
      mutate(tissue = sapply(strsplit(sample,'\\_'), "[", 3)) %>%
      mutate(bird = sapply(strsplit(sample,'\\_'), "[", 1)) %>%
      select(sample, bird, sex, tissue, treatment)

    colData.gon <- colData.gon %>%  
      mutate(tissue = sapply(strsplit(sample,'\\_'), "[", 3)) %>%
      mutate(bird = sapply(strsplit(sample,'\\_'), "[", 1)) %>%
      select(sample, bird, sex, tissue, treatment)

selecting candidate genes counts from the hypothalamus
======================================================

    geneinfo <- read_csv("../metadata/00_geneinfo.csv") %>%  dplyr::select(Name, geneid, entrezid) %>% arrange(Name)

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_character(),
    ##   row.names = col_double(),
    ##   Name = col_character(),
    ##   geneid = col_double(),
    ##   entrezid = col_character()
    ## )

    head(geneinfo)

    ## # A tibble: 6 x 3
    ##   Name      geneid entrezid      
    ##   <chr>      <dbl> <chr>         
    ## 1 A2ML1     418254 XP_015148230.1
    ## 2 A2ML2     427942 XP_004938161.2
    ## 3 A2ML3  100857394 XP_015148584.1
    ## 4 A2ML4  100858010 XP_015154891.1
    ## 5 A4GALT    418223 XP_015145932.1
    ## 6 A4GNT     429136 XP_426692.3

    candidategenes <- c("PRL","PRLR", "CRY1")

    candidates.hyp <- selectcandidatevsds(candidategenes, vsd.hyp, colData.hyp)

    ## [1] "PRL"  "PRLR" "CRY1"
    ## [1] "NP_989576.1"    "NP_990797.2"    "XP_015132722.1"

    ## Joining, by = "entrezid"

    ## Joining, by = "sample"

    candidates.pit <- selectcandidatevsds(candidategenes, vsd.pit, colData.pit)

    ## [1] "PRL"  "PRLR" "CRY1"
    ## [1] "NP_989576.1"    "NP_990797.2"    "XP_015132722.1"

    ## Joining, by = "entrezid"
    ## Joining, by = "sample"

    candidates.gon <- selectcandidatevsds(candidategenes, vsd.gon, colData.gon)

    ## [1] "PRL"  "PRLR" "CRY1"
    ## [1] "NP_989576.1"    "NP_990797.2"    "XP_015132722.1"

    ## Joining, by = "entrezid"
    ## Joining, by = "sample"

    plotcanddateexpression <- function(candidateexpression,  mysubtitle, whichgene, myylab){
      
      ggplot(candidateexpression, aes(x = as.numeric(treatment), y = whichgene)) + 
            geom_smooth(aes(colour = sex)) +
        geom_boxplot(aes(fill = treatment, alpha = sex, color = sex)) + 
        scale_alpha_manual(values = c(0.75,1)) +
         theme_B3() +
        theme(legend.position = "none") +
        theme(axis.title.y=element_text(face="italic"),
              axis.title.x = element_blank(),
              axis.text.x = element_blank()) +
        scale_color_manual(values = c("female" = "#969696", "male" = "#525252")) +
        labs(subtitle = mysubtitle, y = myylab)
      
    }

    a <- plotcanddateexpression(candidates.hyp,  "hypothalamus", candidates.hyp$PRL, "PRL")
    b <- plotcanddateexpression(candidates.pit, "pituitary", candidates.pit$PRL, "PRL")
    c <- plotcanddateexpression(candidates.gon,  "gonad", candidates.gon$PRL, "PRL")

    d <- plotcanddateexpression(candidates.hyp, "hypothalamus", candidates.hyp$PRLR, "PRLR")
    e <- plotcanddateexpression(candidates.pit, "pituitary", candidates.pit$PRLR, "PRLR")
    f <- plotcanddateexpression(candidates.gon, "gonad", candidates.gon$PRLR, "PRLR")


    g <- plotcanddateexpression(candidates.hyp, "hypothalamus", candidates.hyp$CRY1, "CRY1")
    h <- plotcanddateexpression(candidates.pit, "pituitary", candidates.pit$CRY1, "CRY1")
    i <- plotcanddateexpression(candidates.gon, "gonad", candidates.gon$CRY1, "CRY1")

    abc <- plot_grid(a, d, g, nrow = 1)

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

    def <- plot_grid(b, e, h, nrow = 1)

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

    ghi <- plot_grid(c, f, i, nrow = 1)

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

    plot_grid(abc, def,ghi,nrow = 3)

![](../figures/specificgenes/candidategenes-1.png)

    PRLpit <- candidates.pit %>%
      select(bird, sex,treatment,tissue, PRL)
    head(PRLpit)

    ##      bird    sex treatment    tissue      PRL
    ## 1 L.Blu13   male   control pituitary 16.75792
    ## 2  L.G107   male   control pituitary 17.87347
    ## 3  L.G118 female   control pituitary 17.91121
    ## 4    L.R3   male   control pituitary 16.51080
    ## 5    L.R8   male   control pituitary 17.30663
    ## 6   L.W33   male   control pituitary 19.93211

    write.csv(PRLpit, "../results/10_PRLpit.csv", row.names = F)

    pitcandidates <- c("BMP6", "CGA", "GAL", "XBP1")
    goncandidates <- c("CCND1", "NME1", "PRL", "APLP2", "CYP11A1", "CREBRF", "DBH")

    candidates.pit <- selectcandidatevsds(pitcandidates, vsd.pit, colData.pit) %>% select(-sample) %>% filter(treatment %in% c("bldg", "lay"))

    ## [1] "BMP6" "CGA"  "GAL"  "XBP1"
    ## [1] "XP_015131483.1" "NP_001264950.1" "NP_001138861.1" "NP_001006192.2"

    ## Joining, by = "entrezid"

    ## Joining, by = "sample"

    candidates.gon <- selectcandidatevsds(goncandidates, vsd.gon, colData.gon) %>% select(-sample) %>% filter(treatment %in% c("bldg", "lay"))

    ## [1] "CCND1"   "NME1"    "PRL"     "APLP2"   "CYP11A1" "CREBRF"  "DBH"    
    ## [1] "NP_001006317.2" "NP_990712.1"    "XP_001231574.1" "NP_001001756.1"
    ## [5] "XP_415429.5"    "NP_001191690.1" "NP_990797.2"

    ## Joining, by = "entrezid"
    ## Joining, by = "sample"

    candidates.pit <- candidates.pit %>% rename_at(vars(-bird,-sex, -tissue, -treatment), function(x) paste0("pit_", x)) %>% select(-tissue)
    candidates.gon <- candidates.gon %>% rename_at(vars(-bird,-sex, -tissue, -treatment), function(x) paste0("gon_", x)) %>% select(-tissue)

    counts <- full_join(candidates.pit, candidates.gon) %>% drop_na() %>% select(-bird, -sex, -treatment)

    ## Joining, by = c("bird", "sex", "treatment")

    cols <- full_join(candidates.pit, candidates.gon) %>% drop_na() %>% select(bird, sex, treatment)

    ## Joining, by = c("bird", "sex", "treatment")

    x <- correlate(counts)

    ## 
    ## Correlation method: 'pearson'
    ## Missing treated using: 'pairwise.complete.obs'

    rplot(x, colors = c("skyblue1", "white", "indianred2"))  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(subtitle = "bldg to lay")

    ## Don't know how to automatically pick scale for object of type noquote. Defaulting to continuous.

![](../figures/specificgenes/correlations-bldg-lay-1.png)

    x %>% network_plot.cor_df(min_cor = .6 , colors = c("skyblue1", "white", "indianred2")) + 
      labs(subtitle = "bldg to lay")

![](../figures/specificgenes/correlations-bldg-lay-2.png)

    pitcandidates <- c("PTEN", "CGA", "GAL", "XBP1")
    goncandidates <- c(
      "AVPR1A", "C1QTNF1", "CRY1", "CRY2", "F2R", "GJA1", "INHBA", "RHOA", "SMAD4",
      "ATP2B2", "ATP7B", "CAV1", "EIF2AK3", "ERBB4", "GJA1", "HIF1A", "MED1", "NME1", 
      "PRLR", "VEGFA", "XBP1",  "ZBTB7B", "APLP2", "NCOA1", "NCOA2", "THRB",
      "CREBRF", "GNAQ", "NR3C1")

    candidates.pit <- selectcandidatevsds(pitcandidates, vsd.pit, colData.pit) %>% select(-sample) %>% filter(treatment %in% c("inc.d3", "lay"))

    ## [1] "PTEN" "CGA"  "GAL"  "XBP1"
    ## [1] "NP_001264950.1" "NP_001138861.1" "XP_015134187.1" "NP_001006192.2"

    ## Joining, by = "entrezid"

    ## Joining, by = "sample"

    candidates.gon <- selectcandidatevsds(goncandidates, vsd.gon, colData.gon) %>% select(-sample) %>% filter(treatment %in% c("inc.d3", "lay"))

    ##  [1] "AVPR1A"  "C1QTNF1" "CRY1"    "CRY2"    "F2R"     "GJA1"    "INHBA"  
    ##  [8] "RHOA"    "SMAD4"   "ATP2B2"  "ATP7B"   "CAV1"    "EIF2AK3" "ERBB4"  
    ## [15] "GJA1"    "HIF1A"   "MED1"    "NME1"    "PRLR"    "VEGFA"   "XBP1"   
    ## [22] "ZBTB7B"  "APLP2"   "NCOA1"   "NCOA2"   "THRB"    "CREBRF"  "GNAQ"   
    ## [29] "NR3C1"  
    ##  [1] "NP_001006317.2" "XP_015148963.1" "XP_015131533.1" "NP_001103908.1"
    ##  [5] "XP_001231907.1" "NP_001099134.1" "XP_001231574.1" "NP_989576.1"   
    ##  [9] "NP_989575.1"    "XP_420868.2"    "XP_015144603.1" "XP_004937394.1"
    ## [13] "NP_989917.1"    "NP_001026598.1" "NP_989628.1"    "XP_015136806.1"
    ## [17] "XP_418125.3"    "XP_015140591.1" "XP_015138381.1" "NP_001191690.1"
    ## [21] "XP_015149519.1" "XP_015132722.1" "NP_990035.1"    "XP_015154049.1"
    ## [25] "NP_001239150.1" "NP_990373.1"    "NP_001006192.2" "XP_015153981.1"

    ## Joining, by = "entrezid"
    ## Joining, by = "sample"

    candidates.pit <- candidates.pit %>% rename_at(vars(-bird,-sex, -tissue, -treatment), function(x) paste0("pit_", x)) %>% select(-tissue)
    candidates.gon <- candidates.gon %>% rename_at(vars(-bird,-sex, -tissue, -treatment), function(x) paste0("gon_", x)) %>% select(-tissue)

    counts <- full_join(candidates.pit, candidates.gon) %>% drop_na() %>% select(-bird, -sex, -treatment)

    ## Joining, by = c("bird", "sex", "treatment")

    cols <- full_join(candidates.pit, candidates.gon) %>% drop_na() %>% select(bird, sex, treatment)

    ## Joining, by = c("bird", "sex", "treatment")

    x <- correlate(counts)

    ## 
    ## Correlation method: 'pearson'
    ## Missing treated using: 'pairwise.complete.obs'

    rplot(x, colors = c("skyblue1", "white", "indianred2")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(subtitle = "lay to inc.d3")

    ## Don't know how to automatically pick scale for object of type noquote. Defaulting to continuous.

![](../figures/specificgenes/correlations-lay-incd3-1.png)

    x %>% network_plot.cor_df(min_cor = .6 , colors = c("skyblue1", "white", "indianred2"))  + 
      labs(subtitle = "lay to inc.d3")

![](../figures/specificgenes/correlations-lay-incd3-2.png)

    pitcandidates <- c("ATP6AP2", "CRY1", "TMF1",
                       "EIF2AK3", "PRL", "SOCS2", "STAT5B", "XBP1",
                       "APP", "HDAC4", "NCOA1")
    goncandidates <- c("AVPR2", "BMP6", "CRHR2", "DRD5", "TAC1",
                       "SOCS2", "ZBTB7B",
                       "DRD5", "NCOA2", "BRINP1", "NR3C1")

    candidates.pit <- selectcandidatevsds(pitcandidates, vsd.pit, colData.pit) %>% select(-sample) %>% filter(treatment %in% c("inc.d9", "inc.d17"))

    ##  [1] "ATP6AP2" "CRY1"    "TMF1"    "EIF2AK3" "PRL"     "SOCS2"   "STAT5B" 
    ##  [8] "XBP1"    "APP"     "HDAC4"   "NCOA1"  
    ##  [1] "XP_015154580.1" "NP_001025972.2" "NP_989576.1"    "XP_420868.2"   
    ##  [5] "XP_015144578.1" "XP_015140591.1" "NP_990797.2"    "NP_989871.1"   
    ##  [9] "XP_015155078.1" "XP_423749.2"    "NP_001006192.2"

    ## Joining, by = "entrezid"

    ## Joining, by = "sample"

    candidates.gon <- selectcandidatevsds(goncandidates, vsd.gon, colData.gon) %>% select(-sample) %>% filter(treatment %in% c("inc.d9", "inc.d17"))

    ##  [1] "AVPR2"  "BMP6"   "CRHR2"  "DRD5"   "TAC1"   "SOCS2"  "ZBTB7B"
    ##  [8] "DRD5"   "NCOA2"  "BRINP1" "NR3C1" 
    ##  [1] "NP_001026650.1" "XP_015131483.1" "NP_989780.1"    "NP_989785.1"   
    ##  [5] "XP_015141299.1" "XP_015138381.1" "XP_015149519.1" "NP_989871.1"   
    ##  [9] "XP_004939375.1" "XP_015153981.1"

    ## Joining, by = "entrezid"
    ## Joining, by = "sample"

    candidates.pit <- candidates.pit %>% rename_at(vars(-bird,-sex, -tissue, -treatment), function(x) paste0("pit_", x)) %>% select(-tissue)
    candidates.gon <- candidates.gon %>% rename_at(vars(-bird,-sex, -tissue, -treatment), function(x) paste0("gon_", x)) %>% select(-tissue)

    counts <- full_join(candidates.pit, candidates.gon) %>% drop_na() %>% select(-bird, -sex, -treatment)

    ## Joining, by = c("bird", "sex", "treatment")

    cols <- full_join(candidates.pit, candidates.gon) %>% drop_na() %>% select(bird, sex, treatment)

    ## Joining, by = c("bird", "sex", "treatment")

    x <- correlate(counts)

    ## 
    ## Correlation method: 'pearson'
    ## Missing treated using: 'pairwise.complete.obs'

    rplot(x, colors = c("skyblue1", "white", "indianred2"))  + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(subtitletitle = "inc.d9 to inc.d17")

    ## Don't know how to automatically pick scale for object of type noquote. Defaulting to continuous.

![](../figures/specificgenes/correlations-incd9-incd17-1.png)

    x %>% network_plot.cor_df(min_cor = .6 , colors = c("skyblue1", "white", "indianred2")) + 
      labs(subtitletitle = "inc.d9 to inc.d17")

![](../figures/specificgenes/correlations-incd9-incd17-2.png)

    pitcandidates <- c( "ATP6AP2", "AVPR2", "CRY2", "RHOA", 
                        "EIF2AK3", "NME1", "PRL", "XBP1",
                        "APP", "NCOA1")
    hypcandidates <- c("CRHBP", "CRHR2", "DRD3", "EDN1", "FAM129B", "FOXL2", "KRAS",
                       "OPRK1", "POMC", "RAB11FIP5", "RUNX1",
                       "CAV1", "ZBTB7B", "HEXB", "PGR")

    candidates.pit <- selectcandidatevsds(pitcandidates, vsd.pit, colData.pit) %>% select(-sample) %>% filter(treatment %in% c("hatch", "n5"))

    ##  [1] "ATP6AP2" "AVPR2"   "CRY2"    "RHOA"    "EIF2AK3" "NME1"    "PRL"    
    ##  [8] "XBP1"    "APP"     "NCOA1"  
    ##  [1] "XP_015154580.1" "NP_001025972.2" "NP_001026650.1" "NP_989575.1"   
    ##  [5] "XP_420868.2"    "XP_015140591.1" "NP_001191690.1" "NP_990797.2"   
    ##  [9] "NP_990035.1"    "NP_001006192.2"

    ## Joining, by = "entrezid"

    ## Joining, by = "sample"

    candidates.hyp <- selectcandidatevsds(hypcandidates, vsd.hyp, colData.hyp) %>% select(-sample) %>% filter(treatment %in% c("hatch", "n5"))

    ##  [1] "CRHBP"     "CRHR2"     "DRD3"      "EDN1"      "FAM129B"  
    ##  [6] "FOXL2"     "KRAS"      "OPRK1"     "POMC"      "RAB11FIP5"
    ## [11] "RUNX1"     "CAV1"      "ZBTB7B"    "HEXB"      "PGR"      
    ##  [1] "NP_001099134.1" "XP_003643006.2" "NP_989785.1"    "XP_004938262.1"
    ##  [5] "XP_418943.2"    "XP_015135170.1" "NP_001012630.1" "XP_424791.3"   
    ##  [9] "XP_015145758.1" "XP_426087.2"    "NP_990593.1"    "NP_001026269.1"
    ## [13] "XP_420890.5"    "XP_015154785.1" "XP_015153981.1"

    ## Joining, by = "entrezid"
    ## Joining, by = "sample"

    candidates.pit <- candidates.pit %>% rename_at(vars(-bird,-sex, -tissue, -treatment), function(x) paste0("pit_", x)) %>% select(-tissue)
    candidates.hyp <- candidates.hyp %>% rename_at(vars(-bird,-sex, -tissue, -treatment), function(x) paste0("hyp_", x)) %>% select(-tissue)

    counts <- full_join(candidates.pit, candidates.hyp) %>% drop_na() %>% select(-bird, -sex, -treatment)

    ## Joining, by = c("bird", "sex", "treatment")

    cols <- full_join(candidates.pit, candidates.hyp) %>% drop_na() %>% select(bird, sex, treatment)

    ## Joining, by = c("bird", "sex", "treatment")

    x <- correlate(counts)

    ## 
    ## Correlation method: 'pearson'
    ## Missing treated using: 'pairwise.complete.obs'

    rplot(x,  colors = c("skyblue1", "white", "indianred2")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      labs(subtitle = "hatch to n5")

    ## Don't know how to automatically pick scale for object of type noquote. Defaulting to continuous.

![](../figures/specificgenes/correlations-hatch-n5-1.png)

    x %>% network_plot.cor_df(min_cor = .6, colors = c("skyblue1", "white", "indianred2"))  + 
      labs(subtitle = "hatch to n5")

![](../figures/specificgenes/correlations-hatch-n5-2.png)

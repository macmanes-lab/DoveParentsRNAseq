selecting candidate genes counts from the hypothalamus
======================================================

    # import "colData" which contains sample information and "countData" which contains read counts
    c.colData <- read.csv("../metadata/00_colData_characterization.csv", header = T, row.names = 1)
    c.countData <- read.csv("../results/00_countData_characterization.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)

    # select only hypothalamus samples
    hypsamples <- c.colData %>% 
      filter(tissue == "hypothalamus") %>%
      select(-study) %>%
      droplevels()
    countstokeep <- as.character(hypsamples$V1)
    hypcounts <- c.countData %>% 
      mutate(genes = row.names(c.countData)) %>% 
      select(genes, countstokeep)

    # select genes to keep 
    # CISH, SOCS1, SOCS2, SOCS2, SOCS4, SOCS5, SOCS6
    candidates <- c("NP_989957.1", "NP_001131120.1", "NP_989871.1", "NP_001186037.1", "NP_001120786.1")
    hypcounts <- hypcounts %>% 
      filter(genes %in% candidates)


    # transform and combine into 1 dataframe
    hypcounts <- as.data.frame(hypcounts)
    row.names(hypcounts) <- hypcounts$genes
    hypcounts$genes <- NULL
    hypcounts <- as.data.frame(t(hypcounts))
    hypcounts$V1 <- row.names(hypcounts)
    hypcandidatecounts <- left_join(hypsamples, hypcounts)

    ## Joining, by = "V1"

    ## Warning: Column `V1` joining factor and character vector, coercing into
    ## character vector

    # rename columsn
    colnames(hypcandidatecounts)[colnames(hypcandidatecounts)=="NP_989957.1"] <- "CISH"
    colnames(hypcandidatecounts)[colnames(hypcandidatecounts)=="NP_001131120.1"] <- "SOCS1"
    colnames(hypcandidatecounts)[colnames(hypcandidatecounts)=="NP_989871.1"] <- "SOCS2"
    colnames(hypcandidatecounts)[colnames(hypcandidatecounts)=="NP_001186037.1"] <- "SOCS4"
    colnames(hypcandidatecounts)[colnames(hypcandidatecounts)=="NP_001120786.1"] <- "SOCS5"
    colnames(hypcandidatecounts)[colnames(hypcandidatecounts)=="NP_001120784.1"] <- "SOCS6"

    # set rownames
    row.names(hypcandidatecounts) <- hypcandidatecounts$V1
    hypcandidatecounts$V1 <- NULL

    # preview and write
    head(hypcandidatecounts)

    ##                                            bird    sex       tissue
    ## L.Blu13_male_hypothalamus_control.NYNO  L.Blu13   male hypothalamus
    ## L.G107_male_hypothalamus_control         L.G107   male hypothalamus
    ## L.G118_female_hypothalamus_control.NYNO  L.G118 female hypothalamus
    ## L.R3_male_hypothalamus_control             L.R3   male hypothalamus
    ## L.R8_male_hypothalamus_control             L.R8   male hypothalamus
    ## L.W33_male_hypothalamus_control.NYNO      L.W33   male hypothalamus
    ##                                         treatment
    ## L.Blu13_male_hypothalamus_control.NYNO    control
    ## L.G107_male_hypothalamus_control          control
    ## L.G118_female_hypothalamus_control.NYNO   control
    ## L.R3_male_hypothalamus_control            control
    ## L.R8_male_hypothalamus_control            control
    ## L.W33_male_hypothalamus_control.NYNO      control
    ##                                                               group SOCS5
    ## L.Blu13_male_hypothalamus_control.NYNO    male.hypothalamus.control    92
    ## L.G107_male_hypothalamus_control          male.hypothalamus.control   196
    ## L.G118_female_hypothalamus_control.NYNO female.hypothalamus.control   138
    ## L.R3_male_hypothalamus_control            male.hypothalamus.control   333
    ## L.R8_male_hypothalamus_control            male.hypothalamus.control    23
    ## L.W33_male_hypothalamus_control.NYNO      male.hypothalamus.control   325
    ##                                         SOCS1 SOCS4 SOCS2 CISH
    ## L.Blu13_male_hypothalamus_control.NYNO      1     7     6   26
    ## L.G107_male_hypothalamus_control            6    18    16   88
    ## L.G118_female_hypothalamus_control.NYNO     1    11    10    0
    ## L.R3_male_hypothalamus_control              5    34    17   92
    ## L.R8_male_hypothalamus_control              1     4     3    7
    ## L.W33_male_hypothalamus_control.NYNO        1    23    22   79

    write.csv(hypcandidatecounts, "../results/10_hypcandidatecounts.csv")

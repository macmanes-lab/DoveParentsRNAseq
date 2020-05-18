    # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

    library(tidyverse)

    source("../R/themes.R")

    knitr::opts_chunk$set(message=F, warning=F)

Data Wrangling
==============

import Kallisto transcript data, make gene info file
----------------------------------------------------

    # import count data, set rows at entreziz
    kallistodata <- read.table("../results/kallistocounts.txt", 
                            sep = ",", row.names = NULL) %>%
      rename("NCBI" = "entrezid",
             "gene" = "Name")
    #head(kallistodata)

    ## save gene information
    geneinfo <- kallistodata %>%
      select(gene, geneid, NCBI) 
    head(geneinfo)

    ##       gene geneid           NCBI
    ## 1    EDNRB 408082 NP_001001127.1
    ## 2  CYP26A1 408183 NP_001001129.1
    ## 3    CFDP1 374073 NP_001001189.1
    ## 4    AvBD7 407777 NP_001001194.1
    ## 5     KRT5 407779 NP_001001195.1
    ## 6 HSD11B1L 408034 NP_001001201.1

\# check for genes that have multple transcripts expressed
----------------------------------------------------------

    isoforms <- kallistodata %>%
      group_by(gene) %>%
      summarize(n = n()) %>%
      filter(n > 1) %>%
      group_by(n) %>% 
      summarize(genes = str_c(gene, collapse = ", ")) %>%
      arrange(desc(n)) %>%
      rename("n(isoforms)" = "n") %>% 
      mutate("n(counts)" = str_count(genes, ",") + 1)
    head(isoforms)

    ## # A tibble: 6 x 3
    ##   `n(isoforms)` genes                                           `n(counts)`
    ##           <int> <chr>                                                 <dbl>
    ## 1            57 SSPO                                                      1
    ## 2             6 CEND1, LOC101747901, TIRAP                                3
    ## 3             5 AIP, CCL5, GINS3, SORBS2, TCF7L2                          5
    ## 4             4 ATG13, BPTF, DOCK3, DTNA, EIF4G3, ELAVL2, EPB4…          14
    ## 5             3 AGAP3, ARFIP1, ARHGAP17, BZRAP1, C2CD5, CACNA1…          70
    ## 6             2 ABCB10, ABCC3, ABCC6, ABCC9, ABCG8, ABI1, ABI2…         697

    # for example, these gene have 2 and 3 isoforms
    geneinfo %>% filter(gene == "GRIN1")

    ##    gene geneid           NCBI
    ## 1 GRIN1 404296    NP_996862.1
    ## 2 GRIN1 404296 XP_015134834.1

    geneinfo %>% filter(gene == "CACNA1C")

    ##      gene geneid           NCBI
    ## 1 CACNA1C 395891 XP_015141927.1
    ## 2 CACNA1C 395891 XP_015141993.1
    ## 3 CACNA1C 395891 XP_015142129.1

group data by gene
------------------

    #head(kallistodata)
    # aggregate transcript counts to gene counts
    countData <- kallistodata %>% 
      select(-row.names, -geneid, -NCBI) %>% 
      pivot_longer(-gene, names_to = "samples", values_to = "counts") %>%  
      pivot_wider(
        names_from = samples, 
        values_from = counts,
        values_fn = list(counts = sum))  %>% 
      arrange(gene) %>% 
      filter(gene != "")
    countData <- as.data.frame(countData)
    row.names(countData) <- countData$gene ## make gene the row name
    countData[1] <- NULL ## make gene the row name
    countData <- round(countData) #round all value to nearest 1s place
    head(countData[13:15])

    ##        L.R8_male_gonad_control L.R8_male_hypothalamus_control
    ## A2ML1                       79                              7
    ## A2ML2                        2                              0
    ## A2ML3                      233                            211
    ## A2ML4                       12                              0
    ## A4GALT                       9                              0
    ## A4GNT                        0                              0
    ##        L.R8_male_pituitary_control
    ## A2ML1                           38
    ## A2ML2                            0
    ## A2ML3                           73
    ## A2ML4                            7
    ## A4GALT                          63
    ## A4GNT                            0

    # print tolal num of genes and samples
    dim(countData)

    ## [1] 13966   987

wrangle colData
---------------

    colData <- read.table(file.path( "../metadata/kallistosamples.txt"),
                          header = F, stringsAsFactors = F) %>%
      # use strsplit to cut the filename into meaningful columns
      mutate(bird = sapply(strsplit(V1,'\\_'), "[", 1),
             sex = sapply(strsplit(V1,'\\_'), "[", 2),
             tissue = sapply(strsplit(V1,'\\_'), "[", 3),
             temp = sapply(strsplit(V1,'\\_'), "[", 4)) %>%
      mutate(treatmenttemp = sapply(strsplit(temp,'\\.'), "[", 1),
             NYNO = sapply(strsplit(temp,'\\.'), "[", 2)) %>%
     
      # rename variables
      mutate(treatment = ifelse(grepl("extend-hatch", treatmenttemp), "extend",
                                ifelse(grepl("inc-prolong", treatmenttemp), "prolong",
                                       ifelse(grepl("m.hatch", treatmenttemp), "m.n2",
                                              ifelse(grepl("m.inc.d8", treatmenttemp), "early",
                                       treatmenttemp))))) %>%
       select(-temp, -NYNO, -treatmenttemp ) %>%
      # replace dashes with periods (to make deseq happy)
      mutate(bird = gsub("-", ".", bird),
             treatment = gsub("-", ".", treatment),
             V1 = gsub("-", ".", V1))
    # specify levels that will be factors
    cols <- c("sex", "treatment", "tissue")
    colData[cols] <- lapply(colData[cols], factor) 
    colData <- colData %>%
      mutate(tissue = factor(tissue, levels = c("hypothalamus", "pituitary", "gonad"))) %>% 
      mutate(tissue = fct_recode(tissue, "gonads" = "gonad")) %>%
      mutate(group = paste(sex, tissue, treatment, sep = "."),
             study = ifelse(grepl("m.|extend|prolong", treatment), 
                            "manipulation", "charcterization")) %>%
      mutate(study = factor(study, levels = c("charcterization", "manipulation")))
    colData <- as.data.frame(colData)
    row.names(colData) <- colData$V1
    str(colData)

    ## 'data.frame':    987 obs. of  7 variables:
    ##  $ V1       : chr  "L.Blu13_male_gonad_control.NYNO" "L.Blu13_male_hypothalamus_control.NYNO" "L.Blu13_male_pituitary_control.NYNO" "L.G107_male_gonad_control" ...
    ##  $ bird     : chr  "L.Blu13" "L.Blu13" "L.Blu13" "L.G107" ...
    ##  $ sex      : Factor w/ 2 levels "female","male": 2 2 2 2 2 2 1 1 1 2 ...
    ##  $ tissue   : Factor w/ 3 levels "hypothalamus",..: 3 1 2 3 1 2 3 1 2 3 ...
    ##  $ treatment: Factor w/ 16 levels "bldg","control",..: 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ group    : chr  "male.gonads.control" "male.hypothalamus.control" "male.pituitary.control" "male.gonads.control" ...
    ##  $ study    : Factor w/ 2 levels "charcterization",..: 1 1 1 1 1 1 1 1 1 1 ...

    ## check that rownames and colnames match for DESeq
    ncol(countData) == nrow(colData)

    ## [1] TRUE

summarize data
--------------

    colData %>% select(sex, tissue, treatment, study)  %>%  summary()

    ##      sex               tissue        treatment               study    
    ##  female:497   hypothalamus:327   control  : 73   charcterization:636  
    ##  male  :490   pituitary   :330   inc.d9   : 71   manipulation   :351  
    ##               gonads      :330   inc.d17  : 66                        
    ##                                  n9       : 66                        
    ##                                  m.inc.d17: 63                        
    ##                                  bldg     : 60                        
    ##                                  (Other)  :588

    table(colData$sex, colData$treatment, colData$tissue)

    ## , ,  = hypothalamus
    ## 
    ##         
    ##          bldg control early extend hatch inc.d17 inc.d3 inc.d9 lay
    ##   female   10      11    10     10    10      11     10     12  10
    ##   male     10      11    10     10    10      11     10     11  10
    ##         
    ##          m.inc.d17 m.inc.d3 m.inc.d9 m.n2 n5 n9 prolong
    ##   female        11       10        9   10 10 11      10
    ##   male          10       10        8   10 10 11      10
    ## 
    ## , ,  = pituitary
    ## 
    ##         
    ##          bldg control early extend hatch inc.d17 inc.d3 inc.d9 lay
    ##   female   10      11    10     10    10      11     10     13  10
    ##   male     10      14    10     10    10      11     10     11  10
    ##         
    ##          m.inc.d17 m.inc.d3 m.inc.d9 m.n2 n5 n9 prolong
    ##   female        11       10        8   10 10 11      10
    ##   male          10       10        8   10 10 11      10
    ## 
    ## , ,  = gonads
    ## 
    ##         
    ##          bldg control early extend hatch inc.d17 inc.d3 inc.d9 lay
    ##   female   10      13    10     10    10      11     10     13  10
    ##   male     10      13    10     10    10      11     10     11  10
    ##         
    ##          m.inc.d17 m.inc.d3 m.inc.d9 m.n2 n5 n9 prolong
    ##   female        11       10        8   10 10 11      10
    ##   male          10       10        8    9 10 11      10

    length(unique(colData$bird))

    ## [1] 332

    table(colData$sex, colData$tissue)

    ##         
    ##          hypothalamus pituitary gonads
    ##   female          165       165    167
    ##   male            162       165    163

save bird data
--------------

    birds <- colData %>%
      select(bird, sex, treatment) %>%
      distinct()
    head(birds)

    ##      bird    sex treatment
    ## 1 L.Blu13   male   control
    ## 2  L.G107   male   control
    ## 3  L.G118 female   control
    ## 4    L.R3   male   control
    ## 5    L.R8   male   control
    ## 6   L.W33   male   control

save files for downstream use
-----------------------------

    write.csv(countData, "../results/00_counts.csv")
    write.csv(isoforms, "../results/00_geneswithisoforms.csv", row.names = T)

    write.csv(colData, "../metadata/00_colData.csv")
    write.csv(birds, "../metadata/00_birds.csv")
    write.csv(geneinfo, "../metadata/00_geneinfo.csv", row.names = TRUE)

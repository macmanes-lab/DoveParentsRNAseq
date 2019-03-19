Data Wrangling
==============

countData and geneinfo
----------------------

    # import count data, set rows at entreziz
    countData <- read.table("../results/kallistocounts.txt", 
                            sep = ",", row.names = NULL)

    ## set row names as entrezid for deseq2
    row.names(countData) <- countData$entrezid

    # replace dashes with periods (to make deseq happy)
    colnames(countData) <- gsub('-', '.', colnames(countData))

    # create df with gene info and remove this from countData
    geneinfo <- countData %>%
      select(row.names, Name, geneid, entrezid)
    head(geneinfo)

    ##                row.names     Name geneid       entrezid
    ## NP_001001127.1    408082    EDNRB 408082 NP_001001127.1
    ## NP_001001129.1    408183  CYP26A1 408183 NP_001001129.1
    ## NP_001001189.1    374073    CFDP1 374073 NP_001001189.1
    ## NP_001001194.1    407777    AvBD7 407777 NP_001001194.1
    ## NP_001001195.1    407779     KRT5 407779 NP_001001195.1
    ## NP_001001201.1    408034 HSD11B1L 408034 NP_001001201.1

    countData <- countData %>%
      select (-c(row.names, Name, geneid, entrezid))

    # because DESeq require integers
    countData <- round(countData)

    # print tolal num of genes and samples
    dim(countData)

    ## [1] 14937   987

    # make a tibble for the sole purpose of previewing this huge df
    countDataTibble <- as.tbl(countData)
    print(countDataTibble, n_extra = 0)

    ## # A tibble: 14,937 x 987
    ##    L.Blu13_male_go… L.Blu13_male_hy… L.Blu13_male_pi… L.G107_male_gon…
    ##               <dbl>            <dbl>            <dbl>            <dbl>
    ##  1                6               37               46                7
    ##  2                1                0               47                2
    ##  3              408               72              197              728
    ##  4                1                0                0                1
    ##  5                0                0                0                2
    ##  6               65                3               30              147
    ##  7               12                1                1                0
    ##  8               38               21                1               97
    ##  9                2              142               12                0
    ## 10                3                7                4                4
    ## # … with 14,927 more rows, and 983 more variables

wrangle colData
---------------

    colData <- read.table(file.path( "../metadata/kallistosamples.txt"), header = F, stringsAsFactors = F)

    # use strsplit to cut the filename into meaningful columns
    colData$bird <- sapply(strsplit(as.character(colData$V1),'\\_'), "[", 1)
    colData$sex <- sapply(strsplit(as.character(colData$V1),'\\_'), "[", 2)
    colData$tissue <- sapply(strsplit(as.character(colData$V1),'\\_'), "[", 3)

    # to create columns with variable names in the filename 
    colData$temp <- sapply(strsplit(as.character(colData$V1),'\\_'), "[", 4)
    colData$treatmentTemp <- sapply(strsplit(as.character(colData$temp),'\\.'), "[", 1)
    colData$NYNO <- sapply(strsplit(as.character(colData$temp),'\\.'), "[", 2)

    # relabel incorrectly named variables 
    unique(colData$treatmentTemp)

    ##  [1] "control"      "prolong"      "m-hatch"      "m-inc-d3"    
    ##  [5] "inc-d9"       "m-inc-d8"     "m-n2"         "bldg"        
    ##  [9] "n5"           "inc-d17"      "extend"       "hatch"       
    ## [13] "n9"           "inc-d3"       "inc-prolong"  "m-inc-d9"    
    ## [17] "extend-hatch" "m-inc-d17"    "lay"

    colData$treatment <- ifelse(grepl("extend-hatch", colData$treatmentTemp), "extend", colData$treatmentTemp) 

    # drop no-longer-needed cols
    colData <- colData %>% dplyr::select(-c(temp, treatmentTemp))

    # replace dashes with periods and make rowname, (to make deseq happy)
    colData$bird <- gsub("-", ".", colData$bird)
    colData$treatment <- gsub("-", ".", colData$treatment)
    colData$V1 <- gsub("-", ".", colData$V1)
    row.names(colData) <- colData$V1

    # make variables factors and specify level
    cols <- c("sex", "treatment", "tissue")
    colData[cols] <- lapply(colData[cols], factor) 
    colData$tissue <- factor(colData$tissue, 
                             levels = c("hypothalamus", "pituitary", "gonad"))

    # create a combinatorial variable
    colData$group <- paste(colData$sex, colData$tissue, colData$treatment, sep=".")

    ## view colData
    str(colData)

    ## 'data.frame':    987 obs. of  7 variables:
    ##  $ V1       : chr  "L.Blu13_male_gonad_control.NYNO" "L.Blu13_male_hypothalamus_control.NYNO" "L.Blu13_male_pituitary_control.NYNO" "L.G107_male_gonad_control" ...
    ##  $ bird     : chr  "L.Blu13" "L.Blu13" "L.Blu13" "L.G107" ...
    ##  $ sex      : Factor w/ 2 levels "female","male": 2 2 2 2 2 2 1 1 1 2 ...
    ##  $ tissue   : Factor w/ 3 levels "hypothalamus",..: 3 1 2 3 1 2 3 1 2 3 ...
    ##  $ NYNO     : chr  "NYNO" "NYNO" "NYNO" NA ...
    ##  $ treatment: Factor w/ 18 levels "bldg","control",..: 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ group    : chr  "male.gonad.control" "male.hypothalamus.control" "male.pituitary.control" "male.gonad.control" ...

    head(colData, 5)

    ##                                                                            V1
    ## L.Blu13_male_gonad_control.NYNO               L.Blu13_male_gonad_control.NYNO
    ## L.Blu13_male_hypothalamus_control.NYNO L.Blu13_male_hypothalamus_control.NYNO
    ## L.Blu13_male_pituitary_control.NYNO       L.Blu13_male_pituitary_control.NYNO
    ## L.G107_male_gonad_control                           L.G107_male_gonad_control
    ## L.G107_male_hypothalamus_control             L.G107_male_hypothalamus_control
    ##                                           bird  sex       tissue NYNO
    ## L.Blu13_male_gonad_control.NYNO        L.Blu13 male        gonad NYNO
    ## L.Blu13_male_hypothalamus_control.NYNO L.Blu13 male hypothalamus NYNO
    ## L.Blu13_male_pituitary_control.NYNO    L.Blu13 male    pituitary NYNO
    ## L.G107_male_gonad_control               L.G107 male        gonad <NA>
    ## L.G107_male_hypothalamus_control        L.G107 male hypothalamus <NA>
    ##                                        treatment                     group
    ## L.Blu13_male_gonad_control.NYNO          control        male.gonad.control
    ## L.Blu13_male_hypothalamus_control.NYNO   control male.hypothalamus.control
    ## L.Blu13_male_pituitary_control.NYNO      control    male.pituitary.control
    ## L.G107_male_gonad_control                control        male.gonad.control
    ## L.G107_male_hypothalamus_control         control male.hypothalamus.control

    # check that rownames and colnames match for DESeq
    ncol(countData) == nrow(colData)

    ## [1] TRUE

subset for manipulation study
-----------------------------

    # create a new dfs for just the manipulation study
    colData_manipluation <- colData %>%
      dplyr::filter(grepl('m.|extend|prolong', treatment)) %>%
      droplevels()
    row.names(colData_manipluation) <- colData_manipluation$V1

    # print sample sizes
    colData_manipluation %>% select(sex,treatment, tissue)  %>%  summary()

    ##      sex          treatment           tissue   
    ##  female:208   m.inc.d17:63   hypothalamus:138  
    ##  male  :203   extend   :60   pituitary   :137  
    ##               m.inc.d3 :60   gonad       :136  
    ##               m.inc.d8 :60                     
    ##               m.inc.d9 :49                     
    ##               m.hatch  :47                     
    ##               (Other)  :72

    # create list with filenames and select countData columns matching the file names
    savecols <- as.character(colData_manipluation$V1) 
    savecols <- as.vector(savecols) 
    countData_manipluation <- countData %>% dplyr::select(one_of(savecols)) 

    # confirm that row and col lenghts are equal
    ncol(countData_manipluation) == nrow(colData_manipluation)

    ## [1] TRUE

subset for characterization study
---------------------------------

    # filter out data related to the manipulation study
    colData_characterization <- colData %>%
      dplyr::filter(!grepl('m.|extend|prolong', treatment)) %>%
      droplevels()
    row.names(colData_characterization) <- colData_characterization$V1

    # print sample sizes
    colData_characterization %>% select(sex,treatment, tissue)  %>%  summary()

    ##      sex        treatment            tissue   
    ##  female:289   control: 73   hypothalamus:189  
    ##  male  :287   inc.d9 : 71   pituitary   :193  
    ##               inc.d17: 66   gonad       :194  
    ##               n9     : 66                     
    ##               bldg   : 60                     
    ##               hatch  : 60                     
    ##               (Other):180

    savecols <- as.character(colData_characterization$V1) 
    savecols <- as.vector(savecols) 
    countData_characterization <- countData %>% dplyr::select(one_of(savecols)) 

    # check that row and col lenghts are equal
    ncol(countData_characterization) == nrow(colData_characterization)

    ## [1] TRUE

save files for downstream use
-----------------------------

    write.csv(colData_manipluation, "../results/00_colData_manipluation.csv", row.names = TRUE) 
    write.csv(countData_manipluation, "../results/00_countData_manipluation.csv", row.names = TRUE) 
    write.csv(colData_characterization, "../results/00_colData_characterization.csv", row.names = TRUE) 
    write.csv(countData_characterization, "../results/00_countData_characterization.csv", row.names = TRUE) 
    write.csv(geneinfo, "../results/00_geneinfo.csv", row.names = TRUE)

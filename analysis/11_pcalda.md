pca analysis
------------

    pseudocounts <- read_csv("../results/01_pseudo.counts.csv")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )

    ## See spec(...) for full column specifications.

    colData <- read.csv("../metadata/00_colData_characterization.csv", header = T, row.names = 1)
    colData$treatment <- factor(colData$treatment, levels = 
                                  c("control",  "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", 
                                    "hatch", "n5", "n9"))

    geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)

    subsetcounts <- function(colData, countData){
      
      # save counts that match colData
      savecols <- as.character(colData$V1) 
      savecols <- as.vector(savecols) 
      countData <- countData %>% 
        dplyr::select(X1,savecols) %>% 
        dplyr::rename(entrezid = X1)  
      
      countData <- left_join(geneinfo, countData)
      countData <- countData %>% distinct(Name, .keep_all = TRUE) %>%
        select(-row.names)  %>% select(-geneid) %>% select(-entrezid) 
      
      
      row.names(countData) <- countData$Name
      countData$Name <- NULL
      countData <- as.data.frame(t(countData))
      
      return(countData)
    }

    countData <- subsetcounts(colData, pseudocounts)

    ## Joining, by = "entrezid"

    ## Warning: Column `entrezid` joining factor and character vector, coercing
    ## into character vector

    # confirm
    row.names(countData) == row.names(colData)

    ##   [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [15] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [29] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [43] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [57] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [71] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [85] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [99] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [113] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [127] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [141] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [155] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [169] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [183] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [197] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [211] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [225] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [239] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [253] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [267] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [281] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [295] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [309] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [323] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [337] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [351] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [365] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [379] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [393] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [407] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [421] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [435] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [449] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [463] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [477] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [491] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [505] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [519] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [533] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [547] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [561] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [575] TRUE TRUE

    # make pca
    mypca <- prcomp(countData)

    # percent conribution
    eig.val <- get_eigenvalue(mypca)
    head(eig.val)

    ##         eigenvalue variance.percent cumulative.variance.percent
    ## Dim.1 301234419601       35.9445757                    35.94458
    ## Dim.2 272266589259       32.4880106                    68.43259
    ## Dim.3 213391804875       25.4628203                    93.89541
    ## Dim.4  24543770401        2.9286674                    96.82407
    ## Dim.5  10627670724        1.2681390                    98.09221
    ## Dim.6   3598570939        0.4293968                    98.52161

    mypcadf <- data.frame(PC1 = mypca$x[, 1], PC2 = mypca$x[, 2], PC3 = mypca$x[, 3], 
                      PC4 = mypca$x[, 4],PC5 = mypca$x[, 5],PC6 = mypca$x[, 6],
                      ID = row.names(countData))
    head(mypcadf)

    ##                                              PC1        PC2        PC3
    ## L.Blu13_male_gonad_control.NYNO        -289031.0    9583.47  418914.02
    ## L.Blu13_male_hypothalamus_control.NYNO -885455.5 -657013.51 -557113.05
    ## L.Blu13_male_pituitary_control.NYNO    -294197.4  150695.51  211668.84
    ## L.G107_male_gonad_control              -225279.4  -37945.28  363813.73
    ## L.G107_male_hypothalamus_control       -484957.1 -344082.98  -89470.83
    ## L.G107_male_pituitary_control          -352850.2  276251.11  129694.05
    ##                                              PC4       PC5          PC6
    ## L.Blu13_male_gonad_control.NYNO        -96018.57 -33516.68   55973.7678
    ## L.Blu13_male_hypothalamus_control.NYNO -56372.39  -8591.00 -103116.0847
    ## L.Blu13_male_pituitary_control.NYNO    416519.89  14124.25    9801.6432
    ## L.G107_male_gonad_control              -94072.20 -31764.49   64140.8075
    ## L.G107_male_hypothalamus_control       -65247.42 -16252.64    4662.0743
    ## L.G107_male_pituitary_control          403886.93  13368.66    -237.8222
    ##                                                                            ID
    ## L.Blu13_male_gonad_control.NYNO               L.Blu13_male_gonad_control.NYNO
    ## L.Blu13_male_hypothalamus_control.NYNO L.Blu13_male_hypothalamus_control.NYNO
    ## L.Blu13_male_pituitary_control.NYNO       L.Blu13_male_pituitary_control.NYNO
    ## L.G107_male_gonad_control                           L.G107_male_gonad_control
    ## L.G107_male_hypothalamus_control             L.G107_male_hypothalamus_control
    ## L.G107_male_pituitary_control                   L.G107_male_pituitary_control

    mypcadf$V1 <- row.names(mypcadf)
    mypcadf <- left_join(colData, mypcadf)

    ## Joining, by = "V1"

    ## Warning: Column `V1` joining factor and character vector, coercing into
    ## character vector

    head(mypcadf)

    ##                                       V1    bird  sex       tissue
    ## 1        L.Blu13_male_gonad_control.NYNO L.Blu13 male        gonad
    ## 2 L.Blu13_male_hypothalamus_control.NYNO L.Blu13 male hypothalamus
    ## 3    L.Blu13_male_pituitary_control.NYNO L.Blu13 male    pituitary
    ## 4              L.G107_male_gonad_control  L.G107 male        gonad
    ## 5       L.G107_male_hypothalamus_control  L.G107 male hypothalamus
    ## 6          L.G107_male_pituitary_control  L.G107 male    pituitary
    ##   treatment                     group           study       PC1        PC2
    ## 1   control        male.gonad.control charcterization -289031.0    9583.47
    ## 2   control male.hypothalamus.control charcterization -885455.5 -657013.51
    ## 3   control    male.pituitary.control charcterization -294197.4  150695.51
    ## 4   control        male.gonad.control charcterization -225279.4  -37945.28
    ## 5   control male.hypothalamus.control charcterization -484957.1 -344082.98
    ## 6   control    male.pituitary.control charcterization -352850.2  276251.11
    ##          PC3       PC4       PC5          PC6
    ## 1  418914.02 -96018.57 -33516.68   55973.7678
    ## 2 -557113.05 -56372.39  -8591.00 -103116.0847
    ## 3  211668.84 416519.89  14124.25    9801.6432
    ## 4  363813.73 -94072.20 -31764.49   64140.8075
    ## 5  -89470.83 -65247.42 -16252.64    4662.0743
    ## 6  129694.05 403886.93  13368.66    -237.8222
    ##                                       ID
    ## 1        L.Blu13_male_gonad_control.NYNO
    ## 2 L.Blu13_male_hypothalamus_control.NYNO
    ## 3    L.Blu13_male_pituitary_control.NYNO
    ## 4              L.G107_male_gonad_control
    ## 5       L.G107_male_hypothalamus_control
    ## 6          L.G107_male_pituitary_control

    a <- ggplot(mypcadf, aes(x = PC1, y = PC2, color = colData$treatment)) +
      geom_point() + labs(x = "PC1: 36%", y = "PC2: 32%")
    b <- ggplot(mypcadf, aes(x = PC3, y = PC4, color = colData$treatment)) +
      geom_point() + theme(legend.position = "bottom", legend.title = element_blank()) + labs(x = "PC3: 25%", y = "PC4: 3%")

    mylegend <- get_legend(b)

    ab <- plot_grid(a + theme(legend.position = "none") ,
                    b + theme(legend.position = "none"))

    plot_grid(ab,mylegend, ncol = 1, rel_heights = c(0.9,0.1))

![](../figures/11_PCA.LDA/pca-1.png)

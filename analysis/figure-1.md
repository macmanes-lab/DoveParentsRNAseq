Figure 1. Experimental design and tSNE analysis
===============================================

wrange data
-----------

Note, the input file is too big for storage on github :(

    pseudocounts <- read_csv("../results/01_pseudo.counts.csv")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )

    ## See spec(...) for full column specifications.

    head(pseudocounts[1:3])

    ## # A tibble: 6 x 3
    ##   X1     L.Blu13_male_gonad_control.NY… L.Blu13_male_hypothalamus_control.…
    ##   <chr>                           <dbl>                               <dbl>
    ## 1 A2ML1                          42.7                               201.   
    ## 2 A2ML2                           4.44                                4.86 
    ## 3 A2ML3                         212.                               6920.   
    ## 4 A2ML4                           6.91                                0.810
    ## 5 A4GALT                         14.3                                 8.53 
    ## 6 A4GNT                           0.206                               0.105

    pseudocounts <- as.data.frame(pseudocounts)
    row.names(pseudocounts) <- pseudocounts$X1
    pseudocounts$X1 <- NULL
    # prep count data for all samples
    countData <- as.data.frame(t(pseudocounts))
    head(countData[1:3])

    ##                                            A2ML1     A2ML2      A2ML3
    ## L.Blu13_male_gonad_control.NYNO         42.65677 4.4397552  211.96191
    ## L.Blu13_male_hypothalamus_control.NYNO 201.26331 4.8633048 6919.87335
    ## L.Blu13_male_pituitary_control.NYNO    161.22614 0.3158851  212.10845
    ## L.G107_male_gonad_control               43.22441 2.0404458  203.74048
    ## L.G107_male_hypothalamus_control       382.33084 4.9190817 9531.52550
    ## L.G107_male_pituitary_control           85.34910 0.3577761   69.02124

    # prep col data for all samples
    colData <- read.csv("../metadata/00_samples.csv", header = T, row.names = 1)
    colData$treatment <- factor(colData$treatment, levels = alllevels)
    colData <- colData %>% mutate(tissue = fct_recode(tissue, "gonads" = "gonad"))
    colData$tissue <- factor(colData$tissue, levels = tissuelevels)
    row.names(colData) <- colData$V1

    # check ready for analysis
    #row.names(countData) == row.names(colData)

    head(colData)

    ##                                                                            V1
    ## L.Blu13_male_gonad_control.NYNO               L.Blu13_male_gonad_control.NYNO
    ## L.Blu13_male_hypothalamus_control.NYNO L.Blu13_male_hypothalamus_control.NYNO
    ## L.Blu13_male_pituitary_control.NYNO       L.Blu13_male_pituitary_control.NYNO
    ## L.G107_male_gonad_control                           L.G107_male_gonad_control
    ## L.G107_male_hypothalamus_control             L.G107_male_hypothalamus_control
    ## L.G107_male_pituitary_control                   L.G107_male_pituitary_control
    ##                                           bird  sex       tissue treatment
    ## L.Blu13_male_gonad_control.NYNO        L.Blu13 male       gonads   control
    ## L.Blu13_male_hypothalamus_control.NYNO L.Blu13 male hypothalamus   control
    ## L.Blu13_male_pituitary_control.NYNO    L.Blu13 male    pituitary   control
    ## L.G107_male_gonad_control               L.G107 male       gonads   control
    ## L.G107_male_hypothalamus_control        L.G107 male hypothalamus   control
    ## L.G107_male_pituitary_control           L.G107 male    pituitary   control
    ##                                                            group
    ## L.Blu13_male_gonad_control.NYNO               male.gonad.control
    ## L.Blu13_male_hypothalamus_control.NYNO male.hypothalamus.control
    ## L.Blu13_male_pituitary_control.NYNO       male.pituitary.control
    ## L.G107_male_gonad_control                     male.gonad.control
    ## L.G107_male_hypothalamus_control       male.hypothalamus.control
    ## L.G107_male_pituitary_control             male.pituitary.control
    ##                                                  study
    ## L.Blu13_male_gonad_control.NYNO        charcterization
    ## L.Blu13_male_hypothalamus_control.NYNO charcterization
    ## L.Blu13_male_pituitary_control.NYNO    charcterization
    ## L.G107_male_gonad_control              charcterization
    ## L.G107_male_hypothalamus_control       charcterization
    ## L.G107_male_pituitary_control          charcterization

    # prep for tsne

    subsetmaketsne <- function(whichtissue, whichtreatment, whichsex){

      colData <- colData %>%
        dplyr::filter(tissue %in% whichtissue,
                      treatment %in% whichtreatment,
                      sex %in% whichsex) 
      row.names(colData) <- colData$V1

      # save counts that match colData
      savecols <- as.character(colData$V1) 
      savecols <- as.vector(savecols) 
      
      countData <- as.data.frame(t(countData))
      countData <- countData %>% dplyr::select(one_of(savecols)) 
      countData <- as.data.frame(t(countData))

      euclidist <- dist(countData) # euclidean distances between the rows

      tsne_model <- Rtsne(euclidist, check_duplicates=FALSE, pca=TRUE, perplexity=10, theta=0.5, dims=2)
      tsne_df = as.data.frame(tsne_model$Y) 

      # prep for adding columns
      colData2 <- colData 
      colData2$V1 <- NULL
      tsne_df_cols <- cbind(colData2, tsne_df)
      return(tsne_df_cols)
    }

    chartsne <- subsetmaketsne(tissuelevels, charlevels, sexlevels)
    hyptsne <- subsetmaketsne("hypothalamus", charlevels, sexlevels)
    pittsne <- subsetmaketsne("pituitary", charlevels, sexlevels)
    gontsne <- subsetmaketsne("gonads", charlevels, sexlevels)

    ftsne <-  subsetmaketsne(tissuelevels, charlevels, "female")
    mtsne <-  subsetmaketsne(tissuelevels, charlevels, "male")

make figure
-----------

    plottsneelipse <- function(tsnedf, whichfactor, whichcolors){
      p <- ggplot(tsnedf, aes(x = V1, y = V2)) +
        geom_point(size = 1, aes(color = whichfactor)) +
        theme_B3() +
        labs(x = "tSNE 1", y = "tSNE 2") +
        scale_color_manual(values = whichcolors) +
        theme(legend.position = "none",
              axis.text = element_blank()) +
        stat_ellipse(linetype = 2, aes(color = tissue )) 
      return(p)
    }


    a <- plottsneelipse(chartsne, chartsne$tissue, allcolors) + labs(subtitle = "tissue\n")    
    b <- plottsneelipse(chartsne, chartsne$sex, allcolors)   + labs(y = " ", subtitle = "tissue * sex\n")    
    c <- plottsneelipse(ftsne, ftsne$treatment, allcolors ) + labs(y = " ", subtitle = "females\ntissue * treatment")
    d <- plottsneelipse(mtsne, mtsne$treatment, allcolors ) + labs(y = " ", subtitle = "males\ntissue * treatment") 

    abcd <- plot_grid(a,b,c,d, nrow = 1, labels = c("b", "c", "d", "e"), label_size = 12 )

    expdesign <- png::readPNG("../figures/images/fig1a_fig1a.png")
    expdesign <- ggdraw() +  draw_image(expdesign, scale = 1)

    fig1 <- plot_grid(expdesign, abcd, nrow = 2, labels = c("a", "b"), label_size = 12, rel_heights = c(0.5,1))
    fig1

![](../figures/fig1-1.png)

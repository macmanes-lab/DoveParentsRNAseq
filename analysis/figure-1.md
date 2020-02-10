Figure 1. Experimental design, tSNE analysis, PCA, and prolactin
================================================================

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

    plottsneelipse <- function(tsnedf, pointcolor, whichcolors){
      p <- ggplot(tsnedf, aes(x = V1, y = V2)) +
        geom_point(size = 1, aes(color = pointcolor)) +
        theme_B3() +
        labs(x = "tSNE 1", y = "tSNE 2") +
        scale_color_manual(values = whichcolors) +
        theme(legend.position = "none",
              axis.text = element_blank()) +
        stat_ellipse(linetype = 1, aes(color = tissue )) 
      return(p)
    }


    a <- plottsneelipse(chartsne, chartsne$tissue, allcolors) + labs(subtitle = "tissue\n")    
    b <- plottsneelipse(chartsne, chartsne$sex, allcolors)   + labs(y = " ", subtitle = "tissue * sex\n")    
    c <- plottsneelipse(ftsne, ftsne$treatment, allcolors ) + labs(y = " ", subtitle = "treatment * tissue\nfemales")
    d <- plottsneelipse(mtsne, mtsne$treatment, allcolors ) + labs(y = " ", subtitle = "treatment * tissue\nmales") 

    abcd <- plot_grid(a,b,c,d, nrow = 1, labels = c("b", "c", "d", "e"), label_size = 12 )

    expdesign <- png::readPNG("../figures/images/fig1a_fig1a.png")
    expdesign <- ggdraw() +  draw_image(expdesign, scale = 1)

    abcde <- plot_grid(expdesign, abcd, nrow = 2, labels = c("a", "b"), label_size = 12, rel_heights = c(0.5,1))



    plottsneelipsev2 <- function(tsnedf, pointcolor, whichcolors){
      p <- ggplot(tsnedf, aes(x = V1, y = V2)) +
        geom_point(size = 1, aes(color = pointcolor)) +
        theme_B3() +
        labs(x = "tSNE 1", y = "tSNE 2") +
        scale_color_manual(values = whichcolors) +
        theme(legend.position = "none",
              axis.text = element_blank()) +
        stat_ellipse(linetype = 1, aes(color = treatment)) 
      return(p)
    }

    h <- plottsneelipsev2(hyptsne, hyptsne$treatment, allcolors) + labs(subtitle = "hypothalamus")  + facet_wrap(~sex)
    i <- plottsneelipsev2(pittsne, pittsne$treatment, allcolors ) + labs(subtitle = "pituitary", y = NULL) + facet_wrap(~sex) 
    j <- plottsneelipsev2(gontsne, gontsne$treatment, allcolors ) + labs(subtitle = "gonads", y = NULL)  + facet_wrap(~sex)

    hij <- plot_grid(h,i,j, nrow = 1, labels = c("h", "i", "j"), label_size = 12)

    fig1 <- plot_grid(abcde, hij, nrow = 2, rel_heights = c(0.6,0.4))
    fig1

![](../figures/fig1-1.png)

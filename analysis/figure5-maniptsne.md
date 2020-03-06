Figure 1. Experimental design, tSNE analysis, PCA, and prolactin
================================================================

    # prep for tsne

    hyptsne <- subsetmaketsne("hypothalamus", allmaniplevels, sexlevels)
    pittsne <- subsetmaketsne("pituitary", allmaniplevels, sexlevels)
    gontsne <- subsetmaketsne("gonads", allmaniplevels, sexlevels)

    addgroupings <- function(df){
      
      df <- df %>% mutate(hiloPRL = fct_collapse(treatment, 
                                      lo = c("m.inc.d8", "m.inc.d3", "m.inc.d9", "inc.d3", "inc.d9", "lay"),
                                      hi = c("inc.d17",   "m.inc.d17", "prolong" ,  "hatch"  ,   "m.n2", "extend", "n5")),
                    extint = fct_collapse(treatment, 
                                      eggs = c("inc.d3", "inc.d9", "inc.d17", "prolong", "lay"),
                                      chicks = c( "m.inc.d8", "hatch", "extend", "n5"),
                                      loss = c("m.inc.d3", "m.inc.d9", "m.n2", "m.inc.d17")))
      df$extint <- factor(df$extint, levels = levelsextint)
      return(df)
    }
      

    hyptsne <- addgroupings(hyptsne)
    pittsne <- addgroupings(pittsne)
    gontsne <- addgroupings(gontsne)

make figure
-----------

    expdesign2 <- png::readPNG("../figures/images/DoveParentsRNAseq_hypothesis.png")
    expdesign2 <- ggdraw() +  draw_image(expdesign2, scale = 1)


    b <- plottsneelipsev2(hyptsne, hyptsne$treatment, allcolors) + labs(subtitle = "hypothalamus ~ parental stage", x = NULL)  + facet_wrap(~sex, scales = "free")
    c <- plottsneelipsev2(pittsne, pittsne$treatment, allcolors ) + labs(subtitle = "pituitary ~ parental stage", x = NULL) + facet_wrap(~sex, scales = "free") 
    d <- plottsneelipsev2(gontsne, gontsne$treatment, allcolors ) + labs(subtitle = "gonads ~ parental stage")  + facet_wrap(~sex, scales = "free")  

    bcd <- plot_grid(b,c,d, nrow = 3, labels = c("b", "c", "d"), label_size = 8, rel_heights = c(1,1,1))

    ## Warning in MASS::cov.trob(data[, vars]): Probable convergence failure

    e <- plottsneelipsev2(hyptsne, hyptsne$hiloPRL, allcolors ) + labs(subtitle = "hypothalamus ~ estimated prolactin", x = NULL) + facet_wrap(~sex, scales = "free")  
    f <- plottsneelipsev2(pittsne, pittsne$hiloPRL, allcolors ) + labs(subtitle = "pituitary ~ estimated prolactin", x = NULL) + facet_wrap(~sex, scales = "free")   
    g <- plottsneelipsev2(gontsne, gontsne$hiloPRL, allcolors ) + labs(subtitle = "gonads ~ estimated prolactin" ) + facet_wrap(~sex, scales = "free")   

    efg <- plot_grid(e,f,g, nrow = 3, labels = c("h", "i", "j"), label_size = 8, rel_heights = c(1,1,1))

    h <- plottsneelipsev2(hyptsne, hyptsne$extint, allcolors ) + labs(subtitle = "hypothalamus ~ external enviornment", x = NULL) + facet_wrap(~sex, scales = "free")  
    i <- plottsneelipsev2(pittsne, pittsne$extint, allcolors ) + labs(subtitle = "pituitary ~ external enviornment", x = NULL) + facet_wrap(~sex, scales = "free")   
    j <- plottsneelipsev2(gontsne, gontsne$extint, allcolors ) + labs(subtitle = "gonads ~ external enviornment" ) + facet_wrap(~sex, scales = "free")   

    hij <- plot_grid(h,i,j, nrow = 3, labels = c("e", "f", "g"), label_size = 8, rel_heights = c(1,1,1))

    bcdhijefg <- plot_grid(bcd, hij,  efg, nrow = 1)


    fig6 <- plot_grid(expdesign2, bcdhijefg, nrow = 2, labels = c("e"), label_size = 8, rel_heights = c(0.4,1))
    fig6

![](../figures/fig5-1.png)

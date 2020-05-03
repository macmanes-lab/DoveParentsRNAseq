Music example
-------------

    song <- c("Oh,", "sin-", "ner", "man,", "where", "_", "he", "gon'", "run", "_", "to?")
    stage <- c( "bldg", "lay", "inc.d3", "m.inc.d3" ,
                "inc.d9", "m.inc.d9"  , 
                "inc.d17", "m.inc.d17" ,
                "hatch", "m.n2" ,
                "n5")
    position <- c(1,2,3,4,5,6,7,8,9,10,11)
    gene1 <- c(4,4,4,4,6,6,5,4,4,4,4)
    gene2 <- gene1 + 2
    gene3 <- gene2 + 2

    sinnerman <- as.data.frame(t(rbind(position, gene1,gene2,gene3)))
    sinnerman$song <- song
    sinnerman$stage <- stage
    sinnerman$image <- sample(c("https://cdn.pixabay.com/photo/2016/03/23/20/49/music-note-1275650_960_720.png"), 
                            size=11, replace = TRUE)
    sinnerman <- sinnerman %>%
      pivot_longer(gene1:gene3, names_to = "music", values_to = "genes")  

    sinnermansong <- ggplot(sinnerman, aes(x = position, y = genes)) +
      #geom_line(aes(color = music)) +
      geom_image(aes(image=image), size = 0.15) +
      scale_x_continuous(breaks = position,
                        labels = song,
                         name = element_blank()) +
      scale_y_continuous(breaks = c(3, 5, 7, 9),
                         labels = c( "F",   "A",   "C",  "E"),
                         limits = c(-1,12)) +
      geom_hline(yintercept=c(2.15,4.15,6.15,8.15,10.15))   +
      theme_music() +
      theme(axis.ticks = element_blank()) +
      labs(subtitle = "Sinnerman by Nina Simone", y = NULL)
    sinnermansong

![](../figures/favegenes/transcriptionalsymphony-1.png)

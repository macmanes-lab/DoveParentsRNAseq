R Markdown
----------

    song <- c("What", "you", "want,", "ba-", "by", "I", "got.")
    stage <- c( "bldg", "lay", "inc.d9", "inc.d17", "hatch", "n5", "n9")
    time <- c(1,2,3,4,5,6,7)
    position <- c(1,2,3,4,5,6,7)
    gene1 <- c(11,9,7,8,8,8,6)
    gene2 <- c(7, 9,5,6,8,8,4)
    gene3 <- c(4, 9,2,4,8,8,1)

    respect <- as.data.frame(t(rbind(position, gene1,gene2,gene3)))
    respect$song <- song
    respect$stage <- stage
    respect$image <- sample(c("https://cdn.pixabay.com/photo/2016/03/23/20/49/music-note-1275650_960_720.png"), 
                            size=7, replace = TRUE)
    respect <- respect %>%
      pivot_longer(gene1:gene3, names_to = "music", values_to = "genes")  %>% 
      mutate(songstage = paste(song,stage, sep = "\n\n"))

    songstageslim <- unique(respect$songstage)

    respectsong <- ggplot(respect, aes(x = position, y = genes)) +
      geom_line(aes(color = music)) +
      geom_image(aes(image=image), size = 0.15) +
      scale_x_continuous(breaks = position,
                         labels = songstageslim,
                         name = element_blank()) +
      scale_y_continuous(breaks = c(3, 5, 7, 9),
                         labels = c( "F",   "A",   "C",  "E"),
                         limits = c(-1,12)) +
      geom_hline(yintercept=c(2.15,4.15,6.15,8.15,10.15))   +
      theme_B3() +
      theme(axis.ticks = element_blank()) +
      labs(subtitle = "RESPECT by Otis Redding, Jr.", y = NULL)
    respectsong

![](../figures/favegenes/respect-1.png)

    song <- c("Birds", "fly-", "ing", "high,", "you", "know", "how", "I", "feel.")
    stage <- charlevels
    position <- c(1,2,3,4,5,6,7,8,9)
    gene1 <- c(-2,0,1,2,2,1,1,-2,0)
    gene2 <- gene1 - 1
    gene1 <- gene2
    gene3 <- gene2

    feelinggood <- as.data.frame(t(rbind(position, gene1,gene2,gene3)))
    feelinggood$song <- song
    feelinggood$stage <- stage
    feelinggood$image <- sample(c("https://cdn.pixabay.com/photo/2016/03/23/20/49/music-note-1275650_960_720.png"), 
                            size=9, replace = TRUE)
    feelinggood <- feelinggood %>%
      pivot_longer(gene1:gene3, names_to = "music", values_to = "genes")   %>%
      mutate(songstage = paste(song,stage, sep = "\n\n"))
    feelinggood

    ## # A tibble: 27 x 7
    ##    position song  stage  image                      music genes songstage  
    ##       <dbl> <chr> <chr>  <chr>                      <chr> <dbl> <chr>      
    ##  1        1 Birds contr… https://cdn.pixabay.com/p… gene1    -3 "Birds\n\n…
    ##  2        1 Birds contr… https://cdn.pixabay.com/p… gene2    -3 "Birds\n\n…
    ##  3        1 Birds contr… https://cdn.pixabay.com/p… gene3    -3 "Birds\n\n…
    ##  4        2 fly-  bldg   https://cdn.pixabay.com/p… gene1    -1 "fly-\n\nb…
    ##  5        2 fly-  bldg   https://cdn.pixabay.com/p… gene2    -1 "fly-\n\nb…
    ##  6        2 fly-  bldg   https://cdn.pixabay.com/p… gene3    -1 "fly-\n\nb…
    ##  7        3 ing   lay    https://cdn.pixabay.com/p… gene1     0 "ing\n\nla…
    ##  8        3 ing   lay    https://cdn.pixabay.com/p… gene2     0 "ing\n\nla…
    ##  9        3 ing   lay    https://cdn.pixabay.com/p… gene3     0 "ing\n\nla…
    ## 10        4 high, inc.d3 https://cdn.pixabay.com/p… gene1     1 "high,\n\n…
    ## # … with 17 more rows

    songstageslim <- unique(feelinggood$songstage)

    feelinggoodsong <- ggplot(feelinggood, aes(x = position, y = genes)) +
      geom_line(aes(color = music)) +
      geom_image(aes(image=image), size = 0.15) +
      scale_x_continuous(breaks = position,
                        labels = songstageslim,
                         name = element_blank()) +
      scale_y_continuous(breaks = c(3, 5, 7, 9),
                         labels = c( "F",   "A",   "C",  "E"),
                         limits = c(-6,12)) +
      geom_hline(yintercept=c(2.15,4.15,6.15,8.15,10.15))   +
      theme_B3() +
      theme(axis.ticks = element_blank()) +
      labs(subtitle = "Feeling Good by L. Bricusse & A. Newley", y = NULL)
    feelinggoodsong

![](../figures/favegenes/feelinggood-1.png)

    song <- c("Oh,", "sin-", "ner", "man,,", "where", "_", "he", "gon'", "run", "_", "to?")
    stage <- c( "bldg", "lay", "inc.d3", "m.inc.d3" ,
                "inc.d9", "m.inc.d9"  , 
                "inc.d17", "m.inc.d17" ,
                "hatch", "m.n2" ,
                "n5")
    position <- c(1,2,3,4,5,6,7,8,9,10,11)
    gene1 <- c(0,0,0,0,2,2,1,0,0,0,0)
    gene2 <- gene1 + 2
    gene3 <- gene2 + 2

    sinnerman <- as.data.frame(t(rbind(position, gene1,gene2,gene3)))
    sinnerman$song <- song
    sinnerman$stage <- stage
    sinnerman$image <- sample(c("https://cdn.pixabay.com/photo/2016/03/23/20/49/music-note-1275650_960_720.png"), 
                            size=11, replace = TRUE)
    sinnerman <- sinnerman %>%
      pivot_longer(gene1:gene3, names_to = "music", values_to = "genes")   %>%
      mutate(songstage = paste(song,stage, sep = "\n\n"))
    sinnerman

    ## # A tibble: 33 x 7
    ##    position song  stage  image                     music genes songstage   
    ##       <dbl> <chr> <chr>  <chr>                     <chr> <dbl> <chr>       
    ##  1        1 Oh,   bldg   https://cdn.pixabay.com/… gene1     0 "Oh,\n\nbld…
    ##  2        1 Oh,   bldg   https://cdn.pixabay.com/… gene2     2 "Oh,\n\nbld…
    ##  3        1 Oh,   bldg   https://cdn.pixabay.com/… gene3     4 "Oh,\n\nbld…
    ##  4        2 sin-  lay    https://cdn.pixabay.com/… gene1     0 "sin-\n\nla…
    ##  5        2 sin-  lay    https://cdn.pixabay.com/… gene2     2 "sin-\n\nla…
    ##  6        2 sin-  lay    https://cdn.pixabay.com/… gene3     4 "sin-\n\nla…
    ##  7        3 ner   inc.d3 https://cdn.pixabay.com/… gene1     0 "ner\n\ninc…
    ##  8        3 ner   inc.d3 https://cdn.pixabay.com/… gene2     2 "ner\n\ninc…
    ##  9        3 ner   inc.d3 https://cdn.pixabay.com/… gene3     4 "ner\n\ninc…
    ## 10        4 man,, m.inc… https://cdn.pixabay.com/… gene1     0 "man,,\n\nm…
    ## # … with 23 more rows

    songstageslim <- unique(sinnerman$songstage)

    sinnermansong <- ggplot(sinnerman, aes(x = position, y = genes)) +
      geom_line(aes(color = music)) +
      geom_image(aes(image=image), size = 0.1) +
      scale_x_continuous(breaks = position,
                        labels = songstageslim,
                         name = element_blank()) +
      scale_y_continuous(breaks = c(3, 5, 7, 9),
                         labels = c( "F",   "A",   "C",  "E"),
                         limits = c(-1,12)) +
      geom_hline(yintercept=c(2.15,4.15,6.15,8.15,10.15))   +
      theme_B3() +
      theme(axis.ticks = element_blank()) +
      labs(subtitle = "Sinnerman by Nina Simone", y = NULL)
    sinnermansong

![](../figures/favegenes/sinnerman-1.png)

    p1 <- plot_grid(respectsong + theme(legend.position = "none"), feelinggoodsong + theme(legend.position = "none"), rel_widths = c(0.45,0.55))

    p3 <- plot_grid(p1, sinnermansong, nrow = 2)
    p3

![](../figures/favegenes/soulmusic-1.png)

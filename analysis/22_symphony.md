R Markdown
----------

    song <- c("what", "you", "want", "ba-", "by", "I", "got.")
    position <- c(1,2,3,4,5,6,7)
    note1 <- c(11,9,7,8,8,8,6)
    note2 <- c(7, 9,5,6,8,8,4)
    note3 <- c(4, 9,2,4,8,8,1)

    respect <- as.data.frame(t(rbind(position, note1,note2,note3)))
    respect$song <- song
    respect$image <- sample(c("https://cdn.pixabay.com/photo/2016/03/23/20/49/music-note-1275650_960_720.png"), 
                            size=7, replace = TRUE)
    respect <- respect %>%
      pivot_longer(note1:note3, names_to = "music", values_to = "notes")  %>% 
      mutate(position = as.numeric(position),
             notes = as.numeric(notes))

    song <- ggplot(respect, aes(x = position, y = notes)) +
      geom_line(aes(color = music)) +
      geom_image(aes(image=image), size = 0.15) +
      scale_x_continuous(breaks = position,
                         labels = song,
                         name = element_blank()) +
      scale_y_continuous(breaks = c(3, 5, 7, 9),
                         labels = c( "F",   "A",   "C",  "E"),
                         limits = c(-1,12)) +
      geom_hline(yintercept=c(2.15,4.15,6.15,8.15,10.15))   +
      theme_B3() +
      theme(axis.ticks = element_blank())

    stage <- c( "bldg", "lay", "inc.d9", "inc.d17", "hatch", "n5", "n9")
    time <- c(1,2,3,4,5,6,7)
    gene1 <- c(11,9,7,8,8,8,6)
    gene2 <- c(7, 9,5,6,8,8,4)
    gene3 <- c(4, 9,2,4,8,8,1)

    hypothesisgenes <- as.data.frame(t(rbind(time, gene1,gene2,gene3)))
    hypothesisgenes$stage <- stage
    hypothesisgenes$image <- sample(c("https://cdn.pixabay.com/photo/2016/03/23/20/49/music-note-1275650_960_720.png"), 
                            size=7, replace = TRUE)
    hypothesisgenes <- hypothesisgenes %>%
      pivot_longer(gene1:gene3, names_to = "data", values_to = "expression")  %>% 
      mutate( time = as.numeric(time))

    theory <- ggplot(hypothesisgenes, aes(x = time, y = expression)) +
      geom_line(aes(color = data)) +
      geom_image(aes(image=image), size = 0.15) +
      scale_x_continuous(breaks = time,
                         labels = stage,
                         name = element_blank()) +
      scale_y_continuous(breaks = c(2, 4, 6, 8, 10),
                         limits = c(-1,12)) +
      geom_hline(yintercept=c(2.15,4.15,6.15,8.15,10.15)) +
      theme_B3() +
      theme(axis.ticks = element_blank())

    plot_grid(song , theory, nrow = 2)

![](../figures/favegenes/respect-1.png)

    hypothesisgenes2 <- hypothesisgenes %>% select(-image)

    respect <- respect %>% select(-image)

    write.csv(respect, "../results/22_respect.csv")
    write.csv(hypothesisgenes2, "../results/22_respectgenes.csv" )

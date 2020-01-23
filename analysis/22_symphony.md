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

    modules <- read_csv("../results/08_genes_modules.csv") %>%
      rename(modulecolor = `net$colors`)  

    ## Parsed with column specification:
    ## cols(
    ##   `net$colors` = col_character(),
    ##   gene = col_character()
    ## )

    candidategenes <- c("PRL", "PRLR", 
                     "VIP", "VIPR1", "VIPR2", 
                     "OXT", "AVP", "AVPR1A", "AVPR1B", 
                     "GNRH1","GNRHR", "NPVF",
                     "NR3C1", "NR3C2",
                     "ESR1", "ESR2", "AR",
                     "DIO2","LEPR", "DIO3", "DIO1","CYP19A1",
                     "HSPA14", "HSPA12A",
                     "PTGES3", "HSD11B2","DRD5", "DRD1", "DRD2","PGE1", "PGF")

    modulescandidates <- modules %>% filter(gene %in% candidategenes)  %>%
      group_by(modulecolor)  %>%
      summarize(gene = str_c(gene, collapse = ", "))

    modulescandidates

    ## # A tibble: 8 x 2
    ##   modulecolor gene                                                         
    ##   <chr>       <chr>                                                        
    ## 1 black       DRD5                                                         
    ## 2 green       PTGES3                                                       
    ## 3 greenyellow AVP, GNRH1, NPVF, OXT                                        
    ## 4 grey        DIO1, DIO3                                                   
    ## 5 magenta     PRLR                                                         
    ## 6 red         PRL                                                          
    ## 7 salmon      VIP                                                          
    ## 8 turquoise   AR, AVPR1A, AVPR1B, CYP19A1, DIO2, DRD1, ESR1, ESR2, GNRHR, …

    datapath <- "../results/"   # path to the data
    datafiles <- dir(datapath, pattern = "*allvsd.csv") # get file names
    datapathfile <- paste0(datapath, datafiles)

    df <- datapathfile %>%
      setNames(nm = .) %>%
      map_df(~read_csv(.x, col_types = cols(), col_names = T), .id = "filename") %>% 
      mutate(tissue = sapply(strsplit(as.character(filename),'../results/06_'), "[", 2)) %>% 
      mutate(tissue = sapply(strsplit(as.character(tissue),'allvsd.csv'), "[", 1))  %>% 
      select(tissue, X1, everything()) 

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Warning: Missing column names filled in: 'X1' [1]

    df2 <-  df  %>%
      filter(X1 %in% c(candidategenes)) %>%
      pivot_longer(cols = L.Blu13_male_gonad_control.NYNO:y98.o50.x_male_pituitary_inc.d3, 
                   names_to = "sample", values_to = "vsd") %>%
       mutate(sex = sapply(strsplit(as.character(sample),'_'), "[", 2)) %>%
       mutate(treatment = sapply(strsplit(as.character(sample),'_'), "[", 4))  %>%
       mutate(treatment = sapply(strsplit(as.character(treatment),'.NYNO'), "[", 1)) %>%
      mutate(bird = sapply(strsplit(as.character(sample),'_'), "[", 1)) %>%
      filter(treatment %in% charlevels) %>%
      select(bird, sex, treatment, tissue, X1, vsd) %>%
      mutate(tissue = fct_recode(tissue, "hypothalamus" = "hyp",
                        "pituitary" = "pit",
                        "gonads" = "gon"
                        )) %>%
      rename(gene = X1) %>%
      drop_na() %>%
      droplevels()
    head(df2)

    ## # A tibble: 6 x 6
    ##   bird    sex    treatment tissue gene    vsd
    ##   <chr>   <chr>  <chr>     <fct>  <chr> <dbl>
    ## 1 L.Blu13 male   control   gonads AR     7.64
    ## 2 L.G107  male   control   gonads AR     7.61
    ## 3 L.G118  female control   gonads AR     8.93
    ## 4 L.R3    male   control   gonads AR     7.85
    ## 5 L.R8    male   control   gonads AR     7.35
    ## 6 L.W33   male   control   gonads AR     7.57

    df2$treatment <- factor(df2$treatment, levels = alllevels)
    df2$tissue <- factor(df2$tissue, levels = tissuelevels)
    df2$gene <- factor(df2$gene)


    df3 <- df2 %>% 
      mutate(treatment = fct_relevel(treatment, charlevels)) %>% 
      group_by(treatment, tissue, gene)  %>% 
      summarize(m = mean(vsd, na.rm = T), se = sd(vsd,  na.rm = T)/sqrt(length(vsd))) %>%
      mutate(image = "../figures/images/DoveParentsRNAseq_note.png") 
    head(df3)  

    ## # A tibble: 6 x 6
    ## # Groups:   treatment, tissue [1]
    ##   treatment tissue      gene       m     se image                          
    ##   <fct>     <fct>       <fct>  <dbl>  <dbl> <chr>                          
    ## 1 control   hypothalam… AR      7.15 0.0800 ../figures/images/DoveParentsR…
    ## 2 control   hypothalam… AVP    10.0  0.313  ../figures/images/DoveParentsR…
    ## 3 control   hypothalam… AVPR1A  8.32 0.0993 ../figures/images/DoveParentsR…
    ## 4 control   hypothalam… AVPR1B  5.59 0.0462 ../figures/images/DoveParentsR…
    ## 5 control   hypothalam… CYP19…  9.06 0.188  ../figures/images/DoveParentsR…
    ## 6 control   hypothalam… DIO1    5.61 0.0541 ../figures/images/DoveParentsR…

    d4 <- left_join(df3, modules, by = "gene")

    ## Warning: Column `gene` joining factor and character vector, coercing into
    ## character vector

    for (i in levels(d4$tissue)) {
      p <-  d4 %>%
        filter(tissue == i) %>%
        ggplot(aes(x = treatment, y = m)) +
        geom_errorbar(aes(ymin=m-se, ymax=m+se, color=gene), width=.1) +
        geom_point(size = 1, aes(color = gene)) +
        geom_line(aes(x = as.numeric(treatment), y = m, color = gene)) +
        scale_alpha_manual(values = c(0.5, 1)) +
        labs(subtitle = i, y = "average expression", x = "parental stage") +
        facet_wrap(~modulecolor, nrow = 2) +
        theme_B3() +
        theme(legend.position = "bottom")
     print(p)
    }

![](../figures/favegenes/symphonymodules-1.png)![](../figures/favegenes/symphonymodules-2.png)![](../figures/favegenes/symphonymodules-3.png)

    d4 %>%
      filter(modulecolor != "turquoise",
             tissue == "pituitary") %>%
      droplevels() %>% 
        ggplot(aes(x = treatment, y = m)) +
        geom_errorbar(aes(ymin=m-se, ymax=m+se, color=gene), width=.1) +
        geom_point(size = 1, aes(color = gene)) +
        geom_line(aes(x = as.numeric(treatment), y = m, color = gene)) +
        scale_alpha_manual(values = c(0.5, 1)) +
        labs(subtitle = "WGCNA + canddiate genes in the pituitary", y = "average expression", x = "parental stage") +
        facet_wrap(~modulecolor, scales = "free_y", nrow = 2) +
        theme_B3() +
        theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
        guides(color = guide_legend(nrow = 2, byrow = T)) 

![](../figures/favegenes/symphonymodules-4.png)

    d4 %>%
      filter( modulecolor != "turquoise",
              # modulecolor != "greenyellow",
              tissue == "pituitary") %>%
      droplevels() %>% 
        ggplot(aes(x = treatment, y = m)) +
      geom_image(aes(image=image), size = 0.15) +
      facet_wrap(~modulecolor, scales = "free_y", nrow = 2) +
      labs(subtitle = "WGCNA + canddiate genes in the pituitary", y = "gene expression", x = "parental stage") +
        facet_wrap(~gene, scales = "free_y", nrow = ) +
        theme_B3() +
        theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) 

![](../figures/favegenes/symphonymodules-5.png)

    d4 %>%
      filter( modulecolor != "turquoise",
              # modulecolor != "greenyellow",
              tissue == "pituitary") %>%
      droplevels() %>% 
        ggplot(aes(x = treatment, y = m)) +
      geom_image(aes(image=image), size = 0.1) +
      labs(subtitle = "WGCNA + canddiate genes in the pituitary", y = "gene expression", x = "parental stage") +
        theme_B3() +
        theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_y_log10()

![](../figures/favegenes/symphonymodules-6.png)

outline
-------

-   data-driven
-   hypothesis-drive
-   sonically-driven
-   visually-driven
-   statistcally-driven

Data and hypothesis driven research for the visually, sonically, and
statistcially inclined biologists

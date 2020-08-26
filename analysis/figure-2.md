Fig 2
=====

    library(tidyverse)

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(cowplot)

    ## 
    ## ********************************************************

    ## Note: As of version 1.0.0, cowplot does not change the

    ##   default ggplot2 theme anymore. To recover the previous

    ##   behavior, execute:
    ##   theme_set(theme_cowplot())

    ## ********************************************************

    library(ggsignif)

    source("../R/themes.R")
    source("../R/functions.R")
    source("../R/genelists.R")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_double(),
    ##   geneid = col_double(),
    ##   NCBI = col_character(),
    ##   literature = col_character(),
    ##   GO = col_character(),
    ##   gene = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   x = col_character()
    ## )

    knitr::opts_chunk$set(echo = TRUE, fig.path = '../figures/')

Treatment specific DEGs
-----------------------

    allDEG <- read_csv("../results/03_allDEG.csv") %>%
      mutate(tissue = factor(tissue, levels = tissuelevel),
             direction = factor(direction, levels = alllevels),
             label = gsub("_","\nvs. ", comparison)) %>%
      filter(lfc > 1.1 | lfc < -1.1 )

    ## Parsed with column specification:
    ## cols(
    ##   sex = col_character(),
    ##   tissue = col_character(),
    ##   comparison = col_character(),
    ##   direction = col_character(),
    ##   gene = col_character(),
    ##   lfc = col_double(),
    ##   padj = col_double(),
    ##   logpadj = col_double()
    ## )

    head(allDEG)

    ## # A tibble: 6 x 9
    ##   sex    tissue comparison  direction gene       lfc     padj logpadj label     
    ##   <chr>  <fct>  <chr>       <fct>     <chr>    <dbl>    <dbl>   <dbl> <chr>     
    ## 1 female gonad  bldg_extend extend    CDK3     19.5  2.69e-20   19.6  "bldg\nvs…
    ## 2 female gonad  bldg_extend extend    CRISP2    6.52 3.48e- 3    2.46 "bldg\nvs…
    ## 3 female gonad  bldg_extend extend    KRT20     5.47 3.67e- 4    3.44 "bldg\nvs…
    ## 4 female gonad  bldg_extend extend    CLDN34    5.01 5.14e- 3    2.29 "bldg\nvs…
    ## 5 female gonad  bldg_extend extend    LOC1070…  4.89 7.51e- 2    1.12 "bldg\nvs…
    ## 6 female gonad  bldg_extend extend    OMD       3.53 5.38e- 4    3.27 "bldg\nvs…

    # for suppl figures
    DEGcontrol <- allDEG %>% 
      filter(grepl("control", comparison),
             !grepl("m.|early|extend|prolong", comparison))  %>%
      mutate(comparison = factor(comparison, levels = comparisonlevelscontrol)) %>%
      group_by(sex, tissue, comparison, direction, label) %>%
        summarise(n = n()) %>%
        mutate(n = ifelse(direction == "control", n*-1, n*1 ))


    DEGbldg <- allDEG %>% 
      filter(grepl("bldg", comparison),
             !grepl("m.|early|extend|prolong|control", comparison))  %>%
      drop_na() %>%
      mutate(comparison = factor(comparison, levels = comparisonlevelsbldg))  %>%
      group_by(sex, tissue, comparison, direction, label) %>%
        summarise(n = n()) %>%
        mutate(n = ifelse(direction == "bldg", n*-1, n*1 ))

    DEGchar <- allDEG %>% 
      filter(comparison %in% comparisonlevelschar,
             !grepl("control|bldg", comparison)) %>%
      mutate(comparison = factor(comparison, levels = comparisonlevelschar)) %>%
      mutate(updown = ifelse(lfc > 0, 1, -1)) %>%
      group_by(sex, tissue, comparison, direction, label, updown) %>%
      summarise(n = n()) %>%
      mutate(n = n*updown ) %>%
      select(-updown)

    DEGremove <- allDEG %>% 
      filter(comparison %in% comparisonlevelsremoval) %>%
      mutate(comparison = factor(comparison, levels = comparisonlevelsremoval)) %>%
      mutate(updown = ifelse(lfc > 0, 1, -1)) %>%
      group_by(sex, tissue, comparison, direction, label, updown) %>%
      summarise(n = n()) %>%
      mutate(n = n*updown ) %>%
      select(-updown)
    DEGremove

    ## # A tibble: 41 x 6
    ## # Groups:   sex, tissue, comparison, direction, label [41]
    ##    sex    tissue       comparison        direction label                       n
    ##    <chr>  <fct>        <fct>             <fct>     <chr>                   <dbl>
    ##  1 female hypothalamus inc.d3_m.inc.d3   inc.d3    "inc.d3\nvs. m.inc.d3"    -20
    ##  2 female hypothalamus inc.d3_m.inc.d3   m.inc.d3  "inc.d3\nvs. m.inc.d3"      9
    ##  3 female hypothalamus inc.d9_m.inc.d9   inc.d9    "inc.d9\nvs. m.inc.d9"     -1
    ##  4 female hypothalamus inc.d9_m.inc.d9   m.inc.d9  "inc.d9\nvs. m.inc.d9"     19
    ##  5 female hypothalamus inc.d17_m.inc.d17 inc.d17   "inc.d17\nvs. m.inc.d1…   -87
    ##  6 female hypothalamus inc.d17_m.inc.d17 m.inc.d17 "inc.d17\nvs. m.inc.d1…   109
    ##  7 female hypothalamus hatch_m.n2        hatch     "hatch\nvs. m.n2"        -157
    ##  8 female hypothalamus hatch_m.n2        m.n2      "hatch\nvs. m.n2"         156
    ##  9 female pituitary    inc.d3_m.inc.d3   inc.d3    "inc.d3\nvs. m.inc.d3"    -42
    ## 10 female pituitary    inc.d3_m.inc.d3   m.inc.d3  "inc.d3\nvs. m.inc.d3"     14
    ## # … with 31 more rows

    DEGreplace <- allDEG %>% 
      filter(comparison %in% comparisonlevelsreplace) %>%
      mutate(comparison = factor(comparison, levels = comparisonlevelsreplace)) %>%
      mutate(updown = ifelse(lfc > 0, 1, -1)) %>%
      group_by(sex, tissue, comparison, direction, label, updown) %>%
      summarise(n = n()) %>%
      mutate(n = n*updown ) %>%
      select(-updown)

    ## Warning: Factor `tissue` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`

    ## Warning: Factor `tissue` contains implicit NA, consider using
    ## `forcats::fct_explicit_na`

    DEGreplace

    ## # A tibble: 59 x 6
    ## # Groups:   sex, tissue, comparison, direction, label [59]
    ##    sex    tissue       comparison      direction label                      n
    ##    <chr>  <fct>        <fct>           <fct>     <chr>                  <dbl>
    ##  1 female hypothalamus hatch_early     early     "hatch\nvs. early"       127
    ##  2 female hypothalamus hatch_early     hatch     "hatch\nvs. early"      -149
    ##  3 female hypothalamus inc.d17_prolong inc.d17   "inc.d17\nvs. prolong"   -72
    ##  4 female hypothalamus inc.d17_prolong prolong   "inc.d17\nvs. prolong"    82
    ##  5 female hypothalamus hatch_prolong   prolong   "hatch\nvs. prolong"     168
    ##  6 female hypothalamus hatch_prolong   hatch     "hatch\nvs. prolong"    -182
    ##  7 female hypothalamus extend_hatch    hatch     "extend\nvs. hatch"      160
    ##  8 female hypothalamus extend_hatch    extend    "extend\nvs. hatch"     -149
    ##  9 female pituitary    hatch_early     early     "hatch\nvs. early"        74
    ## 10 female pituitary    hatch_early     hatch     "hatch\nvs. early"      -223
    ## # … with 49 more rows

Variance stabilized data
------------------------

    countData <- read_csv("../results/01_limma.csv") %>%
      column_to_rownames(var = "X1")

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )

    ## See spec(...) for full column specifications.

    countData <- as.data.frame(t(countData))
    head(countData[1:3])

    ##                                            A2ML1     A2ML2      A2ML3
    ## L.Blu13_male_gonad_control.NYNO         47.78082 4.7790384  236.34992
    ## L.Blu13_male_hypothalamus_control.NYNO 184.83464 4.4706096 6385.50191
    ## L.Blu13_male_pituitary_control.NYNO    157.88698 0.3126307  207.56568
    ## L.G107_male_gonad_control               48.12817 2.2530380  225.98641
    ## L.G107_male_hypothalamus_control       324.49068 4.3381750 8072.76202
    ## L.G107_male_pituitary_control           83.47082 0.3541458   67.54887

    colData <- read_csv("../metadata/00_colData.csv") %>%
      mutate(treatment = factor(treatment, levels = alllevels),
             tissue = factor(tissue, levels = tissuelevels)) %>% 
      column_to_rownames(var = "X1") 

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   X1 = col_character(),
    ##   V1 = col_character(),
    ##   bird = col_character(),
    ##   sex = col_character(),
    ##   tissue = col_character(),
    ##   treatment = col_character(),
    ##   external = col_character(),
    ##   internal = col_character(),
    ##   group = col_character(),
    ##   study = col_character()
    ## )

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
    ##                                        external internal
    ## L.Blu13_male_gonad_control.NYNO         control  control
    ## L.Blu13_male_hypothalamus_control.NYNO  control  control
    ## L.Blu13_male_pituitary_control.NYNO     control  control
    ## L.G107_male_gonad_control               control  control
    ## L.G107_male_hypothalamus_control        control  control
    ## L.G107_male_pituitary_control           control  control
    ##                                                            group
    ## L.Blu13_male_gonad_control.NYNO              male.gonads.control
    ## L.Blu13_male_hypothalamus_control.NYNO male.hypothalamus.control
    ## L.Blu13_male_pituitary_control.NYNO       male.pituitary.control
    ## L.G107_male_gonad_control                    male.gonads.control
    ## L.G107_male_hypothalamus_control       male.hypothalamus.control
    ## L.G107_male_pituitary_control             male.pituitary.control
    ##                                                  study
    ## L.Blu13_male_gonad_control.NYNO        charcterization
    ## L.Blu13_male_hypothalamus_control.NYNO charcterization
    ## L.Blu13_male_pituitary_control.NYNO    charcterization
    ## L.G107_male_gonad_control              charcterization
    ## L.G107_male_hypothalamus_control       charcterization
    ## L.G107_male_pituitary_control          charcterization

    # check ready for analysis
    # row.names(countData) == row.names(colData)
    head(row.names(countData) == row.names(colData))

    ## [1] TRUE TRUE TRUE TRUE TRUE TRUE

    vsd_path <- "../results/DEseq2/treatment/"   # path to the data
    vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
    vsd_pathfiles <- paste0(vsd_path, vsd_files)
    vsd_files

    ## [1] "female_gonad_vsd.csv"        "female_gonads_vsd.csv"      
    ## [3] "female_hypothalamus_vsd.csv" "female_pituitary_vsd.csv"   
    ## [5] "male_gonad_vsd.csv"          "male_gonads_vsd.csv"        
    ## [7] "male_hypothalamus_vsd.csv"   "male_pituitary_vsd.csv"

    allvsd <- vsd_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name")  %>% 
      dplyr::rename("gene" = "X1") 

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## Warning: Missing column names filled in: 'X1' [1]

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   X1 = col_character()
    ## )
    ## See spec(...) for full column specifications.

    ## before pivoting, check names of df with
    ## head(names(allvsd)) and tail(names(allvsd))
    allvsd <- allvsd %>% 
      pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
                   names_to = "samples", values_to = "counts")  %>% 
        drop_na() 
    allvsd <- as_tibble(allvsd)
    head(allvsd[2:4])

    ## # A tibble: 6 x 3
    ##   gene  samples                            counts
    ##   <chr> <chr>                               <dbl>
    ## 1 A2ML1 L.G118_female_gonad_control         14.7 
    ## 2 A2ML1 R.G106_female_gonad_control          7.79
    ## 3 A2ML1 R.R20_female_gonad_control          14.9 
    ## 4 A2ML1 R.R9_female_gonad_control           11.7 
    ## 5 A2ML1 R.W44_female_gonad_control          13.2 
    ## 6 A2ML1 blk.s031.pu.d_female_gonad_prolong   7.42

    tail(allvsd[2:4])

    ## # A tibble: 6 x 3
    ##   gene  samples                          counts
    ##   <chr> <chr>                             <dbl>
    ## 1 ZZZ3  y18.x_male_pituitary_m.inc.d3     10.0 
    ## 2 ZZZ3  y4.x_male_pituitary_m.inc.d17      9.88
    ## 3 ZZZ3  y55.x_male_pituitary_m.inc.d8      9.91
    ## 4 ZZZ3  y63.x_male_pituitary_m.inc.d9      9.92
    ## 5 ZZZ3  y95.g131.x_male_pituitary_inc.d9   9.79
    ## 6 ZZZ3  y98.o50.x_male_pituitary_inc.d3    9.89

    # for plotting gene expression over time
    getcandidatevsdmanip <- function(whichgenes, whichtissue){
      
      candidates  <- allvsd %>%
        filter(gene %in% whichgenes) %>%
        dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1),
                      sextissue = sapply(strsplit(sextissue, '../results/DEseq2/treatment/'), "[", 2),
                      sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                    sex = factor(sex),
                      tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                    treatment = sapply(strsplit(samples, '\\_'), "[", 4),
                    treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1),
                    treatment = recode(treatment, "m.hatch" = "m.n2", "extend.hatch" = "m.n2",
                                                       "inc.prolong" = "prolong", "m.inc.d8" = "early"),
                    treatment =  factor(treatment, levels = alllevels),
                    earlylate = fct_collapse(treatment, 
                                      early = c("early", "m.inc.d3", "m.inc.d9",
                                               "lay", "inc.d3", "inc.d9" ),
                                      
                                      late = c("inc.d17", "m.inc.d17", "prolong" ,  "hatch" , 
                                             "m.n2", "extend", "n5", "n9"),
                                      reference = c("control", "bldg")),
                    extint = fct_collapse(treatment, 
                                      loss = c("m.inc.d3", "m.inc.d9", "m.n2", "m.inc.d17"),
                                      eggs = c("lay","inc.d3", "inc.d9", "inc.d17", "prolong"),
                                      chicks = c("early", "hatch", "extend", "n5", "n9"),
                                      reference = c("control", "bldg")))  %>%
        dplyr::select(gene, counts, sex, tissue, treatment, earlylate, extint, samples) 
      return(candidates)
    }

    candidatevsd <- getcandidatevsdmanip(candidategenes)

    hypvsd <- candidatevsd %>% filter(tissue %in% "hypothalamus")  
    pitvsd <- candidatevsd %>% filter(tissue %in% "pituitary") 
    gonvsd <- candidatevsd %>% filter(tissue %in% "gonad") 

    hypvsdf <- candidatevsd %>% filter(tissue %in% "hypothalamus", sex == "female")  
    pitvsdf <- candidatevsd %>% filter(tissue %in% "pituitary", sex == "female")
    gonvsdf <- candidatevsd %>% filter(tissue %in% "gonad", sex == "female")  
    hypvsdm <- candidatevsd %>% filter(tissue %in% "hypothalamus", sex == "male")  
    pitvsdm <- candidatevsd %>% filter(tissue %in% "pituitary", sex == "male") 
    gonvsdm <- candidatevsd %>% filter(tissue %in% "gonad", sex == "male") 

    summary(hypvsd)

    ##      gene               counts           sex           tissue         
    ##  Length:37908       Min.   : 5.547   female:18954   Length:37908      
    ##  Class :character   1st Qu.: 6.408   male  :18954   Class :character  
    ##  Mode  :character   Median : 7.861                  Mode  :character  
    ##                     Mean   : 7.957                                    
    ##                     3rd Qu.: 9.187                                    
    ##                     Max.   :15.891                                    
    ##                                                                       
    ##    treatment         earlylate           extint        samples         
    ##  control: 2574   reference: 4914   reference: 4914   Length:37908      
    ##  inc.d9 : 2574   early    :13806   eggs     :12168   Class :character  
    ##  inc.d17: 2574   late     :19188   loss     : 9126   Mode  :character  
    ##  m.n2   : 2574                     chicks   :11700                     
    ##  n9     : 2574                                                         
    ##  bldg   : 2340                                                         
    ##  (Other):22698

    summary(pitvsd)

    ##      gene               counts           sex           tissue         
    ##  Length:38280       Min.   : 5.096   female:18975   Length:38280      
    ##  Class :character   1st Qu.: 6.599   male  :19305   Class :character  
    ##  Mode  :character   Median : 7.937                  Mode  :character  
    ##                     Mean   : 8.207                                    
    ##                     3rd Qu.: 9.356                                    
    ##                     Max.   :21.505                                    
    ##                                                                       
    ##      treatment         earlylate           extint        samples         
    ##  control  : 2903   reference: 5223   reference: 5223   Length:38280      
    ##  inc.d9   : 2782   early    :13918   eggs     :12294   Class :character  
    ##  inc.d17  : 2552   late     :19139   loss     : 9163   Mode  :character  
    ##  m.n2     : 2552                     chicks   :11600                     
    ##  n9       : 2552                                                         
    ##  m.inc.d17: 2435                                                         
    ##  (Other)  :22504

    g <- plotcandidatechar(hypvsdf, "HTR2C") + labs(subtitle = "Hypothalamus")
    h <- plotcandidatechar(hypvsdm, "HTR2C") + labs(y = NULL, x = "") + labs(subtitle = " ")
    i <- plotremoval(hypvsdf, "HTR2C")+ labs(y = NULL) + labs(subtitle = " ")
    j <- plotremoval(hypvsdm, "HTR2C") + labs(y = NULL, x = "")+ labs(subtitle = " ")
    k <- plotreplacement(hypvsdf, "HTR2C") + labs(y = NULL)+ labs(subtitle = " ")
    l <- plotreplacement(hypvsdm, "HTR2C") + labs(y = NULL, x = "")+ labs(subtitle = " ")

    a <- plot.volcano("hypothalamus", sexlevels,  "control_bldg") + 
      facet_wrap(~sex) 

    b <- makebargraphv3(DEGcontrol, "hypothalamus","No. of DEGs\n(-) decreased  increased (+)", comparisonlabelscontrol) +
      labs(x = "Control versus all other reproductive and parental stages") +
      geom_rect(mapping=aes(xmin=0.5, xmax=1.5, ymin=-1000, ymax = 350, fill = F), color="black", alpha=0.5) +
      annotate("text", x = 1, y = -1075, label = "D", size = 2.5)   

    c <- makebargraphv3(DEGbldg, "hypothalamus", "No. of DEGs\n (-)decreased  increased (+)", comparisonlabelsbldg) +
      labs(x = "Nest-building versus all other parental stages")
    d <- makebargraphv3(DEGchar, "hypothalamus", NULL,  comparisonlabelscharnobldg) +
      labs(x = "Comparison of sequential parental stages")

    e <- makebargraphv3(DEGremove, "hypothalamus", "No. of DEGs\n (-)decreased  increased (+)", comparisonlevelsremoval) +
      labs(x = "Offspring removal versus temporal control")
    f <- makebargraphv3(DEGreplace, "hypothalamus", NULL,  comparisonlevelsreplace) +
      labs(x = "Offspring removal versus temporal or external control")


    m <- png::readPNG("../figures/venn-eggs-hyp.png")
    m <- ggdraw() +  draw_image(m, scale = 1)

    n <- png::readPNG("../figures/venn-chicks-hyp.png")
    n <- ggdraw() +  draw_image(n, scale = 1)

    ghi <- plot_grid(g,h,k,l,nrow = 1,
                    labels = c("A", "", "B"), label_size = 8, hjust = 0,
                    rel_widths = c(9,9,6,6))

    ab <- plot_grid(a,b,rel_widths = c(1,2.5), 
                    labels = c("C", "D"), label_size = 8, hjust = 0)

    cd <- plot_grid(c,d,rel_widths = c(1.1,1), align = "h",
                    labels = c("E", "F"), label_size = 8, hjust = 0)

    ef <- plot_grid(f,m,n, rel_widths = c(2,1,1), nrow =1,
                    labels = c("G", "H", "I"), label_size = 8, hjust = 0)


    fig2 <- plot_grid(ghi,ab,cd,ef,  ncol = 1)
    fig2

![](../figures/fig2-1.png)

    a <- plot.volcano("pituitary", sexlevels,  "control_bldg") + 
      facet_wrap(~sex) 

    b <- makebargraphv3(DEGcontrol, "pituitary","No. of DEGs\n(-) decreased  increased (+)", comparisonlabelscontrol) +
      labs(x = "Control versus all other reproductive and parental stages") 

    c <- makebargraphv3(DEGbldg, "pituitary", "No. of DEGs\n (-)decreased  increased (+)", comparisonlabelsbldg) +
      labs(x = "Nest-building versus all other parental stages")
    d <- makebargraphv3(DEGchar, "pituitary", NULL,  comparisonlabelscharnobldg) +
      labs(x = "Comparison of sequential parental stages")

    e <- makebargraphv3(DEGremove, "pituitary", "No. of DEGs\n (-)decreased  increased (+)", comparisonlevelsremoval) +
      labs(x = "Offspring removal versus temporal control")
    f <- makebargraphv3(DEGreplace, "pituitary", NULL,  comparisonlablelssreplace) +
      labs(x = "Offspring removal versus temporal or external control")

    g <- plotcandidatechar(pitvsdf, "PRL") + labs(subtitle = "Pituitary")
    h <- plotcandidatechar(pitvsdm, "PRL") + labs(y = NULL, x = "") + labs(subtitle = " ")
    i <- plotremoval(pitvsdf, "PRL")+ labs(y = NULL) + labs(subtitle = " ")
    j <- plotremoval(pitvsdm, "PRL") + labs(y = NULL, x = "")+ labs(subtitle = " ")
    k <- plotreplacement(pitvsdf, "PRL") + labs(y = NULL)+ labs(subtitle = " ")
    l <- plotreplacement(pitvsdm, "PRL") + labs(y = NULL, x = "")+ labs(subtitle = " ")

    m <- png::readPNG("../figures/venn-eggs-pit.png")
    m <- ggdraw() +  draw_image(m, scale = 1)

    n <- png::readPNG("../figures/venn-chicks-pit.png")
    n <- ggdraw() +  draw_image(n, scale = 1)

    ghi <- plot_grid(g,h,k,l,nrow = 1,
                    labels = c("A", "", "B"), label_size = 8, hjust = 0,
                    rel_widths = c(9,9,6,6))

    ab <- plot_grid(a,b,rel_widths = c(1,2.5), 
                    labels = c("C", "D"), label_size = 8, hjust = 0)

    cd <- plot_grid(c,d,rel_widths = c(1.1,1), align = "h",
                    labels = c("E", "F"), label_size = 8, hjust = 0)

    ef <- plot_grid(f,m,n, rel_widths = c(2,1,1), nrow =1,
                    labels = c("G", "H", "I"), label_size = 8, hjust = 0)


    fig3 <- plot_grid(ghi,ab,cd,ef,  ncol = 1)
    fig3

![](../figures/fig3-1.png)

    a <- plot.volcano("gonad", sexlevels,  "control_bldg") + 
      facet_wrap(~sex) 

    b <- makebargraphv3(DEGcontrol, "gonad","No. of DEGs\n(-) decreased  increased (+)", comparisonlabelscontrol) +
      labs(x = "Control versus all other reproductive and parental stages") 

    c <- makebargraphv3(DEGbldg, "gonad", "No. of DEGs\n (-)decreased  increased (+)", comparisonlabelsbldg) +
      labs(x = "Nest-building versus all other parental stages")
    d <- makebargraphv3(DEGchar, "gonad", NULL,  comparisonlabelscharnobldg) +
      labs(x = "Comparison of sequential parental stages")

    e <- makebargraphv3(DEGremove, "gonad", "No. of DEGs\n (-)decreased  increased (+)", comparisonlevelsremoval) +
      labs(x = "Offspring removal versus temporal control")
    f <- makebargraphv3(DEGreplace, "gonad", NULL,  comparisonlablelssreplace) +
      labs(x = "Offspring removal versus temporal or external control")

    g <- plotcandidatechar(pitvsdf, "ESR1") + labs(subtitle = "Gonads")
    h <- plotcandidatechar(pitvsdm, "ESR1") + labs(y = NULL, x = "") + labs(subtitle = " ")
    i <- plotremoval(pitvsdf, "ESR1")+ labs(y = NULL) + labs(subtitle = " ")
    j <- plotremoval(pitvsdm, "ESR1") + labs(y = NULL, x = "")+ labs(subtitle = " ")
    k <- plotreplacement(pitvsdf, "ESR1") + labs(y = NULL)+ labs(subtitle = " ")
    l <- plotreplacement(pitvsdm, "ESR1") + labs(y = NULL, x = "")+ labs(subtitle = " ")

    m <- png::readPNG("../figures/venn-eggs-gon.png")
    m <- ggdraw() +  draw_image(m, scale = 1)

    n <- png::readPNG("../figures/venn-chicks-gon.png")
    n <- ggdraw() +  draw_image(n, scale = 1)

    ghi <- plot_grid(g,h,k,l,nrow = 1,
                    labels = c("A", "", "B"), label_size = 8, hjust = 0,
                    rel_widths = c(9,9,6,6))

    ab <- plot_grid(a,b,rel_widths = c(1,2.5), 
                    labels = c("C", "D"), label_size = 8, hjust = 0)

    cd <- plot_grid(c,d,rel_widths = c(1.1,1), align = "h",
                    labels = c("E", "F"), label_size = 8, hjust = 0)

    ef <- plot_grid(f,m,n, rel_widths = c(2,1,1), nrow =1,
                    labels = c("G", "H", "I"), label_size = 8, hjust = 0)


    fig4 <- plot_grid(ghi,ab,cd,ef,  ncol = 1)
    fig4

![](../figures/fig4-1.png)

Save files
----------

    pdf(file="../figures/fig2-1.pdf", width=7, height=7)
    plot(fig2)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    png("../figures/fig2-1.png", width = 7, height = 7, 
        units = 'in', res = 300)
    plot(fig2) 
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/fig3-1.pdf", width=7, height=7)
    plot(fig3)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    png("../figures/fig3-1.png", width = 7, height = 7, 
        units = 'in', res = 300)
    plot(fig3) 
    dev.off()

    ## quartz_off_screen 
    ##                 2

    pdf(file="../figures/fig4-1.pdf", width=7, height=7)
    plot(fig4)
    dev.off()

    ## quartz_off_screen 
    ##                 2

    png("../figures/fig4-1.png", width = 7, height = 7, 
        units = 'in', res = 300)
    plot(fig4) 
    dev.off()

    ## quartz_off_screen 
    ##                 2

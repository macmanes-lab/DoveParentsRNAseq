    library(tidyverse)

    ## ── Attaching packages ──────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.2.1     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.3
    ## ✓ tidyr   1.0.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(cowplot)

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:ggplot2':
    ## 
    ##     ggsave

    library(readxl)
    library(modelr)
    library(lubridate)

    ## 
    ## Attaching package: 'lubridate'

    ## The following object is masked from 'package:base':
    ## 
    ##     date

    library(ggsignif)
    library(apaTables)

    source("../R/themes.R")  # load custom themes and color palletes
    source("../R/icons.R")

    ## Warning: Column `icons` joining factor and character vector, coercing into
    ## character vector

    knitr::opts_chunk$set(fig.path = '../figures/hormones/',message=F, warning=FALSE, cache = T)

    colData  <- read_csv("../metadata/00_birds.csv") %>%
      mutate(RNAseq = "RNAseq",
             bird_id = bird)  %>%
       select(-X1, -bird)
    colData

    ## # A tibble: 334 x 4
    ##    sex    treatment RNAseq bird_id
    ##    <chr>  <chr>     <chr>  <chr>  
    ##  1 male   control   RNAseq L.Blu13
    ##  2 male   control   RNAseq L.G107 
    ##  3 female control   RNAseq L.G118 
    ##  4 male   control   RNAseq L.R3   
    ##  5 male   control   RNAseq L.R8   
    ##  6 male   control   RNAseq L.W33  
    ##  7 male   control   RNAseq L.W3   
    ##  8 male   control   RNAseq L.W4   
    ##  9 female control   RNAseq R.G106 
    ## 10 female control   RNAseq R.R20  
    ## # … with 324 more rows

    prolactin <- read_excel("../results/Pigeon prolactin concentrations juil 2018.xlsx", sheet = 1) %>% 
      filter(Study %in% c("Baseline", "ParentalCare")) %>%
        dplyr::mutate(sex = fct_recode(Sex,
                                "female" = "f",
                                "male" = "m"),
               treatment = fct_recode(Treatment,
                                "hatch" = "Hatch",
                                "inc.d17" = "Inc_d17",
                                "inc.d17" = "inc_d17",
                                "inc.d3" = "Inc_d3",
                                 "inc.d3" = "inc_d3",
                                "inc.d9" = "Inc_d9",
                                 "inc.d9" = "inc_d9",
                                "m.inc.d9" = "M_Inc9",
                                "m.inc.d9" = "M_inc9",
                                "m.inc.d3" = "M_Inc3",
                                "m.inc.d8" = "M_Inc8",
                                "m.inc.d8" = "M_inc8",
                                "m.inc.d17" = "M_Inc17",
                                "m.n2" = "M_hatch",
                                "control" = "baseline",
                                "n5" = "N5", 
                                "n9" = "N9"),
               study = fct_collapse(treatment,
                                     characterization = charlevels,
                                     manipulation = maniplevels1)) %>%
              dplyr::rename("plasma_conc" = "Prolactin ng/mL") %>%
              mutate(bird_id = gsub("[[:punct:]]", "." , ColorBands)) %>% 
              dplyr::mutate(hormone = "prolactin") %>% 
              dplyr::select(study, treatment, sex, bird_id, hormone, plasma_conc)  %>% 
              drop_na()
    head(prolactin)

    ## # A tibble: 6 x 6
    ##   study            treatment sex    bird_id   hormone   plasma_conc
    ##   <fct>            <fct>     <fct>  <chr>     <chr>           <dbl>
    ## 1 characterization control   male   x.g       prolactin        3.83
    ## 2 characterization control   male   x.g.g     prolactin        3.28
    ## 3 characterization control   male   x.blk.blk prolactin        4.15
    ## 4 characterization control   male   x.g.g.g   prolactin       25.3 
    ## 5 characterization control   female x.g.g.f   prolactin       21.5 
    ## 6 characterization control   male   x.blu.o   prolactin       14.9

    PETC <- read_excel("../results/hormones.xlsx", sheet = 1)  %>% 
                  dplyr::rename(corticosterone = cort, progesterone = p4,
                                estradiol = e2,  testosterone = t) %>% 
                   dplyr::mutate(treatment = fct_recode(treatment...3,
                                "inc.d17" = "incd17",
                                "inc.d3" = "incd3",
                                "inc.d9" = "incd9",
                                "m.inc.d9" = "minc9",
                                "m.inc.d3" = "minc3",
                                "m.inc.d8" = "minc8",
                                "m.inc.d17" = "minc17",
                                "m.n2" = "m hatch",
                                "control" = "baseline")) %>%
                    dplyr::filter(treatment %in% alllevels2) %>%   
                    dplyr::mutate(sex = fct_recode(sex,
                                "female" = "f",
                                "male" = "m"),
                           study = fct_collapse(treatment,
                                    characterization = charlevels,
                                    manipulation = maniplevels1),
                           moltbin = as.integer(moltbin)) %>% 
                  dplyr::mutate(bird_id = gsub("[[:punct:]]", "." , id)) %>% 
                  dplyr::select(bird_id, treatment, sex, study, 
                                corticosterone, estradiol, testosterone, progesterone, moltbin) %>% 
                  pivot_longer(cols = corticosterone:moltbin,
                               names_to = "hormone", values_to = "plasma_conc",
                               values_drop_na = TRUE) %>% 
                  drop_na() %>% droplevels()  
    head(PETC)

    ## # A tibble: 6 x 6
    ##   bird_id   treatment sex   study            hormone        plasma_conc
    ##   <chr>     <fct>     <fct> <fct>            <chr>                <dbl>
    ## 1 blk.d.g.s n9        male  characterization corticosterone        1.31
    ## 2 blk.d.g.s n9        male  characterization testosterone          4.32
    ## 3 blk.d.g.s n9        male  characterization moltbin               1   
    ## 4 blk.s.o.g prolong   male  manipulation     corticosterone        1.75
    ## 5 blk.s.o.g prolong   male  manipulation     testosterone          1.10
    ## 6 blk.s.o.g prolong   male  manipulation     moltbin               1

    hormones <- rbind(prolactin, PETC)

    hormones$okay <- ifelse(hormones$hormone == "corticosterone" & hormones$plasma_conc > 30, "bad",
                        ifelse(hormones$hormone == "progesterone" & hormones$plasma_conc > 5, "bad", 
                               ifelse(hormones$hormone == "prolactin" & hormones$plasma_conc > 150, "bad", 
                            ifelse(hormones$hormone == "testosterone" & hormones$sex == "female", "bad",
                                   ifelse(hormones$hormone == "estradiol" & hormones$sex == "male", "bad", "okay")))))
    hormones <- hormones %>% filter(okay == "okay") %>% droplevels() %>% select(-okay)


    # make a winder one for correlations
    hormoneswide <- hormones %>% pivot_wider(names_from = "hormone", values_from = "plasma_conc", 
                                             values_fn = list(plasma_conc = mean))
    hormoneswide$treatment <- factor(hormoneswide$treatment, levels = alllevels)
    head(hormoneswide)

    ## # A tibble: 6 x 10
    ##   study treatment sex   bird_id prolactin corticosterone testosterone
    ##   <fct> <fct>     <fct> <chr>       <dbl>          <dbl>        <dbl>
    ## 1 char… control   male  x.g          3.83          0.590        0.187
    ## 2 char… control   male  x.g.g        3.28          0.817        0.661
    ## 3 char… control   male  x.blk.…      4.15          0.850        1.00 
    ## 4 char… control   male  x.g.g.g     25.3           1.03         1.74 
    ## 5 char… control   fema… x.g.g.f     21.5           1.15        NA    
    ## 6 char… control   male  x.blu.o     14.9           1.44         3.05 
    ## # … with 3 more variables: moltbin <dbl>, estradiol <dbl>,
    ## #   progesterone <dbl>

    # for faceting
    hormones <- hormones %>% filter(hormone != "moltbin")
    hormones$treatment <- factor(hormones$treatment, levels = alllevels)
    hormones$hormone <- factor(hormones$hormone, 
                               levels = c("prolactin", "corticosterone", "progesterone", 
                                          "estradiol", "testosterone"))
    # for stats
    hormones$logconc <- log10(hormones$plasma_conc)
    head(hormones)

    ## # A tibble: 6 x 7
    ##   study            treatment sex    bird_id   hormone   plasma_conc logconc
    ##   <fct>            <fct>     <fct>  <chr>     <fct>           <dbl>   <dbl>
    ## 1 characterization control   male   x.g       prolactin        3.83   0.583
    ## 2 characterization control   male   x.g.g     prolactin        3.28   0.516
    ## 3 characterization control   male   x.blk.blk prolactin        4.15   0.618
    ## 4 characterization control   male   x.g.g.g   prolactin       25.3    1.40 
    ## 5 characterization control   female x.g.g.f   prolactin       21.5    1.33 
    ## 6 characterization control   male   x.blu.o   prolactin       14.9    1.17

do control bird with high prolactin hormone have high PRL expression in the pituitary? yes.
===========================================================================================

    PRLpit <- read_csv("../results/10_PRLpit.csv") %>% 
      filter(treatment == "control") %>% 
      arrange(desc(PRL))
    head(PRLpit,2)

    ## # A tibble: 2 x 5
    ##   bird          sex    treatment tissue      PRL
    ##   <chr>         <chr>  <chr>     <chr>     <dbl>
    ## 1 blu.o.x.ATLAS female control   pituitary  20.9
    ## 2 L.W33         male   control   pituitary  19.9

summary for owen
----------------

    meanTE <- hormoneswide %>% 
      dplyr::group_by(treatment) %>%
      dplyr::summarise(meanTestosterone = (mean(testosterone, na.rm = TRUE)), 
                       meanEstradiol = mean(estradiol, na.rm = TRUE)) %>%
      dplyr::mutate(meanTestosterone = round(meanTestosterone,2)) %>%
      dplyr::mutate(meanEstradiol = round(meanEstradiol,2)) %>%
      dplyr::filter(treatment %in% charlevels) %>%
      droplevels() %>%
      dplyr::mutate(timepoint = as.numeric(treatment)) %>%
      select(timepoint, treatment, meanTestosterone, meanEstradiol)
    meanTE

    ## # A tibble: 9 x 4
    ##   timepoint treatment meanTestosterone meanEstradiol
    ##       <dbl> <fct>                <dbl>         <dbl>
    ## 1         1 control               0.96          0.12
    ## 2         2 bldg                  1.54          0.34
    ## 3         3 lay                   0.66          0.3 
    ## 4         4 inc.d3                0.59          0.33
    ## 5         5 inc.d9                0.81          0.15
    ## 6         6 inc.d17               1.39          0.25
    ## 7         7 hatch                 1.23          0.31
    ## 8         8 n5                    0.51          0.22
    ## 9         9 n9                    2.33          0.54

    write.csv(meanTE, "../results/07_meanTE.csv", row.names = F)
    write.csv(hormones, "../results/07_hormones.csv", row.names = F)
    write.csv(hormoneswide, "../results/07_hormoneswide.csv", row.names = F)

    #write.csv(hormones, "../../parentalhormones/data/hormones.csv", row.names = F)
    #write.csv(hormoneswide, "../../parentalhormones/data/hormoneswide.csv", row.names = F)

    library(tidyverse)

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    library(readr)
    library(corrr)
    library(forcats)

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

Wrangle data for hyotheses testing
----------------------------------

    # col data with internal (hiloPRL) and external  hypothesese
    colData <- read_csv("../metadata/04_colData.csv") %>%
      column_to_rownames(var = "V1")

    ## Parsed with column specification:
    ## cols(
    ##   V1 = col_character(),
    ##   bird = col_character(),
    ##   sex = col_character(),
    ##   tissue = col_character(),
    ##   treatment = col_character(),
    ##   group = col_character(),
    ##   study = col_character(),
    ##   sextissue = col_character(),
    ##   samples = col_character(),
    ##   hiloPRL = col_character(),
    ##   external = col_character()
    ## )

    tail(colData)

    ##                                          bird    sex       tissue
    ## y98.g54_female_gonad_m.hatch          y98.g54 female       gonads
    ## y98.g54_female_hypothalamus_m.hatch   y98.g54 female hypothalamus
    ## y98.g54_female_pituitary_m.hatch      y98.g54 female    pituitary
    ## y98.o50.x_male_gonad_inc.d3         y98.o50.x   male       gonads
    ## y98.o50.x_male_hypothalamus_inc.d3  y98.o50.x   male hypothalamus
    ## y98.o50.x_male_pituitary_inc.d3     y98.o50.x   male    pituitary
    ##                                     treatment                    group
    ## y98.g54_female_gonad_m.hatch             m.n2       female.gonads.m.n2
    ## y98.g54_female_hypothalamus_m.hatch      m.n2 female.hypothalamus.m.n2
    ## y98.g54_female_pituitary_m.hatch         m.n2    female.pituitary.m.n2
    ## y98.o50.x_male_gonad_inc.d3            inc.d3       male.gonads.inc.d3
    ## y98.o50.x_male_hypothalamus_inc.d3     inc.d3 male.hypothalamus.inc.d3
    ## y98.o50.x_male_pituitary_inc.d3        inc.d3    male.pituitary.inc.d3
    ##                                               study           sextissue
    ## y98.g54_female_gonad_m.hatch           manipulation       female_gonads
    ## y98.g54_female_hypothalamus_m.hatch    manipulation female_hypothalamus
    ## y98.g54_female_pituitary_m.hatch       manipulation    female_pituitary
    ## y98.o50.x_male_gonad_inc.d3         charcterization         male_gonads
    ## y98.o50.x_male_hypothalamus_inc.d3  charcterization   male_hypothalamus
    ## y98.o50.x_male_pituitary_inc.d3     charcterization      male_pituitary
    ##                                                                 samples
    ## y98.g54_female_gonad_m.hatch               y98.g54_female_gonad_m.hatch
    ## y98.g54_female_hypothalamus_m.hatch y98.g54_female_hypothalamus_m.hatch
    ## y98.g54_female_pituitary_m.hatch       y98.g54_female_pituitary_m.hatch
    ## y98.o50.x_male_gonad_inc.d3                 y98.o50.x_male_gonad_inc.d3
    ## y98.o50.x_male_hypothalamus_inc.d3   y98.o50.x_male_hypothalamus_inc.d3
    ## y98.o50.x_male_pituitary_inc.d3         y98.o50.x_male_pituitary_inc.d3
    ##                                     hiloPRL external
    ## y98.g54_female_gonad_m.hatch             hi     nest
    ## y98.g54_female_hypothalamus_m.hatch      hi     nest
    ## y98.g54_female_pituitary_m.hatch         hi     nest
    ## y98.o50.x_male_gonad_inc.d3              lo     eggs
    ## y98.o50.x_male_hypothalamus_inc.d3       lo     eggs
    ## y98.o50.x_male_pituitary_inc.d3          lo     eggs

### variance stabilized gene expression (vsd)

    vsd_path <- "../results/DEseq2/treatment/"   # path to the data
    vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
    vsd_pathfiles <- paste0(vsd_path, vsd_files)
    vsd_files

    ## [1] "female_gonad_vsd.csv"        "female_hypothalamus_vsd.csv"
    ## [3] "female_pituitary_vsd.csv"    "male_gonad_vsd.csv"         
    ## [5] "male_hypothalamus_vsd.csv"   "male_pituitary_vsd.csv"

    ## before pivoting, check names of df with
    ## head(names(allvsd)) and tail(names(allvsd))

    allvsd <- vsd_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name")  %>% 
      dplyr::rename("gene" = "X1") %>% 
      pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
                   names_to = "samples", values_to = "counts") 

    allvsd %>%  select(-file_name) %>% head()

    ## # A tibble: 6 x 3
    ##   gene  samples                            counts
    ##   <chr> <chr>                               <dbl>
    ## 1 A2ML1 L.G118_female_gonad_control         14.7 
    ## 2 A2ML1 R.G106_female_gonad_control          7.79
    ## 3 A2ML1 R.R20_female_gonad_control          14.9 
    ## 4 A2ML1 R.R9_female_gonad_control           11.7 
    ## 5 A2ML1 R.W44_female_gonad_control          13.2 
    ## 6 A2ML1 blk.s031.pu.d_female_gonad_prolong   7.42

    colData %>% filter(tissue == "gonads", treatment == "bldg")  

    ##               bird    sex tissue treatment              group
    ## 1          blk11.x female gonads      bldg female.gonads.bldg
    ## 2  blu114.r38.w198   male gonads      bldg   male.gonads.bldg
    ## 3      blu33.y88.x   male gonads      bldg   male.gonads.bldg
    ## 4     blu38.g135.x female gonads      bldg female.gonads.bldg
    ## 5       g104.w82.x   male gonads      bldg   male.gonads.bldg
    ## 6     g141.blu27.x female gonads      bldg female.gonads.bldg
    ## 7        g52.blu58   male gonads      bldg   male.gonads.bldg
    ## 8     l.s120.y.blk female gonads      bldg female.gonads.bldg
    ## 9      o38.blu29.x female gonads      bldg female.gonads.bldg
    ## 10      r.s056.g.o female gonads      bldg female.gonads.bldg
    ## 11      r.s059.d.o   male gonads      bldg   male.gonads.bldg
    ## 12         r83.g45 female gonads      bldg female.gonads.bldg
    ## 13   s.pu148.blk.r   male gonads      bldg   male.gonads.bldg
    ## 14    s063.d.blk.l female gonads      bldg female.gonads.bldg
    ## 15      s065.l.d.o   male gonads      bldg   male.gonads.bldg
    ## 16      s066.l.d.r   male gonads      bldg   male.gonads.bldg
    ## 17    s092.blk.r.o female gonads      bldg female.gonads.bldg
    ## 18      x.o30.g134   male gonads      bldg   male.gonads.bldg
    ## 19       x.r39.g10 female gonads      bldg female.gonads.bldg
    ## 20     x.r67.blu35   male gonads      bldg   male.gonads.bldg
    ##              study     sextissue                         samples hiloPRL
    ## 1  charcterization female_gonads       blk11.x_female_gonad_bldg      lo
    ## 2  charcterization   male_gonads blu114.r38.w198_male_gonad_bldg      lo
    ## 3  charcterization   male_gonads     blu33.y88.x_male_gonad_bldg      lo
    ## 4  charcterization female_gonads  blu38.g135.x_female_gonad_bldg      lo
    ## 5  charcterization   male_gonads      g104.w82.x_male_gonad_bldg      lo
    ## 6  charcterization female_gonads  g141.blu27.x_female_gonad_bldg      lo
    ## 7  charcterization   male_gonads       g52.blu58_male_gonad_bldg      lo
    ## 8  charcterization female_gonads  l.s120.y.blk_female_gonad_bldg      lo
    ## 9  charcterization female_gonads   o38.blu29.x_female_gonad_bldg      lo
    ## 10 charcterization female_gonads    r.s056.g.o_female_gonad_bldg      lo
    ## 11 charcterization   male_gonads      r.s059.d.o_male_gonad_bldg      lo
    ## 12 charcterization female_gonads       r83.g45_female_gonad_bldg      lo
    ## 13 charcterization   male_gonads   s.pu148.blk.r_male_gonad_bldg      lo
    ## 14 charcterization female_gonads  s063.d.blk.l_female_gonad_bldg      lo
    ## 15 charcterization   male_gonads      s065.l.d.o_male_gonad_bldg      lo
    ## 16 charcterization   male_gonads      s066.l.d.r_male_gonad_bldg      lo
    ## 17 charcterization female_gonads  s092.blk.r.o_female_gonad_bldg      lo
    ## 18 charcterization   male_gonads      x.o30.g134_male_gonad_bldg      lo
    ## 19 charcterization female_gonads     x.r39.g10_female_gonad_bldg      lo
    ## 20 charcterization   male_gonads     x.r67.blu35_male_gonad_bldg      lo
    ##    external
    ## 1      nest
    ## 2      nest
    ## 3      nest
    ## 4      nest
    ## 5      nest
    ## 6      nest
    ## 7      nest
    ## 8      nest
    ## 9      nest
    ## 10     nest
    ## 11     nest
    ## 12     nest
    ## 13     nest
    ## 14     nest
    ## 15     nest
    ## 16     nest
    ## 17     nest
    ## 18     nest
    ## 19     nest
    ## 20     nest

    allvsd %>% filter(samples == "blk11.x_female_gonad_bldg")

    ## # A tibble: 81,947 x 4
    ##    file_name                             gene    samples             counts
    ##    <chr>                                 <chr>   <chr>                <dbl>
    ##  1 ../results/DEseq2/treatment/female_g… A2ML1   blk11.x_female_gon…   7.60
    ##  2 ../results/DEseq2/treatment/female_g… A2ML2   blk11.x_female_gon…   5.33
    ##  3 ../results/DEseq2/treatment/female_g… A2ML3   blk11.x_female_gon…  10.2 
    ##  4 ../results/DEseq2/treatment/female_g… A2ML4   blk11.x_female_gon…   4.69
    ##  5 ../results/DEseq2/treatment/female_g… A4GALT  blk11.x_female_gon…   6.37
    ##  6 ../results/DEseq2/treatment/female_g… A4GNT   blk11.x_female_gon…   5.59
    ##  7 ../results/DEseq2/treatment/female_g… AAAS    blk11.x_female_gon…   8.55
    ##  8 ../results/DEseq2/treatment/female_g… AACS    blk11.x_female_gon…  10.4 
    ##  9 ../results/DEseq2/treatment/female_g… AADACL4 blk11.x_female_gon…   5.69
    ## 10 ../results/DEseq2/treatment/female_g… AADAT   blk11.x_female_gon…   7.52
    ## # … with 81,937 more rows

All DEGs (but currently only char DEGs)
---------------------------------------

    DEG_path <- "../results/DEseq2/hypothesis/"   # path to the data
    DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
    DEG_pathfiles <- paste0(DEG_path, DEG_files)
    #DEG_files

    allDEG <- DEG_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name") %>% 
      mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/hypothesis/'), "[", 2))  %>% 
      mutate(DEG = sapply(strsplit(as.character(DEG),'_DEGs.csv'), "[", 1))  %>% 
      mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
      mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%
      select(DEG, sex, tissue, direction, everything())   %>%
      mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%
      mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%
      mutate(comparison = paste(down,up, sep = "_")) %>%
      mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
      mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%
      dplyr::select(sex,tissue,comparison, direction, gene, lfc, padj, logpadj)  %>%
      mutate(tissue = factor(tissue)) %>% 
      drop_na()

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_character(),
    ##   logpadj = col_character(),
    ##   lfc = col_character(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_character(),
    ##   logpadj = col_character(),
    ##   lfc = col_character(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_character(),
    ##   logpadj = col_character(),
    ##   lfc = col_character(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_character(),
    ##   logpadj = col_character(),
    ##   lfc = col_character(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )
    ## Parsed with column specification:
    ## cols(
    ##   gene = col_character(),
    ##   padj = col_double(),
    ##   logpadj = col_double(),
    ##   lfc = col_double(),
    ##   sextissue = col_character(),
    ##   direction = col_character()
    ## )

    head(allDEG)

    ## # A tibble: 6 x 8
    ##   sex    tissue comparison  direction gene           lfc     padj logpadj
    ##   <chr>  <fct>  <chr>       <chr>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  eggs_chicks chicks    OVALX         6.17 7.05e- 3    2.15
    ## 2 female gonad  eggs_chicks chicks    LOC101749216  5.19 2.45e- 6    5.61
    ## 3 female gonad  eggs_chicks chicks    BPIFB2        4.95 5.23e- 3    2.28
    ## 4 female gonad  eggs_chicks chicks    OvoDA1        4.90 2.69e- 4    3.57
    ## 5 female gonad  eggs_chicks chicks    CDC20B        4.07 2.99e-10    9.52
    ## 6 female gonad  eggs_chicks chicks    WFIKKN2       4.01 1.90e- 2    1.72

gene lists
----------

    ## "internal clock" PRL genes

    hilogenes <- allDEG %>%
      filter(comparison == "lo_hi")  %>%
      drop_na() %>%
      arrange(padj) %>% 
      top_n(20) %>% pull(gene)

    ## Selecting by logpadj

    ## external environment gens 
    eggchickgenes <- allDEG %>%
      filter(comparison != "lo_hi")  %>%
      drop_na() %>%
      arrange(padj)  %>% 
      top_n(20) %>% pull(gene)

    ## Selecting by logpadj

### vsd for all and candidate, hypothesis, and data-driven genes

    # all list of candidate gens accumlated here

    candidategenes <- c(parentalcaregenes, WGCNAgenes, hilogenes, eggchickgenes, 
                        suszynskaagenes, shaidgenes)

    getcandidatevsd2 <- function(whichgenes, whichtissue, whichsex){
      candidates  <- allvsd %>%
        filter(gene %in% whichgenes) %>%
        dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
        dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/treatment/'), "[", 2)) %>%
        dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                      tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                      treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
        dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
        dplyr::select(sex, tissue, treatment, gene, samples, counts) %>%
        filter(tissue == whichtissue, sex %in% whichsex)  %>%
        drop_na()
      #candidates$treatment <- factor(candidates$treatment, levels = alllevels)
      return(candidates)
    }

    hypvsd <- getcandidatevsd2(candidategenes, "hypothalamus", sexlevels)
    pitvsd <- getcandidatevsd2(candidategenes, "pituitary", sexlevels)
    gonvsd <- getcandidatevsd2(candidategenes, "gonad", sexlevels)

    candidatevsd <- rbind(hypvsd, pitvsd)
    candidatevsd <- rbind(candidatevsd, gonvsd) 
    candidatevsd <- candidatevsd %>%
      mutate(treatment = factor(treatment)) %>%
      mutate(treatment = fct_recode(treatment, "early" = "m.inc.d8" )) %>%
      mutate(treatment = factor(treatment, alllevels))
    tail(candidatevsd)

    ## # A tibble: 6 x 6
    ##   sex   tissue treatment gene  samples                      counts
    ##   <chr> <chr>  <fct>     <chr> <chr>                         <dbl>
    ## 1 male  gonad  m.inc.d3  ZFX   y18.x_male_gonad_m.inc.d3      10.5
    ## 2 male  gonad  m.inc.d17 ZFX   y4.x_male_gonad_m.inc.d17      10.5
    ## 3 male  gonad  early     ZFX   y55.x_male_gonad_m.inc.d8      10.3
    ## 4 male  gonad  m.inc.d9  ZFX   y63.x_male_gonad_m.inc.d9      10.3
    ## 5 male  gonad  inc.d9    ZFX   y95.g131.x_male_gonad_inc.d9   10.5
    ## 6 male  gonad  inc.d3    ZFX   y98.o50.x_male_gonad_inc.d3    10.2

    write.csv(candidatevsd, "../results/06_candidatevsd.csv", row.names = F)
    write.csv(allDEG, "../results/06_allDEG.csv", row.names = F)

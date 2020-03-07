DEGs
----

    DEG_path <- "../results/DEseq2/"   # path to the data
    DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
    DEG_pathfiles <- paste0(DEG_path, DEG_files)
    #DEG_files

    allDEG <- DEG_pathfiles %>%
      setNames(nm = .) %>% 
      map_df(~read_csv(.x), .id = "file_name") %>% 
      mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/'), "[", 2))  %>% 
      mutate(DEG = sapply(strsplit(as.character(DEG),'_diffexp.csv'), "[", 1))  %>% 
      mutate(tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1)) %>%
      mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%
      mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%
      mutate(comparison = paste(down,up, sep = "_")) %>%
      mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
      mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%
    dplyr::select(sex,tissue,comparison, direction, gene, lfc, padj, logpadj) 
    head(allDEG)

    ## # A tibble: 6 x 8
    ##   sex    tissue comparison direction gene           lfc     padj logpadj
    ##   <chr>  <chr>  <chr>      <chr>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  bldg_lay   lay       LOC107053414  9.72 4.76e- 4    3.32
    ## 2 female gonad  bldg_lay   lay       MUC           5.82 2.72e- 3    2.57
    ## 3 female gonad  bldg_lay   lay       OVSTL         5.45 9.58e-10    9.02
    ## 4 female gonad  bldg_lay   lay       AOC1          4.41 2.75e- 3    2.56
    ## 5 female gonad  bldg_lay   lay       ETNPPL        4.25 4.76e- 5    4.32
    ## 6 female gonad  bldg_lay   lay       GKN2          3.99 1.63e- 2    1.79

    allDEG$tissue <- factor(allDEG$tissue , levels = tissuelevel)
    allDEG$comparison <- factor(allDEG$comparison , levels = comparisonlevels)
    allDEG$direction <- factor(allDEG$direction, levels = charlevels)

    makebargraph <- function(whichtissue, myylab, lowlim, higherlim){
      p <- allDEG %>%
        filter(tissue == whichtissue,
               comparison != "control_bldg") %>%
      ggplot(aes(x = comparison,  fill = direction)) +
        geom_bar(position = "dodge") +
        facet_grid(tissue~sex) +
        theme_B3() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")  +
        guides(fill = guide_legend(nrow = 1)) +
        labs(x = "Sequential parental care stage comparisons", 
             y = myylab,
             subtitle = " ") +
      scale_fill_manual(values = allcolors,
                           name = " ",
                           drop = FALSE) +
      scale_color_manual(values = allcolors) +
      geom_text(stat='count', aes(label=..count..), vjust =-0.5, 
                position = position_dodge(width = 1),
                size = 2, color = "black")  +
      ylim(lowlim, higherlim)
      return(p)
    }


    # hyp
    p1 <- makebargraph("hypothalamus","DEGs", 0, 1250) + theme(axis.text.x = element_blank(), 
                                                                   axis.title.x = element_blank())
    # pit
    p2 <- makebargraph("pituitary","DEGs", 0, 1250)  +  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
                                                                  strip.text.x = element_blank())
    # gon
    p3 <- makebargraph("gonad","DEGs", 0, 1250) +  theme(strip.text.x = element_blank())
    bcd <- plot_grid(p1,p2,p3, nrow = 3, rel_heights = c(1.2,1,1.5), labels = c("b", "c", "d"), label_size = 8)

    expdesign <- png::readPNG("../figures/images/fig_fig2a.png")
    expdesign <- ggdraw() +  draw_image(expdesign, scale = 1)

    plot_grid(expdesign, bcd, labels = c("a", " "), label_size = 8, nrow = 2, rel_heights = c(0.25,1))

![](../figures/fig2-1.png)

total degs
----------

    allDEG %>%
      group_by(sex, tissue, comparison) %>%
      summarize(totalDEGs = n()) %>%
      arrange(tissue, comparison)

    ## # A tibble: 37 x 4
    ## # Groups:   sex, tissue [6]
    ##    sex    tissue       comparison     totalDEGs
    ##    <chr>  <fct>        <fct>              <int>
    ##  1 female hypothalamus control_bldg        5683
    ##  2 male   hypothalamus control_bldg        6683
    ##  3 female hypothalamus bldg_lay               1
    ##  4 male   hypothalamus bldg_lay               1
    ##  5 female hypothalamus inc.d3_inc.d9          1
    ##  6 female hypothalamus inc.d9_inc.d17         5
    ##  7 male   hypothalamus inc.d9_inc.d17        87
    ##  8 female hypothalamus inc.d17_hatch          3
    ##  9 male   hypothalamus inc.d17_hatch          6
    ## 10 female hypothalamus hatch_n5            1927
    ## # … with 27 more rows

candidate genes
---------------

    geneids <- read_csv("../metadata/00_geneinfo.csv")

    ## Warning: Missing column names filled in: 'X1' [1]

    candidategenes <- c("OXT", "AVP", "GNRH1", "GNRHR", "CGNRH-R",
                        "AR", "POMC", "AGRP", 
                           "CRH", "AVPR1A", "AVPR1B", "AVPR2","VIP",
                           "CYP19A1", "DRD1", "DRD2", "PRL", "PRLR", "SOX9", 
                        "ESR1","ESR2", "LBH", "CDK1", "BRCA1",
                        "PTEN", "CREBBP", "FOS", "JUN", "EGR1",
                         "BDNF", "GRM2","GRIA1",
                        "KCNJ5", "CISH", "PTGER3", "CEBPD", "ZBTB16", 
                        "DIO3", "DIO2", "DIO1") 

    table1 <- allDEG %>%
      filter(gene %in% candidategenes,
             comparison != "control_bldg") %>%
        arrange(gene) %>%
      group_by(sex, tissue, comparison) %>%
      summarize(genes = str_c(gene, collapse = " ")) %>%
      pivot_wider(names_from = comparison, values_from = genes ) %>%
      select(sex, tissue, bldg_lay, lay_inc.d3, inc.d3_inc.d9,
            inc.d9_inc.d17,inc.d17_hatch, hatch_n5, n5_n9)  %>%
      arrange( tissue, sex)
    kable(table1)

<table>
<thead>
<tr>
<th style="text-align:left;">
sex
</th>
<th style="text-align:left;">
tissue
</th>
<th style="text-align:left;">
bldg\_lay
</th>
<th style="text-align:left;">
lay\_inc.d3
</th>
<th style="text-align:left;">
inc.d3\_inc.d9
</th>
<th style="text-align:left;">
inc.d9\_inc.d17
</th>
<th style="text-align:left;">
inc.d17\_hatch
</th>
<th style="text-align:left;">
hatch\_n5
</th>
<th style="text-align:left;">
n5\_n9
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
hypothalamus
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
BDNF BRCA1 CISH CYP19A1 DRD1 EGR1 GRIA1 POMC
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
male
</td>
<td style="text-align:left;">
hypothalamus
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
AR
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
pituitary
</td>
<td style="text-align:left;">
DIO3 ESR1 GNRHR
</td>
<td style="text-align:left;">
ESR1 GNRHR PTEN ZBTB16
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
BRCA1 CDK1 DIO2 KCNJ5 LBH PRL
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
AVPR2 BRCA1 CDK1 GRIA1 LBH PRL
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
male
</td>
<td style="text-align:left;">
pituitary
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
BRCA1 CDK1 CISH GRM2 LBH PRL VIP
</td>
<td style="text-align:left;">
GRM2
</td>
<td style="text-align:left;">
CEBPD ZBTB16
</td>
<td style="text-align:left;">
BRCA1 CDK1 CEBPD
</td>
</tr>
<tr>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
gonad
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
AGRP AVPR1A DIO2 EGR1 FOS PRLR PTGER3
</td>
<td style="text-align:left;">
AVPR1A
</td>
<td style="text-align:left;">
SOX9
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
male
</td>
<td style="text-align:left;">
gonad
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
SOX9
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
</tbody>
</table>

    write_csv(table1, "../results/table1.csv")

supplemental tables, all DEGS
-----------------------------

    suppletable1 <- allDEG %>%
      filter(comparison != "control_bldg") %>%
      group_by(sex, tissue, comparison) %>%
      arrange( tissue, sex, direction, gene)
    head(suppletable1)

    ## # A tibble: 6 x 8
    ## # Groups:   sex, tissue, comparison [3]
    ##   sex    tissue    comparison   direction gene         lfc     padj logpadj
    ##   <chr>  <fct>     <fct>        <fct>     <chr>      <dbl>    <dbl>   <dbl>
    ## 1 female hypothal… bldg_lay     bldg      HEMGN     -1.37  1.93e- 2    1.72
    ## 2 female hypothal… inc.d3_inc.… inc.d3    LOC1070… -17.2   9.56e-16   15.0 
    ## 3 female hypothal… inc.d9_inc.… inc.d9    CFAP44    -0.708 4.54e- 2    1.34
    ## 4 female hypothal… inc.d9_inc.… inc.d9    GMNN      -0.512 4.54e- 2    1.34
    ## 5 female hypothal… inc.d9_inc.… inc.d17   IGLL1      4.20  2.45e- 2    1.61
    ## 6 female hypothal… inc.d9_inc.… inc.d17   LOC1070…  17.8   9.74e-19   18.0

    suppletable1 %>%
      group_by(tissue) %>%
      summarize(totalDEGs = n())

    ## # A tibble: 3 x 2
    ##   tissue       totalDEGs
    ##   <fct>            <int>
    ## 1 hypothalamus      2032
    ## 2 pituitary         4440
    ## 3 gonad             3770

    write_csv(suppletable1, "../results/suppletable1.csv")


    head(allDEG) 

    ## # A tibble: 6 x 8
    ##   sex    tissue comparison direction gene           lfc     padj logpadj
    ##   <chr>  <fct>  <fct>      <fct>     <chr>        <dbl>    <dbl>   <dbl>
    ## 1 female gonad  bldg_lay   lay       LOC107053414  9.72 4.76e- 4    3.32
    ## 2 female gonad  bldg_lay   lay       MUC           5.82 2.72e- 3    2.57
    ## 3 female gonad  bldg_lay   lay       OVSTL         5.45 9.58e-10    9.02
    ## 4 female gonad  bldg_lay   lay       AOC1          4.41 2.75e- 3    2.56
    ## 5 female gonad  bldg_lay   lay       ETNPPL        4.25 4.76e- 5    4.32
    ## 6 female gonad  bldg_lay   lay       GKN2          3.99 1.63e- 2    1.79

    allDEG %>%
      drop_na()

    ## # A tibble: 46,764 x 8
    ##    sex    tissue comparison direction gene           lfc     padj logpadj
    ##    <chr>  <fct>  <fct>      <fct>     <chr>        <dbl>    <dbl>   <dbl>
    ##  1 female gonad  bldg_lay   lay       LOC107053414  9.72 4.76e- 4    3.32
    ##  2 female gonad  bldg_lay   lay       MUC           5.82 2.72e- 3    2.57
    ##  3 female gonad  bldg_lay   lay       OVSTL         5.45 9.58e-10    9.02
    ##  4 female gonad  bldg_lay   lay       AOC1          4.41 2.75e- 3    2.56
    ##  5 female gonad  bldg_lay   lay       ETNPPL        4.25 4.76e- 5    4.32
    ##  6 female gonad  bldg_lay   lay       GKN2          3.99 1.63e- 2    1.79
    ##  7 female gonad  bldg_lay   lay       COL10A1       3.69 4.29e- 2    1.37
    ##  8 female gonad  bldg_lay   lay       CA4           3.47 3.89e- 4    3.41
    ##  9 female gonad  bldg_lay   lay       GAL3ST2       3.37 1.16e- 2    1.93
    ## 10 female gonad  bldg_lay   lay       NEU4          3.12 8.98e- 3    2.05
    ## # … with 46,754 more rows

    write.csv(allDEG, "../../musicalgenes/data/allDEG.csv")

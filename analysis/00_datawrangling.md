    # https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
    library(tidyverse)

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0.9000     ✓ purrr   0.3.3     
    ## ✓ tibble  2.1.3          ✓ dplyr   0.8.3     
    ## ✓ tidyr   1.0.0          ✓ stringr 1.4.0     
    ## ✓ readr   1.3.1          ✓ forcats 0.4.0

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

    source("../R/themes.R")

Data Wrangling
==============

import Kallisto transcript data, make gene info file
----------------------------------------------------

    # import count data, set rows at entreziz
    kallistodata <- read.table("../results/kallistocounts.txt", 
                            sep = ",", row.names = NULL) %>%
      rename("NCBI" = "entrezid",
             "gene" = "Name")
    head(kallistodata)

    ##   row.names     gene geneid           NCBI L.Blu13_male_gonad_control.NYNO
    ## 1    408082    EDNRB 408082 NP_001001127.1                               6
    ## 2    408183  CYP26A1 408183 NP_001001129.1                               1
    ## 3    374073    CFDP1 374073 NP_001001189.1                             408
    ## 4    407777    AvBD7 407777 NP_001001194.1                               1
    ## 5    407779     KRT5 407779 NP_001001195.1                               0
    ## 6    408034 HSD11B1L 408034 NP_001001201.1                              65
    ##   L.Blu13_male_hypothalamus_control.NYNO
    ## 1                               37.00000
    ## 2                                0.00000
    ## 3                               72.00000
    ## 4                                0.00000
    ## 5                                0.00000
    ## 6                                2.72853
    ##   L.Blu13_male_pituitary_control.NYNO L.G107_male_gonad_control
    ## 1                                  46                     7.000
    ## 2                                  47                     2.000
    ## 3                                 197                   728.000
    ## 4                                   0                     1.000
    ## 5                                   0                     2.000
    ## 6                                  30                   146.948
    ##   L.G107_male_hypothalamus_control L.G107_male_pituitary_control
    ## 1                               43                       60.0000
    ## 2                                1                       25.0000
    ## 3                              159                      186.0000
    ## 4                                0                        0.0000
    ## 5                                0                        0.0000
    ## 6                                8                       18.5133
    ##   L.G118_female_gonad_control L.G118_female_hypothalamus_control.NYNO
    ## 1                          33                                      60
    ## 2                          94                                       0
    ## 3                         733                                      66
    ## 4                          96                                       0
    ## 5                           2                                       0
    ## 6                          53                                       3
    ##   L.G118_female_pituitary_control.NYNO L.R3_male_gonad_control.NYNO
    ## 1                             45.00000                            4
    ## 2                             14.00000                            0
    ## 3                            298.00000                          273
    ## 4                              0.00000                            0
    ## 5                              1.05712                            0
    ## 6                             38.00000                           40
    ##   L.R3_male_hypothalamus_control L.R3_male_pituitary_control.NYNO
    ## 1                            113                              114
    ## 2                              0                               37
    ## 3                            244                              249
    ## 4                              0                                0
    ## 5                              0                                0
    ## 6                             15                               26
    ##   L.R8_male_gonad_control L.R8_male_hypothalamus_control
    ## 1                  15.000                              6
    ## 2                   1.000                              0
    ## 3                1814.000                             23
    ## 4                   0.000                              0
    ## 5                   5.000                              0
    ## 6                 241.948                              0
    ##   L.R8_male_pituitary_control L.W33_male_gonad_control
    ## 1                          87                        3
    ## 2                          58                        1
    ## 3                         226                      301
    ## 4                           0                        0
    ## 5                           1                        1
    ## 6                          44                       37
    ##   L.W33_male_hypothalamus_control.NYNO L.W33_male_pituitary_control
    ## 1                            112.00000                     117.0000
    ## 2                              0.00000                      54.0000
    ## 3                            200.00000                     252.0000
    ## 4                              0.00000                       0.0000
    ## 5                              0.00000                       1.0000
    ## 6                              7.91449                      27.7011
    ##   L.W3_male_gonad_control.NYNO L.W3_male_hypothalamus_control
    ## 1                            1                       11.00000
    ## 2                            0                        0.00000
    ## 3                          247                       24.00000
    ## 4                            1                        0.00000
    ## 5                            1                        0.00000
    ## 6                           38                        1.13262
    ##   L.W3_male_pituitary_control L.W4_male_gonad_control.NYNO
    ## 1                          71                           10
    ## 2                          33                            1
    ## 3                         301                          345
    ## 4                           0                            0
    ## 5                           1                            1
    ## 6                          49                           57
    ##   L.W4_male_hypothalamus_control L.W4_male_pituitary_control
    ## 1                       111.0000                          68
    ## 2                         0.0000                           7
    ## 3                       427.0000                         246
    ## 4                         0.0000                           0
    ## 5                         0.0000                           0
    ## 6                        14.1608                          39
    ##   R.G106_female_gonad_control R.G106_female_hypothalamus_control
    ## 1                          95                           180.0000
    ## 2                         175                             2.0000
    ## 3                         823                           499.0000
    ## 4                          58                             0.0000
    ## 5                           2                             0.0000
    ## 6                          60                            27.2655
    ##   R.G106_female_pituitary_control R.R20_female_gonad_control
    ## 1                         98.0000                         61
    ## 2                          4.0000                         71
    ## 3                        176.0000                        320
    ## 4                          0.0000                         54
    ## 5                          0.0000                          0
    ## 6                         40.9282                         25
    ##   R.R20_female_hypothalamus_control.NYNO R.R20_female_pituitary_control
    ## 1                                     53                             34
    ## 2                                      0                              7
    ## 3                                    204                            226
    ## 4                                      0                              1
    ## 5                                      0                              2
    ## 6                                      6                             22
    ##   R.R9_female_gonad_control R.R9_female_hypothalamus_control
    ## 1                  39.10632                         270.0000
    ## 2                  76.00000                           1.0000
    ## 3                 327.00000                         448.0000
    ## 4                  43.00000                           0.0000
    ## 5                   1.00000                           2.0000
    ## 6                  12.20620                          16.1167
    ##   R.R9_female_pituitary_control.NYNO R.W44_female_gonad_control
    ## 1                             120.00                  159.40370
    ## 2                               2.00                   63.00000
    ## 3                             264.00                  775.00000
    ## 4                               0.00                  100.00000
    ## 5                               1.00                    6.60068
    ## 6                              17.21                   55.00000
    ##   R.W44_female_hypothalamus_control R.W44_female_pituitary_control.NYNO
    ## 1                          540.0000                                  48
    ## 2                            4.0000                                  33
    ## 3                          886.0000                                 141
    ## 4                            1.0000                                   0
    ## 5                            0.0000                                   1
    ## 6                           34.4786                                  19
    ##   R.Y108.W29_male_gonad_control R.Y108.W29_male_hypothalamus_control.NYNO
    ## 1                            10                                   49.0000
    ## 2                             0                                    1.0000
    ## 3                           954                                   69.0000
    ## 4                             0                                    0.0000
    ## 5                             5                                    0.0000
    ## 6                           153                                    6.0712
    ##   R.Y108.W29_male_pituitary_control blk.s030.o.g_male_gonad_prolong
    ## 1                                93                          52.000
    ## 2                                 7                           0.000
    ## 3                               256                        1127.000
    ## 4                                 0                           2.000
    ## 5                                 0                           7.000
    ## 6                                49                         349.972
    ##   blk.s030.o.g_male_hypothalamus_prolong
    ## 1                                77.0000
    ## 2                                 1.0000
    ## 3                               271.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                15.0411
    ##   blk.s030.o.g_male_pituitary_prolong blk.s031.pu.d_female_gonad_prolong
    ## 1                           306.00000                                268
    ## 2                            56.00000                                312
    ## 3                           550.00000                               1430
    ## 4                             0.00000                                 44
    ## 5                             3.70651                                 11
    ## 6                            96.27890                                178
    ##   blk.s031.pu.d_female_hypothalamus_prolong
    ## 1                                   86.0000
    ## 2                                    0.0000
    ## 3                                  253.0000
    ## 4                                    0.0000
    ## 5                                    0.0000
    ## 6                                   14.3261
    ##   blk.s031.pu.d_female_pituitary_prolong blk.s032.g.w_female_gonad_m.hatch
    ## 1                                412.000                            47.000
    ## 2                                134.000                           248.000
    ## 3                                417.000                          1017.000
    ## 4                                  1.000                            67.000
    ## 5                                  1.000                             4.000
    ## 6                                185.615                            95.941
    ##   blk.s032.g.w_female_hypothalamus_m.hatch
    ## 1                                 268.0000
    ## 2                                   2.0000
    ## 3                                 549.0000
    ## 4                                   2.0000
    ## 5                                   1.0000
    ## 6                                  43.4202
    ##   blk.s032.g.w_female_pituitary_m.hatch blk.s049.y.g_female_gonad_m.inc.d3
    ## 1                               336.000                                487
    ## 2                                12.000                                828
    ## 3                               669.000                               1609
    ## 4                                 0.000                                264
    ## 5                                 3.000                                 21
    ## 6                               288.056                                172
    ##   blk.s049.y.g_female_hypothalamus_m.inc.d3
    ## 1                                   78.0000
    ## 2                                    2.0000
    ## 3                                  241.0000
    ## 4                                    0.0000
    ## 5                                    0.0000
    ## 6                                   13.1815
    ##   blk.s049.y.g_female_pituitary_m.inc.d3
    ## 1                                267.000
    ## 2                                129.000
    ## 3                                494.000
    ## 4                                  0.000
    ## 5                                  1.000
    ## 6                                120.259
    ##   blk.s060.pu.w_female_gonad_m.inc.d3
    ## 1                                  68
    ## 2                                 149
    ## 3                                1079
    ## 4                                   7
    ## 5                                   7
    ## 6                                 106
    ##   blk.s060.pu.w_female_hypothalamus_m.inc.d3
    ## 1                                   215.4358
    ## 2                                     3.0000
    ## 3                                   545.0000
    ## 4                                     0.0000
    ## 5                                     1.0000
    ## 6                                    40.1182
    ##   blk.s060.pu.w_female_pituitary_m.inc.d3.NYNO
    ## 1                                      170.000
    ## 2                                      105.000
    ## 3                                      418.000
    ## 4                                        0.000
    ## 5                                        0.000
    ## 6                                      132.596
    ##   blk.s061.pu.y_female_gonad_inc.d9
    ## 1                                45
    ## 2                                77
    ## 3                               233
    ## 4                               154
    ## 5                                 1
    ## 6                                26
    ##   blk.s061.pu.y_female_hypothalamus_inc.d9
    ## 1                                 391.0000
    ## 2                                   4.0000
    ## 3                                 583.0000
    ## 4                                   0.0000
    ## 5                                   2.0000
    ## 6                                  55.2234
    ##   blk.s061.pu.y_female_pituitary_inc.d9 blk.y.l.s109_female_gonad_m.inc.d8
    ## 1                               85.2333                           190.8635
    ## 2                               61.0000                           329.0000
    ## 3                              189.0000                           691.0000
    ## 4                                0.0000                           151.0000
    ## 5                                0.0000                             5.0000
    ## 6                               54.9364                            75.0000
    ##   blk.y.l.s109_female_hypothalamus_m.inc.d8
    ## 1                                  285.0000
    ## 2                                    2.0000
    ## 3                                  528.0000
    ## 4                                    0.0000
    ## 5                                    1.0000
    ## 6                                   31.6834
    ##   blk.y.l.s109_female_pituitary_m.inc.d8 blk0.x_female_gonad_m.n2
    ## 1                              323.00000                       32
    ## 2                               27.00000                       38
    ## 3                              663.00000                      263
    ## 4                                0.00000                        1
    ## 5                                7.09367                        0
    ## 6                              138.42200                       37
    ##   blk0.x_female_hypothalamus_m.n2 blk0.x_female_pituitary_m.n2
    ## 1                         169.000                      89.0000
    ## 2                           1.000                      24.0000
    ## 3                         234.000                     261.0000
    ## 4                           0.000                       0.0000
    ## 5                           0.000                       1.0000
    ## 6                          17.462                      69.8865
    ##   blk11.x_female_gonad_bldg blk11.x_female_hypothalamus_bldg
    ## 1                        50                         334.0000
    ## 2                        32                           3.0000
    ## 3                       222                         594.0000
    ## 4                        38                           0.0000
    ## 5                         0                           1.0000
    ## 6                        22                          38.1538
    ##   blk11.x_female_pituitary_bldg blk12.x_male_gonad_n5
    ## 1                        206.00                 6.000
    ## 2                        137.00                 0.000
    ## 3                        701.00               355.000
    ## 4                          0.00                 0.000
    ## 5                          3.00                 1.000
    ## 6                        135.76               132.874
    ##   blk12.x_male_hypothalamus_n5.NYNO blk12.x_male_pituitary_n5
    ## 1                          171.0000                   256.000
    ## 2                            1.0000                    82.000
    ## 3                          567.0000                   757.000
    ## 4                            0.0000                     0.000
    ## 5                            0.0000                     1.000
    ## 6                           42.2418                   254.513
    ##   blk17.x_male_gonad_inc.d17 blk17.x_male_hypothalamus_inc.d17
    ## 1                          9                          249.0000
    ## 2                          3                            4.0000
    ## 3                        401                          681.0000
    ## 4                          0                            1.0000
    ## 5                          5                            2.0000
    ## 6                         83                           56.1321
    ##   blk17.x_male_pituitary_inc.d17 blk19.x_female_gonad_extend
    ## 1                        293.000                   171.00000
    ## 2                         58.000                   161.00000
    ## 3                        703.000                   510.00000
    ## 4                          0.000                    21.00000
    ## 5                          1.000                     8.32329
    ## 6                        162.923                    57.00000
    ##   blk19.x_female_hypothalamus_extend blk19.x_female_pituitary_extend
    ## 1                           245.0000                         123.000
    ## 2                             2.0000                           8.000
    ## 3                           746.0000                         286.000
    ## 4                             0.0000                           0.000
    ## 5                             1.0000                           4.000
    ## 6                            91.5648                         121.688
    ##   blk21.x_female_gonad_hatch blk21.x_female_hypothalamus_hatch
    ## 1                         39                               282
    ## 2                        135                                 0
    ## 3                        324                              1038
    ## 4                          2                                 0
    ## 5                          2                                 0
    ## 6                         27                                93
    ##   blk21.x_female_pituitary_hatch blk4.x_female_gonad_n9
    ## 1                        92.5414                     74
    ## 2                        16.0000                      9
    ## 3                       231.0000                    245
    ## 4                         0.0000                      1
    ## 5                         0.0000                      2
    ## 6                        70.8845                     95
    ##   blk4.x_female_hypothalamus_n9 blk4.x_female_pituitary_n9
    ## 1                       112.000                  183.00000
    ## 2                         3.000                   23.00000
    ## 3                      1148.000                  859.00000
    ## 4                         0.000                    0.00000
    ## 5                         0.000                    3.80532
    ## 6                        45.113                  206.00000
    ##   blk5.x_male_gonad_m.inc.d3 blk5.x_male_hypothalamus_m.inc.d3
    ## 1                     34.000                           95.0000
    ## 2                      4.000                            1.0000
    ## 3                   1375.000                          241.0000
    ## 4                      0.000                            0.0000
    ## 5                      4.000                            0.0000
    ## 6                    735.984                           12.0309
    ##   blk5.x_male_pituitary_m.inc.d3 blu.o.x.ATLAS_female_gonad_control
    ## 1                        270.000                                 45
    ## 2                         72.000                                123
    ## 3                        522.000                                602
    ## 4                          0.000                                  5
    ## 5                          1.000                                  2
    ## 6                        157.821                                 32
    ##   blu.o.x.ATLAS_female_hypothalamus_control
    ## 1                                        22
    ## 2                                         0
    ## 3                                        60
    ## 4                                         0
    ## 5                                         0
    ## 6                                         4
    ##   blu.o.x.ATLAS_female_pituitary_control blu10.w26.x_male_gonad_m.hatch
    ## 1                                142.000                       58.00000
    ## 2                                 34.000                        1.00000
    ## 3                                450.000                      680.00000
    ## 4                                  0.000                        5.00000
    ## 5                                  1.000                        9.37313
    ## 6                                 57.062                      342.00000
    ##   blu10.w26.x_male_hypothalamus_m.hatch blu10.w26.x_male_pituitary_m.hatch
    ## 1                              251.0000                             170.00
    ## 2                                3.0000                              59.00
    ## 3                              546.0000                             353.00
    ## 4                                0.0000                               0.00
    ## 5                                0.0000                               0.00
    ## 6                               42.4404                             102.69
    ##   blu103.x_female_gonad_hatch.NYNO blu103.x_female_hypothalamus_hatch
    ## 1                          101.000                            361.000
    ## 2                           35.000                              2.000
    ## 3                          805.000                            661.000
    ## 4                           91.000                              0.000
    ## 5                            6.000                              0.000
    ## 6                          171.977                             56.517
    ##   blu103.x_female_pituitary_hatch.NYNO blu104.w120.x_male_gonad_hatch
    ## 1                               234.00                         6.0000
    ## 2                                43.00                         0.0000
    ## 3                               482.00                       410.0000
    ## 4                                 0.00                         0.0000
    ## 5                                 0.00                         2.0000
    ## 6                               158.13                        84.8101
    ##   blu104.w120.x_male_hypothalamus_hatch
    ## 1                                   461
    ## 2                                     7
    ## 3                                   654
    ## 4                                     0
    ## 5                                     2
    ## 6                                    58
    ##   blu104.w120.x_male_pituitary_hatch.NYNO
    ## 1                                      96
    ## 2                                      19
    ## 3                                     426
    ## 4                                       0
    ## 5                                       0
    ## 6                                      84
    ##   blu108.w40.o158_male_gonad_inc.d9
    ## 1                            11.000
    ## 2                             0.000
    ## 3                           417.000
    ## 4                             0.000
    ## 5                             2.000
    ## 6                           237.985
    ##   blu108.w40.o158_male_hypothalamus_inc.d9
    ## 1                                  452.000
    ## 2                                    4.000
    ## 3                                  636.000
    ## 4                                    0.000
    ## 5                                    1.000
    ## 6                                   33.126
    ##   blu108.w40.o158_male_pituitary_inc.d9 blu111.w113.x_male_gonad_inc.d3
    ## 1                              108.1732                              11
    ## 2                                9.0000                               0
    ## 3                              203.0000                             415
    ## 4                                0.0000                               0
    ## 5                                1.0000                               0
    ## 6                               72.5762                             121
    ##   blu111.w113.x_male_hypothalamus_inc.d3
    ## 1                               163.0000
    ## 2                                 0.0000
    ## 3                               370.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                23.2068
    ##   blu111.w113.x_male_pituitary_inc.d3 blu113.w124.x_male_gonad_inc.d17
    ## 1                                 155                                6
    ## 2                                  68                                1
    ## 3                                 889                              380
    ## 4                                   0                                0
    ## 5                                   4                                0
    ## 6                                 149                              102
    ##   blu113.w124.x_male_hypothalamus_inc.d17
    ## 1                                728.0000
    ## 2                                  2.0000
    ## 3                                822.0000
    ## 4                                  0.0000
    ## 5                                  5.0000
    ## 6                                 46.7947
    ##   blu113.w124.x_male_pituitary_inc.d17.NYNO
    ## 1                                   219.000
    ## 2                                    22.000
    ## 3                                   535.000
    ## 4                                     0.000
    ## 5                                     4.000
    ## 6                                   140.482
    ##   blu114.r38.w198_male_gonad_bldg blu114.r38.w198_male_hypothalamus_bldg
    ## 1                           3.000                               312.6289
    ## 2                           0.000                                 4.0000
    ## 3                         527.000                               461.0000
    ## 4                           2.000                                 0.0000
    ## 5                           1.000                                 0.0000
    ## 6                         153.888                                30.1273
    ##   blu114.r38.w198_male_pituitary_bldg
    ## 1                                  95
    ## 2                                  93
    ## 3                                 281
    ## 4                                   0
    ## 5                                   0
    ## 6                                 109
    ##   blu115.y150.x_female_gonad_inc.prolong
    ## 1                                     31
    ## 2                                     48
    ## 3                                    331
    ## 4                                      4
    ## 5                                      1
    ## 6                                     66
    ##   blu115.y150.x_female_hypothalamus_inc.prolong
    ## 1                                      239.0000
    ## 2                                        6.0000
    ## 3                                      596.0000
    ## 4                                        1.0000
    ## 5                                        1.0000
    ## 6                                       47.3668
    ##   blu115.y150.x_female_pituitary_inc.prolong
    ## 1                                         71
    ## 2                                          9
    ## 3                                        206
    ## 4                                          0
    ## 5                                          0
    ## 6                                         61
    ##   blu119.w84.x_female_gonad_m.inc.d8
    ## 1                                 95
    ## 2                                 38
    ## 3                                843
    ## 4                                 57
    ## 5                                  2
    ## 6                                 94
    ##   blu119.w84.x_female_hypothalamus_m.inc.d8
    ## 1                                  314.0000
    ## 2                                    9.0000
    ## 3                                  999.0000
    ## 4                                    0.0000
    ## 5                                    4.0000
    ## 6                                   89.3176
    ##   blu119.w84.x_female_pituitary_m.inc.d8 blu121.w91.x_male_gonad_inc.d17
    ## 1                                221.000                           8.000
    ## 2                                 42.000                           0.000
    ## 3                                461.000                         395.000
    ## 4                                  0.000                           1.000
    ## 5                                  0.000                           1.000
    ## 6                                184.234                         117.874
    ##   blu121.w91.x_male_hypothalamus_inc.d17
    ## 1                               467.0000
    ## 2                                 1.0000
    ## 3                               893.0000
    ## 4                                 0.0000
    ## 5                                 3.0000
    ## 6                                66.1144
    ##   blu121.w91.x_male_pituitary_inc.d17 blu124.w180.x_female_gonad_hatch
    ## 1                             74.0000                               53
    ## 2                              5.0000                               70
    ## 3                            229.0000                              245
    ## 4                              0.0000                               15
    ## 5                              0.0000                                0
    ## 6                             45.4034                               28
    ##   blu124.w180.x_female_hypothalamus_hatch
    ## 1                                 471.000
    ## 2                                   0.000
    ## 3                                 714.000
    ## 4                                   0.000
    ## 5                                   2.000
    ## 6                                  51.193
    ##   blu124.w180.x_female_pituitary_hatch blu33.y88.x_male_gonad_bldg
    ## 1                              69.0000                     8.09308
    ## 2                               7.0000                     3.00000
    ## 3                             165.0000                   631.00000
    ## 4                               0.0000                     0.00000
    ## 5                               0.0000                     6.00000
    ## 6                              74.0788                   243.76100
    ##   blu33.y88.x_male_hypothalamus_bldg blu33.y88.x_male_pituitary_bldg
    ## 1                                179                         474.000
    ## 2                                  1                          34.000
    ## 3                                425                         528.000
    ## 4                                  0                           0.000
    ## 5                                  0                           2.000
    ## 6                                 29                         130.845
    ##   blu36.w16_female_gonad_n9 blu36.w16_female_hypothalamus_n9
    ## 1                        64                         109.0000
    ## 2                        72                           1.0000
    ## 3                       373                         592.0000
    ## 4                        10                           0.0000
    ## 5                         1                           1.0000
    ## 6                        69                          20.2337
    ##   blu36.w16_female_pituitary_n9 blu37.r65.x_male_gonad_n5
    ## 1                         245.0                     5.000
    ## 2                          67.0                     0.000
    ## 3                         739.0                   355.000
    ## 4                           0.0                     0.000
    ## 5                           6.0                     0.000
    ## 6                         171.6                   131.811
    ##   blu37.r65.x_male_hypothalamus_n5 blu37.r65.x_male_pituitary_n5
    ## 1                         291.0000                       223.000
    ## 2                           1.0000                        45.000
    ## 3                         753.0000                       620.000
    ## 4                           0.0000                         0.000
    ## 5                           0.0000                         6.000
    ## 6                          72.7684                       184.411
    ##   blu38.g135.x_female_gonad_bldg blu38.g135.x_female_hypothalamus_bldg
    ## 1                              6                              256.0000
    ## 2                              3                                4.0000
    ## 3                            262                              450.0000
    ## 4                              0                                0.0000
    ## 5                              1                                0.0000
    ## 6                             83                               46.1741
    ##   blu38.g135.x_female_pituitary_bldg blu39.o26.x_female_gonad_inc.d3
    ## 1                                 63                              31
    ## 2                                 25                              10
    ## 3                                239                             360
    ## 4                                  0                              29
    ## 5                                  3                               0
    ## 6                                 88                              15
    ##   blu39.o26.x_female_hypothalamus_inc.d3.NYNO
    ## 1                                    305.0000
    ## 2                                      4.0000
    ## 3                                    516.0000
    ## 4                                      0.0000
    ## 5                                      0.0000
    ## 6                                     28.4879
    ##   blu39.o26.x_female_pituitary_inc.d3.NYNO blu41.y100.x_male_gonad_n5
    ## 1                                263.00000                          2
    ## 2                                 52.00000                          1
    ## 3                                976.00000                        479
    ## 4                                  0.00000                          1
    ## 5                                  2.02492                          2
    ## 6                                316.27000                        104
    ##   blu41.y100.x_male_hypothalamus_n5.NYNO blu41.y100.x_male_pituitary_n5
    ## 1                               149.0000                        261.000
    ## 2                                 4.0000                         31.000
    ## 3                               451.0000                        594.000
    ## 4                                 0.0000                          0.000
    ## 5                                 1.0000                          2.000
    ## 6                                41.3472                        214.119
    ##   blu44.y102_female_gonad_extend blu44.y102_female_hypothalamus_extend
    ## 1                      101.00000                              471.0000
    ## 2                       96.00000                                6.0000
    ## 3                      629.00000                              919.0000
    ## 4                        7.00000                                1.0000
    ## 5                        7.83441                                0.0000
    ## 6                      177.00000                               54.1446
    ##   blu44.y102_female_pituitary_extend blu47.y96.x_female_gonad_inc.d9
    ## 1                            268.000                              73
    ## 2                             29.000                             120
    ## 3                            356.000                             329
    ## 4                              0.000                              43
    ## 5                              1.000                               3
    ## 6                            111.138                              47
    ##   blu47.y96.x_female_hypothalamus_inc.d9
    ## 1                                     93
    ## 2                                      0
    ## 3                                    247
    ## 4                                      0
    ## 5                                      0
    ## 6                                     10
    ##   blu47.y96.x_female_pituitary_inc.d9 blu55.g51_female_gonad_n5
    ## 1                                  85                        19
    ## 2                                  37                        13
    ## 3                                 195                       251
    ## 4                                   0                         7
    ## 5                                   2                         1
    ## 6                                  64                        41
    ##   blu55.g51_female_hypothalamus_n5 blu55.g51_female_pituitary_n5
    ## 1                         461.4700                       62.0000
    ## 2                           3.0000                       11.0000
    ## 3                         772.0000                      177.0000
    ## 4                           0.0000                        0.0000
    ## 5                           5.0000                        2.0000
    ## 6                          73.6446                       42.7425
    ##   blu56.o53_female_gonad_m.inc.d3 blu56.o53_female_hypothalamus_m.inc.d3
    ## 1                              97                               118.0000
    ## 2                             137                                 1.0000
    ## 3                             995                               367.0000
    ## 4                              79                                 0.0000
    ## 5                               1                                 0.0000
    ## 6                              93                                26.1381
    ##   blu56.o53_female_pituitary_m.inc.d3 blu63.g62_female_gonad_m.inc.d9
    ## 1                             260.000                         29.0000
    ## 2                              53.000                         27.0000
    ## 3                             752.000                        254.0000
    ## 4                               0.000                         20.0000
    ## 5                               2.000                          0.0000
    ## 6                             276.256                         58.1971
    ##   blu63.g62_female_hypothalamus_m.inc.d9
    ## 1                               93.00000
    ## 2                                1.00000
    ## 3                              237.00000
    ## 4                                0.00000
    ## 5                                0.00000
    ## 6                                8.12378
    ##   blu63.g62_female_pituitary_m.inc.d9 blu80.r97_female_gonad_m.inc.d8
    ## 1                                  74                              67
    ## 2                                  13                             206
    ## 3                                 214                            1083
    ## 4                                   0                              16
    ## 5                                   1                               4
    ## 6                                  50                             105
    ##   blu80.r97_female_hypothalamus_m.inc.d8
    ## 1                               243.0000
    ## 2                                 0.0000
    ## 3                               588.0000
    ## 4                                 0.0000
    ## 5                                 1.0000
    ## 6                                38.1029
    ##   blu80.r97_female_pituitary_m.inc.d8 blu81.r88_male_gonad_n9
    ## 1                             164.000                 33.0000
    ## 2                              23.000                  5.0000
    ## 3                             610.000               1516.0000
    ## 4                               0.000                  0.0000
    ## 5                               1.000                  6.1045
    ## 6                             155.062                524.9910
    ##   blu81.r88_male_hypothalamus_n9 blu81.r88_male_pituitary_n9
    ## 1                       210.0000                          74
    ## 2                         6.0000                          19
    ## 3                       792.0000                         300
    ## 4                         0.0000                           0
    ## 5                         0.0000                           0
    ## 6                        57.2532                          46
    ##   blu84.x_male_gonad_extend.hatch blu84.x_male_hypothalamus_extend.hatch
    ## 1                           17.00                               225.0000
    ## 2                            0.00                                 3.0000
    ## 3                          184.00                               432.0000
    ## 4                           18.00                                 0.0000
    ## 5                            0.00                                 1.0000
    ## 6                          272.73                                31.2052
    ##   blu84.x_male_pituitary_extend.hatch d.r.blk.s159_female_gonad_m.inc.d9
    ## 1                             556.000                          238.14460
    ## 2                              72.000                          542.00000
    ## 3                             757.000                         1028.00000
    ## 4                               0.000                           22.00000
    ## 5                               1.000                           10.98146
    ## 6                             270.416                          119.96700
    ##   d.r.blk.s159_female_hypothalamus_m.inc.d9
    ## 1                                  401.0000
    ## 2                                    2.0000
    ## 3                                  893.0000
    ## 4                                    0.0000
    ## 5                                    2.0000
    ## 6                                   85.4804
    ##   d.r.blk.s159_female_pituitary_m.inc.d9 d.s008.y.blk_male_gonad_n5
    ## 1                               233.4187                          3
    ## 2                                38.0000                          0
    ## 3                               412.0000                        308
    ## 4                                 0.0000                          1
    ## 5                                 3.0000                          5
    ## 6                               120.8720                        100
    ##   d.s008.y.blk_male_hypothalamus_n5 d.s008.y.blk_male_pituitary_n5
    ## 1                               360                        91.0000
    ## 2                                 9                         3.0000
    ## 3                               777                       173.0000
    ## 4                                 0                         0.0000
    ## 5                                 1                         3.0000
    ## 6                                34                        32.3868
    ##   d.s047.blk.o_male_gonad_n5 d.s047.blk.o_male_hypothalamus_n5
    ## 1                     10.000                          270.0000
    ## 2                      4.000                            5.0000
    ## 3                   1368.000                          595.0000
    ## 4                      1.000                            0.0000
    ## 5                      2.000                            0.0000
    ## 6                    496.665                           35.3459
    ##   d.s047.blk.o_male_pituitary_n5 d.s110.g.blk_male_gonad_m.inc.d3
    ## 1                         85.000                           17.000
    ## 2                         41.000                            4.000
    ## 3                        202.000                         1328.000
    ## 4                          0.000                            2.000
    ## 5                          0.000                            3.000
    ## 6                         62.738                          498.963
    ##   d.s110.g.blk_male_hypothalamus_m.inc.d3
    ## 1                                     128
    ## 2                                       1
    ## 3                                     290
    ## 4                                       0
    ## 5                                       0
    ## 6                                      25
    ##   d.s110.g.blk_male_pituitary_m.inc.d3 d.s112.blk.w_female_gonad_m.inc.d17
    ## 1                              240.000                                  52
    ## 2                               41.000                                  94
    ## 3                              504.000                                 438
    ## 4                                0.000                                  35
    ## 5                                3.000                                  12
    ## 6                              183.331                                  89
    ##   d.s112.blk.w_female_hypothalamus_m.inc.d17
    ## 1                                    71.0000
    ## 2                                     3.0000
    ## 3                                   258.0000
    ## 4                                     0.0000
    ## 5                                     0.0000
    ## 6                                    15.4873
    ##   d.s112.blk.w_female_pituitary_m.inc.d17
    ## 1                                  385.00
    ## 2                                  198.00
    ## 3                                  858.00
    ## 4                                    0.00
    ## 5                                    6.00
    ## 6                                  273.03
    ##   d.s177.blk.r_female_gonad_m.inc.d3
    ## 1                                248
    ## 2                                447
    ## 3                               1095
    ## 4                                118
    ## 5                                  9
    ## 6                                248
    ##   d.s177.blk.r_female_hypothalamus_m.inc.d3
    ## 1                                  257.0000
    ## 2                                    3.0000
    ## 3                                  434.0000
    ## 4                                    0.0000
    ## 5                                    2.0000
    ## 6                                   51.2969
    ##   d.s177.blk.r_female_pituitary_m.inc.d3 g.blk.s004.pk_female_gonad_lay
    ## 1                                    305                            438
    ## 2                                    113                            112
    ## 3                                    463                            731
    ## 4                                      0                            262
    ## 5                                      0                              0
    ## 6                                    183                            120
    ##   g.blk.s004.pk_female_hypothalamus_lay g.blk.s004.pk_female_pituitary_lay
    ## 1                              230.8799                            182.000
    ## 2                                2.0000                             10.000
    ## 3                              778.0000                            521.000
    ## 4                                0.0000                              1.000
    ## 5                                2.0000                              1.000
    ## 6                               71.5386                            279.425
    ##   g.blk.s041.r_male_gonad_m.inc.d3 g.blk.s041.r_male_hypothalamus_m.inc.d3
    ## 1                               12                                216.0000
    ## 2                                2                                  2.0000
    ## 3                              853                                613.0000
    ## 4                                0                                  0.0000
    ## 5                                9                                  5.0000
    ## 6                              442                                 49.1011
    ##   g.blk.s041.r_male_pituitary_m.inc.d3 g.o.y.s037_male_gonad_m.inc.d17
    ## 1                              296.000                               7
    ## 2                               85.000                               2
    ## 3                              448.000                            1061
    ## 4                                0.000                               1
    ## 5                                2.000                               6
    ## 6                              173.263                             389
    ##   g.o.y.s037_male_hypothalamus_m.inc.d17
    ## 1                               210.0000
    ## 2                                 1.0000
    ## 3                               407.0000
    ## 4                                 0.0000
    ## 5                                 1.0000
    ## 6                                31.2125
    ##   g.o.y.s037_male_pituitary_m.inc.d17 g.s.blk.d_male_gonad_n9
    ## 1                              275.00                  21.000
    ## 2                               31.00                  12.000
    ## 3                              737.00                1244.000
    ## 4                                0.00                   1.000
    ## 5                                5.00                   9.000
    ## 6                              242.46                 556.928
    ##   g.s.blk.d_male_hypothalamus_n9 g.s.blk.d_male_pituitary_n9
    ## 1                       128.0000                     205.000
    ## 2                         2.0000                      28.000
    ## 3                       294.0000                     351.000
    ## 4                         0.0000                       0.000
    ## 5                         1.0000                       1.000
    ## 6                        20.1458                     110.423
    ##   g.s.blk.y_male_gonad_lay g.s.blk.y_male_hypothalamus_lay
    ## 1                   14.000                        110.0000
    ## 2                    4.000                          0.0000
    ## 3                  816.000                       1427.0000
    ## 4                    0.000                          0.0000
    ## 5                    0.000                          0.0000
    ## 6                  381.932                         27.1544
    ##   g.s.blk.y_male_pituitary_lay g.s043.pu.blk_male_gonad_lay
    ## 1                      245.000                            2
    ## 2                       88.000                            0
    ## 3                      493.000                           54
    ## 4                        0.000                            0
    ## 5                        0.000                            0
    ## 6                      179.105                           20
    ##   g.s043.pu.blk_male_hypothalamus_lay g.s043.pu.blk_male_pituitary_lay
    ## 1                            190.0000                         243.1996
    ## 2                              1.0000                         117.0000
    ## 3                            489.0000                         694.0000
    ## 4                              0.0000                           1.0000
    ## 5                              1.0000                           2.0000
    ## 6                             41.1324                         163.3460
    ##   g.s075.pk.pu_male_gonad_m.hatch g.s075.pk.pu_male_hypothalamus_m.hatch
    ## 1                          6.0000                               139.0000
    ## 2                          3.0000                                 0.0000
    ## 3                        712.0000                               490.0000
    ## 4                          1.0000                                 0.0000
    ## 5                         14.5647                                 0.0000
    ## 6                        331.9850                                48.2458
    ##   g.s075.pk.pu_male_pituitary_m.hatch g.s076.pk.r_female_gonad_m.hatch
    ## 1                             398.000                               77
    ## 2                             223.000                              131
    ## 3                             604.000                              549
    ## 4                               0.000                               16
    ## 5                               1.000                                9
    ## 6                             202.075                              151
    ##   g.s076.pk.r_female_hypothalamus_m.hatch
    ## 1                                140.0000
    ## 2                                  0.0000
    ## 3                                444.0000
    ## 4                                  0.0000
    ## 5                                  2.0000
    ## 6                                 31.2192
    ##   g.s076.pk.r_female_pituitary_m.hatch g.s078.blk.o_female_gonad_lay
    ## 1                                  235                     248.00000
    ## 2                                   68                     508.00000
    ## 3                                  272                    1629.00000
    ## 4                                    0                     376.00000
    ## 5                                    2                       2.00738
    ## 6                                  109                     176.00000
    ##   g.s078.blk.o_female_hypothalamus_lay g.s078.blk.o_female_pituitary_lay
    ## 1                             148.0000                          330.2731
    ## 2                               4.0000                           42.0000
    ## 3                             479.0000                          631.0000
    ## 4                               0.0000                            0.0000
    ## 5                               1.0000                            3.9136
    ## 6                              24.1861                          173.4140
    ##   g.s111.r.blk_male_gonad_m.inc.d8 g.s111.r.blk_male_hypothalamus_m.inc.d8
    ## 1                            8.000                                420.0000
    ## 2                            6.000                                 14.0000
    ## 3                          772.000                                987.0000
    ## 4                            6.000                                  0.0000
    ## 5                            3.000                                  2.0000
    ## 6                          236.814                                 87.2718
    ##   g.s111.r.blk_male_pituitary_m.inc.d8 g.s179.o.pk_male_gonad_m.inc.d8
    ## 1                             180.0000                              17
    ## 2                              34.0000                               2
    ## 3                             463.0000                            1141
    ## 4                               0.0000                               2
    ## 5                               2.0000                              10
    ## 6                              76.7042                             404
    ##   g.s179.o.pk_male_hypothalamus_m.inc.d8
    ## 1                             261.000000
    ## 2                               7.000000
    ## 3                             533.000000
    ## 4                               0.000000
    ## 5                               0.377751
    ## 6                              57.145900
    ##   g.s179.o.pk_male_pituitary_m.inc.d8 g.s351.pk.w_male_gonad_extend
    ## 1                             370.000                            45
    ## 2                              15.000                            11
    ## 3                            1158.000                          1842
    ## 4                               0.000                             0
    ## 5                               2.000                             8
    ## 6                             207.395                           581
    ##   g.s351.pk.w_male_hypothalamus_extend g.s351.pk.w_male_pituitary_extend
    ## 1                             154.0000                               172
    ## 2                               0.0000                               103
    ## 3                             283.0000                               408
    ## 4                               0.0000                                 0
    ## 5                               1.0000                                 2
    ## 6                              24.7213                               135
    ##   g.x.ATLAS_female_gonad_control g.y.blk.s006_female_gonad_m.inc.d17
    ## 1                            100                            60.00000
    ## 2                            127                           189.00000
    ## 3                            650                           958.00000
    ## 4                            968                             0.00000
    ## 5                              3                             4.00688
    ## 6                             48                           218.00000
    ##   g.y.blk.s006_female_hypothalamus_m.inc.d17
    ## 1                                    66.0000
    ## 2                                     0.0000
    ## 3                                   247.0000
    ## 4                                     0.0000
    ## 5                                     2.0000
    ## 6                                    22.3853
    ##   g.y.blk.s006_female_pituitary_m.inc.d17 g.y.o.s_male_gonad_prolong
    ## 1                                 226.000                     97.000
    ## 2                                 103.000                      2.000
    ## 3                                 584.000                   1068.000
    ## 4                                   0.000                      7.000
    ## 5                                   6.000                      9.000
    ## 6                                 255.062                    395.959
    ##   g.y.o.s_male_hypothalamus_prolong g.y.o.s_male_pituitary_prolong
    ## 1                           342.000                        299.000
    ## 2                             3.000                         19.000
    ## 3                           616.000                        549.000
    ## 4                             0.000                          0.000
    ## 5                             0.000                          3.000
    ## 6                            55.274                        172.128
    ##   g104.w82.x_male_gonad_bldg g104.w82.x_male_hypothalamus_bldg
    ## 1                      7.000                           64.0000
    ## 2                      0.000                            0.0000
    ## 3                    341.000                          199.0000
    ## 4                      0.000                            0.0000
    ## 5                      0.000                            1.0000
    ## 6                    112.909                           10.1845
    ##   g104.w82.x_male_pituitary_bldg g114.w83.x_male_gonad_hatch.NYNO
    ## 1                        213.000                          25.0883
    ## 2                         42.000                           0.0000
    ## 3                        717.000                        1026.0000
    ## 4                          1.000                           0.0000
    ## 5                          0.000                           2.0000
    ## 6                        153.026                         390.9390
    ##   g114.w83.x_male_hypothalamus_hatch g114.w83.x_male_pituitary_hatch.NYNO
    ## 1                           266.0000                             200.0000
    ## 2                             3.0000                              25.0000
    ## 3                           584.0000                             429.0000
    ## 4                             0.0000                               0.0000
    ## 5                             0.0000                               0.0000
    ## 6                            60.6454                              98.7944
    ##   g130.y81.x_male_gonad_inc.d17 g130.y81.x_male_hypothalamus_inc.d17
    ## 1                         1.000                             179.0000
    ## 2                         1.000                               2.0000
    ## 3                       576.000                             525.0000
    ## 4                         2.000                               0.0000
    ## 5                         3.000                               0.0000
    ## 6                       225.923                              28.1043
    ##   g130.y81.x_male_pituitary_inc.d17 g137.r24.w5_male_gonad_m.inc.d8
    ## 1                           246.000                          29.000
    ## 2                            29.000                           4.000
    ## 3                           641.000                        1627.000
    ## 4                             0.000                           2.000
    ## 5                             1.000                           6.000
    ## 6                           206.651                         398.817
    ##   g137.r24.w5_male_hypothalamus_m.inc.d8
    ## 1                                 229.00
    ## 2                                   0.00
    ## 3                                 552.00
    ## 4                                   0.00
    ## 5                                   0.00
    ## 6                                  56.52
    ##   g137.r24.w5_male_pituitary_m.inc.d8 g141.blu27.x_female_gonad_bldg
    ## 1                             447.348                             54
    ## 2                              74.000                             13
    ## 3                             846.000                            212
    ## 4                               0.000                            140
    ## 5                               2.000                              0
    ## 6                             297.533                             61
    ##   g141.blu27.x_female_hypothalamus_bldg g141.blu27.x_female_pituitary_bldg
    ## 1                               445.000                            203.000
    ## 2                                 2.000                             41.000
    ## 3                               631.000                            807.000
    ## 4                                 0.000                              0.000
    ## 5                                 0.000                              1.000
    ## 6                                35.303                            227.845
    ##   g142.r40.x_female_gonad_inc.d17 g142.r40.x_female_hypothalamus_inc.d17
    ## 1                              70                               135.0000
    ## 2                              74                                 0.0000
    ## 3                             326                               611.0000
    ## 4                              11                                 0.0000
    ## 5                               2                                 0.0000
    ## 6                              37                                49.4282
    ##   g142.r40.x_female_pituitary_inc.d17 g143.blu32.x_male_gonad_inc.d17
    ## 1                                  46                           5.000
    ## 2                                   5                           2.000
    ## 3                                 130                         451.000
    ## 4                                   0                           0.000
    ## 5                                   0                           1.000
    ## 6                                  41                         211.936
    ##   g143.blu32.x_male_hypothalamus_inc.d17
    ## 1                               190.0000
    ## 2                                 2.0000
    ## 3                               424.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                14.1323
    ##   g143.blu32.x_male_pituitary_inc.d17 g144.r54.x_female_gonad_m.inc.d3
    ## 1                            122.0000                          106.000
    ## 2                             14.0000                          223.000
    ## 3                            278.0000                          729.000
    ## 4                              0.0000                            6.000
    ## 5                              1.0000                            2.000
    ## 6                             29.6071                           78.943
    ##   g144.r54.x_female_hypothalamus_m.inc.d3
    ## 1                                118.0000
    ## 2                                  1.0000
    ## 3                                255.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 14.0372
    ##   g144.r54.x_female_pituitary_m.inc.d3 g146.blu51_male_gonad_inc.d3
    ## 1                              284.000                            3
    ## 2                               30.000                            0
    ## 3                              531.000                          562
    ## 4                                0.000                            0
    ## 5                                3.000                            1
    ## 6                              144.855                          113
    ##   g146.blu51_male_hypothalamus_inc.d3.NYNO
    ## 1                                 223.0000
    ## 2                                   6.0000
    ## 3                                 778.0000
    ## 4                                   0.0000
    ## 5                                   0.0000
    ## 6                                  25.2943
    ##   g146.blu51_male_pituitary_inc.d3 g17.w108.x_female_gonad_extend
    ## 1                              104                            516
    ## 2                               25                             53
    ## 3                              296                            611
    ## 4                                0                             11
    ## 5                                0                              4
    ## 6                               60                            108
    ##   g17.w108.x_female_hypothalamus_extend g17.w108.x_female_pituitary_extend
    ## 1                              418.0000                                339
    ## 2                                3.0000                                 11
    ## 3                              824.0000                                784
    ## 4                                0.0000                                  0
    ## 5                                4.0000                                  4
    ## 6                               62.6417                                158
    ##   g20.w106.x_male_gonad_inc.d3 g20.w106.x_male_hypothalamus_inc.d3
    ## 1                           10                                  97
    ## 2                            1                                   1
    ## 3                          561                                 137
    ## 4                            2                                   0
    ## 5                            3                                   0
    ## 6                          197                                  17
    ##   g20.w106.x_male_pituitary_inc.d3 g22.blu118_female_gonad_extend
    ## 1                         132.0000                      200.00000
    ## 2                           4.0000                      102.00000
    ## 3                         332.0000                      630.00000
    ## 4                           1.0000                       25.00000
    ## 5                           0.0000                        5.94182
    ## 6                          56.2266                       70.30920
    ##   g22.blu118_female_hypothalamus_extend g22.blu118_female_pituitary_extend
    ## 1                              255.5420                            224.000
    ## 2                                2.0000                             15.000
    ## 3                              533.0000                            401.000
    ## 4                                0.0000                              5.000
    ## 5                                1.0000                              1.000
    ## 6                               48.2148                            107.239
    ##   g3.g119.w20_male_gonad_extend g3.g119.w20_male_hypothalamus_extend
    ## 1                        32.000                                  126
    ## 2                         5.000                                    1
    ## 3                      1740.000                                  337
    ## 4                         2.000                                    0
    ## 5                         7.000                                    0
    ## 6                       592.975                                   30
    ##   g3.g119.w20_male_pituitary_extend g32.blu79_male_gonad_m.inc.d17
    ## 1                           202.000                       33.00000
    ## 2                            57.000                        6.00000
    ## 3                           436.000                      543.00000
    ## 4                             0.000                        0.00000
    ## 5                             2.000                       12.37906
    ## 6                           140.306                      400.00000
    ##   g32.blu79_male_hypothalamus_m.inc.d17 g32.blu79_male_pituitary_m.inc.d17
    ## 1                             202.00000                           344.0043
    ## 2                               4.00000                            28.0000
    ## 3                             736.00000                           605.0000
    ## 4                               0.00000                             0.0000
    ## 5                               2.56112                             3.0000
    ## 6                              39.04830                           148.1250
    ##   g34.x_male_gonad_m.hatch.NYNO g34.x_male_hypothalamus_m.hatch
    ## 1                            14                        283.0000
    ## 2                             3                          1.0000
    ## 3                           737                        536.0000
    ## 4                             1                          0.0000
    ## 5                             9                          0.0000
    ## 6                           302                         64.1078
    ##   g34.x_male_pituitary_m.hatch g38.x_male_gonad_inc.prolong
    ## 1                      275.000                        6.000
    ## 2                       64.000                        1.000
    ## 3                      392.000                      502.000
    ## 4                        0.000                        2.000
    ## 5                        2.000                        3.000
    ## 6                      142.112                      186.925
    ##   g38.x_male_hypothalamus_inc.prolong g38.x_male_pituitary_inc.prolong
    ## 1                             361.000                               51
    ## 2                               4.000                               13
    ## 3                             845.000                              225
    ## 4                               0.000                                0
    ## 5                               0.000                                0
    ## 6                              48.479                               40
    ##   g52.blu58_male_gonad_bldg g52.blu58_male_hypothalamus_bldg
    ## 1                         9                         262.0000
    ## 2                         3                           0.0000
    ## 3                       470                         370.0000
    ## 4                         0                           0.0000
    ## 5                         1                           0.0000
    ## 6                       144                          27.2265
    ##   g52.blu58_male_pituitary_bldg g53.y84_male_gonad_hatch
    ## 1                       142.000                    8.000
    ## 2                        58.000                    2.000
    ## 3                       419.000                  553.000
    ## 4                         0.000                    0.000
    ## 5                         3.000                    0.000
    ## 6                       119.003                  172.985
    ##   g53.y84_male_hypothalamus_hatch g53.y84_male_pituitary_hatch
    ## 1                        72.00000                           52
    ## 2                         4.00000                           15
    ## 3                       231.00000                          178
    ## 4                         0.00000                            0
    ## 5                         0.00000                            0
    ## 6                         8.06621                           47
    ##   g6.w197.x_female_gonad_inc.d3 g6.w197.x_female_hypothalamus_inc.d3
    ## 1                            40                             359.0000
    ## 2                            19                               4.0000
    ## 3                           208                            1078.0000
    ## 4                             8                               0.0000
    ## 5                             1                               0.0000
    ## 6                            31                              89.3317
    ##   g6.w197.x_female_pituitary_inc.d3 g63.blu65_female_gonad_m.inc.d17
    ## 1                           229.000                              153
    ## 2                            82.000                              174
    ## 3                           766.000                              836
    ## 4                             0.000                               37
    ## 5                             3.000                               15
    ## 6                           180.453                              177
    ##   g63.blu65_female_hypothalamus_m.inc.d17
    ## 1                                263.0000
    ## 2                                  3.0000
    ## 3                                551.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 53.5458
    ##   g63.blu65_female_pituitary_m.inc.d17 g73.x_female_gonad_m.inc.d9
    ## 1                              253.000                          32
    ## 2                               27.000                          85
    ## 3                              660.000                         181
    ## 4                                1.000                          30
    ## 5                                0.000                           2
    ## 6                              248.763                          24
    ##   g73.x_female_hypothalamus_m.inc.d9 g73.x_female_pituitary_m.inc.d9
    ## 1                           348.0000                        122.0000
    ## 2                             3.0000                          4.0000
    ## 3                           725.0000                        286.0000
    ## 4                             0.0000                          0.0000
    ## 5                             2.0000                          0.0000
    ## 6                            63.0935                         79.6828
    ##   g75.x_female_gonad_inc.d9 g75.x_female_hypothalamus_inc.d9
    ## 1                        37                         603.0000
    ## 2                        59                           7.0000
    ## 3                       309                        1367.0000
    ## 4                        14                           0.0000
    ## 5                         0                           2.0000
    ## 6                        22                          90.4118
    ##   g75.x_female_pituitary_inc.d9 g8.y197_male_gonad_extend
    ## 1                       222.000                    16.000
    ## 2                        44.000                     7.000
    ## 3                       569.000                  1360.000
    ## 4                         0.000                     1.000
    ## 5                         1.000                     9.000
    ## 6                       134.608                   516.808
    ##   g8.y197_male_hypothalamus_extend g8.y197_male_pituitary_extend
    ## 1                         218.0000                     175.00000
    ## 2                           0.0000                      71.00000
    ## 3                         426.0000                     377.00000
    ## 4                           0.0000                       0.00000
    ## 5                           1.0000                       1.64771
    ## 6                          34.6207                     155.89600
    ##   l.s.o.blk_male_gonad_extend l.s.o.blk_male_hypothalamus_extend
    ## 1                      16.000                           361.0000
    ## 2                       6.000                             6.0000
    ## 3                     725.000                          1008.0000
    ## 4                       5.000                             0.0000
    ## 5                       5.000                             1.0000
    ## 6                     336.986                            69.4175
    ##   l.s.o.blk_male_pituitary_extend l.s.w.d_female_gonad_m.hatch
    ## 1                         512.000                      182.000
    ## 2                         119.000                      229.000
    ## 3                         650.000                      701.000
    ## 4                           0.000                       28.000
    ## 5                           5.000                        1.000
    ## 6                         266.052                      125.902
    ##   l.s.w.d_female_hypothalamus_m.hatch l.s.w.d_female_pituitary_m.hatch
    ## 1                            178.0000                          484.000
    ## 2                              2.0000                           45.000
    ## 3                            397.0000                          636.000
    ## 4                              0.0000                            2.000
    ## 5                              1.0000                            2.000
    ## 6                             23.0574                          279.357
    ##   l.s024.y.g_male_gonad_m.inc.d17 l.s024.y.g_male_hypothalamus_m.inc.d17
    ## 1                          54.000                                85.0000
    ## 2                           0.000                                 0.0000
    ## 3                         802.000                               346.0000
    ## 4                           2.000                                 0.0000
    ## 5                           3.000                                 1.0000
    ## 6                         591.981                                19.2081
    ##   l.s024.y.g_male_pituitary_m.inc.d17 l.s052.pk.r_female_gonad_prolong
    ## 1                             338.000                               91
    ## 2                              24.000                               54
    ## 3                             467.000                              149
    ## 4                               0.000                               23
    ## 5                               0.000                                0
    ## 6                             146.359                               78
    ##   l.s052.pk.r_female_hypothalamus_prolong
    ## 1                                116.0000
    ## 2                                  0.0000
    ## 3                                315.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 33.3344
    ##   l.s052.pk.r_female_pituitary_prolong l.s080.blk.r_male_gonad_prolong
    ## 1                            443.00000                          22.000
    ## 2                            110.00000                           7.000
    ## 3                            755.00000                        1076.000
    ## 4                              0.00000                           1.000
    ## 5                              1.65596                           2.000
    ## 6                            506.38900                         450.953
    ##   l.s080.blk.r_male_hypothalamus_prolong
    ## 1                                    123
    ## 2                                      1
    ## 3                                    278
    ## 4                                      0
    ## 5                                      0
    ## 6                                     16
    ##   l.s080.blk.r_male_pituitary_prolong l.s120.y.blk_female_gonad_bldg
    ## 1                             348.000                        142.000
    ## 2                              42.000                        392.000
    ## 3                             615.000                        895.000
    ## 4                              19.000                         18.000
    ## 5                               5.000                          1.000
    ## 6                             157.277                        116.804
    ##   l.s120.y.blk_female_hypothalamus_bldg l.s120.y.blk_female_pituitary_bldg
    ## 1                                   264                          357.02200
    ## 2                                     4                           39.00000
    ## 3                                   325                          633.00000
    ## 4                                     0                            0.00000
    ## 5                                     0                            4.07053
    ## 6                                    26                          273.84700
    ##   l.s166.o.w_female_gonad_prolong l.s166.o.w_female_hypothalamus_prolong
    ## 1                        160.0000                               135.0000
    ## 2                        132.0000                                 2.0000
    ## 3                        656.0000                               303.0000
    ## 4                         35.0000                                 0.0000
    ## 5                          8.0473                                 1.0000
    ## 6                        139.0000                                21.7751
    ##   l.s166.o.w_female_pituitary_prolong l.s280.g.blk_male_gonad_m.inc.d9
    ## 1                             295.000                               18
    ## 2                              15.000                                4
    ## 3                             402.000                             1466
    ## 4                               1.000                                4
    ## 5                               5.000                                8
    ## 6                             174.643                              440
    ##   l.s280.g.blk_male_hypothalamus_m.inc.d9
    ## 1                                389.0000
    ## 2                                  7.0000
    ## 3                                660.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 57.2579
    ##   l.s280.g.blk_male_pituitary_m.inc.d9 o.d.s009_male_gonad_m.inc.d9
    ## 1                              323.000                      12.0895
    ## 2                               40.000                       2.0000
    ## 3                              401.000                     322.0000
    ## 4                                0.000                       0.0000
    ## 5                                3.000                       4.0000
    ## 6                              137.913                     115.0000
    ##   o.d.s009_male_hypothalamus_m.inc.d9.NYNO
    ## 1                                 147.0000
    ## 2                                   1.0000
    ## 3                                 845.0000
    ## 4                                   0.0000
    ## 5                                   3.0000
    ## 6                                  54.4877
    ##   o.d.s009_male_pituitary_m.inc.d9 o.s.w.r_male_gonad_lay
    ## 1                          94.0000                     10
    ## 2                          30.0000                      4
    ## 3                         248.0000                    651
    ## 4                           0.0000                      2
    ## 5                           0.0000                      2
    ## 6                          65.8021                    267
    ##   o.s.w.r_male_hypothalamus_lay o.s.w.r_male_pituitary_lay
    ## 1                           257                    446.212
    ## 2                             0                    100.000
    ## 3                           603                    611.000
    ## 4                             0                      0.000
    ## 5                             1                      2.000
    ## 6                            32                    178.323
    ##   o.s010.r.blk_male_gonad_m.inc.d3 o.s010.r.blk_male_hypothalamus_m.inc.d3
    ## 1                         17.00000                                 118.000
    ## 2                          3.00000                                   4.000
    ## 3                        789.00000                                 299.000
    ## 4                          2.00000                                   0.000
    ## 5                          3.80319                                   0.000
    ## 6                        452.81300                                  16.341
    ##   o.s010.r.blk_male_pituitary_m.inc.d3 o.s084.w.blk_female_gonad_m.inc.d17
    ## 1                                  234                                 283
    ## 2                                  177                                 448
    ## 3                                  490                                1087
    ## 4                                    0                                  45
    ## 5                                    3                                   7
    ## 6                                  135                                 260
    ##   o.s084.w.blk_female_hypothalamus_m.inc.d17
    ## 1                                   133.0000
    ## 2                                     0.0000
    ## 3                                   312.0000
    ## 4                                     0.0000
    ## 5                                     0.0000
    ## 6                                    25.0745
    ##   o.s084.w.blk_female_pituitary_m.inc.d17 o114.blu9_male_gonad_inc.prolong
    ## 1                                 441.000                               13
    ## 2                                  82.000                                3
    ## 3                                 521.000                              591
    ## 4                                   0.000                                1
    ## 5                                   0.000                                0
    ## 6                                 198.724                              219
    ##   o114.blu9_male_hypothalamus_inc.prolong
    ## 1                                252.0000
    ## 2                                  1.0000
    ## 3                                732.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 44.3434
    ##   o114.blu9_male_pituitary_inc.prolong o152.o120.w42_male_gonad_n5
    ## 1                              98.0000                    80.00000
    ## 2                              84.0000                     4.00000
    ## 3                             246.0000                  1134.00000
    ## 4                               0.0000                     2.00000
    ## 5                               0.0000                     3.79367
    ## 6                              93.7873                   692.93800
    ##   o152.o120.w42_male_hypothalamus_n5 o152.o120.w42_male_pituitary_n5
    ## 1                           241.0000                          323.00
    ## 2                             3.0000                          104.00
    ## 3                           508.0000                          750.00
    ## 4                             0.0000                            0.00
    ## 5                             0.0000                            0.00
    ## 6                            54.4265                          280.78
    ##   o156.w80.x_female_gonad_inc.d3 o156.w80.x_female_hypothalamus_inc.d3
    ## 1                             30                              372.0000
    ## 2                             53                                5.0000
    ## 3                            243                              567.0000
    ## 4                              7                                0.0000
    ## 5                              0                                1.0000
    ## 6                             18                               60.3712
    ##   o156.w80.x_female_pituitary_inc.d3 o165.w122.x_female_gonad_inc.d3.NYNO
    ## 1                            275.000                                  121
    ## 2                            159.000                                  201
    ## 3                            491.000                                  780
    ## 4                              0.000                                  186
    ## 5                              1.000                                    0
    ## 6                            183.222                                   57
    ##   o165.w122.x_female_hypothalamus_inc.d3
    ## 1                               333.0000
    ## 2                                 1.0000
    ## 3                               804.0000
    ## 4                                 0.0000
    ## 5                                 2.0000
    ## 6                                71.2884
    ##   o165.w122.x_female_pituitary_inc.d3.NYNO
    ## 1                                  348.000
    ## 2                                  113.000
    ## 3                                  599.000
    ## 4                                    0.000
    ## 5                                    2.000
    ## 6                                  120.699
    ##   o169.r28.x_female_gonad_m.inc.d17
    ## 1                          53.00000
    ## 2                         135.00000
    ## 3                         979.00000
    ## 4                          14.00000
    ## 5                           5.01482
    ## 6                         153.00000
    ##   o169.r28.x_female_hypothalamus_m.inc.d17
    ## 1                                  93.0000
    ## 2                                   0.0000
    ## 3                                 279.0000
    ## 4                                   0.0000
    ## 5                                   0.0000
    ## 6                                  18.1848
    ##   o169.r28.x_female_pituitary_m.inc.d17 o172.w115.x_female_gonad_hatch
    ## 1                             314.04870                             52
    ## 2                              48.00000                             21
    ## 3                             716.00000                            285
    ## 4                               0.00000                             64
    ## 5                               4.03228                              0
    ## 6                             259.22500                             18
    ##   o172.w115.x_female_hypothalamus_hatch
    ## 1                              587.0000
    ## 2                                4.0000
    ## 3                              735.0000
    ## 4                                1.0000
    ## 5                                2.0000
    ## 6                               56.5556
    ##   o172.w115.x_female_pituitary_hatch.NYNO o173.w179.x_female_gonad_inc.d3
    ## 1                                139.0000                              39
    ## 2                                 18.0000                             103
    ## 3                                339.0000                             339
    ## 4                                  0.0000                               1
    ## 5                                  0.0000                               1
    ## 6                                 99.8199                              14
    ##   o173.w179.x_female_hypothalamus_inc.d3
    ## 1                               166.0000
    ## 2                                 1.0000
    ## 3                               517.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                28.0671
    ##   o173.w179.x_female_pituitary_inc.d3 o3.x_male_gonad_m.n2
    ## 1                             99.0000               15.000
    ## 2                             19.0000                1.000
    ## 3                            189.0000              410.000
    ## 4                              0.0000                1.000
    ## 5                              1.0000                1.000
    ## 6                             48.3061              170.904
    ##   o3.x_male_hypothalamus_m.n2 o3.x_male_pituitary_m.n2
    ## 1                    217.0000                  83.0000
    ## 2                      7.0000                  27.0000
    ## 3                    489.0000                 204.0000
    ## 4                      0.0000                   0.0000
    ## 5                      1.0000                   0.0000
    ## 6                     38.0962                  64.9791
    ##   o35.r51.x_female_gonad_inc.d17 o35.r51.x_female_hypothalamus_inc.d17
    ## 1                            158                              435.0000
    ## 2                            170                                2.0000
    ## 3                            771                              743.0000
    ## 4                            260                                0.0000
    ## 5                              3                                1.0000
    ## 6                            123                               39.5935
    ##   o35.r51.x_female_pituitary_inc.d17 o36.r62.x_female_gonad_m.inc.d9
    ## 1                                109                              49
    ## 2                                 30                              84
    ## 3                                286                             285
    ## 4                                  0                               2
    ## 5                                  1                               2
    ## 6                                 88                              23
    ##   o36.r62.x_female_hypothalamus_m.inc.d9
    ## 1                               278.0000
    ## 2                                 6.0000
    ## 3                               549.0000
    ## 4                                 0.0000
    ## 5                                 2.0000
    ## 6                                26.0506
    ##   o36.r62.x_female_pituitary_m.inc.d9 o38.blu29.x_female_gonad_bldg
    ## 1                                  93                            61
    ## 2                                  10                           173
    ## 3                                 241                           497
    ## 4                                   0                           135
    ## 5                                   0                             3
    ## 6                                  36                            82
    ##   o38.blu29.x_female_hypothalamus_bldg o38.blu29.x_female_pituitary_bldg
    ## 1                             205.0000                          143.0000
    ## 2                               2.0000                           19.0000
    ## 3                             582.0000                          263.0000
    ## 4                               0.0000                            0.0000
    ## 5                               0.0000                            3.0000
    ## 6                              45.2665                           85.1016
    ##   o39.y77.x_male_gonad_hatch o39.y77.x_male_hypothalamus_hatch
    ## 1                   14.00000                           463.791
    ## 2                    0.00000                             5.000
    ## 3                  623.00000                           631.000
    ## 4                    0.00000                             0.000
    ## 5                    6.75478                             0.000
    ## 6                  153.76000                            39.713
    ##   o39.y77.x_male_pituitary_hatch o44.blu26.x_male_gonad_hatch
    ## 1                             93                            7
    ## 2                             61                            1
    ## 3                            226                          365
    ## 4                              0                            1
    ## 5                              0                            1
    ## 6                             92                          126
    ##   o44.blu26.x_male_hypothalamus_hatch o44.blu26.x_male_pituitary_hatch
    ## 1                            334.0000                          311.000
    ## 2                              5.0000                           25.000
    ## 3                            811.0000                          417.000
    ## 4                              0.0000                            0.000
    ## 5                              1.0000                            0.000
    ## 6                             66.3776                          170.714
    ##   o45.g128.x_female_gonad_m.inc.d9.NYNO
    ## 1                                    51
    ## 2                                   133
    ## 3                                   599
    ## 4                                     5
    ## 5                                     9
    ## 6                                    85
    ##   o45.g128.x_female_hypothalamus_m.inc.d9
    ## 1                                229.0000
    ## 2                                  2.0000
    ## 3                                651.0000
    ## 4                                  0.0000
    ## 5                                  3.0000
    ## 6                                 35.2079
    ##   o45.g128.x_female_pituitary_m.inc.d9 o48.r197.x_male_gonad_inc.d3
    ## 1                              312.000                            7
    ## 2                              113.000                            2
    ## 3                              965.000                          281
    ## 4                                3.000                            3
    ## 5                                4.000                            3
    ## 6                              253.029                          110
    ##   o48.r197.x_male_hypothalamus_inc.d3 o48.r197.x_male_pituitary_inc.d3
    ## 1                            238.0000                            189.0
    ## 2                              8.0000                            249.0
    ## 3                            565.0000                            767.0
    ## 4                              0.0000                              0.0
    ## 5                              0.0000                              2.0
    ## 6                             48.0802                            176.2
    ##   o49.x_male_gonad_inc.d9 o49.x_male_hypothalamus_inc.d9
    ## 1                  10.000                            142
    ## 2                   3.000                              3
    ## 3                 300.000                            231
    ## 4                   0.000                              0
    ## 5                   3.000                              0
    ## 6                 197.982                             17
    ##   o49.x_male_pituitary_inc.d9 o52.blu53_female_gonad_inc.d17
    ## 1                     88.1613                            113
    ## 2                     52.0000                            225
    ## 3                    218.0000                            703
    ## 4                      0.0000                            142
    ## 5                      2.0000                              4
    ## 6                     65.2782                             60
    ##   o52.blu53_female_hypothalamus_inc.d17 o52.blu53_female_pituitary_inc.d17
    ## 1                              323.0000                           174.0000
    ## 2                                1.0000                            12.0000
    ## 3                              644.0000                           356.0000
    ## 4                                0.0000                             0.0000
    ## 5                                5.0000                             3.0000
    ## 6                               37.8189                            69.9528
    ##   o57.g59_male_gonad_inc.d9 o57.g59_male_hypothalamus_inc.d9
    ## 1                     6.000                         126.0000
    ## 2                     1.000                           2.0000
    ## 3                   565.000                         205.0000
    ## 4                     0.000                           0.0000
    ## 5                     2.000                           1.0000
    ## 6                   218.655                          14.4038
    ##   o57.g59_male_pituitary_inc.d9 o59.blu64_male_gonad_m.inc.d17
    ## 1                           104                       18.00000
    ## 2                            28                        4.00000
    ## 3                           259                      960.00000
    ## 4                             0                        1.00000
    ## 5                             1                        8.66485
    ## 6                            70                      515.00000
    ##   o59.blu64_male_hypothalamus_m.inc.d17 o59.blu64_male_pituitary_m.inc.d17
    ## 1                              339.4690                           267.2821
    ## 2                                0.0000                           110.0000
    ## 3                              714.0000                           557.0000
    ## 4                                0.0000                             0.0000
    ## 5                                2.0000                             0.0000
    ## 6                               48.3326                           206.2510
    ##   o73.x_female_gonad_inc.d9 o73.x_female_hypothalamus_inc.d9
    ## 1                        35                              129
    ## 2                        73                                1
    ## 3                       366                              359
    ## 4                        63                                0
    ## 5                         1                                1
    ## 6                        47                               24
    ##   o73.x_female_pituitary_inc.d9 p.g.blk.s040_female_gonad_m.inc.d8
    ## 1                           148                                176
    ## 2                            15                                182
    ## 3                           171                                611
    ## 4                             0                                110
    ## 5                             0                                  2
    ## 6                            46                                 42
    ##   p.g.blk.s040_female_hypothalamus_m.inc.d8
    ## 1                                   221.000
    ## 2                                     5.000
    ## 3                                   591.000
    ## 4                                     0.000
    ## 5                                     1.000
    ## 6                                    36.132
    ##   p.g.blk.s040_female_pituitary_m.inc.d8 pk.s.d.g_male_gonad_prolong
    ## 1                                292.544                         108
    ## 2                                 26.000                         215
    ## 3                                474.000                         831
    ## 4                                  0.000                           8
    ## 5                                  0.000                          13
    ## 6                                130.283                         234
    ##   pk.s.d.g_male_hypothalamus_prolong pk.s.d.g_male_pituitary_prolong
    ## 1                            97.0000                         224.000
    ## 2                             1.0000                         108.000
    ## 3                           364.0000                         473.000
    ## 4                             0.0000                           0.000
    ## 5                             1.0000                           1.000
    ## 6                            34.1803                         172.298
    ##   pk.s011.o.y_female_gonad_m.inc.d17
    ## 1                           283.1745
    ## 2                           597.0000
    ## 3                          1838.0000
    ## 4                            30.0000
    ## 5                             4.0000
    ## 6                           188.0000
    ##   pk.s011.o.y_female_hypothalamus_m.inc.d17
    ## 1                                       133
    ## 2                                         5
    ## 3                                       393
    ## 4                                         0
    ## 5                                         0
    ## 6                                        15
    ##   pk.s011.o.y_female_pituitary_m.inc.d17 pk.s054.d.g_female_gonad_m.hatch
    ## 1                                232.000                               71
    ## 2                                 51.000                              192
    ## 3                                416.000                              463
    ## 4                                  0.000                               14
    ## 5                                  8.000                               10
    ## 6                                136.775                              102
    ##   pk.s054.d.g_female_hypothalamus_m.hatch
    ## 1                                217.0000
    ## 2                                  1.0000
    ## 3                                603.0000
    ## 4                                  0.0000
    ## 5                                  1.0000
    ## 6                                 39.1933
    ##   pk.s054.d.g_female_pituitary_m.hatch pk.s055.d.l_female_gonad_m.inc.d8
    ## 1                              559.000                          64.00000
    ## 2                               57.000                         126.00000
    ## 3                              889.000                         925.00000
    ## 4                                0.000                         423.00000
    ## 5                                5.000                           4.02273
    ## 6                              290.393                          99.00000
    ##   pk.s055.d.l_female_hypothalamus_m.inc.d8
    ## 1                                 239.0000
    ## 2                                   1.0000
    ## 3                                 547.0000
    ## 4                                   0.0000
    ## 5                                   1.0000
    ## 6                                  56.5185
    ##   pk.s055.d.l_female_pituitary_m.inc.d8 pk.s238.blk.w_male_gonad_lay
    ## 1                               273.000                     69.00000
    ## 2                                38.000                     14.00000
    ## 3                               768.000                   2576.00000
    ## 4                                 0.000                      6.00000
    ## 5                                 2.000                     16.43527
    ## 6                               347.064                    890.20200
    ##   pk.s238.blk.w_male_hypothalamus_lay pk.s238.blk.w_male_pituitary_lay
    ## 1                            206.0000                        305.00000
    ## 2                              3.0000                         55.00000
    ## 3                            720.0000                        540.00000
    ## 4                              0.0000                          1.00000
    ## 5                              1.0000                          1.06018
    ## 6                             44.4001                        128.69600
    ##   pk.w.s141.o_male_gonad_lay pk.w.s141.o_male_hypothalamus_lay
    ## 1                          9                          348.0000
    ## 2                          7                            3.0000
    ## 3                       2365                          689.0000
    ## 4                          5                            1.0000
    ## 5                          5                            0.0000
    ## 6                        481                           47.6348
    ##   pk.w.s141.o_male_pituitary_lay pu.blk.s102.y_male_gonad_m.inc.d8
    ## 1                        345.000                            22.000
    ## 2                        106.000                             3.000
    ## 3                        723.000                          1722.000
    ## 4                          0.000                             4.000
    ## 5                          4.000                             9.000
    ## 6                        244.974                           648.811
    ##   pu.blk.s102.y_male_hypothalamus_m.inc.d8
    ## 1                                 178.0000
    ## 2                                   0.0000
    ## 3                                 457.0000
    ## 4                                   0.0000
    ## 5                                   1.0000
    ## 6                                  26.0551
    ##   pu.blk.s102.y_male_pituitary_m.inc.d8 pu.s.o.r_male_hypothalamus_m.hatch
    ## 1                               204.000                           213.0000
    ## 2                                46.000                             3.0000
    ## 3                               503.000                           554.0000
    ## 4                                 1.000                             0.0000
    ## 5                                 3.000                             2.0000
    ## 6                               192.009                            50.2942
    ##   pu.s.o.r_male_pituitary_m.hatch r.r.x.ATLAS.R2XR_female_gonad_control
    ## 1                         359.000                                    16
    ## 2                          25.000                                    22
    ## 3                         502.000                                   259
    ## 4                           0.000                                    17
    ## 5                           0.000                                     0
    ## 6                         211.369                                    23
    ##   r.r.x.ATLAS.R2XR_female_hypothalamus_control
    ## 1                                      48.0000
    ## 2                                       0.0000
    ## 3                                     163.0000
    ## 4                                       0.0000
    ## 5                                       0.0000
    ## 6                                      19.2357
    ##   r.r.x.ATLAS.R2XR_female_pituitary_control
    ## 1                                   88.0000
    ## 2                                   14.0000
    ## 3                                  228.0000
    ## 4                                    0.0000
    ## 5                                    1.0000
    ## 6                                   48.2428
    ##   r.r.x.ATLAS_female_gonad_control r.r.x.ATLAS_female_hypothalamus_control
    ## 1                              134                                 48.0000
    ## 2                              154                                  0.0000
    ## 3                              773                                163.0000
    ## 4                              806                                  0.0000
    ## 5                                0                                  0.0000
    ## 6                               28                                 19.2353
    ##   r.r.x.ATLAS_female_pituitary_control r.s005.pk.blk_male_gonad_lay
    ## 1                                   23                           44
    ## 2                                    5                            2
    ## 3                                   71                          938
    ## 4                                    0                            0
    ## 5                                    0                            8
    ## 6                                   12                          692
    ##   r.s005.pk.blk_male_hypothalamus_lay r.s005.pk.blk_male_pituitary_lay
    ## 1                            288.0000                          146.000
    ## 2                              5.0000                           12.000
    ## 3                            459.0000                          636.000
    ## 4                              0.0000                            0.000
    ## 5                              1.0000                            3.000
    ## 6                             32.2842                          143.958
    ##   r.s035.y.blk_female_gonad_m.inc.d3
    ## 1                          219.00000
    ## 2                          150.00000
    ## 3                          852.00000
    ## 4                           79.00000
    ## 5                            5.02095
    ## 6                          158.00000
    ##   r.s035.y.blk_female_hypothalamus_m.inc.d3
    ## 1                                  321.0000
    ## 2                                    3.0000
    ## 3                                  669.0000
    ## 4                                    0.0000
    ## 5                                    1.0000
    ## 6                                   54.1273
    ##   r.s035.y.blk_female_pituitary_m.inc.d3 r.s056.g.o_female_gonad_bldg
    ## 1                                497.000                           15
    ## 2                                 21.000                           49
    ## 3                                530.000                          265
    ## 4                                  0.000                           90
    ## 5                                  0.000                            0
    ## 6                                176.551                            8
    ##   r.s056.g.o_female_hypothalamus_bldg r.s056.g.o_female_pituitary_bldg
    ## 1                            154.0000                          94.0000
    ## 2                              0.0000                           7.0000
    ## 3                            237.0000                         388.0000
    ## 4                              0.0000                           0.0000
    ## 5                              0.0000                           0.0000
    ## 6                             16.1403                          55.2826
    ##   r.s057.g.pk_male_gonad_extend r.s057.g.pk_male_hypothalamus_extend
    ## 1                            26                             152.0000
    ## 2                             3                               1.0000
    ## 3                           610                             576.0000
    ## 4                             0                               0.0000
    ## 5                             2                               1.0000
    ## 6                           291                              35.1602
    ##   r.s057.g.pk_male_pituitary_extend r.s058.d.l_male_gonad_m.inc.d8
    ## 1                           265.000                         22.000
    ## 2                           122.000                          2.000
    ## 3                           451.000                        579.000
    ## 4                             0.000                          3.000
    ## 5                             3.000                          2.000
    ## 6                           169.844                        388.958
    ##   r.s058.d.l_male_hypothalamus_m.inc.d8 r.s058.d.l_male_pituitary_m.inc.d8
    ## 1                              231.0000                            254.000
    ## 2                                0.0000                             38.000
    ## 3                              639.0000                            556.000
    ## 4                                0.0000                              0.000
    ## 5                                0.0000                              0.000
    ## 6                               24.2426                            205.905
    ##   r.s059.d.o_male_gonad_bldg r.s059.d.o_male_hypothalamus_bldg
    ## 1                      5.000                               364
    ## 2                      1.000                                 5
    ## 3                    471.000                               637
    ## 4                      0.000                                 0
    ## 5                      2.000                                 0
    ## 6                    192.961                                49
    ##   r.s059.d.o_male_pituitary_bldg r.s086.l.blk_male_gonad_extend
    ## 1                       164.0000                         12.000
    ## 2                         9.0000                          1.000
    ## 3                       245.0000                       1486.000
    ## 4                         0.0000                          3.000
    ## 5                         0.0000                          3.000
    ## 6                        70.5407                        489.739
    ##   r.s086.l.blk_male_hypothalamus_extend r.s086.l.blk_male_pituitary_extend
    ## 1                              218.0000                             160.00
    ## 2                                4.0000                              41.00
    ## 3                              680.0000                             344.00
    ## 4                                0.0000                               0.00
    ## 5                                2.0000                               1.00
    ## 6                               61.1766                             207.88
    ##   r.s116.blk.pu_male_gonad_lay r.s116.blk.pu_male_hypothalamus_lay
    ## 1                            9                            229.0000
    ## 2                            2                              0.0000
    ## 3                          805                            558.0000
    ## 4                            0                              1.0000
    ## 5                            4                              0.0000
    ## 6                          304                             40.1807
    ##   r.s116.blk.pu_male_pituitary_lay r.s131.o.d_female_gonad_extend
    ## 1                          441.000                        195.000
    ## 2                           46.000                        305.000
    ## 3                          780.000                        686.000
    ## 4                            0.000                         31.000
    ## 5                            3.000                          2.000
    ## 6                          205.162                        109.928
    ##   r.s131.o.d_female_hypothalamus_extend r.s131.o.d_female_pituitary_extend
    ## 1                              496.0000                            524.000
    ## 2                                3.0000                             22.000
    ## 3                              949.0000                            982.000
    ## 4                                0.0000                              0.000
    ## 5                                6.0000                              6.000
    ## 6                               79.8977                            343.488
    ##   r.s171.l.w_female_gonad_n9 r.s171.l.w_female_hypothalamus_n9
    ## 1                    213.000                          280.0000
    ## 2                    214.000                            0.0000
    ## 3                   1060.000                          709.0000
    ## 4                     11.000                            0.0000
    ## 5                      8.000                            0.0000
    ## 6                    326.992                           55.3649
    ##   r.s171.l.w_female_pituitary_n9 r.s172.l.y_male_gonad_extend
    ## 1                        208.000                       46.000
    ## 2                         51.000                        8.000
    ## 3                        496.000                      620.000
    ## 4                          0.000                        0.000
    ## 5                          1.000                        4.000
    ## 6                        155.265                      246.953
    ##   r.s172.l.y_male_hypothalamus_extend r.s172.l.y_male_pituitary_extend
    ## 1                            346.0000                          167.000
    ## 2                              7.0000                           84.000
    ## 3                            867.0000                          409.000
    ## 4                              0.0000                            0.000
    ## 5                              1.0000                            2.000
    ## 6                             74.0479                          134.466
    ##   r.y.s007.blk_male_gonad_n9 r.y.s007.blk_male_hypothalamus_n9
    ## 1                     14.000                               102
    ## 2                      2.000                                 1
    ## 3                    771.000                               229
    ## 4                      1.000                                 0
    ## 5                      5.000                                 0
    ## 6                    252.869                                 9
    ##   r.y.s007.blk_male_pituitary_n9 r176.blu54_male_gonad_inc.d17
    ## 1                       123.0000                         3.000
    ## 2                        15.0000                         2.000
    ## 3                       218.0000                       886.000
    ## 4                         0.0000                         0.000
    ## 5                         0.0000                         1.000
    ## 6                        90.2129                       236.929
    ##   r176.blu54_male_hypothalamus_inc.d17 r176.blu54_male_pituitary_inc.d17
    ## 1                            118.00000                           156.000
    ## 2                              0.00000                            14.000
    ## 3                            209.00000                           338.000
    ## 4                              0.00000                             0.000
    ## 5                              0.00000                             0.000
    ## 6                              7.13669                           104.282
    ##   r183.o22_female_gonad_hatch r183.o22_female_hypothalamus_hatch
    ## 1                          60                           345.6826
    ## 2                          16                             5.0000
    ## 3                         196                           670.0000
    ## 4                          11                             0.0000
    ## 5                           0                             0.0000
    ## 6                          23                            46.5800
    ##   r183.o22_female_pituitary_hatch r190.o43.x_male_gonad_lay
    ## 1                              76                     3.000
    ## 2                              10                     3.000
    ## 3                             303                   906.000
    ## 4                               0                     2.000
    ## 5                               0                     8.000
    ## 6                              70                   234.871
    ##   r190.o43.x_male_hypothalamus_lay r190.o43.x_male_pituitary_lay
    ## 1                          273.000                      130.0000
    ## 2                            7.000                       51.0000
    ## 3                          696.000                      528.0000
    ## 4                            0.000                        0.0000
    ## 5                            2.000                        2.0000
    ## 6                           38.167                       69.0822
    ##   r194.x_female_gonad_prolong r194.x_female_hypothalamus_prolong
    ## 1                         163                             92.000
    ## 2                         268                              0.000
    ## 3                         989                            283.000
    ## 4                          19                              0.000
    ## 5                          12                              0.000
    ## 6                         105                             22.054
    ##   r194.x_female_pituitary_prolong r195.x_male_gonad_n9
    ## 1                        332.3145                 51.0
    ## 2                         21.0000                  7.0
    ## 3                        541.0000               1209.0
    ## 4                          0.0000                  1.0
    ## 5                          2.0000                 11.0
    ## 6                        268.7010                539.5
    ##   r195.x_male_hypothalamus_n9 r195.x_male_pituitary_n9
    ## 1                    283.1267                  244.000
    ## 2                      6.0000                   69.000
    ## 3                    427.0000                  476.000
    ## 4                      0.0000                    0.000
    ## 5                      0.0000                    4.000
    ## 6                     24.0000                  135.442
    ##   r196.x_male_gonad_m.inc.d9 r196.x_male_hypothalamus_m.inc.d9
    ## 1                          9                          251.0000
    ## 2                          1                            2.0000
    ## 3                        479                          683.0000
    ## 4                          1                            0.0000
    ## 5                          1                            1.0000
    ## 6                        249                           59.3666
    ##   r196.x_male_pituitary_m.inc.d9 r27.w111.blu125_female_gonad_inc.d3
    ## 1                       126.0000                            102.2163
    ## 2                        64.0000                            149.0000
    ## 3                       249.0000                            503.0000
    ## 4                         0.0000                              9.0000
    ## 5                         0.0000                              2.0000
    ## 6                        46.3643                             69.9099
    ##   r27.w111.blu125_female_hypothalamus_inc.d3
    ## 1                                   255.0000
    ## 2                                     5.0000
    ## 3                                   640.0000
    ## 4                                     0.0000
    ## 5                                     1.0000
    ## 6                                    69.4284
    ##   r27.w111.blu125_female_pituitary_inc.d3 r30.w112.r46_female_gonad_inc.d9
    ## 1                                  86.000                               44
    ## 2                                  21.000                               78
    ## 3                                 828.000                              317
    ## 4                                   0.000                                7
    ## 5                                   0.000                                1
    ## 6                                 123.205                               43
    ##   r30.w112.r46_female_hypothalamus_inc.d9
    ## 1                                 258.000
    ## 2                                   9.000
    ## 3                                 863.000
    ## 4                                   0.000
    ## 5                                   0.000
    ## 6                                  81.642
    ##   r30.w112.r46_female_pituitary_inc.d9 r36.w184.x_female_gonad_inc.d9
    ## 1                              370.000                             33
    ## 2                               11.000                             25
    ## 3                              611.000                            186
    ## 4                                0.000                             13
    ## 5                                7.000                              0
    ## 6                              242.469                             16
    ##   r36.w184.x_female_hypothalamus_inc.d9 r36.w184.x_female_pituitary_inc.d9
    ## 1                              335.5032                            64.0000
    ## 2                                3.0000                             8.0000
    ## 3                              634.0000                           145.0000
    ## 4                                0.0000                             0.0000
    ## 5                                2.0000                             0.0000
    ## 6                               50.0884                            29.1238
    ##   r37.w100.x_male_gonad_n9 r37.w100.x_male_hypothalamus_n9
    ## 1                        7                       241.00000
    ## 2                        3                         4.00000
    ## 3                      570                       805.00000
    ## 4                        0                         0.00000
    ## 5                        1                         7.12279
    ## 6                      234                        69.09420
    ##   r37.w100.x_male_pituitary_n9 r4.x_male_gonad_m.inc.d8
    ## 1                      56.0000                    20.00
    ## 2                      50.0000                     6.00
    ## 3                     176.0000                   568.00
    ## 4                       0.0000                     3.00
    ## 5                       1.0000                     6.00
    ## 6                      60.1219                   264.74
    ##   r4.x_male_hypothalamus_m.inc.d8 r4.x_male_pituitary_m.inc.d8
    ## 1                        162.0000                     361.4383
    ## 2                          2.0000                      61.0000
    ## 3                        483.0000                    1063.0000
    ## 4                          0.0000                       2.0000
    ## 5                          0.0000                       1.0000
    ## 6                         36.3035                     374.2990
    ##   r41.w99.x_male_gonad_hatch r41.w99.x_male_hypothalamus_hatch
    ## 1                          7                          253.0000
    ## 2                          4                            2.0000
    ## 3                        493                          555.0000
    ## 4                          1                            0.0000
    ## 5                          2                            0.0000
    ## 6                        142                           53.2058
    ##   r41.w99.x_male_pituitary_hatch r45.X_male_gonad_inc.d9
    ## 1                        216.000                   9.000
    ## 2                         53.000                   0.000
    ## 3                        490.000                 481.000
    ## 4                          0.000                   0.000
    ## 5                          0.000                   0.000
    ## 6                        199.285                 281.956
    ##   r45.X_male_pituitary_inc.d9 r45.x_male_hypothalamus_inc.d9
    ## 1                    143.0000                       145.4928
    ## 2                     28.0000                         1.0000
    ## 3                    339.0000                       389.0000
    ## 4                      0.0000                         0.0000
    ## 5                      1.0000                         0.0000
    ## 6                     61.5609                         8.0000
    ##   r49.w189.x_female_gonad_inc.d17 r49.w189.x_female_hypothalamus_inc.d17
    ## 1                              71                                    150
    ## 2                             132                                      4
    ## 3                             306                                    362
    ## 4                              58                                      0
    ## 5                               3                                      1
    ## 6                              55                                     33
    ##   r49.w189.x_female_pituitary_inc.d17 r6.x_female_gonad_control.NYNO
    ## 1                             193.000                             96
    ## 2                              15.000                            101
    ## 3                             489.000                            603
    ## 4                               0.000                             37
    ## 5                               5.000                              1
    ## 6                             258.941                             34
    ##   r6.x_female_hypothalamus_control.NYNO r6.x_female_pituitary_control
    ## 1                                   226                            15
    ## 2                                     4                             3
    ## 3                                   474                            41
    ## 4                                     0                             0
    ## 5                                     0                             1
    ## 6                                    24                            14
    ##   r72.y83.x_male_gonad_hatch r72.y83.x_male_hypothalamus_hatch
    ## 1                          2                          488.4638
    ## 2                          0                            3.0000
    ## 3                        521                          589.0000
    ## 4                         19                            0.0000
    ## 5                          3                            2.0000
    ## 6                        181                           70.9725
    ##   r72.y83.x_male_pituitary_hatch r73.g127.x_female_gonad_inc.d3
    ## 1                        132.000                            252
    ## 2                         31.000                            286
    ## 3                        666.000                           2155
    ## 4                          1.000                            177
    ## 5                          2.000                              2
    ## 6                        145.901                            202
    ##   r73.g127.x_female_hypothalamus_inc.d3 r73.g127.x_female_pituitary_inc.d3
    ## 1                              229.0000                                103
    ## 2                                2.0000                                 55
    ## 3                              563.0000                                226
    ## 4                                0.0000                                  0
    ## 5                                3.0000                                  1
    ## 6                               56.2033                                 86
    ##   r81.x_male_gonad_inc.prolong r81.x_male_hypothalamus_inc.prolong
    ## 1                        3.000                            240.0000
    ## 2                        0.000                              5.0000
    ## 3                      437.000                            474.0000
    ## 4                        0.000                              0.0000
    ## 5                        0.000                              0.0000
    ## 6                      115.944                             21.1408
    ##   r81.x_male_pituitary_inc.prolong r83.g45_female_gonad_bldg
    ## 1                         112.0000                        22
    ## 2                           9.0000                         5
    ## 3                         194.0000                       245
    ## 4                           0.0000                         1
    ## 5                           1.0000                         0
    ## 6                          76.7789                        36
    ##   r83.g45_female_hypothalamus_bldg r83.g45_female_pituitary_bldg
    ## 1                              101                       288.000
    ## 2                                0                        42.000
    ## 3                              225                       498.000
    ## 4                                0                         0.000
    ## 5                                0                         4.000
    ## 6                               15                       128.952
    ##   r84.x_male_gonad_m.inc.d9 r84.x_male_hypothalamus_m.inc.d9
    ## 1                     5.000                              102
    ## 2                     2.000                                2
    ## 3                   511.000                              309
    ## 4                     0.000                                0
    ## 5                     1.000                                0
    ## 6                   124.916                               22
    ##   r84.x_male_pituitary_m.inc.d9 r85.g39_male_gonad_m.inc.d9.NYNO
    ## 1                            81                         14.00000
    ## 2                            68                          3.00000
    ## 3                           205                        844.00000
    ## 4                             0                          1.00000
    ## 5                             1                          1.02191
    ## 6                            49                        335.00000
    ##   r85.g39_male_hypothalamus_m.inc.d9 r85.g39_male_pituitary_m.inc.d9.NYNO
    ## 1                            72.0000                              283.320
    ## 2                             0.0000                              138.000
    ## 3                           247.0000                             1034.000
    ## 4                             0.0000                                0.000
    ## 5                             1.0000                                1.000
    ## 6                            11.0597                              290.658
    ##   r90.x_male_gonad_inc.prolong r90.x_male_hypothalamus_inc.prolong
    ## 1                            3                            299.0000
    ## 2                            0                              1.0000
    ## 3                          399                            687.0000
    ## 4                            0                              0.0000
    ## 5                            2                              1.0000
    ## 6                          148                             49.0618
    ##   r90.x_male_pituitary_inc.prolong r95.blu99_female_gonad_n9
    ## 1                         153.0000                       198
    ## 2                          26.0000                       400
    ## 3                         322.0000                       777
    ## 4                           0.0000                        32
    ## 5                           0.0000                         7
    ## 6                          77.5533                        80
    ##   r95.blu99_female_hypothalamus_n9 r95.blu99_female_pituitary_n9
    ## 1                        370.00000                           116
    ## 2                          0.00000                            12
    ## 3                        817.00000                           273
    ## 4                          0.00000                             0
    ## 5                          3.79307                             2
    ## 6                         45.30540                            61
    ##   s.o.pk_female_gonad_lay s.o.pk_female_hypothalamus_lay
    ## 1                     501                       289.0000
    ## 2                     338                         2.0000
    ## 3                     587                       521.0000
    ## 4                     132                         0.0000
    ## 5                       2                         0.0000
    ## 6                      57                        28.3354
    ##   s.o.pk_female_pituitary_lay s.pu.pk_female_gonad_prolong
    ## 1                     282.000                     137.8081
    ## 2                      42.000                     176.0000
    ## 3                     615.000                     602.0000
    ## 4                       6.000                      17.0000
    ## 5                       4.000                       1.0000
    ## 6                     147.841                      93.0000
    ##   s.pu.pk_female_hypothalamus_prolong s.pu.pk_female_pituitary_prolong
    ## 1                            239.0000                          249.000
    ## 2                              3.0000                           65.000
    ## 3                            710.0000                          659.000
    ## 4                              0.0000                            1.000
    ## 5                              1.0000                            0.000
    ## 6                             65.7126                          243.232
    ##   s.pu148.blk.r_male_gonad_bldg s.pu148.blk.r_male_hypothalamus_bldg
    ## 1                         4.000                             138.0000
    ## 2                         1.000                               2.0000
    ## 3                       357.000                             241.0000
    ## 4                         0.000                               0.0000
    ## 5                         0.000                               0.0000
    ## 6                       107.891                              16.0806
    ##   s.pu148.blk.r_male_pituitary_bldg s.x.ATLAS_female_gonad_control
    ## 1                           70.0000                             68
    ## 2                           47.0000                            105
    ## 3                          255.0000                            441
    ## 4                            0.0000                            546
    ## 5                            1.0000                              9
    ## 6                           50.4172                             25
    ##   s.x.ATLAS_female_hypothalamus_control s.x.ATLAS_female_pituitary_control
    ## 1                                    74                                 15
    ## 2                                     2                                  0
    ## 3                                   277                                110
    ## 4                                     0                                  0
    ## 5                                     0                                  1
    ## 6                                    26                                 13
    ##   s002.x_female_gonad_m.inc.d8 s002.x_female_hypothalamus_m.inc.d8
    ## 1                     214.1315                            172.5488
    ## 2                     392.0000                              1.0000
    ## 3                    1263.0000                            433.0000
    ## 4                     105.0000                              0.0000
    ## 5                       0.0000                              0.0000
    ## 6                     287.9320                             22.0661
    ##   s002.x_female_pituitary_m.inc.d8 s038.g.d.blk_female_gonad_m.inc.d17
    ## 1                          205.000                           189.45850
    ## 2                           20.000                           415.00000
    ## 3                          591.000                          1451.00000
    ## 4                            0.000                            50.00000
    ## 5                            7.000                            11.92053
    ## 6                          201.708                           150.00000
    ##   s038.g.d.blk_female_hypothalamus_m.inc.d17
    ## 1                                        137
    ## 2                                          1
    ## 3                                        287
    ## 4                                          0
    ## 5                                          1
    ## 6                                         25
    ##   s038.g.d.blk_female_pituitary_m.inc.d17 s044.blk.d.r_male_gonad_m.inc.d8
    ## 1                                 454.281                           26.000
    ## 2                                 200.000                            5.000
    ## 3                                 711.000                         1552.000
    ## 4                                   0.000                            8.000
    ## 5                                   3.000                            5.000
    ## 6                                 239.919                          753.463
    ##   s044.blk.d.r_male_hypothalamus_m.inc.d8
    ## 1                                115.0000
    ## 2                                  2.0000
    ## 3                                212.0000
    ## 4                                  0.0000
    ## 5                                  2.0000
    ## 6                                 21.2082
    ##   s044.blk.d.r_male_pituitary_m.inc.d8 s062.d.blk.g_male_gonad_m.hatch
    ## 1                            259.00000                           32.00
    ## 2                             69.00000                            4.00
    ## 3                            586.00000                          941.00
    ## 4                              0.00000                            0.00
    ## 5                              5.92508                            4.00
    ## 6                            131.74600                          299.97
    ##   s062.d.blk.g_male_hypothalamus_m.hatch
    ## 1                               164.0000
    ## 2                                 1.0000
    ## 3                               546.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                11.1595
    ##   s062.d.blk.g_male_pituitary_m.hatch s063.d.blk.l_female_gonad_bldg
    ## 1                             650.000                       67.00000
    ## 2                             151.000                       71.00000
    ## 3                            1055.000                      303.00000
    ## 4                               0.000                      102.00000
    ## 5                               7.000                        2.01395
    ## 6                             321.659                      129.82500
    ##   s063.d.blk.l_female_hypothalamus_bldg s063.d.blk.l_female_pituitary_bldg
    ## 1                              150.0000                            69.0000
    ## 2                                3.0000                            18.0000
    ## 3                              290.0000                           272.0000
    ## 4                                0.0000                             0.0000
    ## 5                                0.0000                             0.0000
    ## 6                               16.0577                            80.9741
    ##   s064.g.blk.pu_male_gonad_m.inc.d17
    ## 1                             36.000
    ## 2                              4.000
    ## 3                           1334.000
    ## 4                              6.000
    ## 5                              5.000
    ## 6                            479.961
    ##   s064.g.blk.pu_male_hypothalamus_m.inc.d17
    ## 1                                  118.0000
    ## 2                                    0.0000
    ## 3                                  402.0000
    ## 4                                    0.0000
    ## 5                                    1.0000
    ## 6                                   13.0267
    ##   s064.g.blk.pu_male_pituitary_m.inc.d17 s065.l.d.o_male_gonad_bldg
    ## 1                                375.369                         11
    ## 2                                 38.000                          4
    ## 3                                611.000                        462
    ## 4                                  0.000                          0
    ## 5                                  3.000                          7
    ## 6                                216.276                        539
    ##   s065.l.d.o_male_hypothalamus_bldg s065.l.d.o_male_pituitary_bldg
    ## 1                               167                        351.000
    ## 2                                 2                         35.000
    ## 3                               450                        967.000
    ## 4                                 0                          0.000
    ## 5                                 1                          9.000
    ## 6                                17                        177.632
    ##   s066.l.d.r_male_gonad_bldg s066.l.d.r_male_hypothalamus_bldg
    ## 1                         12                               102
    ## 2                          0                                 3
    ## 3                        515                               193
    ## 4                          0                                 0
    ## 5                          2                                 0
    ## 6                        172                                16
    ##   s066.l.d.r_male_pituitary_bldg s067.o.l.y_male_gonad_m.inc.d8
    ## 1                        84.0000                         19.000
    ## 2                        10.0000                          0.000
    ## 3                       255.0000                       1193.000
    ## 4                         0.0000                          1.000
    ## 5                         0.0000                          4.000
    ## 6                        63.0944                        422.807
    ##   s067.o.l.y_male_hypothalamus_m.inc.d8 s067.o.l.y_male_pituitary_m.inc.d8
    ## 1                              186.0000                            371.000
    ## 2                                2.0000                             50.000
    ## 3                              372.0000                            620.000
    ## 4                                0.0000                              0.000
    ## 5                                1.0000                             13.000
    ## 6                               37.5402                            209.668
    ##   s069.pk.pu.g_female_gonad_m.inc.d3
    ## 1                           644.3560
    ## 2                           212.0000
    ## 3                           914.0000
    ## 4                            44.0000
    ## 5                            20.0000
    ## 6                            62.8809
    ##   s069.pk.pu.g_female_hypothalamus_m.inc.d3
    ## 1                                   193.000
    ## 2                                     1.000
    ## 3                                   451.000
    ## 4                                     0.000
    ## 5                                     0.000
    ## 6                                    25.141
    ##   s069.pk.pu.g_female_pituitary_m.inc.d3 s071.pu.g.pk_male_gonad_m.inc.d3
    ## 1                                213.000                          19.0324
    ## 2                                 69.000                           3.0000
    ## 3                                597.000                         942.0000
    ## 4                                  0.000                           2.0000
    ## 5                                  4.000                           1.0000
    ## 6                                201.018                         443.0000
    ##   s071.pu.g.pk_male_hypothalamus_m.inc.d3
    ## 1                                113.0000
    ## 2                                  0.0000
    ## 3                                324.0000
    ## 4                                  0.0000
    ## 5                                  4.0000
    ## 6                                 17.0477
    ##   s071.pu.g.pk_male_pituitary_m.inc.d3 s089.blk.pk.pu_female_gonad_extend
    ## 1                           421.000000                                244
    ## 2                             5.000000                                222
    ## 3                           622.000000                                557
    ## 4                             0.000000                                 21
    ## 5                             0.997325                                  1
    ## 6                           166.419000                                 77
    ##   s089.blk.pk.pu_female_hypothalamus_extend
    ## 1                                  588.0000
    ## 2                                    2.0000
    ## 3                                 1094.0000
    ## 4                                    0.0000
    ## 5                                    1.0000
    ## 6                                   78.3128
    ##   s089.blk.pk.pu_female_pituitary_extend s090.blk.pk.w_male_gonad_m.hatch
    ## 1                              328.00000                               14
    ## 2                               30.00000                                2
    ## 3                              996.00000                              915
    ## 4                                0.00000                               10
    ## 5                                5.09787                                5
    ## 6                              321.05000                              368
    ##   s090.blk.pk.w_male_hypothalamus_m.hatch
    ## 1                                  260.00
    ## 2                                    4.00
    ## 3                                  787.00
    ## 4                                    0.00
    ## 5                                    3.00
    ## 6                                   65.82
    ##   s090.blk.pk.w_male_pituitary_m.hatch s091.blk.r.g_female_gonad_m.inc.d8
    ## 1                              207.000                             103.00
    ## 2                               99.000                             406.00
    ## 3                              458.000                            3222.47
    ## 4                                0.000                             299.00
    ## 5                                0.000                               6.00
    ## 6                              225.737                             111.00
    ##   s091.blk.r.g_female_hypothalamus_m.inc.d8
    ## 1                                   486.000
    ## 2                                     5.000
    ## 3                                   883.000
    ## 4                                     7.000
    ## 5                                     1.000
    ## 6                                    48.186
    ##   s091.blk.r.g_female_pituitary_m.inc.d8 s092.blk.r.o_female_gonad_bldg
    ## 1                                379.000                        46.0000
    ## 2                                 30.000                       140.0000
    ## 3                                834.000                       518.0000
    ## 4                                  0.000                        32.0000
    ## 5                                  3.000                         3.0000
    ## 6                                256.668                        76.8905
    ##   s092.blk.r.o_female_hypothalamus_bldg s092.blk.r.o_female_pituitary_bldg
    ## 1                                   514                           103.0000
    ## 2                                     6                             8.0000
    ## 3                                   643                           218.0000
    ## 4                                     0                             0.0000
    ## 5                                     3                             0.0000
    ## 6                                    55                            62.0045
    ##   s093.blk.y.d_female_gonad_m.inc.d3
    ## 1                            267.000
    ## 2                            624.000
    ## 3                           1529.000
    ## 4                            586.000
    ## 5                              7.000
    ## 6                            242.925
    ##   s093.blk.y.d_female_hypothalamus_m.inc.d3
    ## 1                                   87.0000
    ## 2                                    1.0000
    ## 3                                  292.0000
    ## 4                                    0.0000
    ## 5                                    0.0000
    ## 6                                   15.0684
    ##   s093.blk.y.d_female_pituitary_m.inc.d3 s095.g.blk.o_female_gonad_lay
    ## 1                                356.000                     133.00000
    ## 2                                 28.000                     282.00000
    ## 3                                640.000                     572.00000
    ## 4                                  0.000                      71.00000
    ## 5                                  2.000                       5.08058
    ## 6                                189.069                      56.00000
    ##   s095.g.blk.o_female_hypothalamus_lay s095.g.blk.o_female_pituitary_lay
    ## 1                             506.0000                           238.000
    ## 2                               5.0000                            54.000
    ## 3                             985.0000                           699.000
    ## 4                               1.0000                             0.000
    ## 5                               1.0000                             3.000
    ## 6                              86.7824                           150.369
    ##   s096.g.blk.pk_male_gonad_m.inc.d17
    ## 1                             80.000
    ## 2                              2.000
    ## 3                            958.000
    ## 4                              1.000
    ## 5                              3.000
    ## 6                            616.645
    ##   s096.g.blk.pk_male_hypothalamus_m.inc.d17
    ## 1                                 87.000000
    ## 2                                  2.000000
    ## 3                                302.000000
    ## 4                                  0.000000
    ## 5                                  0.574588
    ## 6                                 19.034600
    ##   s096.g.blk.pk_male_pituitary_m.inc.d17
    ## 1                              234.00000
    ## 2                               19.00000
    ## 3                              476.00000
    ## 4                                0.00000
    ## 5                                3.04165
    ## 6                              175.44100
    ##   s100.l.pk.blk_female_gonad_m.inc.d17
    ## 1                             89.00000
    ## 2                            106.00000
    ## 3                            529.00000
    ## 4                             55.00000
    ## 5                              4.01436
    ## 6                             98.00000
    ##   s100.l.pk.blk_female_hypothalamus_m.inc.d17
    ## 1                                         106
    ## 2                                           0
    ## 3                                         295
    ## 4                                           0
    ## 5                                           2
    ## 6                                          17
    ##   s100.l.pk.blk_female_pituitary_m.inc.d17
    ## 1                                 377.0000
    ## 2                                  90.0000
    ## 3                                 698.0000
    ## 4                                   0.0000
    ## 5                                   3.0488
    ## 6                                 322.4510
    ##   s103.y.pk.blk_female_gonad_m.inc.d3
    ## 1                            146.0000
    ## 2                             38.0000
    ## 3                            509.0000
    ## 4                             30.0000
    ## 5                              5.0000
    ## 6                             49.8351
    ##   s103.y.pk.blk_female_hypothalamus_m.inc.d3
    ## 1                                   112.0000
    ## 2                                     1.0000
    ## 3                                   285.0000
    ## 4                                     0.0000
    ## 5                                     1.0000
    ## 6                                    17.0905
    ##   s103.y.pk.blk_female_pituitary_m.inc.d3 s136.d.w.o_female_gonad_lay
    ## 1                                 162.000                         279
    ## 2                                  10.000                          78
    ## 3                                 395.000                         381
    ## 4                                   2.000                         150
    ## 5                                   2.000                           2
    ## 6                                  99.456                          57
    ##   s136.d.w.o_female_hypothalamus_lay s136.d.w.o_female_pituitary_lay
    ## 1                            307.000                        277.5261
    ## 2                              7.000                         14.0000
    ## 3                            537.000                        515.0000
    ## 4                              0.000                          0.0000
    ## 5                              0.000                          4.0000
    ## 6                             60.513                        227.4950
    ##   s137.g.blk.y_female_gonad_extend s137.g.blk.y_female_hypothalamus_extend
    ## 1                              268                                493.0000
    ## 2                              101                                  3.0000
    ## 3                              490                                735.0000
    ## 4                                3                                  0.0000
    ## 5                                1                                  0.0000
    ## 6                               89                                 43.1931
    ##   s137.g.blk.y_female_pituitary_extend s139.l.blk.w_male_gonad_m.inc.d3
    ## 1                             404.0919                            35.00
    ## 2                              29.0000                            11.00
    ## 3                             904.0000                          1245.00
    ## 4                               0.0000                             0.00
    ## 5                               4.0000                             5.00
    ## 6                             190.7940                           455.87
    ##   s139.l.blk.w_male_hypothalamus_m.inc.d3
    ## 1                                 87.0000
    ## 2                                  1.0000
    ## 3                                198.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 14.0403
    ##   s139.l.blk.w_male_pituitary_m.inc.d3 s142.o.pk.pu_female_gonad_lay
    ## 1                              370.000                            48
    ## 2                               28.000                             1
    ## 3                              401.000                           616
    ## 4                                0.000                             7
    ## 5                                2.000                             0
    ## 6                              110.345                            23
    ##   s142.o.pk.pu_female_hypothalamus_lay s142.o.pk.pu_female_pituitary_lay
    ## 1                              445.438                          298.0856
    ## 2                                2.000                           32.0000
    ## 3                              844.000                          401.0000
    ## 4                                0.000                            0.0000
    ## 5                                2.000                            2.0000
    ## 6                               67.000                          164.9370
    ##   s150.w.g.blk_male_gonad_lay s150.w.g.blk_male_hypothalamus_lay
    ## 1                   14.000000                           285.0000
    ## 2                    3.000000                             2.0000
    ## 3                  647.000000                           832.0000
    ## 4                    0.000000                             0.0000
    ## 5                    0.457772                             0.0000
    ## 6                  259.000000                            71.1904
    ##   s150.w.g.blk_male_pituitary_lay s175.blk.pu.pk_male_gonad_m.inc.d3
    ## 1                         363.000                             16.000
    ## 2                         101.000                              1.000
    ## 3                         661.000                           1029.000
    ## 4                           0.000                              0.000
    ## 5                           0.000                              8.000
    ## 6                         165.657                            362.863
    ##   s175.blk.pu.pk_male_hypothalamus_m.inc.d3
    ## 1                                   98.0000
    ## 2                                    1.0000
    ## 3                                  230.0000
    ## 4                                    0.0000
    ## 5                                    0.0000
    ## 6                                   16.0568
    ##   s175.blk.pu.pk_male_pituitary_m.inc.d3 s176.blk.pu.r_female_gonad_lay
    ## 1                                248.000                        66.1214
    ## 2                                 56.000                        77.0000
    ## 3                                465.000                       417.0000
    ## 4                                  0.000                        23.0000
    ## 5                                  0.000                         0.0000
    ## 6                                212.525                        32.0000
    ##   s176.blk.pu.r_female_hypothalamus_lay s176.blk.pu.r_female_pituitary_lay
    ## 1                               353.442                            387.322
    ## 2                                 1.000                             45.000
    ## 3                               522.000                            819.000
    ## 4                                 0.000                              0.000
    ## 5                                 1.000                              1.000
    ## 6                                44.398                            237.888
    ##   s186.l.o.pk_female_gonad_extend s186.l.o.pk_female_hypothalamus_extend
    ## 1                          125.00                               207.0000
    ## 2                           47.00                                 0.0000
    ## 3                          378.00                               561.0000
    ## 4                           18.00                                 0.0000
    ## 5                            1.00                                 1.0000
    ## 6                          123.85                                26.1634
    ##   s186.l.o.pk_female_pituitary_extend s187.l.o.r_male_gonad_n9
    ## 1                             457.000                  41.1015
    ## 2                              62.000                   6.0000
    ## 3                             699.000                1173.0000
    ## 4                               3.000                   0.0000
    ## 5                               4.000                   6.0000
    ## 6                             184.162                 630.9350
    ##   s187.l.o.r_male_hypothalamus_n9 s187.l.o.r_male_pituitary_n9
    ## 1                        218.0000                     246.0000
    ## 2                          1.0000                      51.0000
    ## 3                        728.0000                     476.0000
    ## 4                          0.0000                       0.0000
    ## 5                          0.0000                       1.0000
    ## 6                         18.0298                      90.5408
    ##   s243.blk.pk.r_male_gonad_lay s243.blk.pk.r_male_hypothalamus_lay
    ## 1                       18.000                            312.0000
    ## 2                        1.000                              5.0000
    ## 3                     1037.000                            670.0000
    ## 4                        0.000                              0.0000
    ## 5                       10.000                              2.0000
    ## 6                      571.952                             54.2584
    ##   s243.blk.pk.r_male_pituitary_lay s333.y.blk.pk_female_gonad_m.hatch
    ## 1                          335.000                                 73
    ## 2                           82.000                                138
    ## 3                         1083.000                                412
    ## 4                            0.000                                  3
    ## 5                            0.000                                  6
    ## 6                          263.717                                 63
    ##   s333.y.blk.pk_female_hypothalamus_m.hatch
    ## 1                                  252.0000
    ## 2                                    6.0000
    ## 3                                  594.0000
    ## 4                                    0.0000
    ## 5                                    3.0000
    ## 6                                   41.2738
    ##   s333.y.blk.pk_female_pituitary_m.hatch w191.r1_female_gonad_control
    ## 1                               194.2892                           11
    ## 2                                32.0000                           23
    ## 3                               338.0000                           72
    ## 4                                 0.0000                            8
    ## 5                                 1.0000                            1
    ## 6                               120.9890                            7
    ##   w191.r1_female_hypothalamus_control w191.r1_female_pituitary_control
    ## 1                                  37                               37
    ## 2                                   0                                7
    ## 3                                  35                               43
    ## 4                                   0                                0
    ## 5                                   0                                1
    ## 6                                   2                               14
    ##   w34.x_male_gonad_inc.d9 w34.x_male_hypothalamus_inc.d9
    ## 1                   6.000                        215.000
    ## 2                   3.000                          2.000
    ## 3                 475.000                        671.000
    ## 4                   0.000                          0.000
    ## 5                   3.000                          1.000
    ## 6                 216.756                         40.241
    ##   w34.x_male_pituitary_inc.d9 x.blk.blk.ATLAS_male_gonad_control
    ## 1                     50.0000                            24.0747
    ## 2                     31.0000                             2.0000
    ## 3                    235.0000                           932.0000
    ## 4                      2.0000                             0.0000
    ## 5                      1.0000                             6.0000
    ## 6                     54.6582                           212.9450
    ##   x.blk.blk.ATLAS_male_hypothalamus_control
    ## 1                                 147.00000
    ## 2                                   2.00000
    ## 3                                 129.00000
    ## 4                                   0.00000
    ## 5                                   0.00000
    ## 6                                   8.05413
    ##   x.blk.blk.ATLAS_male_pituitary_control x.blk16_male_gonad_n9.NYNO
    ## 1                               30.00000                    43.3122
    ## 2                               25.00000                     4.0000
    ## 3                               54.00000                   868.0000
    ## 4                                0.00000                     0.0000
    ## 5                                1.00000                     6.0000
    ## 6                                6.29407                   287.9820
    ##   x.blk16_male_hypothalamus_n9.NYNO x.blk16_male_pituitary_n9
    ## 1                          157.5762                   222.000
    ## 2                            3.0000                    75.000
    ## 3                          560.0000                   719.000
    ## 4                            0.0000                     0.000
    ## 5                            0.0000                     4.000
    ## 6                           33.0000                   221.994
    ##   x.blk20_female_gonad_prolong x.blk20_female_hypothalamus_prolong
    ## 1                       21.000                            102.0000
    ## 2                        2.000                              2.0000
    ## 3                      756.000                            292.0000
    ## 4                       30.000                              0.0000
    ## 5                        2.000                              1.0000
    ## 6                      425.959                             45.1597
    ##   x.blk20_female_pituitary_prolong x.blk23_male_gonad_m.n2
    ## 1                          165.000                       5
    ## 2                           78.000                       0
    ## 3                          508.000                     371
    ## 4                            0.000                       0
    ## 5                            3.000                       0
    ## 6                          270.421                     134
    ##   x.blk23_male_hypothalamus_m.n2 x.blk23_male_pituitary_m.n2
    ## 1                       131.0000                    135.0000
    ## 2                         1.0000                    106.0000
    ## 3                       267.0000                    297.0000
    ## 4                         0.0000                      0.0000
    ## 5                         0.0000                      4.0000
    ## 6                        21.0632                     65.7529
    ##   x.blu.o.ATLAS_male_pituitary_control x.blu101.w43_female_gonad_inc.d9
    ## 1                                  119                               27
    ## 2                                   36                               78
    ## 3                                  290                              228
    ## 4                                    0                               28
    ## 5                                    2                                1
    ## 6                                   32                               27
    ##   x.blu101.w43_female_hypothalamus_inc.d9
    ## 1                                 81.0000
    ## 2                                  1.0000
    ## 3                                163.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 15.0698
    ##   x.blu101.w43_female_pituitary_inc.d9 x.blu102.w105_female_gonad_inc.d3
    ## 1                              292.385                                28
    ## 2                               26.000                                61
    ## 3                              736.000                               387
    ## 4                                0.000                                52
    ## 5                                0.000                                 1
    ## 6                              162.969                                20
    ##   x.blu102.w105_female_hypothalamus_inc.d3
    ## 1                                 331.6843
    ## 2                                   6.0000
    ## 3                                 666.0000
    ## 4                                   0.0000
    ## 5                                   1.0000
    ## 6                                  55.3019
    ##   x.blu102.w105_female_pituitary_inc.d3
    ## 1                               184.000
    ## 2                                22.000
    ## 3                               892.000
    ## 4                                 0.000
    ## 5                                 5.000
    ## 6                               160.137
    ##   x.blu106.o153_male_gonad_inc.d9.NYNO
    ## 1                               18.000
    ## 2                                3.000
    ## 3                              954.000
    ## 4                                1.000
    ## 5                                7.000
    ## 6                              274.692
    ##   x.blu106.o153_male_hypothalamus_inc.d9
    ## 1                               223.3344
    ## 2                                 4.0000
    ## 3                               514.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                32.1358
    ##   x.blu106.o153_male_pituitary_inc.d9.NYNO x.blu109.w121_female_gonad_n5
    ## 1                                 319.0000                            70
    ## 2                                  18.0000                            33
    ## 3                                 552.0000                           279
    ## 4                                   0.0000                            21
    ## 5                                   0.0000                             0
    ## 6                                  95.3834                            43
    ##   x.blu109.w121_female_hypothalamus_n5 x.blu109.w121_female_pituitary_n5
    ## 1                             297.0000                           112.000
    ## 2                               7.0000                            25.000
    ## 3                             536.0000                           230.000
    ## 4                               0.0000                             0.000
    ## 5                               3.0000                             1.000
    ## 6                              29.2823                            65.915
    ##   x.blu116.w107_female_gonad_inc.d17
    ## 1                                161
    ## 2                                 29
    ## 3                                258
    ## 4                                 13
    ## 5                                  0
    ## 6                                  7
    ##   x.blu116.w107_female_hypothalamus_inc.d17
    ## 1                                  315.0000
    ## 2                                    4.0000
    ## 3                                  529.0000
    ## 4                                    0.0000
    ## 5                                    0.0000
    ## 6                                   42.2818
    ##   x.blu116.w107_female_pituitary_inc.d17.NYNO
    ## 1                                     160.000
    ## 2                                      38.000
    ## 3                                     487.000
    ## 4                                       0.000
    ## 5                                       0.000
    ## 6                                     135.431
    ##   x.blu117.w89_male_gonad_inc.d17 x.blu117.w89_male_hypothalamus_inc.d17
    ## 1                               4                               263.0000
    ## 2                               1                                 0.0000
    ## 3                             441                               346.0000
    ## 4                               2                                 0.0000
    ## 5                               1                                 3.0000
    ## 6                             134                                23.1035
    ##   x.blu117.w89_male_pituitary_inc.d17 x.blu122.r66_female_gonad_inc.d9
    ## 1                                  87                               80
    ## 2                                  12                               93
    ## 3                                 190                              758
    ## 4                                   0                                8
    ## 5                                   0                                1
    ## 6                                  52                               50
    ##   x.blu122.r66_female_hypothalamus_inc.d9
    ## 1                                285.0000
    ## 2                                  1.0000
    ## 3                                753.0000
    ## 4                                  0.0000
    ## 5                                  1.0000
    ## 6                                 28.0435
    ##   x.blu122.r66_female_pituitary_inc.d9 x.blu23.w14_male_gonad_n9
    ## 1                             158.0000                     11.00
    ## 2                               6.0000                      3.00
    ## 3                             274.0000                    351.00
    ## 4                               0.0000                      1.00
    ## 5                               0.0000                      0.00
    ## 6                              75.5002                    123.42
    ##   x.blu23.w14_male_hypothalamus_n9 x.blu23.w14_male_pituitary_n9
    ## 1                         115.0000                       198.000
    ## 2                           0.0000                        59.000
    ## 3                         255.0000                       623.000
    ## 4                           0.0000                         0.000
    ## 5                           1.0000                         3.000
    ## 6                          11.2852                       190.169
    ##   x.blu25_female_gonad_m.inc.d9 x.blu25_female_hypothalamus_m.inc.d9
    ## 1                            12                                  331
    ## 2                             6                                    1
    ## 3                           247                                  604
    ## 4                             0                                    0
    ## 5                             1                                    1
    ## 6                            30                                   60
    ##   x.blu25_female_pituitary_m.inc.d9 x.blu30_male_gonad_n5
    ## 1                           45.0000                 20.00
    ## 2                           29.0000                  0.00
    ## 3                          235.0000               1211.00
    ## 4                            0.0000                  0.00
    ## 5                            0.0000                 11.00
    ## 6                           98.5696                623.71
    ##   x.blu30_male_hypothalamus_n5 x.blu30_male_pituitary_n5
    ## 1                     115.0000                  134.2518
    ## 2                       1.0000                   36.0000
    ## 3                     349.0000                  423.0000
    ## 4                       0.0000                    0.0000
    ## 5                       1.0000                    3.0000
    ## 6                      35.2691                  262.6560
    ##   x.blu42.o28_male_gonad_inc.d3 x.blu42.o28_male_hypothalamus_inc.d3.NYNO
    ## 1                             3                                       233
    ## 2                             1                                         0
    ## 3                           464                                       534
    ## 4                             0                                         0
    ## 5                             2                                         0
    ## 6                           136                                        43
    ##   x.blu42.o28_male_pituitary_inc.d3 x.blu43.g132_female_gonad_n9
    ## 1                           302.000                           29
    ## 2                            55.000                           36
    ## 3                           918.000                          341
    ## 4                             0.000                           16
    ## 5                             1.000                            2
    ## 6                           238.091                           81
    ##   x.blu43.g132_female_hypothalamus_n9 x.blu43.g132_female_pituitary_n9
    ## 1                            231.0000                               56
    ## 2                              3.0000                               27
    ## 3                            711.0000                              169
    ## 4                              0.0000                                0
    ## 5                              2.0000                                4
    ## 6                             66.8668                               56
    ##   x.blu57_male_gonad_prolong x.blu57_male_hypothalamus_prolong
    ## 1                     17.000                          292.0000
    ## 2                      2.000                            3.0000
    ## 3                   1114.000                          254.0000
    ## 4                      0.000                            0.0000
    ## 5                      6.000                            0.0000
    ## 6                    437.963                           40.5562
    ##   x.blu57_male_pituitary_prolong x.blu6.y80_female_gonad_lay
    ## 1                      277.00000                     80.0000
    ## 2                       27.00000                    148.0000
    ## 3                      379.00000                    603.0000
    ## 4                        0.00000                      2.0000
    ## 5                        2.80812                      6.0000
    ## 6                      193.26500                     46.3754
    ##   x.blu6.y80_female_hypothalamus_lay x.blu6.y80_female_pituitary_lay
    ## 1                           406.0000                              89
    ## 2                             7.0000                              19
    ## 3                           646.0000                             321
    ## 4                             0.0000                               0
    ## 5                             2.0000                               0
    ## 6                            49.4926                              37
    ##   x.blu69_male_gonad_extend x.blu69_male_hypothalamus_extend
    ## 1                     9.000                          233.000
    ## 2                     3.000                            4.000
    ## 3                  1038.000                          561.000
    ## 4                     0.000                            0.000
    ## 5                     9.000                            0.000
    ## 6                   358.874                           43.096
    ##   x.blu69_male_pituitary_extend x.blu73_male_gonad_extend
    ## 1                       295.000                    14.000
    ## 2                        80.000                     4.000
    ## 3                       642.000                  1967.000
    ## 4                         0.000                     4.000
    ## 5                         0.000                    10.000
    ## 6                       270.535                   534.797
    ##   x.blu73_male_hypothalamus_extend x.blu73_male_pituitary_extend
    ## 1                         218.0000                       210.000
    ## 2                           4.0000                        52.000
    ## 3                         622.0000                       362.000
    ## 4                           0.0000                         0.000
    ## 5                           4.0000                         5.000
    ## 6                          63.3916                       148.857
    ##   x.blu76_male_gonad_m.hatch x.blu76_male_hypothalamus_m.hatch
    ## 1                   71.00000                          238.4507
    ## 2                    6.00000                            1.0000
    ## 3                 1851.00000                          629.0000
    ## 4                    1.00000                            0.0000
    ## 5                    6.38926                            0.0000
    ## 6                  463.46000                           69.2807
    ##   x.blu76_male_pituitary_m.hatch x.blu94_female_gonad_m.inc.d17
    ## 1                        205.000                      113.00000
    ## 2                         34.000                      207.00000
    ## 3                        325.000                      829.00000
    ## 4                          0.000                        7.00000
    ## 5                          0.000                       11.88084
    ## 6                        139.204                       70.86650
    ##   x.blu94_female_hypothalamus_m.inc.d17 x.blu94_female_pituitary_m.inc.d17
    ## 1                                    83                            283.000
    ## 2                                     0                            124.000
    ## 3                                   184                            633.000
    ## 4                                     0                              0.000
    ## 5                                     2                              2.000
    ## 6                                    13                            139.258
    ##   x.blu96_female_gonad_m.n2 x.blu96_female_hypothalamus_m.n2
    ## 1                        89                          248.000
    ## 2                       421                            7.000
    ## 3                       820                          524.000
    ## 4                         8                            0.000
    ## 5                         0                            0.000
    ## 6                        85                           41.083
    ##   x.blu96_female_pituitary_m.n2 x.g.ATLAS_male_gonad_control
    ## 1                            76                            0
    ## 2                            10                            0
    ## 3                           203                          212
    ## 4                             0                            0
    ## 5                             1                            0
    ## 6                            49                           19
    ##   x.g.ATLAS_male_hypothalamus_control x.g.ATLAS_male_pituitary_control
    ## 1                            93.00000                         159.0000
    ## 2                             1.00000                          12.0000
    ## 3                           149.00000                         297.0000
    ## 4                             1.00000                           0.0000
    ## 5                             0.00000                           0.0000
    ## 6                             8.05182                          21.7122
    ##   x.g.g.ATLAS_female_gonad_control x.g.g.ATLAS_male_gonad_control
    ## 1                               80                         14.000
    ## 2                               49                          3.000
    ## 3                              542                       1126.000
    ## 4                               14                          0.000
    ## 5                                2                          2.000
    ## 6                               19                        285.787
    ##   x.g.g.ATLAS_male_hypothalamus_control x.g.g.ATLAS_male_pituitary_control
    ## 1                                    16                                 13
    ## 2                                     2                                 58
    ## 3                                   199                                240
    ## 4                                     0                                  0
    ## 5                                     0                                  0
    ## 6                                     3                                 21
    ##   x.g.g.g.ATLAS_male_gonad_control x.g.g.g.ATLAS_male_pituitary_control
    ## 1                            9.000                                  126
    ## 2                            2.000                                   52
    ## 3                         1142.260                                  214
    ## 4                            0.000                                    0
    ## 5                           10.000                                    0
    ## 6                          222.832                                   14
    ##   x.g13.w109_male_gonad_inc.d9 x.g13.w109_male_hypothalamus_inc.d9
    ## 1                            4                            243.0000
    ## 2                            0                              0.0000
    ## 3                          573                            502.0000
    ## 4                            1                              0.0000
    ## 5                            1                              1.0000
    ## 6                          149                             45.5407
    ##   x.g13.w109_male_pituitary_inc.d9 x.g14.w199_male_gonad_inc.d17
    ## 1                         133.0000                             5
    ## 2                           9.0000                             1
    ## 3                         297.0000                           344
    ## 4                           0.0000                             2
    ## 5                           1.0000                             5
    ## 6                          64.2791                           121
    ##   x.g14.w199_male_hypothalamus_inc.d17 x.g14.w199_male_pituitary_inc.d17
    ## 1                             361.0000                          329.3536
    ## 2                               2.0000                           54.0000
    ## 3                             648.0000                          738.0000
    ## 4                               0.0000                            0.0000
    ## 5                               0.0000                            5.0000
    ## 6                              35.8286                          162.6550
    ##   x.g147.blu28_male_gonad_inc.d3 x.g147.blu28_male_hypothalamus_inc.d3
    ## 1                          9.000                              261.0000
    ## 2                          3.000                                5.0000
    ## 3                       1562.000                              985.0000
    ## 4                        124.000                                0.0000
    ## 5                          3.000                                2.0000
    ## 6                        553.985                               43.1875
    ##   x.g147.blu28_male_pituitary_inc.d3 x.g26_female_gonad_m.inc.d8
    ## 1                            106.000                      130.00
    ## 2                             13.000                       84.00
    ## 3                            236.000                      949.00
    ## 4                              0.000                      231.00
    ## 5                              1.000                        3.00
    ## 6                             54.323                      194.43
    ##   x.g26_female_hypothalamus_m.inc.d8 x.g26_female_pituitary_m.inc.d8
    ## 1                            328.000                         296.000
    ## 2                              8.000                          81.000
    ## 3                            592.000                         567.000
    ## 4                              0.000                           0.000
    ## 5                              1.000                           5.000
    ## 6                             46.645                         146.854
    ##   x.g37_female_gonad_n5 x.g37_female_hypothalamus_n5
    ## 1                    42                          261
    ## 2                    52                            2
    ## 3                   205                          648
    ## 4                     2                            0
    ## 5                     1                            3
    ## 6                    41                           43
    ##   x.g37_female_pituitary_n5 x.g4.w50_female_gonad_n9
    ## 1                   80.0000                       55
    ## 2                   37.0000                       38
    ## 3                  225.0000                      261
    ## 4                    0.0000                        4
    ## 5                    0.0000                        2
    ## 6                   61.1526                       61
    ##   x.g4.w50_female_hypothalamus_n9 x.g4.w50_female_pituitary_n9
    ## 1                             155                           85
    ## 2                               1                            9
    ## 3                             482                          195
    ## 4                               0                            0
    ## 5                               0                            2
    ## 6                              36                           43
    ##   x.g43_female_gonad_n5 x.g43_female_hypothalamus_n5
    ## 1                    85                          148
    ## 2                    60                            2
    ## 3                   202                          363
    ## 4                   105                            1
    ## 5                     0                            0
    ## 6                    18                           21
    ##   x.g43_female_pituitary_n5 x.g49_female_gonad_n5
    ## 1                   74.0000                   141
    ## 2                    7.0000                   132
    ## 3                  184.0000                   729
    ## 4                    0.0000                     0
    ## 5                    4.0000                     5
    ## 6                   64.7108                    93
    ##   x.g49_female_hypothalamus_n5 x.g49_female_pituitary_n5.NYNO
    ## 1                     172.0000                        100.000
    ## 2                       2.0000                         11.000
    ## 3                     475.0000                        391.000
    ## 4                       4.0000                          0.000
    ## 5                       4.0000                          2.000
    ## 6                      36.0876                        130.174
    ##   x.g5.w79_male_gonad_m.inc.d3 x.g5.w79_male_hypothalamus_m.inc.d3
    ## 1                        9.000                             98.0000
    ## 2                        0.000                              3.0000
    ## 3                     2103.000                            329.0000
    ## 4                        0.000                              0.0000
    ## 5                        4.000                              0.0000
    ## 6                      479.937                             16.0705
    ##   x.g5.w79_male_pituitary_m.inc.d3 x.g50_female_gonad_inc.prolong
    ## 1                          225.000                            125
    ## 2                           32.000                            138
    ## 3                          519.000                            420
    ## 4                            0.000                             18
    ## 5                            2.000                              1
    ## 6                          138.927                             39
    ##   x.g50_female_hypothalamus_inc.prolong x.g50_female_pituitary_inc.prolong
    ## 1                             234.00000                           101.0000
    ## 2                               1.00000                            20.0000
    ## 3                             606.00000                           313.0000
    ## 4                               0.00000                             0.0000
    ## 5                               0.64101                             0.0000
    ## 6                              63.33810                            62.5335
    ##   x.g70_male_gonad_hatch x.g70_male_hypothalamus_hatch
    ## 1                     10                      252.0000
    ## 2                      0                        4.0000
    ## 3                    440                      623.0000
    ## 4                      0                        0.0000
    ## 5                      4                        0.0000
    ## 6                    151                       31.3563
    ##   x.g70_male_pituitary_hatch x.g9.o166_female_gonad_inc.d9.NYNO
    ## 1                    95.0000                                 79
    ## 2                     9.0000                                 58
    ## 3                   209.0000                                609
    ## 4                     0.0000                                 31
    ## 5                     0.0000                                  1
    ## 6                    74.9382                                 36
    ##   x.g9.o166_female_hypothalamus_inc.d9
    ## 1                                  360
    ## 2                                    4
    ## 3                                  418
    ## 4                                    0
    ## 5                                    0
    ## 6                                   38
    ##   x.g9.o166_female_pituitary_inc.d9.NYNO x.o117_male_gonad_m.inc.d9
    ## 1                               190.0000                    6.00000
    ## 2                                 8.0000                    1.00000
    ## 3                               440.0000                  564.00000
    ## 4                                 0.0000                    0.00000
    ## 5                                 1.0000                    3.66482
    ## 6                                95.7844                  185.00000
    ##   x.o117_male_hypothalamus_m.inc.d9 x.o117_male_pituitary_m.inc.d9
    ## 1                                96                             78
    ## 2                                 1                             74
    ## 3                               243                            226
    ## 4                                 0                              0
    ## 5                                 0                              1
    ## 6                                16                             57
    ##   x.o159.w90_female_gonad_inc.d17 x.o159.w90_female_hypothalamus_inc.d17
    ## 1                              45                                432.000
    ## 2                               9                                  2.000
    ## 3                             220                                549.000
    ## 4                              22                                  0.000
    ## 5                               0                                  1.000
    ## 6                              26                                 74.282
    ##   x.o159.w90_female_pituitary_inc.d17 x.o160.w102_male_gonad_hatch
    ## 1                              62.000                           12
    ## 2                               9.000                            3
    ## 3                             185.000                          363
    ## 4                               0.000                            0
    ## 5                               0.000                            1
    ## 6                              62.627                          151
    ##   x.o160.w102_male_hypothalamus_hatch x.o160.w102_male_pituitary_hatch
    ## 1                                 329                              117
    ## 2                                   3                               20
    ## 3                                 610                              188
    ## 4                                   0                                0
    ## 5                                   0                                0
    ## 6                                  36                               63
    ##   x.o163.w101_male_gonad_inc.d3 x.o163.w101_male_hypothalamus_inc.d3
    ## 1                             1                             360.0000
    ## 2                             0                               8.0000
    ## 3                           435                             664.0000
    ## 4                             0                               0.0000
    ## 5                             1                               2.0000
    ## 6                            96                              47.3857
    ##   x.o163.w101_male_pituitary_inc.d3.NYNO x.o164.w123_male_gonad_n5
    ## 1                                244.000                  28.00000
    ## 2                                 12.000                   6.00000
    ## 3                                540.000                1226.00000
    ## 4                                  1.000                   6.00000
    ## 5                                  1.000                   9.66689
    ## 6                                136.983                 508.00000
    ##   x.o164.w123_male_hypothalamus_n5 x.o164.w123_male_pituitary_n5.NYNO
    ## 1                         282.0000                           162.0000
    ## 2                           5.0000                            24.0000
    ## 3                         521.0000                           398.0000
    ## 4                           0.0000                             0.0000
    ## 5                           1.0000                             2.0000
    ## 6                          37.5231                            96.5617
    ##   x.o171.w45_female_gonad_m.hatch x.o171.w45_female_hypothalamus_m.hatch
    ## 1                        80.00000                               295.0000
    ## 2                       152.00000                                 3.0000
    ## 3                       493.00000                               564.0000
    ## 4                         8.00000                                 0.0000
    ## 5                        12.85527                                 0.0000
    ## 6                        70.00000                                36.1629
    ##   x.o171.w45_female_pituitary_m.hatch x.o175.g21_female_gonad_n5
    ## 1                             271.000                        137
    ## 2                              69.000                        119
    ## 3                             600.000                        648
    ## 4                               0.000                          5
    ## 5                               1.000                          1
    ## 6                             180.606                         99
    ##   x.o175.g21_female_hypothalamus_n5 x.o175.g21_female_pituitary_n5
    ## 1                          170.0000                        316.000
    ## 2                            2.0000                         35.000
    ## 3                          470.0000                        989.000
    ## 4                            0.0000                          0.000
    ## 5                            0.0000                          3.000
    ## 6                           24.0849                        321.965
    ##   x.o2_male_gonad_n9 x.o2_male_hypothalamus_n9 x.o2_male_pituitary_n9
    ## 1                 17                       224                187.000
    ## 2                  0                         5                 31.000
    ## 3               1960                       397                456.000
    ## 4                  3                         0                  0.000
    ## 5                  5                         1                  0.000
    ## 6                352                        30                133.531
    ##   x.o30.g134_male_gonad_bldg x.o30.g134_male_hypothalamus_bldg
    ## 1                   16.00000                               177
    ## 2                    6.00000                                 1
    ## 3                 2493.00000                               365
    ## 4                    4.00000                                 0
    ## 5                   13.02614                                 0
    ## 6                  647.75300                                30
    ##   x.o30.g134_male_pituitary_bldg x.o37.blu50_female_gonad_hatch
    ## 1                         90.000                             40
    ## 2                         36.000                             38
    ## 3                        190.000                            210
    ## 4                          0.000                             10
    ## 5                          0.000                              1
    ## 6                         47.377                             28
    ##   x.o37.blu50_female_hypothalamus_hatch.NYNO
    ## 1                                   445.0000
    ## 2                                     6.0000
    ## 3                                   532.0000
    ## 4                                     0.0000
    ## 5                                     2.0000
    ## 6                                    45.1149
    ##   x.o37.blu50_female_pituitary_hatch x.o40.r70_male_gonad_m.inc.d17
    ## 1                                 94                         38.000
    ## 2                                  8                          7.000
    ## 3                                318                       1515.000
    ## 4                                  0                          8.000
    ## 5                                  0                          9.000
    ## 6                                 85                        608.791
    ##   x.o40.r70_male_hypothalamus_m.inc.d17 x.o40.r70_male_pituitary_m.inc.d17
    ## 1                               545.000                            263.000
    ## 2                                 2.000                            200.000
    ## 3                              1544.070                            559.000
    ## 4                                 0.000                              0.000
    ## 5                                 2.000                              2.000
    ## 6                                76.367                            152.371
    ##   x.o47.y82_male_gonad_inc.d9 x.o47.y82_male_hypothalamus_inc.d9
    ## 1                           3                                187
    ## 2                           1                                  1
    ## 3                         445                                261
    ## 4                           2                                  0
    ## 5                           4                                  0
    ## 6                         120                                 17
    ##   x.o47.y82_male_pituitary_inc.d9 x.o4_male_gonad_m.hatch
    ## 1                         232.000                  18.000
    ## 2                          39.000                   4.000
    ## 3                         620.000                1559.000
    ## 4                           0.000                   1.000
    ## 5                           0.000                   7.000
    ## 6                         161.228                 574.817
    ##   x.o4_male_hypothalamus_m.hatch x.o4_male_pituitary_m.hatch
    ## 1                       274.0000                     426.000
    ## 2                         6.0000                     273.000
    ## 3                       769.0000                     914.000
    ## 4                         0.0000                       0.000
    ## 5                         3.0000                       4.000
    ## 6                        49.1308                     259.865
    ##   x.o61_female_gonad_extend.hatch x.o61_female_hypothalamus_extend.hatch
    ## 1                          92.000                               229.0000
    ## 2                         141.000                                 3.0000
    ## 3                        1154.000                               482.0000
    ## 4                          57.000                                 0.0000
    ## 5                           0.000                                 1.0000
    ## 6                         308.829                                34.8294
    ##   x.o61_female_pituitary_extend.hatch x.o68_male_gonad_n5
    ## 1                             210.000              66.000
    ## 2                              33.000               1.000
    ## 3                             474.000            1465.000
    ## 4                               0.000               1.000
    ## 5                               0.000               3.000
    ## 6                             126.266             454.917
    ##   x.o68_male_hypothalamus_n5 x.o68_male_pituitary_n5 x.o70_female_gonad_n5
    ## 1                   206.0000                 387.000                    46
    ## 2                     4.0000                  46.000                    78
    ## 3                   421.0000                 837.000                   279
    ## 4                     0.0000                   0.000                     2
    ## 5                     1.0000                   3.000                     1
    ## 6                    21.0568                 320.081                    44
    ##   x.o70_female_hypothalamus_n5.NYNO x.o70_female_pituitary_n5
    ## 1                          160.0000                   162.000
    ## 2                            5.0000                    12.000
    ## 3                          453.0000                   714.000
    ## 4                            0.0000                     0.000
    ## 5                            2.0000                     4.000
    ## 6                           41.1347                   205.473
    ##   x.r10.w18_male_gonad_m.inc.d17 x.r10.w18_male_hypothalamus_m.inc.d17
    ## 1                       82.00000                               70.0000
    ## 2                        9.00000                                1.0000
    ## 3                     1939.00000                              288.0000
    ## 4                        0.00000                                0.0000
    ## 5                       10.43282                                0.0000
    ## 6                      752.89100                               21.3397
    ##   x.r10.w18_male_pituitary_m.inc.d17 x.r178_male_gonad_hatch
    ## 1                            205.000                30.00000
    ## 2                            308.000                 3.00000
    ## 3                            406.000              1919.00000
    ## 4                              0.000                 0.00000
    ## 5                              0.000                14.59269
    ## 6                            141.758               679.00000
    ##   x.r178_male_hypothalamus_hatch x.r178_male_pituitary_hatch
    ## 1                       425.0000                     86.0000
    ## 2                         2.0000                     12.0000
    ## 3                       581.0000                    199.0000
    ## 4                         0.0000                      1.0000
    ## 5                         0.0000                      0.0000
    ## 6                        52.3416                     84.8959
    ##   x.r180_female_gonad_m.inc.d3 x.r180_female_hypothalamus_m.inc.d3
    ## 1                      255.000                                  97
    ## 2                      359.000                                   0
    ## 3                     1033.000                                 312
    ## 4                       53.000                                   0
    ## 5                        9.000                                   0
    ## 6                      236.976                                  14
    ##   x.r180_female_pituitary_m.inc.d3 x.r181_male_gonad_n5
    ## 1                          180.000                   11
    ## 2                          155.000                    1
    ## 3                          616.000                  417
    ## 4                            0.000                    0
    ## 5                            4.000                    1
    ## 6                          163.813                  121
    ##   x.r181_male_hypothalamus_n5 x.r181_male_pituitary_n5
    ## 1                    409.2810                       60
    ## 2                      7.0000                        3
    ## 3                   1063.0000                      161
    ## 4                      0.0000                        0
    ## 5                      1.0000                        0
    ## 6                     48.3819                       59
    ##   x.r185_male_gonad_m.inc.d3 x.r185_male_hypothalamus_m.inc.d3
    ## 1                   11.00000                          407.0000
    ## 2                    2.00000                            5.0000
    ## 3                 1062.00000                          637.0000
    ## 4                    0.00000                            0.0000
    ## 5                    6.17202                            3.0000
    ## 6                  459.00000                           39.3343
    ##   x.r185_male_pituitary_m.inc.d3 x.r29.w96_male_gonad_inc.d17
    ## 1                        535.000                            8
    ## 2                        366.000                            0
    ## 3                       1494.000                          371
    ## 4                          0.000                            0
    ## 5                          0.000                            0
    ## 6                        357.461                          210
    ##   x.r29.w96_male_hypothalamus_inc.d17 x.r29.w96_male_pituitary_inc.d17
    ## 1                            450.0000                               81
    ## 2                              8.0000                               14
    ## 3                           1060.0000                              173
    ## 4                              0.0000                                0
    ## 5                              0.0000                                0
    ## 6                             49.2768                               37
    ##   x.r33.w183_female_gonad_inc.d3 x.r33.w183_female_hypothalamus_inc.d3
    ## 1                             49                              198.0000
    ## 2                            126                                1.0000
    ## 3                            347                              537.0000
    ## 4                            319                                0.0000
    ## 5                              1                                0.0000
    ## 6                             42                               35.7369
    ##   x.r33.w183_female_pituitary_inc.d3 x.r39.g10_female_gonad_bldg
    ## 1                            97.0000                          28
    ## 2                            15.0000                           8
    ## 3                           281.0000                         324
    ## 4                             0.0000                          61
    ## 5                             2.0000                           2
    ## 6                            67.5426                          46
    ##   x.r39.g10_female_hypothalamus_bldg x.r39.g10_female_pituitary_bldg
    ## 1                           111.0000                         196.000
    ## 2                             1.0000                          71.000
    ## 3                           428.0000                         606.000
    ## 4                             0.0000                          11.000
    ## 5                             0.0000                           2.000
    ## 6                            22.3529                         195.592
    ##   x.r44.w95_female_gonad_hatch x.r44.w95_female_hypothalamus_hatch
    ## 1                      74.1788                            246.0000
    ## 2                      92.0000                              2.0000
    ## 3                     338.0000                            931.0000
    ## 4                      12.0000                              0.0000
    ## 5                       6.0000                              5.0000
    ## 6                      34.0000                             58.0675
    ##   x.r44.w95_female_pituitary_hatch x.r48.y139_female_gonad_inc.d17
    ## 1                               35                              23
    ## 2                                4                              68
    ## 3                              285                             166
    ## 4                                1                              36
    ## 5                                0                               2
    ## 6                               67                              13
    ##   x.r48.y139_female_hypothalamus_inc.d17
    ## 1                               271.0000
    ## 2                                 1.0000
    ## 3                               633.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                49.1383
    ##   x.r48.y139_female_pituitary_inc.d17.NYNO x.r50.w97_female_gonad_n5
    ## 1                                  128.000                       101
    ## 2                                   29.000                       189
    ## 3                                  822.000                       969
    ## 4                                    0.000                         5
    ## 5                                    0.000                        10
    ## 6                                  166.599                       362
    ##   x.r50.w97_female_hypothalamus_n5 x.r50.w97_female_pituitary_n5
    ## 1                         464.0000                           358
    ## 2                           7.0000                            60
    ## 3                         972.0000                           793
    ## 4                           0.0000                             0
    ## 5                           1.0000                             0
    ## 6                          49.4283                           183
    ##   x.r64.g140_male_gonad_inc.d3 x.r64.g140_male_hypothalamus_inc.d3
    ## 1                       21.000                           334.00000
    ## 2                        6.000                             8.00000
    ## 3                     1120.000                           676.00000
    ## 4                        4.000                             0.00000
    ## 5                        1.000                             8.86496
    ## 6                      336.924                            55.36700
    ##   x.r64.g140_male_pituitary_inc.d3 x.r67.blu35_male_gonad_bldg
    ## 1                          101.000                           1
    ## 2                           59.000                           3
    ## 3                          362.000                         460
    ## 4                            0.000                           2
    ## 5                            3.000                           4
    ## 6                          116.705                         183
    ##   x.r67.blu35_male_hypothalamus_bldg x.r67.blu35_male_pituitary_bldg.NYNO
    ## 1                             96.000                            195.00000
    ## 2                              3.000                             45.00000
    ## 3                            616.000                            465.00000
    ## 4                              0.000                              0.00000
    ## 5                              0.000                              4.10243
    ## 6                             53.623                            147.86600
    ##   x.r74.o29_female_gonad_m.hatch x.r74.o29_female_hypothalamus_m.hatch
    ## 1                             75                                332.00
    ## 2                             82                                  7.00
    ## 3                            607                                566.00
    ## 4                             14                                  0.00
    ## 5                              4                                  1.00
    ## 6                            117                                 48.23
    ##   x.r74.o29_female_pituitary_m.hatch x.r76_female_gonad_m.inc.d9
    ## 1                            444.553                          63
    ## 2                             68.000                         142
    ## 3                           1003.000                         344
    ## 4                              0.000                          10
    ## 5                              0.000                           4
    ## 6                            331.455                          32
    ##   x.r76_female_hypothalamus_m.inc.d9 x.r76_female_pituitary_m.inc.d9
    ## 1                           239.0000                         53.0000
    ## 2                             0.0000                         13.0000
    ## 3                           610.0000                        178.0000
    ## 4                             0.0000                          0.0000
    ## 5                             0.0000                          2.0000
    ## 6                            26.1712                         18.2918
    ##   x.s188_female_gonad_m.inc.d8 x.s188_female_hypothalamus_m.inc.d8
    ## 1                     241.0000                                 244
    ## 2                     324.0000                                   5
    ## 3                     705.0000                                 557
    ## 4                      35.0000                                   0
    ## 5                       3.0000                                   1
    ## 6                      87.9361                                  38
    ##   x.s188_female_pituitary_m.inc.d8 x.w178_female_gonad_n9
    ## 1                        374.00000                60.0000
    ## 2                        136.00000                64.0000
    ## 3                        966.00000               344.0000
    ## 4                          0.00000                 0.0000
    ## 5                          4.03775                 2.0000
    ## 6                        352.99000                64.9316
    ##   x.w178_female_hypothalamus_n9 x.w178_female_pituitary_n9
    ## 1                      336.0000                        100
    ## 2                        4.0000                         30
    ## 3                      915.0000                        275
    ## 4                        0.0000                          0
    ## 5                        0.0000                          0
    ## 6                       54.5463                         82
    ##   x.w192.o157_male_gonad_inc.d9 x.w192.o157_male_hypothalamus_inc.d9
    ## 1                       5.00000                             388.0000
    ## 2                       1.00000                               7.0000
    ## 3                     498.00000                             658.0000
    ## 4                       4.00000                              37.0000
    ## 5                       6.88128                               1.0000
    ## 6                     172.00000                              65.3199
    ##   x.w192.o157_male_pituitary_inc.d9 x.w51_female_gonad_lay
    ## 1                           169.000                    175
    ## 2                           141.000                    177
    ## 3                           554.000                    730
    ## 4                             0.000                     47
    ## 5                             1.000                      5
    ## 6                           226.246                     81
    ##   x.w51_female_hypothalamus_lay x.w51_female_pituitary_lay
    ## 1                      433.3630                    104.000
    ## 2                       18.0000                     16.000
    ## 3                      541.0000                    530.000
    ## 4                        0.0000                      0.000
    ## 5                        0.0000                      1.000
    ## 6                       44.2758                    228.802
    ##   x.w6_female_gonad_n9 x.w6_female_hypothalamus_n9
    ## 1                  165                    153.0000
    ## 2                   35                      0.0000
    ## 3                  798                    449.0000
    ## 4                   21                      0.0000
    ## 5                    8                      1.0000
    ## 6                  188                     20.0494
    ##   x.w6_female_pituitary_n9 x.y.s.ATLAS_male_gonad_control
    ## 1                  160.000                         13.000
    ## 2                   36.000                          0.000
    ## 3                  476.000                        772.000
    ## 4                    0.000                         15.000
    ## 5                    2.000                          3.000
    ## 6                  159.063                        135.693
    ##   x.y.s.ATLAS_male_pituitary_control x.y109_female_gonad_inc.d9
    ## 1                                 27                         27
    ## 2                                 55                         64
    ## 3                                139                        310
    ## 4                                  0                         14
    ## 5                                  0                          2
    ## 6                                 10                         28
    ##   x.y109_female_hypothalamus_inc.d9 x.y109_female_pituitary_inc.d9
    ## 1                          115.0000                             89
    ## 2                            0.0000                             13
    ## 3                          385.0000                            216
    ## 4                            0.0000                              0
    ## 5                            0.0000                              0
    ## 6                           21.1942                             35
    ##   x.y119.w11_female_gonad_m.inc.d17
    ## 1                                55
    ## 2                               108
    ## 3                               349
    ## 4                                 5
    ## 5                                 2
    ## 6                                75
    ##   x.y119.w11_female_hypothalamus_m.inc.d17
    ## 1                                 106.0000
    ## 2                                   2.0000
    ## 3                                 255.0000
    ## 4                                   0.0000
    ## 5                                   0.0000
    ## 6                                  14.0527
    ##   x.y119.w11_female_pituitary_m.inc.d17 x.y132.w76_male_gonad_inc.d17
    ## 1                               223.000                            12
    ## 2                               118.000                             2
    ## 3                               481.000                           661
    ## 4                                 4.000                             0
    ## 5                                 5.000                             2
    ## 6                               137.661                           234
    ##   x.y132.w76_male_hypothalamus_inc.d17 x.y132.w76_male_pituitary_inc.d17
    ## 1                             213.0000                                83
    ## 2                               3.0000                                53
    ## 3                             470.0000                               325
    ## 4                               0.0000                                 1
    ## 5                               0.0000                                 0
    ## 6                              50.1584                                86
    ##   x.y138.w176_female_gonad_n9 x.y138.w176_female_hypothalamus_n9
    ## 1                          12                            261.000
    ## 2                          50                              6.000
    ## 3                         426                            573.000
    ## 4                           0                              0.000
    ## 5                           0                              1.000
    ## 6                          90                             48.274
    ##   x.y138.w176_female_pituitary_n9 x.y141.w116_male_gonad_inc.d9
    ## 1                       92.000000                        20.000
    ## 2                       38.000000                         0.000
    ## 3                      437.000000                      1705.000
    ## 4                        0.000000                         7.000
    ## 5                        2.976245                         3.000
    ## 6                       69.697900                       434.976
    ##   x.y141.w116_male_hypothalamus_inc.d9 x.y141.w116_male_pituitary_inc.d9
    ## 1                             166.0000                           358.000
    ## 2                               3.0000                            71.000
    ## 3                             372.0000                           707.000
    ## 4                               0.0000                             0.000
    ## 5                               0.0000                             1.000
    ## 6                              25.3426                           166.149
    ##   x.y145.r55_female_gonad_inc.prolong
    ## 1                                  50
    ## 2                                  62
    ## 3                                 281
    ## 4                                   3
    ## 5                                   1
    ## 6                                  41
    ##   x.y145.r55_female_hypothalamus_inc.prolong
    ## 1                                   300.0000
    ## 2                                     0.0000
    ## 3                                   739.0000
    ## 4                                     0.0000
    ## 5                                     1.0000
    ## 6                                    65.8116
    ##   x.y145.r55_female_pituitary_inc.prolong x.y51_male_gonad_m.inc.d8
    ## 1                                      75                         8
    ## 2                                      20                         3
    ## 3                                     242                       607
    ## 4                                       0                         0
    ## 5                                       1                         2
    ## 6                                      51                       435
    ##   x.y51_male_hypothalamus_m.inc.d8 x.y51_male_pituitary_m.inc.d8
    ## 1                         257.0000                      422.0000
    ## 2                           2.0000                       85.0000
    ## 3                         830.0000                      464.0000
    ## 4                           0.0000                        0.0000
    ## 5                           2.0000                        2.8266
    ## 6                          88.3668                      262.1060
    ##   x.y5_female_gonad_m.inc.d8 x.y5_female_hypothalamus_m.inc.d8
    ## 1                        101                          249.0000
    ## 2                        313                           14.0000
    ## 3                       1243                          479.0000
    ## 4                         12                            0.0000
    ## 5                          2                            0.0000
    ## 6                         67                           23.2768
    ##   x.y5_female_pituitary_m.inc.d8 x.y64_female_gonad_extend
    ## 1                       309.3845                       138
    ## 2                        29.0000                       255
    ## 3                       630.0000                       642
    ## 4                         0.0000                        24
    ## 5                         5.0000                         0
    ## 6                       135.6580                        93
    ##   x.y64_female_hypothalamus_extend x.y64_female_pituitary_extend
    ## 1                         336.0000                         249.0
    ## 2                           3.0000                          13.0
    ## 3                         557.0000                         572.0
    ## 4                           0.0000                           0.0
    ## 5                           1.0000                           2.0
    ## 6                          48.2873                         206.6
    ##   x.y79_male_gonad_prolong x.y79_male_hypothalamus_prolong
    ## 1                   32.000                        385.0000
    ## 2                    1.000                          3.0000
    ## 3                  378.000                        649.0000
    ## 4                    0.000                          0.0000
    ## 5                    5.000                          0.0000
    ## 6                  377.955                         64.2577
    ##   x.y79_male_pituitary_prolong x.y90_female_gonad_hatch
    ## 1                     261.3429                       36
    ## 2                     122.0000                      106
    ## 3                     579.0000                      338
    ## 4                       1.0000                        7
    ## 5                       2.0000                        0
    ## 6                     246.1290                       20
    ##   x.y90_female_hypothalamus_hatch x.y90_female_pituitary_hatch
    ## 1                              87                      144.000
    ## 2                               0                       12.000
    ## 3                             197                      678.000
    ## 4                               0                        0.000
    ## 5                               0                        1.000
    ## 6                               5                      244.345
    ##   x.y93.g126_female_gonad_inc.d9 x.y93.g126_female_hypothalamus_inc.d9
    ## 1                        33.0000                                   122
    ## 2                       117.0000                                     0
    ## 3                       292.0000                                   407
    ## 4                        18.0000                                     0
    ## 5                         0.0000                                     1
    ## 6                        39.8289                                    37
    ##   x.y93.g126_female_pituitary_inc.d9 x.y9_female_gonad_n9
    ## 1                            260.000                   85
    ## 2                             57.000                   55
    ## 3                            627.000                  602
    ## 4                              0.000                   14
    ## 5                              2.000                    0
    ## 6                            220.045                  332
    ##   x.y9_female_hypothalamus_n9 x.y9_female_pituitary_n9
    ## 1                    179.0000                 299.4387
    ## 2                      0.0000                  26.0000
    ## 3                    386.0000                 414.0000
    ## 4                      0.0000                   0.0000
    ## 5                      3.0000                   4.0000
    ## 6                     27.4368                 166.2900
    ##   y.s156.o.r_female_gonad_lay y.s156.o.r_female_hypothalamus_lay
    ## 1                   586.00000                           848.6660
    ## 2                    71.00000                             5.0000
    ## 3                   427.00000                          1314.0000
    ## 4                    64.00000                             0.0000
    ## 5                    13.03804                             3.0000
    ## 6                    50.00000                            66.7481
    ##   y.s156.o.r_female_pituitary_lay y.w.blk_male_gonad_m.inc.d9
    ## 1                         419.583                     5.03668
    ## 2                          30.000                     3.00000
    ## 3                         772.000                   610.00000
    ## 4                           0.000                     0.00000
    ## 5                           1.000                     5.00000
    ## 6                         201.084                   175.00000
    ##   y.w.blk_male_hypothalamus_m.inc.d9 y.w.blk_male_pituitary_m.inc.d9
    ## 1                           293.0000                              76
    ## 2                             2.0000                              13
    ## 3                           638.0000                             210
    ## 4                             0.0000                               0
    ## 5                             1.0000                               0
    ## 6                            71.7467                              50
    ##   y12.x_female_gonad_m.inc.d9 y12.x_female_hypothalamus_m.inc.d9
    ## 1                          33                           201.0000
    ## 2                          43                             1.0000
    ## 3                         225                           289.0000
    ## 4                          21                             0.0000
    ## 5                           2                             0.0000
    ## 6                          25                            23.2119
    ##   y12.x_female_pituitary_m.inc.d9 y126.w92.x_female_gonad_inc.d17
    ## 1                         70.0000                              30
    ## 2                          5.0000                              80
    ## 3                        281.0000                             289
    ## 4                          0.0000                              23
    ## 5                          1.0000                               1
    ## 6                         49.0215                              37
    ##   y126.w92.x_female_hypothalamus_inc.d17
    ## 1                               445.0000
    ## 2                                 7.0000
    ## 3                               456.0000
    ## 4                                 0.0000
    ## 5                                 3.0000
    ## 6                                46.4953
    ##   y126.w92.x_female_pituitary_inc.d17 y128.g23.x_female_gonad_inc.d9
    ## 1                             72.0000                        41.0000
    ## 2                              3.0000                         5.0000
    ## 3                            173.0000                       303.0000
    ## 4                              0.0000                         7.0000
    ## 5                              0.0000                         1.0000
    ## 6                             66.9048                        79.9014
    ##   y128.g23.x_female_hypothalamus_m.inc.d9
    ## 1                                283.0000
    ## 2                                  6.0000
    ## 3                                639.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 40.4052
    ##   y128.g23.x_female_pituitary_inc.d9 y129.x_male_gonad_n9
    ## 1                            70.2786                7.000
    ## 2                            22.0000                1.000
    ## 3                           236.0000              686.000
    ## 4                             0.0000                2.000
    ## 5                             2.0000                2.000
    ## 6                           150.5000              246.735
    ##   y129.x_male_hypothalamus_n9 y129.x_male_pituitary_n9
    ## 1                     145.000                  77.0000
    ## 2                       0.000                  59.0000
    ## 3                     519.000                 143.0000
    ## 4                       0.000                   0.0000
    ## 5                       1.000                   0.0000
    ## 6                      33.115                  50.6307
    ##   y13.x_female_gonad_inc.d3 y13.x_female_hypothalamus_inc.d3
    ## 1                        63                         183.0000
    ## 2                        23                           0.0000
    ## 3                       280                         430.0000
    ## 4                         2                           0.0000
    ## 5                         0                           0.0000
    ## 6                        19                          38.7971
    ##   y13.x_female_pituitary_inc.d3 y130.o170.x_female_gonad_inc.d17
    ## 1                            86                               44
    ## 2                            27                              119
    ## 3                           254                              171
    ## 4                             0                               10
    ## 5                             0                                0
    ## 6                            42                               33
    ##   y130.o170.x_female_hypothalamus_inc.d17
    ## 1                                287.0000
    ## 2                                  0.0000
    ## 3                                577.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 65.3311
    ##   y130.o170.x_female_pituitary_inc.d17 y131.w185.x_male_gonad_n9
    ## 1                                   91                         3
    ## 2                                   15                         1
    ## 3                                  206                       504
    ## 4                                    0                         1
    ## 5                                    1                         4
    ## 6                                   73                       130
    ##   y131.w185.x_male_hypothalamus_n9 y131.w185.x_male_pituitary_n9
    ## 1                         294.0000                       204.000
    ## 2                           4.0000                       183.000
    ## 3                         704.0000                       662.000
    ## 4                           0.0000                         0.000
    ## 5                           1.0000                         0.000
    ## 6                          39.1707                       152.338
    ##   y133.w77.r58_male_gonad_inc.d17 y133.w77.r58_male_hypothalamus_inc.d17
    ## 1                              12                               296.0000
    ## 2                               0                                 8.0000
    ## 3                             514                               582.0000
    ## 4                               0                                 0.0000
    ## 5                               0                                 2.7485
    ## 6                             137                                50.7092
    ##   y133.w77.r58_male_pituitary_inc.d17 y135.blu107.x_female_gonad_inc.d17
    ## 1                           102.00000                                 35
    ## 2                            21.00000                                107
    ## 3                           206.00000                                291
    ## 4                             0.00000                                 23
    ## 5                             1.08287                                  1
    ## 6                            75.00000                                 32
    ##   y135.blu107.x_female_hypothalamus_inc.d17
    ## 1                                  388.0000
    ## 2                                    5.0000
    ## 3                                  646.0000
    ## 4                                    0.0000
    ## 5                                    2.0000
    ## 6                                   55.3267
    ##   y135.blu107.x_female_pituitary_inc.d17.NYNO y136.x_female_gonad_inc.d17
    ## 1                                      209.00                          33
    ## 2                                        9.00                          65
    ## 3                                      444.00                         192
    ## 4                                        0.00                          57
    ## 5                                        0.00                           1
    ## 6                                      163.22                          20
    ##   y136.x_female_hypothalamus_inc.d17 y136.x_female_pituitary_inc.d17
    ## 1                           253.0000                        106.0000
    ## 2                             1.0000                          7.0000
    ## 3                           595.0000                        271.0000
    ## 4                             0.0000                          0.0000
    ## 5                             0.0000                          0.0000
    ## 6                            38.0986                         44.2683
    ##   y140.w119.x_female_gonad_inc.d9 y140.w119.x_female_hypothalamus_inc.d9
    ## 1                              66                                    228
    ## 2                             106                                      1
    ## 3                             296                                    553
    ## 4                              77                                      0
    ## 5                               6                                      0
    ## 6                              31                                     31
    ##   y140.w119.x_female_pituitary_inc.d9 y146.r32.x_female_gonad_inc.prolong
    ## 1                              203.00                                  44
    ## 2                               77.00                                  71
    ## 3                              679.00                                 335
    ## 4                                0.00                                   6
    ## 5                                0.00                                   2
    ## 6                              111.51                                  59
    ##   y146.r32.x_female_hypothalamus_inc.prolong
    ## 1                                   243.0000
    ## 2                                     2.0000
    ## 3                                   717.0000
    ## 4                                     0.0000
    ## 5                                     1.0000
    ## 6                                    85.9667
    ##   y146.r32.x_female_pituitary_inc.prolong
    ## 1                                 65.0000
    ## 2                                 31.0000
    ## 3                                250.0000
    ## 4                                  0.0000
    ## 5                                  2.0000
    ## 6                                 94.4045
    ##   y148.w81.g145_male_gonad_m.inc.d17
    ## 1                            58.0000
    ## 2                             2.0000
    ## 3                          1702.0000
    ## 4                             7.0000
    ## 5                            14.8705
    ## 6                           795.8680
    ##   y148.w81.g145_male_hypothalamus_m.inc.d17
    ## 1                                   69.0000
    ## 2                                    2.0000
    ## 3                                  230.0000
    ## 4                                    0.0000
    ## 5                                    0.0000
    ## 6                                   16.0543
    ##   y148.w81.g145_male_pituitary_m.inc.d17 y149.r52.x_male_gonad_inc.d3
    ## 1                                547.000                            8
    ## 2                                153.000                            4
    ## 3                               1001.000                          444
    ## 4                                  0.000                            0
    ## 5                                  1.000                            2
    ## 6                                389.522                          137
    ##   y149.r52.x_male_hypothalamus_inc.d3 y149.r52.x_male_pituitary_inc.d3
    ## 1                             407.000                              135
    ## 2                               3.000                               11
    ## 3                             797.000                              222
    ## 4                               0.000                                0
    ## 5                               1.000                                3
    ## 6                              74.137                               55
    ##   y15.x_female_gonad_hatch y15.x_female_hypothalamus_hatch
    ## 1                       63                             307
    ## 2                       93                               7
    ## 3                      353                             581
    ## 4                        3                               0
    ## 5                        4                               3
    ## 6                       42                              68
    ##   y15.x_female_pituitary_hatch y18.x_male_gonad_m.inc.d3
    ## 1                      77.2156                   44.0000
    ## 2                      17.0000                    7.0000
    ## 3                     213.0000                 1732.0000
    ## 4                       0.0000                    0.0000
    ## 5                       3.0000                   15.8283
    ## 6                      59.7258                  897.8130
    ##   y18.x_male_hypothalamus_m.inc.d3 y18.x_male_pituitary_m.inc.d3
    ## 1                         128.0000                           480
    ## 2                           1.0000                            66
    ## 3                         357.0000                           725
    ## 4                           0.0000                             0
    ## 5                           0.0000                             5
    ## 6                          19.3011                           153
    ##   y4.x_male_gonad_m.inc.d17 y4.x_male_hypothalamus_m.inc.d17
    ## 1                         3                         127.0000
    ## 2                         0                           1.0000
    ## 3                       514                         253.0000
    ## 4                         0                           0.0000
    ## 5                         1                           0.0000
    ## 6                       134                          19.1706
    ##   y4.x_male_pituitary_m.inc.d17 y55.x_male_gonad_m.inc.d8
    ## 1                      100.0000                    46.000
    ## 2                       14.0000                     0.000
    ## 3                      191.0000                  1756.000
    ## 4                        0.0000                     0.000
    ## 5                        0.0000                    20.000
    ## 6                       67.7619                   496.855
    ##   y55.x_male_hypothalamus_m.inc.d8 y55.x_male_pituitary_m.inc.d8
    ## 1                         228.2791                       182.250
    ## 2                           3.0000                        90.000
    ## 3                         733.0000                       487.000
    ## 4                           0.0000                         0.000
    ## 5                           1.0000                         1.000
    ## 6                          46.1482                       142.488
    ##   y6.o54_female_gonad_n5 y6.o54_female_hypothalamus_n5
    ## 1                    198                      214.0000
    ## 2                    277                        5.0000
    ## 3                    674                      635.0000
    ## 4                     20                        0.0000
    ## 5                     10                        0.0000
    ## 6                    130                       56.2059
    ##   y6.o54_female_pituitary_n5 y63.x_male_gonad_m.inc.d9
    ## 1                    272.000                         7
    ## 2                     39.000                         3
    ## 3                    599.000                       497
    ## 4                      0.000                         0
    ## 5                      2.000                         1
    ## 6                    310.354                       158
    ##   y63.x_male_hypothalamus_m.inc.d9 y63.x_male_pituitary_m.inc.d9
    ## 1                          67.6528                      124.0000
    ## 2                           1.0000                       27.0000
    ## 3                         294.0000                      361.0000
    ## 4                           0.0000                        2.0000
    ## 5                           0.0000                        0.0000
    ## 6                           9.0000                       68.5814
    ##   y7.g58_female_gonad_hatch y7.g58_female_hypothalamus_hatch
    ## 1                        32                           241.00
    ## 2                        68                             2.00
    ## 3                       178                           516.00
    ## 4                         3                             0.00
    ## 5                         1                             3.00
    ## 6                        25                            49.09
    ##   y7.g58_female_pituitary_hatch y85.r71.x_female_gonad_m.inc.d17
    ## 1                       187.000                               44
    ## 2                        34.000                               26
    ## 3                       511.000                              332
    ## 4                         0.000                                9
    ## 5                         3.000                                2
    ## 6                       130.283                               46
    ##   y85.r71.x_female_hypothalamus_m.inc.d17
    ## 1                                313.0000
    ## 2                                  4.0000
    ## 3                                607.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 62.5207
    ##   y85.r71.x_female_pituitary_m.inc.d17 y94.g133.x_female_gonad_n5
    ## 1                                  122                         72
    ## 2                                   28                         61
    ## 3                                  273                        178
    ## 4                                    0                         13
    ## 5                                    1                          1
    ## 6                                   82                         27
    ##   y94.g133.x_female_hypothalamus_n5.NYNO y94.g133.x_female_pituitary_n5
    ## 1                               175.0000                        300.000
    ## 2                                 3.0000                         44.000
    ## 3                               528.0000                        612.000
    ## 4                                 0.0000                          0.000
    ## 5                                 1.0000                          2.000
    ## 6                                31.4769                        142.049
    ##   y95.g131.x_male_gonad_inc.d9 y95.g131.x_male_hypothalamus_inc.d9
    ## 1                        5.000                            743.0000
    ## 2                        1.000                             10.0000
    ## 3                      560.000                           1540.0000
    ## 4                        0.000                              0.0000
    ## 5                        6.000                              1.0000
    ## 6                      391.944                             87.5236
    ##   y95.g131.x_male_pituitary_inc.d9 y97.x_female_gonad_n9
    ## 1                               91             455.09400
    ## 2                                7             422.00000
    ## 3                              145             717.00000
    ## 4                                0               8.00000
    ## 5                                0              13.64535
    ## 6                               42             185.00000
    ##   y97.x_female_hypothalamus_n9 y97.x_female_pituitary_n9
    ## 1                     88.00000                   187.000
    ## 2                      0.00000                    75.000
    ## 3                    242.00000                   536.000
    ## 4                      0.00000                     0.000
    ## 5                      0.00000                     3.000
    ## 6                      6.04743                   129.126
    ##   y98.g54_female_gonad_m.hatch y98.g54_female_hypothalamus_m.hatch
    ## 1                           26                            241.0000
    ## 2                           35                             11.0000
    ## 3                          380                            564.0000
    ## 4                            1                              0.0000
    ## 5                           12                              1.0000
    ## 6                          109                             43.7713
    ##   y98.g54_female_pituitary_m.hatch y98.o50.x_male_gonad_inc.d3
    ## 1                          270.000                           7
    ## 2                          222.000                           1
    ## 3                          378.000                         462
    ## 4                            0.000                           0
    ## 5                            0.000                           3
    ## 6                          111.913                         174
    ##   y98.o50.x_male_hypothalamus_inc.d3 y98.o50.x_male_pituitary_inc.d3
    ## 1                           245.0000                        130.0000
    ## 2                             3.0000                         11.0000
    ## 3                           567.0000                        223.0000
    ## 4                             0.0000                          0.0000
    ## 5                             2.0000                          2.0000
    ## 6                            46.1529                         67.6398

    ## save gene information
    geneinfo <- kallistodata %>%
      select(gene, geneid, NCBI) 
    head(geneinfo)

    ##       gene geneid           NCBI
    ## 1    EDNRB 408082 NP_001001127.1
    ## 2  CYP26A1 408183 NP_001001129.1
    ## 3    CFDP1 374073 NP_001001189.1
    ## 4    AvBD7 407777 NP_001001194.1
    ## 5     KRT5 407779 NP_001001195.1
    ## 6 HSD11B1L 408034 NP_001001201.1

\# check for genes that have multple transcripts expressed
----------------------------------------------------------

    isoforms <- kallistodata %>%
      group_by(gene) %>%
      summarize(n = n()) %>%
      filter(n > 1) %>%
      group_by(n) %>% 
      summarize(genes = str_c(gene, collapse = ", ")) %>%
      arrange(desc(n)) %>%
      rename("n(isoforms)" = "n") %>% 
      mutate("n(counts)" = str_count(genes, ",") + 1)
    head(isoforms)

    ## # A tibble: 6 x 3
    ##   `n(isoforms)` genes                                           `n(counts)`
    ##           <int> <chr>                                                 <dbl>
    ## 1            57 SSPO                                                      1
    ## 2             6 CEND1, LOC101747901, TIRAP                                3
    ## 3             5 AIP, CCL5, GINS3, SORBS2, TCF7L2                          5
    ## 4             4 ATG13, BPTF, DOCK3, DTNA, EIF4G3, ELAVL2, EPB4…          14
    ## 5             3 AGAP3, ARFIP1, ARHGAP17, BZRAP1, C2CD5, CACNA1…          70
    ## 6             2 ABCB10, ABCC3, ABCC6, ABCC9, ABCG8, ABI1, ABI2…         697

    # for example, these gene have 2 and 3 isoforms
    geneinfo %>% filter(gene == "GRIN1")

    ##    gene geneid           NCBI
    ## 1 GRIN1 404296    NP_996862.1
    ## 2 GRIN1 404296 XP_015134834.1

    geneinfo %>% filter(gene == "CACNA1C")

    ##      gene geneid           NCBI
    ## 1 CACNA1C 395891 XP_015141927.1
    ## 2 CACNA1C 395891 XP_015141993.1
    ## 3 CACNA1C 395891 XP_015142129.1

group data by gene
------------------

    head(kallistodata)

    ##   row.names     gene geneid           NCBI L.Blu13_male_gonad_control.NYNO
    ## 1    408082    EDNRB 408082 NP_001001127.1                               6
    ## 2    408183  CYP26A1 408183 NP_001001129.1                               1
    ## 3    374073    CFDP1 374073 NP_001001189.1                             408
    ## 4    407777    AvBD7 407777 NP_001001194.1                               1
    ## 5    407779     KRT5 407779 NP_001001195.1                               0
    ## 6    408034 HSD11B1L 408034 NP_001001201.1                              65
    ##   L.Blu13_male_hypothalamus_control.NYNO
    ## 1                               37.00000
    ## 2                                0.00000
    ## 3                               72.00000
    ## 4                                0.00000
    ## 5                                0.00000
    ## 6                                2.72853
    ##   L.Blu13_male_pituitary_control.NYNO L.G107_male_gonad_control
    ## 1                                  46                     7.000
    ## 2                                  47                     2.000
    ## 3                                 197                   728.000
    ## 4                                   0                     1.000
    ## 5                                   0                     2.000
    ## 6                                  30                   146.948
    ##   L.G107_male_hypothalamus_control L.G107_male_pituitary_control
    ## 1                               43                       60.0000
    ## 2                                1                       25.0000
    ## 3                              159                      186.0000
    ## 4                                0                        0.0000
    ## 5                                0                        0.0000
    ## 6                                8                       18.5133
    ##   L.G118_female_gonad_control L.G118_female_hypothalamus_control.NYNO
    ## 1                          33                                      60
    ## 2                          94                                       0
    ## 3                         733                                      66
    ## 4                          96                                       0
    ## 5                           2                                       0
    ## 6                          53                                       3
    ##   L.G118_female_pituitary_control.NYNO L.R3_male_gonad_control.NYNO
    ## 1                             45.00000                            4
    ## 2                             14.00000                            0
    ## 3                            298.00000                          273
    ## 4                              0.00000                            0
    ## 5                              1.05712                            0
    ## 6                             38.00000                           40
    ##   L.R3_male_hypothalamus_control L.R3_male_pituitary_control.NYNO
    ## 1                            113                              114
    ## 2                              0                               37
    ## 3                            244                              249
    ## 4                              0                                0
    ## 5                              0                                0
    ## 6                             15                               26
    ##   L.R8_male_gonad_control L.R8_male_hypothalamus_control
    ## 1                  15.000                              6
    ## 2                   1.000                              0
    ## 3                1814.000                             23
    ## 4                   0.000                              0
    ## 5                   5.000                              0
    ## 6                 241.948                              0
    ##   L.R8_male_pituitary_control L.W33_male_gonad_control
    ## 1                          87                        3
    ## 2                          58                        1
    ## 3                         226                      301
    ## 4                           0                        0
    ## 5                           1                        1
    ## 6                          44                       37
    ##   L.W33_male_hypothalamus_control.NYNO L.W33_male_pituitary_control
    ## 1                            112.00000                     117.0000
    ## 2                              0.00000                      54.0000
    ## 3                            200.00000                     252.0000
    ## 4                              0.00000                       0.0000
    ## 5                              0.00000                       1.0000
    ## 6                              7.91449                      27.7011
    ##   L.W3_male_gonad_control.NYNO L.W3_male_hypothalamus_control
    ## 1                            1                       11.00000
    ## 2                            0                        0.00000
    ## 3                          247                       24.00000
    ## 4                            1                        0.00000
    ## 5                            1                        0.00000
    ## 6                           38                        1.13262
    ##   L.W3_male_pituitary_control L.W4_male_gonad_control.NYNO
    ## 1                          71                           10
    ## 2                          33                            1
    ## 3                         301                          345
    ## 4                           0                            0
    ## 5                           1                            1
    ## 6                          49                           57
    ##   L.W4_male_hypothalamus_control L.W4_male_pituitary_control
    ## 1                       111.0000                          68
    ## 2                         0.0000                           7
    ## 3                       427.0000                         246
    ## 4                         0.0000                           0
    ## 5                         0.0000                           0
    ## 6                        14.1608                          39
    ##   R.G106_female_gonad_control R.G106_female_hypothalamus_control
    ## 1                          95                           180.0000
    ## 2                         175                             2.0000
    ## 3                         823                           499.0000
    ## 4                          58                             0.0000
    ## 5                           2                             0.0000
    ## 6                          60                            27.2655
    ##   R.G106_female_pituitary_control R.R20_female_gonad_control
    ## 1                         98.0000                         61
    ## 2                          4.0000                         71
    ## 3                        176.0000                        320
    ## 4                          0.0000                         54
    ## 5                          0.0000                          0
    ## 6                         40.9282                         25
    ##   R.R20_female_hypothalamus_control.NYNO R.R20_female_pituitary_control
    ## 1                                     53                             34
    ## 2                                      0                              7
    ## 3                                    204                            226
    ## 4                                      0                              1
    ## 5                                      0                              2
    ## 6                                      6                             22
    ##   R.R9_female_gonad_control R.R9_female_hypothalamus_control
    ## 1                  39.10632                         270.0000
    ## 2                  76.00000                           1.0000
    ## 3                 327.00000                         448.0000
    ## 4                  43.00000                           0.0000
    ## 5                   1.00000                           2.0000
    ## 6                  12.20620                          16.1167
    ##   R.R9_female_pituitary_control.NYNO R.W44_female_gonad_control
    ## 1                             120.00                  159.40370
    ## 2                               2.00                   63.00000
    ## 3                             264.00                  775.00000
    ## 4                               0.00                  100.00000
    ## 5                               1.00                    6.60068
    ## 6                              17.21                   55.00000
    ##   R.W44_female_hypothalamus_control R.W44_female_pituitary_control.NYNO
    ## 1                          540.0000                                  48
    ## 2                            4.0000                                  33
    ## 3                          886.0000                                 141
    ## 4                            1.0000                                   0
    ## 5                            0.0000                                   1
    ## 6                           34.4786                                  19
    ##   R.Y108.W29_male_gonad_control R.Y108.W29_male_hypothalamus_control.NYNO
    ## 1                            10                                   49.0000
    ## 2                             0                                    1.0000
    ## 3                           954                                   69.0000
    ## 4                             0                                    0.0000
    ## 5                             5                                    0.0000
    ## 6                           153                                    6.0712
    ##   R.Y108.W29_male_pituitary_control blk.s030.o.g_male_gonad_prolong
    ## 1                                93                          52.000
    ## 2                                 7                           0.000
    ## 3                               256                        1127.000
    ## 4                                 0                           2.000
    ## 5                                 0                           7.000
    ## 6                                49                         349.972
    ##   blk.s030.o.g_male_hypothalamus_prolong
    ## 1                                77.0000
    ## 2                                 1.0000
    ## 3                               271.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                15.0411
    ##   blk.s030.o.g_male_pituitary_prolong blk.s031.pu.d_female_gonad_prolong
    ## 1                           306.00000                                268
    ## 2                            56.00000                                312
    ## 3                           550.00000                               1430
    ## 4                             0.00000                                 44
    ## 5                             3.70651                                 11
    ## 6                            96.27890                                178
    ##   blk.s031.pu.d_female_hypothalamus_prolong
    ## 1                                   86.0000
    ## 2                                    0.0000
    ## 3                                  253.0000
    ## 4                                    0.0000
    ## 5                                    0.0000
    ## 6                                   14.3261
    ##   blk.s031.pu.d_female_pituitary_prolong blk.s032.g.w_female_gonad_m.hatch
    ## 1                                412.000                            47.000
    ## 2                                134.000                           248.000
    ## 3                                417.000                          1017.000
    ## 4                                  1.000                            67.000
    ## 5                                  1.000                             4.000
    ## 6                                185.615                            95.941
    ##   blk.s032.g.w_female_hypothalamus_m.hatch
    ## 1                                 268.0000
    ## 2                                   2.0000
    ## 3                                 549.0000
    ## 4                                   2.0000
    ## 5                                   1.0000
    ## 6                                  43.4202
    ##   blk.s032.g.w_female_pituitary_m.hatch blk.s049.y.g_female_gonad_m.inc.d3
    ## 1                               336.000                                487
    ## 2                                12.000                                828
    ## 3                               669.000                               1609
    ## 4                                 0.000                                264
    ## 5                                 3.000                                 21
    ## 6                               288.056                                172
    ##   blk.s049.y.g_female_hypothalamus_m.inc.d3
    ## 1                                   78.0000
    ## 2                                    2.0000
    ## 3                                  241.0000
    ## 4                                    0.0000
    ## 5                                    0.0000
    ## 6                                   13.1815
    ##   blk.s049.y.g_female_pituitary_m.inc.d3
    ## 1                                267.000
    ## 2                                129.000
    ## 3                                494.000
    ## 4                                  0.000
    ## 5                                  1.000
    ## 6                                120.259
    ##   blk.s060.pu.w_female_gonad_m.inc.d3
    ## 1                                  68
    ## 2                                 149
    ## 3                                1079
    ## 4                                   7
    ## 5                                   7
    ## 6                                 106
    ##   blk.s060.pu.w_female_hypothalamus_m.inc.d3
    ## 1                                   215.4358
    ## 2                                     3.0000
    ## 3                                   545.0000
    ## 4                                     0.0000
    ## 5                                     1.0000
    ## 6                                    40.1182
    ##   blk.s060.pu.w_female_pituitary_m.inc.d3.NYNO
    ## 1                                      170.000
    ## 2                                      105.000
    ## 3                                      418.000
    ## 4                                        0.000
    ## 5                                        0.000
    ## 6                                      132.596
    ##   blk.s061.pu.y_female_gonad_inc.d9
    ## 1                                45
    ## 2                                77
    ## 3                               233
    ## 4                               154
    ## 5                                 1
    ## 6                                26
    ##   blk.s061.pu.y_female_hypothalamus_inc.d9
    ## 1                                 391.0000
    ## 2                                   4.0000
    ## 3                                 583.0000
    ## 4                                   0.0000
    ## 5                                   2.0000
    ## 6                                  55.2234
    ##   blk.s061.pu.y_female_pituitary_inc.d9 blk.y.l.s109_female_gonad_m.inc.d8
    ## 1                               85.2333                           190.8635
    ## 2                               61.0000                           329.0000
    ## 3                              189.0000                           691.0000
    ## 4                                0.0000                           151.0000
    ## 5                                0.0000                             5.0000
    ## 6                               54.9364                            75.0000
    ##   blk.y.l.s109_female_hypothalamus_m.inc.d8
    ## 1                                  285.0000
    ## 2                                    2.0000
    ## 3                                  528.0000
    ## 4                                    0.0000
    ## 5                                    1.0000
    ## 6                                   31.6834
    ##   blk.y.l.s109_female_pituitary_m.inc.d8 blk0.x_female_gonad_m.n2
    ## 1                              323.00000                       32
    ## 2                               27.00000                       38
    ## 3                              663.00000                      263
    ## 4                                0.00000                        1
    ## 5                                7.09367                        0
    ## 6                              138.42200                       37
    ##   blk0.x_female_hypothalamus_m.n2 blk0.x_female_pituitary_m.n2
    ## 1                         169.000                      89.0000
    ## 2                           1.000                      24.0000
    ## 3                         234.000                     261.0000
    ## 4                           0.000                       0.0000
    ## 5                           0.000                       1.0000
    ## 6                          17.462                      69.8865
    ##   blk11.x_female_gonad_bldg blk11.x_female_hypothalamus_bldg
    ## 1                        50                         334.0000
    ## 2                        32                           3.0000
    ## 3                       222                         594.0000
    ## 4                        38                           0.0000
    ## 5                         0                           1.0000
    ## 6                        22                          38.1538
    ##   blk11.x_female_pituitary_bldg blk12.x_male_gonad_n5
    ## 1                        206.00                 6.000
    ## 2                        137.00                 0.000
    ## 3                        701.00               355.000
    ## 4                          0.00                 0.000
    ## 5                          3.00                 1.000
    ## 6                        135.76               132.874
    ##   blk12.x_male_hypothalamus_n5.NYNO blk12.x_male_pituitary_n5
    ## 1                          171.0000                   256.000
    ## 2                            1.0000                    82.000
    ## 3                          567.0000                   757.000
    ## 4                            0.0000                     0.000
    ## 5                            0.0000                     1.000
    ## 6                           42.2418                   254.513
    ##   blk17.x_male_gonad_inc.d17 blk17.x_male_hypothalamus_inc.d17
    ## 1                          9                          249.0000
    ## 2                          3                            4.0000
    ## 3                        401                          681.0000
    ## 4                          0                            1.0000
    ## 5                          5                            2.0000
    ## 6                         83                           56.1321
    ##   blk17.x_male_pituitary_inc.d17 blk19.x_female_gonad_extend
    ## 1                        293.000                   171.00000
    ## 2                         58.000                   161.00000
    ## 3                        703.000                   510.00000
    ## 4                          0.000                    21.00000
    ## 5                          1.000                     8.32329
    ## 6                        162.923                    57.00000
    ##   blk19.x_female_hypothalamus_extend blk19.x_female_pituitary_extend
    ## 1                           245.0000                         123.000
    ## 2                             2.0000                           8.000
    ## 3                           746.0000                         286.000
    ## 4                             0.0000                           0.000
    ## 5                             1.0000                           4.000
    ## 6                            91.5648                         121.688
    ##   blk21.x_female_gonad_hatch blk21.x_female_hypothalamus_hatch
    ## 1                         39                               282
    ## 2                        135                                 0
    ## 3                        324                              1038
    ## 4                          2                                 0
    ## 5                          2                                 0
    ## 6                         27                                93
    ##   blk21.x_female_pituitary_hatch blk4.x_female_gonad_n9
    ## 1                        92.5414                     74
    ## 2                        16.0000                      9
    ## 3                       231.0000                    245
    ## 4                         0.0000                      1
    ## 5                         0.0000                      2
    ## 6                        70.8845                     95
    ##   blk4.x_female_hypothalamus_n9 blk4.x_female_pituitary_n9
    ## 1                       112.000                  183.00000
    ## 2                         3.000                   23.00000
    ## 3                      1148.000                  859.00000
    ## 4                         0.000                    0.00000
    ## 5                         0.000                    3.80532
    ## 6                        45.113                  206.00000
    ##   blk5.x_male_gonad_m.inc.d3 blk5.x_male_hypothalamus_m.inc.d3
    ## 1                     34.000                           95.0000
    ## 2                      4.000                            1.0000
    ## 3                   1375.000                          241.0000
    ## 4                      0.000                            0.0000
    ## 5                      4.000                            0.0000
    ## 6                    735.984                           12.0309
    ##   blk5.x_male_pituitary_m.inc.d3 blu.o.x.ATLAS_female_gonad_control
    ## 1                        270.000                                 45
    ## 2                         72.000                                123
    ## 3                        522.000                                602
    ## 4                          0.000                                  5
    ## 5                          1.000                                  2
    ## 6                        157.821                                 32
    ##   blu.o.x.ATLAS_female_hypothalamus_control
    ## 1                                        22
    ## 2                                         0
    ## 3                                        60
    ## 4                                         0
    ## 5                                         0
    ## 6                                         4
    ##   blu.o.x.ATLAS_female_pituitary_control blu10.w26.x_male_gonad_m.hatch
    ## 1                                142.000                       58.00000
    ## 2                                 34.000                        1.00000
    ## 3                                450.000                      680.00000
    ## 4                                  0.000                        5.00000
    ## 5                                  1.000                        9.37313
    ## 6                                 57.062                      342.00000
    ##   blu10.w26.x_male_hypothalamus_m.hatch blu10.w26.x_male_pituitary_m.hatch
    ## 1                              251.0000                             170.00
    ## 2                                3.0000                              59.00
    ## 3                              546.0000                             353.00
    ## 4                                0.0000                               0.00
    ## 5                                0.0000                               0.00
    ## 6                               42.4404                             102.69
    ##   blu103.x_female_gonad_hatch.NYNO blu103.x_female_hypothalamus_hatch
    ## 1                          101.000                            361.000
    ## 2                           35.000                              2.000
    ## 3                          805.000                            661.000
    ## 4                           91.000                              0.000
    ## 5                            6.000                              0.000
    ## 6                          171.977                             56.517
    ##   blu103.x_female_pituitary_hatch.NYNO blu104.w120.x_male_gonad_hatch
    ## 1                               234.00                         6.0000
    ## 2                                43.00                         0.0000
    ## 3                               482.00                       410.0000
    ## 4                                 0.00                         0.0000
    ## 5                                 0.00                         2.0000
    ## 6                               158.13                        84.8101
    ##   blu104.w120.x_male_hypothalamus_hatch
    ## 1                                   461
    ## 2                                     7
    ## 3                                   654
    ## 4                                     0
    ## 5                                     2
    ## 6                                    58
    ##   blu104.w120.x_male_pituitary_hatch.NYNO
    ## 1                                      96
    ## 2                                      19
    ## 3                                     426
    ## 4                                       0
    ## 5                                       0
    ## 6                                      84
    ##   blu108.w40.o158_male_gonad_inc.d9
    ## 1                            11.000
    ## 2                             0.000
    ## 3                           417.000
    ## 4                             0.000
    ## 5                             2.000
    ## 6                           237.985
    ##   blu108.w40.o158_male_hypothalamus_inc.d9
    ## 1                                  452.000
    ## 2                                    4.000
    ## 3                                  636.000
    ## 4                                    0.000
    ## 5                                    1.000
    ## 6                                   33.126
    ##   blu108.w40.o158_male_pituitary_inc.d9 blu111.w113.x_male_gonad_inc.d3
    ## 1                              108.1732                              11
    ## 2                                9.0000                               0
    ## 3                              203.0000                             415
    ## 4                                0.0000                               0
    ## 5                                1.0000                               0
    ## 6                               72.5762                             121
    ##   blu111.w113.x_male_hypothalamus_inc.d3
    ## 1                               163.0000
    ## 2                                 0.0000
    ## 3                               370.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                23.2068
    ##   blu111.w113.x_male_pituitary_inc.d3 blu113.w124.x_male_gonad_inc.d17
    ## 1                                 155                                6
    ## 2                                  68                                1
    ## 3                                 889                              380
    ## 4                                   0                                0
    ## 5                                   4                                0
    ## 6                                 149                              102
    ##   blu113.w124.x_male_hypothalamus_inc.d17
    ## 1                                728.0000
    ## 2                                  2.0000
    ## 3                                822.0000
    ## 4                                  0.0000
    ## 5                                  5.0000
    ## 6                                 46.7947
    ##   blu113.w124.x_male_pituitary_inc.d17.NYNO
    ## 1                                   219.000
    ## 2                                    22.000
    ## 3                                   535.000
    ## 4                                     0.000
    ## 5                                     4.000
    ## 6                                   140.482
    ##   blu114.r38.w198_male_gonad_bldg blu114.r38.w198_male_hypothalamus_bldg
    ## 1                           3.000                               312.6289
    ## 2                           0.000                                 4.0000
    ## 3                         527.000                               461.0000
    ## 4                           2.000                                 0.0000
    ## 5                           1.000                                 0.0000
    ## 6                         153.888                                30.1273
    ##   blu114.r38.w198_male_pituitary_bldg
    ## 1                                  95
    ## 2                                  93
    ## 3                                 281
    ## 4                                   0
    ## 5                                   0
    ## 6                                 109
    ##   blu115.y150.x_female_gonad_inc.prolong
    ## 1                                     31
    ## 2                                     48
    ## 3                                    331
    ## 4                                      4
    ## 5                                      1
    ## 6                                     66
    ##   blu115.y150.x_female_hypothalamus_inc.prolong
    ## 1                                      239.0000
    ## 2                                        6.0000
    ## 3                                      596.0000
    ## 4                                        1.0000
    ## 5                                        1.0000
    ## 6                                       47.3668
    ##   blu115.y150.x_female_pituitary_inc.prolong
    ## 1                                         71
    ## 2                                          9
    ## 3                                        206
    ## 4                                          0
    ## 5                                          0
    ## 6                                         61
    ##   blu119.w84.x_female_gonad_m.inc.d8
    ## 1                                 95
    ## 2                                 38
    ## 3                                843
    ## 4                                 57
    ## 5                                  2
    ## 6                                 94
    ##   blu119.w84.x_female_hypothalamus_m.inc.d8
    ## 1                                  314.0000
    ## 2                                    9.0000
    ## 3                                  999.0000
    ## 4                                    0.0000
    ## 5                                    4.0000
    ## 6                                   89.3176
    ##   blu119.w84.x_female_pituitary_m.inc.d8 blu121.w91.x_male_gonad_inc.d17
    ## 1                                221.000                           8.000
    ## 2                                 42.000                           0.000
    ## 3                                461.000                         395.000
    ## 4                                  0.000                           1.000
    ## 5                                  0.000                           1.000
    ## 6                                184.234                         117.874
    ##   blu121.w91.x_male_hypothalamus_inc.d17
    ## 1                               467.0000
    ## 2                                 1.0000
    ## 3                               893.0000
    ## 4                                 0.0000
    ## 5                                 3.0000
    ## 6                                66.1144
    ##   blu121.w91.x_male_pituitary_inc.d17 blu124.w180.x_female_gonad_hatch
    ## 1                             74.0000                               53
    ## 2                              5.0000                               70
    ## 3                            229.0000                              245
    ## 4                              0.0000                               15
    ## 5                              0.0000                                0
    ## 6                             45.4034                               28
    ##   blu124.w180.x_female_hypothalamus_hatch
    ## 1                                 471.000
    ## 2                                   0.000
    ## 3                                 714.000
    ## 4                                   0.000
    ## 5                                   2.000
    ## 6                                  51.193
    ##   blu124.w180.x_female_pituitary_hatch blu33.y88.x_male_gonad_bldg
    ## 1                              69.0000                     8.09308
    ## 2                               7.0000                     3.00000
    ## 3                             165.0000                   631.00000
    ## 4                               0.0000                     0.00000
    ## 5                               0.0000                     6.00000
    ## 6                              74.0788                   243.76100
    ##   blu33.y88.x_male_hypothalamus_bldg blu33.y88.x_male_pituitary_bldg
    ## 1                                179                         474.000
    ## 2                                  1                          34.000
    ## 3                                425                         528.000
    ## 4                                  0                           0.000
    ## 5                                  0                           2.000
    ## 6                                 29                         130.845
    ##   blu36.w16_female_gonad_n9 blu36.w16_female_hypothalamus_n9
    ## 1                        64                         109.0000
    ## 2                        72                           1.0000
    ## 3                       373                         592.0000
    ## 4                        10                           0.0000
    ## 5                         1                           1.0000
    ## 6                        69                          20.2337
    ##   blu36.w16_female_pituitary_n9 blu37.r65.x_male_gonad_n5
    ## 1                         245.0                     5.000
    ## 2                          67.0                     0.000
    ## 3                         739.0                   355.000
    ## 4                           0.0                     0.000
    ## 5                           6.0                     0.000
    ## 6                         171.6                   131.811
    ##   blu37.r65.x_male_hypothalamus_n5 blu37.r65.x_male_pituitary_n5
    ## 1                         291.0000                       223.000
    ## 2                           1.0000                        45.000
    ## 3                         753.0000                       620.000
    ## 4                           0.0000                         0.000
    ## 5                           0.0000                         6.000
    ## 6                          72.7684                       184.411
    ##   blu38.g135.x_female_gonad_bldg blu38.g135.x_female_hypothalamus_bldg
    ## 1                              6                              256.0000
    ## 2                              3                                4.0000
    ## 3                            262                              450.0000
    ## 4                              0                                0.0000
    ## 5                              1                                0.0000
    ## 6                             83                               46.1741
    ##   blu38.g135.x_female_pituitary_bldg blu39.o26.x_female_gonad_inc.d3
    ## 1                                 63                              31
    ## 2                                 25                              10
    ## 3                                239                             360
    ## 4                                  0                              29
    ## 5                                  3                               0
    ## 6                                 88                              15
    ##   blu39.o26.x_female_hypothalamus_inc.d3.NYNO
    ## 1                                    305.0000
    ## 2                                      4.0000
    ## 3                                    516.0000
    ## 4                                      0.0000
    ## 5                                      0.0000
    ## 6                                     28.4879
    ##   blu39.o26.x_female_pituitary_inc.d3.NYNO blu41.y100.x_male_gonad_n5
    ## 1                                263.00000                          2
    ## 2                                 52.00000                          1
    ## 3                                976.00000                        479
    ## 4                                  0.00000                          1
    ## 5                                  2.02492                          2
    ## 6                                316.27000                        104
    ##   blu41.y100.x_male_hypothalamus_n5.NYNO blu41.y100.x_male_pituitary_n5
    ## 1                               149.0000                        261.000
    ## 2                                 4.0000                         31.000
    ## 3                               451.0000                        594.000
    ## 4                                 0.0000                          0.000
    ## 5                                 1.0000                          2.000
    ## 6                                41.3472                        214.119
    ##   blu44.y102_female_gonad_extend blu44.y102_female_hypothalamus_extend
    ## 1                      101.00000                              471.0000
    ## 2                       96.00000                                6.0000
    ## 3                      629.00000                              919.0000
    ## 4                        7.00000                                1.0000
    ## 5                        7.83441                                0.0000
    ## 6                      177.00000                               54.1446
    ##   blu44.y102_female_pituitary_extend blu47.y96.x_female_gonad_inc.d9
    ## 1                            268.000                              73
    ## 2                             29.000                             120
    ## 3                            356.000                             329
    ## 4                              0.000                              43
    ## 5                              1.000                               3
    ## 6                            111.138                              47
    ##   blu47.y96.x_female_hypothalamus_inc.d9
    ## 1                                     93
    ## 2                                      0
    ## 3                                    247
    ## 4                                      0
    ## 5                                      0
    ## 6                                     10
    ##   blu47.y96.x_female_pituitary_inc.d9 blu55.g51_female_gonad_n5
    ## 1                                  85                        19
    ## 2                                  37                        13
    ## 3                                 195                       251
    ## 4                                   0                         7
    ## 5                                   2                         1
    ## 6                                  64                        41
    ##   blu55.g51_female_hypothalamus_n5 blu55.g51_female_pituitary_n5
    ## 1                         461.4700                       62.0000
    ## 2                           3.0000                       11.0000
    ## 3                         772.0000                      177.0000
    ## 4                           0.0000                        0.0000
    ## 5                           5.0000                        2.0000
    ## 6                          73.6446                       42.7425
    ##   blu56.o53_female_gonad_m.inc.d3 blu56.o53_female_hypothalamus_m.inc.d3
    ## 1                              97                               118.0000
    ## 2                             137                                 1.0000
    ## 3                             995                               367.0000
    ## 4                              79                                 0.0000
    ## 5                               1                                 0.0000
    ## 6                              93                                26.1381
    ##   blu56.o53_female_pituitary_m.inc.d3 blu63.g62_female_gonad_m.inc.d9
    ## 1                             260.000                         29.0000
    ## 2                              53.000                         27.0000
    ## 3                             752.000                        254.0000
    ## 4                               0.000                         20.0000
    ## 5                               2.000                          0.0000
    ## 6                             276.256                         58.1971
    ##   blu63.g62_female_hypothalamus_m.inc.d9
    ## 1                               93.00000
    ## 2                                1.00000
    ## 3                              237.00000
    ## 4                                0.00000
    ## 5                                0.00000
    ## 6                                8.12378
    ##   blu63.g62_female_pituitary_m.inc.d9 blu80.r97_female_gonad_m.inc.d8
    ## 1                                  74                              67
    ## 2                                  13                             206
    ## 3                                 214                            1083
    ## 4                                   0                              16
    ## 5                                   1                               4
    ## 6                                  50                             105
    ##   blu80.r97_female_hypothalamus_m.inc.d8
    ## 1                               243.0000
    ## 2                                 0.0000
    ## 3                               588.0000
    ## 4                                 0.0000
    ## 5                                 1.0000
    ## 6                                38.1029
    ##   blu80.r97_female_pituitary_m.inc.d8 blu81.r88_male_gonad_n9
    ## 1                             164.000                 33.0000
    ## 2                              23.000                  5.0000
    ## 3                             610.000               1516.0000
    ## 4                               0.000                  0.0000
    ## 5                               1.000                  6.1045
    ## 6                             155.062                524.9910
    ##   blu81.r88_male_hypothalamus_n9 blu81.r88_male_pituitary_n9
    ## 1                       210.0000                          74
    ## 2                         6.0000                          19
    ## 3                       792.0000                         300
    ## 4                         0.0000                           0
    ## 5                         0.0000                           0
    ## 6                        57.2532                          46
    ##   blu84.x_male_gonad_extend.hatch blu84.x_male_hypothalamus_extend.hatch
    ## 1                           17.00                               225.0000
    ## 2                            0.00                                 3.0000
    ## 3                          184.00                               432.0000
    ## 4                           18.00                                 0.0000
    ## 5                            0.00                                 1.0000
    ## 6                          272.73                                31.2052
    ##   blu84.x_male_pituitary_extend.hatch d.r.blk.s159_female_gonad_m.inc.d9
    ## 1                             556.000                          238.14460
    ## 2                              72.000                          542.00000
    ## 3                             757.000                         1028.00000
    ## 4                               0.000                           22.00000
    ## 5                               1.000                           10.98146
    ## 6                             270.416                          119.96700
    ##   d.r.blk.s159_female_hypothalamus_m.inc.d9
    ## 1                                  401.0000
    ## 2                                    2.0000
    ## 3                                  893.0000
    ## 4                                    0.0000
    ## 5                                    2.0000
    ## 6                                   85.4804
    ##   d.r.blk.s159_female_pituitary_m.inc.d9 d.s008.y.blk_male_gonad_n5
    ## 1                               233.4187                          3
    ## 2                                38.0000                          0
    ## 3                               412.0000                        308
    ## 4                                 0.0000                          1
    ## 5                                 3.0000                          5
    ## 6                               120.8720                        100
    ##   d.s008.y.blk_male_hypothalamus_n5 d.s008.y.blk_male_pituitary_n5
    ## 1                               360                        91.0000
    ## 2                                 9                         3.0000
    ## 3                               777                       173.0000
    ## 4                                 0                         0.0000
    ## 5                                 1                         3.0000
    ## 6                                34                        32.3868
    ##   d.s047.blk.o_male_gonad_n5 d.s047.blk.o_male_hypothalamus_n5
    ## 1                     10.000                          270.0000
    ## 2                      4.000                            5.0000
    ## 3                   1368.000                          595.0000
    ## 4                      1.000                            0.0000
    ## 5                      2.000                            0.0000
    ## 6                    496.665                           35.3459
    ##   d.s047.blk.o_male_pituitary_n5 d.s110.g.blk_male_gonad_m.inc.d3
    ## 1                         85.000                           17.000
    ## 2                         41.000                            4.000
    ## 3                        202.000                         1328.000
    ## 4                          0.000                            2.000
    ## 5                          0.000                            3.000
    ## 6                         62.738                          498.963
    ##   d.s110.g.blk_male_hypothalamus_m.inc.d3
    ## 1                                     128
    ## 2                                       1
    ## 3                                     290
    ## 4                                       0
    ## 5                                       0
    ## 6                                      25
    ##   d.s110.g.blk_male_pituitary_m.inc.d3 d.s112.blk.w_female_gonad_m.inc.d17
    ## 1                              240.000                                  52
    ## 2                               41.000                                  94
    ## 3                              504.000                                 438
    ## 4                                0.000                                  35
    ## 5                                3.000                                  12
    ## 6                              183.331                                  89
    ##   d.s112.blk.w_female_hypothalamus_m.inc.d17
    ## 1                                    71.0000
    ## 2                                     3.0000
    ## 3                                   258.0000
    ## 4                                     0.0000
    ## 5                                     0.0000
    ## 6                                    15.4873
    ##   d.s112.blk.w_female_pituitary_m.inc.d17
    ## 1                                  385.00
    ## 2                                  198.00
    ## 3                                  858.00
    ## 4                                    0.00
    ## 5                                    6.00
    ## 6                                  273.03
    ##   d.s177.blk.r_female_gonad_m.inc.d3
    ## 1                                248
    ## 2                                447
    ## 3                               1095
    ## 4                                118
    ## 5                                  9
    ## 6                                248
    ##   d.s177.blk.r_female_hypothalamus_m.inc.d3
    ## 1                                  257.0000
    ## 2                                    3.0000
    ## 3                                  434.0000
    ## 4                                    0.0000
    ## 5                                    2.0000
    ## 6                                   51.2969
    ##   d.s177.blk.r_female_pituitary_m.inc.d3 g.blk.s004.pk_female_gonad_lay
    ## 1                                    305                            438
    ## 2                                    113                            112
    ## 3                                    463                            731
    ## 4                                      0                            262
    ## 5                                      0                              0
    ## 6                                    183                            120
    ##   g.blk.s004.pk_female_hypothalamus_lay g.blk.s004.pk_female_pituitary_lay
    ## 1                              230.8799                            182.000
    ## 2                                2.0000                             10.000
    ## 3                              778.0000                            521.000
    ## 4                                0.0000                              1.000
    ## 5                                2.0000                              1.000
    ## 6                               71.5386                            279.425
    ##   g.blk.s041.r_male_gonad_m.inc.d3 g.blk.s041.r_male_hypothalamus_m.inc.d3
    ## 1                               12                                216.0000
    ## 2                                2                                  2.0000
    ## 3                              853                                613.0000
    ## 4                                0                                  0.0000
    ## 5                                9                                  5.0000
    ## 6                              442                                 49.1011
    ##   g.blk.s041.r_male_pituitary_m.inc.d3 g.o.y.s037_male_gonad_m.inc.d17
    ## 1                              296.000                               7
    ## 2                               85.000                               2
    ## 3                              448.000                            1061
    ## 4                                0.000                               1
    ## 5                                2.000                               6
    ## 6                              173.263                             389
    ##   g.o.y.s037_male_hypothalamus_m.inc.d17
    ## 1                               210.0000
    ## 2                                 1.0000
    ## 3                               407.0000
    ## 4                                 0.0000
    ## 5                                 1.0000
    ## 6                                31.2125
    ##   g.o.y.s037_male_pituitary_m.inc.d17 g.s.blk.d_male_gonad_n9
    ## 1                              275.00                  21.000
    ## 2                               31.00                  12.000
    ## 3                              737.00                1244.000
    ## 4                                0.00                   1.000
    ## 5                                5.00                   9.000
    ## 6                              242.46                 556.928
    ##   g.s.blk.d_male_hypothalamus_n9 g.s.blk.d_male_pituitary_n9
    ## 1                       128.0000                     205.000
    ## 2                         2.0000                      28.000
    ## 3                       294.0000                     351.000
    ## 4                         0.0000                       0.000
    ## 5                         1.0000                       1.000
    ## 6                        20.1458                     110.423
    ##   g.s.blk.y_male_gonad_lay g.s.blk.y_male_hypothalamus_lay
    ## 1                   14.000                        110.0000
    ## 2                    4.000                          0.0000
    ## 3                  816.000                       1427.0000
    ## 4                    0.000                          0.0000
    ## 5                    0.000                          0.0000
    ## 6                  381.932                         27.1544
    ##   g.s.blk.y_male_pituitary_lay g.s043.pu.blk_male_gonad_lay
    ## 1                      245.000                            2
    ## 2                       88.000                            0
    ## 3                      493.000                           54
    ## 4                        0.000                            0
    ## 5                        0.000                            0
    ## 6                      179.105                           20
    ##   g.s043.pu.blk_male_hypothalamus_lay g.s043.pu.blk_male_pituitary_lay
    ## 1                            190.0000                         243.1996
    ## 2                              1.0000                         117.0000
    ## 3                            489.0000                         694.0000
    ## 4                              0.0000                           1.0000
    ## 5                              1.0000                           2.0000
    ## 6                             41.1324                         163.3460
    ##   g.s075.pk.pu_male_gonad_m.hatch g.s075.pk.pu_male_hypothalamus_m.hatch
    ## 1                          6.0000                               139.0000
    ## 2                          3.0000                                 0.0000
    ## 3                        712.0000                               490.0000
    ## 4                          1.0000                                 0.0000
    ## 5                         14.5647                                 0.0000
    ## 6                        331.9850                                48.2458
    ##   g.s075.pk.pu_male_pituitary_m.hatch g.s076.pk.r_female_gonad_m.hatch
    ## 1                             398.000                               77
    ## 2                             223.000                              131
    ## 3                             604.000                              549
    ## 4                               0.000                               16
    ## 5                               1.000                                9
    ## 6                             202.075                              151
    ##   g.s076.pk.r_female_hypothalamus_m.hatch
    ## 1                                140.0000
    ## 2                                  0.0000
    ## 3                                444.0000
    ## 4                                  0.0000
    ## 5                                  2.0000
    ## 6                                 31.2192
    ##   g.s076.pk.r_female_pituitary_m.hatch g.s078.blk.o_female_gonad_lay
    ## 1                                  235                     248.00000
    ## 2                                   68                     508.00000
    ## 3                                  272                    1629.00000
    ## 4                                    0                     376.00000
    ## 5                                    2                       2.00738
    ## 6                                  109                     176.00000
    ##   g.s078.blk.o_female_hypothalamus_lay g.s078.blk.o_female_pituitary_lay
    ## 1                             148.0000                          330.2731
    ## 2                               4.0000                           42.0000
    ## 3                             479.0000                          631.0000
    ## 4                               0.0000                            0.0000
    ## 5                               1.0000                            3.9136
    ## 6                              24.1861                          173.4140
    ##   g.s111.r.blk_male_gonad_m.inc.d8 g.s111.r.blk_male_hypothalamus_m.inc.d8
    ## 1                            8.000                                420.0000
    ## 2                            6.000                                 14.0000
    ## 3                          772.000                                987.0000
    ## 4                            6.000                                  0.0000
    ## 5                            3.000                                  2.0000
    ## 6                          236.814                                 87.2718
    ##   g.s111.r.blk_male_pituitary_m.inc.d8 g.s179.o.pk_male_gonad_m.inc.d8
    ## 1                             180.0000                              17
    ## 2                              34.0000                               2
    ## 3                             463.0000                            1141
    ## 4                               0.0000                               2
    ## 5                               2.0000                              10
    ## 6                              76.7042                             404
    ##   g.s179.o.pk_male_hypothalamus_m.inc.d8
    ## 1                             261.000000
    ## 2                               7.000000
    ## 3                             533.000000
    ## 4                               0.000000
    ## 5                               0.377751
    ## 6                              57.145900
    ##   g.s179.o.pk_male_pituitary_m.inc.d8 g.s351.pk.w_male_gonad_extend
    ## 1                             370.000                            45
    ## 2                              15.000                            11
    ## 3                            1158.000                          1842
    ## 4                               0.000                             0
    ## 5                               2.000                             8
    ## 6                             207.395                           581
    ##   g.s351.pk.w_male_hypothalamus_extend g.s351.pk.w_male_pituitary_extend
    ## 1                             154.0000                               172
    ## 2                               0.0000                               103
    ## 3                             283.0000                               408
    ## 4                               0.0000                                 0
    ## 5                               1.0000                                 2
    ## 6                              24.7213                               135
    ##   g.x.ATLAS_female_gonad_control g.y.blk.s006_female_gonad_m.inc.d17
    ## 1                            100                            60.00000
    ## 2                            127                           189.00000
    ## 3                            650                           958.00000
    ## 4                            968                             0.00000
    ## 5                              3                             4.00688
    ## 6                             48                           218.00000
    ##   g.y.blk.s006_female_hypothalamus_m.inc.d17
    ## 1                                    66.0000
    ## 2                                     0.0000
    ## 3                                   247.0000
    ## 4                                     0.0000
    ## 5                                     2.0000
    ## 6                                    22.3853
    ##   g.y.blk.s006_female_pituitary_m.inc.d17 g.y.o.s_male_gonad_prolong
    ## 1                                 226.000                     97.000
    ## 2                                 103.000                      2.000
    ## 3                                 584.000                   1068.000
    ## 4                                   0.000                      7.000
    ## 5                                   6.000                      9.000
    ## 6                                 255.062                    395.959
    ##   g.y.o.s_male_hypothalamus_prolong g.y.o.s_male_pituitary_prolong
    ## 1                           342.000                        299.000
    ## 2                             3.000                         19.000
    ## 3                           616.000                        549.000
    ## 4                             0.000                          0.000
    ## 5                             0.000                          3.000
    ## 6                            55.274                        172.128
    ##   g104.w82.x_male_gonad_bldg g104.w82.x_male_hypothalamus_bldg
    ## 1                      7.000                           64.0000
    ## 2                      0.000                            0.0000
    ## 3                    341.000                          199.0000
    ## 4                      0.000                            0.0000
    ## 5                      0.000                            1.0000
    ## 6                    112.909                           10.1845
    ##   g104.w82.x_male_pituitary_bldg g114.w83.x_male_gonad_hatch.NYNO
    ## 1                        213.000                          25.0883
    ## 2                         42.000                           0.0000
    ## 3                        717.000                        1026.0000
    ## 4                          1.000                           0.0000
    ## 5                          0.000                           2.0000
    ## 6                        153.026                         390.9390
    ##   g114.w83.x_male_hypothalamus_hatch g114.w83.x_male_pituitary_hatch.NYNO
    ## 1                           266.0000                             200.0000
    ## 2                             3.0000                              25.0000
    ## 3                           584.0000                             429.0000
    ## 4                             0.0000                               0.0000
    ## 5                             0.0000                               0.0000
    ## 6                            60.6454                              98.7944
    ##   g130.y81.x_male_gonad_inc.d17 g130.y81.x_male_hypothalamus_inc.d17
    ## 1                         1.000                             179.0000
    ## 2                         1.000                               2.0000
    ## 3                       576.000                             525.0000
    ## 4                         2.000                               0.0000
    ## 5                         3.000                               0.0000
    ## 6                       225.923                              28.1043
    ##   g130.y81.x_male_pituitary_inc.d17 g137.r24.w5_male_gonad_m.inc.d8
    ## 1                           246.000                          29.000
    ## 2                            29.000                           4.000
    ## 3                           641.000                        1627.000
    ## 4                             0.000                           2.000
    ## 5                             1.000                           6.000
    ## 6                           206.651                         398.817
    ##   g137.r24.w5_male_hypothalamus_m.inc.d8
    ## 1                                 229.00
    ## 2                                   0.00
    ## 3                                 552.00
    ## 4                                   0.00
    ## 5                                   0.00
    ## 6                                  56.52
    ##   g137.r24.w5_male_pituitary_m.inc.d8 g141.blu27.x_female_gonad_bldg
    ## 1                             447.348                             54
    ## 2                              74.000                             13
    ## 3                             846.000                            212
    ## 4                               0.000                            140
    ## 5                               2.000                              0
    ## 6                             297.533                             61
    ##   g141.blu27.x_female_hypothalamus_bldg g141.blu27.x_female_pituitary_bldg
    ## 1                               445.000                            203.000
    ## 2                                 2.000                             41.000
    ## 3                               631.000                            807.000
    ## 4                                 0.000                              0.000
    ## 5                                 0.000                              1.000
    ## 6                                35.303                            227.845
    ##   g142.r40.x_female_gonad_inc.d17 g142.r40.x_female_hypothalamus_inc.d17
    ## 1                              70                               135.0000
    ## 2                              74                                 0.0000
    ## 3                             326                               611.0000
    ## 4                              11                                 0.0000
    ## 5                               2                                 0.0000
    ## 6                              37                                49.4282
    ##   g142.r40.x_female_pituitary_inc.d17 g143.blu32.x_male_gonad_inc.d17
    ## 1                                  46                           5.000
    ## 2                                   5                           2.000
    ## 3                                 130                         451.000
    ## 4                                   0                           0.000
    ## 5                                   0                           1.000
    ## 6                                  41                         211.936
    ##   g143.blu32.x_male_hypothalamus_inc.d17
    ## 1                               190.0000
    ## 2                                 2.0000
    ## 3                               424.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                14.1323
    ##   g143.blu32.x_male_pituitary_inc.d17 g144.r54.x_female_gonad_m.inc.d3
    ## 1                            122.0000                          106.000
    ## 2                             14.0000                          223.000
    ## 3                            278.0000                          729.000
    ## 4                              0.0000                            6.000
    ## 5                              1.0000                            2.000
    ## 6                             29.6071                           78.943
    ##   g144.r54.x_female_hypothalamus_m.inc.d3
    ## 1                                118.0000
    ## 2                                  1.0000
    ## 3                                255.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 14.0372
    ##   g144.r54.x_female_pituitary_m.inc.d3 g146.blu51_male_gonad_inc.d3
    ## 1                              284.000                            3
    ## 2                               30.000                            0
    ## 3                              531.000                          562
    ## 4                                0.000                            0
    ## 5                                3.000                            1
    ## 6                              144.855                          113
    ##   g146.blu51_male_hypothalamus_inc.d3.NYNO
    ## 1                                 223.0000
    ## 2                                   6.0000
    ## 3                                 778.0000
    ## 4                                   0.0000
    ## 5                                   0.0000
    ## 6                                  25.2943
    ##   g146.blu51_male_pituitary_inc.d3 g17.w108.x_female_gonad_extend
    ## 1                              104                            516
    ## 2                               25                             53
    ## 3                              296                            611
    ## 4                                0                             11
    ## 5                                0                              4
    ## 6                               60                            108
    ##   g17.w108.x_female_hypothalamus_extend g17.w108.x_female_pituitary_extend
    ## 1                              418.0000                                339
    ## 2                                3.0000                                 11
    ## 3                              824.0000                                784
    ## 4                                0.0000                                  0
    ## 5                                4.0000                                  4
    ## 6                               62.6417                                158
    ##   g20.w106.x_male_gonad_inc.d3 g20.w106.x_male_hypothalamus_inc.d3
    ## 1                           10                                  97
    ## 2                            1                                   1
    ## 3                          561                                 137
    ## 4                            2                                   0
    ## 5                            3                                   0
    ## 6                          197                                  17
    ##   g20.w106.x_male_pituitary_inc.d3 g22.blu118_female_gonad_extend
    ## 1                         132.0000                      200.00000
    ## 2                           4.0000                      102.00000
    ## 3                         332.0000                      630.00000
    ## 4                           1.0000                       25.00000
    ## 5                           0.0000                        5.94182
    ## 6                          56.2266                       70.30920
    ##   g22.blu118_female_hypothalamus_extend g22.blu118_female_pituitary_extend
    ## 1                              255.5420                            224.000
    ## 2                                2.0000                             15.000
    ## 3                              533.0000                            401.000
    ## 4                                0.0000                              5.000
    ## 5                                1.0000                              1.000
    ## 6                               48.2148                            107.239
    ##   g3.g119.w20_male_gonad_extend g3.g119.w20_male_hypothalamus_extend
    ## 1                        32.000                                  126
    ## 2                         5.000                                    1
    ## 3                      1740.000                                  337
    ## 4                         2.000                                    0
    ## 5                         7.000                                    0
    ## 6                       592.975                                   30
    ##   g3.g119.w20_male_pituitary_extend g32.blu79_male_gonad_m.inc.d17
    ## 1                           202.000                       33.00000
    ## 2                            57.000                        6.00000
    ## 3                           436.000                      543.00000
    ## 4                             0.000                        0.00000
    ## 5                             2.000                       12.37906
    ## 6                           140.306                      400.00000
    ##   g32.blu79_male_hypothalamus_m.inc.d17 g32.blu79_male_pituitary_m.inc.d17
    ## 1                             202.00000                           344.0043
    ## 2                               4.00000                            28.0000
    ## 3                             736.00000                           605.0000
    ## 4                               0.00000                             0.0000
    ## 5                               2.56112                             3.0000
    ## 6                              39.04830                           148.1250
    ##   g34.x_male_gonad_m.hatch.NYNO g34.x_male_hypothalamus_m.hatch
    ## 1                            14                        283.0000
    ## 2                             3                          1.0000
    ## 3                           737                        536.0000
    ## 4                             1                          0.0000
    ## 5                             9                          0.0000
    ## 6                           302                         64.1078
    ##   g34.x_male_pituitary_m.hatch g38.x_male_gonad_inc.prolong
    ## 1                      275.000                        6.000
    ## 2                       64.000                        1.000
    ## 3                      392.000                      502.000
    ## 4                        0.000                        2.000
    ## 5                        2.000                        3.000
    ## 6                      142.112                      186.925
    ##   g38.x_male_hypothalamus_inc.prolong g38.x_male_pituitary_inc.prolong
    ## 1                             361.000                               51
    ## 2                               4.000                               13
    ## 3                             845.000                              225
    ## 4                               0.000                                0
    ## 5                               0.000                                0
    ## 6                              48.479                               40
    ##   g52.blu58_male_gonad_bldg g52.blu58_male_hypothalamus_bldg
    ## 1                         9                         262.0000
    ## 2                         3                           0.0000
    ## 3                       470                         370.0000
    ## 4                         0                           0.0000
    ## 5                         1                           0.0000
    ## 6                       144                          27.2265
    ##   g52.blu58_male_pituitary_bldg g53.y84_male_gonad_hatch
    ## 1                       142.000                    8.000
    ## 2                        58.000                    2.000
    ## 3                       419.000                  553.000
    ## 4                         0.000                    0.000
    ## 5                         3.000                    0.000
    ## 6                       119.003                  172.985
    ##   g53.y84_male_hypothalamus_hatch g53.y84_male_pituitary_hatch
    ## 1                        72.00000                           52
    ## 2                         4.00000                           15
    ## 3                       231.00000                          178
    ## 4                         0.00000                            0
    ## 5                         0.00000                            0
    ## 6                         8.06621                           47
    ##   g6.w197.x_female_gonad_inc.d3 g6.w197.x_female_hypothalamus_inc.d3
    ## 1                            40                             359.0000
    ## 2                            19                               4.0000
    ## 3                           208                            1078.0000
    ## 4                             8                               0.0000
    ## 5                             1                               0.0000
    ## 6                            31                              89.3317
    ##   g6.w197.x_female_pituitary_inc.d3 g63.blu65_female_gonad_m.inc.d17
    ## 1                           229.000                              153
    ## 2                            82.000                              174
    ## 3                           766.000                              836
    ## 4                             0.000                               37
    ## 5                             3.000                               15
    ## 6                           180.453                              177
    ##   g63.blu65_female_hypothalamus_m.inc.d17
    ## 1                                263.0000
    ## 2                                  3.0000
    ## 3                                551.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 53.5458
    ##   g63.blu65_female_pituitary_m.inc.d17 g73.x_female_gonad_m.inc.d9
    ## 1                              253.000                          32
    ## 2                               27.000                          85
    ## 3                              660.000                         181
    ## 4                                1.000                          30
    ## 5                                0.000                           2
    ## 6                              248.763                          24
    ##   g73.x_female_hypothalamus_m.inc.d9 g73.x_female_pituitary_m.inc.d9
    ## 1                           348.0000                        122.0000
    ## 2                             3.0000                          4.0000
    ## 3                           725.0000                        286.0000
    ## 4                             0.0000                          0.0000
    ## 5                             2.0000                          0.0000
    ## 6                            63.0935                         79.6828
    ##   g75.x_female_gonad_inc.d9 g75.x_female_hypothalamus_inc.d9
    ## 1                        37                         603.0000
    ## 2                        59                           7.0000
    ## 3                       309                        1367.0000
    ## 4                        14                           0.0000
    ## 5                         0                           2.0000
    ## 6                        22                          90.4118
    ##   g75.x_female_pituitary_inc.d9 g8.y197_male_gonad_extend
    ## 1                       222.000                    16.000
    ## 2                        44.000                     7.000
    ## 3                       569.000                  1360.000
    ## 4                         0.000                     1.000
    ## 5                         1.000                     9.000
    ## 6                       134.608                   516.808
    ##   g8.y197_male_hypothalamus_extend g8.y197_male_pituitary_extend
    ## 1                         218.0000                     175.00000
    ## 2                           0.0000                      71.00000
    ## 3                         426.0000                     377.00000
    ## 4                           0.0000                       0.00000
    ## 5                           1.0000                       1.64771
    ## 6                          34.6207                     155.89600
    ##   l.s.o.blk_male_gonad_extend l.s.o.blk_male_hypothalamus_extend
    ## 1                      16.000                           361.0000
    ## 2                       6.000                             6.0000
    ## 3                     725.000                          1008.0000
    ## 4                       5.000                             0.0000
    ## 5                       5.000                             1.0000
    ## 6                     336.986                            69.4175
    ##   l.s.o.blk_male_pituitary_extend l.s.w.d_female_gonad_m.hatch
    ## 1                         512.000                      182.000
    ## 2                         119.000                      229.000
    ## 3                         650.000                      701.000
    ## 4                           0.000                       28.000
    ## 5                           5.000                        1.000
    ## 6                         266.052                      125.902
    ##   l.s.w.d_female_hypothalamus_m.hatch l.s.w.d_female_pituitary_m.hatch
    ## 1                            178.0000                          484.000
    ## 2                              2.0000                           45.000
    ## 3                            397.0000                          636.000
    ## 4                              0.0000                            2.000
    ## 5                              1.0000                            2.000
    ## 6                             23.0574                          279.357
    ##   l.s024.y.g_male_gonad_m.inc.d17 l.s024.y.g_male_hypothalamus_m.inc.d17
    ## 1                          54.000                                85.0000
    ## 2                           0.000                                 0.0000
    ## 3                         802.000                               346.0000
    ## 4                           2.000                                 0.0000
    ## 5                           3.000                                 1.0000
    ## 6                         591.981                                19.2081
    ##   l.s024.y.g_male_pituitary_m.inc.d17 l.s052.pk.r_female_gonad_prolong
    ## 1                             338.000                               91
    ## 2                              24.000                               54
    ## 3                             467.000                              149
    ## 4                               0.000                               23
    ## 5                               0.000                                0
    ## 6                             146.359                               78
    ##   l.s052.pk.r_female_hypothalamus_prolong
    ## 1                                116.0000
    ## 2                                  0.0000
    ## 3                                315.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 33.3344
    ##   l.s052.pk.r_female_pituitary_prolong l.s080.blk.r_male_gonad_prolong
    ## 1                            443.00000                          22.000
    ## 2                            110.00000                           7.000
    ## 3                            755.00000                        1076.000
    ## 4                              0.00000                           1.000
    ## 5                              1.65596                           2.000
    ## 6                            506.38900                         450.953
    ##   l.s080.blk.r_male_hypothalamus_prolong
    ## 1                                    123
    ## 2                                      1
    ## 3                                    278
    ## 4                                      0
    ## 5                                      0
    ## 6                                     16
    ##   l.s080.blk.r_male_pituitary_prolong l.s120.y.blk_female_gonad_bldg
    ## 1                             348.000                        142.000
    ## 2                              42.000                        392.000
    ## 3                             615.000                        895.000
    ## 4                              19.000                         18.000
    ## 5                               5.000                          1.000
    ## 6                             157.277                        116.804
    ##   l.s120.y.blk_female_hypothalamus_bldg l.s120.y.blk_female_pituitary_bldg
    ## 1                                   264                          357.02200
    ## 2                                     4                           39.00000
    ## 3                                   325                          633.00000
    ## 4                                     0                            0.00000
    ## 5                                     0                            4.07053
    ## 6                                    26                          273.84700
    ##   l.s166.o.w_female_gonad_prolong l.s166.o.w_female_hypothalamus_prolong
    ## 1                        160.0000                               135.0000
    ## 2                        132.0000                                 2.0000
    ## 3                        656.0000                               303.0000
    ## 4                         35.0000                                 0.0000
    ## 5                          8.0473                                 1.0000
    ## 6                        139.0000                                21.7751
    ##   l.s166.o.w_female_pituitary_prolong l.s280.g.blk_male_gonad_m.inc.d9
    ## 1                             295.000                               18
    ## 2                              15.000                                4
    ## 3                             402.000                             1466
    ## 4                               1.000                                4
    ## 5                               5.000                                8
    ## 6                             174.643                              440
    ##   l.s280.g.blk_male_hypothalamus_m.inc.d9
    ## 1                                389.0000
    ## 2                                  7.0000
    ## 3                                660.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 57.2579
    ##   l.s280.g.blk_male_pituitary_m.inc.d9 o.d.s009_male_gonad_m.inc.d9
    ## 1                              323.000                      12.0895
    ## 2                               40.000                       2.0000
    ## 3                              401.000                     322.0000
    ## 4                                0.000                       0.0000
    ## 5                                3.000                       4.0000
    ## 6                              137.913                     115.0000
    ##   o.d.s009_male_hypothalamus_m.inc.d9.NYNO
    ## 1                                 147.0000
    ## 2                                   1.0000
    ## 3                                 845.0000
    ## 4                                   0.0000
    ## 5                                   3.0000
    ## 6                                  54.4877
    ##   o.d.s009_male_pituitary_m.inc.d9 o.s.w.r_male_gonad_lay
    ## 1                          94.0000                     10
    ## 2                          30.0000                      4
    ## 3                         248.0000                    651
    ## 4                           0.0000                      2
    ## 5                           0.0000                      2
    ## 6                          65.8021                    267
    ##   o.s.w.r_male_hypothalamus_lay o.s.w.r_male_pituitary_lay
    ## 1                           257                    446.212
    ## 2                             0                    100.000
    ## 3                           603                    611.000
    ## 4                             0                      0.000
    ## 5                             1                      2.000
    ## 6                            32                    178.323
    ##   o.s010.r.blk_male_gonad_m.inc.d3 o.s010.r.blk_male_hypothalamus_m.inc.d3
    ## 1                         17.00000                                 118.000
    ## 2                          3.00000                                   4.000
    ## 3                        789.00000                                 299.000
    ## 4                          2.00000                                   0.000
    ## 5                          3.80319                                   0.000
    ## 6                        452.81300                                  16.341
    ##   o.s010.r.blk_male_pituitary_m.inc.d3 o.s084.w.blk_female_gonad_m.inc.d17
    ## 1                                  234                                 283
    ## 2                                  177                                 448
    ## 3                                  490                                1087
    ## 4                                    0                                  45
    ## 5                                    3                                   7
    ## 6                                  135                                 260
    ##   o.s084.w.blk_female_hypothalamus_m.inc.d17
    ## 1                                   133.0000
    ## 2                                     0.0000
    ## 3                                   312.0000
    ## 4                                     0.0000
    ## 5                                     0.0000
    ## 6                                    25.0745
    ##   o.s084.w.blk_female_pituitary_m.inc.d17 o114.blu9_male_gonad_inc.prolong
    ## 1                                 441.000                               13
    ## 2                                  82.000                                3
    ## 3                                 521.000                              591
    ## 4                                   0.000                                1
    ## 5                                   0.000                                0
    ## 6                                 198.724                              219
    ##   o114.blu9_male_hypothalamus_inc.prolong
    ## 1                                252.0000
    ## 2                                  1.0000
    ## 3                                732.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 44.3434
    ##   o114.blu9_male_pituitary_inc.prolong o152.o120.w42_male_gonad_n5
    ## 1                              98.0000                    80.00000
    ## 2                              84.0000                     4.00000
    ## 3                             246.0000                  1134.00000
    ## 4                               0.0000                     2.00000
    ## 5                               0.0000                     3.79367
    ## 6                              93.7873                   692.93800
    ##   o152.o120.w42_male_hypothalamus_n5 o152.o120.w42_male_pituitary_n5
    ## 1                           241.0000                          323.00
    ## 2                             3.0000                          104.00
    ## 3                           508.0000                          750.00
    ## 4                             0.0000                            0.00
    ## 5                             0.0000                            0.00
    ## 6                            54.4265                          280.78
    ##   o156.w80.x_female_gonad_inc.d3 o156.w80.x_female_hypothalamus_inc.d3
    ## 1                             30                              372.0000
    ## 2                             53                                5.0000
    ## 3                            243                              567.0000
    ## 4                              7                                0.0000
    ## 5                              0                                1.0000
    ## 6                             18                               60.3712
    ##   o156.w80.x_female_pituitary_inc.d3 o165.w122.x_female_gonad_inc.d3.NYNO
    ## 1                            275.000                                  121
    ## 2                            159.000                                  201
    ## 3                            491.000                                  780
    ## 4                              0.000                                  186
    ## 5                              1.000                                    0
    ## 6                            183.222                                   57
    ##   o165.w122.x_female_hypothalamus_inc.d3
    ## 1                               333.0000
    ## 2                                 1.0000
    ## 3                               804.0000
    ## 4                                 0.0000
    ## 5                                 2.0000
    ## 6                                71.2884
    ##   o165.w122.x_female_pituitary_inc.d3.NYNO
    ## 1                                  348.000
    ## 2                                  113.000
    ## 3                                  599.000
    ## 4                                    0.000
    ## 5                                    2.000
    ## 6                                  120.699
    ##   o169.r28.x_female_gonad_m.inc.d17
    ## 1                          53.00000
    ## 2                         135.00000
    ## 3                         979.00000
    ## 4                          14.00000
    ## 5                           5.01482
    ## 6                         153.00000
    ##   o169.r28.x_female_hypothalamus_m.inc.d17
    ## 1                                  93.0000
    ## 2                                   0.0000
    ## 3                                 279.0000
    ## 4                                   0.0000
    ## 5                                   0.0000
    ## 6                                  18.1848
    ##   o169.r28.x_female_pituitary_m.inc.d17 o172.w115.x_female_gonad_hatch
    ## 1                             314.04870                             52
    ## 2                              48.00000                             21
    ## 3                             716.00000                            285
    ## 4                               0.00000                             64
    ## 5                               4.03228                              0
    ## 6                             259.22500                             18
    ##   o172.w115.x_female_hypothalamus_hatch
    ## 1                              587.0000
    ## 2                                4.0000
    ## 3                              735.0000
    ## 4                                1.0000
    ## 5                                2.0000
    ## 6                               56.5556
    ##   o172.w115.x_female_pituitary_hatch.NYNO o173.w179.x_female_gonad_inc.d3
    ## 1                                139.0000                              39
    ## 2                                 18.0000                             103
    ## 3                                339.0000                             339
    ## 4                                  0.0000                               1
    ## 5                                  0.0000                               1
    ## 6                                 99.8199                              14
    ##   o173.w179.x_female_hypothalamus_inc.d3
    ## 1                               166.0000
    ## 2                                 1.0000
    ## 3                               517.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                28.0671
    ##   o173.w179.x_female_pituitary_inc.d3 o3.x_male_gonad_m.n2
    ## 1                             99.0000               15.000
    ## 2                             19.0000                1.000
    ## 3                            189.0000              410.000
    ## 4                              0.0000                1.000
    ## 5                              1.0000                1.000
    ## 6                             48.3061              170.904
    ##   o3.x_male_hypothalamus_m.n2 o3.x_male_pituitary_m.n2
    ## 1                    217.0000                  83.0000
    ## 2                      7.0000                  27.0000
    ## 3                    489.0000                 204.0000
    ## 4                      0.0000                   0.0000
    ## 5                      1.0000                   0.0000
    ## 6                     38.0962                  64.9791
    ##   o35.r51.x_female_gonad_inc.d17 o35.r51.x_female_hypothalamus_inc.d17
    ## 1                            158                              435.0000
    ## 2                            170                                2.0000
    ## 3                            771                              743.0000
    ## 4                            260                                0.0000
    ## 5                              3                                1.0000
    ## 6                            123                               39.5935
    ##   o35.r51.x_female_pituitary_inc.d17 o36.r62.x_female_gonad_m.inc.d9
    ## 1                                109                              49
    ## 2                                 30                              84
    ## 3                                286                             285
    ## 4                                  0                               2
    ## 5                                  1                               2
    ## 6                                 88                              23
    ##   o36.r62.x_female_hypothalamus_m.inc.d9
    ## 1                               278.0000
    ## 2                                 6.0000
    ## 3                               549.0000
    ## 4                                 0.0000
    ## 5                                 2.0000
    ## 6                                26.0506
    ##   o36.r62.x_female_pituitary_m.inc.d9 o38.blu29.x_female_gonad_bldg
    ## 1                                  93                            61
    ## 2                                  10                           173
    ## 3                                 241                           497
    ## 4                                   0                           135
    ## 5                                   0                             3
    ## 6                                  36                            82
    ##   o38.blu29.x_female_hypothalamus_bldg o38.blu29.x_female_pituitary_bldg
    ## 1                             205.0000                          143.0000
    ## 2                               2.0000                           19.0000
    ## 3                             582.0000                          263.0000
    ## 4                               0.0000                            0.0000
    ## 5                               0.0000                            3.0000
    ## 6                              45.2665                           85.1016
    ##   o39.y77.x_male_gonad_hatch o39.y77.x_male_hypothalamus_hatch
    ## 1                   14.00000                           463.791
    ## 2                    0.00000                             5.000
    ## 3                  623.00000                           631.000
    ## 4                    0.00000                             0.000
    ## 5                    6.75478                             0.000
    ## 6                  153.76000                            39.713
    ##   o39.y77.x_male_pituitary_hatch o44.blu26.x_male_gonad_hatch
    ## 1                             93                            7
    ## 2                             61                            1
    ## 3                            226                          365
    ## 4                              0                            1
    ## 5                              0                            1
    ## 6                             92                          126
    ##   o44.blu26.x_male_hypothalamus_hatch o44.blu26.x_male_pituitary_hatch
    ## 1                            334.0000                          311.000
    ## 2                              5.0000                           25.000
    ## 3                            811.0000                          417.000
    ## 4                              0.0000                            0.000
    ## 5                              1.0000                            0.000
    ## 6                             66.3776                          170.714
    ##   o45.g128.x_female_gonad_m.inc.d9.NYNO
    ## 1                                    51
    ## 2                                   133
    ## 3                                   599
    ## 4                                     5
    ## 5                                     9
    ## 6                                    85
    ##   o45.g128.x_female_hypothalamus_m.inc.d9
    ## 1                                229.0000
    ## 2                                  2.0000
    ## 3                                651.0000
    ## 4                                  0.0000
    ## 5                                  3.0000
    ## 6                                 35.2079
    ##   o45.g128.x_female_pituitary_m.inc.d9 o48.r197.x_male_gonad_inc.d3
    ## 1                              312.000                            7
    ## 2                              113.000                            2
    ## 3                              965.000                          281
    ## 4                                3.000                            3
    ## 5                                4.000                            3
    ## 6                              253.029                          110
    ##   o48.r197.x_male_hypothalamus_inc.d3 o48.r197.x_male_pituitary_inc.d3
    ## 1                            238.0000                            189.0
    ## 2                              8.0000                            249.0
    ## 3                            565.0000                            767.0
    ## 4                              0.0000                              0.0
    ## 5                              0.0000                              2.0
    ## 6                             48.0802                            176.2
    ##   o49.x_male_gonad_inc.d9 o49.x_male_hypothalamus_inc.d9
    ## 1                  10.000                            142
    ## 2                   3.000                              3
    ## 3                 300.000                            231
    ## 4                   0.000                              0
    ## 5                   3.000                              0
    ## 6                 197.982                             17
    ##   o49.x_male_pituitary_inc.d9 o52.blu53_female_gonad_inc.d17
    ## 1                     88.1613                            113
    ## 2                     52.0000                            225
    ## 3                    218.0000                            703
    ## 4                      0.0000                            142
    ## 5                      2.0000                              4
    ## 6                     65.2782                             60
    ##   o52.blu53_female_hypothalamus_inc.d17 o52.blu53_female_pituitary_inc.d17
    ## 1                              323.0000                           174.0000
    ## 2                                1.0000                            12.0000
    ## 3                              644.0000                           356.0000
    ## 4                                0.0000                             0.0000
    ## 5                                5.0000                             3.0000
    ## 6                               37.8189                            69.9528
    ##   o57.g59_male_gonad_inc.d9 o57.g59_male_hypothalamus_inc.d9
    ## 1                     6.000                         126.0000
    ## 2                     1.000                           2.0000
    ## 3                   565.000                         205.0000
    ## 4                     0.000                           0.0000
    ## 5                     2.000                           1.0000
    ## 6                   218.655                          14.4038
    ##   o57.g59_male_pituitary_inc.d9 o59.blu64_male_gonad_m.inc.d17
    ## 1                           104                       18.00000
    ## 2                            28                        4.00000
    ## 3                           259                      960.00000
    ## 4                             0                        1.00000
    ## 5                             1                        8.66485
    ## 6                            70                      515.00000
    ##   o59.blu64_male_hypothalamus_m.inc.d17 o59.blu64_male_pituitary_m.inc.d17
    ## 1                              339.4690                           267.2821
    ## 2                                0.0000                           110.0000
    ## 3                              714.0000                           557.0000
    ## 4                                0.0000                             0.0000
    ## 5                                2.0000                             0.0000
    ## 6                               48.3326                           206.2510
    ##   o73.x_female_gonad_inc.d9 o73.x_female_hypothalamus_inc.d9
    ## 1                        35                              129
    ## 2                        73                                1
    ## 3                       366                              359
    ## 4                        63                                0
    ## 5                         1                                1
    ## 6                        47                               24
    ##   o73.x_female_pituitary_inc.d9 p.g.blk.s040_female_gonad_m.inc.d8
    ## 1                           148                                176
    ## 2                            15                                182
    ## 3                           171                                611
    ## 4                             0                                110
    ## 5                             0                                  2
    ## 6                            46                                 42
    ##   p.g.blk.s040_female_hypothalamus_m.inc.d8
    ## 1                                   221.000
    ## 2                                     5.000
    ## 3                                   591.000
    ## 4                                     0.000
    ## 5                                     1.000
    ## 6                                    36.132
    ##   p.g.blk.s040_female_pituitary_m.inc.d8 pk.s.d.g_male_gonad_prolong
    ## 1                                292.544                         108
    ## 2                                 26.000                         215
    ## 3                                474.000                         831
    ## 4                                  0.000                           8
    ## 5                                  0.000                          13
    ## 6                                130.283                         234
    ##   pk.s.d.g_male_hypothalamus_prolong pk.s.d.g_male_pituitary_prolong
    ## 1                            97.0000                         224.000
    ## 2                             1.0000                         108.000
    ## 3                           364.0000                         473.000
    ## 4                             0.0000                           0.000
    ## 5                             1.0000                           1.000
    ## 6                            34.1803                         172.298
    ##   pk.s011.o.y_female_gonad_m.inc.d17
    ## 1                           283.1745
    ## 2                           597.0000
    ## 3                          1838.0000
    ## 4                            30.0000
    ## 5                             4.0000
    ## 6                           188.0000
    ##   pk.s011.o.y_female_hypothalamus_m.inc.d17
    ## 1                                       133
    ## 2                                         5
    ## 3                                       393
    ## 4                                         0
    ## 5                                         0
    ## 6                                        15
    ##   pk.s011.o.y_female_pituitary_m.inc.d17 pk.s054.d.g_female_gonad_m.hatch
    ## 1                                232.000                               71
    ## 2                                 51.000                              192
    ## 3                                416.000                              463
    ## 4                                  0.000                               14
    ## 5                                  8.000                               10
    ## 6                                136.775                              102
    ##   pk.s054.d.g_female_hypothalamus_m.hatch
    ## 1                                217.0000
    ## 2                                  1.0000
    ## 3                                603.0000
    ## 4                                  0.0000
    ## 5                                  1.0000
    ## 6                                 39.1933
    ##   pk.s054.d.g_female_pituitary_m.hatch pk.s055.d.l_female_gonad_m.inc.d8
    ## 1                              559.000                          64.00000
    ## 2                               57.000                         126.00000
    ## 3                              889.000                         925.00000
    ## 4                                0.000                         423.00000
    ## 5                                5.000                           4.02273
    ## 6                              290.393                          99.00000
    ##   pk.s055.d.l_female_hypothalamus_m.inc.d8
    ## 1                                 239.0000
    ## 2                                   1.0000
    ## 3                                 547.0000
    ## 4                                   0.0000
    ## 5                                   1.0000
    ## 6                                  56.5185
    ##   pk.s055.d.l_female_pituitary_m.inc.d8 pk.s238.blk.w_male_gonad_lay
    ## 1                               273.000                     69.00000
    ## 2                                38.000                     14.00000
    ## 3                               768.000                   2576.00000
    ## 4                                 0.000                      6.00000
    ## 5                                 2.000                     16.43527
    ## 6                               347.064                    890.20200
    ##   pk.s238.blk.w_male_hypothalamus_lay pk.s238.blk.w_male_pituitary_lay
    ## 1                            206.0000                        305.00000
    ## 2                              3.0000                         55.00000
    ## 3                            720.0000                        540.00000
    ## 4                              0.0000                          1.00000
    ## 5                              1.0000                          1.06018
    ## 6                             44.4001                        128.69600
    ##   pk.w.s141.o_male_gonad_lay pk.w.s141.o_male_hypothalamus_lay
    ## 1                          9                          348.0000
    ## 2                          7                            3.0000
    ## 3                       2365                          689.0000
    ## 4                          5                            1.0000
    ## 5                          5                            0.0000
    ## 6                        481                           47.6348
    ##   pk.w.s141.o_male_pituitary_lay pu.blk.s102.y_male_gonad_m.inc.d8
    ## 1                        345.000                            22.000
    ## 2                        106.000                             3.000
    ## 3                        723.000                          1722.000
    ## 4                          0.000                             4.000
    ## 5                          4.000                             9.000
    ## 6                        244.974                           648.811
    ##   pu.blk.s102.y_male_hypothalamus_m.inc.d8
    ## 1                                 178.0000
    ## 2                                   0.0000
    ## 3                                 457.0000
    ## 4                                   0.0000
    ## 5                                   1.0000
    ## 6                                  26.0551
    ##   pu.blk.s102.y_male_pituitary_m.inc.d8 pu.s.o.r_male_hypothalamus_m.hatch
    ## 1                               204.000                           213.0000
    ## 2                                46.000                             3.0000
    ## 3                               503.000                           554.0000
    ## 4                                 1.000                             0.0000
    ## 5                                 3.000                             2.0000
    ## 6                               192.009                            50.2942
    ##   pu.s.o.r_male_pituitary_m.hatch r.r.x.ATLAS.R2XR_female_gonad_control
    ## 1                         359.000                                    16
    ## 2                          25.000                                    22
    ## 3                         502.000                                   259
    ## 4                           0.000                                    17
    ## 5                           0.000                                     0
    ## 6                         211.369                                    23
    ##   r.r.x.ATLAS.R2XR_female_hypothalamus_control
    ## 1                                      48.0000
    ## 2                                       0.0000
    ## 3                                     163.0000
    ## 4                                       0.0000
    ## 5                                       0.0000
    ## 6                                      19.2357
    ##   r.r.x.ATLAS.R2XR_female_pituitary_control
    ## 1                                   88.0000
    ## 2                                   14.0000
    ## 3                                  228.0000
    ## 4                                    0.0000
    ## 5                                    1.0000
    ## 6                                   48.2428
    ##   r.r.x.ATLAS_female_gonad_control r.r.x.ATLAS_female_hypothalamus_control
    ## 1                              134                                 48.0000
    ## 2                              154                                  0.0000
    ## 3                              773                                163.0000
    ## 4                              806                                  0.0000
    ## 5                                0                                  0.0000
    ## 6                               28                                 19.2353
    ##   r.r.x.ATLAS_female_pituitary_control r.s005.pk.blk_male_gonad_lay
    ## 1                                   23                           44
    ## 2                                    5                            2
    ## 3                                   71                          938
    ## 4                                    0                            0
    ## 5                                    0                            8
    ## 6                                   12                          692
    ##   r.s005.pk.blk_male_hypothalamus_lay r.s005.pk.blk_male_pituitary_lay
    ## 1                            288.0000                          146.000
    ## 2                              5.0000                           12.000
    ## 3                            459.0000                          636.000
    ## 4                              0.0000                            0.000
    ## 5                              1.0000                            3.000
    ## 6                             32.2842                          143.958
    ##   r.s035.y.blk_female_gonad_m.inc.d3
    ## 1                          219.00000
    ## 2                          150.00000
    ## 3                          852.00000
    ## 4                           79.00000
    ## 5                            5.02095
    ## 6                          158.00000
    ##   r.s035.y.blk_female_hypothalamus_m.inc.d3
    ## 1                                  321.0000
    ## 2                                    3.0000
    ## 3                                  669.0000
    ## 4                                    0.0000
    ## 5                                    1.0000
    ## 6                                   54.1273
    ##   r.s035.y.blk_female_pituitary_m.inc.d3 r.s056.g.o_female_gonad_bldg
    ## 1                                497.000                           15
    ## 2                                 21.000                           49
    ## 3                                530.000                          265
    ## 4                                  0.000                           90
    ## 5                                  0.000                            0
    ## 6                                176.551                            8
    ##   r.s056.g.o_female_hypothalamus_bldg r.s056.g.o_female_pituitary_bldg
    ## 1                            154.0000                          94.0000
    ## 2                              0.0000                           7.0000
    ## 3                            237.0000                         388.0000
    ## 4                              0.0000                           0.0000
    ## 5                              0.0000                           0.0000
    ## 6                             16.1403                          55.2826
    ##   r.s057.g.pk_male_gonad_extend r.s057.g.pk_male_hypothalamus_extend
    ## 1                            26                             152.0000
    ## 2                             3                               1.0000
    ## 3                           610                             576.0000
    ## 4                             0                               0.0000
    ## 5                             2                               1.0000
    ## 6                           291                              35.1602
    ##   r.s057.g.pk_male_pituitary_extend r.s058.d.l_male_gonad_m.inc.d8
    ## 1                           265.000                         22.000
    ## 2                           122.000                          2.000
    ## 3                           451.000                        579.000
    ## 4                             0.000                          3.000
    ## 5                             3.000                          2.000
    ## 6                           169.844                        388.958
    ##   r.s058.d.l_male_hypothalamus_m.inc.d8 r.s058.d.l_male_pituitary_m.inc.d8
    ## 1                              231.0000                            254.000
    ## 2                                0.0000                             38.000
    ## 3                              639.0000                            556.000
    ## 4                                0.0000                              0.000
    ## 5                                0.0000                              0.000
    ## 6                               24.2426                            205.905
    ##   r.s059.d.o_male_gonad_bldg r.s059.d.o_male_hypothalamus_bldg
    ## 1                      5.000                               364
    ## 2                      1.000                                 5
    ## 3                    471.000                               637
    ## 4                      0.000                                 0
    ## 5                      2.000                                 0
    ## 6                    192.961                                49
    ##   r.s059.d.o_male_pituitary_bldg r.s086.l.blk_male_gonad_extend
    ## 1                       164.0000                         12.000
    ## 2                         9.0000                          1.000
    ## 3                       245.0000                       1486.000
    ## 4                         0.0000                          3.000
    ## 5                         0.0000                          3.000
    ## 6                        70.5407                        489.739
    ##   r.s086.l.blk_male_hypothalamus_extend r.s086.l.blk_male_pituitary_extend
    ## 1                              218.0000                             160.00
    ## 2                                4.0000                              41.00
    ## 3                              680.0000                             344.00
    ## 4                                0.0000                               0.00
    ## 5                                2.0000                               1.00
    ## 6                               61.1766                             207.88
    ##   r.s116.blk.pu_male_gonad_lay r.s116.blk.pu_male_hypothalamus_lay
    ## 1                            9                            229.0000
    ## 2                            2                              0.0000
    ## 3                          805                            558.0000
    ## 4                            0                              1.0000
    ## 5                            4                              0.0000
    ## 6                          304                             40.1807
    ##   r.s116.blk.pu_male_pituitary_lay r.s131.o.d_female_gonad_extend
    ## 1                          441.000                        195.000
    ## 2                           46.000                        305.000
    ## 3                          780.000                        686.000
    ## 4                            0.000                         31.000
    ## 5                            3.000                          2.000
    ## 6                          205.162                        109.928
    ##   r.s131.o.d_female_hypothalamus_extend r.s131.o.d_female_pituitary_extend
    ## 1                              496.0000                            524.000
    ## 2                                3.0000                             22.000
    ## 3                              949.0000                            982.000
    ## 4                                0.0000                              0.000
    ## 5                                6.0000                              6.000
    ## 6                               79.8977                            343.488
    ##   r.s171.l.w_female_gonad_n9 r.s171.l.w_female_hypothalamus_n9
    ## 1                    213.000                          280.0000
    ## 2                    214.000                            0.0000
    ## 3                   1060.000                          709.0000
    ## 4                     11.000                            0.0000
    ## 5                      8.000                            0.0000
    ## 6                    326.992                           55.3649
    ##   r.s171.l.w_female_pituitary_n9 r.s172.l.y_male_gonad_extend
    ## 1                        208.000                       46.000
    ## 2                         51.000                        8.000
    ## 3                        496.000                      620.000
    ## 4                          0.000                        0.000
    ## 5                          1.000                        4.000
    ## 6                        155.265                      246.953
    ##   r.s172.l.y_male_hypothalamus_extend r.s172.l.y_male_pituitary_extend
    ## 1                            346.0000                          167.000
    ## 2                              7.0000                           84.000
    ## 3                            867.0000                          409.000
    ## 4                              0.0000                            0.000
    ## 5                              1.0000                            2.000
    ## 6                             74.0479                          134.466
    ##   r.y.s007.blk_male_gonad_n9 r.y.s007.blk_male_hypothalamus_n9
    ## 1                     14.000                               102
    ## 2                      2.000                                 1
    ## 3                    771.000                               229
    ## 4                      1.000                                 0
    ## 5                      5.000                                 0
    ## 6                    252.869                                 9
    ##   r.y.s007.blk_male_pituitary_n9 r176.blu54_male_gonad_inc.d17
    ## 1                       123.0000                         3.000
    ## 2                        15.0000                         2.000
    ## 3                       218.0000                       886.000
    ## 4                         0.0000                         0.000
    ## 5                         0.0000                         1.000
    ## 6                        90.2129                       236.929
    ##   r176.blu54_male_hypothalamus_inc.d17 r176.blu54_male_pituitary_inc.d17
    ## 1                            118.00000                           156.000
    ## 2                              0.00000                            14.000
    ## 3                            209.00000                           338.000
    ## 4                              0.00000                             0.000
    ## 5                              0.00000                             0.000
    ## 6                              7.13669                           104.282
    ##   r183.o22_female_gonad_hatch r183.o22_female_hypothalamus_hatch
    ## 1                          60                           345.6826
    ## 2                          16                             5.0000
    ## 3                         196                           670.0000
    ## 4                          11                             0.0000
    ## 5                           0                             0.0000
    ## 6                          23                            46.5800
    ##   r183.o22_female_pituitary_hatch r190.o43.x_male_gonad_lay
    ## 1                              76                     3.000
    ## 2                              10                     3.000
    ## 3                             303                   906.000
    ## 4                               0                     2.000
    ## 5                               0                     8.000
    ## 6                              70                   234.871
    ##   r190.o43.x_male_hypothalamus_lay r190.o43.x_male_pituitary_lay
    ## 1                          273.000                      130.0000
    ## 2                            7.000                       51.0000
    ## 3                          696.000                      528.0000
    ## 4                            0.000                        0.0000
    ## 5                            2.000                        2.0000
    ## 6                           38.167                       69.0822
    ##   r194.x_female_gonad_prolong r194.x_female_hypothalamus_prolong
    ## 1                         163                             92.000
    ## 2                         268                              0.000
    ## 3                         989                            283.000
    ## 4                          19                              0.000
    ## 5                          12                              0.000
    ## 6                         105                             22.054
    ##   r194.x_female_pituitary_prolong r195.x_male_gonad_n9
    ## 1                        332.3145                 51.0
    ## 2                         21.0000                  7.0
    ## 3                        541.0000               1209.0
    ## 4                          0.0000                  1.0
    ## 5                          2.0000                 11.0
    ## 6                        268.7010                539.5
    ##   r195.x_male_hypothalamus_n9 r195.x_male_pituitary_n9
    ## 1                    283.1267                  244.000
    ## 2                      6.0000                   69.000
    ## 3                    427.0000                  476.000
    ## 4                      0.0000                    0.000
    ## 5                      0.0000                    4.000
    ## 6                     24.0000                  135.442
    ##   r196.x_male_gonad_m.inc.d9 r196.x_male_hypothalamus_m.inc.d9
    ## 1                          9                          251.0000
    ## 2                          1                            2.0000
    ## 3                        479                          683.0000
    ## 4                          1                            0.0000
    ## 5                          1                            1.0000
    ## 6                        249                           59.3666
    ##   r196.x_male_pituitary_m.inc.d9 r27.w111.blu125_female_gonad_inc.d3
    ## 1                       126.0000                            102.2163
    ## 2                        64.0000                            149.0000
    ## 3                       249.0000                            503.0000
    ## 4                         0.0000                              9.0000
    ## 5                         0.0000                              2.0000
    ## 6                        46.3643                             69.9099
    ##   r27.w111.blu125_female_hypothalamus_inc.d3
    ## 1                                   255.0000
    ## 2                                     5.0000
    ## 3                                   640.0000
    ## 4                                     0.0000
    ## 5                                     1.0000
    ## 6                                    69.4284
    ##   r27.w111.blu125_female_pituitary_inc.d3 r30.w112.r46_female_gonad_inc.d9
    ## 1                                  86.000                               44
    ## 2                                  21.000                               78
    ## 3                                 828.000                              317
    ## 4                                   0.000                                7
    ## 5                                   0.000                                1
    ## 6                                 123.205                               43
    ##   r30.w112.r46_female_hypothalamus_inc.d9
    ## 1                                 258.000
    ## 2                                   9.000
    ## 3                                 863.000
    ## 4                                   0.000
    ## 5                                   0.000
    ## 6                                  81.642
    ##   r30.w112.r46_female_pituitary_inc.d9 r36.w184.x_female_gonad_inc.d9
    ## 1                              370.000                             33
    ## 2                               11.000                             25
    ## 3                              611.000                            186
    ## 4                                0.000                             13
    ## 5                                7.000                              0
    ## 6                              242.469                             16
    ##   r36.w184.x_female_hypothalamus_inc.d9 r36.w184.x_female_pituitary_inc.d9
    ## 1                              335.5032                            64.0000
    ## 2                                3.0000                             8.0000
    ## 3                              634.0000                           145.0000
    ## 4                                0.0000                             0.0000
    ## 5                                2.0000                             0.0000
    ## 6                               50.0884                            29.1238
    ##   r37.w100.x_male_gonad_n9 r37.w100.x_male_hypothalamus_n9
    ## 1                        7                       241.00000
    ## 2                        3                         4.00000
    ## 3                      570                       805.00000
    ## 4                        0                         0.00000
    ## 5                        1                         7.12279
    ## 6                      234                        69.09420
    ##   r37.w100.x_male_pituitary_n9 r4.x_male_gonad_m.inc.d8
    ## 1                      56.0000                    20.00
    ## 2                      50.0000                     6.00
    ## 3                     176.0000                   568.00
    ## 4                       0.0000                     3.00
    ## 5                       1.0000                     6.00
    ## 6                      60.1219                   264.74
    ##   r4.x_male_hypothalamus_m.inc.d8 r4.x_male_pituitary_m.inc.d8
    ## 1                        162.0000                     361.4383
    ## 2                          2.0000                      61.0000
    ## 3                        483.0000                    1063.0000
    ## 4                          0.0000                       2.0000
    ## 5                          0.0000                       1.0000
    ## 6                         36.3035                     374.2990
    ##   r41.w99.x_male_gonad_hatch r41.w99.x_male_hypothalamus_hatch
    ## 1                          7                          253.0000
    ## 2                          4                            2.0000
    ## 3                        493                          555.0000
    ## 4                          1                            0.0000
    ## 5                          2                            0.0000
    ## 6                        142                           53.2058
    ##   r41.w99.x_male_pituitary_hatch r45.X_male_gonad_inc.d9
    ## 1                        216.000                   9.000
    ## 2                         53.000                   0.000
    ## 3                        490.000                 481.000
    ## 4                          0.000                   0.000
    ## 5                          0.000                   0.000
    ## 6                        199.285                 281.956
    ##   r45.X_male_pituitary_inc.d9 r45.x_male_hypothalamus_inc.d9
    ## 1                    143.0000                       145.4928
    ## 2                     28.0000                         1.0000
    ## 3                    339.0000                       389.0000
    ## 4                      0.0000                         0.0000
    ## 5                      1.0000                         0.0000
    ## 6                     61.5609                         8.0000
    ##   r49.w189.x_female_gonad_inc.d17 r49.w189.x_female_hypothalamus_inc.d17
    ## 1                              71                                    150
    ## 2                             132                                      4
    ## 3                             306                                    362
    ## 4                              58                                      0
    ## 5                               3                                      1
    ## 6                              55                                     33
    ##   r49.w189.x_female_pituitary_inc.d17 r6.x_female_gonad_control.NYNO
    ## 1                             193.000                             96
    ## 2                              15.000                            101
    ## 3                             489.000                            603
    ## 4                               0.000                             37
    ## 5                               5.000                              1
    ## 6                             258.941                             34
    ##   r6.x_female_hypothalamus_control.NYNO r6.x_female_pituitary_control
    ## 1                                   226                            15
    ## 2                                     4                             3
    ## 3                                   474                            41
    ## 4                                     0                             0
    ## 5                                     0                             1
    ## 6                                    24                            14
    ##   r72.y83.x_male_gonad_hatch r72.y83.x_male_hypothalamus_hatch
    ## 1                          2                          488.4638
    ## 2                          0                            3.0000
    ## 3                        521                          589.0000
    ## 4                         19                            0.0000
    ## 5                          3                            2.0000
    ## 6                        181                           70.9725
    ##   r72.y83.x_male_pituitary_hatch r73.g127.x_female_gonad_inc.d3
    ## 1                        132.000                            252
    ## 2                         31.000                            286
    ## 3                        666.000                           2155
    ## 4                          1.000                            177
    ## 5                          2.000                              2
    ## 6                        145.901                            202
    ##   r73.g127.x_female_hypothalamus_inc.d3 r73.g127.x_female_pituitary_inc.d3
    ## 1                              229.0000                                103
    ## 2                                2.0000                                 55
    ## 3                              563.0000                                226
    ## 4                                0.0000                                  0
    ## 5                                3.0000                                  1
    ## 6                               56.2033                                 86
    ##   r81.x_male_gonad_inc.prolong r81.x_male_hypothalamus_inc.prolong
    ## 1                        3.000                            240.0000
    ## 2                        0.000                              5.0000
    ## 3                      437.000                            474.0000
    ## 4                        0.000                              0.0000
    ## 5                        0.000                              0.0000
    ## 6                      115.944                             21.1408
    ##   r81.x_male_pituitary_inc.prolong r83.g45_female_gonad_bldg
    ## 1                         112.0000                        22
    ## 2                           9.0000                         5
    ## 3                         194.0000                       245
    ## 4                           0.0000                         1
    ## 5                           1.0000                         0
    ## 6                          76.7789                        36
    ##   r83.g45_female_hypothalamus_bldg r83.g45_female_pituitary_bldg
    ## 1                              101                       288.000
    ## 2                                0                        42.000
    ## 3                              225                       498.000
    ## 4                                0                         0.000
    ## 5                                0                         4.000
    ## 6                               15                       128.952
    ##   r84.x_male_gonad_m.inc.d9 r84.x_male_hypothalamus_m.inc.d9
    ## 1                     5.000                              102
    ## 2                     2.000                                2
    ## 3                   511.000                              309
    ## 4                     0.000                                0
    ## 5                     1.000                                0
    ## 6                   124.916                               22
    ##   r84.x_male_pituitary_m.inc.d9 r85.g39_male_gonad_m.inc.d9.NYNO
    ## 1                            81                         14.00000
    ## 2                            68                          3.00000
    ## 3                           205                        844.00000
    ## 4                             0                          1.00000
    ## 5                             1                          1.02191
    ## 6                            49                        335.00000
    ##   r85.g39_male_hypothalamus_m.inc.d9 r85.g39_male_pituitary_m.inc.d9.NYNO
    ## 1                            72.0000                              283.320
    ## 2                             0.0000                              138.000
    ## 3                           247.0000                             1034.000
    ## 4                             0.0000                                0.000
    ## 5                             1.0000                                1.000
    ## 6                            11.0597                              290.658
    ##   r90.x_male_gonad_inc.prolong r90.x_male_hypothalamus_inc.prolong
    ## 1                            3                            299.0000
    ## 2                            0                              1.0000
    ## 3                          399                            687.0000
    ## 4                            0                              0.0000
    ## 5                            2                              1.0000
    ## 6                          148                             49.0618
    ##   r90.x_male_pituitary_inc.prolong r95.blu99_female_gonad_n9
    ## 1                         153.0000                       198
    ## 2                          26.0000                       400
    ## 3                         322.0000                       777
    ## 4                           0.0000                        32
    ## 5                           0.0000                         7
    ## 6                          77.5533                        80
    ##   r95.blu99_female_hypothalamus_n9 r95.blu99_female_pituitary_n9
    ## 1                        370.00000                           116
    ## 2                          0.00000                            12
    ## 3                        817.00000                           273
    ## 4                          0.00000                             0
    ## 5                          3.79307                             2
    ## 6                         45.30540                            61
    ##   s.o.pk_female_gonad_lay s.o.pk_female_hypothalamus_lay
    ## 1                     501                       289.0000
    ## 2                     338                         2.0000
    ## 3                     587                       521.0000
    ## 4                     132                         0.0000
    ## 5                       2                         0.0000
    ## 6                      57                        28.3354
    ##   s.o.pk_female_pituitary_lay s.pu.pk_female_gonad_prolong
    ## 1                     282.000                     137.8081
    ## 2                      42.000                     176.0000
    ## 3                     615.000                     602.0000
    ## 4                       6.000                      17.0000
    ## 5                       4.000                       1.0000
    ## 6                     147.841                      93.0000
    ##   s.pu.pk_female_hypothalamus_prolong s.pu.pk_female_pituitary_prolong
    ## 1                            239.0000                          249.000
    ## 2                              3.0000                           65.000
    ## 3                            710.0000                          659.000
    ## 4                              0.0000                            1.000
    ## 5                              1.0000                            0.000
    ## 6                             65.7126                          243.232
    ##   s.pu148.blk.r_male_gonad_bldg s.pu148.blk.r_male_hypothalamus_bldg
    ## 1                         4.000                             138.0000
    ## 2                         1.000                               2.0000
    ## 3                       357.000                             241.0000
    ## 4                         0.000                               0.0000
    ## 5                         0.000                               0.0000
    ## 6                       107.891                              16.0806
    ##   s.pu148.blk.r_male_pituitary_bldg s.x.ATLAS_female_gonad_control
    ## 1                           70.0000                             68
    ## 2                           47.0000                            105
    ## 3                          255.0000                            441
    ## 4                            0.0000                            546
    ## 5                            1.0000                              9
    ## 6                           50.4172                             25
    ##   s.x.ATLAS_female_hypothalamus_control s.x.ATLAS_female_pituitary_control
    ## 1                                    74                                 15
    ## 2                                     2                                  0
    ## 3                                   277                                110
    ## 4                                     0                                  0
    ## 5                                     0                                  1
    ## 6                                    26                                 13
    ##   s002.x_female_gonad_m.inc.d8 s002.x_female_hypothalamus_m.inc.d8
    ## 1                     214.1315                            172.5488
    ## 2                     392.0000                              1.0000
    ## 3                    1263.0000                            433.0000
    ## 4                     105.0000                              0.0000
    ## 5                       0.0000                              0.0000
    ## 6                     287.9320                             22.0661
    ##   s002.x_female_pituitary_m.inc.d8 s038.g.d.blk_female_gonad_m.inc.d17
    ## 1                          205.000                           189.45850
    ## 2                           20.000                           415.00000
    ## 3                          591.000                          1451.00000
    ## 4                            0.000                            50.00000
    ## 5                            7.000                            11.92053
    ## 6                          201.708                           150.00000
    ##   s038.g.d.blk_female_hypothalamus_m.inc.d17
    ## 1                                        137
    ## 2                                          1
    ## 3                                        287
    ## 4                                          0
    ## 5                                          1
    ## 6                                         25
    ##   s038.g.d.blk_female_pituitary_m.inc.d17 s044.blk.d.r_male_gonad_m.inc.d8
    ## 1                                 454.281                           26.000
    ## 2                                 200.000                            5.000
    ## 3                                 711.000                         1552.000
    ## 4                                   0.000                            8.000
    ## 5                                   3.000                            5.000
    ## 6                                 239.919                          753.463
    ##   s044.blk.d.r_male_hypothalamus_m.inc.d8
    ## 1                                115.0000
    ## 2                                  2.0000
    ## 3                                212.0000
    ## 4                                  0.0000
    ## 5                                  2.0000
    ## 6                                 21.2082
    ##   s044.blk.d.r_male_pituitary_m.inc.d8 s062.d.blk.g_male_gonad_m.hatch
    ## 1                            259.00000                           32.00
    ## 2                             69.00000                            4.00
    ## 3                            586.00000                          941.00
    ## 4                              0.00000                            0.00
    ## 5                              5.92508                            4.00
    ## 6                            131.74600                          299.97
    ##   s062.d.blk.g_male_hypothalamus_m.hatch
    ## 1                               164.0000
    ## 2                                 1.0000
    ## 3                               546.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                11.1595
    ##   s062.d.blk.g_male_pituitary_m.hatch s063.d.blk.l_female_gonad_bldg
    ## 1                             650.000                       67.00000
    ## 2                             151.000                       71.00000
    ## 3                            1055.000                      303.00000
    ## 4                               0.000                      102.00000
    ## 5                               7.000                        2.01395
    ## 6                             321.659                      129.82500
    ##   s063.d.blk.l_female_hypothalamus_bldg s063.d.blk.l_female_pituitary_bldg
    ## 1                              150.0000                            69.0000
    ## 2                                3.0000                            18.0000
    ## 3                              290.0000                           272.0000
    ## 4                                0.0000                             0.0000
    ## 5                                0.0000                             0.0000
    ## 6                               16.0577                            80.9741
    ##   s064.g.blk.pu_male_gonad_m.inc.d17
    ## 1                             36.000
    ## 2                              4.000
    ## 3                           1334.000
    ## 4                              6.000
    ## 5                              5.000
    ## 6                            479.961
    ##   s064.g.blk.pu_male_hypothalamus_m.inc.d17
    ## 1                                  118.0000
    ## 2                                    0.0000
    ## 3                                  402.0000
    ## 4                                    0.0000
    ## 5                                    1.0000
    ## 6                                   13.0267
    ##   s064.g.blk.pu_male_pituitary_m.inc.d17 s065.l.d.o_male_gonad_bldg
    ## 1                                375.369                         11
    ## 2                                 38.000                          4
    ## 3                                611.000                        462
    ## 4                                  0.000                          0
    ## 5                                  3.000                          7
    ## 6                                216.276                        539
    ##   s065.l.d.o_male_hypothalamus_bldg s065.l.d.o_male_pituitary_bldg
    ## 1                               167                        351.000
    ## 2                                 2                         35.000
    ## 3                               450                        967.000
    ## 4                                 0                          0.000
    ## 5                                 1                          9.000
    ## 6                                17                        177.632
    ##   s066.l.d.r_male_gonad_bldg s066.l.d.r_male_hypothalamus_bldg
    ## 1                         12                               102
    ## 2                          0                                 3
    ## 3                        515                               193
    ## 4                          0                                 0
    ## 5                          2                                 0
    ## 6                        172                                16
    ##   s066.l.d.r_male_pituitary_bldg s067.o.l.y_male_gonad_m.inc.d8
    ## 1                        84.0000                         19.000
    ## 2                        10.0000                          0.000
    ## 3                       255.0000                       1193.000
    ## 4                         0.0000                          1.000
    ## 5                         0.0000                          4.000
    ## 6                        63.0944                        422.807
    ##   s067.o.l.y_male_hypothalamus_m.inc.d8 s067.o.l.y_male_pituitary_m.inc.d8
    ## 1                              186.0000                            371.000
    ## 2                                2.0000                             50.000
    ## 3                              372.0000                            620.000
    ## 4                                0.0000                              0.000
    ## 5                                1.0000                             13.000
    ## 6                               37.5402                            209.668
    ##   s069.pk.pu.g_female_gonad_m.inc.d3
    ## 1                           644.3560
    ## 2                           212.0000
    ## 3                           914.0000
    ## 4                            44.0000
    ## 5                            20.0000
    ## 6                            62.8809
    ##   s069.pk.pu.g_female_hypothalamus_m.inc.d3
    ## 1                                   193.000
    ## 2                                     1.000
    ## 3                                   451.000
    ## 4                                     0.000
    ## 5                                     0.000
    ## 6                                    25.141
    ##   s069.pk.pu.g_female_pituitary_m.inc.d3 s071.pu.g.pk_male_gonad_m.inc.d3
    ## 1                                213.000                          19.0324
    ## 2                                 69.000                           3.0000
    ## 3                                597.000                         942.0000
    ## 4                                  0.000                           2.0000
    ## 5                                  4.000                           1.0000
    ## 6                                201.018                         443.0000
    ##   s071.pu.g.pk_male_hypothalamus_m.inc.d3
    ## 1                                113.0000
    ## 2                                  0.0000
    ## 3                                324.0000
    ## 4                                  0.0000
    ## 5                                  4.0000
    ## 6                                 17.0477
    ##   s071.pu.g.pk_male_pituitary_m.inc.d3 s089.blk.pk.pu_female_gonad_extend
    ## 1                           421.000000                                244
    ## 2                             5.000000                                222
    ## 3                           622.000000                                557
    ## 4                             0.000000                                 21
    ## 5                             0.997325                                  1
    ## 6                           166.419000                                 77
    ##   s089.blk.pk.pu_female_hypothalamus_extend
    ## 1                                  588.0000
    ## 2                                    2.0000
    ## 3                                 1094.0000
    ## 4                                    0.0000
    ## 5                                    1.0000
    ## 6                                   78.3128
    ##   s089.blk.pk.pu_female_pituitary_extend s090.blk.pk.w_male_gonad_m.hatch
    ## 1                              328.00000                               14
    ## 2                               30.00000                                2
    ## 3                              996.00000                              915
    ## 4                                0.00000                               10
    ## 5                                5.09787                                5
    ## 6                              321.05000                              368
    ##   s090.blk.pk.w_male_hypothalamus_m.hatch
    ## 1                                  260.00
    ## 2                                    4.00
    ## 3                                  787.00
    ## 4                                    0.00
    ## 5                                    3.00
    ## 6                                   65.82
    ##   s090.blk.pk.w_male_pituitary_m.hatch s091.blk.r.g_female_gonad_m.inc.d8
    ## 1                              207.000                             103.00
    ## 2                               99.000                             406.00
    ## 3                              458.000                            3222.47
    ## 4                                0.000                             299.00
    ## 5                                0.000                               6.00
    ## 6                              225.737                             111.00
    ##   s091.blk.r.g_female_hypothalamus_m.inc.d8
    ## 1                                   486.000
    ## 2                                     5.000
    ## 3                                   883.000
    ## 4                                     7.000
    ## 5                                     1.000
    ## 6                                    48.186
    ##   s091.blk.r.g_female_pituitary_m.inc.d8 s092.blk.r.o_female_gonad_bldg
    ## 1                                379.000                        46.0000
    ## 2                                 30.000                       140.0000
    ## 3                                834.000                       518.0000
    ## 4                                  0.000                        32.0000
    ## 5                                  3.000                         3.0000
    ## 6                                256.668                        76.8905
    ##   s092.blk.r.o_female_hypothalamus_bldg s092.blk.r.o_female_pituitary_bldg
    ## 1                                   514                           103.0000
    ## 2                                     6                             8.0000
    ## 3                                   643                           218.0000
    ## 4                                     0                             0.0000
    ## 5                                     3                             0.0000
    ## 6                                    55                            62.0045
    ##   s093.blk.y.d_female_gonad_m.inc.d3
    ## 1                            267.000
    ## 2                            624.000
    ## 3                           1529.000
    ## 4                            586.000
    ## 5                              7.000
    ## 6                            242.925
    ##   s093.blk.y.d_female_hypothalamus_m.inc.d3
    ## 1                                   87.0000
    ## 2                                    1.0000
    ## 3                                  292.0000
    ## 4                                    0.0000
    ## 5                                    0.0000
    ## 6                                   15.0684
    ##   s093.blk.y.d_female_pituitary_m.inc.d3 s095.g.blk.o_female_gonad_lay
    ## 1                                356.000                     133.00000
    ## 2                                 28.000                     282.00000
    ## 3                                640.000                     572.00000
    ## 4                                  0.000                      71.00000
    ## 5                                  2.000                       5.08058
    ## 6                                189.069                      56.00000
    ##   s095.g.blk.o_female_hypothalamus_lay s095.g.blk.o_female_pituitary_lay
    ## 1                             506.0000                           238.000
    ## 2                               5.0000                            54.000
    ## 3                             985.0000                           699.000
    ## 4                               1.0000                             0.000
    ## 5                               1.0000                             3.000
    ## 6                              86.7824                           150.369
    ##   s096.g.blk.pk_male_gonad_m.inc.d17
    ## 1                             80.000
    ## 2                              2.000
    ## 3                            958.000
    ## 4                              1.000
    ## 5                              3.000
    ## 6                            616.645
    ##   s096.g.blk.pk_male_hypothalamus_m.inc.d17
    ## 1                                 87.000000
    ## 2                                  2.000000
    ## 3                                302.000000
    ## 4                                  0.000000
    ## 5                                  0.574588
    ## 6                                 19.034600
    ##   s096.g.blk.pk_male_pituitary_m.inc.d17
    ## 1                              234.00000
    ## 2                               19.00000
    ## 3                              476.00000
    ## 4                                0.00000
    ## 5                                3.04165
    ## 6                              175.44100
    ##   s100.l.pk.blk_female_gonad_m.inc.d17
    ## 1                             89.00000
    ## 2                            106.00000
    ## 3                            529.00000
    ## 4                             55.00000
    ## 5                              4.01436
    ## 6                             98.00000
    ##   s100.l.pk.blk_female_hypothalamus_m.inc.d17
    ## 1                                         106
    ## 2                                           0
    ## 3                                         295
    ## 4                                           0
    ## 5                                           2
    ## 6                                          17
    ##   s100.l.pk.blk_female_pituitary_m.inc.d17
    ## 1                                 377.0000
    ## 2                                  90.0000
    ## 3                                 698.0000
    ## 4                                   0.0000
    ## 5                                   3.0488
    ## 6                                 322.4510
    ##   s103.y.pk.blk_female_gonad_m.inc.d3
    ## 1                            146.0000
    ## 2                             38.0000
    ## 3                            509.0000
    ## 4                             30.0000
    ## 5                              5.0000
    ## 6                             49.8351
    ##   s103.y.pk.blk_female_hypothalamus_m.inc.d3
    ## 1                                   112.0000
    ## 2                                     1.0000
    ## 3                                   285.0000
    ## 4                                     0.0000
    ## 5                                     1.0000
    ## 6                                    17.0905
    ##   s103.y.pk.blk_female_pituitary_m.inc.d3 s136.d.w.o_female_gonad_lay
    ## 1                                 162.000                         279
    ## 2                                  10.000                          78
    ## 3                                 395.000                         381
    ## 4                                   2.000                         150
    ## 5                                   2.000                           2
    ## 6                                  99.456                          57
    ##   s136.d.w.o_female_hypothalamus_lay s136.d.w.o_female_pituitary_lay
    ## 1                            307.000                        277.5261
    ## 2                              7.000                         14.0000
    ## 3                            537.000                        515.0000
    ## 4                              0.000                          0.0000
    ## 5                              0.000                          4.0000
    ## 6                             60.513                        227.4950
    ##   s137.g.blk.y_female_gonad_extend s137.g.blk.y_female_hypothalamus_extend
    ## 1                              268                                493.0000
    ## 2                              101                                  3.0000
    ## 3                              490                                735.0000
    ## 4                                3                                  0.0000
    ## 5                                1                                  0.0000
    ## 6                               89                                 43.1931
    ##   s137.g.blk.y_female_pituitary_extend s139.l.blk.w_male_gonad_m.inc.d3
    ## 1                             404.0919                            35.00
    ## 2                              29.0000                            11.00
    ## 3                             904.0000                          1245.00
    ## 4                               0.0000                             0.00
    ## 5                               4.0000                             5.00
    ## 6                             190.7940                           455.87
    ##   s139.l.blk.w_male_hypothalamus_m.inc.d3
    ## 1                                 87.0000
    ## 2                                  1.0000
    ## 3                                198.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 14.0403
    ##   s139.l.blk.w_male_pituitary_m.inc.d3 s142.o.pk.pu_female_gonad_lay
    ## 1                              370.000                            48
    ## 2                               28.000                             1
    ## 3                              401.000                           616
    ## 4                                0.000                             7
    ## 5                                2.000                             0
    ## 6                              110.345                            23
    ##   s142.o.pk.pu_female_hypothalamus_lay s142.o.pk.pu_female_pituitary_lay
    ## 1                              445.438                          298.0856
    ## 2                                2.000                           32.0000
    ## 3                              844.000                          401.0000
    ## 4                                0.000                            0.0000
    ## 5                                2.000                            2.0000
    ## 6                               67.000                          164.9370
    ##   s150.w.g.blk_male_gonad_lay s150.w.g.blk_male_hypothalamus_lay
    ## 1                   14.000000                           285.0000
    ## 2                    3.000000                             2.0000
    ## 3                  647.000000                           832.0000
    ## 4                    0.000000                             0.0000
    ## 5                    0.457772                             0.0000
    ## 6                  259.000000                            71.1904
    ##   s150.w.g.blk_male_pituitary_lay s175.blk.pu.pk_male_gonad_m.inc.d3
    ## 1                         363.000                             16.000
    ## 2                         101.000                              1.000
    ## 3                         661.000                           1029.000
    ## 4                           0.000                              0.000
    ## 5                           0.000                              8.000
    ## 6                         165.657                            362.863
    ##   s175.blk.pu.pk_male_hypothalamus_m.inc.d3
    ## 1                                   98.0000
    ## 2                                    1.0000
    ## 3                                  230.0000
    ## 4                                    0.0000
    ## 5                                    0.0000
    ## 6                                   16.0568
    ##   s175.blk.pu.pk_male_pituitary_m.inc.d3 s176.blk.pu.r_female_gonad_lay
    ## 1                                248.000                        66.1214
    ## 2                                 56.000                        77.0000
    ## 3                                465.000                       417.0000
    ## 4                                  0.000                        23.0000
    ## 5                                  0.000                         0.0000
    ## 6                                212.525                        32.0000
    ##   s176.blk.pu.r_female_hypothalamus_lay s176.blk.pu.r_female_pituitary_lay
    ## 1                               353.442                            387.322
    ## 2                                 1.000                             45.000
    ## 3                               522.000                            819.000
    ## 4                                 0.000                              0.000
    ## 5                                 1.000                              1.000
    ## 6                                44.398                            237.888
    ##   s186.l.o.pk_female_gonad_extend s186.l.o.pk_female_hypothalamus_extend
    ## 1                          125.00                               207.0000
    ## 2                           47.00                                 0.0000
    ## 3                          378.00                               561.0000
    ## 4                           18.00                                 0.0000
    ## 5                            1.00                                 1.0000
    ## 6                          123.85                                26.1634
    ##   s186.l.o.pk_female_pituitary_extend s187.l.o.r_male_gonad_n9
    ## 1                             457.000                  41.1015
    ## 2                              62.000                   6.0000
    ## 3                             699.000                1173.0000
    ## 4                               3.000                   0.0000
    ## 5                               4.000                   6.0000
    ## 6                             184.162                 630.9350
    ##   s187.l.o.r_male_hypothalamus_n9 s187.l.o.r_male_pituitary_n9
    ## 1                        218.0000                     246.0000
    ## 2                          1.0000                      51.0000
    ## 3                        728.0000                     476.0000
    ## 4                          0.0000                       0.0000
    ## 5                          0.0000                       1.0000
    ## 6                         18.0298                      90.5408
    ##   s243.blk.pk.r_male_gonad_lay s243.blk.pk.r_male_hypothalamus_lay
    ## 1                       18.000                            312.0000
    ## 2                        1.000                              5.0000
    ## 3                     1037.000                            670.0000
    ## 4                        0.000                              0.0000
    ## 5                       10.000                              2.0000
    ## 6                      571.952                             54.2584
    ##   s243.blk.pk.r_male_pituitary_lay s333.y.blk.pk_female_gonad_m.hatch
    ## 1                          335.000                                 73
    ## 2                           82.000                                138
    ## 3                         1083.000                                412
    ## 4                            0.000                                  3
    ## 5                            0.000                                  6
    ## 6                          263.717                                 63
    ##   s333.y.blk.pk_female_hypothalamus_m.hatch
    ## 1                                  252.0000
    ## 2                                    6.0000
    ## 3                                  594.0000
    ## 4                                    0.0000
    ## 5                                    3.0000
    ## 6                                   41.2738
    ##   s333.y.blk.pk_female_pituitary_m.hatch w191.r1_female_gonad_control
    ## 1                               194.2892                           11
    ## 2                                32.0000                           23
    ## 3                               338.0000                           72
    ## 4                                 0.0000                            8
    ## 5                                 1.0000                            1
    ## 6                               120.9890                            7
    ##   w191.r1_female_hypothalamus_control w191.r1_female_pituitary_control
    ## 1                                  37                               37
    ## 2                                   0                                7
    ## 3                                  35                               43
    ## 4                                   0                                0
    ## 5                                   0                                1
    ## 6                                   2                               14
    ##   w34.x_male_gonad_inc.d9 w34.x_male_hypothalamus_inc.d9
    ## 1                   6.000                        215.000
    ## 2                   3.000                          2.000
    ## 3                 475.000                        671.000
    ## 4                   0.000                          0.000
    ## 5                   3.000                          1.000
    ## 6                 216.756                         40.241
    ##   w34.x_male_pituitary_inc.d9 x.blk.blk.ATLAS_male_gonad_control
    ## 1                     50.0000                            24.0747
    ## 2                     31.0000                             2.0000
    ## 3                    235.0000                           932.0000
    ## 4                      2.0000                             0.0000
    ## 5                      1.0000                             6.0000
    ## 6                     54.6582                           212.9450
    ##   x.blk.blk.ATLAS_male_hypothalamus_control
    ## 1                                 147.00000
    ## 2                                   2.00000
    ## 3                                 129.00000
    ## 4                                   0.00000
    ## 5                                   0.00000
    ## 6                                   8.05413
    ##   x.blk.blk.ATLAS_male_pituitary_control x.blk16_male_gonad_n9.NYNO
    ## 1                               30.00000                    43.3122
    ## 2                               25.00000                     4.0000
    ## 3                               54.00000                   868.0000
    ## 4                                0.00000                     0.0000
    ## 5                                1.00000                     6.0000
    ## 6                                6.29407                   287.9820
    ##   x.blk16_male_hypothalamus_n9.NYNO x.blk16_male_pituitary_n9
    ## 1                          157.5762                   222.000
    ## 2                            3.0000                    75.000
    ## 3                          560.0000                   719.000
    ## 4                            0.0000                     0.000
    ## 5                            0.0000                     4.000
    ## 6                           33.0000                   221.994
    ##   x.blk20_female_gonad_prolong x.blk20_female_hypothalamus_prolong
    ## 1                       21.000                            102.0000
    ## 2                        2.000                              2.0000
    ## 3                      756.000                            292.0000
    ## 4                       30.000                              0.0000
    ## 5                        2.000                              1.0000
    ## 6                      425.959                             45.1597
    ##   x.blk20_female_pituitary_prolong x.blk23_male_gonad_m.n2
    ## 1                          165.000                       5
    ## 2                           78.000                       0
    ## 3                          508.000                     371
    ## 4                            0.000                       0
    ## 5                            3.000                       0
    ## 6                          270.421                     134
    ##   x.blk23_male_hypothalamus_m.n2 x.blk23_male_pituitary_m.n2
    ## 1                       131.0000                    135.0000
    ## 2                         1.0000                    106.0000
    ## 3                       267.0000                    297.0000
    ## 4                         0.0000                      0.0000
    ## 5                         0.0000                      4.0000
    ## 6                        21.0632                     65.7529
    ##   x.blu.o.ATLAS_male_pituitary_control x.blu101.w43_female_gonad_inc.d9
    ## 1                                  119                               27
    ## 2                                   36                               78
    ## 3                                  290                              228
    ## 4                                    0                               28
    ## 5                                    2                                1
    ## 6                                   32                               27
    ##   x.blu101.w43_female_hypothalamus_inc.d9
    ## 1                                 81.0000
    ## 2                                  1.0000
    ## 3                                163.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 15.0698
    ##   x.blu101.w43_female_pituitary_inc.d9 x.blu102.w105_female_gonad_inc.d3
    ## 1                              292.385                                28
    ## 2                               26.000                                61
    ## 3                              736.000                               387
    ## 4                                0.000                                52
    ## 5                                0.000                                 1
    ## 6                              162.969                                20
    ##   x.blu102.w105_female_hypothalamus_inc.d3
    ## 1                                 331.6843
    ## 2                                   6.0000
    ## 3                                 666.0000
    ## 4                                   0.0000
    ## 5                                   1.0000
    ## 6                                  55.3019
    ##   x.blu102.w105_female_pituitary_inc.d3
    ## 1                               184.000
    ## 2                                22.000
    ## 3                               892.000
    ## 4                                 0.000
    ## 5                                 5.000
    ## 6                               160.137
    ##   x.blu106.o153_male_gonad_inc.d9.NYNO
    ## 1                               18.000
    ## 2                                3.000
    ## 3                              954.000
    ## 4                                1.000
    ## 5                                7.000
    ## 6                              274.692
    ##   x.blu106.o153_male_hypothalamus_inc.d9
    ## 1                               223.3344
    ## 2                                 4.0000
    ## 3                               514.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                32.1358
    ##   x.blu106.o153_male_pituitary_inc.d9.NYNO x.blu109.w121_female_gonad_n5
    ## 1                                 319.0000                            70
    ## 2                                  18.0000                            33
    ## 3                                 552.0000                           279
    ## 4                                   0.0000                            21
    ## 5                                   0.0000                             0
    ## 6                                  95.3834                            43
    ##   x.blu109.w121_female_hypothalamus_n5 x.blu109.w121_female_pituitary_n5
    ## 1                             297.0000                           112.000
    ## 2                               7.0000                            25.000
    ## 3                             536.0000                           230.000
    ## 4                               0.0000                             0.000
    ## 5                               3.0000                             1.000
    ## 6                              29.2823                            65.915
    ##   x.blu116.w107_female_gonad_inc.d17
    ## 1                                161
    ## 2                                 29
    ## 3                                258
    ## 4                                 13
    ## 5                                  0
    ## 6                                  7
    ##   x.blu116.w107_female_hypothalamus_inc.d17
    ## 1                                  315.0000
    ## 2                                    4.0000
    ## 3                                  529.0000
    ## 4                                    0.0000
    ## 5                                    0.0000
    ## 6                                   42.2818
    ##   x.blu116.w107_female_pituitary_inc.d17.NYNO
    ## 1                                     160.000
    ## 2                                      38.000
    ## 3                                     487.000
    ## 4                                       0.000
    ## 5                                       0.000
    ## 6                                     135.431
    ##   x.blu117.w89_male_gonad_inc.d17 x.blu117.w89_male_hypothalamus_inc.d17
    ## 1                               4                               263.0000
    ## 2                               1                                 0.0000
    ## 3                             441                               346.0000
    ## 4                               2                                 0.0000
    ## 5                               1                                 3.0000
    ## 6                             134                                23.1035
    ##   x.blu117.w89_male_pituitary_inc.d17 x.blu122.r66_female_gonad_inc.d9
    ## 1                                  87                               80
    ## 2                                  12                               93
    ## 3                                 190                              758
    ## 4                                   0                                8
    ## 5                                   0                                1
    ## 6                                  52                               50
    ##   x.blu122.r66_female_hypothalamus_inc.d9
    ## 1                                285.0000
    ## 2                                  1.0000
    ## 3                                753.0000
    ## 4                                  0.0000
    ## 5                                  1.0000
    ## 6                                 28.0435
    ##   x.blu122.r66_female_pituitary_inc.d9 x.blu23.w14_male_gonad_n9
    ## 1                             158.0000                     11.00
    ## 2                               6.0000                      3.00
    ## 3                             274.0000                    351.00
    ## 4                               0.0000                      1.00
    ## 5                               0.0000                      0.00
    ## 6                              75.5002                    123.42
    ##   x.blu23.w14_male_hypothalamus_n9 x.blu23.w14_male_pituitary_n9
    ## 1                         115.0000                       198.000
    ## 2                           0.0000                        59.000
    ## 3                         255.0000                       623.000
    ## 4                           0.0000                         0.000
    ## 5                           1.0000                         3.000
    ## 6                          11.2852                       190.169
    ##   x.blu25_female_gonad_m.inc.d9 x.blu25_female_hypothalamus_m.inc.d9
    ## 1                            12                                  331
    ## 2                             6                                    1
    ## 3                           247                                  604
    ## 4                             0                                    0
    ## 5                             1                                    1
    ## 6                            30                                   60
    ##   x.blu25_female_pituitary_m.inc.d9 x.blu30_male_gonad_n5
    ## 1                           45.0000                 20.00
    ## 2                           29.0000                  0.00
    ## 3                          235.0000               1211.00
    ## 4                            0.0000                  0.00
    ## 5                            0.0000                 11.00
    ## 6                           98.5696                623.71
    ##   x.blu30_male_hypothalamus_n5 x.blu30_male_pituitary_n5
    ## 1                     115.0000                  134.2518
    ## 2                       1.0000                   36.0000
    ## 3                     349.0000                  423.0000
    ## 4                       0.0000                    0.0000
    ## 5                       1.0000                    3.0000
    ## 6                      35.2691                  262.6560
    ##   x.blu42.o28_male_gonad_inc.d3 x.blu42.o28_male_hypothalamus_inc.d3.NYNO
    ## 1                             3                                       233
    ## 2                             1                                         0
    ## 3                           464                                       534
    ## 4                             0                                         0
    ## 5                             2                                         0
    ## 6                           136                                        43
    ##   x.blu42.o28_male_pituitary_inc.d3 x.blu43.g132_female_gonad_n9
    ## 1                           302.000                           29
    ## 2                            55.000                           36
    ## 3                           918.000                          341
    ## 4                             0.000                           16
    ## 5                             1.000                            2
    ## 6                           238.091                           81
    ##   x.blu43.g132_female_hypothalamus_n9 x.blu43.g132_female_pituitary_n9
    ## 1                            231.0000                               56
    ## 2                              3.0000                               27
    ## 3                            711.0000                              169
    ## 4                              0.0000                                0
    ## 5                              2.0000                                4
    ## 6                             66.8668                               56
    ##   x.blu57_male_gonad_prolong x.blu57_male_hypothalamus_prolong
    ## 1                     17.000                          292.0000
    ## 2                      2.000                            3.0000
    ## 3                   1114.000                          254.0000
    ## 4                      0.000                            0.0000
    ## 5                      6.000                            0.0000
    ## 6                    437.963                           40.5562
    ##   x.blu57_male_pituitary_prolong x.blu6.y80_female_gonad_lay
    ## 1                      277.00000                     80.0000
    ## 2                       27.00000                    148.0000
    ## 3                      379.00000                    603.0000
    ## 4                        0.00000                      2.0000
    ## 5                        2.80812                      6.0000
    ## 6                      193.26500                     46.3754
    ##   x.blu6.y80_female_hypothalamus_lay x.blu6.y80_female_pituitary_lay
    ## 1                           406.0000                              89
    ## 2                             7.0000                              19
    ## 3                           646.0000                             321
    ## 4                             0.0000                               0
    ## 5                             2.0000                               0
    ## 6                            49.4926                              37
    ##   x.blu69_male_gonad_extend x.blu69_male_hypothalamus_extend
    ## 1                     9.000                          233.000
    ## 2                     3.000                            4.000
    ## 3                  1038.000                          561.000
    ## 4                     0.000                            0.000
    ## 5                     9.000                            0.000
    ## 6                   358.874                           43.096
    ##   x.blu69_male_pituitary_extend x.blu73_male_gonad_extend
    ## 1                       295.000                    14.000
    ## 2                        80.000                     4.000
    ## 3                       642.000                  1967.000
    ## 4                         0.000                     4.000
    ## 5                         0.000                    10.000
    ## 6                       270.535                   534.797
    ##   x.blu73_male_hypothalamus_extend x.blu73_male_pituitary_extend
    ## 1                         218.0000                       210.000
    ## 2                           4.0000                        52.000
    ## 3                         622.0000                       362.000
    ## 4                           0.0000                         0.000
    ## 5                           4.0000                         5.000
    ## 6                          63.3916                       148.857
    ##   x.blu76_male_gonad_m.hatch x.blu76_male_hypothalamus_m.hatch
    ## 1                   71.00000                          238.4507
    ## 2                    6.00000                            1.0000
    ## 3                 1851.00000                          629.0000
    ## 4                    1.00000                            0.0000
    ## 5                    6.38926                            0.0000
    ## 6                  463.46000                           69.2807
    ##   x.blu76_male_pituitary_m.hatch x.blu94_female_gonad_m.inc.d17
    ## 1                        205.000                      113.00000
    ## 2                         34.000                      207.00000
    ## 3                        325.000                      829.00000
    ## 4                          0.000                        7.00000
    ## 5                          0.000                       11.88084
    ## 6                        139.204                       70.86650
    ##   x.blu94_female_hypothalamus_m.inc.d17 x.blu94_female_pituitary_m.inc.d17
    ## 1                                    83                            283.000
    ## 2                                     0                            124.000
    ## 3                                   184                            633.000
    ## 4                                     0                              0.000
    ## 5                                     2                              2.000
    ## 6                                    13                            139.258
    ##   x.blu96_female_gonad_m.n2 x.blu96_female_hypothalamus_m.n2
    ## 1                        89                          248.000
    ## 2                       421                            7.000
    ## 3                       820                          524.000
    ## 4                         8                            0.000
    ## 5                         0                            0.000
    ## 6                        85                           41.083
    ##   x.blu96_female_pituitary_m.n2 x.g.ATLAS_male_gonad_control
    ## 1                            76                            0
    ## 2                            10                            0
    ## 3                           203                          212
    ## 4                             0                            0
    ## 5                             1                            0
    ## 6                            49                           19
    ##   x.g.ATLAS_male_hypothalamus_control x.g.ATLAS_male_pituitary_control
    ## 1                            93.00000                         159.0000
    ## 2                             1.00000                          12.0000
    ## 3                           149.00000                         297.0000
    ## 4                             1.00000                           0.0000
    ## 5                             0.00000                           0.0000
    ## 6                             8.05182                          21.7122
    ##   x.g.g.ATLAS_female_gonad_control x.g.g.ATLAS_male_gonad_control
    ## 1                               80                         14.000
    ## 2                               49                          3.000
    ## 3                              542                       1126.000
    ## 4                               14                          0.000
    ## 5                                2                          2.000
    ## 6                               19                        285.787
    ##   x.g.g.ATLAS_male_hypothalamus_control x.g.g.ATLAS_male_pituitary_control
    ## 1                                    16                                 13
    ## 2                                     2                                 58
    ## 3                                   199                                240
    ## 4                                     0                                  0
    ## 5                                     0                                  0
    ## 6                                     3                                 21
    ##   x.g.g.g.ATLAS_male_gonad_control x.g.g.g.ATLAS_male_pituitary_control
    ## 1                            9.000                                  126
    ## 2                            2.000                                   52
    ## 3                         1142.260                                  214
    ## 4                            0.000                                    0
    ## 5                           10.000                                    0
    ## 6                          222.832                                   14
    ##   x.g13.w109_male_gonad_inc.d9 x.g13.w109_male_hypothalamus_inc.d9
    ## 1                            4                            243.0000
    ## 2                            0                              0.0000
    ## 3                          573                            502.0000
    ## 4                            1                              0.0000
    ## 5                            1                              1.0000
    ## 6                          149                             45.5407
    ##   x.g13.w109_male_pituitary_inc.d9 x.g14.w199_male_gonad_inc.d17
    ## 1                         133.0000                             5
    ## 2                           9.0000                             1
    ## 3                         297.0000                           344
    ## 4                           0.0000                             2
    ## 5                           1.0000                             5
    ## 6                          64.2791                           121
    ##   x.g14.w199_male_hypothalamus_inc.d17 x.g14.w199_male_pituitary_inc.d17
    ## 1                             361.0000                          329.3536
    ## 2                               2.0000                           54.0000
    ## 3                             648.0000                          738.0000
    ## 4                               0.0000                            0.0000
    ## 5                               0.0000                            5.0000
    ## 6                              35.8286                          162.6550
    ##   x.g147.blu28_male_gonad_inc.d3 x.g147.blu28_male_hypothalamus_inc.d3
    ## 1                          9.000                              261.0000
    ## 2                          3.000                                5.0000
    ## 3                       1562.000                              985.0000
    ## 4                        124.000                                0.0000
    ## 5                          3.000                                2.0000
    ## 6                        553.985                               43.1875
    ##   x.g147.blu28_male_pituitary_inc.d3 x.g26_female_gonad_m.inc.d8
    ## 1                            106.000                      130.00
    ## 2                             13.000                       84.00
    ## 3                            236.000                      949.00
    ## 4                              0.000                      231.00
    ## 5                              1.000                        3.00
    ## 6                             54.323                      194.43
    ##   x.g26_female_hypothalamus_m.inc.d8 x.g26_female_pituitary_m.inc.d8
    ## 1                            328.000                         296.000
    ## 2                              8.000                          81.000
    ## 3                            592.000                         567.000
    ## 4                              0.000                           0.000
    ## 5                              1.000                           5.000
    ## 6                             46.645                         146.854
    ##   x.g37_female_gonad_n5 x.g37_female_hypothalamus_n5
    ## 1                    42                          261
    ## 2                    52                            2
    ## 3                   205                          648
    ## 4                     2                            0
    ## 5                     1                            3
    ## 6                    41                           43
    ##   x.g37_female_pituitary_n5 x.g4.w50_female_gonad_n9
    ## 1                   80.0000                       55
    ## 2                   37.0000                       38
    ## 3                  225.0000                      261
    ## 4                    0.0000                        4
    ## 5                    0.0000                        2
    ## 6                   61.1526                       61
    ##   x.g4.w50_female_hypothalamus_n9 x.g4.w50_female_pituitary_n9
    ## 1                             155                           85
    ## 2                               1                            9
    ## 3                             482                          195
    ## 4                               0                            0
    ## 5                               0                            2
    ## 6                              36                           43
    ##   x.g43_female_gonad_n5 x.g43_female_hypothalamus_n5
    ## 1                    85                          148
    ## 2                    60                            2
    ## 3                   202                          363
    ## 4                   105                            1
    ## 5                     0                            0
    ## 6                    18                           21
    ##   x.g43_female_pituitary_n5 x.g49_female_gonad_n5
    ## 1                   74.0000                   141
    ## 2                    7.0000                   132
    ## 3                  184.0000                   729
    ## 4                    0.0000                     0
    ## 5                    4.0000                     5
    ## 6                   64.7108                    93
    ##   x.g49_female_hypothalamus_n5 x.g49_female_pituitary_n5.NYNO
    ## 1                     172.0000                        100.000
    ## 2                       2.0000                         11.000
    ## 3                     475.0000                        391.000
    ## 4                       4.0000                          0.000
    ## 5                       4.0000                          2.000
    ## 6                      36.0876                        130.174
    ##   x.g5.w79_male_gonad_m.inc.d3 x.g5.w79_male_hypothalamus_m.inc.d3
    ## 1                        9.000                             98.0000
    ## 2                        0.000                              3.0000
    ## 3                     2103.000                            329.0000
    ## 4                        0.000                              0.0000
    ## 5                        4.000                              0.0000
    ## 6                      479.937                             16.0705
    ##   x.g5.w79_male_pituitary_m.inc.d3 x.g50_female_gonad_inc.prolong
    ## 1                          225.000                            125
    ## 2                           32.000                            138
    ## 3                          519.000                            420
    ## 4                            0.000                             18
    ## 5                            2.000                              1
    ## 6                          138.927                             39
    ##   x.g50_female_hypothalamus_inc.prolong x.g50_female_pituitary_inc.prolong
    ## 1                             234.00000                           101.0000
    ## 2                               1.00000                            20.0000
    ## 3                             606.00000                           313.0000
    ## 4                               0.00000                             0.0000
    ## 5                               0.64101                             0.0000
    ## 6                              63.33810                            62.5335
    ##   x.g70_male_gonad_hatch x.g70_male_hypothalamus_hatch
    ## 1                     10                      252.0000
    ## 2                      0                        4.0000
    ## 3                    440                      623.0000
    ## 4                      0                        0.0000
    ## 5                      4                        0.0000
    ## 6                    151                       31.3563
    ##   x.g70_male_pituitary_hatch x.g9.o166_female_gonad_inc.d9.NYNO
    ## 1                    95.0000                                 79
    ## 2                     9.0000                                 58
    ## 3                   209.0000                                609
    ## 4                     0.0000                                 31
    ## 5                     0.0000                                  1
    ## 6                    74.9382                                 36
    ##   x.g9.o166_female_hypothalamus_inc.d9
    ## 1                                  360
    ## 2                                    4
    ## 3                                  418
    ## 4                                    0
    ## 5                                    0
    ## 6                                   38
    ##   x.g9.o166_female_pituitary_inc.d9.NYNO x.o117_male_gonad_m.inc.d9
    ## 1                               190.0000                    6.00000
    ## 2                                 8.0000                    1.00000
    ## 3                               440.0000                  564.00000
    ## 4                                 0.0000                    0.00000
    ## 5                                 1.0000                    3.66482
    ## 6                                95.7844                  185.00000
    ##   x.o117_male_hypothalamus_m.inc.d9 x.o117_male_pituitary_m.inc.d9
    ## 1                                96                             78
    ## 2                                 1                             74
    ## 3                               243                            226
    ## 4                                 0                              0
    ## 5                                 0                              1
    ## 6                                16                             57
    ##   x.o159.w90_female_gonad_inc.d17 x.o159.w90_female_hypothalamus_inc.d17
    ## 1                              45                                432.000
    ## 2                               9                                  2.000
    ## 3                             220                                549.000
    ## 4                              22                                  0.000
    ## 5                               0                                  1.000
    ## 6                              26                                 74.282
    ##   x.o159.w90_female_pituitary_inc.d17 x.o160.w102_male_gonad_hatch
    ## 1                              62.000                           12
    ## 2                               9.000                            3
    ## 3                             185.000                          363
    ## 4                               0.000                            0
    ## 5                               0.000                            1
    ## 6                              62.627                          151
    ##   x.o160.w102_male_hypothalamus_hatch x.o160.w102_male_pituitary_hatch
    ## 1                                 329                              117
    ## 2                                   3                               20
    ## 3                                 610                              188
    ## 4                                   0                                0
    ## 5                                   0                                0
    ## 6                                  36                               63
    ##   x.o163.w101_male_gonad_inc.d3 x.o163.w101_male_hypothalamus_inc.d3
    ## 1                             1                             360.0000
    ## 2                             0                               8.0000
    ## 3                           435                             664.0000
    ## 4                             0                               0.0000
    ## 5                             1                               2.0000
    ## 6                            96                              47.3857
    ##   x.o163.w101_male_pituitary_inc.d3.NYNO x.o164.w123_male_gonad_n5
    ## 1                                244.000                  28.00000
    ## 2                                 12.000                   6.00000
    ## 3                                540.000                1226.00000
    ## 4                                  1.000                   6.00000
    ## 5                                  1.000                   9.66689
    ## 6                                136.983                 508.00000
    ##   x.o164.w123_male_hypothalamus_n5 x.o164.w123_male_pituitary_n5.NYNO
    ## 1                         282.0000                           162.0000
    ## 2                           5.0000                            24.0000
    ## 3                         521.0000                           398.0000
    ## 4                           0.0000                             0.0000
    ## 5                           1.0000                             2.0000
    ## 6                          37.5231                            96.5617
    ##   x.o171.w45_female_gonad_m.hatch x.o171.w45_female_hypothalamus_m.hatch
    ## 1                        80.00000                               295.0000
    ## 2                       152.00000                                 3.0000
    ## 3                       493.00000                               564.0000
    ## 4                         8.00000                                 0.0000
    ## 5                        12.85527                                 0.0000
    ## 6                        70.00000                                36.1629
    ##   x.o171.w45_female_pituitary_m.hatch x.o175.g21_female_gonad_n5
    ## 1                             271.000                        137
    ## 2                              69.000                        119
    ## 3                             600.000                        648
    ## 4                               0.000                          5
    ## 5                               1.000                          1
    ## 6                             180.606                         99
    ##   x.o175.g21_female_hypothalamus_n5 x.o175.g21_female_pituitary_n5
    ## 1                          170.0000                        316.000
    ## 2                            2.0000                         35.000
    ## 3                          470.0000                        989.000
    ## 4                            0.0000                          0.000
    ## 5                            0.0000                          3.000
    ## 6                           24.0849                        321.965
    ##   x.o2_male_gonad_n9 x.o2_male_hypothalamus_n9 x.o2_male_pituitary_n9
    ## 1                 17                       224                187.000
    ## 2                  0                         5                 31.000
    ## 3               1960                       397                456.000
    ## 4                  3                         0                  0.000
    ## 5                  5                         1                  0.000
    ## 6                352                        30                133.531
    ##   x.o30.g134_male_gonad_bldg x.o30.g134_male_hypothalamus_bldg
    ## 1                   16.00000                               177
    ## 2                    6.00000                                 1
    ## 3                 2493.00000                               365
    ## 4                    4.00000                                 0
    ## 5                   13.02614                                 0
    ## 6                  647.75300                                30
    ##   x.o30.g134_male_pituitary_bldg x.o37.blu50_female_gonad_hatch
    ## 1                         90.000                             40
    ## 2                         36.000                             38
    ## 3                        190.000                            210
    ## 4                          0.000                             10
    ## 5                          0.000                              1
    ## 6                         47.377                             28
    ##   x.o37.blu50_female_hypothalamus_hatch.NYNO
    ## 1                                   445.0000
    ## 2                                     6.0000
    ## 3                                   532.0000
    ## 4                                     0.0000
    ## 5                                     2.0000
    ## 6                                    45.1149
    ##   x.o37.blu50_female_pituitary_hatch x.o40.r70_male_gonad_m.inc.d17
    ## 1                                 94                         38.000
    ## 2                                  8                          7.000
    ## 3                                318                       1515.000
    ## 4                                  0                          8.000
    ## 5                                  0                          9.000
    ## 6                                 85                        608.791
    ##   x.o40.r70_male_hypothalamus_m.inc.d17 x.o40.r70_male_pituitary_m.inc.d17
    ## 1                               545.000                            263.000
    ## 2                                 2.000                            200.000
    ## 3                              1544.070                            559.000
    ## 4                                 0.000                              0.000
    ## 5                                 2.000                              2.000
    ## 6                                76.367                            152.371
    ##   x.o47.y82_male_gonad_inc.d9 x.o47.y82_male_hypothalamus_inc.d9
    ## 1                           3                                187
    ## 2                           1                                  1
    ## 3                         445                                261
    ## 4                           2                                  0
    ## 5                           4                                  0
    ## 6                         120                                 17
    ##   x.o47.y82_male_pituitary_inc.d9 x.o4_male_gonad_m.hatch
    ## 1                         232.000                  18.000
    ## 2                          39.000                   4.000
    ## 3                         620.000                1559.000
    ## 4                           0.000                   1.000
    ## 5                           0.000                   7.000
    ## 6                         161.228                 574.817
    ##   x.o4_male_hypothalamus_m.hatch x.o4_male_pituitary_m.hatch
    ## 1                       274.0000                     426.000
    ## 2                         6.0000                     273.000
    ## 3                       769.0000                     914.000
    ## 4                         0.0000                       0.000
    ## 5                         3.0000                       4.000
    ## 6                        49.1308                     259.865
    ##   x.o61_female_gonad_extend.hatch x.o61_female_hypothalamus_extend.hatch
    ## 1                          92.000                               229.0000
    ## 2                         141.000                                 3.0000
    ## 3                        1154.000                               482.0000
    ## 4                          57.000                                 0.0000
    ## 5                           0.000                                 1.0000
    ## 6                         308.829                                34.8294
    ##   x.o61_female_pituitary_extend.hatch x.o68_male_gonad_n5
    ## 1                             210.000              66.000
    ## 2                              33.000               1.000
    ## 3                             474.000            1465.000
    ## 4                               0.000               1.000
    ## 5                               0.000               3.000
    ## 6                             126.266             454.917
    ##   x.o68_male_hypothalamus_n5 x.o68_male_pituitary_n5 x.o70_female_gonad_n5
    ## 1                   206.0000                 387.000                    46
    ## 2                     4.0000                  46.000                    78
    ## 3                   421.0000                 837.000                   279
    ## 4                     0.0000                   0.000                     2
    ## 5                     1.0000                   3.000                     1
    ## 6                    21.0568                 320.081                    44
    ##   x.o70_female_hypothalamus_n5.NYNO x.o70_female_pituitary_n5
    ## 1                          160.0000                   162.000
    ## 2                            5.0000                    12.000
    ## 3                          453.0000                   714.000
    ## 4                            0.0000                     0.000
    ## 5                            2.0000                     4.000
    ## 6                           41.1347                   205.473
    ##   x.r10.w18_male_gonad_m.inc.d17 x.r10.w18_male_hypothalamus_m.inc.d17
    ## 1                       82.00000                               70.0000
    ## 2                        9.00000                                1.0000
    ## 3                     1939.00000                              288.0000
    ## 4                        0.00000                                0.0000
    ## 5                       10.43282                                0.0000
    ## 6                      752.89100                               21.3397
    ##   x.r10.w18_male_pituitary_m.inc.d17 x.r178_male_gonad_hatch
    ## 1                            205.000                30.00000
    ## 2                            308.000                 3.00000
    ## 3                            406.000              1919.00000
    ## 4                              0.000                 0.00000
    ## 5                              0.000                14.59269
    ## 6                            141.758               679.00000
    ##   x.r178_male_hypothalamus_hatch x.r178_male_pituitary_hatch
    ## 1                       425.0000                     86.0000
    ## 2                         2.0000                     12.0000
    ## 3                       581.0000                    199.0000
    ## 4                         0.0000                      1.0000
    ## 5                         0.0000                      0.0000
    ## 6                        52.3416                     84.8959
    ##   x.r180_female_gonad_m.inc.d3 x.r180_female_hypothalamus_m.inc.d3
    ## 1                      255.000                                  97
    ## 2                      359.000                                   0
    ## 3                     1033.000                                 312
    ## 4                       53.000                                   0
    ## 5                        9.000                                   0
    ## 6                      236.976                                  14
    ##   x.r180_female_pituitary_m.inc.d3 x.r181_male_gonad_n5
    ## 1                          180.000                   11
    ## 2                          155.000                    1
    ## 3                          616.000                  417
    ## 4                            0.000                    0
    ## 5                            4.000                    1
    ## 6                          163.813                  121
    ##   x.r181_male_hypothalamus_n5 x.r181_male_pituitary_n5
    ## 1                    409.2810                       60
    ## 2                      7.0000                        3
    ## 3                   1063.0000                      161
    ## 4                      0.0000                        0
    ## 5                      1.0000                        0
    ## 6                     48.3819                       59
    ##   x.r185_male_gonad_m.inc.d3 x.r185_male_hypothalamus_m.inc.d3
    ## 1                   11.00000                          407.0000
    ## 2                    2.00000                            5.0000
    ## 3                 1062.00000                          637.0000
    ## 4                    0.00000                            0.0000
    ## 5                    6.17202                            3.0000
    ## 6                  459.00000                           39.3343
    ##   x.r185_male_pituitary_m.inc.d3 x.r29.w96_male_gonad_inc.d17
    ## 1                        535.000                            8
    ## 2                        366.000                            0
    ## 3                       1494.000                          371
    ## 4                          0.000                            0
    ## 5                          0.000                            0
    ## 6                        357.461                          210
    ##   x.r29.w96_male_hypothalamus_inc.d17 x.r29.w96_male_pituitary_inc.d17
    ## 1                            450.0000                               81
    ## 2                              8.0000                               14
    ## 3                           1060.0000                              173
    ## 4                              0.0000                                0
    ## 5                              0.0000                                0
    ## 6                             49.2768                               37
    ##   x.r33.w183_female_gonad_inc.d3 x.r33.w183_female_hypothalamus_inc.d3
    ## 1                             49                              198.0000
    ## 2                            126                                1.0000
    ## 3                            347                              537.0000
    ## 4                            319                                0.0000
    ## 5                              1                                0.0000
    ## 6                             42                               35.7369
    ##   x.r33.w183_female_pituitary_inc.d3 x.r39.g10_female_gonad_bldg
    ## 1                            97.0000                          28
    ## 2                            15.0000                           8
    ## 3                           281.0000                         324
    ## 4                             0.0000                          61
    ## 5                             2.0000                           2
    ## 6                            67.5426                          46
    ##   x.r39.g10_female_hypothalamus_bldg x.r39.g10_female_pituitary_bldg
    ## 1                           111.0000                         196.000
    ## 2                             1.0000                          71.000
    ## 3                           428.0000                         606.000
    ## 4                             0.0000                          11.000
    ## 5                             0.0000                           2.000
    ## 6                            22.3529                         195.592
    ##   x.r44.w95_female_gonad_hatch x.r44.w95_female_hypothalamus_hatch
    ## 1                      74.1788                            246.0000
    ## 2                      92.0000                              2.0000
    ## 3                     338.0000                            931.0000
    ## 4                      12.0000                              0.0000
    ## 5                       6.0000                              5.0000
    ## 6                      34.0000                             58.0675
    ##   x.r44.w95_female_pituitary_hatch x.r48.y139_female_gonad_inc.d17
    ## 1                               35                              23
    ## 2                                4                              68
    ## 3                              285                             166
    ## 4                                1                              36
    ## 5                                0                               2
    ## 6                               67                              13
    ##   x.r48.y139_female_hypothalamus_inc.d17
    ## 1                               271.0000
    ## 2                                 1.0000
    ## 3                               633.0000
    ## 4                                 0.0000
    ## 5                                 0.0000
    ## 6                                49.1383
    ##   x.r48.y139_female_pituitary_inc.d17.NYNO x.r50.w97_female_gonad_n5
    ## 1                                  128.000                       101
    ## 2                                   29.000                       189
    ## 3                                  822.000                       969
    ## 4                                    0.000                         5
    ## 5                                    0.000                        10
    ## 6                                  166.599                       362
    ##   x.r50.w97_female_hypothalamus_n5 x.r50.w97_female_pituitary_n5
    ## 1                         464.0000                           358
    ## 2                           7.0000                            60
    ## 3                         972.0000                           793
    ## 4                           0.0000                             0
    ## 5                           1.0000                             0
    ## 6                          49.4283                           183
    ##   x.r64.g140_male_gonad_inc.d3 x.r64.g140_male_hypothalamus_inc.d3
    ## 1                       21.000                           334.00000
    ## 2                        6.000                             8.00000
    ## 3                     1120.000                           676.00000
    ## 4                        4.000                             0.00000
    ## 5                        1.000                             8.86496
    ## 6                      336.924                            55.36700
    ##   x.r64.g140_male_pituitary_inc.d3 x.r67.blu35_male_gonad_bldg
    ## 1                          101.000                           1
    ## 2                           59.000                           3
    ## 3                          362.000                         460
    ## 4                            0.000                           2
    ## 5                            3.000                           4
    ## 6                          116.705                         183
    ##   x.r67.blu35_male_hypothalamus_bldg x.r67.blu35_male_pituitary_bldg.NYNO
    ## 1                             96.000                            195.00000
    ## 2                              3.000                             45.00000
    ## 3                            616.000                            465.00000
    ## 4                              0.000                              0.00000
    ## 5                              0.000                              4.10243
    ## 6                             53.623                            147.86600
    ##   x.r74.o29_female_gonad_m.hatch x.r74.o29_female_hypothalamus_m.hatch
    ## 1                             75                                332.00
    ## 2                             82                                  7.00
    ## 3                            607                                566.00
    ## 4                             14                                  0.00
    ## 5                              4                                  1.00
    ## 6                            117                                 48.23
    ##   x.r74.o29_female_pituitary_m.hatch x.r76_female_gonad_m.inc.d9
    ## 1                            444.553                          63
    ## 2                             68.000                         142
    ## 3                           1003.000                         344
    ## 4                              0.000                          10
    ## 5                              0.000                           4
    ## 6                            331.455                          32
    ##   x.r76_female_hypothalamus_m.inc.d9 x.r76_female_pituitary_m.inc.d9
    ## 1                           239.0000                         53.0000
    ## 2                             0.0000                         13.0000
    ## 3                           610.0000                        178.0000
    ## 4                             0.0000                          0.0000
    ## 5                             0.0000                          2.0000
    ## 6                            26.1712                         18.2918
    ##   x.s188_female_gonad_m.inc.d8 x.s188_female_hypothalamus_m.inc.d8
    ## 1                     241.0000                                 244
    ## 2                     324.0000                                   5
    ## 3                     705.0000                                 557
    ## 4                      35.0000                                   0
    ## 5                       3.0000                                   1
    ## 6                      87.9361                                  38
    ##   x.s188_female_pituitary_m.inc.d8 x.w178_female_gonad_n9
    ## 1                        374.00000                60.0000
    ## 2                        136.00000                64.0000
    ## 3                        966.00000               344.0000
    ## 4                          0.00000                 0.0000
    ## 5                          4.03775                 2.0000
    ## 6                        352.99000                64.9316
    ##   x.w178_female_hypothalamus_n9 x.w178_female_pituitary_n9
    ## 1                      336.0000                        100
    ## 2                        4.0000                         30
    ## 3                      915.0000                        275
    ## 4                        0.0000                          0
    ## 5                        0.0000                          0
    ## 6                       54.5463                         82
    ##   x.w192.o157_male_gonad_inc.d9 x.w192.o157_male_hypothalamus_inc.d9
    ## 1                       5.00000                             388.0000
    ## 2                       1.00000                               7.0000
    ## 3                     498.00000                             658.0000
    ## 4                       4.00000                              37.0000
    ## 5                       6.88128                               1.0000
    ## 6                     172.00000                              65.3199
    ##   x.w192.o157_male_pituitary_inc.d9 x.w51_female_gonad_lay
    ## 1                           169.000                    175
    ## 2                           141.000                    177
    ## 3                           554.000                    730
    ## 4                             0.000                     47
    ## 5                             1.000                      5
    ## 6                           226.246                     81
    ##   x.w51_female_hypothalamus_lay x.w51_female_pituitary_lay
    ## 1                      433.3630                    104.000
    ## 2                       18.0000                     16.000
    ## 3                      541.0000                    530.000
    ## 4                        0.0000                      0.000
    ## 5                        0.0000                      1.000
    ## 6                       44.2758                    228.802
    ##   x.w6_female_gonad_n9 x.w6_female_hypothalamus_n9
    ## 1                  165                    153.0000
    ## 2                   35                      0.0000
    ## 3                  798                    449.0000
    ## 4                   21                      0.0000
    ## 5                    8                      1.0000
    ## 6                  188                     20.0494
    ##   x.w6_female_pituitary_n9 x.y.s.ATLAS_male_gonad_control
    ## 1                  160.000                         13.000
    ## 2                   36.000                          0.000
    ## 3                  476.000                        772.000
    ## 4                    0.000                         15.000
    ## 5                    2.000                          3.000
    ## 6                  159.063                        135.693
    ##   x.y.s.ATLAS_male_pituitary_control x.y109_female_gonad_inc.d9
    ## 1                                 27                         27
    ## 2                                 55                         64
    ## 3                                139                        310
    ## 4                                  0                         14
    ## 5                                  0                          2
    ## 6                                 10                         28
    ##   x.y109_female_hypothalamus_inc.d9 x.y109_female_pituitary_inc.d9
    ## 1                          115.0000                             89
    ## 2                            0.0000                             13
    ## 3                          385.0000                            216
    ## 4                            0.0000                              0
    ## 5                            0.0000                              0
    ## 6                           21.1942                             35
    ##   x.y119.w11_female_gonad_m.inc.d17
    ## 1                                55
    ## 2                               108
    ## 3                               349
    ## 4                                 5
    ## 5                                 2
    ## 6                                75
    ##   x.y119.w11_female_hypothalamus_m.inc.d17
    ## 1                                 106.0000
    ## 2                                   2.0000
    ## 3                                 255.0000
    ## 4                                   0.0000
    ## 5                                   0.0000
    ## 6                                  14.0527
    ##   x.y119.w11_female_pituitary_m.inc.d17 x.y132.w76_male_gonad_inc.d17
    ## 1                               223.000                            12
    ## 2                               118.000                             2
    ## 3                               481.000                           661
    ## 4                                 4.000                             0
    ## 5                                 5.000                             2
    ## 6                               137.661                           234
    ##   x.y132.w76_male_hypothalamus_inc.d17 x.y132.w76_male_pituitary_inc.d17
    ## 1                             213.0000                                83
    ## 2                               3.0000                                53
    ## 3                             470.0000                               325
    ## 4                               0.0000                                 1
    ## 5                               0.0000                                 0
    ## 6                              50.1584                                86
    ##   x.y138.w176_female_gonad_n9 x.y138.w176_female_hypothalamus_n9
    ## 1                          12                            261.000
    ## 2                          50                              6.000
    ## 3                         426                            573.000
    ## 4                           0                              0.000
    ## 5                           0                              1.000
    ## 6                          90                             48.274
    ##   x.y138.w176_female_pituitary_n9 x.y141.w116_male_gonad_inc.d9
    ## 1                       92.000000                        20.000
    ## 2                       38.000000                         0.000
    ## 3                      437.000000                      1705.000
    ## 4                        0.000000                         7.000
    ## 5                        2.976245                         3.000
    ## 6                       69.697900                       434.976
    ##   x.y141.w116_male_hypothalamus_inc.d9 x.y141.w116_male_pituitary_inc.d9
    ## 1                             166.0000                           358.000
    ## 2                               3.0000                            71.000
    ## 3                             372.0000                           707.000
    ## 4                               0.0000                             0.000
    ## 5                               0.0000                             1.000
    ## 6                              25.3426                           166.149
    ##   x.y145.r55_female_gonad_inc.prolong
    ## 1                                  50
    ## 2                                  62
    ## 3                                 281
    ## 4                                   3
    ## 5                                   1
    ## 6                                  41
    ##   x.y145.r55_female_hypothalamus_inc.prolong
    ## 1                                   300.0000
    ## 2                                     0.0000
    ## 3                                   739.0000
    ## 4                                     0.0000
    ## 5                                     1.0000
    ## 6                                    65.8116
    ##   x.y145.r55_female_pituitary_inc.prolong x.y51_male_gonad_m.inc.d8
    ## 1                                      75                         8
    ## 2                                      20                         3
    ## 3                                     242                       607
    ## 4                                       0                         0
    ## 5                                       1                         2
    ## 6                                      51                       435
    ##   x.y51_male_hypothalamus_m.inc.d8 x.y51_male_pituitary_m.inc.d8
    ## 1                         257.0000                      422.0000
    ## 2                           2.0000                       85.0000
    ## 3                         830.0000                      464.0000
    ## 4                           0.0000                        0.0000
    ## 5                           2.0000                        2.8266
    ## 6                          88.3668                      262.1060
    ##   x.y5_female_gonad_m.inc.d8 x.y5_female_hypothalamus_m.inc.d8
    ## 1                        101                          249.0000
    ## 2                        313                           14.0000
    ## 3                       1243                          479.0000
    ## 4                         12                            0.0000
    ## 5                          2                            0.0000
    ## 6                         67                           23.2768
    ##   x.y5_female_pituitary_m.inc.d8 x.y64_female_gonad_extend
    ## 1                       309.3845                       138
    ## 2                        29.0000                       255
    ## 3                       630.0000                       642
    ## 4                         0.0000                        24
    ## 5                         5.0000                         0
    ## 6                       135.6580                        93
    ##   x.y64_female_hypothalamus_extend x.y64_female_pituitary_extend
    ## 1                         336.0000                         249.0
    ## 2                           3.0000                          13.0
    ## 3                         557.0000                         572.0
    ## 4                           0.0000                           0.0
    ## 5                           1.0000                           2.0
    ## 6                          48.2873                         206.6
    ##   x.y79_male_gonad_prolong x.y79_male_hypothalamus_prolong
    ## 1                   32.000                        385.0000
    ## 2                    1.000                          3.0000
    ## 3                  378.000                        649.0000
    ## 4                    0.000                          0.0000
    ## 5                    5.000                          0.0000
    ## 6                  377.955                         64.2577
    ##   x.y79_male_pituitary_prolong x.y90_female_gonad_hatch
    ## 1                     261.3429                       36
    ## 2                     122.0000                      106
    ## 3                     579.0000                      338
    ## 4                       1.0000                        7
    ## 5                       2.0000                        0
    ## 6                     246.1290                       20
    ##   x.y90_female_hypothalamus_hatch x.y90_female_pituitary_hatch
    ## 1                              87                      144.000
    ## 2                               0                       12.000
    ## 3                             197                      678.000
    ## 4                               0                        0.000
    ## 5                               0                        1.000
    ## 6                               5                      244.345
    ##   x.y93.g126_female_gonad_inc.d9 x.y93.g126_female_hypothalamus_inc.d9
    ## 1                        33.0000                                   122
    ## 2                       117.0000                                     0
    ## 3                       292.0000                                   407
    ## 4                        18.0000                                     0
    ## 5                         0.0000                                     1
    ## 6                        39.8289                                    37
    ##   x.y93.g126_female_pituitary_inc.d9 x.y9_female_gonad_n9
    ## 1                            260.000                   85
    ## 2                             57.000                   55
    ## 3                            627.000                  602
    ## 4                              0.000                   14
    ## 5                              2.000                    0
    ## 6                            220.045                  332
    ##   x.y9_female_hypothalamus_n9 x.y9_female_pituitary_n9
    ## 1                    179.0000                 299.4387
    ## 2                      0.0000                  26.0000
    ## 3                    386.0000                 414.0000
    ## 4                      0.0000                   0.0000
    ## 5                      3.0000                   4.0000
    ## 6                     27.4368                 166.2900
    ##   y.s156.o.r_female_gonad_lay y.s156.o.r_female_hypothalamus_lay
    ## 1                   586.00000                           848.6660
    ## 2                    71.00000                             5.0000
    ## 3                   427.00000                          1314.0000
    ## 4                    64.00000                             0.0000
    ## 5                    13.03804                             3.0000
    ## 6                    50.00000                            66.7481
    ##   y.s156.o.r_female_pituitary_lay y.w.blk_male_gonad_m.inc.d9
    ## 1                         419.583                     5.03668
    ## 2                          30.000                     3.00000
    ## 3                         772.000                   610.00000
    ## 4                           0.000                     0.00000
    ## 5                           1.000                     5.00000
    ## 6                         201.084                   175.00000
    ##   y.w.blk_male_hypothalamus_m.inc.d9 y.w.blk_male_pituitary_m.inc.d9
    ## 1                           293.0000                              76
    ## 2                             2.0000                              13
    ## 3                           638.0000                             210
    ## 4                             0.0000                               0
    ## 5                             1.0000                               0
    ## 6                            71.7467                              50
    ##   y12.x_female_gonad_m.inc.d9 y12.x_female_hypothalamus_m.inc.d9
    ## 1                          33                           201.0000
    ## 2                          43                             1.0000
    ## 3                         225                           289.0000
    ## 4                          21                             0.0000
    ## 5                           2                             0.0000
    ## 6                          25                            23.2119
    ##   y12.x_female_pituitary_m.inc.d9 y126.w92.x_female_gonad_inc.d17
    ## 1                         70.0000                              30
    ## 2                          5.0000                              80
    ## 3                        281.0000                             289
    ## 4                          0.0000                              23
    ## 5                          1.0000                               1
    ## 6                         49.0215                              37
    ##   y126.w92.x_female_hypothalamus_inc.d17
    ## 1                               445.0000
    ## 2                                 7.0000
    ## 3                               456.0000
    ## 4                                 0.0000
    ## 5                                 3.0000
    ## 6                                46.4953
    ##   y126.w92.x_female_pituitary_inc.d17 y128.g23.x_female_gonad_inc.d9
    ## 1                             72.0000                        41.0000
    ## 2                              3.0000                         5.0000
    ## 3                            173.0000                       303.0000
    ## 4                              0.0000                         7.0000
    ## 5                              0.0000                         1.0000
    ## 6                             66.9048                        79.9014
    ##   y128.g23.x_female_hypothalamus_m.inc.d9
    ## 1                                283.0000
    ## 2                                  6.0000
    ## 3                                639.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 40.4052
    ##   y128.g23.x_female_pituitary_inc.d9 y129.x_male_gonad_n9
    ## 1                            70.2786                7.000
    ## 2                            22.0000                1.000
    ## 3                           236.0000              686.000
    ## 4                             0.0000                2.000
    ## 5                             2.0000                2.000
    ## 6                           150.5000              246.735
    ##   y129.x_male_hypothalamus_n9 y129.x_male_pituitary_n9
    ## 1                     145.000                  77.0000
    ## 2                       0.000                  59.0000
    ## 3                     519.000                 143.0000
    ## 4                       0.000                   0.0000
    ## 5                       1.000                   0.0000
    ## 6                      33.115                  50.6307
    ##   y13.x_female_gonad_inc.d3 y13.x_female_hypothalamus_inc.d3
    ## 1                        63                         183.0000
    ## 2                        23                           0.0000
    ## 3                       280                         430.0000
    ## 4                         2                           0.0000
    ## 5                         0                           0.0000
    ## 6                        19                          38.7971
    ##   y13.x_female_pituitary_inc.d3 y130.o170.x_female_gonad_inc.d17
    ## 1                            86                               44
    ## 2                            27                              119
    ## 3                           254                              171
    ## 4                             0                               10
    ## 5                             0                                0
    ## 6                            42                               33
    ##   y130.o170.x_female_hypothalamus_inc.d17
    ## 1                                287.0000
    ## 2                                  0.0000
    ## 3                                577.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 65.3311
    ##   y130.o170.x_female_pituitary_inc.d17 y131.w185.x_male_gonad_n9
    ## 1                                   91                         3
    ## 2                                   15                         1
    ## 3                                  206                       504
    ## 4                                    0                         1
    ## 5                                    1                         4
    ## 6                                   73                       130
    ##   y131.w185.x_male_hypothalamus_n9 y131.w185.x_male_pituitary_n9
    ## 1                         294.0000                       204.000
    ## 2                           4.0000                       183.000
    ## 3                         704.0000                       662.000
    ## 4                           0.0000                         0.000
    ## 5                           1.0000                         0.000
    ## 6                          39.1707                       152.338
    ##   y133.w77.r58_male_gonad_inc.d17 y133.w77.r58_male_hypothalamus_inc.d17
    ## 1                              12                               296.0000
    ## 2                               0                                 8.0000
    ## 3                             514                               582.0000
    ## 4                               0                                 0.0000
    ## 5                               0                                 2.7485
    ## 6                             137                                50.7092
    ##   y133.w77.r58_male_pituitary_inc.d17 y135.blu107.x_female_gonad_inc.d17
    ## 1                           102.00000                                 35
    ## 2                            21.00000                                107
    ## 3                           206.00000                                291
    ## 4                             0.00000                                 23
    ## 5                             1.08287                                  1
    ## 6                            75.00000                                 32
    ##   y135.blu107.x_female_hypothalamus_inc.d17
    ## 1                                  388.0000
    ## 2                                    5.0000
    ## 3                                  646.0000
    ## 4                                    0.0000
    ## 5                                    2.0000
    ## 6                                   55.3267
    ##   y135.blu107.x_female_pituitary_inc.d17.NYNO y136.x_female_gonad_inc.d17
    ## 1                                      209.00                          33
    ## 2                                        9.00                          65
    ## 3                                      444.00                         192
    ## 4                                        0.00                          57
    ## 5                                        0.00                           1
    ## 6                                      163.22                          20
    ##   y136.x_female_hypothalamus_inc.d17 y136.x_female_pituitary_inc.d17
    ## 1                           253.0000                        106.0000
    ## 2                             1.0000                          7.0000
    ## 3                           595.0000                        271.0000
    ## 4                             0.0000                          0.0000
    ## 5                             0.0000                          0.0000
    ## 6                            38.0986                         44.2683
    ##   y140.w119.x_female_gonad_inc.d9 y140.w119.x_female_hypothalamus_inc.d9
    ## 1                              66                                    228
    ## 2                             106                                      1
    ## 3                             296                                    553
    ## 4                              77                                      0
    ## 5                               6                                      0
    ## 6                              31                                     31
    ##   y140.w119.x_female_pituitary_inc.d9 y146.r32.x_female_gonad_inc.prolong
    ## 1                              203.00                                  44
    ## 2                               77.00                                  71
    ## 3                              679.00                                 335
    ## 4                                0.00                                   6
    ## 5                                0.00                                   2
    ## 6                              111.51                                  59
    ##   y146.r32.x_female_hypothalamus_inc.prolong
    ## 1                                   243.0000
    ## 2                                     2.0000
    ## 3                                   717.0000
    ## 4                                     0.0000
    ## 5                                     1.0000
    ## 6                                    85.9667
    ##   y146.r32.x_female_pituitary_inc.prolong
    ## 1                                 65.0000
    ## 2                                 31.0000
    ## 3                                250.0000
    ## 4                                  0.0000
    ## 5                                  2.0000
    ## 6                                 94.4045
    ##   y148.w81.g145_male_gonad_m.inc.d17
    ## 1                            58.0000
    ## 2                             2.0000
    ## 3                          1702.0000
    ## 4                             7.0000
    ## 5                            14.8705
    ## 6                           795.8680
    ##   y148.w81.g145_male_hypothalamus_m.inc.d17
    ## 1                                   69.0000
    ## 2                                    2.0000
    ## 3                                  230.0000
    ## 4                                    0.0000
    ## 5                                    0.0000
    ## 6                                   16.0543
    ##   y148.w81.g145_male_pituitary_m.inc.d17 y149.r52.x_male_gonad_inc.d3
    ## 1                                547.000                            8
    ## 2                                153.000                            4
    ## 3                               1001.000                          444
    ## 4                                  0.000                            0
    ## 5                                  1.000                            2
    ## 6                                389.522                          137
    ##   y149.r52.x_male_hypothalamus_inc.d3 y149.r52.x_male_pituitary_inc.d3
    ## 1                             407.000                              135
    ## 2                               3.000                               11
    ## 3                             797.000                              222
    ## 4                               0.000                                0
    ## 5                               1.000                                3
    ## 6                              74.137                               55
    ##   y15.x_female_gonad_hatch y15.x_female_hypothalamus_hatch
    ## 1                       63                             307
    ## 2                       93                               7
    ## 3                      353                             581
    ## 4                        3                               0
    ## 5                        4                               3
    ## 6                       42                              68
    ##   y15.x_female_pituitary_hatch y18.x_male_gonad_m.inc.d3
    ## 1                      77.2156                   44.0000
    ## 2                      17.0000                    7.0000
    ## 3                     213.0000                 1732.0000
    ## 4                       0.0000                    0.0000
    ## 5                       3.0000                   15.8283
    ## 6                      59.7258                  897.8130
    ##   y18.x_male_hypothalamus_m.inc.d3 y18.x_male_pituitary_m.inc.d3
    ## 1                         128.0000                           480
    ## 2                           1.0000                            66
    ## 3                         357.0000                           725
    ## 4                           0.0000                             0
    ## 5                           0.0000                             5
    ## 6                          19.3011                           153
    ##   y4.x_male_gonad_m.inc.d17 y4.x_male_hypothalamus_m.inc.d17
    ## 1                         3                         127.0000
    ## 2                         0                           1.0000
    ## 3                       514                         253.0000
    ## 4                         0                           0.0000
    ## 5                         1                           0.0000
    ## 6                       134                          19.1706
    ##   y4.x_male_pituitary_m.inc.d17 y55.x_male_gonad_m.inc.d8
    ## 1                      100.0000                    46.000
    ## 2                       14.0000                     0.000
    ## 3                      191.0000                  1756.000
    ## 4                        0.0000                     0.000
    ## 5                        0.0000                    20.000
    ## 6                       67.7619                   496.855
    ##   y55.x_male_hypothalamus_m.inc.d8 y55.x_male_pituitary_m.inc.d8
    ## 1                         228.2791                       182.250
    ## 2                           3.0000                        90.000
    ## 3                         733.0000                       487.000
    ## 4                           0.0000                         0.000
    ## 5                           1.0000                         1.000
    ## 6                          46.1482                       142.488
    ##   y6.o54_female_gonad_n5 y6.o54_female_hypothalamus_n5
    ## 1                    198                      214.0000
    ## 2                    277                        5.0000
    ## 3                    674                      635.0000
    ## 4                     20                        0.0000
    ## 5                     10                        0.0000
    ## 6                    130                       56.2059
    ##   y6.o54_female_pituitary_n5 y63.x_male_gonad_m.inc.d9
    ## 1                    272.000                         7
    ## 2                     39.000                         3
    ## 3                    599.000                       497
    ## 4                      0.000                         0
    ## 5                      2.000                         1
    ## 6                    310.354                       158
    ##   y63.x_male_hypothalamus_m.inc.d9 y63.x_male_pituitary_m.inc.d9
    ## 1                          67.6528                      124.0000
    ## 2                           1.0000                       27.0000
    ## 3                         294.0000                      361.0000
    ## 4                           0.0000                        2.0000
    ## 5                           0.0000                        0.0000
    ## 6                           9.0000                       68.5814
    ##   y7.g58_female_gonad_hatch y7.g58_female_hypothalamus_hatch
    ## 1                        32                           241.00
    ## 2                        68                             2.00
    ## 3                       178                           516.00
    ## 4                         3                             0.00
    ## 5                         1                             3.00
    ## 6                        25                            49.09
    ##   y7.g58_female_pituitary_hatch y85.r71.x_female_gonad_m.inc.d17
    ## 1                       187.000                               44
    ## 2                        34.000                               26
    ## 3                       511.000                              332
    ## 4                         0.000                                9
    ## 5                         3.000                                2
    ## 6                       130.283                               46
    ##   y85.r71.x_female_hypothalamus_m.inc.d17
    ## 1                                313.0000
    ## 2                                  4.0000
    ## 3                                607.0000
    ## 4                                  0.0000
    ## 5                                  0.0000
    ## 6                                 62.5207
    ##   y85.r71.x_female_pituitary_m.inc.d17 y94.g133.x_female_gonad_n5
    ## 1                                  122                         72
    ## 2                                   28                         61
    ## 3                                  273                        178
    ## 4                                    0                         13
    ## 5                                    1                          1
    ## 6                                   82                         27
    ##   y94.g133.x_female_hypothalamus_n5.NYNO y94.g133.x_female_pituitary_n5
    ## 1                               175.0000                        300.000
    ## 2                                 3.0000                         44.000
    ## 3                               528.0000                        612.000
    ## 4                                 0.0000                          0.000
    ## 5                                 1.0000                          2.000
    ## 6                                31.4769                        142.049
    ##   y95.g131.x_male_gonad_inc.d9 y95.g131.x_male_hypothalamus_inc.d9
    ## 1                        5.000                            743.0000
    ## 2                        1.000                             10.0000
    ## 3                      560.000                           1540.0000
    ## 4                        0.000                              0.0000
    ## 5                        6.000                              1.0000
    ## 6                      391.944                             87.5236
    ##   y95.g131.x_male_pituitary_inc.d9 y97.x_female_gonad_n9
    ## 1                               91             455.09400
    ## 2                                7             422.00000
    ## 3                              145             717.00000
    ## 4                                0               8.00000
    ## 5                                0              13.64535
    ## 6                               42             185.00000
    ##   y97.x_female_hypothalamus_n9 y97.x_female_pituitary_n9
    ## 1                     88.00000                   187.000
    ## 2                      0.00000                    75.000
    ## 3                    242.00000                   536.000
    ## 4                      0.00000                     0.000
    ## 5                      0.00000                     3.000
    ## 6                      6.04743                   129.126
    ##   y98.g54_female_gonad_m.hatch y98.g54_female_hypothalamus_m.hatch
    ## 1                           26                            241.0000
    ## 2                           35                             11.0000
    ## 3                          380                            564.0000
    ## 4                            1                              0.0000
    ## 5                           12                              1.0000
    ## 6                          109                             43.7713
    ##   y98.g54_female_pituitary_m.hatch y98.o50.x_male_gonad_inc.d3
    ## 1                          270.000                           7
    ## 2                          222.000                           1
    ## 3                          378.000                         462
    ## 4                            0.000                           0
    ## 5                            0.000                           3
    ## 6                          111.913                         174
    ##   y98.o50.x_male_hypothalamus_inc.d3 y98.o50.x_male_pituitary_inc.d3
    ## 1                           245.0000                        130.0000
    ## 2                             3.0000                         11.0000
    ## 3                           567.0000                        223.0000
    ## 4                             0.0000                          0.0000
    ## 5                             2.0000                          2.0000
    ## 6                            46.1529                         67.6398

    # aggregate transcript counts to gene counts
    countData <- kallistodata %>% 
      select(-row.names, -geneid, -NCBI) %>% 
      pivot_longer(-gene, names_to = "samples", values_to = "counts") %>%  
      pivot_wider(
        names_from = samples, 
        values_from = counts,
        values_fn = list(counts = sum))  %>% 
      arrange(gene) %>% 
      filter(gene != "")
    countData <- as.data.frame(countData)
    row.names(countData) <- countData$gene ## make gene the row name
    countData[1] <- NULL ## make gene the row name
    countData <- round(countData) #round all value to nearest 1s place
    head(countData[13:15])

    ##        L.R8_male_gonad_control L.R8_male_hypothalamus_control
    ## A2ML1                       79                              7
    ## A2ML2                        2                              0
    ## A2ML3                      233                            211
    ## A2ML4                       12                              0
    ## A4GALT                       9                              0
    ## A4GNT                        0                              0
    ##        L.R8_male_pituitary_control
    ## A2ML1                           38
    ## A2ML2                            0
    ## A2ML3                           73
    ## A2ML4                            7
    ## A4GALT                          63
    ## A4GNT                            0

    # print tolal num of genes and samples
    dim(countData)

    ## [1] 13966   987

wrangle colData
---------------

    colData <- read.table(file.path( "../metadata/kallistosamples.txt"),
                          header = F, stringsAsFactors = F) %>%
      # use strsplit to cut the filename into meaningful columns
      mutate(bird = sapply(strsplit(V1,'\\_'), "[", 1),
             sex = sapply(strsplit(V1,'\\_'), "[", 2),
             tissue = sapply(strsplit(V1,'\\_'), "[", 3),
             temp = sapply(strsplit(V1,'\\_'), "[", 4)) %>%
      mutate(treatmenttemp = sapply(strsplit(temp,'\\.'), "[", 1),
             NYNO = sapply(strsplit(temp,'\\.'), "[", 2)) %>%
     
      # rename variables
      mutate(treatment = ifelse(grepl("extend-hatch", treatmenttemp), "extend",
                                ifelse(grepl("inc-prolong", treatmenttemp), "prolong",
                                       ifelse(grepl("m.hatch", treatmenttemp), "m.n2",
                                              ifelse(grepl("m.inc.d8", treatmenttemp), "early",
                                       treatmenttemp))))) %>%
       select(-temp, -NYNO, -treatmenttemp ) %>%
      # replace dashes with periods (to make deseq happy)
      mutate(bird = gsub("-", ".", bird),
             treatment = gsub("-", ".", treatment),
             V1 = gsub("-", ".", V1))
    # specify levels that will be factors
    cols <- c("sex", "treatment", "tissue")
    colData[cols] <- lapply(colData[cols], factor) 
    colData <- colData %>%
      mutate(tissue = factor(tissue, levels = c("hypothalamus", "pituitary", "gonad"))) %>% 
      mutate(tissue = fct_recode(tissue, "gonads" = "gonad")) %>%
      mutate(group = paste(sex, tissue, treatment, sep = "."),
             study = ifelse(grepl("m.|extend|prolong", treatment), 
                            "manipulation", "charcterization")) %>%
      mutate(study = factor(study, levels = c("charcterization", "manipulation")))
    colData <- as.data.frame(colData)
    row.names(colData) <- colData$V1
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

check that rownames and colnames match for DESeq
------------------------------------------------

    str(colData)

    ## 'data.frame':    987 obs. of  7 variables:
    ##  $ V1       : chr  "L.Blu13_male_gonad_control.NYNO" "L.Blu13_male_hypothalamus_control.NYNO" "L.Blu13_male_pituitary_control.NYNO" "L.G107_male_gonad_control" ...
    ##  $ bird     : chr  "L.Blu13" "L.Blu13" "L.Blu13" "L.G107" ...
    ##  $ sex      : Factor w/ 2 levels "female","male": 2 2 2 2 2 2 1 1 1 2 ...
    ##  $ tissue   : Factor w/ 3 levels "hypothalamus",..: 3 1 2 3 1 2 3 1 2 3 ...
    ##  $ treatment: Factor w/ 16 levels "bldg","control",..: 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ group    : chr  "male.gonads.control" "male.hypothalamus.control" "male.pituitary.control" "male.gonads.control" ...
    ##  $ study    : Factor w/ 2 levels "charcterization",..: 1 1 1 1 1 1 1 1 1 1 ...

    ncol(countData) == nrow(colData)

    ## [1] TRUE

    colData %>% select(sex, tissue, treatment, study)  %>%  summary()

    ##      sex               tissue        treatment               study    
    ##  female:497   hypothalamus:327   control  : 73   charcterization:636  
    ##  male  :490   pituitary   :330   inc.d9   : 71   manipulation   :351  
    ##               gonads      :330   inc.d17  : 66                        
    ##                                  n9       : 66                        
    ##                                  m.inc.d17: 63                        
    ##                                  bldg     : 60                        
    ##                                  (Other)  :588

    table(colData$sex, colData$treatment, colData$tissue)

    ## , ,  = hypothalamus
    ## 
    ##         
    ##          bldg control early extend hatch inc.d17 inc.d3 inc.d9 lay
    ##   female   10      11    10     10    10      11     10     12  10
    ##   male     10      11    10     10    10      11     10     11  10
    ##         
    ##          m.inc.d17 m.inc.d3 m.inc.d9 m.n2 n5 n9 prolong
    ##   female        11       10        9   10 10 11      10
    ##   male          10       10        8   10 10 11      10
    ## 
    ## , ,  = pituitary
    ## 
    ##         
    ##          bldg control early extend hatch inc.d17 inc.d3 inc.d9 lay
    ##   female   10      11    10     10    10      11     10     13  10
    ##   male     10      14    10     10    10      11     10     11  10
    ##         
    ##          m.inc.d17 m.inc.d3 m.inc.d9 m.n2 n5 n9 prolong
    ##   female        11       10        8   10 10 11      10
    ##   male          10       10        8   10 10 11      10
    ## 
    ## , ,  = gonads
    ## 
    ##         
    ##          bldg control early extend hatch inc.d17 inc.d3 inc.d9 lay
    ##   female   10      13    10     10    10      11     10     13  10
    ##   male     10      13    10     10    10      11     10     11  10
    ##         
    ##          m.inc.d17 m.inc.d3 m.inc.d9 m.n2 n5 n9 prolong
    ##   female        11       10        8   10 10 11      10
    ##   male          10       10        8    9 10 11      10

save bird data
--------------

    birds <- colData %>%
      select(bird, sex, treatment) %>%
      distinct()
    birds

    ##                 bird    sex treatment
    ## 1            L.Blu13   male   control
    ## 2             L.G107   male   control
    ## 3             L.G118 female   control
    ## 4               L.R3   male   control
    ## 5               L.R8   male   control
    ## 6              L.W33   male   control
    ## 7               L.W3   male   control
    ## 8               L.W4   male   control
    ## 9             R.G106 female   control
    ## 10             R.R20 female   control
    ## 11              R.R9 female   control
    ## 12             R.W44 female   control
    ## 13        R.Y108.W29   male   control
    ## 14      blk.s030.o.g   male   prolong
    ## 15     blk.s031.pu.d female   prolong
    ## 16      blk.s032.g.w female      m.n2
    ## 17      blk.s049.y.g female  m.inc.d3
    ## 18     blk.s060.pu.w female  m.inc.d3
    ## 19     blk.s061.pu.y female    inc.d9
    ## 20      blk.y.l.s109 female     early
    ## 21            blk0.x female      m.n2
    ## 22           blk11.x female      bldg
    ## 23           blk12.x   male        n5
    ## 24           blk17.x   male   inc.d17
    ## 25           blk19.x female    extend
    ## 26           blk21.x female     hatch
    ## 27            blk4.x female        n9
    ## 28            blk5.x   male  m.inc.d3
    ## 29     blu.o.x.ATLAS female   control
    ## 30       blu10.w26.x   male      m.n2
    ## 31          blu103.x female     hatch
    ## 32     blu104.w120.x   male     hatch
    ## 33   blu108.w40.o158   male    inc.d9
    ## 34     blu111.w113.x   male    inc.d3
    ## 35     blu113.w124.x   male   inc.d17
    ## 36   blu114.r38.w198   male      bldg
    ## 37     blu115.y150.x female   prolong
    ## 38      blu119.w84.x female     early
    ## 39      blu121.w91.x   male   inc.d17
    ## 40     blu124.w180.x female     hatch
    ## 41       blu33.y88.x   male      bldg
    ## 42         blu36.w16 female        n9
    ## 43       blu37.r65.x   male        n5
    ## 44      blu38.g135.x female      bldg
    ## 45       blu39.o26.x female    inc.d3
    ## 46      blu41.y100.x   male        n5
    ## 47        blu44.y102 female    extend
    ## 48       blu47.y96.x female    inc.d9
    ## 49         blu55.g51 female        n5
    ## 50         blu56.o53 female  m.inc.d3
    ## 51         blu63.g62 female  m.inc.d9
    ## 52         blu80.r97 female     early
    ## 53         blu81.r88   male        n9
    ## 54           blu84.x   male    extend
    ## 55      d.r.blk.s159 female  m.inc.d9
    ## 56      d.s008.y.blk   male        n5
    ## 57      d.s047.blk.o   male        n5
    ## 58      d.s110.g.blk   male  m.inc.d3
    ## 59      d.s112.blk.w female m.inc.d17
    ## 60      d.s177.blk.r female  m.inc.d3
    ## 61     g.blk.s004.pk female       lay
    ## 62      g.blk.s041.r   male  m.inc.d3
    ## 63        g.o.y.s037   male m.inc.d17
    ## 64         g.s.blk.d   male        n9
    ## 65         g.s.blk.y   male       lay
    ## 66     g.s043.pu.blk   male       lay
    ## 67      g.s075.pk.pu   male      m.n2
    ## 68       g.s076.pk.r female      m.n2
    ## 69      g.s078.blk.o female       lay
    ## 70      g.s111.r.blk   male     early
    ## 71       g.s179.o.pk   male     early
    ## 72       g.s351.pk.w   male    extend
    ## 73         g.x.ATLAS female   control
    ## 74      g.y.blk.s006 female m.inc.d17
    ## 75           g.y.o.s   male   prolong
    ## 76        g104.w82.x   male      bldg
    ## 77        g114.w83.x   male     hatch
    ## 78        g130.y81.x   male   inc.d17
    ## 79       g137.r24.w5   male     early
    ## 80      g141.blu27.x female      bldg
    ## 81        g142.r40.x female   inc.d17
    ## 82      g143.blu32.x   male   inc.d17
    ## 83        g144.r54.x female  m.inc.d3
    ## 84        g146.blu51   male    inc.d3
    ## 85        g17.w108.x female    extend
    ## 86        g20.w106.x   male    inc.d3
    ## 87        g22.blu118 female    extend
    ## 88       g3.g119.w20   male    extend
    ## 89         g32.blu79   male m.inc.d17
    ## 90             g34.x   male      m.n2
    ## 91             g38.x   male   prolong
    ## 92         g52.blu58   male      bldg
    ## 93           g53.y84   male     hatch
    ## 94         g6.w197.x female    inc.d3
    ## 95         g63.blu65 female m.inc.d17
    ## 96             g73.x female  m.inc.d9
    ## 97             g75.x female    inc.d9
    ## 98           g8.y197   male    extend
    ## 99         l.s.o.blk   male    extend
    ## 100          l.s.w.d female      m.n2
    ## 101       l.s024.y.g   male m.inc.d17
    ## 102      l.s052.pk.r female   prolong
    ## 103     l.s080.blk.r   male   prolong
    ## 104     l.s120.y.blk female      bldg
    ## 105       l.s166.o.w female   prolong
    ## 106     l.s280.g.blk   male  m.inc.d9
    ## 107         o.d.s009   male  m.inc.d9
    ## 108          o.s.w.r   male       lay
    ## 109     o.s010.r.blk   male  m.inc.d3
    ## 110     o.s084.w.blk female m.inc.d17
    ## 111        o114.blu9   male   prolong
    ## 112    o152.o120.w42   male        n5
    ## 113       o156.w80.x female    inc.d3
    ## 114      o165.w122.x female    inc.d3
    ## 115       o169.r28.x female m.inc.d17
    ## 116      o172.w115.x female     hatch
    ## 117      o173.w179.x female    inc.d3
    ## 118             o3.x   male      m.n2
    ## 119        o35.r51.x female   inc.d17
    ## 120        o36.r62.x female  m.inc.d9
    ## 121      o38.blu29.x female      bldg
    ## 122        o39.y77.x   male     hatch
    ## 123      o44.blu26.x   male     hatch
    ## 124       o45.g128.x female  m.inc.d9
    ## 125       o48.r197.x   male    inc.d3
    ## 126            o49.x   male    inc.d9
    ## 127        o52.blu53 female   inc.d17
    ## 128          o57.g59   male    inc.d9
    ## 129        o59.blu64   male m.inc.d17
    ## 130            o73.x female    inc.d9
    ## 131     p.g.blk.s040 female     early
    ## 132         pk.s.d.g   male   prolong
    ## 133      pk.s011.o.y female m.inc.d17
    ## 134      pk.s054.d.g female      m.n2
    ## 135      pk.s055.d.l female     early
    ## 136    pk.s238.blk.w   male       lay
    ## 137      pk.w.s141.o   male       lay
    ## 138    pu.blk.s102.y   male     early
    ## 139         pu.s.o.r   male      m.n2
    ## 140 r.r.x.ATLAS.R2XR female   control
    ## 141      r.r.x.ATLAS female   control
    ## 142    r.s005.pk.blk   male       lay
    ## 143     r.s035.y.blk female  m.inc.d3
    ## 144       r.s056.g.o female      bldg
    ## 145      r.s057.g.pk   male    extend
    ## 146       r.s058.d.l   male     early
    ## 147       r.s059.d.o   male      bldg
    ## 148     r.s086.l.blk   male    extend
    ## 149    r.s116.blk.pu   male       lay
    ## 150       r.s131.o.d female    extend
    ## 151       r.s171.l.w female        n9
    ## 152       r.s172.l.y   male    extend
    ## 153     r.y.s007.blk   male        n9
    ## 154       r176.blu54   male   inc.d17
    ## 155         r183.o22 female     hatch
    ## 156       r190.o43.x   male       lay
    ## 157           r194.x female   prolong
    ## 158           r195.x   male        n9
    ## 159           r196.x   male  m.inc.d9
    ## 160  r27.w111.blu125 female    inc.d3
    ## 161     r30.w112.r46 female    inc.d9
    ## 162       r36.w184.x female    inc.d9
    ## 163       r37.w100.x   male        n9
    ## 164             r4.x   male     early
    ## 165        r41.w99.x   male     hatch
    ## 166            r45.X   male    inc.d9
    ## 167            r45.x   male    inc.d9
    ## 168       r49.w189.x female   inc.d17
    ## 169             r6.x female   control
    ## 170        r72.y83.x   male     hatch
    ## 171       r73.g127.x female    inc.d3
    ## 172            r81.x   male   prolong
    ## 173          r83.g45 female      bldg
    ## 174            r84.x   male  m.inc.d9
    ## 175          r85.g39   male  m.inc.d9
    ## 176            r90.x   male   prolong
    ## 177        r95.blu99 female        n9
    ## 178           s.o.pk female       lay
    ## 179          s.pu.pk female   prolong
    ## 180    s.pu148.blk.r   male      bldg
    ## 181        s.x.ATLAS female   control
    ## 182           s002.x female     early
    ## 183     s038.g.d.blk female m.inc.d17
    ## 184     s044.blk.d.r   male     early
    ## 185     s062.d.blk.g   male      m.n2
    ## 186     s063.d.blk.l female      bldg
    ## 187    s064.g.blk.pu   male m.inc.d17
    ## 188       s065.l.d.o   male      bldg
    ## 189       s066.l.d.r   male      bldg
    ## 190       s067.o.l.y   male     early
    ## 191     s069.pk.pu.g female  m.inc.d3
    ## 192     s071.pu.g.pk   male  m.inc.d3
    ## 193   s089.blk.pk.pu female    extend
    ## 194    s090.blk.pk.w   male      m.n2
    ## 195     s091.blk.r.g female     early
    ## 196     s092.blk.r.o female      bldg
    ## 197     s093.blk.y.d female  m.inc.d3
    ## 198     s095.g.blk.o female       lay
    ## 199    s096.g.blk.pk   male m.inc.d17
    ## 200    s100.l.pk.blk female m.inc.d17
    ## 201    s103.y.pk.blk female  m.inc.d3
    ## 202       s136.d.w.o female       lay
    ## 203     s137.g.blk.y female    extend
    ## 204     s139.l.blk.w   male  m.inc.d3
    ## 205     s142.o.pk.pu female       lay
    ## 206     s150.w.g.blk   male       lay
    ## 207   s175.blk.pu.pk   male  m.inc.d3
    ## 208    s176.blk.pu.r female       lay
    ## 209      s186.l.o.pk female    extend
    ## 210       s187.l.o.r   male        n9
    ## 211    s243.blk.pk.r   male       lay
    ## 212    s333.y.blk.pk female      m.n2
    ## 213          w191.r1 female   control
    ## 214            w34.x   male    inc.d9
    ## 215  x.blk.blk.ATLAS   male   control
    ## 216          x.blk16   male        n9
    ## 217          x.blk20 female   prolong
    ## 218          x.blk23   male      m.n2
    ## 219    x.blu.o.ATLAS   male   control
    ## 220     x.blu101.w43 female    inc.d9
    ## 221    x.blu102.w105 female    inc.d3
    ## 222    x.blu106.o153   male    inc.d9
    ## 223    x.blu109.w121 female        n5
    ## 224    x.blu116.w107 female   inc.d17
    ## 225     x.blu117.w89   male   inc.d17
    ## 226     x.blu122.r66 female    inc.d9
    ## 227      x.blu23.w14   male        n9
    ## 228          x.blu25 female  m.inc.d9
    ## 229          x.blu30   male        n5
    ## 230      x.blu42.o28   male    inc.d3
    ## 231     x.blu43.g132 female        n9
    ## 232          x.blu57   male   prolong
    ## 233       x.blu6.y80 female       lay
    ## 234          x.blu69   male    extend
    ## 235          x.blu73   male    extend
    ## 236          x.blu76   male      m.n2
    ## 237          x.blu94 female m.inc.d17
    ## 238          x.blu96 female      m.n2
    ## 239        x.g.ATLAS   male   control
    ## 240      x.g.g.ATLAS female   control
    ## 241      x.g.g.ATLAS   male   control
    ## 242    x.g.g.g.ATLAS   male   control
    ## 243       x.g13.w109   male    inc.d9
    ## 244       x.g14.w199   male   inc.d17
    ## 245     x.g147.blu28   male    inc.d3
    ## 246            x.g26 female     early
    ## 247            x.g37 female        n5
    ## 248         x.g4.w50 female        n9
    ## 249            x.g43 female        n5
    ## 250            x.g49 female        n5
    ## 251         x.g5.w79   male  m.inc.d3
    ## 252            x.g50 female   prolong
    ## 253            x.g70   male     hatch
    ## 254        x.g9.o166 female    inc.d9
    ## 255           x.o117   male  m.inc.d9
    ## 256       x.o159.w90 female   inc.d17
    ## 257      x.o160.w102   male     hatch
    ## 258      x.o163.w101   male    inc.d3
    ## 259      x.o164.w123   male        n5
    ## 260       x.o171.w45 female      m.n2
    ## 261       x.o175.g21 female        n5
    ## 262             x.o2   male        n9
    ## 263       x.o30.g134   male      bldg
    ## 264      x.o37.blu50 female     hatch
    ## 265        x.o40.r70   male m.inc.d17
    ## 266        x.o47.y82   male    inc.d9
    ## 267             x.o4   male      m.n2
    ## 268            x.o61 female    extend
    ## 269            x.o68   male        n5
    ## 270            x.o70 female        n5
    ## 271        x.r10.w18   male m.inc.d17
    ## 272           x.r178   male     hatch
    ## 273           x.r180 female  m.inc.d3
    ## 274           x.r181   male        n5
    ## 275           x.r185   male  m.inc.d3
    ## 276        x.r29.w96   male   inc.d17
    ## 277       x.r33.w183 female    inc.d3
    ## 278        x.r39.g10 female      bldg
    ## 279        x.r44.w95 female     hatch
    ## 280       x.r48.y139 female   inc.d17
    ## 281        x.r50.w97 female        n5
    ## 282       x.r64.g140   male    inc.d3
    ## 283      x.r67.blu35   male      bldg
    ## 284        x.r74.o29 female      m.n2
    ## 285            x.r76 female  m.inc.d9
    ## 286           x.s188 female     early
    ## 287           x.w178 female        n9
    ## 288      x.w192.o157   male    inc.d9
    ## 289            x.w51 female       lay
    ## 290             x.w6 female        n9
    ## 291      x.y.s.ATLAS   male   control
    ## 292           x.y109 female    inc.d9
    ## 293       x.y119.w11 female m.inc.d17
    ## 294       x.y132.w76   male   inc.d17
    ## 295      x.y138.w176 female        n9
    ## 296      x.y141.w116   male    inc.d9
    ## 297       x.y145.r55 female   prolong
    ## 298            x.y51   male     early
    ## 299             x.y5 female     early
    ## 300            x.y64 female    extend
    ## 301            x.y79   male   prolong
    ## 302            x.y90 female     hatch
    ## 303       x.y93.g126 female    inc.d9
    ## 304             x.y9 female        n9
    ## 305       y.s156.o.r female       lay
    ## 306          y.w.blk   male  m.inc.d9
    ## 307            y12.x female  m.inc.d9
    ## 308       y126.w92.x female   inc.d17
    ## 309       y128.g23.x female    inc.d9
    ## 310       y128.g23.x female  m.inc.d9
    ## 311           y129.x   male        n9
    ## 312            y13.x female    inc.d3
    ## 313      y130.o170.x female   inc.d17
    ## 314      y131.w185.x   male        n9
    ## 315     y133.w77.r58   male   inc.d17
    ## 316    y135.blu107.x female   inc.d17
    ## 317           y136.x female   inc.d17
    ## 318      y140.w119.x female    inc.d9
    ## 319       y146.r32.x female   prolong
    ## 320    y148.w81.g145   male m.inc.d17
    ## 321       y149.r52.x   male    inc.d3
    ## 322            y15.x female     hatch
    ## 323            y18.x   male  m.inc.d3
    ## 324             y4.x   male m.inc.d17
    ## 325            y55.x   male     early
    ## 326           y6.o54 female        n5
    ## 327            y63.x   male  m.inc.d9
    ## 328           y7.g58 female     hatch
    ## 329        y85.r71.x female m.inc.d17
    ## 330       y94.g133.x female        n5
    ## 331       y95.g131.x   male    inc.d9
    ## 332            y97.x female        n9
    ## 333          y98.g54 female      m.n2
    ## 334        y98.o50.x   male    inc.d3

save files for downstream use
-----------------------------

    write.csv(countData, "../results/00_counts.csv")
    write.csv(isoforms, "../results/00_geneswithisoforms.csv", row.names = T)

    write.csv(colData, "../metadata/00_colData.csv")
    write.csv(birds, "../metadata/00_birds.csv")
    write.csv(geneinfo, "../metadata/00_geneinfo.csv", row.names = TRUE)

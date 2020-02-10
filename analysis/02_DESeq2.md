Warning: This script takes a very long time to run.
===================================================

DESeq2 was not designed to run on 300 samples. But, I really like it, so
I do it anyways. Some of these commands take like 15 min to run using 6
cores.

DEseq2 on all chaacterization and manipulations
-----------------------------------------------

    # import "colData" which contains sample information and "countData" which contains read counts
    countData <- read.csv("../results/00_countData_characterization.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)
    head(geneinfo)

    FALSE       Name geneid       entrezid
    FALSE 1    EDNRB 408082 NP_001001127.1
    FALSE 2  CYP26A1 408183 NP_001001129.1
    FALSE 3    CFDP1 374073 NP_001001189.1
    FALSE 4    AvBD7 407777 NP_001001194.1
    FALSE 5     KRT5 407779 NP_001001195.1
    FALSE 6 HSD11B1L 408034 NP_001001201.1

    # craete variable that will be critical for subset later on
    colData <- read.csv("../metadata/00_colData_characterization.csv", header = T, row.names = 1)
    colData$sextissue <- as.factor(paste(colData$sex, colData$tissue, sep = "_"))
    colData$treatment <- factor(colData$treatment, levels = charlevels)
    colData$tissue <- factor(colData$tissue, levels = tissuelevel)
    levels(colData$treatment)

    FALSE [1] "control" "bldg"    "lay"     "inc.d3"  "inc.d9"  "inc.d17" "hatch"  
    FALSE [8] "n5"      "n9"

    levels(colData$sex)

    FALSE [1] "female" "male"

    levels(colData$sextissue)

    FALSE [1] "female_gonad"        "female_hypothalamus" "female_pituitary"   
    FALSE [4] "male_gonad"          "male_hypothalamus"   "male_pituitary"

    levels(colData$tissue)

    FALSE [1] "hypothalamus" "pituitary"    "gonad"

    for(i in levels(colData$sextissue)){
      
      newcolData <- subsetcolData2(colData, i)
      
      # save counts that match colData
      savecols <- as.character(newcolData$V1) 
      savecols <- as.vector(savecols) 
      
      newcountData <- countData %>% dplyr::select(one_of(savecols)) 
      
      dds <- DESeqDataSetFromMatrix(countData = newcountData,
                                    colData = newcolData,
                                    design = ~ treatment )
      dds <- dds[rowSums(counts(dds) > 1) >= 10]  # filter more than sample with less 0 counts
      print(dds)
      print(dim(dds))
      dds <- DESeq(dds, parallel = TRUE) # Differential expression analysis
      
      vsd <- as.data.frame(assay(vst(dds, blind=FALSE)))
      
      myfilename = paste0("../results/DEseq2/", i, "_vsd.csv")
      write.csv(vsd, myfilename)
      
      #return(dds)
      #return(vsd)
      print(head(vsd))

      # save differential gene expression results
      control.bldg <- createDEGdfsave("bldg", "control", i)
      bldg.lay <- createDEGdfsave("lay", "bldg", i)
      lay.inc.d3 <- createDEGdfsave("inc.d3", "lay",  i) 
      inc.d3.inc.d9 <- createDEGdfsave("inc.d9", "inc.d3", i) 
      inc.d9.inc.d17 <- createDEGdfsave("inc.d17", "inc.d9", i)
      inc.d17.hatch <- createDEGdfsave( "hatch", "inc.d17", i) 
      hatch.n5 <- createDEGdfsave( "n5", "hatch",  i) 
      n5.n9 <- createDEGdfsave("n9", "n5",  i) 

    }

    FALSE class: DESeqDataSet 
    FALSE dim: 13675 98 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13675): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(98): L.G118_female_gonad_control
    FALSE   R.G106_female_gonad_control ... y94.g133.x_female_gonad_n5
    FALSE   y97.x_female_gonad_n9
    FALSE colData names(8): V1 bird ... study sextissue
    FALSE [1] 13675    98
    FALSE        L.G118_female_gonad_control R.G106_female_gonad_control
    FALSE A2ML1                    14.283612                    7.461476
    FALSE A2ML2                     4.531659                    4.739738
    FALSE A2ML3                     9.530241                    8.134460
    FALSE A2ML4                     5.260736                    5.447566
    FALSE A4GALT                    6.026943                    6.105466
    FALSE A4GNT                     6.118576                    5.906769
    FALSE        R.R20_female_gonad_control R.R9_female_gonad_control
    FALSE A2ML1                   14.506548                 11.380727
    FALSE A2ML2                    4.531659                  4.531659
    FALSE A2ML3                    9.372404                  8.316314
    FALSE A2ML4                    5.254924                  4.531659
    FALSE A4GALT                   6.248399                  5.178447
    FALSE A4GNT                    5.937142                  5.690995
    FALSE        R.W44_female_gonad_control blk.s061.pu.y_female_gonad_inc.d9
    FALSE A2ML1                   12.839160                          6.691119
    FALSE A2ML2                    4.895462                          4.531659
    FALSE A2ML3                   10.508262                          7.875762
    FALSE A2ML4                    4.828964                          4.531659
    FALSE A4GALT                   7.428649                          6.073973
    FALSE A4GNT                    5.934774                          5.138279
    FALSE        blk11.x_female_gonad_bldg blk21.x_female_gonad_hatch
    FALSE A2ML1                   7.280805                   7.959591
    FALSE A2ML2                   5.129171                   5.031630
    FALSE A2ML3                   9.818386                   7.596957
    FALSE A2ML4                   4.531659                   5.142492
    FALSE A4GALT                  6.101179                   5.235274
    FALSE A4GNT                   5.370840                   6.270667
    FALSE        blk4.x_female_gonad_n9 blu.o.x.ATLAS_female_gonad_control
    FALSE A2ML1                7.372385                           7.880658
    FALSE A2ML2                4.531659                           5.117350
    FALSE A2ML3                9.190966                           7.016450
    FALSE A2ML4                5.050916                           4.795024
    FALSE A4GALT               7.230918                           5.313259
    FALSE A4GNT                5.050916                           5.877825
    FALSE        blu103.x_female_gonad_hatch.NYNO blu124.w180.x_female_gonad_hatch
    FALSE A2ML1                          6.553298                         7.221184
    FALSE A2ML2                          4.531659                         5.028212
    FALSE A2ML3                         10.402267                         7.288068
    FALSE A2ML4                          4.531659                         4.531659
    FALSE A4GALT                         6.206199                         5.968551
    FALSE A4GNT                          5.122310                         5.510775
    FALSE        blu36.w16_female_gonad_n9 blu38.g135.x_female_gonad_bldg
    FALSE A2ML1                   7.347426                       7.630063
    FALSE A2ML2                   4.976848                       4.531659
    FALSE A2ML3                   9.890611                      10.158414
    FALSE A2ML4                   5.296813                       4.531659
    FALSE A4GALT                  7.176727                       6.393136
    FALSE A4GNT                   5.463516                       4.531659
    FALSE        blu39.o26.x_female_gonad_inc.d3 blu47.y96.x_female_gonad_inc.d9
    FALSE A2ML1                         6.146270                        7.114329
    FALSE A2ML2                         4.531659                        4.818260
    FALSE A2ML3                         7.714219                        8.089199
    FALSE A2ML4                         5.438711                        4.936312
    FALSE A4GALT                        5.852668                        6.165389
    FALSE A4GNT                         5.110925                        5.467287
    FALSE        blu55.g51_female_gonad_n5 g.blk.s004.pk_female_gonad_lay
    FALSE A2ML1                   6.174115                       7.126960
    FALSE A2ML2                   5.051246                       4.820014
    FALSE A2ML3                   8.779015                       7.864730
    FALSE A2ML4                   4.900050                       5.232254
    FALSE A4GALT                  6.331972                       6.329055
    FALSE A4GNT                   4.900050                       5.385642
    FALSE        g.s078.blk.o_female_gonad_lay g.x.ATLAS_female_gonad_control
    FALSE A2ML1                       6.740466                       6.973951
    FALSE A2ML2                       4.531659                       4.769373
    FALSE A2ML3                       8.114934                       8.373470
    FALSE A2ML4                       5.113854                       5.198815
    FALSE A4GALT                      6.127530                       6.040540
    FALSE A4GNT                       5.519028                       5.494895
    FALSE        g141.blu27.x_female_gonad_bldg g142.r40.x_female_gonad_inc.d17
    FALSE A2ML1                       11.342267                        7.082760
    FALSE A2ML2                        4.911676                        4.969178
    FALSE A2ML3                       10.029122                        7.996074
    FALSE A2ML4                        5.067554                        5.148079
    FALSE A4GALT                       8.241980                        5.871817
    FALSE A4GNT                        4.531659                        5.773439
    FALSE        g6.w197.x_female_gonad_inc.d3 g75.x_female_gonad_inc.d9
    FALSE A2ML1                       7.208665                  6.998577
    FALSE A2ML2                       4.924532                  4.895546
    FALSE A2ML3                       8.151120                  7.946923
    FALSE A2ML4                       6.272538                  5.336998
    FALSE A4GALT                      6.833575                  6.252144
    FALSE A4GNT                       5.310349                  5.158656
    FALSE        l.s120.y.blk_female_gonad_bldg o156.w80.x_female_gonad_inc.d3
    FALSE A2ML1                        7.350783                       6.943549
    FALSE A2ML2                        4.936018                       4.531659
    FALSE A2ML3                        7.938154                       7.501049
    FALSE A2ML4                        5.101668                       4.901892
    FALSE A4GALT                       5.873394                       7.033523
    FALSE A4GNT                        5.683577                       5.619457
    FALSE        o165.w122.x_female_gonad_inc.d3.NYNO o172.w115.x_female_gonad_hatch
    FALSE A2ML1                              6.875564                       9.612026
    FALSE A2ML2                              4.531659                       4.531659
    FALSE A2ML3                              7.351154                       8.719698
    FALSE A2ML4                              4.841947                       4.912748
    FALSE A4GALT                             5.423759                       5.862051
    FALSE A4GNT                              5.448795                       4.912748
    FALSE        o173.w179.x_female_gonad_inc.d3 o35.r51.x_female_gonad_inc.d17
    FALSE A2ML1                         6.929956                       6.868086
    FALSE A2ML2                         4.531659                       4.839667
    FALSE A2ML3                         7.694441                       8.229033
    FALSE A2ML4                         5.076318                       4.839667
    FALSE A4GALT                        6.335917                       5.713287
    FALSE A4GNT                         5.464322                       5.365138
    FALSE        o38.blu29.x_female_gonad_bldg o52.blu53_female_gonad_inc.d17
    FALSE A2ML1                      12.357107                       7.526111
    FALSE A2ML2                       4.764214                       4.531659
    FALSE A2ML3                       8.561910                       7.974061
    FALSE A2ML4                       4.995276                       5.719467
    FALSE A4GALT                      7.436856                       5.646458
    FALSE A4GNT                       5.501039                       5.787389
    FALSE        o73.x_female_gonad_inc.d9 r.r.x.ATLAS.R2XR_female_gonad_control
    FALSE A2ML1                   6.954354                              7.410583
    FALSE A2ML2                   4.822520                              4.531659
    FALSE A2ML3                   7.631282                              7.266141
    FALSE A2ML4                   4.531659                              5.125138
    FALSE A4GALT                  5.823242                              5.365251
    FALSE A4GNT                   5.823242                              6.267240
    FALSE        r.r.x.ATLAS_female_gonad_control r.s056.g.o_female_gonad_bldg
    FALSE A2ML1                          7.731335                     8.137087
    FALSE A2ML2                          4.917632                     4.531659
    FALSE A2ML3                          7.332560                     8.894264
    FALSE A2ML4                          4.847117                     5.240196
    FALSE A4GALT                         5.758686                     5.962141
    FALSE A4GNT                          5.683508                     5.347154
    FALSE        r.s171.l.w_female_gonad_n9 r183.o22_female_gonad_hatch
    FALSE A2ML1                    7.169694                    6.549634
    FALSE A2ML2                    4.531659                    4.531659
    FALSE A2ML3                    9.237050                    8.221550
    FALSE A2ML4                    4.965345                    4.531659
    FALSE A4GALT                   6.668451                    6.277001
    FALSE A4GNT                    5.256757                    5.075079
    FALSE        r27.w111.blu125_female_gonad_inc.d3
    FALSE A2ML1                             7.482238
    FALSE A2ML2                             4.890941
    FALSE A2ML3                             7.881180
    FALSE A2ML4                             5.038465
    FALSE A4GALT                            6.137076
    FALSE A4GNT                             5.719370
    FALSE        r30.w112.r46_female_gonad_inc.d9 r36.w184.x_female_gonad_inc.d9
    FALSE A2ML1                          7.122328                       7.673725
    FALSE A2ML2                          4.531659                       4.955095
    FALSE A2ML3                          7.854510                       8.694257
    FALSE A2ML4                          4.838954                       5.369740
    FALSE A4GALT                         5.891613                       6.533917
    FALSE A4GNT                          5.061927                       5.551237
    FALSE        r49.w189.x_female_gonad_inc.d17 r6.x_female_gonad_control.NYNO
    FALSE A2ML1                         6.850738                       7.888482
    FALSE A2ML2                         5.047763                       5.063920
    FALSE A2ML3                         7.666638                       8.941283
    FALSE A2ML4                         5.195630                       4.944875
    FALSE A4GALT                        6.141433                       6.435741
    FALSE A4GNT                         5.857305                       5.579072
    FALSE        r73.g127.x_female_gonad_inc.d3 r83.g45_female_gonad_bldg
    FALSE A2ML1                        5.961295                  6.983967
    FALSE A2ML2                        4.692123                  4.531659
    FALSE A2ML3                        7.705150                  8.752019
    FALSE A2ML4                        5.148741                  4.531659
    FALSE A4GALT                       5.922137                  6.768087
    FALSE A4GNT                        5.084423                  4.531659
    FALSE        r95.blu99_female_gonad_n9 s.o.pk_female_gonad_lay
    FALSE A2ML1                   7.187948                7.191174
    FALSE A2ML2                   4.531659                4.859152
    FALSE A2ML3                   7.487660                8.546883
    FALSE A2ML4                   5.258512                4.763480
    FALSE A4GALT                  6.194090                5.764871
    FALSE A4GNT                   5.765971                5.471787
    FALSE        s.x.ATLAS_female_gonad_control s063.d.blk.l_female_gonad_bldg
    FALSE A2ML1                        6.796829                       7.300377
    FALSE A2ML2                        4.531659                       4.531659
    FALSE A2ML3                        8.364321                       9.636211
    FALSE A2ML4                        5.081611                       4.860654
    FALSE A4GALT                       6.049832                       6.101587
    FALSE A4GNT                        5.255959                       6.012289
    FALSE        s092.blk.r.o_female_gonad_bldg s095.g.blk.o_female_gonad_lay
    FALSE A2ML1                        6.975940                     14.606315
    FALSE A2ML2                        4.829118                      4.925513
    FALSE A2ML3                        9.968505                      9.272998
    FALSE A2ML4                        4.829118                      5.335688
    FALSE A4GALT                       6.200036                      7.014369
    FALSE A4GNT                        5.988415                      5.639049
    FALSE        s136.d.w.o_female_gonad_lay s142.o.pk.pu_female_gonad_lay
    FALSE A2ML1                     7.344510                     13.065644
    FALSE A2ML2                     4.811453                      4.531659
    FALSE A2ML3                     8.191153                     10.482584
    FALSE A2ML4                     5.211771                      5.343015
    FALSE A4GALT                    6.399669                      6.886832
    FALSE A4GNT                     5.558389                      4.531659
    FALSE        s176.blk.pu.r_female_gonad_lay w191.r1_female_gonad_control
    FALSE A2ML1                        6.470657                     7.463856
    FALSE A2ML2                        5.000485                     5.217890
    FALSE A2ML3                        8.023601                     8.831896
    FALSE A2ML4                        5.562533                     4.531659
    FALSE A4GALT                       6.619311                     5.699347
    FALSE A4GNT                        5.419281                     5.868727
    FALSE        x.blu101.w43_female_gonad_inc.d9 x.blu102.w105_female_gonad_inc.d3
    FALSE A2ML1                          7.547908                          7.160638
    FALSE A2ML2                          4.531659                          4.890275
    FALSE A2ML3                          7.883496                          7.746265
    FALSE A2ML4                          4.899643                          5.528635
    FALSE A4GALT                         5.818917                          6.166986
    FALSE A4GNT                          5.421314                          5.528635
    FALSE        x.blu109.w121_female_gonad_n5 x.blu116.w107_female_gonad_inc.d17
    FALSE A2ML1                       7.119564                           6.341587
    FALSE A2ML2                       4.531659                           4.531659
    FALSE A2ML3                       9.032447                           8.274810
    FALSE A2ML4                       4.995946                           4.531659
    FALSE A4GALT                      6.425291                           5.850647
    FALSE A4GNT                       5.329111                           4.531659
    FALSE        x.blu122.r66_female_gonad_inc.d9 x.blu43.g132_female_gonad_n9
    FALSE A2ML1                          6.465614                     6.521496
    FALSE A2ML2                          4.531659                     4.841583
    FALSE A2ML3                          7.591098                     9.182631
    FALSE A2ML4                          5.508311                     5.342505
    FALSE A4GALT                         5.784950                     6.556930
    FALSE A4GNT                          5.084716                     5.219472
    FALSE        x.blu6.y80_female_gonad_lay x.g.g.ATLAS_female_gonad_control
    FALSE A2ML1                     7.033121                         7.979655
    FALSE A2ML2                     4.531659                         4.971044
    FALSE A2ML3                     7.653809                         7.445875
    FALSE A2ML4                     4.788302                         5.199451
    FALSE A4GALT                    5.652511                         5.199451
    FALSE A4GNT                     5.476233                         6.103192
    FALSE        x.g37_female_gonad_n5 x.g4.w50_female_gonad_n9
    FALSE A2ML1               9.613434                 8.202535
    FALSE A2ML2               4.531659                 4.531659
    FALSE A2ML3               9.271051                 9.904864
    FALSE A2ML4               5.053654                 4.868917
    FALSE A4GALT              7.084068                 6.593792
    FALSE A4GNT               5.350480                 5.912475
    FALSE        x.g43_female_gonad_n5 x.g49_female_gonad_n5
    FALSE A2ML1               7.343343              6.665471
    FALSE A2ML2               4.884286              4.936573
    FALSE A2ML3               7.782633              7.564011
    FALSE A2ML4               4.531659              4.882611
    FALSE A4GALT              6.233743              5.861696
    FALSE A4GNT               5.569642              5.256382
    FALSE        x.g9.o166_female_gonad_inc.d9.NYNO x.o159.w90_female_gonad_inc.d17
    FALSE A2ML1                            6.064444                        6.848258
    FALSE A2ML2                            4.785247                        4.531659
    FALSE A2ML3                            8.274365                        9.390392
    FALSE A2ML4                            4.531659                        5.089783
    FALSE A4GALT                           5.894192                        6.978855
    FALSE A4GNT                            5.197529                        5.406193
    FALSE        x.o175.g21_female_gonad_n5 x.o37.blu50_female_gonad_hatch
    FALSE A2ML1                    7.052315                       7.012922
    FALSE A2ML2                    4.531659                       4.531659
    FALSE A2ML3                    8.400251                       7.650540
    FALSE A2ML4                    5.073364                       4.531659
    FALSE A4GALT                   5.945212                       6.381260
    FALSE A4GNT                    5.559900                       5.666459
    FALSE        x.o70_female_gonad_n5 x.r33.w183_female_gonad_inc.d3
    FALSE A2ML1               7.934962                       7.029733
    FALSE A2ML2               4.531659                       4.531659
    FALSE A2ML3               8.738916                       8.181880
    FALSE A2ML4               4.531659                       5.482322
    FALSE A4GALT              6.489113                       5.853800
    FALSE A4GNT               5.525997                       5.891756
    FALSE        x.r39.g10_female_gonad_bldg x.r44.w95_female_gonad_hatch
    FALSE A2ML1                     6.651175                     7.424698
    FALSE A2ML2                     4.531659                     5.005997
    FALSE A2ML3                     8.430671                     7.446108
    FALSE A2ML4                     5.194202                     5.276728
    FALSE A4GALT                    6.553520                     5.276728
    FALSE A4GNT                     4.531659                     5.834605
    FALSE        x.r48.y139_female_gonad_inc.d17 x.r50.w97_female_gonad_n5
    FALSE A2ML1                         8.158761                  7.373729
    FALSE A2ML2                         4.531659                  4.531659
    FALSE A2ML3                         8.922181                  9.123366
    FALSE A2ML4                         5.249184                  4.531659
    FALSE A4GALT                        6.077643                  6.172592
    FALSE A4GNT                         5.536452                  5.301746
    FALSE        x.w178_female_gonad_n9 x.w51_female_gonad_lay x.w6_female_gonad_n9
    FALSE A2ML1                7.151404               7.410860            13.734145
    FALSE A2ML2                4.951602               4.742396             4.531659
    FALSE A2ML3               10.052793               9.190386             9.553901
    FALSE A2ML4                5.192209               4.896021             5.051059
    FALSE A4GALT               6.866300               6.504533             6.636179
    FALSE A4GNT                5.254023               5.544070             4.744657
    FALSE        x.y109_female_gonad_inc.d9 x.y138.w176_female_gonad_n9
    FALSE A2ML1                    6.883500                    9.542857
    FALSE A2ML2                    4.531659                    4.885100
    FALSE A2ML3                    7.670338                    9.776698
    FALSE A2ML4                    5.191635                    5.140834
    FALSE A4GALT                   6.025356                    6.294788
    FALSE A4GNT                    5.000348                    5.314323
    FALSE        x.y90_female_gonad_hatch x.y93.g126_female_gonad_inc.d9
    FALSE A2ML1                  6.535919                       7.397737
    FALSE A2ML2                  4.982171                       4.845842
    FALSE A2ML3                  7.588067                       8.223347
    FALSE A2ML4                  5.239750                       4.531659
    FALSE A4GALT                 6.263264                       6.284288
    FALSE A4GNT                  5.366233                       5.753987
    FALSE        x.y9_female_gonad_n9 y.s156.o.r_female_gonad_lay
    FALSE A2ML1              6.774206                    6.887652
    FALSE A2ML2              4.783965                    4.784087
    FALSE A2ML3              9.505967                    9.036542
    FALSE A2ML4              4.531659                    4.888193
    FALSE A4GALT             6.293629                    6.336522
    FALSE A4GNT              5.034367                    5.321041
    FALSE        y126.w92.x_female_gonad_inc.d17 y128.g23.x_female_gonad_inc.d9
    FALSE A2ML1                         7.659478                       6.305112
    FALSE A2ML2                         4.531659                       4.864540
    FALSE A2ML3                         8.310981                       8.878365
    FALSE A2ML4                         5.137310                       4.864540
    FALSE A4GALT                        6.104744                       6.280179
    FALSE A4GNT                         5.807645                       4.531659
    FALSE        y13.x_female_gonad_inc.d3 y130.o170.x_female_gonad_inc.d17
    FALSE A2ML1                   8.165926                         7.880891
    FALSE A2ML2                   4.531659                         5.083347
    FALSE A2ML3                   7.570191                         7.826980
    FALSE A2ML4                   5.996375                         4.922942
    FALSE A4GALT                  6.627573                         5.792836
    FALSE A4GNT                   5.349573                         5.792836
    FALSE        y135.blu107.x_female_gonad_inc.d17 y136.x_female_gonad_inc.d17
    FALSE A2ML1                            8.041525                    6.673271
    FALSE A2ML2                            4.926136                    4.531659
    FALSE A2ML3                            7.801461                    7.468737
    FALSE A2ML4                            5.687609                    4.951127
    FALSE A4GALT                           5.483577                    5.619431
    FALSE A4GNT                            5.746759                    5.253217
    FALSE        y140.w119.x_female_gonad_inc.d9 y15.x_female_gonad_hatch
    FALSE A2ML1                         6.932792                 6.912107
    FALSE A2ML2                         4.531659                 4.531659
    FALSE A2ML3                         7.780123                 8.356013
    FALSE A2ML4                         4.845648                 4.828512
    FALSE A4GALT                        5.888115                 5.654816
    FALSE A4GNT                         5.716551                 5.252557
    FALSE        y6.o54_female_gonad_n5 y7.g58_female_gonad_hatch
    FALSE A2ML1                7.443395                  7.368964
    FALSE A2ML2                4.531659                  4.531659
    FALSE A2ML3                8.208459                  7.595324
    FALSE A2ML4                5.171117                  5.291982
    FALSE A4GALT               6.994504                  6.583516
    FALSE A4GNT                5.756778                  5.595162
    FALSE        y94.g133.x_female_gonad_n5 y97.x_female_gonad_n9
    FALSE A2ML1                    7.906523              7.616673
    FALSE A2ML2                    4.531659              4.895544
    FALSE A2ML3                    7.756493              8.533808
    FALSE A2ML4                    5.027514              4.713963
    FALSE A4GALT                   5.764247              5.975021
    FALSE A4GNT                    6.002932              5.615629
    FALSE 'data.frame': 7143 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13675 levels "A2ML1","A2ML2",..: 6146 6342 13302 6797 6076 2351 7900 2994 8621 11293 ...
    FALSE  $ padj     : num  1.06e-09 2.10e-06 7.27e-11 1.71e-06 2.68e-09 ...
    FALSE  $ logpadj  : num  8.98 5.68 10.14 5.77 8.57 ...
    FALSE  $ lfc      : num  8.32 8.29 7.37 7.28 6.53 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_gonad": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "control","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene         padj   logpadj      lfc    sextissue direction
    FALSE 1 LOC101749216 1.055376e-09  8.976593 8.317625 female_gonad      bldg
    FALSE 2 LOC107049904 2.101399e-06  5.677492 8.291206 female_gonad      bldg
    FALSE 3        WFDC2 7.268511e-11 10.138555 7.368025 female_gonad      bldg
    FALSE 4    LOC418356 1.714308e-06  5.765911 7.284289 female_gonad      bldg
    FALSE 5 LOC101747844 2.682055e-09  8.571532 6.532891 female_gonad      bldg
    FALSE 6      COL10A1 2.581711e-07  6.588092 6.178363 female_gonad      bldg
    FALSE 'data.frame': 789 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13675 levels "A2ML1","A2ML2",..: 6566 7955 8737 566 3552 4349 2351 1514 4184 8250 ...
    FALSE  $ padj     : num  4.76e-04 2.72e-03 9.58e-10 2.75e-03 4.76e-05 ...
    FALSE  $ logpadj  : num  3.32 2.57 9.02 2.56 4.32 ...
    FALSE  $ lfc      : num  9.72 5.82 5.45 4.41 4.25 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_gonad": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "bldg","NS","lay": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene         padj  logpadj      lfc    sextissue direction
    FALSE 1 LOC107053414 4.757321e-04 3.322638 9.720236 female_gonad       lay
    FALSE 2          MUC 2.717796e-03 2.565783 5.822901 female_gonad       lay
    FALSE 3        OVSTL 9.578604e-10 9.018698 5.453279 female_gonad       lay
    FALSE 4         AOC1 2.754271e-03 2.559993 4.413449 female_gonad       lay
    FALSE 5       ETNPPL 4.763302e-05 4.322092 4.254100 female_gonad       lay
    FALSE 6         GKN2 1.629372e-02 1.787980 3.988462 female_gonad       lay
    FALSE 'data.frame': 2108 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13675 levels "A2ML1","A2ML2",..: 6277 12368 7718 3682 3924 7717 4736 10895 8918 8653 ...
    FALSE  $ padj     : num  1.43e-06 2.70e-03 8.80e-03 1.21e-02 4.44e-04 ...
    FALSE  $ logpadj  : num  5.85 2.57 2.06 1.92 3.35 ...
    FALSE  $ lfc      : num  3.3 3.13 2.99 2.94 2.63 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_gonad": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "lay","NS","inc.d3": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene         padj  logpadj      lfc    sextissue direction
    FALSE 1 LOC107048991 1.428002e-06 5.845271 3.296910 female_gonad    inc.d3
    FALSE 2      TMEM196 2.697186e-03 2.569089 3.128755 female_gonad    inc.d3
    FALSE 3        MMP13 8.796695e-03 2.055680 2.987979 female_gonad    inc.d3
    FALSE 4      FAM155A 1.210598e-02 1.917000 2.935815 female_gonad    inc.d3
    FALSE 5       FER1L6 4.439508e-04 3.352665 2.633647 female_gonad    inc.d3
    FALSE 6        MMP10 3.023037e-04 3.519557 2.490842 female_gonad    inc.d3
    FALSE 'data.frame': 12 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13675 levels "A2ML1","A2ML2",..: 9122 8593 9531 2638 12876 8712 9651 7901 958 4509 ...
    FALSE  $ padj     : num  0.0431 0.0101 0.021 0.021 0.0101 ...
    FALSE  $ logpadj  : num  1.37 2 1.68 1.68 2 ...
    FALSE  $ lfc      : num  3.405 -0.396 -0.53 -0.849 -0.908 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_gonad": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "inc.d3","NS",..: 3 1 1 1 1 1 1 1 1 1 ...
    FALSE NULL
    FALSE    gene       padj  logpadj        lfc    sextissue direction
    FALSE 1  PI15 0.04307979 1.365726  3.4050247 female_gonad    inc.d9
    FALSE 2 NUTF2 0.01007557 1.996730 -0.3955877 female_gonad    inc.d3
    FALSE 3  PPT1 0.02098266 1.678140 -0.5301224 female_gonad    inc.d3
    FALSE 4 CTSL2 0.02098266 1.678140 -0.8485209 female_gonad    inc.d3
    FALSE 5 TXNIP 0.01007557 1.996730 -0.9084959 female_gonad    inc.d3
    FALSE 6  OSTN 0.01007557 1.996730 -1.2558668 female_gonad    inc.d3
    FALSE 'data.frame': 326 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13675 levels "A2ML1","A2ML2",..: 6566 2225 11823 6797 6444 2715 3148 558 5899 9417 ...
    FALSE  $ padj     : num  0.0207 0.0473 0.0329 0.0968 0.0372 ...
    FALSE  $ logpadj  : num  1.68 1.32 1.48 1.01 1.43 ...
    FALSE  $ lfc      : num  7.13 4.65 4.28 4.16 2.63 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_gonad": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "inc.d9","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene       padj  logpadj      lfc    sextissue direction
    FALSE 1 LOC107053414 0.02067737 1.684505 7.133079 female_gonad   inc.d17
    FALSE 2       CLDN34 0.04734759 1.324702 4.652447 female_gonad   inc.d17
    FALSE 3      SULT1C3 0.03286733 1.483236 4.279182 female_gonad   inc.d17
    FALSE 4    LOC418356 0.09682698 1.014004 4.164914 female_gonad   inc.d17
    FALSE 5 LOC107050953 0.03717379 1.429763 2.629593 female_gonad   inc.d17
    FALSE 6       CYP2W1 0.04507957 1.346020 2.530213 female_gonad   inc.d17
    FALSE 'data.frame': 2 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13675 levels "A2ML1","A2ML2",..: 6146 6566
    FALSE  $ padj     : num  2.24e-02 1.45e-31
    FALSE  $ logpadj  : num  1.65 30.84
    FALSE  $ lfc      : num  5.65 -24.21
    FALSE  $ sextissue: Factor w/ 1 level "female_gonad": 1 1
    FALSE  $ direction: Factor w/ 3 levels "inc.d17","NS",..: 3 1
    FALSE NULL
    FALSE           gene         padj  logpadj        lfc    sextissue direction
    FALSE 1 LOC101749216 2.242538e-02  1.64926   5.652204 female_gonad     hatch
    FALSE 2 LOC107053414 1.452391e-31 30.83792 -24.206587 female_gonad   inc.d17
    FALSE 'data.frame': 1 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13675 levels "A2ML1","A2ML2",..: 6566
    FALSE  $ padj     : num 1.27e-26
    FALSE  $ logpadj  : num 25.9
    FALSE  $ lfc      : num 22.9
    FALSE  $ sextissue: Factor w/ 1 level "female_gonad": 1
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","n5": 3
    FALSE NULL
    FALSE           gene         padj  logpadj      lfc    sextissue direction
    FALSE 1 LOC107053414 1.270497e-26 25.89603 22.85955 female_gonad        n5
    FALSE 'data.frame': 276 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13675 levels "A2ML1","A2ML2",..: 6342 8731 2351 10045 6565 8734 6797 6122 11572 8621 ...
    FALSE  $ padj     : num  5.94e-07 5.62e-02 2.98e-06 5.94e-07 4.36e-02 ...
    FALSE  $ logpadj  : num  6.23 1.25 5.53 6.23 1.36 ...
    FALSE  $ lfc      : num  10.73 7.05 6.89 6.66 5.11 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_gonad": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "n5","NS","n9": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene         padj  logpadj       lfc    sextissue direction
    FALSE 1 LOC107049904 5.938737e-07 6.226306 10.733938 female_gonad        n9
    FALSE 2        OVALX 5.624106e-02 1.249947  7.049173 female_gonad        n9
    FALSE 3      COL10A1 2.976195e-06 5.526339  6.894354 female_gonad        n9
    FALSE 4          RBP 5.938737e-07 6.226306  6.659174 female_gonad        n9
    FALSE 5 LOC107053409 4.364843e-02 1.360031  5.106673 female_gonad        n9
    FALSE 6       OvoDA1 5.874945e-02 1.230996  4.689597 female_gonad        n9
    FALSE class: DESeqDataSet 
    FALSE dim: 13532 95 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13532): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(95): L.G118_female_hypothalamus_control.NYNO
    FALSE   R.G106_female_hypothalamus_control ...
    FALSE   y94.g133.x_female_hypothalamus_n5.NYNO
    FALSE   y97.x_female_hypothalamus_n9
    FALSE colData names(8): V1 bird ... study sextissue
    FALSE [1] 13532    95
    FALSE        L.G118_female_hypothalamus_control.NYNO
    FALSE A2ML1                                 8.910092
    FALSE A2ML2                                 5.459875
    FALSE A2ML3                                13.054198
    FALSE A2ML4                                 5.459875
    FALSE A4GALT                                5.459875
    FALSE A4GNT                                 5.459875
    FALSE        R.G106_female_hypothalamus_control
    FALSE A2ML1                            8.917711
    FALSE A2ML2                            5.840032
    FALSE A2ML3                           13.045704
    FALSE A2ML4                            5.679781
    FALSE A4GALT                           5.459875
    FALSE A4GNT                            5.679781
    FALSE        R.R20_female_hypothalamus_control.NYNO
    FALSE A2ML1                                7.954052
    FALSE A2ML2                                5.873740
    FALSE A2ML3                               13.495181
    FALSE A2ML4                                5.459875
    FALSE A4GALT                               6.043196
    FALSE A4GNT                                5.459875
    FALSE        R.R9_female_hypothalamus_control R.W44_female_hypothalamus_control
    FALSE A2ML1                          9.090823                          7.616815
    FALSE A2ML2                          5.894477                          5.664646
    FALSE A2ML3                         13.094956                         12.977538
    FALSE A2ML4                          5.945320                          5.664646
    FALSE A4GALT                         5.836604                          5.868394
    FALSE A4GNT                          5.459875                          5.459875
    FALSE        blk.s061.pu.y_female_hypothalamus_inc.d9
    FALSE A2ML1                                  8.325287
    FALSE A2ML2                                  5.877348
    FALSE A2ML3                                 12.779245
    FALSE A2ML4                                  5.801135
    FALSE A4GALT                                 6.216050
    FALSE A4GNT                                  5.459875
    FALSE        blk11.x_female_hypothalamus_bldg blk21.x_female_hypothalamus_hatch
    FALSE A2ML1                          8.742237                          7.592432
    FALSE A2ML2                          6.227804                          5.459875
    FALSE A2ML3                         12.714354                         12.933571
    FALSE A2ML4                          5.847224                          5.459875
    FALSE A4GALT                         6.179242                          6.398491
    FALSE A4GNT                          5.459875                          5.459875
    FALSE        blk4.x_female_hypothalamus_n9
    FALSE A2ML1                       8.194487
    FALSE A2ML2                       5.459875
    FALSE A2ML3                      12.247535
    FALSE A2ML4                       5.459875
    FALSE A4GALT                      5.459875
    FALSE A4GNT                       5.459875
    FALSE        blu.o.x.ATLAS_female_hypothalamus_control
    FALSE A2ML1                                   9.892536
    FALSE A2ML2                                   5.459875
    FALSE A2ML3                                  12.910318
    FALSE A2ML4                                   6.267926
    FALSE A4GALT                                  6.443364
    FALSE A4GNT                                   5.459875
    FALSE        blu103.x_female_hypothalamus_hatch
    FALSE A2ML1                            8.346314
    FALSE A2ML2                            5.992508
    FALSE A2ML3                           12.948661
    FALSE A2ML4                            5.930200
    FALSE A4GALT                           6.208973
    FALSE A4GNT                            5.459875
    FALSE        blu124.w180.x_female_hypothalamus_hatch
    FALSE A2ML1                                 7.788130
    FALSE A2ML2                                 5.719754
    FALSE A2ML3                                12.745319
    FALSE A2ML4                                 5.643762
    FALSE A4GALT                                6.008595
    FALSE A4GNT                                 5.459875
    FALSE        blu36.w16_female_hypothalamus_n9
    FALSE A2ML1                          8.321657
    FALSE A2ML2                          5.878822
    FALSE A2ML3                         12.572339
    FALSE A2ML4                          5.702320
    FALSE A4GALT                         5.459875
    FALSE A4GNT                          5.459875
    FALSE        blu38.g135.x_female_hypothalamus_bldg
    FALSE A2ML1                               8.389455
    FALSE A2ML2                               5.829759
    FALSE A2ML3                              12.717841
    FALSE A2ML4                               5.673817
    FALSE A4GALT                              5.981557
    FALSE A4GNT                               5.459875
    FALSE        blu39.o26.x_female_hypothalamus_inc.d3.NYNO
    FALSE A2ML1                                     7.681612
    FALSE A2ML2                                     5.674371
    FALSE A2ML3                                    12.672341
    FALSE A2ML4                                     5.459875
    FALSE A4GALT                                    5.982895
    FALSE A4GNT                                     5.459875
    FALSE        blu47.y96.x_female_hypothalamus_inc.d9
    FALSE A2ML1                                7.763119
    FALSE A2ML2                                6.037268
    FALSE A2ML3                               12.664876
    FALSE A2ML4                                5.794710
    FALSE A4GALT                               6.037268
    FALSE A4GNT                                5.459875
    FALSE        blu55.g51_female_hypothalamus_n5
    FALSE A2ML1                          8.443677
    FALSE A2ML2                          5.986971
    FALSE A2ML3                         12.809212
    FALSE A2ML4                          5.910059
    FALSE A4GALT                         6.094001
    FALSE A4GNT                          5.619604
    FALSE        g.blk.s004.pk_female_hypothalamus_lay
    FALSE A2ML1                               8.816238
    FALSE A2ML2                               6.130263
    FALSE A2ML3                              13.001981
    FALSE A2ML4                               6.008864
    FALSE A4GALT                              6.192960
    FALSE A4GNT                               5.459875
    FALSE        g.s078.blk.o_female_hypothalamus_lay
    FALSE A2ML1                              8.788106
    FALSE A2ML2                              5.843767
    FALSE A2ML3                             12.681579
    FALSE A2ML4                              5.843767
    FALSE A4GALT                             6.189347
    FALSE A4GNT                              5.459875
    FALSE        g141.blu27.x_female_hypothalamus_bldg
    FALSE A2ML1                               8.670336
    FALSE A2ML2                               5.988065
    FALSE A2ML3                              12.605881
    FALSE A2ML4                               5.765957
    FALSE A4GALT                              5.958163
    FALSE A4GNT                               5.459875
    FALSE        g142.r40.x_female_hypothalamus_inc.d17
    FALSE A2ML1                                8.259190
    FALSE A2ML2                                5.785144
    FALSE A2ML3                               12.593884
    FALSE A2ML4                                5.690118
    FALSE A4GALT                               6.020909
    FALSE A4GNT                                5.459875
    FALSE        g6.w197.x_female_hypothalamus_inc.d3
    FALSE A2ML1                              8.569814
    FALSE A2ML2                              5.759420
    FALSE A2ML3                             12.486935
    FALSE A2ML4                              5.826413
    FALSE A4GALT                             5.976863
    FALSE A4GNT                              5.459875
    FALSE        g75.x_female_hypothalamus_inc.d9
    FALSE A2ML1                          8.852503
    FALSE A2ML2                          5.925205
    FALSE A2ML3                         13.143247
    FALSE A2ML4                          5.972174
    FALSE A4GALT                         6.148751
    FALSE A4GNT                          5.584741
    FALSE        l.s120.y.blk_female_hypothalamus_bldg
    FALSE A2ML1                               8.103393
    FALSE A2ML2                               5.912955
    FALSE A2ML3                              12.790205
    FALSE A2ML4                               5.830317
    FALSE A4GALT                              6.194823
    FALSE A4GNT                               5.459875
    FALSE        o156.w80.x_female_hypothalamus_inc.d3
    FALSE A2ML1                               8.101482
    FALSE A2ML2                               6.029043
    FALSE A2ML3                              12.706558
    FALSE A2ML4                               5.758511
    FALSE A4GALT                              6.204477
    FALSE A4GNT                               5.632498
    FALSE        o165.w122.x_female_hypothalamus_inc.d3
    FALSE A2ML1                                8.725135
    FALSE A2ML2                                5.873347
    FALSE A2ML3                               12.286782
    FALSE A2ML4                                5.837536
    FALSE A4GALT                               6.171255
    FALSE A4GNT                                5.459875
    FALSE        o172.w115.x_female_hypothalamus_hatch
    FALSE A2ML1                               7.265392
    FALSE A2ML2                               5.684264
    FALSE A2ML3                              12.631759
    FALSE A2ML4                               5.459875
    FALSE A4GALT                              5.934216
    FALSE A4GNT                               5.459875
    FALSE        o173.w179.x_female_hypothalamus_inc.d3
    FALSE A2ML1                                8.866187
    FALSE A2ML2                                6.040355
    FALSE A2ML3                               12.774482
    FALSE A2ML4                                5.771645
    FALSE A4GALT                               6.216344
    FALSE A4GNT                                5.680544
    FALSE        o35.r51.x_female_hypothalamus_inc.d17
    FALSE A2ML1                               6.673751
    FALSE A2ML2                               5.459875
    FALSE A2ML3                              13.120869
    FALSE A2ML4                               5.459875
    FALSE A4GALT                              5.809092
    FALSE A4GNT                               5.459875
    FALSE        o38.blu29.x_female_hypothalamus_bldg
    FALSE A2ML1                              8.684876
    FALSE A2ML2                              5.740837
    FALSE A2ML3                             12.612957
    FALSE A2ML4                              5.740837
    FALSE A4GALT                             5.740837
    FALSE A4GNT                              5.459875
    FALSE        o52.blu53_female_hypothalamus_inc.d17
    FALSE A2ML1                               8.216164
    FALSE A2ML2                               6.017325
    FALSE A2ML3                              12.844035
    FALSE A2ML4                               5.459875
    FALSE A4GALT                              6.199643
    FALSE A4GNT                               5.915966
    FALSE        o73.x_female_hypothalamus_inc.d9
    FALSE A2ML1                          8.661425
    FALSE A2ML2                          5.963987
    FALSE A2ML3                         13.020591
    FALSE A2ML4                          5.817238
    FALSE A4GALT                         5.963987
    FALSE A4GNT                          5.712891
    FALSE        r.r.x.ATLAS.R2XR_female_hypothalamus_control
    FALSE A2ML1                                      9.031397
    FALSE A2ML2                                      6.396777
    FALSE A2ML3                                     13.529809
    FALSE A2ML4                                      6.229218
    FALSE A4GALT                                     6.007066
    FALSE A4GNT                                      5.459875
    FALSE        r.r.x.ATLAS_female_hypothalamus_control
    FALSE A2ML1                                 9.030580
    FALSE A2ML2                                 6.396474
    FALSE A2ML3                                13.528849
    FALSE A2ML4                                 6.228966
    FALSE A4GALT                                6.006884
    FALSE A4GNT                                 5.459875
    FALSE        r.s056.g.o_female_hypothalamus_bldg
    FALSE A2ML1                             7.898488
    FALSE A2ML2                             5.791355
    FALSE A2ML3                            12.892718
    FALSE A2ML4                             5.927639
    FALSE A4GALT                            5.791355
    FALSE A4GNT                             5.459875
    FALSE        r.s171.l.w_female_hypothalamus_n9
    FALSE A2ML1                           8.566370
    FALSE A2ML2                           6.015017
    FALSE A2ML3                          12.769358
    FALSE A2ML4                           5.950128
    FALSE A4GALT                          5.950128
    FALSE A4GNT                           5.645937
    FALSE        r183.o22_female_hypothalamus_hatch
    FALSE A2ML1                            8.408671
    FALSE A2ML2                            5.652120
    FALSE A2ML3                           12.841899
    FALSE A2ML4                            5.731549
    FALSE A4GALT                           5.843517
    FALSE A4GNT                            5.652120
    FALSE        r27.w111.blu125_female_hypothalamus_inc.d3
    FALSE A2ML1                                    8.598622
    FALSE A2ML2                                    5.836930
    FALSE A2ML3                                   12.501635
    FALSE A2ML4                                    5.921021
    FALSE A4GALT                                   5.881138
    FALSE A4GNT                                    5.459875
    FALSE        r30.w112.r46_female_hypothalamus_inc.d9
    FALSE A2ML1                                 7.652770
    FALSE A2ML2                                 5.647375
    FALSE A2ML3                                12.517328
    FALSE A2ML4                                 5.647375
    FALSE A4GALT                                5.877969
    FALSE A4GNT                                 5.459875
    FALSE        r36.w184.x_female_hypothalamus_inc.d9
    FALSE A2ML1                               8.163416
    FALSE A2ML2                               5.781832
    FALSE A2ML3                              12.880743
    FALSE A2ML4                               5.914256
    FALSE A4GALT                              5.983836
    FALSE A4GNT                               5.459875
    FALSE        r49.w189.x_female_hypothalamus_inc.d17
    FALSE A2ML1                                8.475097
    FALSE A2ML2                                5.874459
    FALSE A2ML3                               12.527533
    FALSE A2ML4                                6.133095
    FALSE A4GALT                               5.993890
    FALSE A4GNT                                5.459875
    FALSE        r6.x_female_hypothalamus_control.NYNO
    FALSE A2ML1                               8.233463
    FALSE A2ML2                               5.459875
    FALSE A2ML3                              13.216136
    FALSE A2ML4                               5.768502
    FALSE A4GALT                              6.110299
    FALSE A4GNT                               5.459875
    FALSE        r73.g127.x_female_hypothalamus_inc.d3
    FALSE A2ML1                               8.601232
    FALSE A2ML2                               6.065995
    FALSE A2ML3                              12.339291
    FALSE A2ML4                               5.827081
    FALSE A4GALT                              6.038168
    FALSE A4GNT                               5.459875
    FALSE        r83.g45_female_hypothalamus_bldg r95.blu99_female_hypothalamus_n9
    FALSE A2ML1                          7.284646                         7.561976
    FALSE A2ML2                          5.459875                         5.459875
    FALSE A2ML3                         12.945319                        12.977434
    FALSE A2ML4                          5.833886                         5.701974
    FALSE A4GALT                         5.833886                         5.841998
    FALSE A4GNT                          5.459875                         5.701974
    FALSE        s.o.pk_female_hypothalamus_lay
    FALSE A2ML1                        7.614873
    FALSE A2ML2                        5.662164
    FALSE A2ML3                       12.779924
    FALSE A2ML4                        5.459875
    FALSE A4GALT                       6.094953
    FALSE A4GNT                        5.459875
    FALSE        s.x.ATLAS_female_hypothalamus_control
    FALSE A2ML1                               8.661674
    FALSE A2ML2                               6.138848
    FALSE A2ML3                              12.883979
    FALSE A2ML4                               5.987723
    FALSE A4GALT                              6.068269
    FALSE A4GNT                               5.459875
    FALSE        s063.d.blk.l_female_hypothalamus_bldg
    FALSE A2ML1                               8.046216
    FALSE A2ML2                               5.459875
    FALSE A2ML3                              12.533519
    FALSE A2ML4                               5.459875
    FALSE A4GALT                              5.755528
    FALSE A4GNT                               5.459875
    FALSE        s092.blk.r.o_female_hypothalamus_bldg
    FALSE A2ML1                               8.261017
    FALSE A2ML2                               5.902769
    FALSE A2ML3                              12.492760
    FALSE A2ML4                               5.641277
    FALSE A4GALT                              6.200159
    FALSE A4GNT                               5.459875
    FALSE        s095.g.blk.o_female_hypothalamus_lay
    FALSE A2ML1                              8.229733
    FALSE A2ML2                              6.010283
    FALSE A2ML3                             12.817899
    FALSE A2ML4                              5.821444
    FALSE A4GALT                             6.116065
    FALSE A4GNT                              5.459875
    FALSE        s136.d.w.o_female_hypothalamus_lay
    FALSE A2ML1                            7.721650
    FALSE A2ML2                            5.749105
    FALSE A2ML3                           12.574648
    FALSE A2ML4                            5.459875
    FALSE A4GALT                           6.102373
    FALSE A4GNT                            5.459875
    FALSE        s142.o.pk.pu_female_hypothalamus_lay
    FALSE A2ML1                              8.668931
    FALSE A2ML2                              5.918663
    FALSE A2ML3                             12.663561
    FALSE A2ML4                              5.918663
    FALSE A4GALT                             6.406476
    FALSE A4GNT                              5.459875
    FALSE        s176.blk.pu.r_female_hypothalamus_lay
    FALSE A2ML1                               8.193546
    FALSE A2ML2                               5.908769
    FALSE A2ML3                              12.963514
    FALSE A2ML4                               5.459875
    FALSE A4GALT                              6.179145
    FALSE A4GNT                               5.459875
    FALSE        w191.r1_female_hypothalamus_control
    FALSE A2ML1                             8.815980
    FALSE A2ML2                             5.459875
    FALSE A2ML3                            13.552992
    FALSE A2ML4                             5.459875
    FALSE A4GALT                            6.387178
    FALSE A4GNT                             5.459875
    FALSE        x.blu101.w43_female_hypothalamus_inc.d9
    FALSE A2ML1                                 8.367972
    FALSE A2ML2                                 5.459875
    FALSE A2ML3                                12.611396
    FALSE A2ML4                                 5.459875
    FALSE A4GALT                                5.965200
    FALSE A4GNT                                 5.818102
    FALSE        x.blu102.w105_female_hypothalamus_inc.d3
    FALSE A2ML1                                  8.667828
    FALSE A2ML2                                  6.019258
    FALSE A2ML3                                 12.738954
    FALSE A2ML4                                  5.993525
    FALSE A4GALT                                 5.966425
    FALSE A4GNT                                  5.459875
    FALSE        x.blu109.w121_female_hypothalamus_n5
    FALSE A2ML1                              8.057967
    FALSE A2ML2                              5.727809
    FALSE A2ML3                             13.099407
    FALSE A2ML4                              5.959355
    FALSE A4GALT                             6.210307
    FALSE A4GNT                              5.459875
    FALSE        x.blu116.w107_female_hypothalamus_inc.d17
    FALSE A2ML1                                   8.551209
    FALSE A2ML2                                   5.963221
    FALSE A2ML3                                  12.591498
    FALSE A2ML4                                   5.816692
    FALSE A4GALT                                  6.074811
    FALSE A4GNT                                   5.459875
    FALSE        x.blu122.r66_female_hypothalamus_inc.d9
    FALSE A2ML1                                 8.466357
    FALSE A2ML2                                 6.018358
    FALSE A2ML3                                12.228792
    FALSE A2ML4                                 5.893546
    FALSE A4GALT                                5.928002
    FALSE A4GNT                                 5.459875
    FALSE        x.blu43.g132_female_hypothalamus_n9
    FALSE A2ML1                             8.720999
    FALSE A2ML2                             6.085773
    FALSE A2ML3                            12.728607
    FALSE A2ML4                             5.762356
    FALSE A4GALT                            5.952335
    FALSE A4GNT                             5.459875
    FALSE        x.blu6.y80_female_hypothalamus_lay x.g37_female_hypothalamus_n5
    FALSE A2ML1                            8.649894                     8.149124
    FALSE A2ML2                            5.642088                     5.999313
    FALSE A2ML3                           12.887041                    12.737075
    FALSE A2ML4                            5.866243                     5.772525
    FALSE A4GALT                           5.972888                     5.862982
    FALSE A4GNT                            5.459875                     5.459875
    FALSE        x.g4.w50_female_hypothalamus_n9 x.g43_female_hypothalamus_n5
    FALSE A2ML1                         8.082881                     8.655364
    FALSE A2ML2                         5.791399                     5.882787
    FALSE A2ML3                        12.724943                    12.839946
    FALSE A2ML4                         5.459875                     5.882787
    FALSE A4GALT                        6.157866                     5.704625
    FALSE A4GNT                         5.459875                     5.459875
    FALSE        x.g49_female_hypothalamus_n5 x.g9.o166_female_hypothalamus_inc.d9
    FALSE A2ML1                      7.774234                             8.663149
    FALSE A2ML2                      5.779302                             5.667433
    FALSE A2ML3                     12.756730                            12.515485
    FALSE A2ML4                      5.459875                             5.753154
    FALSE A4GALT                     6.132752                             6.304292
    FALSE A4GNT                      5.459875                             6.111251
    FALSE        x.o159.w90_female_hypothalamus_inc.d17
    FALSE A2ML1                                8.656477
    FALSE A2ML2                                5.811353
    FALSE A2ML3                               12.533260
    FALSE A2ML4                                5.889819
    FALSE A4GALT                               6.135990
    FALSE A4GNT                                5.459875
    FALSE        x.o175.g21_female_hypothalamus_n5
    FALSE A2ML1                           8.883273
    FALSE A2ML2                           6.000039
    FALSE A2ML3                          12.978024
    FALSE A2ML4                           5.459875
    FALSE A4GALT                          5.802697
    FALSE A4GNT                           5.459875
    FALSE        x.o37.blu50_female_hypothalamus_hatch.NYNO
    FALSE A2ML1                                    7.968115
    FALSE A2ML2                                    5.903259
    FALSE A2ML3                                   12.779596
    FALSE A2ML4                                    5.803857
    FALSE A4GALT                                   5.856759
    FALSE A4GNT                                    5.459875
    FALSE        x.o70_female_hypothalamus_n5.NYNO
    FALSE A2ML1                           8.102717
    FALSE A2ML2                           5.952638
    FALSE A2ML3                          13.089678
    FALSE A2ML4                           5.681101
    FALSE A4GALT                          6.081393
    FALSE A4GNT                           5.459875
    FALSE        x.r33.w183_female_hypothalamus_inc.d3
    FALSE A2ML1                               8.443501
    FALSE A2ML2                               5.994664
    FALSE A2ML3                              12.378111
    FALSE A2ML4                               5.662998
    FALSE A4GALT                              6.157284
    FALSE A4GNT                               5.459875
    FALSE        x.r39.g10_female_hypothalamus_bldg
    FALSE A2ML1                            8.797592
    FALSE A2ML2                            6.114916
    FALSE A2ML3                           12.952328
    FALSE A2ML4                            5.812148
    FALSE A4GALT                           6.014825
    FALSE A4GNT                            5.890789
    FALSE        x.r44.w95_female_hypothalamus_hatch
    FALSE A2ML1                             8.148040
    FALSE A2ML2                             5.933783
    FALSE A2ML3                            12.834130
    FALSE A2ML4                             6.014658
    FALSE A4GALT                            6.273574
    FALSE A4GNT                             5.459875
    FALSE        x.r48.y139_female_hypothalamus_inc.d17
    FALSE A2ML1                                6.922281
    FALSE A2ML2                                5.710404
    FALSE A2ML3                               12.830305
    FALSE A2ML4                                5.766518
    FALSE A4GALT                               6.094281
    FALSE A4GNT                                5.710404
    FALSE        x.r50.w97_female_hypothalamus_n5 x.w178_female_hypothalamus_n9
    FALSE A2ML1                          8.052164                      7.469400
    FALSE A2ML2                          5.783105                      5.868507
    FALSE A2ML3                         12.885750                     12.937074
    FALSE A2ML4                          5.916045                      5.696324
    FALSE A4GALT                         5.916045                      5.793892
    FALSE A4GNT                          5.621743                      5.459875
    FALSE        x.w51_female_hypothalamus_lay x.w6_female_hypothalamus_n9
    FALSE A2ML1                       8.015834                    8.076602
    FALSE A2ML2                       5.869947                    5.859408
    FALSE A2ML3                      12.931066                   13.056756
    FALSE A2ML4                       5.719750                    5.691036
    FALSE A4GALT                      5.977546                    5.920730
    FALSE A4GNT                       5.719750                    5.459875
    FALSE        x.y109_female_hypothalamus_inc.d9
    FALSE A2ML1                           8.683584
    FALSE A2ML2                           6.046741
    FALSE A2ML3                          12.754908
    FALSE A2ML4                           5.915702
    FALSE A4GALT                          6.199220
    FALSE A4GNT                           5.459875
    FALSE        x.y138.w176_female_hypothalamus_n9 x.y90_female_hypothalamus_hatch
    FALSE A2ML1                            8.689669                        7.888223
    FALSE A2ML2                            6.109853                        5.967751
    FALSE A2ML3                           13.259878                       12.983362
    FALSE A2ML4                            6.109853                        5.459875
    FALSE A4GALT                           6.109853                        6.256877
    FALSE A4GNT                            5.801395                        5.459875
    FALSE        x.y93.g126_female_hypothalamus_inc.d9 x.y9_female_hypothalamus_n9
    FALSE A2ML1                               8.411864                    8.324173
    FALSE A2ML2                               5.800123                    5.711991
    FALSE A2ML3                              12.760428                   12.532008
    FALSE A2ML4                               5.459875                    5.711991
    FALSE A4GALT                              5.876114                    6.450073
    FALSE A4GNT                               5.459875                    5.459875
    FALSE        y.s156.o.r_female_hypothalamus_lay
    FALSE A2ML1                            7.806145
    FALSE A2ML2                            5.893824
    FALSE A2ML3                           12.592523
    FALSE A2ML4                            5.767299
    FALSE A4GALT                           6.007559
    FALSE A4GNT                            5.698185
    FALSE        y126.w92.x_female_hypothalamus_inc.d17
    FALSE A2ML1                                8.032921
    FALSE A2ML2                                6.128543
    FALSE A2ML3                               12.730777
    FALSE A2ML4                                5.796446
    FALSE A4GALT                               6.128543
    FALSE A4GNT                                5.848222
    FALSE        y13.x_female_hypothalamus_inc.d3
    FALSE A2ML1                          8.453839
    FALSE A2ML2                          5.998514
    FALSE A2ML3                         12.307559
    FALSE A2ML4                          5.942202
    FALSE A4GALT                         5.942202
    FALSE A4GNT                          5.459875
    FALSE        y130.o170.x_female_hypothalamus_inc.d17
    FALSE A2ML1                                 8.660381
    FALSE A2ML2                                 5.824421
    FALSE A2ML3                                12.746473
    FALSE A2ML4                                 5.717990
    FALSE A4GALT                                6.160782
    FALSE A4GNT                                 5.459875
    FALSE        y135.blu107.x_female_hypothalamus_inc.d17
    FALSE A2ML1                                   7.768251
    FALSE A2ML2                                   5.826855
    FALSE A2ML3                                  12.523206
    FALSE A2ML4                                   5.643736
    FALSE A4GALT                                  5.908732
    FALSE A4GNT                                   5.459875
    FALSE        y136.x_female_hypothalamus_inc.d17
    FALSE A2ML1                            8.827182
    FALSE A2ML2                            5.989519
    FALSE A2ML3                           12.791736
    FALSE A2ML4                            5.785348
    FALSE A4GALT                           5.879470
    FALSE A4GNT                            5.459875
    FALSE        y140.w119.x_female_hypothalamus_inc.d9
    FALSE A2ML1                                8.754093
    FALSE A2ML2                                5.859415
    FALSE A2ML3                               12.078149
    FALSE A2ML4                                5.948437
    FALSE A4GALT                               5.906221
    FALSE A4GNT                                5.459875
    FALSE        y15.x_female_hypothalamus_hatch y6.o54_female_hypothalamus_n5
    FALSE A2ML1                         6.958676                      8.732890
    FALSE A2ML2                         5.649779                      5.645640
    FALSE A2ML3                        12.565456                     12.734341
    FALSE A2ML4                         5.923382                      5.982795
    FALSE A4GALT                        5.923382                      6.043730
    FALSE A4GNT                         5.459875                      5.645640
    FALSE        y7.g58_female_hypothalamus_hatch
    FALSE A2ML1                          8.500365
    FALSE A2ML2                          5.978181
    FALSE A2ML3                         12.851339
    FALSE A2ML4                          5.827353
    FALSE A4GALT                         6.189034
    FALSE A4GNT                          5.459875
    FALSE        y94.g133.x_female_hypothalamus_n5.NYNO y97.x_female_hypothalamus_n9
    FALSE A2ML1                                8.601284                     8.346067
    FALSE A2ML2                                5.935693                     5.459875
    FALSE A2ML3                               12.881522                    12.922851
    FALSE A2ML4                                6.095984                     5.875976
    FALSE A4GALT                               6.129802                     5.459875
    FALSE A4GNT                                5.459875                     5.968621
    FALSE 'data.frame': 5683 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13532 levels "A2ML1","A2ML2",..: 12078 12217 1174 244 7171 9839 3908 3528 3901 4739 ...
    FALSE  $ padj     : num  7.79e-02 1.52e-02 8.61e-07 1.46e-02 3.45e-02 ...
    FALSE  $ logpadj  : num  1.11 1.82 6.07 1.84 1.46 ...
    FALSE  $ lfc      : num  2.52 2.38 2.34 2.24 2.15 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "control","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE      gene         padj  logpadj      lfc           sextissue direction
    FALSE 1 TIMM10B 7.791452e-02 1.108382 2.516672 female_hypothalamus      bldg
    FALSE 2 TMEM173 1.516419e-02 1.819181 2.380516 female_hypothalamus      bldg
    FALSE 3    BRS3 8.606736e-07 6.065161 2.341286 female_hypothalamus      bldg
    FALSE 4  ADGRF4 1.458942e-02 1.835962 2.242221 female_hypothalamus      bldg
    FALSE 5   LYPD2 3.447821e-02 1.462455 2.145038 female_hypothalamus      bldg
    FALSE 6   RAPSN 3.196418e-02 1.495336 2.100362 female_hypothalamus      bldg
    FALSE 'data.frame': 1 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13532 levels "A2ML1","A2ML2",..: 4742
    FALSE  $ padj     : num 0.0193
    FALSE  $ logpadj  : num 1.72
    FALSE  $ lfc      : num -1.37
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1
    FALSE  $ direction: Factor w/ 3 levels "bldg","NS","lay": 1
    FALSE NULL
    FALSE    gene       padj  logpadj       lfc           sextissue direction
    FALSE 1 HEMGN 0.01926955 1.715128 -1.373336 female_hypothalamus      bldg
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13532 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 
    FALSE  $ direction: Factor w/ 3 levels "lay","NS","inc.d3": 
    FALSE NULL
    FALSE [1] gene      padj      logpadj   lfc       sextissue direction
    FALSE <0 rows> (or 0-length row.names)
    FALSE 'data.frame': 1 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13532 levels "A2ML1","A2ML2",..: 6245
    FALSE  $ padj     : num 9.56e-16
    FALSE  $ logpadj  : num 15
    FALSE  $ lfc      : num -17.2
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1
    FALSE  $ direction: Factor w/ 3 levels "inc.d3","NS",..: 1
    FALSE NULL
    FALSE           gene         padj  logpadj       lfc           sextissue
    FALSE 1 LOC107049632 9.560804e-16 15.01951 -17.18283 female_hypothalamus
    FALSE   direction
    FALSE 1    inc.d3
    FALSE 'data.frame': 5 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13532 levels "A2ML1","A2ML2",..: 6245 5058 7507 4352 2045
    FALSE  $ padj     : num  9.74e-19 2.45e-02 4.54e-02 4.54e-02 4.54e-02
    FALSE  $ logpadj  : num  18.01 1.61 1.34 1.34 1.34
    FALSE  $ lfc      : num  17.84 4.201 0.526 -0.512 -0.708
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1 1 1 1 1
    FALSE  $ direction: Factor w/ 3 levels "inc.d9","NS",..: 3 3 3 1 1
    FALSE NULL
    FALSE           gene         padj   logpadj        lfc           sextissue
    FALSE 1 LOC107049632 9.744986e-19 18.011219 17.8400939 female_hypothalamus
    FALSE 2        IGLL1 2.449263e-02  1.610965  4.2012391 female_hypothalamus
    FALSE 3        MFSD5 4.538999e-02  1.343040  0.5260726 female_hypothalamus
    FALSE 4         GMNN 4.538999e-02  1.343040 -0.5122935 female_hypothalamus
    FALSE 5       CFAP44 4.538999e-02  1.343040 -0.7084001 female_hypothalamus
    FALSE   direction
    FALSE 1   inc.d17
    FALSE 2   inc.d17
    FALSE 3   inc.d17
    FALSE 4    inc.d9
    FALSE 5    inc.d9
    FALSE 'data.frame': 3 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13532 levels "A2ML1","A2ML2",..: 37 1594 9571
    FALSE  $ padj     : num  2.41e-02 1.95e-05 6.50e-03
    FALSE  $ logpadj  : num  1.62 4.71 2.19
    FALSE  $ lfc      : num  4.44 2.19 1.06
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1 1 1
    FALSE  $ direction: Factor w/ 3 levels "inc.d17","NS",..: 3 3 3
    FALSE NULL
    FALSE     gene         padj  logpadj      lfc           sextissue direction
    FALSE 1  ABCB5 2.414266e-02 1.617215 4.444561 female_hypothalamus     hatch
    FALSE 2  CAPN2 1.947616e-05 4.710497 2.192592 female_hypothalamus     hatch
    FALSE 3 PSMD14 6.495688e-03 2.187375 1.059568 female_hypothalamus     hatch
    FALSE 'data.frame': 1927 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13532 levels "A2ML1","A2ML2",..: 5057 6213 4983 6945 13052 9001 5776 5838 11185 12029 ...
    FALSE  $ padj     : num  0.0466 0.051 0.0332 0.0612 0.0208 ...
    FALSE  $ logpadj  : num  1.33 1.29 1.48 1.21 1.68 ...
    FALSE  $ lfc      : num  2.47 2.06 1.66 1.54 1.48 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","n5": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene       padj  logpadj      lfc           sextissue direction
    FALSE 1          IGJ 0.04664654 1.331181 2.467559 female_hypothalamus        n5
    FALSE 2 LOC107049139 0.05099722 1.292453 2.059837 female_hypothalamus        n5
    FALSE 3         IAPP 0.03321504 1.478665 1.657833 female_hypothalamus        n5
    FALSE 4    LOC769174 0.06118946 1.213323 1.542581 female_hypothalamus        n5
    FALSE 5        VSIG1 0.02077613 1.682435 1.479388 female_hypothalamus        n5
    FALSE 6         PI15 0.02763121 1.558600 1.467538 female_hypothalamus        n5
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13532 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "female_hypothalamus": 
    FALSE  $ direction: Factor w/ 3 levels "n5","NS","n9": 
    FALSE NULL
    FALSE [1] gene      padj      logpadj   lfc       sextissue direction
    FALSE <0 rows> (or 0-length row.names)
    FALSE class: DESeqDataSet 
    FALSE dim: 13390 96 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13390): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(96): L.G118_female_pituitary_control.NYNO
    FALSE   R.G106_female_pituitary_control ...
    FALSE   y94.g133.x_female_pituitary_n5 y97.x_female_pituitary_n9
    FALSE colData names(8): V1 bird ... study sextissue
    FALSE [1] 13390    96
    FALSE        L.G118_female_pituitary_control.NYNO
    FALSE A2ML1                              6.955325
    FALSE A2ML2                              5.343272
    FALSE A2ML3                              7.983862
    FALSE A2ML4                              5.753528
    FALSE A4GALT                             7.058252
    FALSE AAAS                               7.338964
    FALSE        R.G106_female_pituitary_control R.R20_female_pituitary_control
    FALSE A2ML1                         7.503694                       7.193169
    FALSE A2ML2                         5.343272                       5.343272
    FALSE A2ML3                         7.555459                       8.049306
    FALSE A2ML4                         5.805696                       5.677082
    FALSE A4GALT                        7.314495                       7.238947
    FALSE AAAS                          7.485973                       7.590293
    FALSE        R.R9_female_pituitary_control.NYNO
    FALSE A2ML1                            8.018062
    FALSE A2ML2                            5.817463
    FALSE A2ML3                           11.818065
    FALSE A2ML4                            5.731024
    FALSE A4GALT                           6.593670
    FALSE AAAS                             7.063703
    FALSE        R.W44_female_pituitary_control.NYNO
    FALSE A2ML1                             7.074247
    FALSE A2ML2                             5.343272
    FALSE A2ML3                             7.677627
    FALSE A2ML4                             6.225665
    FALSE A4GALT                            7.396974
    FALSE AAAS                              7.145866
    FALSE        blk.s061.pu.y_female_pituitary_inc.d9 blk11.x_female_pituitary_bldg
    FALSE A2ML1                               7.385607                      7.717752
    FALSE A2ML2                               6.489640                      5.747120
    FALSE A2ML3                               7.632873                      7.543489
    FALSE A2ML4                               5.797157                      6.038302
    FALSE A4GALT                              7.335099                      7.066071
    FALSE AAAS                                7.503830                      7.539149
    FALSE        blk21.x_female_pituitary_hatch blk4.x_female_pituitary_n9
    FALSE A2ML1                        7.129010                   7.077856
    FALSE A2ML2                        5.722566                   5.778584
    FALSE A2ML3                        8.028719                   7.415263
    FALSE A2ML4                        5.343272                   5.875421
    FALSE A4GALT                       6.691326                   7.231669
    FALSE AAAS                         7.715825                   7.699031
    FALSE        blu.o.x.ATLAS_female_pituitary_control
    FALSE A2ML1                                7.447186
    FALSE A2ML2                                5.911211
    FALSE A2ML3                                7.243884
    FALSE A2ML4                                5.835906
    FALSE A4GALT                               6.918132
    FALSE AAAS                                 6.611747
    FALSE        blu103.x_female_pituitary_hatch.NYNO
    FALSE A2ML1                              7.383351
    FALSE A2ML2                              5.536596
    FALSE A2ML3                              7.510928
    FALSE A2ML4                              5.887256
    FALSE A4GALT                             6.965172
    FALSE AAAS                               7.188008
    FALSE        blu124.w180.x_female_pituitary_hatch blu36.w16_female_pituitary_n9
    FALSE A2ML1                              7.223258                      7.056692
    FALSE A2ML2                              5.946667                      5.676864
    FALSE A2ML3                              7.826865                      7.959956
    FALSE A2ML4                              5.946667                      5.783843
    FALSE A4GALT                             6.744949                      7.324639
    FALSE AAAS                               7.636672                      7.532291
    FALSE        blu38.g135.x_female_pituitary_bldg
    FALSE A2ML1                            7.002616
    FALSE A2ML2                            5.770936
    FALSE A2ML3                            7.812426
    FALSE A2ML4                            5.945907
    FALSE A4GALT                           7.997851
    FALSE AAAS                             7.836059
    FALSE        blu39.o26.x_female_pituitary_inc.d3.NYNO
    FALSE A2ML1                                  7.136107
    FALSE A2ML2                                  5.853816
    FALSE A2ML3                                  7.357372
    FALSE A2ML4                                  6.295759
    FALSE A4GALT                                 7.333067
    FALSE AAAS                                   7.821385
    FALSE        blu47.y96.x_female_pituitary_inc.d9 blu55.g51_female_pituitary_n5
    FALSE A2ML1                             7.131080                      6.411637
    FALSE A2ML2                             5.987312                      5.688070
    FALSE A2ML3                             8.142842                      7.543189
    FALSE A2ML4                             5.752615                      5.343272
    FALSE A4GALT                            6.952044                      7.316396
    FALSE AAAS                              7.736735                      7.768764
    FALSE        g.blk.s004.pk_female_pituitary_lay
    FALSE A2ML1                            7.146972
    FALSE A2ML2                            5.885118
    FALSE A2ML3                            7.492003
    FALSE A2ML4                            5.776432
    FALSE A4GALT                           7.380609
    FALSE AAAS                             7.619066
    FALSE        g.s078.blk.o_female_pituitary_lay
    FALSE A2ML1                           7.403130
    FALSE A2ML2                           5.860494
    FALSE A2ML3                           8.594906
    FALSE A2ML4                           5.656228
    FALSE A4GALT                          7.020951
    FALSE AAAS                            7.474813
    FALSE        g141.blu27.x_female_pituitary_bldg
    FALSE A2ML1                            7.129724
    FALSE A2ML2                            5.816841
    FALSE A2ML3                            7.470767
    FALSE A2ML4                            5.958827
    FALSE A4GALT                           7.300950
    FALSE AAAS                             7.579732
    FALSE        g142.r40.x_female_pituitary_inc.d17
    FALSE A2ML1                             7.042172
    FALSE A2ML2                             5.343272
    FALSE A2ML3                             7.273174
    FALSE A2ML4                             5.688622
    FALSE A4GALT                            6.923310
    FALSE AAAS                              7.787831
    FALSE        g6.w197.x_female_pituitary_inc.d3 g75.x_female_pituitary_inc.d9
    FALSE A2ML1                           6.835635                      7.331524
    FALSE A2ML2                           5.718315                      5.827730
    FALSE A2ML3                           7.735114                      7.822747
    FALSE A2ML4                           5.753881                      5.981937
    FALSE A4GALT                          7.196265                      7.033349
    FALSE AAAS                            7.659344                      7.358852
    FALSE        l.s120.y.blk_female_pituitary_bldg
    FALSE A2ML1                            7.477335
    FALSE A2ML2                            5.864523
    FALSE A2ML3                            7.853079
    FALSE A2ML4                            5.712848
    FALSE A4GALT                           7.306367
    FALSE AAAS                             7.806291
    FALSE        o156.w80.x_female_pituitary_inc.d3
    FALSE A2ML1                            7.090146
    FALSE A2ML2                            5.685562
    FALSE A2ML3                            7.911717
    FALSE A2ML4                            5.639877
    FALSE A4GALT                           7.198558
    FALSE AAAS                             7.461973
    FALSE        o165.w122.x_female_pituitary_inc.d3.NYNO
    FALSE A2ML1                                  7.546115
    FALSE A2ML2                                  6.009401
    FALSE A2ML3                                  7.854682
    FALSE A2ML4                                  5.940117
    FALSE A4GALT                                 7.312135
    FALSE AAAS                                   7.306419
    FALSE        o172.w115.x_female_pituitary_hatch.NYNO
    FALSE A2ML1                                 7.270025
    FALSE A2ML2                                 6.094588
    FALSE A2ML3                                 7.648183
    FALSE A2ML4                                 5.596226
    FALSE A4GALT                                6.595100
    FALSE AAAS                                  7.499366
    FALSE        o173.w179.x_female_pituitary_inc.d3
    FALSE A2ML1                             7.410608
    FALSE A2ML2                             5.962222
    FALSE A2ML3                             7.171067
    FALSE A2ML4                             5.736520
    FALSE A4GALT                            7.450904
    FALSE AAAS                              7.600668
    FALSE        o35.r51.x_female_pituitary_inc.d17
    FALSE A2ML1                            7.088127
    FALSE A2ML2                            5.343272
    FALSE A2ML3                            7.722083
    FALSE A2ML4                            5.849769
    FALSE A4GALT                           7.227015
    FALSE AAAS                             7.629262
    FALSE        o38.blu29.x_female_pituitary_bldg
    FALSE A2ML1                           7.245636
    FALSE A2ML2                           5.872310
    FALSE A2ML3                           7.981304
    FALSE A2ML4                           5.718400
    FALSE A4GALT                          7.245636
    FALSE AAAS                            7.714065
    FALSE        o52.blu53_female_pituitary_inc.d17 o73.x_female_pituitary_inc.d9
    FALSE A2ML1                            6.691821                      7.252249
    FALSE A2ML2                            5.740148                      6.037222
    FALSE A2ML3                            7.816141                      7.988090
    FALSE A2ML4                            5.902798                      5.607726
    FALSE A4GALT                           7.146688                      7.409268
    FALSE AAAS                             7.502561                      7.396897
    FALSE        r.r.x.ATLAS.R2XR_female_pituitary_control
    FALSE A2ML1                                   7.106915
    FALSE A2ML2                                   6.257406
    FALSE A2ML3                                   7.629746
    FALSE A2ML4                                   5.771969
    FALSE A4GALT                                  6.879591
    FALSE AAAS                                    7.035893
    FALSE        r.r.x.ATLAS_female_pituitary_control
    FALSE A2ML1                              7.492996
    FALSE A2ML2                              5.956878
    FALSE A2ML3                              7.461160
    FALSE A2ML4                              5.343272
    FALSE A4GALT                             6.470968
    FALSE AAAS                               7.132025
    FALSE        r.s056.g.o_female_pituitary_bldg r.s171.l.w_female_pituitary_n9
    FALSE A2ML1                          7.085584                       7.020751
    FALSE A2ML2                          5.701752                       5.710107
    FALSE A2ML3                          7.978873                       8.001020
    FALSE A2ML4                          6.174597                       5.603012
    FALSE A4GALT                         6.706911                       7.282734
    FALSE AAAS                           7.505554                       7.482620
    FALSE        r183.o22_female_pituitary_hatch
    FALSE A2ML1                         7.362933
    FALSE A2ML2                         5.871152
    FALSE A2ML3                         7.910473
    FALSE A2ML4                         5.801061
    FALSE A4GALT                        7.059808
    FALSE AAAS                          8.070364
    FALSE        r27.w111.blu125_female_pituitary_inc.d3
    FALSE A2ML1                                 6.946516
    FALSE A2ML2                                 5.343272
    FALSE A2ML3                                 7.205188
    FALSE A2ML4                                 5.816414
    FALSE A4GALT                                6.100879
    FALSE AAAS                                  8.123200
    FALSE        r30.w112.r46_female_pituitary_inc.d9
    FALSE A2ML1                              6.980648
    FALSE A2ML2                              5.743433
    FALSE A2ML3                              7.829375
    FALSE A2ML4                              5.574798
    FALSE A4GALT                             6.987303
    FALSE AAAS                               7.621958
    FALSE        r36.w184.x_female_pituitary_inc.d9
    FALSE A2ML1                            7.193105
    FALSE A2ML2                            5.681920
    FALSE A2ML3                            7.824820
    FALSE A2ML4                            5.343272
    FALSE A4GALT                           6.955236
    FALSE AAAS                             7.329396
    FALSE        r49.w189.x_female_pituitary_inc.d17 r6.x_female_pituitary_control
    FALSE A2ML1                             7.189340                      7.351731
    FALSE A2ML2                             5.517699                      5.343272
    FALSE A2ML3                             7.739432                      8.058421
    FALSE A2ML4                             5.769240                      5.942096
    FALSE A4GALT                            7.004072                      7.210545
    FALSE AAAS                              8.121386                      7.133462
    FALSE        r73.g127.x_female_pituitary_inc.d3 r83.g45_female_pituitary_bldg
    FALSE A2ML1                            7.063381                      6.924343
    FALSE A2ML2                            5.792730                      5.508496
    FALSE A2ML3                            7.508071                      7.853705
    FALSE A2ML4                            5.792730                      5.863223
    FALSE A4GALT                           7.109674                      7.300321
    FALSE AAAS                             7.624246                      7.295085
    FALSE        r95.blu99_female_pituitary_n9 s.o.pk_female_pituitary_lay
    FALSE A2ML1                       7.389499                    7.212726
    FALSE A2ML2                       5.586890                    5.870939
    FALSE A2ML3                       7.400061                    7.634791
    FALSE A2ML4                       5.764238                    5.608206
    FALSE A4GALT                      6.978939                    7.122935
    FALSE AAAS                        7.036714                    7.291948
    FALSE        s.x.ATLAS_female_pituitary_control
    FALSE A2ML1                            7.133346
    FALSE A2ML2                            5.343272
    FALSE A2ML3                            7.133346
    FALSE A2ML4                            5.343272
    FALSE A4GALT                           7.029583
    FALSE AAAS                             7.082581
    FALSE        s063.d.blk.l_female_pituitary_bldg
    FALSE A2ML1                            7.211957
    FALSE A2ML2                            5.727437
    FALSE A2ML3                            7.668112
    FALSE A2ML4                            6.004822
    FALSE A4GALT                           7.369211
    FALSE AAAS                             7.657373
    FALSE        s092.blk.r.o_female_pituitary_bldg
    FALSE A2ML1                            7.443766
    FALSE A2ML2                            6.048913
    FALSE A2ML3                            8.572434
    FALSE A2ML4                            5.844696
    FALSE A4GALT                           7.414498
    FALSE AAAS                             7.514169
    FALSE        s095.g.blk.o_female_pituitary_lay s136.d.w.o_female_pituitary_lay
    FALSE A2ML1                           7.080360                        7.445037
    FALSE A2ML2                           5.709849                        5.665220
    FALSE A2ML3                           7.656097                        7.741977
    FALSE A2ML4                           5.642848                        5.571159
    FALSE A4GALT                          7.296741                        7.134590
    FALSE AAAS                            7.550507                        7.413555
    FALSE        s142.o.pk.pu_female_pituitary_lay
    FALSE A2ML1                           7.264707
    FALSE A2ML2                           5.948244
    FALSE A2ML3                           7.925335
    FALSE A2ML4                           5.772607
    FALSE A4GALT                          7.320833
    FALSE AAAS                            7.448938
    FALSE        s176.blk.pu.r_female_pituitary_lay w191.r1_female_pituitary_control
    FALSE A2ML1                            7.345676                         7.249391
    FALSE A2ML2                            5.814997                         5.343272
    FALSE A2ML3                            7.746464                         7.610658
    FALSE A2ML4                            5.616429                         5.887344
    FALSE A4GALT                           7.119348                         7.416736
    FALSE AAAS                             7.504567                         7.307660
    FALSE        x.blu101.w43_female_pituitary_inc.d9
    FALSE A2ML1                              7.017659
    FALSE A2ML2                              5.980605
    FALSE A2ML3                              7.664285
    FALSE A2ML4                              5.612795
    FALSE A4GALT                             7.052238
    FALSE AAAS                               7.816698
    FALSE        x.blu102.w105_female_pituitary_inc.d3
    FALSE A2ML1                               7.451388
    FALSE A2ML2                               5.664774
    FALSE A2ML3                               8.634484
    FALSE A2ML4                               5.850054
    FALSE A4GALT                              7.301938
    FALSE AAAS                                7.791396
    FALSE        x.blu109.w121_female_pituitary_n5
    FALSE A2ML1                           7.227984
    FALSE A2ML2                           5.601768
    FALSE A2ML3                           7.381126
    FALSE A2ML4                           5.601768
    FALSE A4GALT                          7.254835
    FALSE AAAS                            7.727485
    FALSE        x.blu116.w107_female_pituitary_inc.d17.NYNO
    FALSE A2ML1                                     7.214132
    FALSE A2ML2                                     6.014913
    FALSE A2ML3                                     7.537591
    FALSE A2ML4                                     5.733366
    FALSE A4GALT                                    6.892829
    FALSE AAAS                                      7.081829
    FALSE        x.blu122.r66_female_pituitary_inc.d9
    FALSE A2ML1                              7.271713
    FALSE A2ML2                              5.898478
    FALSE A2ML3                              7.555214
    FALSE A2ML4                              5.343272
    FALSE A4GALT                             6.744316
    FALSE AAAS                               7.574658
    FALSE        x.blu43.g132_female_pituitary_n9 x.blu6.y80_female_pituitary_lay
    FALSE A2ML1                          6.998866                        7.391196
    FALSE A2ML2                          6.038517                        5.826382
    FALSE A2ML3                          8.007939                        8.964289
    FALSE A2ML4                          5.656596                        5.762142
    FALSE A4GALT                         7.022737                        6.956671
    FALSE AAAS                           7.472845                        7.293153
    FALSE        x.g37_female_pituitary_n5 x.g4.w50_female_pituitary_n9
    FALSE A2ML1                   6.799779                     6.584282
    FALSE A2ML2                   5.771381                     5.343272
    FALSE A2ML3                   8.003511                     7.110447
    FALSE A2ML4                   6.016533                     5.849149
    FALSE A4GALT                  7.184744                     7.321749
    FALSE AAAS                    7.313399                     8.020980
    FALSE        x.g43_female_pituitary_n5 x.g49_female_pituitary_n5.NYNO
    FALSE A2ML1                   7.173898                       6.899244
    FALSE A2ML2                   5.922520                       5.804896
    FALSE A2ML3                   7.917994                       7.569002
    FALSE A2ML4                   5.634351                       5.848529
    FALSE A4GALT                  7.209362                       7.244587
    FALSE AAAS                    7.490317                       7.428193
    FALSE        x.g9.o166_female_pituitary_inc.d9.NYNO
    FALSE A2ML1                                7.037748
    FALSE A2ML2                                5.703464
    FALSE A2ML3                                8.025444
    FALSE A2ML4                                5.891599
    FALSE A4GALT                               7.210114
    FALSE AAAS                                 7.164623
    FALSE        x.o159.w90_female_pituitary_inc.d17 x.o175.g21_female_pituitary_n5
    FALSE A2ML1                             7.542794                       7.177088
    FALSE A2ML2                             5.973463                       5.785857
    FALSE A2ML3                             8.043441                       7.531404
    FALSE A2ML4                             5.890103                       5.343272
    FALSE A4GALT                            6.646045                       7.138613
    FALSE AAAS                              7.353947                       7.769600
    FALSE        x.o37.blu50_female_pituitary_hatch x.o70_female_pituitary_n5
    FALSE A2ML1                            6.943858                  7.118585
    FALSE A2ML2                            5.920432                  5.795873
    FALSE A2ML3                            7.735014                  7.938951
    FALSE A2ML4                            5.860179                  5.826842
    FALSE A4GALT                           6.908700                  7.157159
    FALSE AAAS                             7.649026                  7.266131
    FALSE        x.r33.w183_female_pituitary_inc.d3 x.r39.g10_female_pituitary_bldg
    FALSE A2ML1                            7.298289                        7.022068
    FALSE A2ML2                            5.787590                        5.842855
    FALSE A2ML3                            7.602283                        7.503927
    FALSE A2ML4                            5.787590                        6.108047
    FALSE A4GALT                           7.045581                        7.368544
    FALSE AAAS                             7.815308                        7.593364
    FALSE        x.r44.w95_female_pituitary_hatch
    FALSE A2ML1                          7.101910
    FALSE A2ML2                          5.801597
    FALSE A2ML3                          7.469674
    FALSE A2ML4                          5.668035
    FALSE A4GALT                         6.950020
    FALSE AAAS                           8.026808
    FALSE        x.r48.y139_female_pituitary_inc.d17.NYNO
    FALSE A2ML1                                  7.049004
    FALSE A2ML2                                  5.641808
    FALSE A2ML3                                  7.980773
    FALSE A2ML4                                  5.764719
    FALSE A4GALT                                 6.417158
    FALSE AAAS                                   7.612968
    FALSE        x.r50.w97_female_pituitary_n5 x.w178_female_pituitary_n9
    FALSE A2ML1                       7.011593                   7.178227
    FALSE A2ML2                       5.608656                   5.878963
    FALSE A2ML3                       9.105131                   7.688267
    FALSE A2ML4                       5.496635                   5.612269
    FALSE A4GALT                      6.971015                   7.066128
    FALSE AAAS                        7.522061                   7.524012
    FALSE        x.w51_female_pituitary_lay x.w6_female_pituitary_n9
    FALSE A2ML1                    7.422162                 7.093516
    FALSE A2ML2                    5.846840                 5.701891
    FALSE A2ML3                    7.350968                 7.622085
    FALSE A2ML4                    5.923775                 6.009967
    FALSE A4GALT                   7.159302                 6.978611
    FALSE AAAS                     7.942733                 7.683322
    FALSE        x.y109_female_pituitary_inc.d9 x.y138.w176_female_pituitary_n9
    FALSE A2ML1                        7.027414                        7.256580
    FALSE A2ML2                        5.631774                        5.658120
    FALSE A2ML3                        7.528711                        8.338060
    FALSE A2ML4                        5.343272                        6.107070
    FALSE A4GALT                       6.923208                        6.956741
    FALSE AAAS                         7.460077                        7.687681
    FALSE        x.y90_female_pituitary_hatch x.y93.g126_female_pituitary_inc.d9
    FALSE A2ML1                      7.380686                           7.376013
    FALSE A2ML2                      5.717024                           5.749078
    FALSE A2ML3                      7.949162                           7.664902
    FALSE A2ML4                      5.785012                           5.749078
    FALSE A4GALT                     7.020994                           7.304352
    FALSE AAAS                       7.788852                           7.731073
    FALSE        x.y9_female_pituitary_n9 y.s156.o.r_female_pituitary_lay
    FALSE A2ML1                  7.346735                        6.984656
    FALSE A2ML2                  5.726664                        5.822938
    FALSE A2ML3                  7.616354                        7.740720
    FALSE A2ML4                  5.586177                        5.777507
    FALSE A4GALT                 7.346735                        7.203783
    FALSE AAAS                   7.722739                        7.794374
    FALSE        y126.w92.x_female_pituitary_inc.d17
    FALSE A2ML1                             7.554412
    FALSE A2ML2                             6.097999
    FALSE A2ML3                             7.803915
    FALSE A2ML4                             5.879951
    FALSE A4GALT                            7.206529
    FALSE AAAS                              7.374845
    FALSE        y128.g23.x_female_pituitary_inc.d9 y13.x_female_pituitary_inc.d3
    FALSE A2ML1                            6.972728                      7.341076
    FALSE A2ML2                            5.852669                      5.343272
    FALSE A2ML3                            7.668880                      7.929950
    FALSE A2ML4                            5.852669                      5.746803
    FALSE A4GALT                           7.168833                      6.865553
    FALSE AAAS                             7.853756                      7.563499
    FALSE        y130.o170.x_female_pituitary_inc.d17
    FALSE A2ML1                              7.436032
    FALSE A2ML2                              5.923491
    FALSE A2ML3                              7.493414
    FALSE A2ML4                              5.990906
    FALSE A4GALT                             7.521183
    FALSE AAAS                               7.329025
    FALSE        y135.blu107.x_female_pituitary_inc.d17.NYNO
    FALSE A2ML1                                     7.301138
    FALSE A2ML2                                     5.737583
    FALSE A2ML3                                     7.877447
    FALSE A2ML4                                     5.540887
    FALSE A4GALT                                    6.928717
    FALSE AAAS                                      7.468971
    FALSE        y136.x_female_pituitary_inc.d17 y140.w119.x_female_pituitary_inc.d9
    FALSE A2ML1                         7.198867                            7.199883
    FALSE A2ML2                         5.623717                            5.522355
    FALSE A2ML3                         7.515573                            6.819538
    FALSE A2ML4                         6.024943                            5.742694
    FALSE A4GALT                        7.078643                            6.916311
    FALSE AAAS                          7.247003                            7.708942
    FALSE        y15.x_female_pituitary_hatch y6.o54_female_pituitary_n5
    FALSE A2ML1                      6.879339                   7.186986
    FALSE A2ML2                      5.779699                   6.135577
    FALSE A2ML3                      8.001796                   7.817287
    FALSE A2ML4                      5.876779                   5.735111
    FALSE A4GALT                     7.092762                   7.399663
    FALSE AAAS                       7.816935                   7.810631
    FALSE        y7.g58_female_pituitary_hatch y94.g133.x_female_pituitary_n5
    FALSE A2ML1                       7.218568                       7.305126
    FALSE A2ML2                       5.728973                       5.845877
    FALSE A2ML3                       7.609924                       7.443190
    FALSE A2ML4                       5.343272                       5.754333
    FALSE A4GALT                      7.266504                       6.853011
    FALSE AAAS                        7.501983                       7.823681
    FALSE        y97.x_female_pituitary_n9
    FALSE A2ML1                   7.041117
    FALSE A2ML2                   5.926072
    FALSE A2ML3                   7.762960
    FALSE A2ML4                   5.706214
    FALSE A4GALT                  7.278895
    FALSE AAAS                    7.373161
    FALSE 'data.frame': 6129 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13390 levels "A2ML1","A2ML2",..: 11894 6271 867 11937 6475 2946 12208 553 6259 3197 ...
    FALSE  $ padj     : num  1.12e-05 5.91e-04 7.04e-09 2.29e-03 4.31e-07 ...
    FALSE  $ logpadj  : num  4.95 3.23 8.15 2.64 6.37 ...
    FALSE  $ lfc      : num  4.77 4.45 4.18 3.93 3.75 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "control","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene         padj  logpadj      lfc        sextissue direction
    FALSE 1           TH 1.117464e-05 4.951766 4.766259 female_pituitary      bldg
    FALSE 2 LOC107050724 5.911085e-04 3.228333 4.451202 female_pituitary      bldg
    FALSE 3       ATP2B4 7.038760e-09 8.152504 4.183579 female_pituitary      bldg
    FALSE 4      TIMM10B 2.287419e-03 2.640654 3.928345 female_pituitary      bldg
    FALSE 5 LOC107055658 4.305160e-07 6.366011 3.749904 female_pituitary      bldg
    FALSE 6         DMB2 1.276116e-07 6.894110 3.339761 female_pituitary      bldg
    FALSE 'data.frame': 279 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13390 levels "A2ML1","A2ML2",..: 9844 9361 9384 5943 6448 2456 324 2935 11315 12226 ...
    FALSE  $ padj     : num  6.35e-08 7.30e-07 7.99e-03 5.27e-02 2.53e-03 ...
    FALSE  $ logpadj  : num  7.2 6.14 2.1 1.28 2.6 ...
    FALSE  $ lfc      : num  6.6 5.87 4.69 3.86 3.11 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "bldg","NS","lay": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene         padj  logpadj      lfc        sextissue direction
    FALSE 1         REG4 6.349233e-08 7.197279 6.602559 female_pituitary       lay
    FALSE 2       PROCA1 7.295224e-07 6.136961 5.870406 female_pituitary       lay
    FALSE 3        PRPH2 7.992796e-03 2.097301 4.693180 female_pituitary       lay
    FALSE 4 LOC101747844 5.273908e-02 1.277867 3.861362 female_pituitary       lay
    FALSE 5 LOC107054855 2.527663e-03 2.597281 3.106098 female_pituitary       lay
    FALSE 6       CRABP1 1.353903e-02 1.868412 2.629124 female_pituitary       lay
    FALSE 'data.frame': 353 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13390 levels "A2ML1","A2ML2",..: 5682 1182 10654 5774 3846 10781 5573 2538 320 395 ...
    FALSE  $ padj     : num  0.024864 0.000772 0.024864 0.037539 0.018043 ...
    FALSE  $ logpadj  : num  1.6 3.11 1.6 1.43 1.74 ...
    FALSE  $ lfc      : num  2.97 2.34 2.2 1.72 1.66 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "lay","NS","inc.d3": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene         padj  logpadj      lfc        sextissue direction
    FALSE 1       LGALS1 0.0248641664 1.604426 2.970741 female_pituitary    inc.d3
    FALSE 2          BTC 0.0007715944 3.112611 2.342521 female_pituitary    inc.d3
    FALSE 3         SGK2 0.0248641664 1.604426 2.203634 female_pituitary    inc.d3
    FALSE 4 LOC100857380 0.0375385799 1.425522 1.724579 female_pituitary    inc.d3
    FALSE 5        FDX1L 0.0180428659 1.743694 1.663060 female_pituitary    inc.d3
    FALSE 6      SLC12A4 0.0240973893 1.618030 1.653605 female_pituitary    inc.d3
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13390 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 
    FALSE  $ direction: Factor w/ 3 levels "inc.d3","NS",..: 
    FALSE NULL
    FALSE [1] gene      padj      logpadj   lfc       sextissue direction
    FALSE <0 rows> (or 0-length row.names)
    FALSE 'data.frame': 2004 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13390 levels "A2ML1","A2ML2",..: 5563 3981 3876 1868 1896 2152 2836 10705 1754 10787 ...
    FALSE  $ padj     : num  6.61e-23 3.47e-06 3.14e-12 8.86e-08 8.50e-24 ...
    FALSE  $ logpadj  : num  22.18 5.46 11.5 7.05 23.07 ...
    FALSE  $ lfc      : num  4.89 4.67 4.67 4.59 4.5 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "inc.d9","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE    gene         padj   logpadj      lfc        sextissue direction
    FALSE 1 KPNA2 6.613542e-23 22.179566 4.886097 female_pituitary   inc.d17
    FALSE 2 FOXM1 3.466167e-06  5.460151 4.674672 female_pituitary   inc.d17
    FALSE 3  FGF6 3.139670e-12 11.503116 4.672478 female_pituitary   inc.d17
    FALSE 4 CDCA3 8.864060e-08  7.052367 4.590127 female_pituitary   inc.d17
    FALSE 5  CDK1 8.502111e-24 23.070473 4.499808 female_pituitary   inc.d17
    FALSE 6 CKAP2 8.046251e-21 20.094406 4.412120 female_pituitary   inc.d17
    FALSE 'data.frame': 5 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13390 levels "A2ML1","A2ML2",..: 1419 6238 4691 5239 12538
    FALSE  $ padj     : num  0.0077 0.0927 0.0927 0.0077 0.0586
    FALSE  $ logpadj  : num  2.11 1.03 1.03 2.11 1.23
    FALSE  $ lfc      : num  5.057 1.663 1.422 0.293 -0.581
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1
    FALSE  $ direction: Factor w/ 3 levels "inc.d17","NS",..: 3 3 3 3 1
    FALSE NULL
    FALSE           gene        padj  logpadj        lfc        sextissue direction
    FALSE 1    C5H11ORF9 0.007698626 2.113587  5.0566225 female_pituitary     hatch
    FALSE 2 LOC107050516 0.092669478 1.033063  1.6633428 female_pituitary     hatch
    FALSE 3        HEBP2 0.092669478 1.033063  1.4222028 female_pituitary     hatch
    FALSE 4     IVNS1ABP 0.007698626 2.113587  0.2927521 female_pituitary     hatch
    FALSE 5          TTK 0.058573209 1.232301 -0.5811019 female_pituitary   inc.d17
    FALSE 'data.frame': 1062 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13390 levels "A2ML1","A2ML2",..: 1741 6568 5682 8222 5897 6060 11412 10654 12551 625 ...
    FALSE  $ padj     : num  0.07185 0.0607 0.00539 0.00182 0.0191 ...
    FALSE  $ logpadj  : num  1.14 1.22 2.27 2.74 1.72 ...
    FALSE  $ lfc      : num  4.12 3.24 3.12 2.56 2.41 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","n5": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene        padj  logpadj      lfc        sextissue direction
    FALSE 1          CCK 0.071853676 1.143551 4.123401 female_pituitary        n5
    FALSE 2 LOC107057630 0.060695999 1.216840 3.239277 female_pituitary        n5
    FALSE 3       LGALS1 0.005387965 2.268575 3.124433 female_pituitary        n5
    FALSE 4        NPTX1 0.001822944 2.739227 2.560311 female_pituitary        n5
    FALSE 5 LOC100859605 0.019099053 1.718988 2.412841 female_pituitary        n5
    FALSE 6 LOC101750367 0.012534506 1.901893 2.287775 female_pituitary        n5
    FALSE 'data.frame': 22 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13390 levels "A2ML1","A2ML2",..: 10496 12923 9722 2150 6762 10390 4503 1590 8285 12399 ...
    FALSE  $ padj     : num  0.0366 0.0228 0.0228 0.0223 0.0387 ...
    FALSE  $ logpadj  : num  1.44 1.64 1.64 1.65 1.41 ...
    FALSE  $ lfc      : num  2.9 2.34 2.25 1.99 1.47 ...
    FALSE  $ sextissue: Factor w/ 1 level "female_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "n5","NS","n9": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE        gene       padj  logpadj      lfc        sextissue direction
    FALSE 1     SEBOX 0.03664755 1.435955 2.895520 female_pituitary        n9
    FALSE 2       VTN 0.02281199 1.641837 2.340757 female_pituitary        n9
    FALSE 3     RASA4 0.02281199 1.641837 2.250830 female_pituitary        n9
    FALSE 4    CITED4 0.02231446 1.651414 1.990386 female_pituitary        n9
    FALSE 5 LOC424214 0.03867513 1.412568 1.473628 female_pituitary        n9
    FALSE 6     SARM1 0.07569012 1.120961 1.462395 female_pituitary        n9
    FALSE class: DESeqDataSet 
    FALSE dim: 13698 96 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13698): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(96): L.Blu13_male_gonad_control.NYNO
    FALSE   L.G107_male_gonad_control ... y95.g131.x_male_gonad_inc.d9
    FALSE   y98.o50.x_male_gonad_inc.d3
    FALSE colData names(8): V1 bird ... study sextissue
    FALSE [1] 13698    96
    FALSE        L.Blu13_male_gonad_control.NYNO L.G107_male_gonad_control
    FALSE A2ML1                         6.019359                  6.040564
    FALSE A2ML2                         5.167541                  4.916132
    FALSE A2ML3                         7.503053                  7.440497
    FALSE A2ML4                         5.167541                  5.318206
    FALSE A4GALT                        5.403130                  5.560384
    FALSE A4GNT                         4.585769                  4.585769
    FALSE        L.R3_male_gonad_control.NYNO L.R8_male_gonad_control
    FALSE A2ML1                      5.642209                6.336031
    FALSE A2ML2                      5.122765                4.881135
    FALSE A2ML3                      7.402096                7.341107
    FALSE A2ML4                      5.340937                5.303117
    FALSE A4GALT                     5.961968                5.208583
    FALSE A4GNT                      5.122765                4.585769
    FALSE        L.W33_male_gonad_control L.W3_male_gonad_control.NYNO
    FALSE A2ML1                  6.297341                     6.687273
    FALSE A2ML2                  4.585769                     4.585769
    FALSE A2ML3                  7.737439                     7.508276
    FALSE A2ML4                  5.068257                     5.482144
    FALSE A4GALT                 5.741944                     5.482144
    FALSE A4GNT                  4.585769                     4.585769
    FALSE        L.W4_male_gonad_control.NYNO R.Y108.W29_male_gonad_control
    FALSE A2ML1                      6.204505                      5.978584
    FALSE A2ML2                      4.585769                      4.585769
    FALSE A2ML3                      7.278444                      7.296223
    FALSE A2ML4                      5.056558                      5.172669
    FALSE A4GALT                     5.800676                      5.620576
    FALSE A4GNT                      5.056558                      4.958485
    FALSE        blk12.x_male_gonad_n5 blk17.x_male_gonad_inc.d17
    FALSE A2ML1               5.941583                   6.406178
    FALSE A2ML2               4.585769                   4.585769
    FALSE A2ML3               7.253837                   7.048093
    FALSE A2ML4               5.190266                   5.206107
    FALSE A4GALT              5.815764                   5.382638
    FALSE A4GNT               4.936470                   4.585769
    FALSE        blu104.w120.x_male_gonad_hatch blu108.w40.o158_male_gonad_inc.d9
    FALSE A2ML1                        6.286342                          6.206925
    FALSE A2ML2                        4.585769                          5.105855
    FALSE A2ML3                        7.689247                          7.504245
    FALSE A2ML4                        5.192962                          5.105855
    FALSE A4GALT                       5.504528                          5.794841
    FALSE A4GNT                        4.585769                          4.887121
    FALSE        blu111.w113.x_male_gonad_inc.d3 blu113.w124.x_male_gonad_inc.d17
    FALSE A2ML1                         6.276823                         6.146349
    FALSE A2ML2                         5.167583                         4.585769
    FALSE A2ML3                         7.277550                         7.229351
    FALSE A2ML4                         5.061880                         5.328991
    FALSE A4GALT                        5.727286                         5.491192
    FALSE A4GNT                         4.923191                         4.585769
    FALSE        blu114.r38.w198_male_gonad_bldg blu121.w91.x_male_gonad_inc.d17
    FALSE A2ML1                         6.275318                        6.488103
    FALSE A2ML2                         4.585769                        4.946019
    FALSE A2ML3                         7.222322                        7.456877
    FALSE A2ML4                         5.177535                        5.093932
    FALSE A4GALT                        5.177535                        5.206561
    FALSE A4GNT                         5.070061                        4.585769
    FALSE        blu33.y88.x_male_gonad_bldg blu37.r65.x_male_gonad_n5
    FALSE A2ML1                     6.147268                  6.268583
    FALSE A2ML2                     4.994115                  5.197729
    FALSE A2ML3                     7.362804                  7.556385
    FALSE A2ML4                     5.343527                  5.086667
    FALSE A4GALT                    5.487176                  5.511608
    FALSE A4GNT                     4.874994                  4.940842
    FALSE        blu41.y100.x_male_gonad_n5 blu81.r88_male_gonad_n9
    FALSE A2ML1                    6.308279                6.156110
    FALSE A2ML2                    4.585769                4.585769
    FALSE A2ML3                    7.448328                7.386146
    FALSE A2ML4                    5.323438                4.973281
    FALSE A4GALT                   5.408335                5.791690
    FALSE A4GNT                    4.957605                4.831294
    FALSE        d.s008.y.blk_male_gonad_n5 d.s047.blk.o_male_gonad_n5
    FALSE A2ML1                    6.273494                   5.707875
    FALSE A2ML2                    4.585769                   4.777396
    FALSE A2ML3                    7.316192                   7.436762
    FALSE A2ML4                    5.212151                   5.187813
    FALSE A4GALT                   5.654622                   5.793031
    FALSE A4GNT                    4.585769                   4.856572
    FALSE        g.s.blk.d_male_gonad_n9 g.s.blk.y_male_gonad_lay
    FALSE A2ML1                 6.093690                 6.143106
    FALSE A2ML2                 4.848324                 5.040616
    FALSE A2ML3                 7.358803                 7.476730
    FALSE A2ML4                 5.321416                 4.938673
    FALSE A4GALT                5.660601                 5.712377
    FALSE A4GNT                 4.848324                 4.585769
    FALSE        g.s043.pu.blk_male_gonad_lay g104.w82.x_male_gonad_bldg
    FALSE A2ML1                      6.880747                   6.027579
    FALSE A2ML2                      4.585769                   5.274755
    FALSE A2ML3                      6.880747                   6.900669
    FALSE A2ML4                      5.380400                   4.585769
    FALSE A4GALT                     5.696192                   5.874385
    FALSE A4GNT                      4.585769                   4.585769
    FALSE        g114.w83.x_male_gonad_hatch.NYNO g130.y81.x_male_gonad_inc.d17
    FALSE A2ML1                          5.856556                      6.350003
    FALSE A2ML2                          5.027668                      4.585769
    FALSE A2ML3                          7.444107                      7.352624
    FALSE A2ML4                          5.079350                      5.731422
    FALSE A4GALT                         5.559182                      5.441689
    FALSE A4GNT                          4.807366                      4.585769
    FALSE        g143.blu32.x_male_gonad_inc.d17 g146.blu51_male_gonad_inc.d3
    FALSE A2ML1                         6.337309                     6.036937
    FALSE A2ML2                         5.248459                     5.074353
    FALSE A2ML3                         7.461428                     7.555819
    FALSE A2ML4                         5.056403                     4.932072
    FALSE A4GALT                        5.800293                     5.925632
    FALSE A4GNT                         4.585769                     4.585769
    FALSE        g20.w106.x_male_gonad_inc.d3 g52.blu58_male_gonad_bldg
    FALSE A2ML1                      5.886343                  6.189829
    FALSE A2ML2                      4.902169                  4.909937
    FALSE A2ML3                      7.191087                  7.322717
    FALSE A2ML4                      5.413177                  5.144919
    FALSE A4GALT                     5.615271                  5.371730
    FALSE A4GNT                      4.902169                  4.585769
    FALSE        g53.y84_male_gonad_hatch o.s.w.r_male_gonad_lay
    FALSE A2ML1                  6.234055               6.173350
    FALSE A2ML2                  4.867491               5.080907
    FALSE A2ML3                  7.000379               7.191781
    FALSE A2ML4                  5.270495               5.238415
    FALSE A4GALT                 5.781220               5.646170
    FALSE A4GNT                  4.585769               5.015093
    FALSE        o152.o120.w42_male_gonad_n5 o39.y77.x_male_gonad_hatch
    FALSE A2ML1                     6.272271                   6.217981
    FALSE A2ML2                     4.997331                   5.087004
    FALSE A2ML3                     7.541205                   7.479798
    FALSE A2ML4                     5.030058                   5.397649
    FALSE A4GALT                    5.724037                   5.752870
    FALSE A4GNT                     4.585769                   4.995707
    FALSE        o44.blu26.x_male_gonad_hatch o48.r197.x_male_gonad_inc.d3
    FALSE A2ML1                      6.629356                     6.225677
    FALSE A2ML2                      4.585769                     4.962138
    FALSE A2ML3                      7.303527                     7.172820
    FALSE A2ML4                      5.269579                     5.116547
    FALSE A4GALT                     5.372975                     5.747672
    FALSE A4GNT                      4.585769                     4.585769
    FALSE        o49.x_male_gonad_inc.d9 o57.g59_male_gonad_inc.d9
    FALSE A2ML1                 6.388085                  5.746470
    FALSE A2ML2                 4.585769                  4.893028
    FALSE A2ML3                 7.401636                  7.294318
    FALSE A2ML4                 5.024968                  5.267753
    FALSE A4GALT                5.866036                  5.629262
    FALSE A4GNT                 5.122652                  4.893028
    FALSE        pk.s238.blk.w_male_gonad_lay pk.w.s141.o_male_gonad_lay
    FALSE A2ML1                      6.687066                   5.608320
    FALSE A2ML2                      5.091657                   4.585769
    FALSE A2ML3                      7.453791                   7.581869
    FALSE A2ML4                      5.080139                   5.084842
    FALSE A4GALT                     5.630497                   5.569717
    FALSE A4GNT                      4.960257                   4.874864
    FALSE        r.s005.pk.blk_male_gonad_lay r.s059.d.o_male_gonad_bldg
    FALSE A2ML1                      6.284581                   6.159828
    FALSE A2ML2                      4.855011                   4.915710
    FALSE A2ML3                      7.331248                   7.230242
    FALSE A2ML4                      5.362463                   5.154795
    FALSE A4GALT                     5.786430                   5.657624
    FALSE A4GNT                      4.585769                   4.585769
    FALSE        r.s116.blk.pu_male_gonad_lay r.y.s007.blk_male_gonad_n9
    FALSE A2ML1                      6.380730                   6.048368
    FALSE A2ML2                      4.810989                   4.585769
    FALSE A2ML3                      7.369286                   7.769750
    FALSE A2ML4                      5.357512                   4.978578
    FALSE A4GALT                     5.858659                   5.706165
    FALSE A4GNT                      4.810989                   4.585769
    FALSE        r176.blu54_male_gonad_inc.d17 r190.o43.x_male_gonad_lay
    FALSE A2ML1                       6.336703                  5.565330
    FALSE A2ML2                       4.840981                  4.835081
    FALSE A2ML3                       7.400262                  7.375459
    FALSE A2ML4                       5.457671                  5.240571
    FALSE A4GALT                      5.492166                  5.594322
    FALSE A4GNT                       4.840981                  4.835081
    FALSE        r195.x_male_gonad_n9 r37.w100.x_male_gonad_n9
    FALSE A2ML1              6.645658                 6.340805
    FALSE A2ML2              4.867456                 4.996352
    FALSE A2ML3              7.385259                 7.400482
    FALSE A2ML4              5.110695                 5.087791
    FALSE A4GALT             5.880195                 5.651352
    FALSE A4GNT              4.867456                 4.585769
    FALSE        r41.w99.x_male_gonad_hatch r45.X_male_gonad_inc.d9
    FALSE A2ML1                    6.388519                5.877571
    FALSE A2ML2                    4.905506                4.918877
    FALSE A2ML3                    7.431216                7.212326
    FALSE A2ML4                    5.137337                4.918877
    FALSE A4GALT                   5.791261                5.838981
    FALSE A4GNT                    5.037029                5.055819
    FALSE        r72.y83.x_male_gonad_hatch s.pu148.blk.r_male_gonad_bldg
    FALSE A2ML1                    6.181928                      5.962161
    FALSE A2ML2                    4.585769                      4.942160
    FALSE A2ML3                    7.224256                      7.705652
    FALSE A2ML4                    5.251379                      5.293253
    FALSE A4GALT                   5.573649                      5.514945
    FALSE A4GNT                    4.920780                      5.088516
    FALSE        s065.l.d.o_male_gonad_bldg s066.l.d.r_male_gonad_bldg
    FALSE A2ML1                    6.639209                   6.505952
    FALSE A2ML2                    4.925243                   5.275327
    FALSE A2ML3                    7.364071                   7.483752
    FALSE A2ML4                    5.001063                   4.896492
    FALSE A4GALT                   6.109852                   5.795325
    FALSE A4GNT                    4.925243                   4.896492
    FALSE        s150.w.g.blk_male_gonad_lay s187.l.o.r_male_gonad_n9
    FALSE A2ML1                     6.307925                 6.424597
    FALSE A2ML2                     4.999689                 4.793968
    FALSE A2ML3                     7.471373                 7.290589
    FALSE A2ML4                     4.999689                 5.206118
    FALSE A4GALT                    5.556086                 5.716791
    FALSE A4GNT                     4.585769                 4.793968
    FALSE        s243.blk.pk.r_male_gonad_lay w34.x_male_gonad_inc.d9
    FALSE A2ML1                      6.058588                6.204249
    FALSE A2ML2                      4.954091                4.913135
    FALSE A2ML3                      7.449724                7.493275
    FALSE A2ML4                      5.400722                5.236380
    FALSE A4GALT                     5.934364                5.694729
    FALSE A4GNT                      4.585769                4.585769
    FALSE        x.blk.blk.ATLAS_male_gonad_control x.blk16_male_gonad_n9.NYNO
    FALSE A2ML1                            6.676337                   6.321121
    FALSE A2ML2                            4.585769                   4.915860
    FALSE A2ML3                            7.464825                   7.114366
    FALSE A2ML4                            5.362457                   5.155050
    FALSE A4GALT                           5.815730                   5.788424
    FALSE A4GNT                            4.585769                   4.585769
    FALSE        x.blu106.o153_male_gonad_inc.d9.NYNO
    FALSE A2ML1                              6.329461
    FALSE A2ML2                              4.816994
    FALSE A2ML3                              7.093050
    FALSE A2ML4                              5.344695
    FALSE A4GALT                             5.796448
    FALSE A4GNT                              4.585769
    FALSE        x.blu117.w89_male_gonad_inc.d17 x.blu23.w14_male_gonad_n9
    FALSE A2ML1                         5.969851                  6.136715
    FALSE A2ML2                         4.585769                  4.585769
    FALSE A2ML3                         7.510585                  7.616270
    FALSE A2ML4                         5.276574                  5.097495
    FALSE A4GALT                        5.713392                  5.531143
    FALSE A4GNT                         4.585769                  4.585769
    FALSE        x.blu30_male_gonad_n5 x.blu42.o28_male_gonad_inc.d3
    FALSE A2ML1               6.255314                      6.280573
    FALSE A2ML2               4.893868                      4.936726
    FALSE A2ML3               7.265350                      7.424327
    FALSE A2ML4               5.145824                      4.936726
    FALSE A4GALT              5.749492                      5.435232
    FALSE A4GNT               5.055218                      4.936726
    FALSE        x.g.ATLAS_male_gonad_control x.g.g.ATLAS_male_gonad_control
    FALSE A2ML1                      5.998806                       6.421553
    FALSE A2ML2                      5.604085                       5.022850
    FALSE A2ML3                      7.781160                       7.500347
    FALSE A2ML4                      5.422799                       5.148622
    FALSE A4GALT                     5.604085                       5.637077
    FALSE A4GNT                      4.585769                       4.838760
    FALSE        x.g.g.g.ATLAS_male_gonad_control x.g13.w109_male_gonad_inc.d9
    FALSE A2ML1                          6.624649                     6.105341
    FALSE A2ML2                          4.585769                     5.052827
    FALSE A2ML3                          7.453107                     7.826637
    FALSE A2ML4                          5.274040                     5.243467
    FALSE A4GALT                         5.444054                     5.612879
    FALSE A4GNT                          4.585769                     4.585769
    FALSE        x.g14.w199_male_gonad_inc.d17 x.g147.blu28_male_gonad_inc.d3
    FALSE A2ML1                       6.487454                       6.434375
    FALSE A2ML2                       4.585769                       4.862278
    FALSE A2ML3                       7.434877                       7.445798
    FALSE A2ML4                       5.086058                       5.101118
    FALSE A4GALT                      5.289828                       5.469262
    FALSE A4GNT                       4.585769                       4.862278
    FALSE        x.g70_male_gonad_hatch x.o160.w102_male_gonad_hatch
    FALSE A2ML1                6.435718                     6.526544
    FALSE A2ML2                4.929511                     5.083020
    FALSE A2ML3                7.354618                     7.193422
    FALSE A2ML4                5.700586                     5.285593
    FALSE A4GALT               5.916330                     5.676971
    FALSE A4GNT                4.585769                     4.938245
    FALSE        x.o163.w101_male_gonad_inc.d3 x.o164.w123_male_gonad_n5
    FALSE A2ML1                       5.751410                  6.082011
    FALSE A2ML2                       4.585769                  4.772597
    FALSE A2ML3                       7.236910                  7.262636
    FALSE A2ML4                       4.983709                  5.228094
    FALSE A4GALT                      5.619690                  5.790593
    FALSE A4GNT                       4.585769                  4.585769
    FALSE        x.o2_male_gonad_n9 x.o30.g134_male_gonad_bldg
    FALSE A2ML1            6.216685                   6.434445
    FALSE A2ML2            4.817332                   4.748256
    FALSE A2ML3            7.501752                   7.598853
    FALSE A2ML4            5.149985                   5.167974
    FALSE A4GALT           5.378766                   5.706905
    FALSE A4GNT            4.912898                   4.748256
    FALSE        x.o47.y82_male_gonad_inc.d9 x.o68_male_gonad_n5
    FALSE A2ML1                     6.393781            6.100262
    FALSE A2ML2                     4.585769            4.842551
    FALSE A2ML3                     7.424317            7.293290
    FALSE A2ML4                     5.117669            5.097326
    FALSE A4GALT                    5.333848            5.637910
    FALSE A4GNT                     4.962938            4.842551
    FALSE        x.r178_male_gonad_hatch x.r181_male_gonad_n5
    FALSE A2ML1                 6.556255             5.925677
    FALSE A2ML2                 4.976897             4.585769
    FALSE A2ML3                 7.546965             7.410825
    FALSE A2ML4                 5.341871             4.585769
    FALSE A4GALT                5.791055             5.729178
    FALSE A4GNT                 4.935817             4.585769
    FALSE        x.r29.w96_male_gonad_inc.d17 x.r64.g140_male_gonad_inc.d3
    FALSE A2ML1                      5.681985                     6.231683
    FALSE A2ML2                      4.585769                     4.955013
    FALSE A2ML3                      7.415313                     7.382694
    FALSE A2ML4                      5.411328                     5.147779
    FALSE A4GALT                     5.738381                     5.347628
    FALSE A4GNT                      4.585769                     4.887529
    FALSE        x.r67.blu35_male_gonad_bldg x.w192.o157_male_gonad_inc.d9
    FALSE A2ML1                     5.876714                      6.114584
    FALSE A2ML2                     4.585769                      4.585769
    FALSE A2ML3                     7.403948                      7.508026
    FALSE A2ML4                     4.585769                      5.066372
    FALSE A4GALT                    5.876714                      5.641679
    FALSE A4GNT                     4.585769                      4.926388
    FALSE        x.y.s.ATLAS_male_gonad_control x.y132.w76_male_gonad_inc.d17
    FALSE A2ML1                        6.171172                      6.295959
    FALSE A2ML2                        4.852226                      4.585769
    FALSE A2ML3                        7.163941                      7.294480
    FALSE A2ML4                        5.045989                      5.364469
    FALSE A4GALT                       5.774639                      5.596925
    FALSE A4GNT                        4.962065                      4.585769
    FALSE        x.y141.w116_male_gonad_inc.d9 y129.x_male_gonad_n9
    FALSE A2ML1                       6.218385             6.287867
    FALSE A2ML2                       4.788440             5.054381
    FALSE A2ML3                       7.127091             7.340209
    FALSE A2ML4                       5.336229             5.245637
    FALSE A4GALT                      5.456637             5.679796
    FALSE A4GNT                       4.788440             4.857113
    FALSE        y131.w185.x_male_gonad_n9 y133.w77.r58_male_gonad_inc.d17
    FALSE A2ML1                   6.204286                        6.079672
    FALSE A2ML2                   4.585769                        4.910687
    FALSE A2ML3                   7.281473                        7.222637
    FALSE A2ML4                   5.462488                        5.146203
    FALSE A4GALT                  5.593174                        5.686784
    FALSE A4GNT                   4.585769                        4.585769
    FALSE        y149.r52.x_male_gonad_inc.d3 y95.g131.x_male_gonad_inc.d9
    FALSE A2ML1                      5.901603                     6.431413
    FALSE A2ML2                      4.585769                     4.859119
    FALSE A2ML3                      7.266361                     7.401964
    FALSE A2ML4                      5.338536                     4.971770
    FALSE A4GALT                     5.531753                     5.554518
    FALSE A4GNT                      4.585769                     4.971770
    FALSE        y98.o50.x_male_gonad_inc.d3
    FALSE A2ML1                     6.108740
    FALSE A2ML2                     4.910595
    FALSE A2ML3                     7.628793
    FALSE A2ML4                     4.910595
    FALSE A4GALT                    5.491509
    FALSE A4GNT                     4.585769
    FALSE 'data.frame': 4224 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13698 levels "A2ML1","A2ML2",..: 10328 12235 6508 10821 2715 3013 1619 9693 10539 3267 ...
    FALSE  $ padj     : num  8.87e-04 4.97e-02 1.32e-07 7.12e-04 1.64e-03 ...
    FALSE  $ logpadj  : num  3.05 1.3 6.88 3.15 2.78 ...
    FALSE  $ lfc      : num  3.87 2.77 2.37 2.36 2.32 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_gonad": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "control","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene         padj  logpadj      lfc  sextissue direction
    FALSE 1        ROCK2 8.872343e-04 3.051962 3.874341 male_gonad      bldg
    FALSE 2      TIMM10B 4.972620e-02 1.303415 2.769143 male_gonad      bldg
    FALSE 3 LOC107051332 1.317653e-07 6.880199 2.367905 male_gonad      bldg
    FALSE 4       SEMA6A 7.117497e-04 3.147673 2.358634 male_gonad      bldg
    FALSE 5      CYP24A1 1.643143e-03 2.784325 2.324644 male_gonad      bldg
    FALSE 6         DMB2 1.048348e-02 1.979494 2.323495 male_gonad      bldg
    FALSE 'data.frame': 118 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13698 levels "A2ML1","A2ML2",..: 8396 9397 12747 1068 764 9055 1153 1081 11446 4414 ...
    FALSE  $ padj     : num  0.0594 0.0722 0.0934 0.0673 0.0564 ...
    FALSE  $ logpadj  : num  1.23 1.14 1.03 1.17 1.25 ...
    FALSE  $ lfc      : num  0.653 0.652 0.521 0.394 0.393 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_gonad": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "bldg","NS","lay": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE     gene       padj  logpadj       lfc  sextissue direction
    FALSE 1   NOD1 0.05938766 1.226304 0.6532077 male_gonad       lay
    FALSE 2   POLQ 0.07217771 1.141597 0.6522237 male_gonad       lay
    FALSE 3  TRPM3 0.09335596 1.029858 0.5213148 male_gonad       lay
    FALSE 4  BCMO1 0.06725603 1.172269 0.3938371 male_gonad       lay
    FALSE 5 ARRDC4 0.05639252 1.248779 0.3927649 male_gonad       lay
    FALSE 6   PEX2 0.05639252 1.248779 0.3858421 male_gonad       lay
    FALSE 'data.frame': 135 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13698 levels "A2ML1","A2ML2",..: 565 884 6651 4860 6146 10930 10469 10494 10471 11021 ...
    FALSE  $ padj     : num  0.0556 0.0859 0.0986 0.0989 0.0986 ...
    FALSE  $ logpadj  : num  1.26 1.07 1.01 1 1.01 ...
    FALSE  $ lfc      : num  1.905 1.697 1.56 1.007 0.996 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_gonad": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "lay","NS","inc.d3": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene       padj  logpadj       lfc  sextissue direction
    FALSE 1        ANXA5 0.05556781 1.255177 1.9052771 male_gonad    inc.d3
    FALSE 2       ATP2B4 0.08590479 1.065983 1.6973052 male_gonad    inc.d3
    FALSE 3 LOC107055658 0.09857356 1.006240 1.5602205 male_gonad    inc.d3
    FALSE 4     HIST2H4B 0.09889027 1.004846 1.0070840 male_gonad    inc.d3
    FALSE 5 LOC101748962 0.09857356 1.006240 0.9958423 male_gonad    inc.d3
    FALSE 6         SGCZ 0.06341736 1.197792 0.9429218 male_gonad    inc.d3
    FALSE 'data.frame': 3 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13698 levels "A2ML1","A2ML2",..: 2725 2780 10658
    FALSE  $ padj     : num  0.0895 0.0806 0.0895
    FALSE  $ logpadj  : num  1.05 1.09 1.05
    FALSE  $ lfc      : num  0.817 -0.337 -0.436
    FALSE  $ sextissue: Factor w/ 1 level "male_gonad": 1 1 1
    FALSE  $ direction: Factor w/ 3 levels "inc.d3","NS",..: 3 1 1
    FALSE NULL
    FALSE       gene       padj  logpadj        lfc  sextissue direction
    FALSE 1 CYP2J2L3 0.08950931 1.048132  0.8169641 male_gonad    inc.d9
    FALSE 2    DAPK2 0.08056183 1.093871 -0.3374675 male_gonad    inc.d3
    FALSE 3   SAP30L 0.08950931 1.048132 -0.4359006 male_gonad    inc.d3
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13698 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "male_gonad": 
    FALSE  $ direction: Factor w/ 3 levels "inc.d9","NS",..: 
    FALSE NULL
    FALSE [1] gene      padj      logpadj   lfc       sextissue direction
    FALSE <0 rows> (or 0-length row.names)
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13698 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "male_gonad": 
    FALSE  $ direction: Factor w/ 3 levels "inc.d17","NS",..: 
    FALSE NULL
    FALSE [1] gene      padj      logpadj   lfc       sextissue direction
    FALSE <0 rows> (or 0-length row.names)
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13698 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "male_gonad": 
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","n5": 
    FALSE NULL
    FALSE [1] gene      padj      logpadj   lfc       sextissue direction
    FALSE <0 rows> (or 0-length row.names)
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13698 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "male_gonad": 
    FALSE  $ direction: Factor w/ 3 levels "n5","NS","n9": 
    FALSE NULL
    FALSE [1] gene      padj      logpadj   lfc       sextissue direction
    FALSE <0 rows> (or 0-length row.names)
    FALSE class: DESeqDataSet 
    FALSE dim: 13484 94 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13484): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(94): L.Blu13_male_hypothalamus_control.NYNO
    FALSE   L.G107_male_hypothalamus_control ...
    FALSE   y95.g131.x_male_hypothalamus_inc.d9
    FALSE   y98.o50.x_male_hypothalamus_inc.d3
    FALSE colData names(8): V1 bird ... study sextissue
    FALSE [1] 13484    94
    FALSE        L.Blu13_male_hypothalamus_control.NYNO
    FALSE A2ML1                                7.863034
    FALSE A2ML2                                5.419856
    FALSE A2ML3                               12.700124
    FALSE A2ML4                                4.643026
    FALSE A4GALT                               5.729120
    FALSE AAAS                                 6.984973
    FALSE        L.G107_male_hypothalamus_control L.R3_male_hypothalamus_control
    FALSE A2ML1                          8.687471                       7.280754
    FALSE A2ML2                          5.338276                       4.643026
    FALSE A2ML3                         13.131036                      12.903700
    FALSE A2ML4                          4.643026                       4.643026
    FALSE A4GALT                         5.338276                       5.469492
    FALSE AAAS                           6.730763                       6.991186
    FALSE        L.R8_male_hypothalamus_control L.W33_male_hypothalamus_control.NYNO
    FALSE A2ML1                        7.971976                             8.340211
    FALSE A2ML2                        4.643026                             5.098345
    FALSE A2ML3                       12.594947                            12.713689
    FALSE A2ML4                        4.643026                             5.284327
    FALSE A4GALT                       4.643026                             5.098345
    FALSE AAAS                         6.672719                             6.665901
    FALSE        L.W3_male_hypothalamus_control L.W4_male_hypothalamus_control
    FALSE A2ML1                        8.396194                       7.418188
    FALSE A2ML2                        4.643026                       5.019213
    FALSE A2ML3                       13.386500                      13.273104
    FALSE A2ML4                        4.643026                       4.643026
    FALSE A4GALT                       4.643026                       5.173549
    FALSE AAAS                         6.063757                       6.840082
    FALSE        R.Y108.W29_male_hypothalamus_control.NYNO
    FALSE A2ML1                                   8.322348
    FALSE A2ML2                                   4.643026
    FALSE A2ML3                                  12.948854
    FALSE A2ML4                                   5.587858
    FALSE A4GALT                                  5.316982
    FALSE AAAS                                    6.632315
    FALSE        blk12.x_male_hypothalamus_n5.NYNO blk17.x_male_hypothalamus_inc.d17
    FALSE A2ML1                           8.043080                          8.373641
    FALSE A2ML2                           5.234242                          5.194166
    FALSE A2ML3                          12.394451                         12.569342
    FALSE A2ML4                           5.234242                          5.194166
    FALSE A4GALT                          5.512855                          5.351732
    FALSE AAAS                            7.101546                          7.183207
    FALSE        blu104.w120.x_male_hypothalamus_hatch
    FALSE A2ML1                               7.063046
    FALSE A2ML2                               4.997466
    FALSE A2ML3                              12.694348
    FALSE A2ML4                               5.051958
    FALSE A4GALT                              5.500721
    FALSE AAAS                                7.187080
    FALSE        blu108.w40.o158_male_hypothalamus_inc.d9
    FALSE A2ML1                                  7.434692
    FALSE A2ML2                                  5.217934
    FALSE A2ML3                                 12.517868
    FALSE A2ML4                                  5.050883
    FALSE A4GALT                                 5.344868
    FALSE AAAS                                   6.817419
    FALSE        blu111.w113.x_male_hypothalamus_inc.d3
    FALSE A2ML1                                8.270320
    FALSE A2ML2                                5.197709
    FALSE A2ML3                               12.244041
    FALSE A2ML4                                5.036446
    FALSE A4GALT                               5.320326
    FALSE AAAS                                 6.751446
    FALSE        blu113.w124.x_male_hypothalamus_inc.d17
    FALSE A2ML1                                 7.092708
    FALSE A2ML2                                 4.996793
    FALSE A2ML3                                12.453744
    FALSE A2ML4                                 4.996793
    FALSE A4GALT                                5.426397
    FALSE AAAS                                  6.683471
    FALSE        blu114.r38.w198_male_hypothalamus_bldg
    FALSE A2ML1                                8.159537
    FALSE A2ML2                                5.334426
    FALSE A2ML3                               12.322238
    FALSE A2ML4                                5.015110
    FALSE A4GALT                               5.424919
    FALSE AAAS                                 6.841970
    FALSE        blu121.w91.x_male_hypothalamus_inc.d17
    FALSE A2ML1                                6.031695
    FALSE A2ML2                                4.868141
    FALSE A2ML3                               12.791926
    FALSE A2ML4                                4.643026
    FALSE A4GALT                               5.191690
    FALSE AAAS                                 7.733652
    FALSE        blu33.y88.x_male_hypothalamus_bldg blu37.r65.x_male_hypothalamus_n5
    FALSE A2ML1                            7.252618                         6.535857
    FALSE A2ML2                            5.338704                         4.969740
    FALSE A2ML3                           12.171387                        12.391502
    FALSE A2ML4                            5.212853                         4.874294
    FALSE A4GALT                           5.393263                         5.435037
    FALSE AAAS                             6.827549                         7.457371
    FALSE        blu41.y100.x_male_hypothalamus_n5.NYNO
    FALSE A2ML1                                7.810432
    FALSE A2ML2                                5.223458
    FALSE A2ML3                               12.459661
    FALSE A2ML4                                4.643026
    FALSE A4GALT                               5.506607
    FALSE AAAS                                 6.862207
    FALSE        blu81.r88_male_hypothalamus_n9 d.s008.y.blk_male_hypothalamus_n5
    FALSE A2ML1                        8.180836                          8.489376
    FALSE A2ML2                        5.075786                          5.613605
    FALSE A2ML3                       12.798433                         12.575551
    FALSE A2ML4                        5.387126                          5.092613
    FALSE A4GALT                       5.549509                          5.464250
    FALSE AAAS                         6.956806                          6.858198
    FALSE        d.s047.blk.o_male_hypothalamus_n5 g.s.blk.d_male_hypothalamus_n9
    FALSE A2ML1                           6.789554                       8.074147
    FALSE A2ML2                           4.874338                       5.438178
    FALSE A2ML3                          12.697067                      12.594865
    FALSE A2ML4                           4.643026                       5.105951
    FALSE A4GALT                          5.206638                       5.294953
    FALSE AAAS                            6.726111                       7.060858
    FALSE        g.s.blk.y_male_hypothalamus_lay g.s043.pu.blk_male_hypothalamus_lay
    FALSE A2ML1                         8.014005                            7.102712
    FALSE A2ML2                         4.643026                            5.118034
    FALSE A2ML3                        11.998781                           12.464753
    FALSE A2ML4                         4.922671                            5.254447
    FALSE A4GALT                        5.515374                            5.190705
    FALSE AAAS                          7.561193                            6.978471
    FALSE        g104.w82.x_male_hypothalamus_bldg
    FALSE A2ML1                           6.675244
    FALSE A2ML2                           5.065759
    FALSE A2ML3                          12.675633
    FALSE A2ML4                           5.065759
    FALSE A4GALT                          5.575362
    FALSE AAAS                            6.798875
    FALSE        g114.w83.x_male_hypothalamus_hatch
    FALSE A2ML1                            7.030941
    FALSE A2ML2                            5.058570
    FALSE A2ML3                           12.863407
    FALSE A2ML4                            4.883492
    FALSE A4GALT                           5.317787
    FALSE AAAS                             6.655802
    FALSE        g130.y81.x_male_hypothalamus_inc.d17
    FALSE A2ML1                              8.303835
    FALSE A2ML2                              5.197435
    FALSE A2ML3                             12.415093
    FALSE A2ML4                              4.643026
    FALSE A4GALT                             5.123889
    FALSE AAAS                               7.157022
    FALSE        g143.blu32.x_male_hypothalamus_inc.d17
    FALSE A2ML1                                8.485939
    FALSE A2ML2                                5.209589
    FALSE A2ML3                               12.034094
    FALSE A2ML4                                5.209589
    FALSE A4GALT                               5.334744
    FALSE AAAS                                 7.045369
    FALSE        g146.blu51_male_hypothalamus_inc.d3.NYNO
    FALSE A2ML1                                  8.124946
    FALSE A2ML2                                  4.993142
    FALSE A2ML3                                 12.193373
    FALSE A2ML4                                  4.643026
    FALSE A4GALT                                 4.993142
    FALSE AAAS                                   7.348389
    FALSE        g20.w106.x_male_hypothalamus_inc.d3
    FALSE A2ML1                             7.653864
    FALSE A2ML2                             5.138125
    FALSE A2ML3                            12.408670
    FALSE A2ML4                             5.138125
    FALSE A4GALT                            5.138125
    FALSE AAAS                              7.071199
    FALSE        g52.blu58_male_hypothalamus_bldg g53.y84_male_hypothalamus_hatch
    FALSE A2ML1                          8.248671                        7.646077
    FALSE A2ML2                          5.111683                        5.188209
    FALSE A2ML3                         12.718247                       12.583425
    FALSE A2ML4                          5.026240                        5.188209
    FALSE A4GALT                         5.183403                        4.643026
    FALSE AAAS                           6.792290                        6.936693
    FALSE        o.s.w.r_male_hypothalamus_lay o152.o120.w42_male_hypothalamus_n5
    FALSE A2ML1                       6.573063                           8.393650
    FALSE A2ML2                       4.643026                           5.194450
    FALSE A2ML3                      12.100042                          12.223893
    FALSE A2ML4                       4.643026                           5.194450
    FALSE A4GALT                      5.202773                           5.418252
    FALSE AAAS                        7.062425                           6.989636
    FALSE        o39.y77.x_male_hypothalamus_hatch
    FALSE A2ML1                           7.202611
    FALSE A2ML2                           5.042557
    FALSE A2ML3                          12.291121
    FALSE A2ML4                           5.206271
    FALSE A4GALT                          5.367174
    FALSE AAAS                            7.000220
    FALSE        o44.blu26.x_male_hypothalamus_hatch
    FALSE A2ML1                             7.809029
    FALSE A2ML2                             5.098227
    FALSE A2ML3                            12.568038
    FALSE A2ML4                             5.050505
    FALSE A4GALT                            5.344227
    FALSE AAAS                              6.894377
    FALSE        o48.r197.x_male_hypothalamus_inc.d3 o49.x_male_hypothalamus_inc.d9
    FALSE A2ML1                             8.480884                       6.744692
    FALSE A2ML2                             5.215012                       4.643026
    FALSE A2ML3                            12.250477                      12.621065
    FALSE A2ML4                             4.643026                       5.212447
    FALSE A4GALT                            5.215012                       5.212447
    FALSE AAAS                              6.771492                       6.852300
    FALSE        o57.g59_male_hypothalamus_inc.d9
    FALSE A2ML1                          7.133694
    FALSE A2ML2                          4.643026
    FALSE A2ML3                         12.521360
    FALSE A2ML4                          5.231693
    FALSE A4GALT                         5.361555
    FALSE AAAS                           7.088188
    FALSE        pk.s238.blk.w_male_hypothalamus_lay
    FALSE A2ML1                             8.356088
    FALSE A2ML2                             5.105439
    FALSE A2ML3                            12.529536
    FALSE A2ML4                             4.643026
    FALSE A4GALT                            5.404276
    FALSE AAAS                              7.073093
    FALSE        pk.w.s141.o_male_hypothalamus_lay
    FALSE A2ML1                           7.871233
    FALSE A2ML2                           5.083772
    FALSE A2ML3                          12.483171
    FALSE A2ML4                           5.181793
    FALSE A4GALT                          5.590245
    FALSE AAAS                            7.039959
    FALSE        r.s005.pk.blk_male_hypothalamus_lay
    FALSE A2ML1                             7.486563
    FALSE A2ML2                             4.984441
    FALSE A2ML3                            12.819952
    FALSE A2ML4                             5.180990
    FALSE A4GALT                            5.435582
    FALSE AAAS                              6.950106
    FALSE        r.s059.d.o_male_hypothalamus_bldg
    FALSE A2ML1                           7.770300
    FALSE A2ML2                           5.190596
    FALSE A2ML3                          12.194832
    FALSE A2ML4                           5.190596
    FALSE A4GALT                          5.697898
    FALSE AAAS                            6.974461
    FALSE        r.s116.blk.pu_male_hypothalamus_lay
    FALSE A2ML1                             7.552892
    FALSE A2ML2                             5.012830
    FALSE A2ML3                            12.198630
    FALSE A2ML4                             5.330246
    FALSE A4GALT                            5.012830
    FALSE AAAS                              7.089304
    FALSE        r.y.s007.blk_male_hypothalamus_n9
    FALSE A2ML1                           7.290425
    FALSE A2ML2                           4.643026
    FALSE A2ML3                          12.431953
    FALSE A2ML4                           5.058451
    FALSE A4GALT                          5.357721
    FALSE AAAS                            7.166620
    FALSE        r176.blu54_male_hypothalamus_inc.d17
    FALSE A2ML1                              7.315059
    FALSE A2ML2                              5.093933
    FALSE A2ML3                             12.621110
    FALSE A2ML4                              4.643026
    FALSE A4GALT                             5.417856
    FALSE AAAS                               7.220444
    FALSE        r190.o43.x_male_hypothalamus_lay r195.x_male_hypothalamus_n9
    FALSE A2ML1                          8.040818                    8.446858
    FALSE A2ML2                          5.104985                    5.214394
    FALSE A2ML3                         12.733512                   12.356514
    FALSE A2ML4                          5.293603                    4.643026
    FALSE A4GALT                         5.436545                    5.280812
    FALSE AAAS                           6.774934                    6.832603
    FALSE        r37.w100.x_male_hypothalamus_n9 r41.w99.x_male_hypothalamus_hatch
    FALSE A2ML1                         8.267524                          8.063920
    FALSE A2ML2                         5.380494                          5.118850
    FALSE A2ML3                        12.406472                         12.674720
    FALSE A2ML4                         5.430222                          5.243294
    FALSE A4GALT                        5.430222                          5.375619
    FALSE AAAS                          7.049551                          6.877959
    FALSE        r45.x_male_hypothalamus_inc.d9 r72.y83.x_male_hypothalamus_hatch
    FALSE A2ML1                        7.567129                          7.363802
    FALSE A2ML2                        4.643026                          4.643026
    FALSE A2ML3                       12.588931                         12.463099
    FALSE A2ML4                        5.172543                          4.949836
    FALSE A4GALT                       5.172543                          4.949836
    FALSE AAAS                         6.858899                          7.173198
    FALSE        s.pu148.blk.r_male_hypothalamus_bldg
    FALSE A2ML1                              7.632429
    FALSE A2ML2                              5.166971
    FALSE A2ML3                             12.701345
    FALSE A2ML4                              4.643026
    FALSE A4GALT                             5.014523
    FALSE AAAS                               6.839031
    FALSE        s065.l.d.o_male_hypothalamus_bldg s066.l.d.r_male_hypothalamus_bldg
    FALSE A2ML1                           8.014537                          7.657472
    FALSE A2ML2                           5.480468                          5.324029
    FALSE A2ML3                          12.389282                         12.506219
    FALSE A2ML4                           4.643026                          5.038623
    FALSE A4GALT                          5.427717                          5.671070
    FALSE AAAS                            6.926939                          7.129654
    FALSE        s150.w.g.blk_male_hypothalamus_lay s187.l.o.r_male_hypothalamus_n9
    FALSE A2ML1                            6.628053                        8.424432
    FALSE A2ML2                            4.643026                        5.563691
    FALSE A2ML3                           12.157983                       11.907269
    FALSE A2ML4                            4.856303                        5.060326
    FALSE A4GALT                           5.118211                        5.124329
    FALSE AAAS                             7.119195                        6.939854
    FALSE        s243.blk.pk.r_male_hypothalamus_lay w34.x_male_hypothalamus_inc.d9
    FALSE A2ML1                             7.592243                       8.121567
    FALSE A2ML2                             4.860959                       5.054702
    FALSE A2ML3                            12.561922                      12.490797
    FALSE A2ML4                             5.128508                       5.223282
    FALSE A4GALT                            5.420113                       5.173312
    FALSE AAAS                              7.356095                       6.913866
    FALSE        x.blk.blk.ATLAS_male_hypothalamus_control
    FALSE A2ML1                                   8.154061
    FALSE A2ML2                                   5.079237
    FALSE A2ML3                                  12.542338
    FALSE A2ML4                                   5.079237
    FALSE A4GALT                                  4.643026
    FALSE AAAS                                    5.257617
    FALSE        x.blk16_male_hypothalamus_n9.NYNO
    FALSE A2ML1                           8.230217
    FALSE A2ML2                           4.643026
    FALSE A2ML3                          12.209832
    FALSE A2ML4                           5.300318
    FALSE A4GALT                          5.399874
    FALSE AAAS                            7.096680
    FALSE        x.blu106.o153_male_hypothalamus_inc.d9
    FALSE A2ML1                                8.200444
    FALSE A2ML2                                5.578479
    FALSE A2ML3                               12.082324
    FALSE A2ML4                                5.149251
    FALSE A4GALT                               5.610153
    FALSE AAAS                                 7.124259
    FALSE        x.blu117.w89_male_hypothalamus_inc.d17
    FALSE A2ML1                                7.340652
    FALSE A2ML2                                5.317821
    FALSE A2ML3                               12.251501
    FALSE A2ML4                                4.643026
    FALSE A4GALT                               5.550043
    FALSE AAAS                                 7.150399
    FALSE        x.blu23.w14_male_hypothalamus_n9 x.blu30_male_hypothalamus_n5
    FALSE A2ML1                          7.625020                     8.359169
    FALSE A2ML2                          4.643026                     5.417567
    FALSE A2ML3                         12.290667                    11.998159
    FALSE A2ML4                          5.035362                     5.277923
    FALSE A4GALT                         5.196189                     5.277923
    FALSE AAAS                           6.848030                     6.919277
    FALSE        x.blu42.o28_male_hypothalamus_inc.d3.NYNO
    FALSE A2ML1                                   7.359571
    FALSE A2ML2                                   5.253488
    FALSE A2ML3                                  12.503712
    FALSE A2ML4                                   4.917659
    FALSE A4GALT                                  5.253488
    FALSE AAAS                                    7.059415
    FALSE        x.g.ATLAS_male_hypothalamus_control
    FALSE A2ML1                             8.208314
    FALSE A2ML2                             5.354327
    FALSE A2ML3                            12.160255
    FALSE A2ML4                             5.509937
    FALSE A4GALT                            4.643026
    FALSE AAAS                              5.942854
    FALSE        x.g.g.ATLAS_male_hypothalamus_control
    FALSE A2ML1                               8.425908
    FALSE A2ML2                               4.643026
    FALSE A2ML3                              11.933571
    FALSE A2ML4                               4.643026
    FALSE A4GALT                              4.643026
    FALSE AAAS                                7.591821
    FALSE        x.g13.w109_male_hypothalamus_inc.d9
    FALSE A2ML1                             8.383485
    FALSE A2ML2                             5.326710
    FALSE A2ML3                            12.040860
    FALSE A2ML4                             5.040200
    FALSE A4GALT                            5.362961
    FALSE AAAS                              6.996153
    FALSE        x.g14.w199_male_hypothalamus_inc.d17
    FALSE A2ML1                              7.243934
    FALSE A2ML2                              4.643026
    FALSE A2ML3                             12.380448
    FALSE A2ML4                              5.091186
    FALSE A4GALT                             5.347463
    FALSE AAAS                               7.032064
    FALSE        x.g147.blu28_male_hypothalamus_inc.d3 x.g70_male_hypothalamus_hatch
    FALSE A2ML1                               7.930698                      7.770444
    FALSE A2ML2                               5.479407                      4.886021
    FALSE A2ML3                              12.496641                     12.176259
    FALSE A2ML4                               5.319069                      5.234786
    FALSE A4GALT                              5.525835                      5.183850
    FALSE AAAS                                7.243236                      6.716114
    FALSE        x.o160.w102_male_hypothalamus_hatch
    FALSE A2ML1                             7.824114
    FALSE A2ML2                             5.267627
    FALSE A2ML3                            12.518693
    FALSE A2ML4                             4.880682
    FALSE A4GALT                            5.267627
    FALSE AAAS                              6.874843
    FALSE        x.o163.w101_male_hypothalamus_inc.d3
    FALSE A2ML1                              7.852115
    FALSE A2ML2                              5.210111
    FALSE A2ML3                             12.269669
    FALSE A2ML4                              5.173907
    FALSE A4GALT                             5.173907
    FALSE AAAS                               6.876383
    FALSE        x.o164.w123_male_hypothalamus_n5 x.o2_male_hypothalamus_n9
    FALSE A2ML1                          7.683803                  8.066473
    FALSE A2ML2                          5.067028                  5.049328
    FALSE A2ML3                         12.329566                 12.544365
    FALSE A2ML4                          4.888410                  4.643026
    FALSE A4GALT                         5.287657                  5.540030
    FALSE AAAS                           6.917510                  7.224266
    FALSE        x.o30.g134_male_hypothalamus_bldg
    FALSE A2ML1                           8.530812
    FALSE A2ML2                           5.167460
    FALSE A2ML3                          12.228472
    FALSE A2ML4                           4.946916
    FALSE A4GALT                          4.643026
    FALSE AAAS                            6.825993
    FALSE        x.o47.y82_male_hypothalamus_inc.d9 x.o68_male_hypothalamus_n5
    FALSE A2ML1                            6.780096                   8.052485
    FALSE A2ML2                            4.643026                   5.155680
    FALSE A2ML3                           12.428847                  12.383841
    FALSE A2ML4                            4.983253                   5.062335
    FALSE A4GALT                           5.396953                   5.420792
    FALSE AAAS                             6.720788                   6.869110
    FALSE        x.r178_male_hypothalamus_hatch x.r181_male_hypothalamus_n5
    FALSE A2ML1                        7.819956                    7.868084
    FALSE A2ML2                        5.058511                    5.195427
    FALSE A2ML3                       12.347107                   12.748520
    FALSE A2ML4                        5.357822                    5.122143
    FALSE A4GALT                       5.357822                    5.259704
    FALSE AAAS                         7.130422                    7.011544
    FALSE        x.r29.w96_male_hypothalamus_inc.d17
    FALSE A2ML1                             6.619089
    FALSE A2ML2                             5.123224
    FALSE A2ML3                            12.337615
    FALSE A2ML4                             5.049399
    FALSE A4GALT                            5.186811
    FALSE AAAS                              7.056191
    FALSE        x.r64.g140_male_hypothalamus_inc.d3
    FALSE A2ML1                             8.213872
    FALSE A2ML2                             5.247253
    FALSE A2ML3                            12.151913
    FALSE A2ML4                             5.176758
    FALSE A4GALT                            5.578180
    FALSE AAAS                              7.174680
    FALSE        x.r67.blu35_male_hypothalamus_bldg
    FALSE A2ML1                            7.928929
    FALSE A2ML2                            4.931378
    FALSE A2ML3                           12.170066
    FALSE A2ML4                            4.931378
    FALSE A4GALT                           5.140824
    FALSE AAAS                             7.103990
    FALSE        x.w192.o157_male_hypothalamus_inc.d9
    FALSE A2ML1                              8.451696
    FALSE A2ML2                              5.392777
    FALSE A2ML3                             12.255085
    FALSE A2ML4                              5.257453
    FALSE A4GALT                             5.328708
    FALSE AAAS                               6.747134
    FALSE        x.y132.w76_male_hypothalamus_inc.d17
    FALSE A2ML1                              8.092496
    FALSE A2ML2                              5.080432
    FALSE A2ML3                             12.674092
    FALSE A2ML4                              4.896206
    FALSE A4GALT                             5.472341
    FALSE AAAS                               7.236002
    FALSE        x.y141.w116_male_hypothalamus_inc.d9 y129.x_male_hypothalamus_n9
    FALSE A2ML1                              8.244622                    8.064276
    FALSE A2ML2                              5.370735                    5.373500
    FALSE A2ML3                             12.201735                   12.276098
    FALSE A2ML4                              5.239265                    5.222796
    FALSE A4GALT                             5.620050                    5.496272
    FALSE AAAS                               6.815110                    7.133455
    FALSE        y131.w185.x_male_hypothalamus_n9
    FALSE A2ML1                          8.031526
    FALSE A2ML2                          5.422234
    FALSE A2ML3                         12.644932
    FALSE A2ML4                          5.241124
    FALSE A4GALT                         5.319848
    FALSE AAAS                           6.818827
    FALSE        y133.w77.r58_male_hypothalamus_inc.d17
    FALSE A2ML1                                8.250452
    FALSE A2ML2                                5.440733
    FALSE A2ML3                               12.492821
    FALSE A2ML4                                4.972139
    FALSE A4GALT                               5.045672
    FALSE AAAS                                 6.888576
    FALSE        y149.r52.x_male_hypothalamus_inc.d3
    FALSE A2ML1                             7.993270
    FALSE A2ML2                             5.157986
    FALSE A2ML3                            12.317356
    FALSE A2ML4                             4.838548
    FALSE A4GALT                            5.367525
    FALSE AAAS                              6.819234
    FALSE        y95.g131.x_male_hypothalamus_inc.d9
    FALSE A2ML1                             7.648417
    FALSE A2ML2                             4.987404
    FALSE A2ML3                            12.571650
    FALSE A2ML4                             5.064305
    FALSE A4GALT                            5.236726
    FALSE AAAS                              6.920141
    FALSE        y98.o50.x_male_hypothalamus_inc.d3
    FALSE A2ML1                            8.187021
    FALSE A2ML2                            5.335176
    FALSE A2ML3                           12.181694
    FALSE A2ML4                            4.971735
    FALSE A4GALT                           4.971735
    FALSE AAAS                             6.809217
    FALSE 'data.frame': 6683 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13484 levels "A2ML1","A2ML2",..: 11340 3565 6303 12401 6511 3212 6573 4229 1721 6076 ...
    FALSE  $ padj     : num  1.21e-09 5.09e-10 3.70e-10 1.49e-03 2.78e-05 ...
    FALSE  $ logpadj  : num  8.92 9.29 9.43 2.83 4.56 ...
    FALSE  $ lfc      : num  2.96 2.79 2.47 2.45 2.41 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "control","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene         padj  logpadj      lfc         sextissue direction
    FALSE 1         SOX4 1.212858e-09 8.916190 2.959749 male_hypothalamus      bldg
    FALSE 2         F8A3 5.086667e-10 9.293567 2.785894 male_hypothalamus      bldg
    FALSE 3 LOC107050672 3.695351e-10 9.432344 2.468940 male_hypothalamus      bldg
    FALSE 4         TPH2 1.492246e-03 2.826160 2.445432 male_hypothalamus      bldg
    FALSE 5 LOC107055647 2.777762e-05 4.556305 2.411410 male_hypothalamus      bldg
    FALSE 6         ECM1 2.233825e-09 8.650951 2.392867 male_hypothalamus      bldg
    FALSE 'data.frame': 1 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13484 levels "A2ML1","A2ML2",..: 12401
    FALSE  $ padj     : num 0.0131
    FALSE  $ logpadj  : num 1.88
    FALSE  $ lfc      : num -3.22
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1
    FALSE  $ direction: Factor w/ 3 levels "bldg","NS","lay": 1
    FALSE NULL
    FALSE   gene       padj  logpadj       lfc         sextissue direction
    FALSE 1 TPH2 0.01310221 1.882655 -3.224611 male_hypothalamus      bldg
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13484 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 
    FALSE  $ direction: Factor w/ 3 levels "lay","NS","inc.d3": 
    FALSE NULL
    FALSE [1] gene      padj      logpadj   lfc       sextissue direction
    FALSE <0 rows> (or 0-length row.names)
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13484 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 
    FALSE  $ direction: Factor w/ 3 levels "inc.d3","NS",..: 
    FALSE NULL
    FALSE [1] gene      padj      logpadj   lfc       sextissue direction
    FALSE <0 rows> (or 0-length row.names)
    FALSE 'data.frame': 87 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13484 levels "A2ML1","A2ML2",..: 5046 5045 2329 1760 3797 1337 2532 6228 5845 3919 ...
    FALSE  $ padj     : num  0.0471 0.0935 0.0196 0.0771 0.0471 ...
    FALSE  $ logpadj  : num  1.33 1.03 1.71 1.11 1.33 ...
    FALSE  $ lfc      : num  3.63 2.65 1.28 1.11 1.02 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "inc.d9","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE     gene       padj  logpadj       lfc         sextissue direction
    FALSE 1  IGLL1 0.04713019 1.326701 3.6319750 male_hypothalamus   inc.d17
    FALSE 2    IGJ 0.09348029 1.029280 2.6518341 male_hypothalamus   inc.d17
    FALSE 3 COL1A1 0.01960242 1.707690 1.2799150 male_hypothalamus   inc.d17
    FALSE 4  CCM2L 0.07706307 1.113154 1.1078543 male_hypothalamus   inc.d17
    FALSE 5 FBLIM1 0.04713019 1.326701 1.0171116 male_hypothalamus   inc.d17
    FALSE 6   C1QB 0.09868602 1.005744 0.9078655 male_hypothalamus   inc.d17
    FALSE 'data.frame': 6 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13484 levels "A2ML1","A2ML2",..: 7902 7861 6012 863 2532 5882
    FALSE  $ padj     : num  0.0961 0.0984 0.0961 0.0836 0.0984 ...
    FALSE  $ logpadj  : num  1.02 1.01 1.02 1.08 1.01 ...
    FALSE  $ lfc      : num  2.949 2.456 2.423 1.951 -0.838 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1 1 1 1 1 1
    FALSE  $ direction: Factor w/ 3 levels "inc.d17","NS",..: 3 3 3 3 1 1
    FALSE NULL
    FALSE           gene       padj  logpadj        lfc         sextissue direction
    FALSE 1        MYOZ1 0.09607232 1.017402  2.9488760 male_hypothalamus     hatch
    FALSE 2        MYL10 0.09837662 1.007108  2.4559578 male_hypothalamus     hatch
    FALSE 3 LOC101748402 0.09607232 1.017402  2.4228230 male_hypothalamus     hatch
    FALSE 4       ATP2A1 0.08357535 1.077922  1.9506186 male_hypothalamus     hatch
    FALSE 5        CSF3R 0.09837662 1.007108 -0.8378164 male_hypothalamus   inc.d17
    FALSE 6 LOC100858707 0.09837662 1.007108 -1.0314972 male_hypothalamus   inc.d17
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13484 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","n5": 
    FALSE NULL
    FALSE [1] gene      padj      logpadj   lfc       sextissue direction
    FALSE <0 rows> (or 0-length row.names)
    FALSE 'data.frame': 1 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13484 levels "A2ML1","A2ML2",..: 5005
    FALSE  $ padj     : num 0.0753
    FALSE  $ logpadj  : num 1.12
    FALSE  $ lfc      : num -3.28
    FALSE  $ sextissue: Factor w/ 1 level "male_hypothalamus": 1
    FALSE  $ direction: Factor w/ 3 levels "n5","NS","n9": 1
    FALSE NULL
    FALSE   gene       padj  logpadj       lfc         sextissue direction
    FALSE 1 IFI6 0.07525269 1.123478 -3.276777 male_hypothalamus        n5
    FALSE class: DESeqDataSet 
    FALSE dim: 13368 97 
    FALSE metadata(1): version
    FALSE assays(1): counts
    FALSE rownames(13368): A2ML1 A2ML2 ... ZYX ZZZ3
    FALSE rowData names(0):
    FALSE colnames(97): L.Blu13_male_pituitary_control.NYNO
    FALSE   L.G107_male_pituitary_control ...
    FALSE   y95.g131.x_male_pituitary_inc.d9 y98.o50.x_male_pituitary_inc.d3
    FALSE colData names(8): V1 bird ... study sextissue
    FALSE [1] 13368    97
    FALSE        L.Blu13_male_pituitary_control.NYNO L.G107_male_pituitary_control
    FALSE A2ML1                             7.975608                      7.455768
    FALSE A2ML2                             5.995081                      5.995081
    FALSE A2ML3                             8.223021                      7.294633
    FALSE A2ML4                             7.518225                      6.500453
    FALSE A4GALT                            7.576528                      7.618167
    FALSE AAAS                              7.815814                      7.798463
    FALSE        L.R3_male_pituitary_control.NYNO L.R8_male_pituitary_control
    FALSE A2ML1                          8.056285                    7.438916
    FALSE A2ML2                          6.661045                    5.995081
    FALSE A2ML3                          8.188504                    7.932039
    FALSE A2ML4                          6.523303                    6.635682
    FALSE A4GALT                         7.547632                    7.810536
    FALSE AAAS                           7.532443                    7.744974
    FALSE        L.W33_male_pituitary_control L.W3_male_pituitary_control
    FALSE A2ML1                      8.089224                    7.777046
    FALSE A2ML2                      6.320518                    5.995081
    FALSE A2ML3                      7.902721                    7.675319
    FALSE A2ML4                      6.845565                    6.799868
    FALSE A4GALT                     7.364632                    6.974645
    FALSE AAAS                       7.745026                    7.440779
    FALSE        L.W4_male_pituitary_control R.Y108.W29_male_pituitary_control
    FALSE A2ML1                     7.734311                          7.755913
    FALSE A2ML2                     5.995081                          5.995081
    FALSE A2ML3                     8.305309                          8.536533
    FALSE A2ML4                     6.830579                          6.623024
    FALSE A4GALT                    7.813524                          7.520980
    FALSE AAAS                      7.875838                          7.778144
    FALSE        blk12.x_male_pituitary_n5 blk17.x_male_pituitary_inc.d17
    FALSE A2ML1                   7.485313                       7.457277
    FALSE A2ML2                   6.212242                       6.223965
    FALSE A2ML3                   7.743318                       7.808488
    FALSE A2ML4                   6.676114                       6.344222
    FALSE A4GALT                  7.498583                       7.506859
    FALSE AAAS                    8.002775                       8.113590
    FALSE        blu104.w120.x_male_pituitary_hatch.NYNO
    FALSE A2ML1                                 7.181305
    FALSE A2ML2                                 5.995081
    FALSE A2ML3                                 8.398549
    FALSE A2ML4                                 7.270779
    FALSE A4GALT                                6.692590
    FALSE AAAS                                  8.286116
    FALSE        blu108.w40.o158_male_pituitary_inc.d9
    FALSE A2ML1                               7.383439
    FALSE A2ML2                               5.995081
    FALSE A2ML3                               7.957373
    FALSE A2ML4                               6.440625
    FALSE A4GALT                              7.332622
    FALSE AAAS                                7.794859
    FALSE        blu111.w113.x_male_pituitary_inc.d3
    FALSE A2ML1                             7.531010
    FALSE A2ML2                             6.290397
    FALSE A2ML3                             7.426577
    FALSE A2ML4                             6.451476
    FALSE A4GALT                            7.462487
    FALSE AAAS                              7.735145
    FALSE        blu113.w124.x_male_pituitary_inc.d17.NYNO
    FALSE A2ML1                                   7.587996
    FALSE A2ML2                                   6.196572
    FALSE A2ML3                                   8.906421
    FALSE A2ML4                                   6.465911
    FALSE A4GALT                                  7.500781
    FALSE AAAS                                    7.506433
    FALSE        blu114.r38.w198_male_pituitary_bldg
    FALSE A2ML1                             7.751318
    FALSE A2ML2                             6.256115
    FALSE A2ML3                             7.577964
    FALSE A2ML4                             6.921537
    FALSE A4GALT                            7.485349
    FALSE AAAS                              7.902880
    FALSE        blu121.w91.x_male_pituitary_inc.d17 blu33.y88.x_male_pituitary_bldg
    FALSE A2ML1                             7.729622                        7.197791
    FALSE A2ML2                             6.257706                        6.207178
    FALSE A2ML3                             7.776221                        8.171864
    FALSE A2ML4                             6.365980                        6.512304
    FALSE A4GALT                            7.140849                        7.623202
    FALSE AAAS                              7.647696                        7.577155
    FALSE        blu37.r65.x_male_pituitary_n5 blu41.y100.x_male_pituitary_n5
    FALSE A2ML1                       7.558549                       7.604341
    FALSE A2ML2                       6.195413                       6.263379
    FALSE A2ML3                       8.350307                       7.613493
    FALSE A2ML4                       6.341512                       6.660870
    FALSE A4GALT                      7.422978                       7.381748
    FALSE AAAS                        8.108753                       7.852256
    FALSE        blu81.r88_male_pituitary_n9 d.s008.y.blk_male_pituitary_n5
    FALSE A2ML1                     7.559509                       7.736650
    FALSE A2ML2                     6.302115                       6.495867
    FALSE A2ML3                     8.212037                       8.086017
    FALSE A2ML4                     6.566826                       6.350074
    FALSE A4GALT                    7.078558                       7.707466
    FALSE AAAS                      7.608530                       7.646868
    FALSE        d.s047.blk.o_male_pituitary_n5 g.s.blk.d_male_pituitary_n9
    FALSE A2ML1                        7.567323                    7.471870
    FALSE A2ML2                        6.434682                    6.335719
    FALSE A2ML3                        8.063528                    7.733052
    FALSE A2ML4                        6.574960                    6.397758
    FALSE A4GALT                       7.431540                    7.395980
    FALSE AAAS                         7.722838                    7.733052
    FALSE        g.s.blk.y_male_pituitary_lay g.s043.pu.blk_male_pituitary_lay
    FALSE A2ML1                      7.689133                         7.166592
    FALSE A2ML2                      6.323572                         6.182976
    FALSE A2ML3                      7.860361                         7.835326
    FALSE A2ML4                      6.512812                         6.715791
    FALSE A4GALT                     7.458036                         7.602834
    FALSE AAAS                       7.727027                         7.793808
    FALSE        g104.w82.x_male_pituitary_bldg g114.w83.x_male_pituitary_hatch.NYNO
    FALSE A2ML1                        7.602575                             7.414819
    FALSE A2ML2                        6.127966                             6.157197
    FALSE A2ML3                        7.680976                             8.093329
    FALSE A2ML4                        6.750072                             6.318804
    FALSE A4GALT                       7.479033                             7.730749
    FALSE AAAS                         7.955214                             7.520519
    FALSE        g130.y81.x_male_pituitary_inc.d17
    FALSE A2ML1                           7.512865
    FALSE A2ML2                           6.361353
    FALSE A2ML3                           7.680616
    FALSE A2ML4                           7.012618
    FALSE A4GALT                          7.309230
    FALSE AAAS                            8.320715
    FALSE        g143.blu32.x_male_pituitary_inc.d17
    FALSE A2ML1                             7.635861
    FALSE A2ML2                             6.365863
    FALSE A2ML3                             7.352120
    FALSE A2ML4                             6.472890
    FALSE A4GALT                            7.463812
    FALSE AAAS                              7.490045
    FALSE        g146.blu51_male_pituitary_inc.d3 g20.w106.x_male_pituitary_inc.d3
    FALSE A2ML1                          7.593350                         7.707669
    FALSE A2ML2                          6.312171                         6.574912
    FALSE A2ML3                          7.484586                         7.716466
    FALSE A2ML4                          6.312171                         6.831476
    FALSE A4GALT                         7.668688                         7.429825
    FALSE AAAS                           7.795518                         7.625424
    FALSE        g52.blu58_male_pituitary_bldg g53.y84_male_pituitary_hatch
    FALSE A2ML1                       7.531154                     7.475689
    FALSE A2ML2                       6.155964                     5.995081
    FALSE A2ML3                       7.843931                     7.496785
    FALSE A2ML4                       6.789826                     6.380545
    FALSE A4GALT                      7.538161                     7.155188
    FALSE AAAS                        7.981877                     7.387027
    FALSE        o.s.w.r_male_pituitary_lay o152.o120.w42_male_pituitary_n5
    FALSE A2ML1                    7.459058                        7.460885
    FALSE A2ML2                    6.411186                        6.269404
    FALSE A2ML3                    7.727487                        7.524437
    FALSE A2ML4                    6.456343                        6.518290
    FALSE A4GALT                   7.435533                        7.632468
    FALSE AAAS                     7.664324                        8.008353
    FALSE        o39.y77.x_male_pituitary_hatch o44.blu26.x_male_pituitary_hatch
    FALSE A2ML1                        7.498122                         7.563133
    FALSE A2ML2                        6.322153                         6.237541
    FALSE A2ML3                        7.612237                         7.797236
    FALSE A2ML4                        6.559200                         6.728750
    FALSE A4GALT                       7.336823                         7.557940
    FALSE AAAS                         7.963111                         7.775780
    FALSE        o48.r197.x_male_pituitary_inc.d3 o49.x_male_pituitary_inc.d9
    FALSE A2ML1                          7.544618                    7.579307
    FALSE A2ML2                          6.334213                    6.481641
    FALSE A2ML3                          7.857511                    7.941162
    FALSE A2ML4                          6.691974                    7.214395
    FALSE A4GALT                         7.666304                    7.265418
    FALSE AAAS                           7.900992                    7.894199
    FALSE        o57.g59_male_pituitary_inc.d9 pk.s238.blk.w_male_pituitary_lay
    FALSE A2ML1                       7.329932                         7.511473
    FALSE A2ML2                       6.521213                         6.245718
    FALSE A2ML3                       7.874016                         7.842295
    FALSE A2ML4                       6.704260                         6.301856
    FALSE A4GALT                      7.110012                         7.703197
    FALSE AAAS                        7.883473                         7.815715
    FALSE        pk.w.s141.o_male_pituitary_lay r.s005.pk.blk_male_pituitary_lay
    FALSE A2ML1                        7.426177                         7.592941
    FALSE A2ML2                        6.210765                         6.241928
    FALSE A2ML3                        8.195584                         7.642883
    FALSE A2ML4                        6.506312                         6.421599
    FALSE A4GALT                       7.467079                         7.557062
    FALSE AAAS                         8.063720                         8.056594
    FALSE        r.s059.d.o_male_pituitary_bldg r.s116.blk.pu_male_pituitary_lay
    FALSE A2ML1                        7.541115                         7.519124
    FALSE A2ML2                        6.405451                         6.617811
    FALSE A2ML3                        8.272585                         7.925980
    FALSE A2ML4                        5.995081                         6.478928
    FALSE A4GALT                       7.433560                         7.493644
    FALSE AAAS                         7.749299                         7.571872
    FALSE        r.y.s007.blk_male_pituitary_n9 r176.blu54_male_pituitary_inc.d17
    FALSE A2ML1                        7.752724                          7.446709
    FALSE A2ML2                        6.469053                          6.317122
    FALSE A2ML3                        8.041988                          8.088448
    FALSE A2ML4                        7.646303                          6.580246
    FALSE A4GALT                       7.465944                          7.328835
    FALSE AAAS                         7.690017                          7.902560
    FALSE        r190.o43.x_male_pituitary_lay r195.x_male_pituitary_n9
    FALSE A2ML1                       7.469683                 7.651797
    FALSE A2ML2                       6.161985                 5.995081
    FALSE A2ML3                       7.806797                 7.667945
    FALSE A2ML4                       6.615143                 6.508200
    FALSE A4GALT                      7.516626                 7.555825
    FALSE AAAS                        7.770302                 8.077592
    FALSE        r37.w100.x_male_pituitary_n9 r41.w99.x_male_pituitary_hatch
    FALSE A2ML1                      7.321838                       7.434163
    FALSE A2ML2                      5.995081                       6.271519
    FALSE A2ML3                      7.786998                       8.166031
    FALSE A2ML4                      6.829623                       6.567649
    FALSE A4GALT                     7.443011                       7.297308
    FALSE AAAS                       7.828923                       7.935692
    FALSE        r45.X_male_pituitary_inc.d9 r72.y83.x_male_pituitary_hatch
    FALSE A2ML1                     7.360836                       7.618982
    FALSE A2ML2                     6.199483                       6.360082
    FALSE A2ML3                     7.758236                       7.941807
    FALSE A2ML4                     6.402869                       6.527783
    FALSE A4GALT                    7.608927                       7.529002
    FALSE AAAS                      8.070912                       8.003009
    FALSE        s.pu148.blk.r_male_pituitary_bldg s065.l.d.o_male_pituitary_bldg
    FALSE A2ML1                           7.723316                       7.421766
    FALSE A2ML2                           6.318512                       6.188884
    FALSE A2ML3                           7.913257                       7.859442
    FALSE A2ML4                           6.746657                       6.454823
    FALSE A4GALT                          7.322806                       7.648614
    FALSE AAAS                            7.540650                       7.743094
    FALSE        s066.l.d.r_male_pituitary_bldg s150.w.g.blk_male_pituitary_lay
    FALSE A2ML1                        7.460976                        7.633980
    FALSE A2ML2                        6.205026                        6.371195
    FALSE A2ML3                        7.670010                        8.040253
    FALSE A2ML4                        6.848954                        6.720520
    FALSE A4GALT                       7.218844                        7.277094
    FALSE AAAS                         7.810267                        7.746178
    FALSE        s187.l.o.r_male_pituitary_n9 s243.blk.pk.r_male_pituitary_lay
    FALSE A2ML1                      7.744048                         7.424135
    FALSE A2ML2                      6.498935                         6.267142
    FALSE A2ML3                      8.754564                         7.856219
    FALSE A2ML4                      6.676809                         6.878026
    FALSE A4GALT                     7.310777                         7.331612
    FALSE AAAS                       7.846717                         8.272852
    FALSE        w34.x_male_pituitary_inc.d9 x.blk.blk.ATLAS_male_pituitary_control
    FALSE A2ML1                     7.746976                               7.704155
    FALSE A2ML2                     6.215789                               5.995081
    FALSE A2ML3                     7.821740                               7.505697
    FALSE A2ML4                     6.435219                               6.945692
    FALSE A4GALT                    7.403876                               7.086503
    FALSE AAAS                      8.132265                               7.086503
    FALSE        x.blk16_male_pituitary_n9 x.blu.o.ATLAS_male_pituitary_control
    FALSE A2ML1                   7.680518                             7.438016
    FALSE A2ML2                   6.236310                             6.177326
    FALSE A2ML3                   7.754702                             7.233292
    FALSE A2ML4                   6.650826                             6.838398
    FALSE A4GALT                  7.605105                             7.005899
    FALSE AAAS                    7.780585                             7.398072
    FALSE        x.blu106.o153_male_pituitary_inc.d9.NYNO
    FALSE A2ML1                                  7.705174
    FALSE A2ML2                                  6.366519
    FALSE A2ML3                                  7.953637
    FALSE A2ML4                                  6.423589
    FALSE A4GALT                                 7.511795
    FALSE AAAS                                   7.793696
    FALSE        x.blu117.w89_male_pituitary_inc.d17 x.blu23.w14_male_pituitary_n9
    FALSE A2ML1                             7.513195                      7.667079
    FALSE A2ML2                             5.995081                      6.281281
    FALSE A2ML3                             8.118326                      7.676952
    FALSE A2ML4                             7.064166                      6.582011
    FALSE A4GALT                            7.064166                      7.574275
    FALSE AAAS                              7.513195                      8.079773
    FALSE        x.blu30_male_pituitary_n5 x.blu42.o28_male_pituitary_inc.d3
    FALSE A2ML1                   7.742855                          7.457350
    FALSE A2ML2                   6.315201                          5.995081
    FALSE A2ML3                   7.475941                          7.491622
    FALSE A2ML4                   7.112354                          6.400184
    FALSE A4GALT                  7.560624                          7.417480
    FALSE AAAS                    8.379024                          8.044500
    FALSE        x.g.ATLAS_male_pituitary_control x.g.g.ATLAS_male_pituitary_control
    FALSE A2ML1                          7.377520                           7.991071
    FALSE A2ML2                          5.995081                           5.995081
    FALSE A2ML3                          7.340502                           7.812799
    FALSE A2ML4                          6.549545                           7.760585
    FALSE A4GALT                         7.447946                           7.541550
    FALSE AAAS                           6.842597                           7.957384
    FALSE        x.g.g.g.ATLAS_male_pituitary_control
    FALSE A2ML1                              7.510354
    FALSE A2ML2                              6.199647
    FALSE A2ML3                              7.400913
    FALSE A2ML4                              6.570347
    FALSE A4GALT                             6.917513
    FALSE AAAS                               6.895985
    FALSE        x.g13.w109_male_pituitary_inc.d9 x.g14.w199_male_pituitary_inc.d17
    FALSE A2ML1                          7.536725                          7.666909
    FALSE A2ML2                          6.218961                          6.302088
    FALSE A2ML3                          7.344848                          7.998013
    FALSE A2ML4                          6.382081                          6.275427
    FALSE A4GALT                         7.466818                          7.571834
    FALSE AAAS                           7.813173                          8.095948
    FALSE        x.g147.blu28_male_pituitary_inc.d3 x.g70_male_pituitary_hatch
    FALSE A2ML1                            7.727083                   7.497262
    FALSE A2ML2                            6.455500                   5.995081
    FALSE A2ML3                            7.909041                   7.525073
    FALSE A2ML4                            6.538944                   6.442137
    FALSE A4GALT                           7.558690                   7.296995
    FALSE AAAS                             7.950647                   7.847284
    FALSE        x.o160.w102_male_pituitary_hatch
    FALSE A2ML1                          7.348064
    FALSE A2ML2                          5.995081
    FALSE A2ML3                          7.240377
    FALSE A2ML4                          6.225506
    FALSE A4GALT                         7.240377
    FALSE AAAS                           7.671050
    FALSE        x.o163.w101_male_pituitary_inc.d3.NYNO
    FALSE A2ML1                                7.585690
    FALSE A2ML2                                6.681918
    FALSE A2ML3                                8.414913
    FALSE A2ML4                                6.384991
    FALSE A4GALT                               7.608081
    FALSE AAAS                                 7.491485
    FALSE        x.o164.w123_male_pituitary_n5.NYNO x.o2_male_pituitary_n9
    FALSE A2ML1                            7.489367               7.036877
    FALSE A2ML2                            6.297006               6.382644
    FALSE A2ML3                            7.873537               7.857782
    FALSE A2ML4                            6.686971               7.424599
    FALSE A4GALT                           7.686081               7.617463
    FALSE AAAS                             7.671509               7.943071
    FALSE        x.o30.g134_male_pituitary_bldg x.o47.y82_male_pituitary_inc.d9
    FALSE A2ML1                        7.440118                        7.692658
    FALSE A2ML2                        5.995081                        6.598696
    FALSE A2ML3                        7.629827                        7.959784
    FALSE A2ML4                        6.323055                        6.192198
    FALSE A4GALT                       7.471460                        7.295783
    FALSE AAAS                         7.998085                        7.906919
    FALSE        x.o68_male_pituitary_n5 x.r178_male_pituitary_hatch
    FALSE A2ML1                 7.409372                    7.397835
    FALSE A2ML2                 6.218625                    5.995081
    FALSE A2ML3                 7.858092                    7.833679
    FALSE A2ML4                 7.044919                    6.570867
    FALSE A4GALT                7.603581                    7.380569
    FALSE AAAS                  7.901201                    7.891181
    FALSE        x.r181_male_pituitary_n5 x.r29.w96_male_pituitary_inc.d17
    FALSE A2ML1                  7.354703                         7.600381
    FALSE A2ML2                  6.256645                         6.358466
    FALSE A2ML3                  8.112961                         8.424290
    FALSE A2ML4                  5.995081                         6.759020
    FALSE A4GALT                 7.418321                         7.289881
    FALSE AAAS                   8.089613                         7.743847
    FALSE        x.r64.g140_male_pituitary_inc.d3
    FALSE A2ML1                          7.618476
    FALSE A2ML2                          6.327653
    FALSE A2ML3                          8.032437
    FALSE A2ML4                          6.378820
    FALSE A4GALT                         7.636954
    FALSE AAAS                           8.206331
    FALSE        x.r67.blu35_male_pituitary_bldg.NYNO
    FALSE A2ML1                              7.605040
    FALSE A2ML2                              6.324632
    FALSE A2ML3                              8.058658
    FALSE A2ML4                              7.478025
    FALSE A4GALT                             7.538411
    FALSE AAAS                               7.835523
    FALSE        x.w192.o157_male_pituitary_inc.d9
    FALSE A2ML1                           7.497151
    FALSE A2ML2                           6.331787
    FALSE A2ML3                           7.744546
    FALSE A2ML4                           6.729063
    FALSE A4GALT                          7.641507
    FALSE AAAS                            7.974839
    FALSE        x.y.s.ATLAS_male_pituitary_control
    FALSE A2ML1                            7.155427
    FALSE A2ML2                            5.995081
    FALSE A2ML3                            7.739670
    FALSE A2ML4                            6.676849
    FALSE A4GALT                           6.755596
    FALSE AAAS                             7.108384
    FALSE        x.y132.w76_male_pituitary_inc.d17 x.y141.w116_male_pituitary_inc.d9
    FALSE A2ML1                           7.656321                          7.734291
    FALSE A2ML2                           6.302269                          6.443319
    FALSE A2ML3                           8.047553                          7.772187
    FALSE A2ML4                           7.098782                          6.361553
    FALSE A4GALT                          7.276321                          7.407793
    FALSE AAAS                            7.933366                          7.942652
    FALSE        y129.x_male_pituitary_n9 y131.w185.x_male_pituitary_n9
    FALSE A2ML1                  7.399777                      7.341287
    FALSE A2ML2                  6.486556                      6.306537
    FALSE A2ML3                  7.506163                      8.036004
    FALSE A2ML4                  6.686851                      6.629492
    FALSE A4GALT                 7.436391                      7.408240
    FALSE AAAS                   7.963062                      7.923059
    FALSE        y133.w77.r58_male_pituitary_inc.d17
    FALSE A2ML1                             7.328795
    FALSE A2ML2                             5.995081
    FALSE A2ML3                             8.100614
    FALSE A2ML4                             7.685700
    FALSE A4GALT                            7.449946
    FALSE AAAS                              7.835981
    FALSE        y149.r52.x_male_pituitary_inc.d3 y95.g131.x_male_pituitary_inc.d9
    FALSE A2ML1                          7.620134                         7.923577
    FALSE A2ML2                          6.487003                         6.417196
    FALSE A2ML3                          7.607798                         7.573919
    FALSE A2ML4                          6.720590                         5.995081
    FALSE A4GALT                         7.595329                         7.634817
    FALSE AAAS                           7.941289                         7.476046
    FALSE        y98.o50.x_male_pituitary_inc.d3
    FALSE A2ML1                         7.503331
    FALSE A2ML2                         6.353018
    FALSE A2ML3                         7.994313
    FALSE A2ML4                         6.456394
    FALSE A4GALT                        7.278933
    FALSE AAAS                          7.842642
    FALSE 'data.frame': 6660 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13368 levels "A2ML1","A2ML2",..: 11780 1047 6275 6445 864 6461 6244 5832 3543 11296 ...
    FALSE  $ padj     : num  7.12e-10 5.59e-06 8.07e-07 8.86e-07 1.73e-13 ...
    FALSE  $ logpadj  : num  9.15 5.25 6.09 6.05 12.76 ...
    FALSE  $ lfc      : num  4.31 3.64 3.49 3.44 3.35 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "control","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene         padj   logpadj      lfc      sextissue direction
    FALSE 1        TCTN3 7.117631e-10  9.147665 4.307650 male_pituitary      bldg
    FALSE 2         BDNF 5.587675e-06  5.252769 3.635870 male_pituitary      bldg
    FALSE 3 LOC107050880 8.071344e-07  6.093054 3.486202 male_pituitary      bldg
    FALSE 4 LOC107055116 8.856257e-07  6.052750 3.444661 male_pituitary      bldg
    FALSE 5       ATP2B4 1.726749e-13 12.762771 3.345182 male_pituitary      bldg
    FALSE 6 LOC107055647 1.417354e-20 19.848522 3.299186 male_pituitary      bldg
    FALSE 'data.frame': 0 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13368 levels "A2ML1","A2ML2",..: 
    FALSE  $ padj     : num 
    FALSE  $ logpadj  : num 
    FALSE  $ lfc      : num 
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 
    FALSE  $ direction: Factor w/ 3 levels "bldg","NS","lay": 
    FALSE NULL
    FALSE [1] gene      padj      logpadj   lfc       sextissue direction
    FALSE <0 rows> (or 0-length row.names)
    FALSE 'data.frame': 5 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13368 levels "A2ML1","A2ML2",..: 6951 13307 1019 10527 1219
    FALSE  $ padj     : num  4.18e-02 4.18e-02 4.18e-02 6.40e-02 1.23e-13
    FALSE  $ logpadj  : num  1.38 1.38 1.38 1.19 12.91
    FALSE  $ lfc      : num  4.159 3.591 0.324 -0.671 -18.1
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1
    FALSE  $ direction: Factor w/ 3 levels "lay","NS","inc.d3": 3 3 3 1 1
    FALSE NULL
    FALSE          gene         padj   logpadj         lfc      sextissue direction
    FALSE 1       LRIT3 4.180219e-02  1.378801   4.1591252 male_pituitary    inc.d3
    FALSE 2      ZNF609 4.180219e-02  1.378801   3.5911575 male_pituitary    inc.d3
    FALSE 3       BCAS2 4.180219e-02  1.378801   0.3239632 male_pituitary    inc.d3
    FALSE 4      SEMA6A 6.399499e-02  1.193854  -0.6712511 male_pituitary       lay
    FALSE 5 C11H19ORF40 1.228278e-13 12.910703 -18.0996348 male_pituitary       lay
    FALSE 'data.frame': 1 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13368 levels "A2ML1","A2ML2",..: 1219
    FALSE  $ padj     : num 7.28e-24
    FALSE  $ logpadj  : num 23.1
    FALSE  $ lfc      : num 22.3
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1
    FALSE  $ direction: Factor w/ 3 levels "inc.d3","NS",..: 3
    FALSE NULL
    FALSE          gene         padj  logpadj      lfc      sextissue direction
    FALSE 1 C11H19ORF40 7.284167e-24 23.13762 22.27883 male_pituitary    inc.d9
    FALSE 'data.frame': 432 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13368 levels "A2ML1","A2ML2",..: 9070 609 4632 5857 6853 5551 529 1865 4202 10688 ...
    FALSE  $ padj     : num  7.96e-02 3.63e-06 1.89e-02 5.49e-02 1.61e-08 ...
    FALSE  $ logpadj  : num  1.1 5.44 1.72 1.26 7.79 ...
    FALSE  $ lfc      : num  6.58 6.1 5.87 5.2 5.14 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "inc.d9","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene         padj   logpadj      lfc      sextissue direction
    FALSE 1         PMP2 7.964975e-02  1.098816 6.578339 male_pituitary   inc.d17
    FALSE 2         APOH 3.634437e-06  5.439563 6.102884 male_pituitary   inc.d17
    FALSE 3       HAPLN2 1.888864e-02  1.723799 5.868075 male_pituitary   inc.d17
    FALSE 4 LOC100859224 5.485345e-02  1.260796 5.196808 male_pituitary   inc.d17
    FALSE 5    LOC769726 1.605500e-08  7.794390 5.135049 male_pituitary   inc.d17
    FALSE 6        KPNA2 2.113680e-17 16.674961 4.517901 male_pituitary   inc.d17
    FALSE 'data.frame': 62 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13368 levels "A2ML1","A2ML2",..: 7775 12200 4025 6069 8404 2974 3896 11092 13010 8307 ...
    FALSE  $ padj     : num  0.06457 0.03605 0.00441 0.08897 0.08933 ...
    FALSE  $ logpadj  : num  1.19 1.44 2.36 1.05 1.05 ...
    FALSE  $ lfc      : num  3.052 0.852 0.745 0.567 0.496 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "inc.d17","NS",..: 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE           gene        padj  logpadj       lfc      sextissue direction
    FALSE 1        MYH11 0.064572388 1.189953 3.0522691 male_pituitary     hatch
    FALSE 2     TNFRSF1B 0.036047241 1.443128 0.8518940 male_pituitary     hatch
    FALSE 3          FTL 0.004406171 2.355939 0.7447844 male_pituitary     hatch
    FALSE 4 LOC101750767 0.088970349 1.050755 0.5669874 male_pituitary     hatch
    FALSE 5        OLFM1 0.089328912 1.049008 0.4963497 male_pituitary     hatch
    FALSE 6      DNAJB11 0.088822923 1.051475 0.4557011 male_pituitary     hatch
    FALSE 'data.frame': 150 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13368 levels "A2ML1","A2ML2",..: 8372 4969 7750 6633 200 13133 1935 10006 265 3231 ...
    FALSE  $ padj     : num  0.0666 0.068 0.0587 0.0863 0.02 ...
    FALSE  $ logpadj  : num  1.18 1.17 1.23 1.06 1.7 ...
    FALSE  $ lfc      : num  1.75 1.63 1.59 1.58 1.42 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "hatch","NS","n5": 3 3 3 3 3 3 3 3 3 3 ...
    FALSE NULL
    FALSE        gene       padj  logpadj      lfc      sextissue direction
    FALSE 1      OASL 0.06655414 1.176825 1.751298 male_pituitary        n5
    FALSE 2     IFIT5 0.06800865 1.167436 1.630534 male_pituitary        n5
    FALSE 3       MX1 0.05869716 1.231383 1.585976 male_pituitary        n5
    FALSE 4 LOC418892 0.08626582 1.064161 1.579719 male_pituitary        n5
    FALSE 5   ADAMTS1 0.01996699 1.699687 1.415549 male_pituitary        n5
    FALSE 6    ZBTB16 0.02386156 1.622301 1.020738 male_pituitary        n5
    FALSE 'data.frame': 65 obs. of  6 variables:
    FALSE  $ gene     : Factor w/ 13368 levels "A2ML1","A2ML2",..: 29 7225 8020 668 5561 8856 12179 1333 13010 9039 ...
    FALSE  $ padj     : num  2.68e-05 6.21e-02 2.86e-02 3.06e-02 1.10e-02 ...
    FALSE  $ logpadj  : num  4.57 1.21 1.54 1.51 1.96 ...
    FALSE  $ lfc      : num  3.3 2 1.42 1.31 1.28 ...
    FALSE  $ sextissue: Factor w/ 1 level "male_pituitary": 1 1 1 1 1 1 1 1 1 1 ...
    FALSE  $ direction: Factor w/ 3 levels "n5","NS","n9": 3 3 3 3 3 3 3 3 1 1 ...
    FALSE NULL
    FALSE       gene         padj  logpadj      lfc      sextissue direction
    FALSE 1    ABCA4 0.0000268105 4.571695 3.300416 male_pituitary        n9
    FALSE 2    MATN1 0.0620882642 1.206990 1.999787 male_pituitary        n9
    FALSE 3    NETO2 0.0285560577 1.544302 1.416218 male_pituitary        n9
    FALSE 4 ARHGAP40 0.0306263540 1.513905 1.305778 male_pituitary        n9
    FALSE 5    KRT14 0.0110268135 1.957550 1.276082 male_pituitary        n9
    FALSE 6   PHLDA2 0.0419550480 1.377216 0.937829 male_pituitary        n9

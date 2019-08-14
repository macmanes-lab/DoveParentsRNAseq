Based on a discussion on Twitter, <a href="https://twitter.com/Jared_Mamrot/status/1161047070821646336" class="uri">https://twitter.com/Jared_Mamrot/status/1161047070821646336</a>, I decided to try out the R package `topconfects` to test the effect size.
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    library(DESeq2)
    library(tidyverse)
    library(topconfects)

    knitr::opts_chunk$set(fig.path = '../figures/topconfects/', cache = TRUE, message = F, echo = T)

This is a demo from the vignette showing how the tool works.
------------------------------------------------------------

    # Confident log2 fold changes based on a DESeq2 analysis

    # Generate some random data
    n <- 20
    folds <- seq(-8,8,length.out=n)
    row_means <- runif(n, min=0, max=5)
    lib_scale <- c(1,2,3,4)
    means <- 2^(outer(folds, c(-0.5,-0.5,0.5,0.5))) *
      row_means * rep(lib_scale,each=n)
    counts <- rnbinom(length(means), mu=means, size=1/0.1)
    dim(counts) <- dim(means)
    group <- factor(c("A","A","B","B"))
    # Apply DESeq2
    dds <- DESeqDataSetFromMatrix(
      countData = counts,
      colData = data.frame(group=group),
      design = ~group)
    dds <- DESeq(dds)
    # Find top confident effect sizes
    deseq2_confects(dds, name="group_B_vs_A", step=0.1)

    ## $table
    ##  rank index confect effect    baseMean  filtered
    ##   1   18     5.6    10.277392 93.588768 FALSE   
    ##   2   19     4.5     7.828128 94.517448 FALSE   
    ##   3   20     4.2     8.699898 31.525759 FALSE   
    ##   4   17     4.2     8.225876 60.944194 FALSE   
    ##   5   15     2.5     5.194451 38.460041 FALSE   
    ##   6   16     2.5     5.954003 25.077792 FALSE   
    ##   7    2    -1.8    -5.896977 19.872687 FALSE   
    ##   8   14     1.2     3.935115 22.443317 FALSE   
    ##   9    3    -0.7    -5.977252  7.319253 FALSE   
    ##  10    4    -0.4    -4.161852 11.908750 FALSE   
    ## ...
    ## 11 of 20 non-zero effect size at FDR 0.05

    #deseq2_confects(dds, name="group_A_vs_B", step=0.1) # Error: subscript contains invalid names

Now I import the gene expression data and order the treatment according the the parentatl care cycle
----------------------------------------------------------------------------------------------------

    # import "colData" which contains sample information and "countData" which contains read counts
    c.colData <- read.csv("../metadata/00_colData_characterization.csv", header = T, row.names = 1)
    c.countData <- read.csv("../results/00_countData_characterization.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)

    # create grouping
    c.colData$sextissue <- as.factor(paste(c.colData$sex, c.colData$tissue, sep = "_"))

    c.colData$treatment <- factor(c.colData$treatment, levels = 
                                  c("control",  "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", "hatch", "n5", "n9"))
    levels(c.colData$treatment)

    ## [1] "control" "bldg"    "lay"     "inc.d3"  "inc.d9"  "inc.d17" "hatch"  
    ## [8] "n5"      "n9"

This is a function to subset the data to look compare one transition (eg. incubation day 17 to lay) in a single “sex tissue” group (e.g. female pitutiary)
----------------------------------------------------------------------------------------------------------------------------------------------------------

    subsetnconfects <- function(mygroup, mytreatments, mycomparison){
      # subset col and count data
      colData <- c.colData %>%
          dplyr::filter(sextissue == mygroup) %>%
          dplyr::filter(treatment %in% mytreatments) %>%
          droplevels()
      row.names(colData) <- colData$V1
      
      # which counts to save
      savecols <- as.character(colData$V1) 
      savecols <- as.vector(savecols) 
      
      # save counts that match colData
      countData <- c.countData %>% dplyr::select(one_of(savecols)) 
      
      # check that row and col lenghts are equal
      #print(ncol(countData) == nrow(colData))

      # run deseq
      dds <- DESeqDataSetFromMatrix(
        countData = countData,
       colData = colData,
       design = ~treatment)
      dds <- DESeq(dds)
      res <- results(dds)
      print(summary(res))
      # Find top confident effect sizes
      deseq2_confects(dds, name=mycomparison, step=0.1)
    }

Number of differentially expressed genes for the building versus lay comparison
-------------------------------------------------------------------------------

    subsetnconfects("female_hypothalamus", c("bldg", "lay"), "treatment_lay_vs_bldg") #1

    ## 
    ## out of 14817 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0067%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    ## $table
    ##  rank index confect effect     baseMean   name           filtered
    ##   1    8051  0       0.5090839  141.65463 XP_015133392.1 FALSE   
    ##   2    8452 NA       0.4517062  665.00496 XP_015135313.1 FALSE   
    ##   3    9246 NA       0.4600936  103.39163 XP_015139349.1 FALSE   
    ##   4   14937 NA      -1.3775434   22.16119 XP_430508.3    FALSE   
    ##   5   13830 NA       0.2317874 2144.35730 XP_419512.5    FALSE   
    ##   6    8040 NA      -0.7834858  129.46254 XP_015133325.1 FALSE   
    ##   7   12389 NA       0.7361647   67.15191 XP_015155243.1 FALSE   
    ##   8    9984 NA       1.0736861   15.87459 XP_015143094.1 FALSE   
    ##   9    4205 NA      -0.8493417 5746.71731 NP_990820.1    FALSE   
    ##  10    8461 NA       0.8839854  163.66689 XP_015135349.1 FALSE   
    ## ...
    ## 1 of 14937 non-zero effect size at FDR 0.05

    subsetnconfects("male_hypothalamus", c("bldg", "lay"), "treatment_lay_vs_bldg") #0

    ## 
    ## out of 14800 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 1, 0.0068%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 2, 0.014%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    ## $table
    ##  rank index confect effect     baseMean    name           filtered
    ##   1   13665 NA      -5.1835827    7.261780 XP_418694.3    FALSE   
    ##   2    7598 NA       0.6378134   94.208775 XP_015131103.1 FALSE   
    ##   3    8030 NA       1.0638550   12.277322 XP_015133278.1 FALSE   
    ##   4   11408 NA       0.5737248   48.594097 XP_015150420.1 FALSE   
    ##   5   12007 NA      -1.0561058   52.829713 XP_015153385.1 FALSE   
    ##   6    5120 NA       1.5248430    5.900277 XP_003642569.2 FALSE   
    ##   7    3814 NA      -1.8412571    6.268012 NP_990213.1    FALSE   
    ##   8    3171 NA       0.4483693   57.851503 NP_001291963.1 FALSE   
    ##   9    4980 NA       0.2469132 1339.068540 XP_003641807.3 FALSE   
    ##  10   10346 NA      -0.2980188  152.930426 XP_015145107.1 FALSE   
    ## ...
    ## 0 of 14937 non-zero effect size at FDR 0.05

    subsetnconfects("female_pituitary", c("bldg", "lay"), "treatment_lay_vs_bldg") #455

    ## 
    ## out of 14797 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 267, 1.8%
    ## LFC < 0 (down)     : 390, 2.6%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 1435, 9.7%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    ## $table
    ##  rank index confect effect    baseMean name           filtered
    ##   1    2672  3.3     6.613285 216.1745 NP_001264456.1 FALSE   
    ##   2   11625  2.2     5.878851 187.6616 XP_015151573.1 FALSE   
    ##   3    4211 -2.2    -5.388463 166.4979 NP_990826.1    FALSE   
    ##   4    3051 -1.6    -3.505873 451.3353 NP_001278710.1 FALSE   
    ##   5    7733 -1.5    -3.084037 145.2674 XP_015131806.1 FALSE   
    ##   6    7112 -1.2    -2.881740 338.9900 XP_015129004.1 FALSE   
    ##   7    3702 -0.9    -2.576023 397.6196 NP_990049.1    FALSE   
    ##   8    1964  0.9     2.465130 161.1175 NP_001135726.1 FALSE   
    ##   9    3937 -0.8    -2.986388 302.5164 NP_990392.1    FALSE   
    ##  10   11767 -0.7    -3.207275 315.1872 XP_015152260.1 FALSE   
    ## ...
    ## 455 of 14937 non-zero effect size at FDR 0.05

    subsetnconfects("male_pituitary", c("bldg", "lay"), "treatment_lay_vs_bldg") #0

    ## 
    ## out of 14801 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 2, 0.014%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 1, 0.0068%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    ## $table
    ##  rank index confect effect     baseMean   name           filtered
    ##   1   13884 NA      -0.3853415  273.60812 XP_419733.1    FALSE   
    ##   2    4132 NA      -0.7369055   63.39692 NP_990696.1    FALSE   
    ##   3    5316 NA      -0.2059543 4984.51144 XP_004934572.1 FALSE   
    ##   4    4111 NA       1.3786842   17.22548 NP_990666.1    FALSE   
    ##   5   12259 NA      -0.9969856  270.15240 XP_015154663.1 FALSE   
    ##   6    2753 NA       0.4334316  161.39445 NP_001264589.1 FALSE   
    ##   7    7681 NA      -0.5595780  115.95246 XP_015131483.1 FALSE   
    ##   8    6221 NA       0.4632017 1470.28165 XP_004942634.1 FALSE   
    ##   9    1274 NA       0.2526556  823.14180 NP_001026358.1 FALSE   
    ##  10    5901 NA      -0.2260285 1834.49773 XP_004939662.1 FALSE   
    ## ...
    ## 0 of 14937 non-zero effect size at FDR 0.05

    subsetnconfects("female_gonad", c("bldg", "lay"), "treatment_lay_vs_bldg") #452

    ## 
    ## out of 14874 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 204, 1.4%
    ## LFC < 0 (down)     : 555, 3.7%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 867, 5.8%
    ## (mean count < 2)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    ## $table
    ##  rank index confect effect    baseMean  name           filtered
    ##   1    9849  1.7     9.751777 222.02716 XP_015142231.1 FALSE   
    ##   2   14465  1.6     5.455372 173.44431 XP_423478.5    FALSE   
    ##   3    5009 -1.2    -3.027448 316.98811 XP_003641938.1 FALSE   
    ##   4    7963 -1.0    -3.658124 335.68037 XP_015132863.1 FALSE   
    ##   5    3436 -0.9    -4.140441  64.75787 NP_989665.1    FALSE   
    ##   6    6437 -0.7    -2.406599  67.29939 XP_004944800.1 FALSE   
    ##   7   11581  0.6     3.767326  33.61142 XP_015151399.1 FALSE   
    ##   8    9015 -0.6    -2.253199 378.46138 XP_015138239.1 FALSE   
    ##   9    8028 -0.4    -1.803896  46.72208 XP_015133274.1 FALSE   
    ##  10    1131 -0.4    -1.859466  99.48012 NP_001026099.1 FALSE   
    ## ...
    ## 452 of 14937 non-zero effect size at FDR 0.05

    subsetnconfects("male_gonad", c("bldg", "lay"), "treatment_lay_vs_bldg") #26

    ## 
    ## out of 14891 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 51, 0.34%
    ## LFC < 0 (down)     : 54, 0.36%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 8082, 54%
    ## (mean count < 183)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    ## $table
    ##  rank index confect effect     baseMean  name           filtered
    ##   1    2066 0       -1.1857223  314.1953 NP_001165148.1 FALSE   
    ##   2    3137 0       -0.5322666  681.6908 NP_001289088.1 FALSE   
    ##   3    8911 0       -0.6878685 2614.1646 XP_015137811.1 FALSE   
    ##   4    3949 0       -0.9612694  957.4182 NP_990411.1    FALSE   
    ##   5    2468 0       -0.6593242  561.5771 NP_001244223.1 FALSE   
    ##   6    2914 0       -0.6662434 1055.6694 NP_001264941.1 FALSE   
    ##   7    3357 0       -0.6557910 1059.5675 NP_989548.1    FALSE   
    ##   8    7241 0       -0.6613696  547.7171 XP_015129713.1 FALSE   
    ##   9     564 0       -0.5653406 1065.4216 NP_001007968.1 FALSE   
    ##  10   10843 0        0.3927802  255.6496 XP_015147761.1 FALSE   
    ## ...
    ## 26 of 14937 non-zero effect size at FDR 0.05

Number of differentially expressed genes for the incubation day 17 versus hatch comparison
------------------------------------------------------------------------------------------

    subsetnconfects("female_hypothalamus", c("inc.d17", "hatch"), "treatment_hatch_vs_inc.d17") #0

    ## 
    ## out of 14838 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 0, 0%
    ## LFC < 0 (down)     : 0, 0%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 3, 0.02%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    ## $table
    ##  rank index confect effect     baseMean    name           filtered
    ##   1   11445 NA       0.2986722   318.71001 XP_015150747.1 FALSE   
    ##   2    4638 NA       0.6241912   147.04403 XP_001234376.2 FALSE   
    ##   3    1317 NA       1.0541375   326.00853 NP_001026427.2 FALSE   
    ##   4      81 NA       0.8197472    70.02199 NP_001004373.1 FALSE   
    ##   5    5853 NA       0.3901112   420.69883 XP_004939153.1 FALSE   
    ##   6   10245 NA       0.5492389    77.54440 XP_015144479.1 FALSE   
    ##   7   13587 NA       0.2630143   499.68449 XP_418192.3    FALSE   
    ##   8    4209 NA       0.4167797 18187.18333 NP_990824.1    FALSE   
    ##   9    4063 NA      -0.7569964   139.93766 NP_990593.1    FALSE   
    ##  10    6043 NA       0.3177381   292.45875 XP_004940989.1 FALSE   
    ## ...
    ## 0 of 14937 non-zero effect size at FDR 0.05

    subsetnconfects("male_hypothalamus", c("inc.d17", "hatch"), "treatment_hatch_vs_inc.d17") #3

    ## 
    ## out of 14821 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 3, 0.02%
    ## LFC < 0 (down)     : 4, 0.027%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 2, 0.013%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    ## $table
    ##  rank index confect effect     baseMean  name           filtered
    ##   1   14421  0      -0.7988301 322.53048 XP_422744.3    FALSE   
    ##   2   10336  0       0.6831391  60.27952 XP_015145039.1 FALSE   
    ##   3    5027  0      -0.6457629 771.60343 XP_003642058.2 FALSE   
    ##   4    6493 NA       0.3542273 818.17672 XP_004945358.2 FALSE   
    ##   5    4334 NA      -1.3293121 158.35937 XP_001231711.1 FALSE   
    ##   6    7061 NA      -0.8847096  79.17451 XP_015128714.1 FALSE   
    ##   7   12582 NA       1.4148943  10.49862 XP_015155945.1 FALSE   
    ##   8   12654 NA      -1.0234398  68.18624 XP_015156815.1 FALSE   
    ##   9    9298 NA       0.2513510 961.75243 XP_015139586.1 FALSE   
    ##  10      58 NA      -0.8447907  21.64548 NP_001001759.1 FALSE   
    ## ...
    ## 3 of 14937 non-zero effect size at FDR 0.05

    subsetnconfects("female_pituitary", c("inc.d17", "hatch"), "treatment_hatch_vs_inc.d17") #6

    ## 
    ## out of 14758 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 3, 0.02%
    ## LFC < 0 (down)     : 10, 0.068%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 6, 0.041%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    ## $table
    ##  rank index confect effect     baseMean   name           filtered
    ##   1   12794  0      -3.9132061   10.29010 XP_015158199.1 FALSE   
    ##   2    4211  0      -3.0705493   45.90387 NP_990826.1    FALSE   
    ##   3    7582  0      -0.9026169 1243.32775 XP_015131019.1 FALSE   
    ##   4   13628  0       0.4082338 3163.39685 XP_418430.2    FALSE   
    ##   5    7147  0      -1.0748571  153.80684 XP_015129182.1 FALSE   
    ##   6    3308  0       0.6250949   56.70568 NP_989480.2    FALSE   
    ##   7   10831 NA      -1.2532044   35.62841 XP_015147738.1 FALSE   
    ##   8    3679 NA      -0.8955422   69.78917 NP_990018.1    FALSE   
    ##   9   11675 NA      -1.1822734   18.71792 XP_015151838.1 FALSE   
    ##  10    6815 NA      -0.8285358 2701.99646 XP_004949112.1 FALSE   
    ## ...
    ## 6 of 14937 non-zero effect size at FDR 0.05

    subsetnconfects("male_pituitary", c("inc.d17", "hatch"), "treatment_hatch_vs_inc.d17") #24

    ## 
    ## out of 14749 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 12, 0.081%
    ## LFC < 0 (down)     : 42, 0.28%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 859, 5.8%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    ## $table
    ##  rank index confect effect    baseMean   name           filtered
    ##   1    8105 -2.1    -7.979165 536.925795 XP_015133687.1 FALSE   
    ##   2    4359 -0.4    -4.840091 161.475315 XP_001231917.1 FALSE   
    ##   3    4072  0.0     3.050371  33.243461 NP_990605.2    FALSE   
    ##   4    2900  0.0    -5.483364  40.864112 NP_001264923.1 FALSE   
    ##   5    6349  0.0    -1.455984  15.172086 XP_004943912.1 FALSE   
    ##   6    8908  0.0    -2.163509   4.100100 XP_015137776.1 FALSE   
    ##   7    7249  0.0    -1.969774  60.775808 XP_015129737.1 FALSE   
    ##   8   12108  0.0    -4.948164   7.574359 XP_015153905.1 FALSE   
    ##   9     889  0.0    -2.780820  55.935341 NP_001025732.1 FALSE   
    ##  10    9868  0.0    -2.996540  39.851586 XP_015142329.1 FALSE   
    ## ...
    ## 24 of 14937 non-zero effect size at FDR 0.05

    subsetnconfects("female_gonad", c("inc.d17", "hatch"), "treatment_hatch_vs_inc.d17") #3

    ## 
    ## out of 14853 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0067%
    ## LFC < 0 (down)     : 5, 0.034%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 1, 0.0067%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    ## $table
    ##  rank index confect effect     baseMean   name           filtered
    ##   1   12794 -0.3    -3.8799714  10.090834 XP_015158199.1 FALSE   
    ##   2    3061 -0.3    -3.0209818 277.182544 NP_001280019.1 FALSE   
    ##   3    3536  0.1     1.4666278  38.829950 NP_989812.1    FALSE   
    ##   4    9663   NA    -2.4884139   4.883027 XP_015141299.1 FALSE   
    ##   5   13497   NA    -2.9601715  13.222331 XP_417666.1    FALSE   
    ##   6    8854   NA    -0.5970754 215.921229 XP_015137498.1 FALSE   
    ##   7    2680   NA    -2.6210892  11.064269 NP_001264473.1 FALSE   
    ##   8     783   NA     5.2687842   3.455082 NP_001012920.1 FALSE   
    ##   9    4884   NA    -1.1374606  11.781145 XP_003641087.2 FALSE   
    ##  10     330   NA    -1.1050912  26.208385 NP_001006368.1 FALSE   
    ## ...
    ## 3 of 14937 non-zero effect size at FDR 0.05

    subsetnconfects("male_gonad", c("inc.d17", "hatch"), "treatment_hatch_vs_inc.d17") #1

    ## 
    ## out of 14875 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1, 0.0067%
    ## LFC < 0 (down)     : 1, 0.0067%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results
    ## 
    ## NULL

    ## $table
    ##  rank index confect effect     baseMean    name           filtered
    ##   1      67 0.1      1.3661866   30.635531 NP_001001777.1 FALSE   
    ##   2    5810  NA     -0.4156240 1003.456522 XP_004938804.1 FALSE   
    ##   3    9618  NA     -1.1100643   21.766250 XP_015141099.1 FALSE   
    ##   4   12191  NA      1.5424894    7.788929 XP_015154284.1 FALSE   
    ##   5    6029  NA     -0.4729438  308.840071 XP_004940903.1 FALSE   
    ##   6    9580  NA     -1.0491934   40.158760 XP_015140894.1 FALSE   
    ##   7   13656  NA     -0.9690712   19.915158 XP_418613.2    FALSE   
    ##   8    3047  NA      0.5787691  137.934279 NP_001278581.1 FALSE   
    ##   9   11038  NA     -0.9618028   41.722236 XP_015148692.1 FALSE   
    ##  10    4709  NA     -0.2649075  428.264959 XP_001235056.1 FALSE   
    ## ...
    ## 1 of 14937 non-zero effect size at FDR 0.05

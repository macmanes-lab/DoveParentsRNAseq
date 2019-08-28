    library(tidyverse)
    library(DESeq2)
    library(cowplot)
    library(RColorBrewer)
    library(pheatmap)
    library(kableExtra)
    library(viridis)
    library(forcats)
    library("wesanderson")
    library(caret) # LDA analysis
    library(MASS) # LDA analysis

    library("BiocParallel")
    register(MulticoreParam(4))

    source("../R/functions.R")  # load custom functions 
    source("../R/themes.R")  # load custom themes and color palletes

    knitr::opts_chunk$set(fig.path = '../figures/sexes/', cache = TRUE)

DEseq2 analysis on characterization. Looking at treatment AND sex for each tissue
---------------------------------------------------------------------------------

    # import "colData" which contains sample information and "countData" which contains read counts
    c.colData <- read.csv("../metadata/00_colData_characterization.csv", header = T, row.names = 1)
    c.countData <- read.csv("../results/00_countData_characterization.csv", header = T, row.names = 1)
    geneinfo <- read.csv("../metadata/00_geneinfo.csv", row.names = 1)

    # set levels
    c.colData$treatment <- factor(c.colData$treatment, levels = 
                                  c("control",  "bldg", "lay", "inc.d3", "inc.d9", "inc.d17", "hatch", "n5", "n9"))
    levels(c.colData$treatment)

    ## [1] "control" "bldg"    "lay"     "inc.d3"  "inc.d9"  "inc.d17" "hatch"  
    ## [8] "n5"      "n9"

    c.colData$sextissue <- as.factor(paste(c.colData$sex, c.colData$tissue, sep = "_"))

    summary(c.colData[c(7,3,4,5,8)])

    ##              study         sex               tissue      treatment  
    ##  charcterization:576   female:289   gonad       :194   control: 73  
    ##                        male  :287   hypothalamus:189   inc.d9 : 71  
    ##                                     pituitary   :193   inc.d17: 66  
    ##                                                        n9     : 66  
    ##                                                        bldg   : 60  
    ##                                                        lay    : 60  
    ##                                                        (Other):180  
    ##                sextissue 
    ##  female_gonad       :98  
    ##  female_hypothalamus:95  
    ##  female_pituitary   :96  
    ##  male_gonad         :96  
    ##  male_hypothalamus  :94  
    ##  male_pituitary     :97  
    ## 

    c.colData <- c.colData %>%
        mutate(hypothesis = fct_recode(treatment,
                                "anticipation" = "control",
                                "anticipation" = "bldg",
                                "incubation" = "lay",
                                "incubation" = "inc.d3",
                                "incubation" = "inc.d9",
                                "incubation" = "inc.d17",
                                "hatchling.care" = "hatch",
                                "hatchling.care" = "n5",
                                "hatchling.care" = "n9"))

    dds.hypothalamus <- subsetDESeq2(c.colData,  c.countData, c("female_hypothalamus","male_hypothalamus") )

    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 189 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(189): L.Blu13_male_hypothalamus_control.NYNO
    ##   L.G107_male_hypothalamus_control ...
    ##   y97.x_female_hypothalamus_n9 y98.o50.x_male_hypothalamus_inc.d3
    ## colData names(8): V1 bird ... study sextissue
    ## [1] 14597   189

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 4 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 4 workers

    ## -- replacing outliers and refitting for 6 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    dds.pituitary <- subsetDESeq2(c.colData,  c.countData, c("female_pituitary", "male_pituitary" ) )

    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 193 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(193): L.Blu13_male_pituitary_control.NYNO
    ##   L.G107_male_pituitary_control ... y97.x_female_pituitary_n9
    ##   y98.o50.x_male_pituitary_inc.d3
    ## colData names(8): V1 bird ... study sextissue
    ## [1] 14484   193

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 4 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 4 workers

    ## -- replacing outliers and refitting for 28 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    dds.gonad <- subsetDESeq2(c.colData,  c.countData, c("female_gonad", "male_gonad") )

    ## [1] TRUE
    ## class: DESeqDataSet 
    ## dim: 14937 194 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(14937): NP_001001127.1 NP_001001129.1 ... XP_430449.2
    ##   XP_430508.3
    ## rowData names(0):
    ## colnames(194): L.Blu13_male_gonad_control.NYNO
    ##   L.G107_male_gonad_control ... y97.x_female_gonad_n9
    ##   y98.o50.x_male_gonad_inc.d3
    ## colData names(8): V1 bird ... study sextissue
    ## [1] 14843   194

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates: 4 workers

    ## mean-dispersion relationship

    ## final dispersion estimates, fitting model and testing: 4 workers

    ## -- replacing outliers and refitting for 8 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

PCA
---

    hypPCA <- returnPCAs2(dds.hypothalamus)

    ## [1] 20 10  7  6  3  3
    ##                Df Sum Sq Mean Sq F value Pr(>F)    
    ## treatment       8   4783   597.8  17.824 <2e-16 ***
    ## sex             1     45    44.7   1.333  0.250    
    ## treatment:sex   8     53     6.6   0.197  0.991    
    ## Residuals     171   5735    33.5                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment       8   1657  207.09  10.524 5.86e-12 ***
    ## sex             1     33   33.14   1.684    0.196    
    ## treatment:sex   8    119   14.93   0.759    0.640    
    ## Residuals     171   3365   19.68                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment       8    185   23.16   1.128  0.347
    ## sex             1     24   23.82   1.160  0.283
    ## treatment:sex   8    147   18.38   0.895  0.522
    ## Residuals     171   3510   20.53               
    ##                Df Sum Sq Mean Sq F value Pr(>F)    
    ## treatment       8  137.2    17.2   1.756 0.0890 .  
    ## sex             1 1064.2  1064.2 108.914 <2e-16 ***
    ## treatment:sex   8  154.8    19.4   1.981 0.0515 .  
    ## Residuals     171 1670.8     9.8                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    plotPC12(hypPCA, "female and male hypothalamus")

![](../figures/sexes/pca-1.png)

    pitPCA <- returnPCAs2(dds.pituitary)

    ## [1] 13  9  7  5  4  4
    ##                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment       8   3851   481.4  49.296  < 2e-16 ***
    ## sex             1    627   626.6  64.172 1.54e-13 ***
    ## treatment:sex   8     60     7.5   0.771    0.628    
    ## Residuals     175   1709     9.8                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value Pr(>F)    
    ## treatment       8 2276.0   284.5  45.986 <2e-16 ***
    ## sex             1 1007.5  1007.5 162.847 <2e-16 ***
    ## treatment:sex   8   91.1    11.4   1.842 0.0723 .  
    ## Residuals     175 1082.7     6.2                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value Pr(>F)    
    ## treatment       8 1123.8   140.5   47.77 <2e-16 ***
    ## sex             1 1979.9  1979.9  673.36 <2e-16 ***
    ## treatment:sex   8   36.0     4.5    1.53   0.15    
    ## Residuals     175  514.6     2.9                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value Pr(>F)  
    ## treatment       8  219.9   27.49   2.278 0.0242 *
    ## sex             1   72.6   72.56   6.013 0.0152 *
    ## treatment:sex   8   86.7   10.84   0.898 0.5192  
    ## Residuals     175 2111.7   12.07                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    plotPC12(pitPCA, "female and male pituitary")

![](../figures/sexes/pca-2.png)

    gonPCA <- returnPCAs2(dds.gonad)

    ## [1] 91  3  1  1  0  0
    ##                Df Sum Sq Mean Sq   F value   Pr(>F)    
    ## treatment       8    314      39     3.398 0.001182 ** 
    ## sex             1 433334  433334 37453.286  < 2e-16 ***
    ## treatment:sex   8    334      42     3.606 0.000662 ***
    ## Residuals     176   2036      12                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## treatment       8   2349  293.60   5.181 8.06e-06 ***
    ## sex             1      7    7.32   0.129     0.72    
    ## treatment:sex   8   2340  292.53   5.162 8.50e-06 ***
    ## Residuals     176   9974   56.67                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value  Pr(>F)    
    ## treatment       8    914  114.31   3.879 0.00031 ***
    ## sex             1      1    0.81   0.028 0.86842    
    ## treatment:sex   8    655   81.88   2.779 0.00643 ** 
    ## Residuals     176   5187   29.47                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ##                Df Sum Sq Mean Sq F value Pr(>F)  
    ## treatment       8  208.4  26.047   1.994 0.0496 *
    ## sex             1    1.5   1.526   0.117 0.7329  
    ## treatment:sex   8  173.2  21.655   1.658 0.1117  
    ## Residuals     176 2298.6  13.060                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    plotPC12(gonPCA, "female and male gonad")

![](../figures/sexes/pca-3.png)

Linear discriminant analysis (LDA)
----------------------------------

<a href="http://www.sthda.com/english/articles/36-classification-methods-essentials/146-discriminant-analysis-essentials-in-r/" class="uri">http://www.sthda.com/english/articles/36-classification-methods-essentials/146-discriminant-analysis-essentials-in-r/</a>

    theme_set(theme_classic())

    LDAanalysis <- function(mydds){
      # data prep
      vsd <- vst(mydds, blind=FALSE) # variance stabilized 
      tvsd <- as.data.frame(t(assay(vsd)))   # get vsd counts and transform

      tvsd <- tvsd[, sample(ncol(tvsd), 20)] # select 20 random genes
      tvsd$V1 <- row.names(tvsd)             # prep for joining

      colDataPit <- subsetcolData2(c.colData, c("female_pituitary", "male_pituitary"))

      mydata <- left_join(tvsd, colDataPit, by = "V1")  # merge counts and variables
      

      # Split the data into training (80%) and test set (20%)
      set.seed(175)

      training.samples <- mydata$treatment %>%
        createDataPartition(p = 0.8, list = FALSE)
      train.data <- mydata[training.samples, ]
      test.data <- mydata[-training.samples, ]

      # Normalize the data. Categorical variables are automatically ignored.
      # Estimate preprocessing parameters
      preproc.param <- train.data %>% 
        preProcess(method = c("center", "scale"))

      # Transform the data using the estimated parameters
      train.transformed <- preproc.param %>% predict(train.data)
      test.transformed <- preproc.param %>% predict(test.data)

      # LDA analysis

      # Fit the model
      print(model <- lda(treatment~ sex, data = train.transformed))
      # Make predictions
      print(predictions <- model %>% predict(test.transformed))
    # Model accuracy
      print(mean(predictions$class==test.transformed$treatment))

      # results
      print(model)

      # make predictions
      predictions <- model %>% predict(test.transformed)
      
      # Predicted classes
      print(head(predictions$class, 6))
      # Predicted probabilities of class memebership.
      print(head(predictions$posterior, 6) )
      # Linear discriminants
      print(head(predictions$x, 3) )

      # LDA plot using ggplot2 
      lda.data <- cbind(train.transformed, predict(model)$x)
      p<-ggplot(data=lda.data, aes(x=treatment, y=LD1, fill = treatment)) +
        geom_bar( stat = "summary", fun.y = "mean")
      plot(p)

      # model accuracy 
      print(mean(predictions$class==test.transformed$treatment))
    }

    LDAanalysis(dds.pituitary)

    ## Warning: Column `V1` joining character vector and factor, coercing into
    ## character vector

    ## Call:
    ## lda(treatment ~ sex, data = train.transformed)
    ## 
    ## Prior probabilities of groups:
    ##   control      bldg       lay    inc.d3    inc.d9   inc.d17     hatch 
    ## 0.1282051 0.1025641 0.1025641 0.1025641 0.1282051 0.1153846 0.1025641 
    ##        n5        n9 
    ## 0.1025641 0.1153846 
    ## 
    ## Group means:
    ##           sexmale
    ## control 0.5000000
    ## bldg    0.5625000
    ## lay     0.5625000
    ## inc.d3  0.5000000
    ## inc.d9  0.4500000
    ## inc.d17 0.4444444
    ## hatch   0.4375000
    ## n5      0.5000000
    ## n9      0.5555556
    ## 
    ## Coefficients of linear discriminants:
    ##              LD1
    ## sexmale 1.950186
    ## $class
    ##  [1] control control inc.d9  inc.d9  control control control inc.d9 
    ##  [9] inc.d9  control control control control inc.d9  inc.d9  control
    ## [17] control control inc.d9  inc.d9  inc.d9  inc.d9  inc.d9  control
    ## [25] control inc.d9  control control inc.d9  inc.d9  control inc.d9 
    ## [33] inc.d9  control inc.d9  control inc.d9 
    ## Levels: control bldg lay inc.d3 inc.d9 inc.d17 hatch n5 n9
    ## 
    ## $posterior
    ##       control       bldg        lay    inc.d3    inc.d9   inc.d17
    ## 5   0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 8   0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 10  0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 15  0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 16  0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 22  0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 25  0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 30  0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 40  0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 42  0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 45  0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 46  0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 50  0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 60  0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 64  0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 66  0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 68  0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 83  0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 88  0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 97  0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 105 0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 106 0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 108 0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 120 0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 126 0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 128 0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 130 0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 131 0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 146 0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 160 0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 164 0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 168 0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 174 0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 178 0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 179 0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 186 0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 187 0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ##          hatch        n5        n9
    ## 5   0.09041847 0.1025887 0.1275217
    ## 8   0.09041847 0.1025887 0.1275217
    ## 10  0.11467544 0.1025841 0.1032292
    ## 15  0.11467544 0.1025841 0.1032292
    ## 16  0.09041847 0.1025887 0.1275217
    ## 22  0.09041847 0.1025887 0.1275217
    ## 25  0.09041847 0.1025887 0.1275217
    ## 30  0.11467544 0.1025841 0.1032292
    ## 40  0.11467544 0.1025841 0.1032292
    ## 42  0.09041847 0.1025887 0.1275217
    ## 45  0.09041847 0.1025887 0.1275217
    ## 46  0.09041847 0.1025887 0.1275217
    ## 50  0.09041847 0.1025887 0.1275217
    ## 60  0.11467544 0.1025841 0.1032292
    ## 64  0.11467544 0.1025841 0.1032292
    ## 66  0.09041847 0.1025887 0.1275217
    ## 68  0.09041847 0.1025887 0.1275217
    ## 83  0.09041847 0.1025887 0.1275217
    ## 88  0.11467544 0.1025841 0.1032292
    ## 97  0.11467544 0.1025841 0.1032292
    ## 105 0.11467544 0.1025841 0.1032292
    ## 106 0.11467544 0.1025841 0.1032292
    ## 108 0.11467544 0.1025841 0.1032292
    ## 120 0.09041847 0.1025887 0.1275217
    ## 126 0.09041847 0.1025887 0.1275217
    ## 128 0.11467544 0.1025841 0.1032292
    ## 130 0.09041847 0.1025887 0.1275217
    ## 131 0.09041847 0.1025887 0.1275217
    ## 146 0.11467544 0.1025841 0.1032292
    ## 160 0.11467544 0.1025841 0.1032292
    ## 164 0.09041847 0.1025887 0.1275217
    ## 168 0.11467544 0.1025841 0.1032292
    ## 174 0.11467544 0.1025841 0.1032292
    ## 178 0.09041847 0.1025887 0.1275217
    ## 179 0.11467544 0.1025841 0.1032292
    ## 186 0.09041847 0.1025887 0.1275217
    ## 187 0.11467544 0.1025841 0.1032292
    ## 
    ## $x
    ##            LD1
    ## 5    0.9750932
    ## 8    0.9750932
    ## 10  -0.9750932
    ## 15  -0.9750932
    ## 16   0.9750932
    ## 22   0.9750932
    ## 25   0.9750932
    ## 30  -0.9750932
    ## 40  -0.9750932
    ## 42   0.9750932
    ## 45   0.9750932
    ## 46   0.9750932
    ## 50   0.9750932
    ## 60  -0.9750932
    ## 64  -0.9750932
    ## 66   0.9750932
    ## 68   0.9750932
    ## 83   0.9750932
    ## 88  -0.9750932
    ## 97  -0.9750932
    ## 105 -0.9750932
    ## 106 -0.9750932
    ## 108 -0.9750932
    ## 120  0.9750932
    ## 126  0.9750932
    ## 128 -0.9750932
    ## 130  0.9750932
    ## 131  0.9750932
    ## 146 -0.9750932
    ## 160 -0.9750932
    ## 164  0.9750932
    ## 168 -0.9750932
    ## 174 -0.9750932
    ## 178  0.9750932
    ## 179 -0.9750932
    ## 186  0.9750932
    ## 187 -0.9750932
    ## 
    ## [1] 0.1621622
    ## Call:
    ## lda(treatment ~ sex, data = train.transformed)
    ## 
    ## Prior probabilities of groups:
    ##   control      bldg       lay    inc.d3    inc.d9   inc.d17     hatch 
    ## 0.1282051 0.1025641 0.1025641 0.1025641 0.1282051 0.1153846 0.1025641 
    ##        n5        n9 
    ## 0.1025641 0.1153846 
    ## 
    ## Group means:
    ##           sexmale
    ## control 0.5000000
    ## bldg    0.5625000
    ## lay     0.5625000
    ## inc.d3  0.5000000
    ## inc.d9  0.4500000
    ## inc.d17 0.4444444
    ## hatch   0.4375000
    ## n5      0.5000000
    ## n9      0.5555556
    ## 
    ## Coefficients of linear discriminants:
    ##              LD1
    ## sexmale 1.950186
    ## [1] control control inc.d9  inc.d9  control control
    ## Levels: control bldg lay inc.d3 inc.d9 inc.d17 hatch n5 n9
    ##      control       bldg        lay    inc.d3    inc.d9   inc.d17
    ## 5  0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 8  0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 10 0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 15 0.1282302 0.09041446 0.09041446 0.1025841 0.1403520 0.1275161
    ## 16 0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ## 22 0.1282358 0.11468053 0.11468053 0.1025887 0.1160518 0.1032338
    ##         hatch        n5        n9
    ## 5  0.09041847 0.1025887 0.1275217
    ## 8  0.09041847 0.1025887 0.1275217
    ## 10 0.11467544 0.1025841 0.1032292
    ## 15 0.11467544 0.1025841 0.1032292
    ## 16 0.09041847 0.1025887 0.1275217
    ## 22 0.09041847 0.1025887 0.1275217
    ##           LD1
    ## 5   0.9750932
    ## 8   0.9750932
    ## 10 -0.9750932

![](../figures/sexes/LDA-1.png)

    ## [1] 0.1621622

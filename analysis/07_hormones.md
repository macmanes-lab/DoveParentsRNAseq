    library(tidyverse)

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.2.1     ✔ purrr   0.3.2
    ## ✔ tibble  2.1.3     ✔ dplyr   0.8.1
    ## ✔ tidyr   0.8.3     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

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

    source("../R/functions.R")  # load custom functions 
    source("../R/themes.R")  # load custom themes and color palletes

    knitr::opts_chunk$set(fig.path = '../figures/hormones/',message=F, warning=FALSE)

    charlevels <- c("control", "bldg", "lay", "inc_d3", "inc_d9", "inc_d17", "hatch", "n5", "n9" )
    removelevels <- c("M_inc3",  "M_inc9", "M_inc17", "M_n2")
    maniplevels <- c("early", "prolong", "extend")
    combolevels <- c("control", "bldg", "lay", "inc_d3", "M_inc3",
                      "inc_d9",  "M_inc9", "early",  "inc_d17","M_inc17",  "hatch",
                     "prolong", "extend", "M_n2",
                      "n5", "n9" )

    prolactin <- read_excel("../results/Pigeon prolactin concentrations juil 2018.xlsx", sheet = 1)

    # keep only samples realated to parental care
    prolactin <- prolactin %>% filter(Study %in% c("Baseline", "ParentalCare"))

    # create lists of factos for characterization and manipluation studies

    # rename some of the treatments levels, store in new columns
    prolactin <- prolactin %>%
        mutate(sex = fct_recode(Sex,
                                "female" = "f",
                                "male" = "m"),
               treatment = fct_recode(Treatment,
                                "hatch" = "Hatch",
                                "inc_d17" = "Inc_d17",
                                "inc_d3" = "Inc_d3",
                                "inc_d9" = "Inc_d9",
                                "M_inc9" = "M_Inc9",
                                "M_inc3" = "M_Inc3",
                                "early" = "M_Inc8",
                                "early" = "M_inc8",
                                "M_inc17" = "M_Inc17",
                                "M_n2" = "M_hatch",
                                "control" = "baseline",
                                "n5" = "N5", 
                                "n9" = "N9"),
               study = fct_collapse(treatment,
                                     characterization = charlevels,
                                     manipulation = maniplevels,
                                    removal = removelevels))

    colnames(prolactin)[colnames(prolactin)=="Prolactin ng/mL"] <- "plasma_conc"
    prolactin$hormone <- "prolactin"
    prolactin <- prolactin %>% select(study, treatment, sex, hormone, plasma_conc)   %>% drop_na()

    summary(prolactin)

    ##               study       treatment       sex        hormone         
    ##  characterization:189   inc_d9 : 24   female:160   Length:325        
    ##  manipulation    : 59   control: 23   male  :165   Class :character  
    ##  removal         : 77   inc_d17: 21                Mode  :character  
    ##                         n5     : 21                                  
    ##                         bldg   : 20                                  
    ##                         extend : 20                                  
    ##                         (Other):196                                  
    ##   plasma_conc    
    ##  Min.   :  2.02  
    ##  1st Qu.: 11.07  
    ##  Median : 23.46  
    ##  Mean   : 34.76  
    ##  3rd Qu.: 56.98  
    ##  Max.   :203.69  
    ## 

    PETC <- read_excel("../results/parental_care_hormone_RIA_data_master.xlsx", sheet = 2)

    PETC <- PETC %>% select(stage, sex, hormone, plasma_conc)  %>%
                    filter(stage %in% combolevels)  %>%
                    droplevels()  
    PETC$stage <- factor(PETC$stage)
    levels(PETC$stage)

    ##  [1] "bldg"    "extend"  "hatch"   "inc_d17" "inc_d3"  "inc_d9"  "lay"    
    ##  [8] "n5"      "n9"      "prolong"

    PETC <- PETC %>%
        mutate(sex = fct_recode(sex,
                                "female" = "f",
                                "male" = "m"),
               study = fct_collapse(stage,
                                     characterization = charlevels,
                                    manipulation = maniplevels,
                                    removal = removelevels)) %>% 
      drop_na()
    colnames(PETC)[colnames(PETC)=="stage"] <- "treatment"
    PETC <- PETC %>% select(study, treatment, sex, hormone, plasma_conc)
    summary(PETC)

    ##               study       treatment       sex        hormone         
    ##  characterization:490   inc_d9 : 88   female:332   Length:603        
    ##  manipulation    :113   inc_d17: 68   male  :271   Class :character  
    ##                         bldg   : 60                Mode  :character  
    ##                         hatch  : 60                                  
    ##                         extend : 59                                  
    ##                         n9     : 58                                  
    ##                         (Other):210                                  
    ##   plasma_conc      
    ##  Min.   :  0.0355  
    ##  1st Qu.:  0.2402  
    ##  Median :  0.8358  
    ##  Mean   :  1.4721  
    ##  3rd Qu.:  1.6961  
    ##  Max.   :117.7708  
    ## 

    hormones <- rbind(prolactin, PETC)
    hormones$treatment <- factor(hormones$treatment, levels = combolevels)

    hormones$okay <- ifelse(hormones$hormone == "cort" & hormones$plasma_conc > 30, "bad",
                        ifelse(hormones$hormone == "progesterone" & hormones$plasma_conc > 5, "bad", 
                               ifelse(hormones$hormone == "prolactin" & hormones$plasma_conc > 150, "bad", 
                                ifelse(hormones$hormone == "prolactin" & hormones$treatment == "control", "bad", 
                            ifelse(hormones$hormone == "testosterone" & hormones$sex == "female", "bad",
                                   ifelse(hormones$hormone == "estradiol" & hormones$sex == "male", "bad", "okay"))))))
    hormones <- hormones %>% filter(okay == "okay") %>% droplevels()


    prl.char <- hormones %>% filter(hormone == "prolactin", study == "characterization")   %>%  droplevels()
    test.char <- hormones %>% filter(hormone == "testosterone", treatment %in% charlevels)   %>%  droplevels()
    est.char <- hormones %>% filter(hormone == "estradiol", treatment %in% charlevels)   %>%  droplevels()
    prog.char <- hormones %>% filter(hormone == "progesterone", treatment %in% charlevels)   %>%  droplevels()
    cort.char <- hormones %>% filter(hormone == "cort", treatment %in% charlevels)   %>%  droplevels()

    hormonecharplot <- function(myhormone, myylab){
      
      mycolors <- c("female" = "#F8766D", "male" = "#00BFC4")
      
      hormones %>% 
        filter(study == "characterization",
               hormone %in% c(myhormone))  %>% 
      ggplot(aes(x = treatment, y = plasma_conc, fill = sex)) +
        geom_boxplot() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none") +
        scale_fill_manual(values = mycolors) +
        labs(y = myylab, x = NULL) 
    }

    a <- hormonecharplot("estradiol", "E (ng/mL)")
    b <- hormonecharplot("testosterone", "T (ng/mL)")
    c <- hormonecharplot("cort", "CORT (ng/mL)")
    d <- hormonecharplot("progesterone", "PROG (ng/mL)")

    e <- hormonecharplot("prolactin", "PRL (ng/mL)")

    # model and plot prolactin data
    mod_prolactin <-  lm(data = prl.char,  plasma_conc ~ treatment + sex)
    grid <- prl.char %>%
      data_grid(treatment, sex) %>%
      add_predictions(mod_prolactin) 

    f <- ggplot(grid, aes(x = treatment, y = pred, color = sex)) +
      geom_point() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            legend.title = element_blank()) +
      labs(y = "predicted PRL (ng/mL)", x = NULL) +
      geom_smooth(se = F, aes(x = as.numeric(treatment))) 

    plot_grid(e,f,  rel_widths = c(0.65,0.35))

![](../figures/hormones/characterization-hormones-1.png)

    plot_grid(d + theme(axis.text.x=element_blank()),
              a + theme(axis.text.x=element_blank()),
              c, b, 
              rel_widths = c(0.65,0.35),
              rel_heights = c(0.45,0.55), ncol = 2)

![](../figures/hormones/characterization-hormones-2.png)

    aov1 <- aov(data = prl.char, plasma_conc ~ treatment + sex)
    aov1

    ## Call:
    ##    aov(formula = plasma_conc ~ treatment + sex, data = prl.char)
    ## 
    ## Terms:
    ##                 treatment      sex Residuals
    ## Sum of Squares   85878.31  1420.71  47900.04
    ## Deg. of Freedom         7        1       156
    ## 
    ## Residual standard error: 17.52289
    ## Estimated effects may be unbalanced

    summary(aov1)

    ##              Df Sum Sq Mean Sq F value Pr(>F)    
    ## treatment     7  85878   12268  39.955 <2e-16 ***
    ## sex           1   1421    1421   4.627  0.033 *  
    ## Residuals   156  47900     307                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov1, which = "treatment")

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = plasma_conc ~ treatment + sex, data = prl.char)
    ## 
    ## $treatment
    ##                       diff        lwr        upr     p adj
    ## lay-bldg         2.8425808 -14.182911  19.868073 0.9995830
    ## inc_d3-bldg     11.2175616  -5.807930  28.243053 0.4689723
    ## inc_d9-bldg     11.4894835  -4.811184  27.790151 0.3783317
    ## inc_d17-bldg    51.0003252  34.178739  67.821911 0.0000000
    ## hatch-bldg      65.3352160  48.309724  82.360708 0.0000000
    ## n5-bldg         45.7474144  28.721923  62.772906 0.0000000
    ## n9-bldg         31.2185488  14.193057  48.244041 0.0000022
    ## inc_d3-lay       8.3749808  -8.650511  25.400473 0.8003128
    ## inc_d9-lay       8.6469027  -7.653765  24.947570 0.7315620
    ## inc_d17-lay     48.1577444  31.336158  64.979331 0.0000000
    ## hatch-lay       62.4926352  45.467143  79.518127 0.0000000
    ## n5-lay          42.9048336  25.879342  59.930325 0.0000000
    ## n9-lay          28.3759680  11.350476  45.401460 0.0000241
    ## inc_d9-inc_d3    0.2719219 -16.028745  16.572589 1.0000000
    ## inc_d17-inc_d3  39.7827636  22.961177  56.604350 0.0000000
    ## hatch-inc_d3    54.1176544  37.092163  71.143146 0.0000000
    ## n5-inc_d3       34.5298528  17.504361  51.555345 0.0000001
    ## n9-inc_d3       20.0009872   2.975495  37.026479 0.0095965
    ## inc_d17-inc_d9  39.5108417  23.423264  55.598419 0.0000000
    ## hatch-inc_d9    53.8457325  37.545065  70.146400 0.0000000
    ## n5-inc_d9       34.2579309  17.957264  50.558598 0.0000000
    ## n9-inc_d9       19.7290653   3.428398  36.029733 0.0066240
    ## hatch-inc_d17   14.3348908  -2.486695  31.156477 0.1572868
    ## n5-inc_d17      -5.2529108 -22.074497  11.568675 0.9792976
    ## n9-inc_d17     -19.7817764 -36.603363  -2.960190 0.0094781
    ## n5-hatch       -19.5878016 -36.613293  -2.562310 0.0122801
    ## n9-hatch       -34.1166672 -51.142159 -17.091175 0.0000002
    ## n9-n5          -14.5288656 -31.554357   2.496626 0.1560036

    aov2 <- aov(data = prog.char, plasma_conc ~ treatment + sex)
    aov2

    ## Call:
    ##    aov(formula = plasma_conc ~ treatment + sex, data = prog.char)
    ## 
    ## Terms:
    ##                 treatment      sex Residuals
    ## Sum of Squares    2.62134  0.08244  67.87352
    ## Deg. of Freedom         7        1       181
    ## 
    ## Residual standard error: 0.6123658
    ## Estimated effects may be unbalanced

    summary(aov2)

    ##              Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment     7   2.62  0.3745   0.999  0.434
    ## sex           1   0.08  0.0824   0.220  0.640
    ## Residuals   181  67.87  0.3750

    aov3 <- aov(data = cort.char, plasma_conc ~ treatment + sex)
    aov3

    ## Call:
    ##    aov(formula = plasma_conc ~ treatment + sex, data = cort.char)
    ## 
    ## Terms:
    ##                 treatment      sex Residuals
    ## Sum of Squares    34.7566   0.1393  368.9301
    ## Deg. of Freedom         7        1       168
    ## 
    ## Residual standard error: 1.481895
    ## Estimated effects may be unbalanced

    summary(aov3)

    ##              Df Sum Sq Mean Sq F value Pr(>F)  
    ## treatment     7   34.8   4.965   2.261 0.0318 *
    ## sex           1    0.1   0.139   0.063 0.8015  
    ## Residuals   168  368.9   2.196                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov3, which = "treatment")

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = plasma_conc ~ treatment + sex, data = cort.char)
    ## 
    ## $treatment
    ##                         diff        lwr        upr     p adj
    ## lay-bldg        0.5878084905 -0.8175481 1.99316509 0.9037028
    ## inc_d3-bldg     0.8797827790 -0.5980619 2.35762741 0.6024523
    ## inc_d9-bldg     0.9196932651 -0.3488813 2.18826782 0.3420874
    ## inc_d17-bldg    0.9714091341 -0.4193232 2.36214144 0.3913452
    ## hatch-bldg      0.9204152143 -0.5368169 2.37764734 0.5264276
    ## n5-bldg         1.5748863713  0.1536853 2.99608740 0.0185863
    ## n9-bldg         0.1833293407 -1.2945153 1.66117398 0.9999436
    ## inc_d3-lay      0.2919742885 -1.1537006 1.73764917 0.9985698
    ## inc_d9-lay      0.3318847747 -0.8990630 1.56283258 0.9913179
    ## inc_d17-lay     0.3836006436 -0.9728976 1.74009887 0.9884590
    ## hatch-lay       0.3326067238 -1.0919903 1.75720370 0.9964109
    ## n5-lay          0.9870778808 -0.4006409 2.37479668 0.3671432
    ## n9-lay         -0.4044791498 -1.8501540 1.04119573 0.9891669
    ## inc_d9-inc_d3   0.0399104862 -1.2731889 1.35300990 1.0000000
    ## inc_d17-inc_d3  0.0916263551 -1.3398362 1.52308890 0.9999994
    ## hatch-inc_d3    0.0406324353 -1.4555207 1.53678555 1.0000000
    ## n5-inc_d3       0.6951035923 -0.7659786 2.15618574 0.8271686
    ## n9-inc_d3      -0.6964534382 -2.2126900 0.81978316 0.8516210
    ## inc_d17-inc_d9  0.0517158689 -1.1625088 1.26594058 1.0000000
    ## hatch-inc_d9    0.0007219491 -1.2891349 1.29057883 1.0000000
    ## n5-inc_d9       0.6551931061 -0.5938136 1.90419979 0.7436121
    ## n9-inc_d9      -0.7363639244 -2.0494633 0.57673549 0.6733815
    ## hatch-inc_d17  -0.0509939198 -1.4611661 1.35917831 1.0000000
    ## n5-inc_d17      0.6034772372 -0.7694294 1.97638387 0.8783490
    ## n9-inc_d17     -0.7880797933 -2.2195423 0.64338275 0.6937894
    ## n5-hatch        0.6544711570 -0.7857586 2.09470091 0.8585405
    ## n9-hatch       -0.7370858735 -2.2332390 0.75906724 0.7999615
    ## n9-n5          -1.3915570305 -2.8526392 0.06952512 0.0743432

    aov4 <- aov(data = est.char, plasma_conc ~ treatment )
    aov4

    ## Call:
    ##    aov(formula = plasma_conc ~ treatment, data = est.char)
    ## 
    ## Terms:
    ##                 treatment Residuals
    ## Sum of Squares   0.812090  4.213773
    ## Deg. of Freedom         7        68
    ## 
    ## Residual standard error: 0.2489322
    ## Estimated effects may be unbalanced

    summary(aov4)

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## treatment    7  0.812 0.11601   1.872 0.0877 .
    ## Residuals   68  4.214 0.06197                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

    TukeyHSD(aov4, which = "treatment")

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = plasma_conc ~ treatment, data = est.char)
    ## 
    ## $treatment
    ##                        diff         lwr       upr     p adj
    ## lay-bldg       -0.045252581 -0.38539108 0.2948859 0.9998904
    ## inc_d3-bldg    -0.008174226 -0.38456016 0.3682117 1.0000000
    ## inc_d9-bldg    -0.197921854 -0.52986304 0.1340193 0.5793020
    ## inc_d17-bldg   -0.095196466 -0.43533497 0.2449420 0.9873323
    ## hatch-bldg     -0.025001550 -0.40138748 0.3513844 0.9999991
    ## n5-bldg        -0.116222267 -0.45636077 0.2239162 0.9611413
    ## n9-bldg         0.165831705 -0.17430680 0.5059702 0.7915163
    ## inc_d3-lay      0.037078355 -0.34655639 0.4207131 0.9999875
    ## inc_d9-lay     -0.152669273 -0.49280777 0.1874692 0.8524877
    ## inc_d17-lay    -0.049943885 -0.39808674 0.2981990 0.9998183
    ## hatch-lay       0.020251031 -0.36338371 0.4038858 0.9999998
    ## n5-lay         -0.070969686 -0.41911254 0.2771732 0.9981889
    ## n9-lay          0.211084286 -0.13705857 0.5592271 0.5583251
    ## inc_d9-inc_d3  -0.189747628 -0.56613356 0.1866383 0.7623475
    ## inc_d17-inc_d3 -0.087022240 -0.47065698 0.2966125 0.9964490
    ## hatch-inc_d3   -0.016827324 -0.43293763 0.3992830 1.0000000
    ## n5-inc_d3      -0.108048041 -0.49168278 0.2755867 0.9868573
    ## n9-inc_d3       0.174005930 -0.20962881 0.5576407 0.8456181
    ## inc_d17-inc_d9  0.102725388 -0.23741311 0.4428639 0.9803468
    ## hatch-inc_d9    0.172920304 -0.20346563 0.5493062 0.8368839
    ## n5-inc_d9       0.081699587 -0.25843892 0.4218381 0.9949411
    ## n9-inc_d9       0.363753558  0.02361506 0.7038921 0.0277306
    ## hatch-inc_d17   0.070194916 -0.31343983 0.4538297 0.9990962
    ## n5-inc_d17     -0.021025801 -0.36916866 0.3271171 0.9999995
    ## n9-inc_d17      0.261028171 -0.08711469 0.6091710 0.2854692
    ## n5-hatch       -0.091220717 -0.47485546 0.2924140 0.9952458
    ## n9-hatch        0.190833254 -0.19280149 0.5744680 0.7743180
    ## n9-n5           0.282053971 -0.06608888 0.6301968 0.1993673

    aov5 <- aov(data = test.char, plasma_conc ~ treatment )
    aov5

    ## Call:
    ##    aov(formula = plasma_conc ~ treatment, data = test.char)
    ## 
    ## Terms:
    ##                 treatment Residuals
    ## Sum of Squares    7.58604  45.89673
    ## Deg. of Freedom         7        33
    ## 
    ## Residual standard error: 1.179326
    ## Estimated effects may be unbalanced

    summary(aov5)

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## treatment    7   7.59   1.084   0.779  0.609
    ## Residuals   33  45.90   1.391

    write.csv(hormones, "../results/hormones.csv", row.names = F)

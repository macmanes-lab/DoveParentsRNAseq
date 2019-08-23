All characterization data
-------------------------

This first image was made using
[`08_WGCNA_1.R`](https://github.com/macmanes-lab/DoveParentsRNAseq/blob/master/analysis/08_WGCNA_1.R),
which was modified from the [first in a series of WGCNA
tutorial](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-DataInput.R).
The first tutorial focuses on data input, cleaning and pre-processing;
outlier removal; clustering of samples by their Euclidean distance; a
comparison of sample cluster to sample meta data.

    ##                  treatment sex tissue
    ## M.gon.control.45   control   M    gon
    ## M.hyp.control.45   control   M    hyp
    ## M.pit.control.45   control   M    pit
    ## M.gon.control.46   control   M    gon
    ## M.hyp.control.46   control   M    hyp
    ## M.pit.control.46   control   M    pit
    ##   treatment sex tissue
    ## 1         2   2      1
    ## 2         2   2      2
    ## 3         2   2      3
    ## 4         2   2      1
    ## 5         2   2      2
    ## 6         2   2      3

![](../figures/wgcna/wgcna-1.png)

Females
-------

    ##                  treatment sex tissue
    ## F.gon.control.47   control   F    gon
    ## F.hyp.control.47   control   F    hyp
    ## F.pit.control.47   control   F    pit
    ## F.gon.control.71   control   F    gon
    ## F.hyp.control.71   control   F    hyp
    ## F.pit.control.71   control   F    pit
    ##   treatment sex tissue
    ## 1         2   1      1
    ## 2         2   1      2
    ## 3         2   1      3
    ## 4         2   1      1
    ## 5         2   1      2
    ## 6         2   1      3

![](../figures/wgcna/wgcna-female-1.png)

Males
-----

    ##                  treatment sex tissue
    ## M.gon.control.45   control   M    gon
    ## M.hyp.control.45   control   M    hyp
    ## M.pit.control.45   control   M    pit
    ## M.gon.control.46   control   M    gon
    ## M.hyp.control.46   control   M    hyp
    ## M.pit.control.46   control   M    pit
    ##   treatment sex tissue
    ## 1         2   1      1
    ## 2         2   1      2
    ## 3         2   1      3
    ## 4         2   1      1
    ## 5         2   1      2
    ## 6         2   1      3

![](../figures/wgcna/wgcna-male-1.png)

hypothalamus
------------

    ##                  treatment sex tissue
    ## M.hyp.control.45   control   M    hyp
    ## M.hyp.control.46   control   M    hyp
    ## F.hyp.control.47   control   F    hyp
    ## M.hyp.control.48   control   M    hyp
    ## M.hyp.control.49   control   M    hyp
    ## M.hyp.control.52   control   M    hyp
    ##   treatment sex tissue
    ## 1         2   2      1
    ## 2         2   2      1
    ## 3         2   1      1
    ## 4         2   2      1
    ## 5         2   2      1
    ## 6         2   2      1

![](../figures/wgcna/wgcna-hyp-1.png)

pituitary
---------

gonads
------

    ##                  treatment sex tissue
    ## M.gon.control.45   control   M    gon
    ## M.gon.control.46   control   M    gon
    ## F.gon.control.47   control   F    gon
    ## M.gon.control.48   control   M    gon
    ## M.gon.control.49   control   M    gon
    ## M.gon.control.52   control   M    gon
    ##   treatment sex tissue
    ## 1         2   2      1
    ## 2         2   2      1
    ## 3         2   1      1
    ## 4         2   2      1
    ## 5         2   2      1
    ## 6         2   2      1

![](../figures/wgcna/wgcna-gon-1.png)

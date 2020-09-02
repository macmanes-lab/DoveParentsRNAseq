[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/macmanes-lab/DoveParentsRNAseq/master?urlpath=rstudio) _Note: Binder not working :(_

[![Shiny](https://img.shields.io/badge/Explore%20the%20data-Open%20and%20R%20Shiny%20app-blue)](https://raynamharris.shinyapps.io/musicalgenes/) To quickly explore the data in an internet browser, use [our R Shiny app](https://raynamharris.shinyapps.io/musicalgenes/). 

#  The reproductive and parental care transcriptome of the rock dove 

The aims of the research associated with this GitHub repository reree-fold: 

1. produce a high quality and reproducible transcriptional characterization of the HPG over reproduction and parental care in male and female rock doves. We were specifically interested in understanding how gene activity in the HPG changed as individuals transitioned from a sexually mature, non-reproductive state into reproduction and parenting.
1. characterize the major patterns of transcriptional changes between the tissues, sexes, and parental stages.
1. identify the drivers underlying changes in transcription. Specifically, we tested whether external or internal factors were regulating gene activity.


## Methods

These data anlyseses are conducted in R and the workflow is automated via Snakemake. The order of operations are descriped in the [Snakefile](https://github.com/macmanes-lab/DoveParentsRNAseq/blob/master/Snakefile).

The data, scripts, and resutls are oranzied in the following manner:

- [R](https://github.com/macmanes-lab/DoveParentsRNAseq/tree/master/R): contains custom functions and themes that are used by multiple anlaysis files
- [analysis](https://github.com/macmanes-lab/DoveParentsRNAseq/tree/master/analysis): contains .R scripts that are executed by snakemake live
- [figures](https://github.com/macmanes-lab/DoveParentsRNAseq/tree/master/figures): contains images and illustrations
- [metadata](https://github.com/macmanes-lab/DoveParentsRNAseq/tree/master/metadata): contains files with biological descriptions of samples or genes
- [results](https://github.com/macmanes-lab/DoveParentsRNAseq/tree/master/results): contains outputs 

### Figures

**_[High quality PDFs of all figures are available here](https://github.com/macmanes-lab/DoveParentsRNAseq/tree/master/figures)_**

![](./figures/fig1-1.png)
**Figure 1 Experimental design.** **A)** We investigated nine timepoints that spanned the majority of reproductive efforts in this species. These time-points consisted of a "control" or non-parenting state (from MacManes et al 2016 and Calisi et al. 2017), nest-building, clutch initiation or the onset of incubation, clutch completion and early incubation, mid-incubation, late incubation, initiation of nestling care, early nestling care, and mid-nestling care. To teste whether external or internal factors were regulating gene activity, we also conducted a series of offspring removal **(B)** and replacement **(C)** manipulations to test whether transcriptional changes were dependent upon offspring presence or were regulated by an internal clock. **(D)** We used t-Distributed Stochastic Neighbor Embedding (t-SNE) to reduce the dimensionality of the transcriptomes from the hypothalamus, pituitary, and gonads of male and female rock. Circles and triangles represent female and male samples, respectively, and points are colored by treatment. 


### Figure 2 

![](./figures/fig2-1.png)
**Figure 2. Hypothesis-driven and data-driven approaches to identifying changes in gene expression.** Box-and-whisker plots show changes in gene expression of serotonin receptor (*HTR2C*) in the hypothalamus **(A)**, prolactin (*PRL*) in the pituitary **(B)**,  and estrogen receptor 1 (*ESR1*) in the gonads.Bar and statistics are shown for one contrast of interest between treatment groups. Boxes are colored by treatment. Volcano plots show the log-fold change (LFC) and -log10(FDR) for all genes that are significantly different differentially between the aforementioned treatment groups. Boxes are colored by the treatment were gene expression was higher.    

### Figure 3 

![](./figures/fig3-1.png)

**Figure 3. Summary of differentially expressed genes.** Bar plots show the total number of significantly differentially expressed genes for all two-way comparisons. Contrasts are made relative to the non-breeding controls, the nest-building group, the previous stage, or the internal or external control group. Bars are colored by the treatment where gene expression was higher. A negative number of differentially expressed genes means that many genes had increase expression in the reference group.

## Tables

- [Table 1](https://github.com/macmanes-lab/DoveParentsRNAseq/blob/master/results/table1.csv)
- [Table 2](https://github.com/macmanes-lab/DoveParentsRNAseq/blob/master/results/table2.csv)
- [Suppl Table 1](https://github.com/macmanes-lab/DoveParentsRNAseq/blob/master/results/suppltable1.csv)
- [Suppl Table 2](https://github.com/macmanes-lab/DoveParentsRNAseq/blob/master/results/suppltable2.csv)
- [All differentially expressed genes](https://github.com/macmanes-lab/DoveParentsRNAseq/blob/master/results/03_allDEG.csv)


## Related documentation 

- Dove Genomics Project website <http://www.dovelovegenomics.org/>
- Shiny app for Data exploration <https://raynamharris.shinyapps.io/musicalgenes/>
- A talk from the Society for Integrative and Comparative Biology 2020 
<https://speakerdeck.com/raynamharris/peaks-and-valleys-of-prolactin-driven-gene-expression-during-parental-care>.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/macmanes-lab/DoveParentsRNAseq/master?urlpath=rstudio)

# Characterizing the neurogenomics of parental care in the rock dove

## Overview

This repository contains the data and analysis for a collaboration between Drs. Rebecca Calisi and Matt MacManes that focuses on one characterizing the neurogenomocs of parental care in the rock dove.

## Main research questions

1. What genes differ in expression between timepoints?
1. What stages are most different, stressful, responsive? 
1. How do genes in the HPG affect genes in other regions?
1. How do male and females differ?
1. How do tissues differ?
1. What is the relationship between genes and hormones across timepoints?
 

## Organization

These repositories is broken down into the following sub-repositories, each with their own unique purpose and structure.

- **analysis**: where the *.Rmd* script and the *.md* outputs live. The prefix corresponds to the order of operation. 
- **figures**: where figure generated by the scripts live. The prefix correspond to the script prefix that created the files
- **metadata**: contains files that describe the sample variables and a rosseta stone for transcripts to gene ids. 
- **results**: where data generated by the scripts live. The prefix correspond to the script prefix that created the files

There are two hidden directories, `kallisto_mappings` and `mapping`, which contain the results of the kallisto and salmon algorithms that transform read counts into gene counts.  

## More information

For more details about the collaboration, visit <http://www.dovelovegenomics.org/>. See also the GitHub repository <https://github.com/macmanes-lab/RockDove> which contains related data and scripts. 


rule all:
  input:
    "results/01_limma.csv",
    "results/DESeq2/treatment/male_hypothalamus_vsd.csv",
    "results/DESeq2/treatment/female_hypothalamus_vsd.csv"

rule wrangledata:
  input:
    "metadata/kallistosamples.txt",
    "results/kallistocounts.txt"
  output:
    "metadata/00_colData.csv",
    "metadata/00_geneinfo.csv",
    "metadata/00_birds.csv",
    "results/00_counts.csv",
    "results/00_geneswithisoforms.csv"
  shell:
    "Rscript analysis/00_datawrangling.R"
    
rule limma:
  input:
    "metadata/00_colData.csv",
    "results/00_counts.csv",
  output:
    "results/01_limma.csv"
  shell:
    "Rscript analysis/01_limma.R"
    
rule deseq2first:
  input:
    "metadata/00_colData.csv",
    "results/00_counts.csv",
  output:
    "results/DESeq2/treatment/male_hypothalamus_vsd.csv",
    "results/DESeq2/treatment/female_hypothalamus_vsd.csv",
    "results/DESeq2/treatment/male_hypothalamus_control_bldg_DEGs.csv"
  shell:
    "Rscript analysis/02_DESeq2.R"
rule wrangledata:
  input:
    "metadata/kallistosamples.txt",
    "results/kallistocounts.txt",
  output:
    "metadata/00_colData.csv",
    "metadata/00_geneinfo.csv",
    "metadata/00_birds.csv",
    "results/00_counts.csv",
    "results/00_geneswithisoforms.csv",
  shell:
    "Rscript analysis/00_datawrangling.R"
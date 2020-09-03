rule wrangle:
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
    "Rscript analysis/00_wrangle.R"
    
rule fig1:
  input:
    "metadata/00_colData.csv",
    "results/00_counts.csv"
  output:
    "figures/fig1-1.pdf"
  shell:
    "Rscript analysis/01_fig1.R"  
    
rule deseq2:
  input:
    "metadata/00_colData.csv",
    "results/00_counts.csv"
  output:
    "results/DESeq2/treatment/female_pituitary_control_bldg_DEGs.csv"
  threads: 6
  shell:
    "Rscript analysis/02_DESeq2.R"
    
rule vsd:
  input:
    "results/DESeq2/treatment/female_hypothalamus_vsd.csv",
    "R/genelists.R"
  output:
    "results/03_hypvsdf.csv",
    "results/03_candidatevsd.csv"
  shell:
    "Rscript analysis/03_vsd.R"    
    
rule degs:
  input:
    "results/DESeq2/treatment/female_gonads_control_early_DEGs.csv",
    "results/DESeq2/treatment/male_gonads_inc.d3_m.inc.d3_DEGs.csv"
  output:
    "results/03_allDEG.csv"
  shell:
    "Rscript analysis/04_DEGs.R"
    

rule figs23:
  input:
    "results/03_allDEG.csv"
  output:
    "figures/fig2-1.pdf",
    "figures/fig3-1.pdf"
  shell:
    "Rscript analysis/05_figs23.R" 
    
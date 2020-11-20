rule wrangle:
  input:
    "metadata/kallistosamples.txt",
    "results/kallistocounts.txt"
  output:
    "metadata/00_colData.csv",
    "results/00_countData.csv"
  shell:
    "Rscript analysis/00_wrangle.R"
    
rule fig1:
  input:
    "metadata/00_colData.csv",
    "results/00_countData.csv"
  output:
    "figures/fig1-1.pdf"
  shell:
    "Rscript analysis/01_fig1.R"  
    
rule deseq2:
  input:
    "metadata/00_colData.csv",
    "results/00_countData.csv"
  output:
    "results/DESeq2/treatment/female_hypothalamus_vsd.csv",
    "results/DESeq2/treatment/male_hypothalamus_vsd.csv"
  threads: 6
  shell:
    "Rscript analysis/02_DESeq2.R"
    
rule vsd:
  input:
    "results/DESeq2/treatment/female_hypothalamus_vsd.csv",
    "R/genelists.R"
  output:
    "results/03_hypvsdf.csv",
    "results/03_candidatevsd.csv",
    "results/03_shinyvsd.csv"
  shell:
    "Rscript analysis/03_vsd.R"    
    
rule degs:
  input:
    "results/DESeq2/treatment/female_gonads_control_early_DEGs.csv"
  output:
    "results/04_allDEG.csv",
    "results/table1a.csv",
    "results/tableS1.csv"
  shell:
    "Rscript analysis/04_DEGs.R"
    

rule figs23:
  input:
    "results/03_hypvsdf.csv",
    "results/04_allDEG.csv"
  output:
    "figures/fig2-1.pdf",
    "figures/fig3-1.pdf"
  shell:
    "Rscript analysis/05_figs23.R" 
    
rule figsup1:
  input:
    "figures/images/fig_supfig1a.png"
  output:
    "figures/figsup1-1.pdf"
  shell:
    "Rscript analysis/06_figsup1.R"     
library(tidyverse)
library(VennDiagram)

d1 <- read_csv("results/DEseq2/treatment/female_hypothalamus_hatch_n5_DEGs.csv") %>% 
  filter(direction != "NS")
d2 <- read_csv("results/DEseq2/treatment/female_hypothalamus_n5_n9_DEGs.csv") %>%
  filter(direction != "NS")
d3 <- read_csv("results/DEseq2/treatment/female_hypothalamus_hatch_prolong_DEGs.csv") %>% 
  filter(direction != "NS")
d4 <- read_csv("results/DEseq2/treatment/female_hypothalamus_hatch_extend_DEGs.csv") %>% 
  filter(direction != "NS")

d1genes <- d1$gene
d2genes <- d2$gene
d3genes <- d3$gene
d4genes <- d4$gene

venn.diagram(
  x = list(d1genes, d2genes, d3genes, d4genes),
  category.names = c("female hyp\nhatch v n5" , "female hyp\n n5 v n9" , 
                     "female hyp\n hatch v prolong" , "female hyp\n hatch v extend" ),
  filename = '~/Desktop/venn_diagram1.png',
  output=TRUE, 
  main = "Number of genes that change expression"
  )


d5 <- read_csv("results/DEseq2/treatment/male_hypothalamus_hatch_n5_DEGs.csv") %>% 
  filter(direction != "NS")
d6 <- read_csv("results/DEseq2/treatment/male_hypothalamus_n5_n9_DEGs.csv") %>%
  filter(direction != "NS")
d7 <- read_csv("results/DEseq2/treatment/male_hypothalamus_hatch_prolong_DEGs.csv") %>% 
  filter(direction != "NS")
d8 <- read_csv("results/DEseq2/treatment/male_hypothalamus_hatch_extend_DEGs.csv") %>% 
  filter(direction != "NS")

d5genes <- d5$gene
d6genes <- d6$gene
d7genes <- d7$gene
d8genes <- d8$gene

venn.diagram(
  x = list(d5genes, d6genes, d7genes, d8genes),
  category.names = c("male hyp\nhatch v n5" , "male hyp\n n5 v n9" , 
                     "male hyp\n hatch v prolong" , "male hyp\n hatch v extend" ),
  filename = '~/Desktop/venn_diagram2.png',
  output=TRUE, 
  main = "Number of genes that change expression"
)



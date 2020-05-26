# genelists

## candidate genes from GO and literature
LitGOdf <- read_csv("../metadata/03_parentalcaregenes.csv") 
curleychampagnegenes <- LitGOdf %>% filter(literature != "NA" ) %>% pull(gene)
GOgenes <- LitGOdf %>% filter(GO != "NA") %>% pull(gene)
parentalcaregenes <- LitGOdf %>% pull(gene)

## genes WGCNA prl module
WGCNAgenes <- read_csv("../results/05_PRLmodule.csv") %>% pull(x)

## candidate genes from breast and ovarian cancers 
## https://www.gynecologiconcology-online.net/article/S0090-8258(19)30069-1/fulltext
suszynskaagenes <- c("BRCA1","BRCA2", "CDKN2A", "PTEN", "PALB2", "TP53", "CDH1", "ATM",
                     "BARD1", "MESH6","MSH2", "BRIP1", "NBN", "FANCC", "FANCM", 
                     "RAD51C","RAD51D")

shaidgenes <- c("GNAS", "USB8", "PIK3CA", "GPR101","RAS","MEN1", "AIP", "DICER1", 
                "PRKAR1A", "PRKACA","SDH", "GPR101")

cancergenes <- c(suszynskaagenes, shaidgenes)

candidategenes <- c(parentalcaregenes, WGCNAgenes, cancergenes)

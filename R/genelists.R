library(tidyverse)

# genelists

## all genes
geneids <- read_csv("metadata/00_geneinfo.csv") %>% select(-X1) 

favoritegenes <- c("AR", "AVP", "AVPR1B", "AVPR2",
                   "CGA", "CRH", "CRHR1",
                   "DIO1", "DIO2", "DIO3", "ESR1",  "ESR2", "FSHB",
                   "FSHR", "GABRQ", "GAL", "GALR1",
                   "GNRH1", "GNRHR",
                   "JAK2", 
                   "LHCGR", 
                   "MC3R", "MC4R", 
                   "NPFFR1", "NPVF", "NPY",
                   "NR3C2", "OXT", "POMC", 
                   "PRL", "PRLH",	"PRLHR", 
                   "SERPINA4", "STAT5B", "VIP", "VIPR1")

## parental care genes 
## genes from Ch 17 Evolution of parental care
curleychampagnegenes <- c("AVP", "AVPR1A", "ADRA2A",    
                          "COMT", "CRH", "CRHBP", "CRHR1", "CRHR2", 
                          "DRD1", "DRD4","ESR1", "ESR2", #ERa ERbeta
                          "FOS", "HTR2C", "NR3C1", #GR
                          "OPRM1", "OXT" , "PGR", "PRL", "PRLR",  "SLC6A4") #5HTT

datadrivengenes <- c("COX1", "COX2", "CYTB", "ATP2B4", "LOC107055658", 
                     "CHGA", "TF", "OVALY", "ANXA5")

literaturegenes <- c(curleychampagnegenes, datadrivengenes,
                     "VIP", "GNIH", "GAL", "TH", "THRB",
                     "GNRHR", "GNRH1", "CGNRH-R", "FSHB", "FSHR")

literaturegenes <- as.data.frame(literaturegenes) %>%
  mutate(GO = "literature", gene = literaturegenes)  %>%
  select(GO, gene)

## candidate gens from parental care GO terms 
GO_path <- "metadata/goterms/"   # path to the data
GO_files <- dir(GO_path, pattern = "*.txt") # get file names
GO_pathfiles <- paste0(GO_path, GO_files)

GOgenes <- GO_pathfiles %>%
  setNames(nm = .) %>% 
  map_df(~read_table(.x, col_types = cols(), col_names = FALSE), .id = "file_name") %>% 
  mutate(GO = sapply(strsplit(as.character(file_name),'metadata/goterms/'), "[", 2)) %>% 
  mutate(GO = sapply(strsplit(as.character(GO),'.txt'), "[", 1)) %>% 
  mutate(gene = sapply(strsplit(as.character(X1), "[\\\\]|[^[:print:]]" ), "[", 2)) %>% 
  select(GO, gene)  %>%
  filter(gene != "Symbol") %>%
  distinct(GO,gene)  %>%
  mutate(gene = toupper(gene))

parentalcaregenes <- GOgenes %>%
  rbind(., literaturegenes) %>%
  arrange(gene) %>%
  left_join(., geneids, by = "gene") %>%
  arrange(gene) %>%
  drop_na()  %>% 
  pivot_wider(
    names_from = GO,
    values_from = gene,
    values_fill = list(gene = "")) %>%
  left_join(geneids, by = c("geneid","NCBI")) %>%
  dplyr::rename("GO"= "parentalbehavior" )
head(parentalcaregenes)

listoparentalcaregenes <- parentalcaregenes %>% pull(gene)

## breast and ovarian cancers  https://www.gynecologiconcology-online.net/article/S0090-8258(19)30069-1/fulltext
breastcancergenes <- c("BRCA1","BRCA2", "CDKN2A", "PTEN", "PALB2", "TP53", "CDH1", "ATM",
                     "BARD1", "MESH6","MSH2", "BRIP1", "NBN", "FANCC", "FANCM", 
                     "RAD51C","RAD51D")
ovariancancergenes <- c("GNAS", "USB8", "PIK3CA", "GPR101","RAS","MEN1", "AIP", "DICER1", 
                "PRKAR1A", "PRKACA","SDH", "GPR101")


# most popular genes in human genome
# https://www.nature.com/articles/d41586-017-07291-9

top10 <- c("EGFR", "VEGFA", "IL6", "TGFB1", "MTHFR", "ESR1", "AKT1")

## stress genes from CORT paper

stressgenes <- c("A2ML1", "ABCA1", "ACOT12", "ANLN",
                 "APCDD1L", "APOD", "APOH", "APOLD1",
                 "CDH19", "CEBPD", "CISH",
                 "CITED2",  "CKMT1A", "CNP", "CREB5")

## genes from github issues 
## https://github.com/raynamharris/musicalgenes/blob/master/.github/ISSUE_TEMPLATE/request-a-gene.md
githubgenes <- c("BRCA1", "NPVF")

#candidategenes <- GOgenesLong %>% distinct(gene) %>% pull(gene)
candidategenes <- c(favoritegenes, curleychampagnegenes)
candidategenes <- unique(candidategenes) %>% sort(.) 
length(candidategenes)


# check that all canididates actually in transcriptome
geneids %>% filter(gene %in% candidategenes) %>%  arrange(gene) %>% pull(gene) %>% unique()

shinygenes <- c(favoritegenes, listoparentalcaregenes, 
                top10, breastcancergenes, githubgenes,
                stressgenes)
shinygenes <- unique(shinygenes) %>% sort(.) 
shinygenes
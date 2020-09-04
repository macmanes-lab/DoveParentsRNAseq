library(tidyverse)

# genelists

## all genes
geneids <- read_csv("metadata/00_geneinfo.csv") %>% select(-X1) 

favoritegenes <- c("AR", "AVP", "AVPR1B", "AVPR2", "CRH", "CRHR1",
                   "DIO1", "DIO2", "DIO3", "ESR2", "FSHB",
                   "FSHR", "GABRQ", "GALR1", "GNRH1", "GNRHR",
                   "JAK2", "LHCGR", "NPFFR1", "NPVF", "NPY",
                   "NR3C2", "OXT", "OXTR", "POMC", "PRL", "PRLH",	
                   "PRLHR", "SERPINA4", "STAT5A", "VIP", "VIPR1")

## parental care genes 
## genes from Ch 17 Evolution of parental care
curleychampagnegenes <- c("AVP", "AVPR1A", "ADRA2A",    
                          "COMT", "CRH", "CRHBP", "CRHR1", "CRHR2", 
                          "DRD1", "DRD4","ESR1", "ESR2", #ERa ERbeta
                          "FOS", "HTR2C", "MEST", "NR3C1", #GR
                          "OPRM1", "OXT" , "PGR", "PRL", "PRLR",  "SLC6A4") #5HTT

datadrivengenes <- c("COX1", "COX2", "CYTB", "ATP2B4", "LOC107055658", 
                     "CHGA", "TF", "OVALY", "ANXA5")

literaturegenes <- c(curleychampagnegenes, datadrivengenes,
                     "VIP", "GNIH", "GAL", "TH", "THRB",
                     "GNRHR", "GNRH1", "CGNRH-R", "FSHB", "FSHR")
literaturegenes

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


#####

## old gene lists


## breast and ovarian cancers  https://www.gynecologiconcology-online.net/article/S0090-8258(19)30069-1/fulltext
suszynskaagenes <- c("BRCA1","BRCA2", "CDKN2A", "PTEN", "PALB2", "TP53", "CDH1", "ATM",
                     "BARD1", "MESH6","MSH2", "BRIP1", "NBN", "FANCC", "FANCM", 
                     "RAD51C","RAD51D")

shaidgenes <- c("GNAS", "USB8", "PIK3CA", "GPR101","RAS","MEN1", "AIP", "DICER1", 
                "PRKAR1A", "PRKACA","SDH", "GPR101")

cancergenes <- c(suszynskaagenes, shaidgenes)

## genes WGCNA prl module
WGCNAgenes <- read_csv("results/05_PRLmodule.csv") %>% pull(x)


## all degs

degs <- read_csv("results/03_allDEG.csv") %>%
  filter(!grepl("LOC", gene),
         !grepl("control|bldg", comparison))  %>%
  distinct(gene) %>%
  arrange(gene) %>%
  pull(gene)
head(degs)
tail(degs)

#candidategenes <- GOgenesLong %>% distinct(gene) %>% pull(gene)
candidategenes <- c(favoritegenes)
shinygenes <- c(favoritegenes, cancergenes, parentalcaregenes)


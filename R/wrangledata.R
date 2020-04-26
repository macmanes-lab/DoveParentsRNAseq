# data wrangle for multiple figures

#### this is for tsne and pca
# Note, pseudocounts is too big for storage on github :(
pseudocounts <- read_csv("../results/01_pseudo.counts.csv")
head(pseudocounts[1:3])
pseudocounts <- as.data.frame(pseudocounts)
row.names(pseudocounts) <- pseudocounts$X1
pseudocounts$X1 <- NULL
# prep count data for all samples
countData <- as.data.frame(t(pseudocounts))
head(countData[1:3])

# this is for tsne and pca
# prep col data for all samples
colData <- read.csv("../metadata/00_samples.csv", header = T, row.names = 1)
colData$treatment <- factor(colData$treatment, levels = alllevels)
colData <- colData %>% mutate(tissue = fct_recode(tissue, "gonads" = "gonad"))
colData$tissue <- factor(colData$tissue, levels = tissuelevels)
row.names(colData) <- colData$V1

# check ready for analysis
#row.names(countData) == row.names(colData)


#### candidate gene anlayses

# all genes in this parental care dataset 
geneids <- read_csv("../metadata/00_geneinfo.csv") %>%
  select(-X1) %>%
  rename("NCBI" = "entrezid",
         "gene" = "Name")

# genes from Ch 17 Evolution of parental care
curleychampagnegenes <- c("AVP", "AVPR1A", "ADRA2A",    
                          "COMT", "CRH", "CRHBP", "CRHR1", "CRHR2", 
                          "DRD1", "DRD4","ESR1", "ESR2", #ERa ERbeta
                          "FOS", "HTR2C", "MEST", "NR3C1", #GR
                          "OPRM1", "OXT" , "PGR", "PRL", "PRLR",  "SLC6A4") #5HTT
curleychampagnegenes <- as.data.frame(curleychampagnegenes) %>%
  mutate(GO = "parentalcare", gene = curleychampagnegenes)  %>%
  select(GO, gene)
curleychampagnegenes


# candidate gens from parental care GO terms 
GO_path <- "../metadata/goterms/"   # path to the data
GO_files <- dir(GO_path, pattern = "*.txt") # get file names
GO_pathfiles <- paste0(GO_path, GO_files)

GOgenesLong <- GO_pathfiles %>%
  setNames(nm = .) %>% 
  map_df(~read_table(.x, col_types = cols(), col_names = FALSE), .id = "file_name") %>% 
  mutate(GO = sapply(strsplit(as.character(file_name),'../metadata/goterms/'), "[", 2)) %>% 
  mutate(GO = sapply(strsplit(as.character(GO),'.txt'), "[", 1)) %>% 
  mutate(gene = sapply(strsplit(as.character(X1), "[\\\\]|[^[:print:]]" ), "[", 2)) %>% 
  select(GO, gene)  %>%
  filter(gene != "Symbol") %>%
  distinct(GO,gene)  %>%
  mutate(gene = toupper(gene)) %>%
  rbind(., curleychampagnegenes) %>%
  arrange(gene) %>%
  left_join(., geneids, by = "gene") %>%
  arrange(gene) %>%
  drop_na() 
GOgenesLong



mamglanddev <- GOgenesLong %>% 
  filter(!GO %in% c("parentalbehavior", "parentalcare")) %>% 
  pivot_wider(
    names_from = GO,
    values_from = gene) %>%
  mutate(numGOs = 4 - rowSums(is.na(.))) %>%
  arrange(desc(numGOs)) %>%
  select(-numGOs) %>%
  left_join(geneids, by = c("geneid","NCBI")) %>%
  select(gene, mamglanddev, prostglanddev,  NCBI) %>%
  arrange(mamglanddev, prostglanddev)
mamglanddev


parentalcaregenes <- GOgenesLong %>% 
  filter(!GO %in% c("mamglanddev", "prostglanddev") ) %>% 
  pivot_wider(
    names_from = GO,
    values_from = gene) %>%
  mutate(numGOs = 4 - rowSums(is.na(.))) %>%
  arrange(desc(numGOs)) %>%
  select(-numGOs) %>%
  left_join(geneids, by = c("geneid","NCBI"))
parentalcaregenes

candidategenes <- parentalcaregenes %>% distinct(gene) %>% pull(gene)

# note: I can't find these genes in dataset: NOS1 or OXTR or FOX1B or HTR5A



#### variance stabilized gene expression  (vsd) for all and candidate genes
vsd_path <- "../results/DEseq2/"   # path to the data
vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
vsd_pathfiles <- paste0(vsd_path, vsd_files)
vsd_files

allvsd <- vsd_pathfiles %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x), .id = "file_name")  %>% 
  dplyr::rename("gene" = "X1") %>% 
  pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
               names_to = "samples", values_to = "counts") 

candidategenes <- GOgenesLong %>% distinct(gene) %>% pull(gene)
hypvsd <- getcandidatevsd(candidategenes, "hypothalamus", sexlevels)
pitvsd <- getcandidatevsd(candidategenes, "pituitary", sexlevels)
gonvsd <- getcandidatevsd(candidategenes, "gonad", sexlevels)
candidatevsd <- rbind(hypvsd, pitvsd)
candidatevsd <- rbind(candidatevsd, gonvsd)

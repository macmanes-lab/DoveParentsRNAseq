library(tidyverse)
library(readr)

source("../R/themes.R")

### candidate gene anlayses

# all genes in this parental care dataset 
geneids <- read_csv("../metadata/00_geneinfo.csv") %>% select(-X1) 

# genes from Ch 17 Evolution of parental care
curleychampagnegenes <- c("AVP", "AVPR1A", "ADRA2A",    
                          "COMT", "CRH", "CRHBP", "CRHR1", "CRHR2", 
                          "DRD1", "DRD4","ESR1", "ESR2", #ERa ERbeta
                          "FOS", "HTR2C", "MEST", "NR3C1", #GR
                          "OPRM1", "OXT" , "PGR", "PRL", "PRLR",  "SLC6A4") #5HTT
curleychampagnegenes
curleychampagnegenes <- as.data.frame(curleychampagnegenes) %>%
  mutate(GO = "parentalcare", gene = curleychampagnegenes)  %>%
  select(GO, gene)

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
head(GOgenesLong)

parentalcaregenes <- GOgenesLong %>% 
  pivot_wider(
    names_from = GO,
    values_from = gene) %>%
  mutate(numGOs = 4 - rowSums(is.na(.))) %>%
  arrange(desc(numGOs)) %>%
  select(-numGOs) %>%
  left_join(geneids, by = c("geneid","NCBI"))
head(parentalcaregenes)

candidategenes <- parentalcaregenes %>% distinct(gene) %>% pull(gene)

# note: I can't find these genes in dataset: NOS1 or OXTR or FOX1B or HTR5A


### variance stabilized gene expression  (vsd) 


vsd_path <- "../results/DEseq2/treatment/"   # path to the data
vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
vsd_pathfiles <- paste0(vsd_path, vsd_files)
vsd_files

## before pivoting, check names of df with
## head(names(allvsd)) and tail(names(allvsd))

allvsd <- vsd_pathfiles %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x), .id = "file_name")  %>% 
  dplyr::rename("gene" = "X1") %>% 
  pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
               names_to = "samples", values_to = "counts") 
allvsd %>% select(-file_name) %>% head()


### vsd for all and candidate genes

candidategenes <- GOgenesLong %>% distinct(gene) %>% pull(gene)

getcandidatevsd <- function(whichgenes, whichtissue, whichsex){
  candidates  <- allvsd %>%
    filter(gene %in% whichgenes) %>%
    dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
    dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/treatment/'), "[", 2)) %>%
    dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                  tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                  treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
    dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
    dplyr::select(sex, tissue, treatment, gene, samples, counts) %>%
    filter(tissue == whichtissue, sex %in% whichsex)  %>%
    drop_na()
  #candidates$treatment <- factor(candidates$treatment, levels = alllevels)
  return(candidates)
}


hypvsd <- getcandidatevsd(candidategenes, "hypothalamus", sexlevels)
pitvsd <- getcandidatevsd(candidategenes, "pituitary", sexlevels)
gonvsd <- getcandidatevsd(candidategenes, "gonad", sexlevels)
candidatevsd <- rbind(hypvsd, pitvsd)
candidatevsd <- rbind(candidatevsd, gonvsd)

head(candidatevsd)


### manip 	



DEG_path <- "../results/DEseq2/hypothesis/"   # path to the data	
DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names	
DEG_pathfiles <- paste0(DEG_path, DEG_files)	
#DEG_files	

allDEG2 <- DEG_pathfiles %>%	
  setNames(nm = .) %>% 	
  map_df(~read_csv(.x), .id = "file_name") %>% 	
  mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/hypothesis/'), "[", 2))  %>% 	
  mutate(DEG = sapply(strsplit(as.character(DEG),'_diffexp.csv'), "[", 1))  %>% 	
  mutate(tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1)) %>%	
  mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%	
  mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%	
  mutate(comparison = paste(down,up, sep = "_")) %>%	
  mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%	
  mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%	
  dplyr::select(sex,tissue,comparison, direction, gene, lfc, padj, logpadj) 	
head(allDEG2)	

allDEG2$tissue <- factor(allDEG2$tissue, levels = tissuelevel)	

allDEG2$comparison <- factor(allDEG2$comparison, levels = c("eggs_chicks", "lo_hi"))	
allDEG2 <- allDEG2 %>% mutate(comparison = fct_recode(comparison, "lo vs. hi PRL   " = "lo_hi",	
                                                      "eggs vs. chicks" = "eggs_chicks"))	
allDEG2$direction <- factor(allDEG2$direction, levels = c("eggs", "chicks", "lo", "hi"))	


PRLDEGs <- allDEG2 %>%	
  filter(tissue == "pituitary", comparison == "lo vs. hi PRL   ",	
         direction == "hi") %>%	
  arrange(desc(logpadj))	

envDEGs <- allDEG2 %>%	
  filter(tissue == "pituitary", comparison == "eggs vs. chicks",	
         direction == "chicks") %>%	
  arrange(desc(logpadj))	

hilogenes <- PRLDEGs %>% head(20) %>% distinct(gene) %>% pull(gene)	
eggchickgenes <- envDEGs %>% head(20) %>% distinct(gene) %>% pull(gene)	

hypothesisDEGs <- c(hilogenes, eggchickgenes) %>% unique()

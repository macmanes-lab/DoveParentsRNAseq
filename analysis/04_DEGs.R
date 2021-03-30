library(tidyverse)
library(readr)

source("R/themes.R")
source("R/genelists.R")

## All DEGs 

DEG_path <- "results/DEseq2/treatment/"   # path to the data
DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
DEG_pathfiles <- paste0(DEG_path, DEG_files)
#DEG_files

allDEG <- DEG_pathfiles %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x), .id = "file_name") %>% 
  mutate(DEG = sapply(strsplit(as.character(file_name),'results/DEseq2/treatment/'), "[", 2))  %>% 
  mutate(DEG = sapply(strsplit(as.character(DEG),'_diffexp.csv'), "[", 1))  %>% 
  mutate(tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1)) %>%
  mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%
  mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%
  mutate(comparison = paste(down,up, sep = "_")) %>%
  mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
  mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%
  dplyr::select(sex,tissue,comparison, gene, lfc, 
                padj, direction, pvalue, direction2,  logpadj)  %>%
  mutate(posneg = ifelse(lfc >= 0, "+", "-")) %>%
  drop_na()
head(allDEG)


## drop pvalues from allDEG

allDEG <- allDEG %>%
  #remove pvalue
  select(-pvalue, -direction2) %>% 
  filter(direction != "NS")

# make df of candidate genes to add NS genes
gene <- candidategenes
candidategenedf <- as.data.frame(gene) 

## subset only candidates
candidateDEGs <- allDEG %>%
  mutate(sex = recode(sex, "female" = "F", "male" = "M" ),
         tissue = recode(tissue, 
                         "hypothalamus" = "H",
                         "pituitary" = "P",
                         "gonad" = "G",
                         "gonads" = "G"),
         group = paste(sex, tissue, sep = "")) %>%
  filter(gene %in% candidategenes) %>%
  mutate(compres = paste(group, posneg, sep = "")) %>%
  group_by(gene, comparison) %>%
  summarize(results = str_c(compres, collapse = "  ")) %>%
  pivot_wider(names_from = comparison, values_from = results) %>%
  full_join(., candidategenedf) %>%
  arrange(gene)
head(candidateDEGs)

# make tables


# helper function to add column of NAs if no DEGs exist
fncols <- function(data, cname) {
  add <-cname[!cname%in%names(data)]
  if(length(add)!=0) data[add] <- NA
  data
}

# make table of candidate pairwise degs
makefinaltables <- function(df, whichlevels){
  table <- fncols(df, whichlevels) %>% 
    select(gene, whichlevels)  
  return(table)
}

table1 <- makefinaltables(candidateDEGs, levelsbldgchar) %>%
  select(gene, bldg_control:bldg_n9)
table2 <- makefinaltables(candidateDEGs, c(levelsrm, levelsreplace))

tableS1 <- makefinaltables(candidateDEGs, levelscontrolcharmanip) %>%
  select(gene, bldg_control:control_n9)
tableS2 <- makefinaltables(candidateDEGs, levelssequential) 
tableS3 <- makefinaltables(candidateDEGs, levelscontrolcharmanip) %>%
  select(gene, control_m.inc.d3:control_extend)


## save files
write.csv(allDEG, "results/04_allDEG.csv", row.names = F)
write.csv(allDEG, "results/04_candidateDEGs.csv", row.names = F)

write_csv(table1, "results/table1.csv")
write_csv(table2, "results/table2.csv")

write_csv(tableS1, "results/tableS1.csv")
write_csv(tableS2, "results/tableS2.csv")
write_csv(tableS3, "results/tableS3.csv")

## sex differences

hypsex <- read_csv("results/DEseq2/sex/pituitary_female_male_DEGs.csv") %>%
  rename(tissue = sextissue)
pitsex <- read_csv("results/DEseq2/sex/hypothalamus_female_male_DEGs.csv") %>%
  rename(tissue = sextissue)
gonsex <- read_csv("results/DEseq2/sex/gonad_female_male_DEGs.csv")

sharedsexdifs <- rbind(hypsex, pitsex) %>%
  mutate(res = paste("higher in the", direction, tissue, sep = " "))  %>%
  filter(gene %in% candidategenes) %>%
  select(gene, res) %>%
  group_by(gene) %>%
  summarize(res = str_c(res, collapse = "; ")) %>%
  group_by(res) %>%
  summarize(genes = str_c(gene, collapse = " ")) %>%
  mutate(n = str_count(genes, " "),
         n = n + 1) %>%
  select(n, res, genes) %>%
  arrange(desc(n))

head(sharedsexdifs) 
write.csv(sharedsexdifs, "results/sharedsexdifs.csv")

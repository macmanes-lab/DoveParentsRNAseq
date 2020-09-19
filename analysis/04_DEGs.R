library(tidyverse)
library(readr)

source("R/themes.R")
source("R/genelists.R")

## functions


# helper function to add column of NAs if no DEGs exist
fncols <- function(data, cname) {
  add <-cname[!cname%in%names(data)]
  if(length(add)!=0) data[add] <- NA
  data
}


# make table of candidate pairwise degs
maketable1 <- function(df, whichlevels){
  table <- fncols(df, whichlevels) %>% 
    select(gene, whichlevels)  %>%
    ungroup() %>%
    mutate(allres = rowSums(is.na(.))) 
  
  numcomparisons <- ncol(table) -2
  
  table <- table %>%
    filter(allres < numcomparisons) %>%
    select(-allres)
  print(table)
  return(table)
  
}

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
  dplyr::select(sex,tissue,comparison, gene, lfc, padj, direction, pvalue, direction2,  logpadj)  %>%
  mutate(posneg = ifelse(lfc >= 0, "+", "-"),
         sex = recode(sex, "female" = "F", "male" = "M" ),
         tissue = recode(tissue, 
                         "hypothalamus" = "H",
                         "pituitary" = "P",
                         "gonad" = "G",
                         "gonads" = "G"),
         group = paste(sex, tissue, sep = "")) %>%
  drop_na()
head(allDEG)

allDEG %>% filter(gene == "AR", tissue == "H", sex == "M")

candidateDEGs <- allDEG %>%
  #remove padj
  select(-padj, -direction) %>%
  filter(gene %in% candidategenes) %>%
  mutate(res = paste("(", posneg, ")", sep = "")) %>%
  mutate(compres = paste(group, res, sep = "")) %>%
  group_by(gene, comparison) %>%
  summarize(results = str_c(compres, collapse = "; ")) %>%
  pivot_wider(names_from = comparison, values_from = results) 

print(candidategenes)
print(head(candidateDEGs))

## drop pvalues from allDEG


allDEG <- allDEG %>%
  #remove pvalue
  select(-pvalue, -direction2) %>% 
  filter(direction != "NS")

# make tables

table1a <- maketable1(candidateDEGs, levelssequential) 
table1b <- maketable1(candidateDEGs, levelsrm) 
table1c <- maketable1(candidateDEGs, levelsreplace) 

tableS1 <- maketable1(candidateDEGs, levelscontrolcharmanip) %>%
  select(gene, control_bldg:control_n9)
tableS3 <- maketable1(candidateDEGs, levelscontrolcharmanip) %>%
  select(gene, control_m.inc.d3:control_extend)

tableS2a <- maketable1(candidateDEGs, levelsbldgcharmanip) %>%
  select(gene, bldg_lay:bldg_n9)
tableS2b <- maketable1(candidateDEGs, levelsbldgcharmanip) %>%
  select(gene, bldg_m.inc.d3:bldg_extend)

## save files
write.csv(allDEG, "results/04_allDEG.csv", row.names = F)
write.csv(candidateDEGs, "results/04_candidateDEGs.csv")

write_csv(table1a, "results/table1a.csv")
write_csv(table1b, "results/table1b.csv")
write_csv(table1c, "results/table1c.csv")

write_csv(tableS1, "results/tableS1.csv")
write_csv(tableS2a, "results/tableS2a.csv")
write_csv(tableS2b, "results/tableS2b.csv")
write_csv(tableS3, "results/tableS3.csv")
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
  dplyr::select(sex,tissue,comparison, direction, gene, lfc, padj, logpadj)  %>%
  mutate(posneg = ifelse(lfc >= 0, "+", "-"),
         sex = recode(sex, "female" = "F", "male" = "M" ),
         tissue = recode(tissue, 
                         "hypothalamus" = "H",
                         "pituitary" = "P",
                         "gonad" = "G",
                         "gonads" = "G"),
         group = paste(sex, tissue, sep = ""))
head(allDEG)

unique(allDEG$comparison)

print(candidategenes)

candidateDEGs <- allDEG %>%
  filter(gene %in% candidategenes) %>%
  mutate(res = paste("(", posneg, ")", sep = "")) %>%
  mutate(compres = paste(group, res, sep = "")) %>%
  group_by(gene, comparison) %>%
  summarize(results = str_c(compres, collapse = "; ")) %>%
  pivot_wider(names_from = comparison, values_from = results) 
head(candidateDEGs)

### manipuluation

makemanipulationDEGtables <- function(whichcomparisons){
  df <- allDEG %>%
    filter(comparison %in% whichcomparisons) %>%
    mutate(res = paste("(", posneg, ")", sep = "")) %>%
    mutate(compres = paste(group, res, sep = "")) %>%
    select(gene, group, compres) %>%
    group_by(group,gene) %>%
    summarize(results = str_c(compres, collapse = "; ")) %>%
    mutate(n = (str_count(results, pattern = ";"))+1) %>%
    arrange(desc(n))  %>%
    filter(n>1) %>%
    group_by(group,results) %>%
    summarize(genes = str_c(gene, collapse = "; "))  %>%
    #arrange(group, results) %>%
    mutate(n = (str_count(genes, pattern = ";"))+1) %>%
    arrange(desc(n))
  print(head(df))
  return(df) 
}

removalcomps <- c("inc.d3_m.inc.d3", "inc.d9_m.inc.d9", "inc.d17_m.inc.d17", "hatch_m.n2")
earlys <- c("inc.d9_early", "inc.d17_early", "hatch_early", "n5_early")
prolongs <- c("inc.d9_prolong","inc.d17_prolong", "hatch_prolong",  "n5_prolong")
extends <- c("inc.d9_extend","inc.d17_extend", "hatch_extend",  "n5_extend")

removalDEGs <- makemanipulationDEGtables(removalcomps)
earlyDEGs <- makemanipulationDEGtables(earlys)
prolongDEGs <- makemanipulationDEGtables(prolongs)
extendDEGs <- makemanipulationDEGtables(extends)

manipDEGs <- rbind(removalDEGs, earlyDEGs) %>%
  rbind(., prolongDEGs) %>%
  rbind(., extendDEGs) %>%
  arrange(desc(n))
head(manipDEGs)


# make tables

# helper function to add column of NAs if no DEGs exist
fncols <- function(data, cname) {
  add <-cname[!cname%in%names(data)]
  if(length(add)!=0) data[add] <- NA
  data
}


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

table1a <- maketables123(candidateDEGs, levelssequential) 
table1b <- maketables123(candidateDEGs, levelsrm) 
table1c <- maketables123(candidateDEGs, levelsreplace) 

## save files
write.csv(allDEG, "results/04_allDEG.csv", row.names = F)
write.csv(manipDEGs, "results/04_manipDEGs.csv")
write.csv(candidateDEGs, "results/04_candidateDEGs.csv")

write_csv(table1a, "results/table1a.csv")
write_csv(table1b, "results/table1b.csv")
write_csv(table1c, "results/table1c.csv")

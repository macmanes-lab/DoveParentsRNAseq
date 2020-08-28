library(tidyverse)
library(readr)
library(VennDiagram)

source("../R/themes.R")
source("../R/genelists.R")

## All DEGs 

DEG_path <- "../results/DEseq2/treatment/"   # path to the data
DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
DEG_pathfiles <- paste0(DEG_path, DEG_files)
#DEG_files

allDEG <- DEG_pathfiles %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x), .id = "file_name") %>% 
  mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/treatment/'), "[", 2))  %>% 
  mutate(DEG = sapply(strsplit(as.character(DEG),'_diffexp.csv'), "[", 1))  %>% 
  mutate(tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1)) %>%
  mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%
  mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%
  mutate(comparison = paste(down,up, sep = "_")) %>%
  mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
  mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)) %>%
  dplyr::select(sex,tissue,comparison, direction, gene, lfc, padj, logpadj) 
head(allDEG)

unique(allDEG$comparison)

tempDEGs <- DEG_pathfiles %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x), .id = "file_name") %>% 
  mutate(DEG = sapply(strsplit(as.character(file_name),'./results/DEseq2/treatment/'), "[", 2))  %>% 
  mutate(DEG = sapply(strsplit(as.character(DEG),'_diffexp.csv'), "[", 1))  %>% 
  mutate(tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1)) %>%
  mutate(down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3)) %>%
  mutate(up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4)) %>%
  mutate(comparison = paste(down,up, sep = "_")) %>%
  mutate(sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1)) %>%
  mutate(tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2))
head(tempDEGs)


# sequential in 2 or more comparisons
wideDEGscandidates <- tempDEGs %>%
  filter(gene %in% candidategenes) %>%
  mutate(posneg = ifelse(lfc >= 0, "+", "-"),
         sex = recode(sex, "female" = "F", "male" = "M" ),
         tissue = recode(tissue, 
                         "hypothalamus" = "H",
                         "pituitary" = "P", "gonad" = "G"),
         group = paste(sex, tissue, sep = "")) %>%
  mutate(res = paste("(", posneg, ")", sep = "")) %>%
  mutate(compres = paste(group, res, sep = "")) %>%
  group_by(gene, comparison) %>%
  summarize(results = str_c(compres, collapse = "; ")) %>%
  pivot_wider(names_from = comparison, values_from = results) 
wideDEGscandidates


### manipuluation


head(allDEG)

filteredDEGs <- allDEG %>%
  filter(lfc > 0.14 | lfc < -0.14 )
dim(allDEG)
dim(filteredDEGs)

makemanipulationDEGtables <- function(whichcomparisons){
  df <- filteredDEGs %>%
    filter(comparison %in% whichcomparisons) %>%
    mutate(posneg = ifelse(lfc >= 0, "+", "-"),
           sex = recode(sex, "female" = "F", "male" = "M" ),
           tissue = recode(tissue, 
                           "hypothalamus" = "H",
                           "pituitary" = "P", 
                           "gonad" = "G",
                           "gonads" = "G"),
           group = paste(sex, tissue, sep = "")) %>%
    mutate(res = paste("(", posneg, ")", sep = "")) %>%
    mutate(compres = paste(comparison, res, sep = "")) %>%
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


# for suppl figures
DEGcontrol <- allDEG %>% 
  filter(grepl("control", comparison),
         !grepl("m.|early|extend|prolong", comparison))  %>%
  mutate(comparison = factor(comparison, levels = comparisonlevelscontrol)) %>%
  group_by(sex, tissue, comparison, direction, label) %>%
  summarise(n = n()) %>%
  mutate(n = ifelse(direction == "control", n*-1, n*1 ))


DEGbldg <- allDEG %>% 
  filter(grepl("bldg", comparison),
         !grepl("m.|early|extend|prolong|control", comparison))  %>%
  drop_na() %>%
  mutate(comparison = factor(comparison, levels = comparisonlevelsbldg))  %>%
  group_by(sex, tissue, comparison, direction, label) %>%
  summarise(n = n()) %>%
  mutate(n = ifelse(direction == "bldg", n*-1, n*1 ))

DEGchar <- allDEG %>% 
  filter(comparison %in% comparisonlevelschar,
         !grepl("control|bldg", comparison)) %>%
  mutate(comparison = factor(comparison, levels = comparisonlevelschar)) %>%
  mutate(updown = ifelse(lfc > 0, 1, -1)) %>%
  group_by(sex, tissue, comparison, direction, label, updown) %>%
  summarise(n = n()) %>%
  mutate(n = n*updown ) %>%
  select(-updown)

DEGremove <- allDEG %>% 
  filter(comparison %in% comparisonlevelsremoval) %>%
  mutate(comparison = factor(comparison, levels = comparisonlevelsremoval)) %>%
  mutate(updown = ifelse(lfc > 0, 1, -1)) %>%
  group_by(sex, tissue, comparison, direction, label, updown) %>%
  summarise(n = n()) %>%
  mutate(n = n*updown ) %>%
  select(-updown)
DEGremove

DEGreplace <- allDEG %>% 
  filter(comparison %in% comparisonlevelsreplace) %>%
  mutate(comparison = factor(comparison, levels = comparisonlevelsreplace)) %>%
  mutate(updown = ifelse(lfc > 0, 1, -1)) %>%
  group_by(sex, tissue, comparison, direction, label, updown) %>%
  summarise(n = n()) %>%
  mutate(n = n*updown ) %>%
  select(-updown)
DEGreplace


DEGcontrolreplace <- allDEG %>% 
  filter(comparison %in% comparisonlevelscontrolreplace)  %>%
  mutate(comparison = factor(comparison, levels = comparisonlevelscontrolreplace)) %>%
  mutate(updown = ifelse(lfc > 0, 1, -1)) %>%
  group_by(sex, tissue, comparison, direction, label, updown) %>%
  summarise(n = n()) %>%
  mutate(n = n*updown ) %>%
  select(-updown)

## save files


write.csv(allDEG, "../results/03_allDEG.csv", row.names = F)
write.csv(candidatevsd, "../results/03_candidatevsd.csv")
write.csv(manipDEGs, "../results/03_manipDEGs.csv")
write.csv(wideDEGscandidates, "../results/03_wideDEGscandidates.csv", row.names = F)



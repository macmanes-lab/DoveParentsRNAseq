library(tidyverse)
source("../../R/themes.R")

# collaborator wants nest info...

# read bird data
birds <- read_csv("../../metadata/00_birds.csv") %>% select(-X1)

#wrangle nest data
nests <- read_csv("../../metadata/nestsummary.csv") %>%
  dplyr::rename("bird" = "AdultId_1BandNum",
                "nestmate" = "AdultId_2BandNum")  %>%
  mutate(pair = as.numeric(row.names(.)),
         bird = gsub("-", ".", bird),
         bird = gsub("/", ".", bird),
         bird = sapply(strsplit(bird, ' '), "[", 1),
         bird = sapply(strsplit(bird, '\\('), "[", 1),
         nestmate = gsub("-", ".", nestmate),
         nestmate = gsub("/", ".", nestmate),
         nestmate = sapply(strsplit(nestmate, ' '), "[", 1),
         nestmate = sapply(strsplit(nestmate, '\\C'), "[", 1)) 
  
neststemp1 <- nests %>%
  select(bird, nestmate, pair, c_BandDate) %>%
  drop_na()  

neststemp2 <- neststemp1 %>%
  select(nestmate, bird, pair, c_BandDate)
names(neststemp2) <-  c("bird", "nestmate", "pair", "c_BandDate") 

nestsnew <- rbind(neststemp1, neststemp2)

## join to get paired nests

nestpairs <- left_join(birds, nestsnew, by = "bird") %>%
  mutate(treatment = factor(treatment, levels = alllevels)) %>%
  arrange(treatment, pair) 
head(nestpairs)


write.csv(nestpairs, "../../metadata/nestpairs.csv", row.names = F)

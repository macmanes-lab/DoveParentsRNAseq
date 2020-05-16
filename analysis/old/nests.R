library(tidyverse)

# collaborator wants nest info...

# read bird data
birds <- read_csv("../../metadata/00_birds.csv") %>% select(-X1)

#wrangle nest data
neststemp1 <- read_csv("../../metadata/nestsummary.csv") %>%
  select(contains("AdultId")) %>%
  drop_na()  %>%
  dplyr::rename("bird" = "AdultId_1BandNum",
                "nestmate" = "AdultId_2BandNum") %>%
  mutate(pair = as.numeric(row.names(.)),
         bird = gsub("-", ".", bird),
         bird = gsub("/", ".", bird),
         bird = sapply(strsplit(bird, ' '), "[", 1),
         bird = sapply(strsplit(bird, '\\('), "[", 1),
         nestmate = gsub("-", ".", nestmate),
         nestmate = gsub("/", ".", nestmate),
         nestmate = sapply(strsplit(nestmate, ' '), "[", 1),
         nestmate = sapply(strsplit(nestmate, '\\C'), "[", 1)) 

neststemp2 <- neststemp1 %>%
  select(nestmate, bird, pair)
names(neststemp2) <-  c("bird", "nestmate", "pair") 

nests <- rbind(neststemp1, neststemp2)

## join to get paired nests

nestpairs <- left_join(birds, nests, by = "bird") %>%
  arrange(treatment, pair)
head(nestpairs)

write.csv(nestpairs, "../../metadata/nestpairs.csv", row.names = F)

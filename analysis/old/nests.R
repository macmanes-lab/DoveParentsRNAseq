library(tidyverse)

# collaborator wants nest info...

# read bird data
birds <- read_csv("../../metadata/00_birds.csv") %>% select(-X1)
head(birds)

#wrangle nest data
neststemp1 <- read_csv("../../metadata/nestsummary.csv") %>%
  select(contains("AdultId")) %>%
  drop_na()  %>%
  dplyr::rename("bird" = "AdultId_1BandNum",
                "nestmate" = "AdultId_2BandNum") %>%
  mutate(pair = row.names(.),
         bird = gsub("-", ".", bird),
         bird = gsub("/", ".", bird),
         nestmate = gsub("-", ".", nestmate),
         nestmate = gsub("/", ".", nestmate)) 

neststemp2 <- neststemp1 %>%
  select(nestmate, bird, pair)
names(neststemp2) <-  c("bird", "nestmate", "pair") 

nests <- rbind(neststemp1, neststemp2)

## join to get paired nests

nestpairs <- left_join(birds, nests, by = "bird")
head(nestpairs)

write.csv(nestpairs, "../../metadata/nestpairs.csv", row.names = F)

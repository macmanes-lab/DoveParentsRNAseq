source("R/themes.R")
library(tidyverse)

cDEGs <- read_csv("results/03_wideDEGscandidates.csv") %>%
  select(-X1)


mycomparsisons <- c(levelssequential)

# helper function to add column of NAs if no DEGs exist
fncols <- function(data, cname) {
  add <-cname[!cname%in%names(data)]
  if(length(add)!=0) data[add] <- NA
  data
}

table1 <- fncols(cDEGs, mycomparsisons) %>%
  select(gene,levelssequential) 
table1

write_csv(table1, "results/table1.csv")

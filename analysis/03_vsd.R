library(tidyverse)
library(corrr)

samples <- read_csv("metadata/00_birds_sachecked.csv") %>%
  mutate(id = bird) %>%
  select(id, horm.id)
head(samples)

hormones <- read_csv("results/AllHorm_02042020_nomiss.csv") %>%
  mutate(id = str_replace_all(id, c("-" = ".", "/" = "."))) %>%
  mutate(prl = as.numeric(prl),
         cort = as.numeric(cort),
         p4 = as.numeric(p4),
         e2 = as.numeric(e2),
         t = as.numeric(t)) %>%
  full_join(., samples) %>%
  select(id, treatment, sex, prl, cort, p4, e2, t) %>%
  filter(treatment != "NA", treatment != ".")
tail(hormones)

femalegenesnhomrmones <- function(vsds){
  df <- vsds %>% 
    inner_join(hormones, ., by = "id") %>%
    filter(sex == "f") %>%
    select(-t) 
  print(head(df))
  return(df)
}

malegenesnhomrmones <- function(vsds){
  df <-  vsds %>% 
    inner_join(hormones, ., by = "id") %>%
    filter(sex == "m") %>%
    select(-e2)
  print(head(df))
  return(df)
}

pitvsdAll <- read_csv("results/03_pitvsdAll.csv")
hypvsdAll <- read_csv("results/03_hypvsdAll.csv")
gonvsdAll <- read_csv("results/03_gonvsdAll.csv")

pitf <- femalegenesnhomrmones(pitvsdAll)
gonf <- femalegenesnhomrmones(hypvsdAll)
hypf <- femalegenesnhomrmones(gonvsdAll) 

pitm <- malegenesnhomrmones(pitvsdAll)
gonm <- malegenesnhomrmones(hypvsdAll)
hypm <- malegenesnhomrmones(gonvsdAll) 




plotcorrelation <- function(df, xvar, yvar, xlab, ylab){
  df %>%
    ggplot(aes(x = xvar, y = yvar, color = sex)) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_y_log10() +
    scale_x_log10() +
    labs(x = xlab, y = ylab) 
}

plotcorrelation(pitm,
                pitm$prl, pitm$PRL,
                "circulating PRL", "PRL expression")

plotcorrelation(pitf,
                pitf$prl, pitf$PRL,
                "circulating PRL", "PRL expression")

calculatecorrs <- function(df){
  df2 <- df %>%
    select(-id, -sex, -treatment) %>%
    correlate(.)
  print(head(df2)[1:6])
  return(df2)
}

pitfcorrs <- calculatecorrs(pitf)
pitmcorrs <- calculatecorrs(pitm)

library(tidyverse)
library(corrr)

# sample information ----
samples <- read_csv("metadata/00_birds_sachecked.csv") %>%
  mutate(id = bird) %>%
  select(id, horm.id)
head(samples)


# hormones ---
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

# vsds (gene expression) ----
pitvsdAll <- read_csv("results/03_pitvsdAll.csv")
hypvsdAll <- read_csv("results/03_hypvsdAll.csv")
gonvsdAll <- read_csv("results/03_gonvsdAll.csv")

# join vsds and hormones ----
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

# genes and hormones by sex and tissue ---
pitf <- femalegenesnhomrmones(pitvsdAll)
gonf <- femalegenesnhomrmones(hypvsdAll)
hypf <- femalegenesnhomrmones(gonvsdAll) 

pitm <- malegenesnhomrmones(pitvsdAll)
gonm <- malegenesnhomrmones(hypvsdAll)
hypm <- malegenesnhomrmones(gonvsdAll) 



# correlations ---
plotcorrelation <- function(df, xvar, yvar, xlab, ylab){
  df %>%
    ggplot(aes(x = xvar, y = yvar, color = sex)) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_y_log10() +
    scale_x_log10() +
    labs(x = xlab, y = ylab) 
}

calculatecorrs <- function(df){
  df2 <- df %>%
    select(-id, -sex, -treatment) %>%
    correlate(.)
  print(head(df2)[1:6])
  return(df2)
}

plotcorrelation(pitm,
                pitm$prl, pitm$PRL,
                "circulating PRL", "PRL expression")

plotcorrelation(pitf,
                pitf$prl, pitf$PRL,
                "circulating PRL", "PRL expression")

plotcorrelation(pitf,
                pitf$cort, pitf$PRL,
                "circulating cort", "PRL expression")


pitfcorrs <- calculatecorrs(pitf)
pitmcorrs <- calculatecorrs(pitm)

pitfcorrs %>% focus(PRL) %>% arrange(PRL)
pitfcorrs %>% focus(PRL) %>% arrange(desc(PRL))
pitmcorrs %>% focus(PRL) %>% arrange(PRL)
pitmcorrs %>% focus(PRL) %>% arrange(desc(PRL))

pitfcorrs %>% focus(prl) %>% arrange(prl)
pitfcorrs %>% focus(prl) %>% arrange(desc(prl))
pitmcorrs %>% focus(prl) %>% arrange(prl)
pitmcorrs %>% focus(prl) %>% arrange(desc(prl))

pitfcorrs %>% focus(cort) %>% arrange(cort)
pitfcorrs %>% focus(cort) %>% arrange(desc(cort))
pitmcorrs %>% focus(cort) %>% arrange(cort)
pitmcorrs %>% focus(cort) %>% arrange(desc(cort))

pitfcorrs %>% focus(e2) %>% arrange(e2)
pitfcorrs %>% focus(e2) %>% arrange(desc(e2))
pitmcorrs %>% focus(t) %>% arrange(t)
pitmcorrs %>% focus(t) %>% arrange(desc(t))

gonfcorrs <- calculatecorrs(gonf)
gonmcorrs <- calculatecorrs(gonm)

gonfcorrs %>% focus(e2) %>% arrange(e2)
gonfcorrs %>% focus(e2) %>% arrange(desc(e2))
gonmcorrs %>% focus(t) %>% arrange(t)
gonmcorrs %>% focus(t) %>% arrange(desc(t))

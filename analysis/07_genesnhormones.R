library(tidyverse)
library(corrr)
library(cowplot)
library(Hmisc)

source("R/genelists.R")
source("R/themes.R")

# sample information ----
samples <- read_csv("metadata/00_birds_sachecked.csv") %>%
  mutate(id = bird) %>%
  select(id, horm.id)
head(samples)



# hormones ---
hormones <- read_csv("results/AllHorm_02042020_nomiss.csv") %>%
  rename("horm.id" = "id") %>%
  mutate(prl = as.numeric(prl),
         cort = as.numeric(cort),
         p4 = as.numeric(p4),
         e2 = as.numeric(e2),
         t = as.numeric(t)) %>%
  full_join(samples, ., by = "horm.id") %>%
  arrange(desc(studytype)) %>%
  filter(treatment != "NA", treatment != ".", treatment != "stress") %>%
  mutate(e2t = ifelse(sex == "f", e2, t)) %>%
  mutate(sex = recode(sex, "f" = "female", "m" = "male" )) %>%
  select(id, treatment, sex, prl, cort, p4, e2t) %>%
  mutate(treatment = ifelse(treatment == "inc17", "inc.d17",
                            ifelse(treatment == "inc9", "inc.d9",
                                   ifelse(treatment == "inc3", "inc.d3",
                                          ifelse(treatment == "m.inc17", "m.inc.d17",
                                                 ifelse(treatment == "m.inc9", "m.inc.d9",
                                                        ifelse(treatment == "m.inc3", "m.inc.d3",
                                          ifelse(grepl("m.hatch", treatment), "m.n2",
                                                 ifelse(treatment == "m.inc8", "early",
                                                  treatment)))))))))



write_csv(hormones, "../musicalgenes/data/hormones.csv")

# vsds (gene expression) ----


vsd_path <- "results/"   # path to the data
vsd_files <- c("03_gonvsdf.csv" ,  "03_gonvsdm.csv" ,
               "03_hypvsdm.csv" ,   "03_hypvsdf.csv",
               "03_pitvsdf.csv" ,     "03_pitvsdm.csv" )
vsd_pathfiles <- paste0(vsd_path, vsd_files)
vsd_pathfiles

candidatevsds <- vsd_pathfiles %>%
  map_dfr(~read_csv(.x)) %>%
  mutate(id = sapply(strsplit(as.character(samples),'\\_'), "[", 1)) %>%
  select(id, sex, tissue, treatment, gene, counts) %>%
  pivot_wider(names_from = "gene", values_from = "counts") 
head(candidatevsds)

unique(candidatevsds$treatment)

# join genes and hormones ----
genesnhomrmones <- full_join(hormones, candidatevsds, by = "id")  %>%
  select(id, sex, tissue, treatment, 
         prl, cort, p4, e2t, everything()) %>%
  mutate(tissue = factor(tissue, levels = tissuelevels))
head(genesnhomrmones)

genesnhomrmones$treatment.x

x <- genesnhomrmones %>% filter_(~sex.x != sex.y) %>%
  select(id, treatment.x, treatment.y, sex.x, sex.y) %>%
  distinct()
x

# correlations ---

plotcorrelation <- function(xvar, xlab, yvar, ylab){
  p <- genesnhomrmones %>%
    ggplot(aes(x = xvar, y = yvar)) +
    geom_point(aes(color = treatment)) +
    geom_smooth(method = "lm", aes(color = sex)) +
    scale_y_log10() +
    scale_x_log10() +
    labs(x = xlab, y = ylab) +
    scale_color_manual(values = allcolors) +
    facet_wrap(~tissue, ncol = 1, scales = "free", switch = T) +
    theme(legend.position = "none")
  p
}


plotcorrelation(genesnhomrmones$prl, "Circulating prolactin (ng/mL)",
                genesnhomrmones$PRL, "Prolactin (PRL) expression")

plotcorrelation(genesnhomrmones$COMT, "COMT",
                genesnhomrmones$DRD1, "DRD1")

plotcorrelation(genesnhomrmones$OPRM1, "OPRM1",
                genesnhomrmones$PGR, "PGR")


cor.test(genesnhomrmones$prl, genesnhomrmones$PRL, 
                method = "pearson")


# get all correlations for each sex tissue without stats

corrs <- genesnhomrmones %>%
  filter(sex == "female", tissue == "hypothalamus") %>%
  select(-id, -sex, -tissue, -treatment, 
         -NPVF, -MC3R, -VIP, -CRH) %>%
  cor(.) %>%  
  as.data.frame(.) %>%
  rownames_to_column(., var = "var1") %>%
  pivot_longer(cols = prl:VIPR1, 
               names_to = "var2", values_to = "R2")  %>%
  mutate(sex = "female", tissue = "hypothalamus") %>%
  filter(R2 != 1) %>%
  arrange(desc(R2))
corrs


# get all correlations for each sex tissue with stats

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}




makecortable <- function(){
  
  d <- c()
  
  for(i in sexlevels){
    for(j in tissuelevels){
      
      res2 <- genesnhomrmones %>%
        filter(sex == i, tissue == j) %>%
        select(-id, -sex, -tissue, -treatment, 
               -NPVF, -MC3R, -VIP, -CRH)  %>%
        as.matrix()
      res2 <- rcorr(res2)
      
      flatres2 <- flattenCorrMatrix(res2$r, res2$P)
      flatres2 <- as.tibble(flatres2) %>%
        arrange(desc(cor,p)) %>%
        #filter(p < 0.01) %>%
        mutate(sex = i, tissue = j)
      flatres2   
      
      d <- rbind(d, flatres2)
    }
  }
  d
}

allcors <- makecortable()  %>%
  mutate(direction = if_else(cor > 0, 
                             "positive", "negative")) %>%
  mutate(tissue = factor(tissue, levels = tissuelevels))

allcors  %>%
  ggplot(aes(x = tissue, y = cor, fill = sex)) +
  geom_boxplot() +
  facet_wrap(~direction)

top50 <- allcors %>%
  filter(cor > 0.5)

bottom50 <- allcors %>%
  filter(cor < -0.5)

ggplot(allcors, aes(x = cor, fill = sex)) +
  geom_histogram() +
  facet_wrap(sex ~tissue)

allcors %>%
  filter( row == "prl") %>%
  filter(p < 0.01)

allcors %>%
  filter( row == "cort") %>%
  filter(p < 0.01)

allcors %>%
  filter( row == "e2t") %>%
  filter(p < 0.01)

allcors %>%
  filter( row == "p4") %>%
  filter(p < 0.01)

plotcorrelation(genesnhomrmones$cort, "Circulating corticosterone (ng/mL)",
                genesnhomrmones$FOS, "FOS expression")


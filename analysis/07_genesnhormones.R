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
glimpse(samples)


# hormones ---
hormones <- read_csv("results/AllHorm_02042020_nomiss.csv") %>%
  mutate(id = str_replace_all(id, c("-" = ".", "/" = "."))) %>%
  mutate(prl = as.numeric(prl),
         cort = as.numeric(cort),
         p4 = as.numeric(p4),
         e2 = as.numeric(e2),
         t = as.numeric(t)) %>%
  full_join(., samples) %>%
  filter(treatment != "NA", treatment != ".") %>%
  mutate(e2t = ifelse(sex == "f", e2, t)) %>%
  mutate(sex = recode(sex, "f" = "female", "m" = "male" )) %>%
  select(id, treatment, sex, prl, cort, p4, e2t)   %>%
  filter(treatment != "control")
head(hormones)

#write_csv(hormones, "../musicalgenes/data/hormones.csv")

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
  filter(id != "x.g.g.ATLAS") %>%
  select(id, tissue, gene, counts) %>%
  pivot_wider(names_from = "gene", values_from = "counts")
head(candidatevsds)

# join genes and hormones ----
genesnhomrmones <- inner_join(hormones, candidatevsds, by = "id")  %>%
  select(id, sex, tissue, treatment, 
         prl, cort, p4, e2t, everything()) %>%
  mutate(tissue = factor(tissue, levels = tissuelevels))
head(genesnhomrmones)

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


mkFrameForLoop <- function(nRow,nCol) {
  d <- c()
  for(i in seq_len(nRow)) {
    ri <- mkRow(nCol)
    di <- data.frame(ri,
                     stringsAsFactors=FALSE)
    d <- rbind(d,di)
  }
  d
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

sigcors <- makecortable()  %>%
  mutate(direction = if_else(cor > 0, 
                             "positive", "negative")) %>%
  mutate(tissue = factor(tissue, levels = tissuelevels))

sigcors  %>%
  ggplot(aes(x = tissue, y = cor, fill = sex)) +
  geom_boxplot() +
  facet_wrap(~direction)

top70 <- sigcors %>%
  filter(cor > 0.7)

bottom70 <- sigcors %>%
  filter(cor < -0.7)

ggplot(sigcors, aes(x = cor, fill = sex)) +
  geom_histogram() +
  facet_wrap(sex ~tissue)

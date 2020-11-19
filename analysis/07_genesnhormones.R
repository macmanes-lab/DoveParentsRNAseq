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
hormones <- read_csv("results/AllHorm_02042020_nomiss.csv",
                    # skip last line  
                    n_max = 351
                    ) %>%
  rename("horm.id" = "id") %>%
  mutate(prl = as.numeric(prl),
         cort = as.numeric(cort),
         p4 = as.numeric(p4),
         e2 = as.numeric(e2),
         t = as.numeric(t)) %>%
  full_join(samples, ., by = "horm.id") %>%
  arrange(desc(studytype)) %>%
  filter(treatment != "NA", treatment != "stress") %>%
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
                                                  treatment)))))))),
         treatment = ifelse(id == "pk.s054.d.g", "m.n2", treatment),
         sex = ifelse(id == "pk.s054.d.g", "female", sex),
         sex = ifelse(id == "blu119.w84.x", "female", sex))  %>%
  drop_na()
head(hormones)

## fix bad samples in hormones
# pk.s054.d.g is a female m.hatch (n2) bad, drop
# blu119.w84.x is male and female in other ds



# fix bad  samples in gene exrpress

# x.o61 is a female extend fixed
# y128.g23.x is a female m.inc.d9 fixed
# blu108.w40.o158 is a male m.inc.d9 fixed
# x.blu43.g132 is a female m.inc.d9 fixed
# r37.w100.x is a male m.inc.d9 fixed

# remember that x.g.g.ATLAS has same name for males and females


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
  #filter(!gene %in% c("CRH" , "MC3R", "NPVF")) %>%
  pivot_wider(names_from = "gene", values_from = "counts")  %>%
  mutate(treatment = ifelse(grepl("x.o61|blu84.x", id), "extend",
                            ifelse(grepl("y128.g23.x|blu108.w40.o158|x.blu43.g132|r37.w100.x", 
                                         id), "m.inc.d9",
                            treatment)))
head(candidatevsds)


# made new file for shiny app

candidatecounts <- candidatevsds %>%
  pivot_longer(cols = ADRA2A:MC3R, 
               names_to = "gene", values_to = "counts") %>%
  mutate(samples = paste(id, sex, tissue, treatment, sep = "_")) %>%
  select(gene, counts, sex, tissue, treatment, samples)
head(candidatecounts)


# join genes and hormones ----

genesnhomrmones <- full_join(hormones, candidatevsds, by = c("id", "sex", "treatment"))  %>%
  select(id, sex, tissue, treatment, 
         prl, cort, p4, e2t, everything()) %>%
  mutate(tissue = factor(tissue, levels = tissuelevels),
         treatment = factor(treatment, levels = alllevels))
head(genesnhomrmones)

samplesfixed <- genesnhomrmones %>%
  select(id, sex, treatment) %>%
  distinct() %>% arrange(id)
head(samplesfixed)

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


## plot correlations

plotcorrelation(genesnhomrmones$prl, "Circulating prolactin (ng/mL)",
                genesnhomrmones$PRL, "Prolactin (PRL) expression")

plotcorrelation(genesnhomrmones$COMT, "COMT",
                genesnhomrmones$DRD1, "DRD1")

plotcorrelation(genesnhomrmones$OPRM1, "OPRM1",
                genesnhomrmones$PGR, "PGR")


a <- plotcorrelation(genesnhomrmones$HTR2C, "HTR2C",
                genesnhomrmones$DRD1, "DRD1")

b <- plotcorrelation(genesnhomrmones$POMC, "POMC",
                genesnhomrmones$DRD1, "DRD1")

plot_grid(a,b)


## sig correlatoins 

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
        select(-id, -sex, -tissue, -treatment)  %>%	
        as.matrix()	
      res2 <- rcorr(res2)	
      
      flatres2 <- flattenCorrMatrix(res2$r, res2$P)	
      flatres2 <- as_tibble(flatres2) %>%	
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

head(allcors)

sigcors <- allcors %>%
  filter(p < 0.01)

sigcors  %>%	
  ggplot(aes(x = tissue, y = cor, fill = sex)) +	
    geom_boxplot() +	
    facet_wrap(~direction)

sigcors %>%
  ggplot(aes(x = cor, fill = sex)) +	
    geom_histogram() +	
    facet_wrap(sex ~tissue)	

top50 <- sigcors %>%	
  filter(cor > 0.5)	

bottom50 <- sigcors %>%	
  filter(cor < -0.5)	

topbottom50 <- rbind(top50,bottom50) %>%
  arrange(row,column)
topbottom50

hormonecorrs <- sigcors %>%
  filter(row %in% c("prl", "cort", "p4", "e2t"))  %>%
  arrange(row,column)

write_csv(topbottom50, "results/07_topsigcorrsover50.csv")
write_csv(hormonecorrs, "results/07_hormonecorrs.csv")


plotcorrelation(genesnhomrmones$cort, "Circulating corticosterone (ng/mL)",	
                genesnhomrmones$FOS, "FOS expression")

## save file for musical genes 

write_csv(samplesfixed, "metadata/00_birds_fixed.csv")
write_csv(hormones, "../musicalgenes/data/hormones.csv")
write_csv(candidatecounts, "../musicalgenes/data/candidatecounts.csv")

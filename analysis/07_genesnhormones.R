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
                                                  treatment))))))))) %>%
  filter(!id %in% c("pk.s054.d.g", "blu119.w84.x", "blu84.x", "x.g.g.ATLAS")) 
head(hormones)

## fix bad samples in hormones
# pk.s054.d.g is a female m.hatch (n2) bad, drop
# blu119.w84.x is male and female in other ds


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
  pivot_wider(names_from = "gene", values_from = "counts")  %>%
  mutate(treatment = ifelse(id == "x.o61", "extend",
                            ifelse(grepl("y128.g23.x|blu108.w40.o158|x.blu43.g132|r37.w100.x", id), "m.inc.d9",
                            treatment))) %>%
  filter(!id %in% c("pk.s054.d.g", "blu119.w84.x", "blu84.x", "x.g.g.ATLAS") ) 
head(candidatevsds)

# fix bad hormoen samples

# x.o61 is a female extend fixed
# y128.g23.x is a female m.inc.d9 fixed
# blu108.w40.o158 is a male m.inc.d9 fixed
# x.blu43.g132 is a female m.inc.d9 fixed
# r37.w100.x is a male m.inc.d9 fixed

# remove x.g.g.ATLAS has same name for males and females
# rmove blu119.w84.x has same name for males and females


# join genes and hormones ----

genesnhomrmones <- inner_join(hormones, candidatevsds, by = "id")  %>%
  rename("treatment" = "treatment.x", "sex" = "sex.x" ) %>%
  select(-treatment.y, sex.y) %>%
  select(id, sex, tissue, treatment, 
         prl, cort, p4, e2t, everything()) %>%
  mutate(tissue = factor(tissue, levels = tissuelevels),
         treatment = factor(treatment, levels = alllevels))
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


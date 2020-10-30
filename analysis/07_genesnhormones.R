library(tidyverse)
library(corrr)
library(cowplot)

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
  select(id, treatment, sex, prl, cort, p4, e2t)  
head(hormones)

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


plotcorrelation(genesnhomrmones$e2t, "Circulating E2 or T (ng/mL)",
                genesnhomrmones$PRLH, "PRLH expression")

plotcorrelation(genesnhomrmones$prl, "Circulating prolactin (ng/mL)",
                genesnhomrmones$PRL, "Prolactin (PRL) expression")

plotcorrelation(genesnhomrmones$prl, "Circulating prolactin (ng/mL)",
                genesnhomrmones$PRL, "Prolactin (PRL) expression")

plot_grid(a,b,c, ncol = 1)



plotcorrrs <- function(whichtissue, whichsex, mytitle){
  p <- genesnhomrmones %>%
    filter(tissue == whichtissue,
           sex %in% whichsex) %>%
    select(-id, -sex, -tissue,-treatment) %>%
    correlate(.) %>% 
    select(rowname, prl, cort, p4, e2t) %>%
    filter(!rowname %in% c("prl", "cort", "p4", "e2t")) %>%
    rplot(colors = 
            c("skyblue1", "white","indianred2")) +
    theme_minimal(base_size = 9) +
    theme(axis.text.y = element_blank(), 
          legend.position = "none") +
    labs(subtitle = mytitle)
  return(p)
}

a <- plotcorrrs("hypothalamus", sexlevels, "Hypothalamus") + 
  theme(axis.text.y = element_text(face = "italic"))
b <- plotcorrrs("pituitary", sexlevels, "Pituitary")
c <- plotcorrrs("gonads", sexlevels, "Gonads") + 
  theme(legend.position = "right")

p <- plot_grid(a,b,c, nrow = 1, rel_widths = c(1.2,1,1.2))
p

e <- plotcorrrs("hypothalamus", "female", "F Hyp") + 
  theme(axis.text.y = element_text(face = "italic"))
f <- plotcorrrs("pituitary", "female", "F Pit")
g <- plotcorrrs("gonads", "female", "F Gon") 

h <- plotcorrrs("hypothalamus", "male", "M Hyp") 
i <- plotcorrrs("pituitary", "male", "M Pit")
j <- plotcorrrs("gonads", "male", "M Gon") + 
  theme(legend.position = "right") 
p2 <- plot_grid(e,h,f,i,g,j, nrow = 1, rel_widths = c(1.2,1,1,1,1,1,1.2))
p2





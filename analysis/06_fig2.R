# Fig 2

library(tidyverse)
library(cowplot)
library(ggsignif)

source("R/themes.R")
source("R/functions.R")
source("R/genelists.R")

wranglevsds <- function(pathtofile){
  df <- read_csv(pathtofile) %>%
    mutate(treatment = factor(treatment, levels = alllevels))
  return(df)
}

hypvsdf <- wranglevsds("results/03_hypvsdf.csv")
pitvsdf <- wranglevsds("results/03_pitvsdf.csv")
gonvsdf <- wranglevsds("results/03_gonvsdf.csv")
hypvsdm <- wranglevsds("results/03_hypvsdm.csv")
pitvsdm <- wranglevsds("results/03_pitvsdm.csv")
gonvsdm <- wranglevsds("results/03_gonvsdm.csv")

print("done with vsds")

allDEG <- read_csv("results/03_allDEG.csv")


DEGcontrol <- allDEG %>% 
  filter(grepl("control", comparison),
         !grepl("m.|early|extend|prolong", comparison))  %>%
  mutate(comparison = factor(comparison, levels = comparisonlevelscontrol)) %>%
  group_by(sex, tissue, comparison, direction) %>%
  summarise(n = n()) %>%
  mutate(n = ifelse(direction == "control", n*-1, n*1 ))


DEGbldg <- allDEG %>% 
  filter(grepl("bldg", comparison),
         !grepl("m.|early|extend|prolong|control", comparison))  %>%
  drop_na() %>%
  mutate(comparison = factor(comparison, levels = comparisonlevelsbldg))  %>%
  group_by(sex, tissue, comparison, direction) %>%
  summarise(n = n()) %>%
  mutate(n = ifelse(direction == "bldg", n*-1, n*1 ))

DEGchar <- allDEG %>% 
  filter(comparison %in% comparisonlevelschar,
         !grepl("control|bldg", comparison)) %>%
  mutate(comparison = factor(comparison, levels = comparisonlevelschar)) %>%
  mutate(updown = ifelse(lfc > 0, 1, -1)) %>%
  group_by(sex, tissue, comparison, direction, updown) %>%
  summarise(n = n()) %>%
  mutate(n = n*updown ) %>%
  select(-updown)

DEGremove <- allDEG %>% 
  filter(comparison %in% comparisonlevelsremoval) %>%
  mutate(comparison = factor(comparison, levels = comparisonlevelsremoval)) %>%
  mutate(updown = ifelse(lfc > 0, 1, -1)) %>%
  group_by(sex, tissue, comparison, direction, updown) %>%
  summarise(n = n()) %>%
  mutate(n = n*updown ) %>%
  select(-updown)
DEGremove

DEGreplace <- allDEG %>% 
  filter(comparison %in% comparisonlevelsreplace) %>%
  mutate(comparison = factor(comparison, levels = comparisonlevelsreplace)) %>%
  mutate(updown = ifelse(lfc > 0, 1, -1)) %>%
  group_by(sex, tissue, comparison, direction, updown) %>%
  summarise(n = n()) %>%
  mutate(n = n*updown ) %>%
  select(-updown)
DEGreplace


DEGcontrolreplace <- allDEG %>% 
  filter(comparison %in% comparisonlevelscontrolreplace)  %>%
  mutate(comparison = factor(comparison, levels = comparisonlevelscontrolreplace)) %>%
  mutate(updown = ifelse(lfc > 0, 1, -1)) %>%
  group_by(sex, tissue, comparison, direction, updown) %>%
  summarise(n = n()) %>%
  mutate(n = n*updown ) %>%
  select(-updown)


a1 <- plotcandidatechar(hypvsdf, "HTR2C") + labs(subtitle = "Hypothalamus")
a2 <- plotcandidatechar(hypvsdm, "HTR2C") + labs(y = NULL, x = "") + labs(subtitle = " ")

b1 <- plot.volcano("hypothalamus", sexlevels,  "control_hatch") + 
  facet_wrap(~sex) 
b2 <- plot.volcano("hypothalamus", sexlevels,  "hatch_n5") + 
  facet_wrap(~sex) 

ab <- plot_grid(a1,a2,b1,b2, nrow = 1,
                hjust = 0,
                 rel_widths = c(2,2,1,1),
                labels = c("A", "B"), label_size = 8)

c <- makebargraphv3(DEGcontrolreplace, "hypothalamus","No. of DEGs\n(-) decreased  increased (+)", comparisonlabelscontrolreaplce, comparisonlevelscontrolreplace)
  labs(x = "Control versus all other reproductive and parental stages") 

d <- makebargraphv3(DEGbldg, "hypothalamus", "No. of DEGs\n (-)decreased  increased (+)", comparisonlabelsbldg, comparisonlevelsbldg) +
  labs(x = "Nest-building versus all other parental stages")
e <- makebargraphv3(DEGchar, "hypothalamus", NULL,  comparisonlabelscharnobldg,
                    comparisonlevelscharnobldg) +
  labs(x = "Comparison of sequential parental stages")

de <- plot_grid(d,e,
                labels = c("D", "E"), label_size = 8)


f1 <- plotremoval(hypvsdf, "HTR2C")+ labs(y = NULL) + labs(subtitle = " ")
f2 <- plotremoval(hypvsdm, "HTR2C") + labs(y = NULL, x = "")+ labs(subtitle = " ")

g1 <- plotreplacement(hypvsdf, "HTR2C") + labs(y = NULL)+ labs(subtitle = " ")
g2 <- plotreplacement(hypvsdm, "HTR2C") + labs(y = NULL, x = "")+ labs(subtitle = " ")


#h <- makebargraphv3(DEGremove, "hypothalamus", "No. of DEGs\n (-)decreased  increased (+)", comparisonlablelsremoval, comparisonlevelsremoval) +
  labs(x = "Offspring removal versus temporal control")

i <- makebargraphv3(DEGreplace, "hypothalamus", NULL, comparisonlablelssreplace, comparisonlevelsreplace) +
  labs(x = "Offspring removal versus temporal or external control")

fghi <- plot_grid(f1,f2,g1,g2,i, nrow = 1, 
          labels = c("F", " ", "G", " ", "H", "I"), label_size = 8)


fig2 <- plot_grid(ab,c,de, fghi, nrow = 4, 
          labels = c(" ", "C"), label_size = 8)




pdf(file="figures/fig2-1.pdf", width=7, height=7)
plot(fig2)
dev.off()

png("figures/fig2-1.png", width = 7, height = 7, 
    units = 'in', res = 300)
plot(fig2) 
dev.off()


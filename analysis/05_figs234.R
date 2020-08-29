# Fig 2, 3, and 4

library(tidyverse)
library(cowplot)
library(ggsignif)

source("R/themes.R")
source("R/functions.R")

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



makefigs234 <- function(tissue, fvsdpath, mvsdpath, 
                        candidategene, comp1, comp2, figno){
  
  dff <- wranglevsds(fvsdpath)
  dfm <- wranglevsds(mvsdpath)
  
  print(head(dff))
  print(head(dfm))

  a <- plotcandidatechar(dff, candidategene) + labs(subtitle = tissue)
  b <- plotremoval(dff, candidategene) + labs(y = NULL) + labs(subtitle = " ")
  c <- plotreplacement(dff, candidategene) + labs(y = NULL) + labs(subtitle = " ")
  d <- plot.volcano(tissue, "female",  comp1) + facet_wrap(~sex)  + labs(subtitle = " ")
  
  abc <- plot_grid(a,b,c, nrow = 1, align = "h")
  abcd <- plot_grid(abc,d, nrow = 1, rel_widths = c(3,1))
  
  e <- makebargraphv4(DEGcontrolreplace, tissue, 
                      "No. of DEGs\n(-) decreased  increased (+)", 
                      comparisonlabelscontrolreaplce, comparisonlevelscontrolreplace,
                      "female") +
    labs(x = "Control versus all other reproductive and parental stages") 
  
  f <- makebargraphv4(DEGbldg, tissue, 
                      "No. of DEGs\n (-)decreased  increased (+)", 
                      comparisonlevelsbldg, comparisonlevelsbldg,
                      "female") +
    labs(x = "Nest-building versus all other parental stages")
  
  g <- makebargraphv4(DEGchar, tissue, 
                      NULL,  
                      comparisonlevelscharnobldg, comparisonlevelscharnobldg,
                      "female") +
    labs(x = "Comparison of sequential parental stages")
  
  efg <- plot_grid(e,f,g, nrow = 1)
  
  
  h <- plotcandidatechar(dfm, candidategene) + labs(subtitle = tissue)
  i <- plotremoval(dfm, candidategene) + labs(y = NULL) + labs(subtitle = " ")
  j <- plotreplacement(dfm, candidategene) + labs(y = NULL) + labs(subtitle = " ")
  l <- plot.volcano(tissue, "male",  comp2) + facet_wrap(~sex)  + labs(subtitle = " ")
  
  hij <- plot_grid(h,i,j, nrow = 1, align = "h")
  hijl <- plot_grid(hij,l, nrow = 1, rel_widths = c(3,1))
  
  l <- makebargraphv4(DEGcontrolreplace, tissue, 
                      "No. of DEGs\n(-) decreased  increased (+)", 
                      comparisonlabelscontrolreaplce, comparisonlevelscontrolreplace,
                      "male") +
    labs(x = "Control versus all other reproductive and parental stages") 
  
  m <- makebargraphv4(DEGbldg, tissue, 
                      "No. of DEGs\n (-)decreased  increased (+)", 
                      comparisonlevelsbldg, comparisonlevelsbldg,
                      "male") +
    labs(x = "Nest-building versus all other parental stages")
  
  n <- makebargraphv4(DEGchar, tissue, 
                      NULL,  
                      comparisonlevelscharnobldg, comparisonlevelscharnobldg,
                      "male") +
    labs(x = "Comparison of sequential parental stages")
  
  lmn <- plot_grid(l,m,n, nrow = 1)
  
  fig <- plot_grid(abcd,efg,hijl, lmn, nrow = 4)

pdffilename <- paste("figures/", figno, "-1.pdf", sep = "")
pngfilename <- paste("figures/", figno, "-1.png", sep = "")

print(fig)

pdf(file = pdffilename, width=7, height=7)
plot(fig)
dev.off()

png(file = pngfilename, width = 7, height = 7, 
    units = 'in', res = 300)
plot(fig) 
dev.off()

}

makefigs234("hypothalamus",
            "results/03_hypvsdf.csv", 
            "results/03_hypvsdm.csv",
            "HTR2C",
            "control_hatch",
            "control_hatch",
            "fig2")

makefigs234("pituitary",
            "results/03_pitvsdf.csv", 
            "results/03_pitvsdm.csv",
            "PRL",
            "inc.d9_inc.d17",
            "inc.d9_inc.d17",
            "fig3")

makefigs234("gonads",
            "results/03_gonvsdf.csv", 
            "results/03_gonvsdm.csv",
            "ESR1",
            "control_lay",
            "control_lay",
            "fig4")
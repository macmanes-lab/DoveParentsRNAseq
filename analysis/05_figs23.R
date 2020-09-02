# Fig 2, 3, and 4

library(tidyverse)
library(cowplot)
library(ggsignif)

source("R/themes.R")
source("R/functions.R")

###### variance statbilized data

hypf <- wranglevsds("results/03_hypvsdf.csv")
hypm <- wranglevsds("results/03_hypvsdm.csv")
pitf <- wranglevsds("results/03_pitvsdf.csv")
pitm <- wranglevsds("results/03_pitvsdm.csv")
gonf <- wranglevsds("results/03_gonvsdf.csv")
gonm <- wranglevsds("results/03_gonvsdm.csv")

## figure 2 candidate genes

makefig2 <- function(tissue, dff, dfm, candidategene, comp1, comp2, label1){
  
  fsubtitle = paste("Female", tissue, sep = " ")
  msubtitle = paste("Male", tissue, sep = " ")

  a <- plotcandidatechar(dff, candidategene) + labs(subtitle = fsubtitle) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  b <- plotremoval(dff, candidategene) + labs(y = NULL) + labs(subtitle = " ")+
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  c <- plotreplacement(dff, candidategene) + labs(y = NULL) + labs(subtitle = " ")+
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  d <- plot.volcano(tissue, "female",  comp1) + labs(subtitle = " ") +
    theme(axis.title.x = element_blank(),
          legend.position = "none")
  d2 <- plot.volcano(tissue, "female",  comp2) + labs(subtitle = " ")  +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
  
  abc <- plot_grid(a,b,c, nrow = 1)
  abcd <- plot_grid(abc,d,d2, nrow = 1, rel_widths = c(3,1,0.9))
  

  h <- plotcandidatechar(dfm, candidategene) + labs(subtitle = msubtitle)
  i <- plotremoval(dfm, candidategene) + labs(y = NULL) + labs(subtitle = " ")
  j <- plotreplacement(dfm, candidategene) + labs(y = NULL) + labs(subtitle = " ")
  l <- plot.volcano(tissue, "male",  comp1)  + labs(subtitle = " ") 
  l2 <- plot.volcano(tissue, "male",  comp2) + labs(subtitle = " ") +
    theme(axis.title.y = element_blank())
  
  hij <- plot_grid(h,i,j, nrow = 1)
  hijl <- plot_grid(hij,l,l2, nrow = 1, rel_widths = c(3,1,1))
  
  fig <- plot_grid(abcd, hijl, nrow = 2,
                   rel_heights = c(1,1.2),
                   labels = c(label1), label_size = 8)
  return(fig)
}

ab <- makefig2("hypothalamus", hypf, hypm, "HTR2C", "control_hatch", "hatch_m.n2", "A")
cd <- makefig2("pituitary", pitf, pitm, "PRL", "inc.d9_inc.d17", "hatch_early", "B")
ef <- makefig2("gonads", gonf, gonm, "ESR1", "control_lay", "lay_inc.d3", "C")

fig2 <- plot_grid(ab,cd,ef, nrow = 3)

pdf(file = "figures/fig2-1.pdf", width=7, height=7)
plot(fig2)
dev.off()

png(file = "figures/fig2-1.png", width = 7, height = 7, 
    units = 'in', res = 300)
plot(fig2) 
dev.off()

###### DEGs

allDEG <- read_csv("results/03_allDEG.csv") %>%
  mutate(sex = factor(sex, levels = sexlevels),
         tissue = factor(tissue, levels = tissuelevels)) 
head(allDEG)


filterDEGs <- function(whichlevels){  
  df <- allDEG %>% 
    filter(comparison %in% whichlevels) %>%
    droplevels() %>%
    mutate(comparison = factor(comparison, levels = whichlevels)) %>%
    mutate(updown = ifelse(lfc > 0, 1, -1)) %>%
    group_by(sex, tissue, comparison, direction, updown, .drop=FALSE) %>%
    summarise(n = n()) %>%
    mutate(n = n*updown ) %>%
    select(-updown) %>%
    dplyr::mutate(n = replace_na(n, 0))
  print(head(df))
  return(df)
}

DEGcharmanip <- filterDEGs(levelscontrolcharmanip)
DEGbldg <- filterDEGs(levelsbldgcharmanip)
DEGsequential <- filterDEGs(levelssequential)
DEGmanip <- filterDEGs(levelsmanip)

### fig 3 DEGs

makefig3 <- function(tissue, label1){
  
  fsubtitle = paste("Female", tissue, sep = " ")
  msubtitle = paste("Male", tissue, sep = " ")
  
  e <- makebargraphv4(DEGcharmanip, tissue, 
                      "No. of DEGs", 
                      levelscontrolcharmanip, labelscontrolcharmanip, 
                      "female") +
    labs(x =  NULL, subtitle = fsubtitle) +
    theme(axis.text.x = element_blank())
  
  f <- makebargraphv4(DEGbldg, tissue, 
                      NULL, 
                      levelsbldgcharmanip, labelsbldgcharmanip, 
                      "female") +
    labs(x = NULL,
         subtitle = " ") +
    theme(axis.text.x = element_blank())
  
  g <- makebargraphv4(DEGsequential, tissue, 
                      NULL,  
                      levelssequential, labelsssequential,
                      "female") +
    labs(x = NULL,
         subtitle = " ") +
    theme(axis.text.x = element_blank())
  
  g2 <- makebargraphv4(DEGmanip, tissue, 
                       NULL,  
                       levelsmanip, labelssmanip,
                       "female") +
    labs(x = NULL,
         subtitle = " ") +
    theme(axis.text.x = element_blank())
  
  efg <- plot_grid(e,f,g,g2, nrow = 1, rel_widths = c(2,2,1,1))

  l <- makebargraphv4(DEGcharmanip, tissue, 
                      "No. of DEGs", 
                      levelscontrolcharmanip, labelscontrolcharmanip, 
                      "male") +
    labs(x = "Relative to... non-breeing controls,",
         subtitle = msubtitle) 
  
  
  m <- makebargraphv4(DEGbldg, tissue, 
                      NULL,
                      levelsbldgcharmanip, labelsbldgcharmanip, 
                      "male") +
    labs(x = "nest-building controls",
         subtitle = " ")
  
  n <- makebargraphv4(DEGsequential, tissue, 
                      NULL,  
                      levelssequential, labelsssequential,
                      "male") +
    labs(x = "previous stages,", 
         subtitle = " ")
  
  n2 <- makebargraphv4(DEGmanip, tissue, 
                       NULL,  
                       levelsmanip, labelssmanip,
                       "male") +
    labs(x = "other manipulations",
         subtitle = " ")
  
  lmn <- plot_grid(l,m,n,n2, nrow = 1, rel_widths = c(2,2,1,1))
  
  fig <- plot_grid(efg,lmn, nrow = 2, rel_heights = c(1,1.2),
                   labels = label1, label_size = 8)
  return(fig)
  
  
  
}

ab <- makefig3("hypothalamus", "A")
cd <- makefig3("pituitary", "B")
ef <- makefig3("gonads", "C")

fig3 <- plot_grid(ab,cd,ef, nrow = 3)

png(file = "figures/fig3-1.png", width = 7, height = 7, 
    units = 'in', res = 300)
plot(fig3) 
dev.off()

pdf(file = "figures/fig3-1.pdf", width=7, height=7)
plot(fig3)
dev.off()

# Fig 2 and 3

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

allDEG <- read_csv("results/04_allDEG.csv")

## figure 2 candidate genes

makefig2 <- function(label1, label2, tissue, dff, dfm, 
                     candidategene, comp1, comp2, 
                     treatment1, treatment2, treatment3, treatment4){
  
  fsubtitle = paste("Female", tissue, sep = " ")
  msubtitle = paste("Male", tissue, sep = " ")

  a <- plotcandidatechar(dff, candidategene, treatment1, treatment2) + labs(subtitle = fsubtitle) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  b <- plotcandidatemanip(dff, candidategene, treatment3, treatment4) + labs(y = NULL) + labs(subtitle = " ") +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  
  c <- plot.volcano(tissue, "female",  comp1) + labs(subtitle = " ") +
    theme(axis.title.x = element_blank(),
          legend.position = "none")
  d <- plot.volcano(tissue, "female",  comp2) + labs(subtitle = " ")  +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none")
  
  abcd <- plot_grid(a,b,c,d,nrow = 1, rel_widths = c(9,12,7.35,7),
                    labels = c(label1, "", label2), label_size = 8)
  

  e <- plotcandidatechar(dfm, candidategene, treatment1, treatment2) + labs(subtitle = msubtitle)
  f <- plotcandidatemanip(dfm, candidategene, treatment3, treatment4) + labs(y = NULL) + labs(subtitle = " ")
  
  g <- plot.volcano(tissue, "male",  comp1)  + labs(subtitle = " ") 
  h <- plot.volcano(tissue, "male",  comp2) + labs(subtitle = " ") +
    theme(axis.title.y = element_blank())
  
  
  efgh <- plot_grid(e,f,g,h,nrow = 1, rel_widths = c(9,12,7.35,7))
  
  fig <- plot_grid(abcd, efgh, nrow = 2,
                   rel_heights = c(1,1.2))
  return(fig)
}

ab <- makefig2("A", "D","hypothalamus", hypf, hypm, "VIPR1", 
               "control_hatch", "hatch_m.n2", 
               "control", "hatch", "hatch", "m.n2")
cd <- makefig2("B", "E",  "pituitary", pitf, pitm, "PRL", 
               "inc.d9_inc.d17", "hatch_early", 
               "inc.d9", "inc.d17",  "hatch", "early")
ef <- makefig2("C", "F", "gonads", gonf, gonm, "ESR2",
               "control_lay", "inc.d17_prolong",
               "control", "lay", "inc.d17", "prolong")

fig2 <- plot_grid(ab,cd,ef, nrow = 3)


png(file = "figures/fig2-1.png", width = 7, height = 7, 
    units = 'in', res = 300)
plot(fig2) 
dev.off()

pdf(file = "figures/fig2-1.pdf", width=7, height=7)
plot(fig2)
dev.off()

###### DEGs


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

DEGcontrol <- filterDEGs(levelscontrolcharmanip)
DEGbldg <- filterDEGs(levelsbldgcharmanip)
DEGsequential <- filterDEGs(levelssequential)

DEGinc9 <- filterDEGs(levelsinc9)
DEGinc17 <- filterDEGs(levelsinc17)
DEGhatch <- filterDEGs(levelshatch)

### fig 3 DEGs

makefig3 <- function(tissue, label1){
  
  fsubtitle = paste("Female", tissue, sep = " ")
  msubtitle = paste("Male", tissue, sep = " ")
  
  a <- makebargraphv4(DEGcontrol, tissue, "No. of DEGs", 
                      levelscontrolcharmanip, labelscontrolcharmanip, "female",
                      -4000, 4000) +
    labs(x =  NULL, subtitle = fsubtitle) +
    theme(axis.text.x = element_blank()) 
  
  b <- makebargraphv4(DEGbldg, tissue, NULL, 
                      levelsbldgcharmanip, labelsbldgcharmanip, "female",
                      -1500, 1500) +
    labs(x = NULL, subtitle = " ") +
    theme(axis.text.x = element_blank()) 
  
  c <- makebargraphv4(DEGsequential, tissue, NULL,  
                      levelssequential, labelsssequential,"female",
                      -1500, 1500) +
    labs(x = NULL, subtitle = " ") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  d <- makebargraphv4(DEGinc9, tissue, NULL,  
                      levelsinc9, labelsinc9, "female",
                      -1500, 1500) +
    labs(x = NULL, subtitle = " ") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  e <- makebargraphv4(DEGinc17, tissue, NULL,  
                      levelsinc17, labelsinc17, "female",
                      -1500, 1500) +
    labs(x = NULL, subtitle = " ") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  f <- makebargraphv4(DEGhatch, tissue, NULL,  
                      levelshatch, labelshatch, "female",
                      -1500, 1500) +
    labs(x = NULL, subtitle = " ") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  abcdef <- plot_grid(a,b,c,d,e,f, nrow = 1, rel_widths = c(26,20,10,4,4,6))

  g <- makebargraphv4(DEGcontrol, tissue, "No. of DEGs", 
                      levelscontrolcharmanip, labelscontrolcharmanip, "male",
                      -4000, 4000) +
    labs(x = "Relative to... non-breeing controls,", subtitle = msubtitle) 
  
  h <- makebargraphv4(DEGbldg, tissue, NULL,
                      levelsbldgcharmanip, labelsbldgcharmanip, "male",
                      -1500, 1500) +
    labs(x = "nest-building controls...", subtitle = " ")
  
  i <- makebargraphv4(DEGsequential, tissue,  NULL,  
                      levelssequential, labelsssequential, "male",
                      -1500, 1500) +
    labs(x = "previous stage...",  subtitle = " ") +
    theme(axis.text.y = element_blank())
  
  j <- makebargraphv4(DEGinc9, tissue, NULL,  
                      levelsinc9, labelsinc9, "male",
                      -1500, 1500) +
    labs(x = "inc.d9...",  subtitle = " ") +
    theme(axis.text.y = element_blank()) 
  
  k<- makebargraphv4(DEGinc17, tissue, NULL,  
                      levelsinc17, labelsinc17, "male",
                     -1500, 1500) +
    labs(x = NULL, subtitle = " ")+
    labs(x = "inc.d17...",  subtitle = " ") +
    theme(axis.text.y = element_blank())
  
  
  l <- makebargraphv4(DEGhatch, tissue, NULL,  
                      levelshatch, labelshatch, "male",
                      -1500, 1500) +
    labs(x = "hatch",  subtitle = " ") +
    theme(axis.text.y = element_blank())
  
  ghijkl <- plot_grid(g,h,i,j,k,l, nrow = 1, rel_widths = c(26,20,10,4,4,6),
                      align = "h")
  
  fig <- plot_grid(abcdef,ghijkl, nrow = 2, rel_heights = c(1,1.2),
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

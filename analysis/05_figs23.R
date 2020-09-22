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

hyp <- rbind(hypf, hypm)
pit <- rbind(pitf, pitm)
gon <- rbind(gonf, gonm)

allDEG <- read_csv("results/04_allDEG.csv")


##  figure 2

a <- plotcandidatechar(hyp, "NPVF") 
d <- plotcandidatemanip(hyp, "NPVF") 
b <- plotcandidatechar(pit, "PRL")
e <- plotcandidatemanip(pit, "PRL")
c <- plotcandidatechar(gon, "ESR1")
f <- plotcandidatemanip(gon, "ESR1")

abc <- plot_grid(a,b,c, nrow = 1)
def <- plot_grid(d,e,f, nrow = 1)

plot_grid(abc,def, nrow  = 2)

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

max(DEGcontrol$n)
min(DEGcontrol$n)

max(DEGbldg$n)
max(DEGsequential$n)
max(DEGinc9$n)
max(DEGinc17$n)
max(DEGhatch$n)

min(DEGbldg$n)
min(DEGsequential$n)
min(DEGinc9$n)
min(DEGinc17$n)
min(DEGhatch$n)



### fig 3 DEGs

makefig3 <- function(tissue, label1){
  
  fsubtitle = paste("Female", tissue, sep = " ")
  msubtitle = paste("Male", tissue, sep = " ")
  
  a <- makebargraphv5(DEGcontrol, tissue, "No. of DEGs", 
                      levelscontrolcharmanip, labelscontrolcharmanip, "female",
                      -4000, 3000) +
    labs(x =  NULL, subtitle = fsubtitle) +
    theme(axis.text.x = element_blank()) 
  
  b <- makebargraphv4(DEGbldg, tissue, NULL, 
                      levelsbldgcharmanip, labelsbldgcharmanip, "female",
                      -1600, 1100) +
    labs(x = NULL, subtitle = " ") +
    theme(axis.text.x = element_blank()) 
  
  c <- makebargraphv4(DEGsequential, tissue, NULL,  
                      levelssequential, labelsssequential,"female",
                      -1600, 1100) +
    labs(x = NULL, subtitle = " ") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  d <- makebargraphv4(DEGinc9, tissue, NULL,  
                      levelsinc9, labelsinc9, "female",
                      -1600, 1100) +
    labs(x = NULL, subtitle = " ") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  e <- makebargraphv4(DEGinc17, tissue, NULL,  
                      levelsinc17, labelsinc17, "female",
                      -1600, 1100) +
    labs(x = NULL, subtitle = " ") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  f <- makebargraphv4(DEGhatch, tissue, NULL,  
                      levelshatch, labelshatch, "female",
                      -1600, 1100) +
    labs(x = NULL, subtitle = " ") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  abcdef <- plot_grid(a,b,c,d,e,f, nrow = 1, rel_widths = c(20,20,10,4,4,6),
                      labels = label1, label_size = 8)

  g <- makebargraphv5(DEGcontrol, tissue, "No. of DEGs", 
                      levelscontrolcharmanip, labelscontrolcharmanip, "male",
                      -4000, 3000) +
    labs(x = "Relative to... non-breeing controls,", subtitle = msubtitle) 
  
  h <- makebargraphv4(DEGbldg, tissue, NULL,
                      levelsbldgcharmanip, labelsbldgcharmanip, "male",
                      -1600, 1100) +
    labs(x = "nest-building controls...", subtitle = " ") 
  
  i <- makebargraphv4(DEGsequential, tissue,  NULL,  
                      levelssequential, labelsssequential, "male",
                      -1700, 1200) +
    labs(x = "previous stage...",  subtitle = " ") +
    theme(axis.text.y = element_blank())
  
  j <- makebargraphv4(DEGinc9, tissue, NULL,  
                      levelsinc9, labelsinc9, "male",
                      -1600, 1100) +
    labs(x = "inc.d9...",  subtitle = " ") +
    theme(axis.text.y = element_blank()) 
  
  k<- makebargraphv4(DEGinc17, tissue, NULL,  
                      levelsinc17, labelsinc17, "male",
                     -1700, 1200) +
    labs(x = NULL, subtitle = " ")+
    labs(x = "inc.d17...",  subtitle = " ") +
    theme(axis.text.y = element_blank())
  
  
  l <- makebargraphv4(DEGhatch, tissue, NULL,  
                      levelshatch, labelshatch, "male",
                      -1600, 1100) +
    labs(x = "and hatch.",  subtitle = " ") +
    theme(axis.text.y = element_blank())
  
  ghijkl <- plot_grid(g,h,i,j,k,l, nrow = 1, rel_widths = c(20,20,10,4,4,6),
                      align = "h")
  
  fig <- plot_grid(abcdef,ghijkl, nrow = 2, rel_heights = c(1,1.2))
  return(fig)
}

ab <- makefig3("hypothalamus", c("A1", "2", "3", "4", "5", "6"))
cd <- makefig3("pituitary", c("B"))
ef <- makefig3("gonads", c("C1"))

fig3 <- plot_grid(ab,cd,ef, nrow = 3)

png(file = "figures/fig3-1.png", width = 7, height = 7, 
    units = 'in', res = 300)
plot(fig3) 
dev.off()

pdf(file = "figures/fig3-1.pdf", width=7, height=7)
plot(fig3)
dev.off()

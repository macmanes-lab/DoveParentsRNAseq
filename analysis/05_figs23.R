# Fig 2 and 3

library(tidyverse)
library(cowplot)
library(ggsignif)
library(stringr)

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

##  figure 2 DEGs
fig2Acomps <- c("bldg_inc.d17", "inc.d9_inc.d17" )
fig2Dcomps <- c("inc.d17_m.inc.d17", "inc.d17_prolong" )

fig2Alabels <- c("inc.d9_inc.d17 vs.\ninc.d9", "inc.d9 vs.\ninc.d17" )
fig2Dlabels <- c("inc.d17vs\nm.inc.d17", "inc.d17vs\nprolong" )


fig2Adegs <- filterDEGs(fig2Acomps)
fig2Ddegs <- filterDEGs(fig2Dcomps)


a <- plotcandidatechar(hyp, "AVP") + labs(subtitle = "Hypothalamus", x = " ", title = "Characterization") 
f <- plotcandidatemanip(hyp, "AVP") + labs(subtitle = "Hypothalamus",  x = " ", title = "Manipulation")

b <- plotcandidatechar(pit, "PRL") + labs(subtitle = "Pituitary")
g <- plotcandidatemanip(pit, "PRL") + labs(subtitle = "Pituitary")

c <- plotcandidatechar(gon, "ESR1") + labs(subtitle = "Gonads", x = " ")
h <- plotcandidatemanip(gon, "ESR1") + labs(subtitle = "Gonads", x = " ")


d <- plot.volcano("pituitary", "bldg_inc.d17") + labs(subtitle = "Pituitary")
e <- plot.volcano("pituitary", "inc.d9_inc.d17") + labs(subtitle = " ", y = NULL) 

e2 <- makebargraphv4(fig2Adegs, "pituitary", "No. of DEGs", 
                     fig2Acomps, fig2Alabels)  + 
  labs(x = "Comparison", subtitle = "Pituitary", title = " ")


i <- plot.volcano("pituitary", "inc.d17_m.inc.d17")  + labs(subtitle = "Pituitary")
j <- plot.volcano("pituitary", "inc.d17_prolong") + labs(subtitle = " ", y = NULL)

j2 <- makebargraphv4(fig2Ddegs, "pituitary", "No. of DEGs", 
                     fig2Dcomps, fig2Dlabels) + 
  labs(x = "Comparison", subtitle = "Pituitary", title = " ")

abcde <- plot_grid(a,b,c, nrow = 1, 
                   labels = c("A", "B", "C" ), label_size = 8,
                   rel_widths = c(1,1,1))
fghij <- plot_grid(f,g,h,nrow = 1, 
                   labels = c("D", "E", "F" ), label_size = 8,
                   rel_widths = c(1,1,1))

fig2 <- plot_grid(abcde,fghij, nrow  = 2)
fig2

png(file = "figures/fig2-1.png", width = 7, height = 7, 
    units = 'in', res = 300)
plot(fig2) 
dev.off()

pdf(file = "figures/fig2-1.pdf", width=7, height=7)
plot(fig2)
dev.off()

### fig 3 DEGs

DEGbldg <- filterDEGs(levelsbldgchar)
DEGsequential <- filterDEGs(levelssequential)
DEGinc9 <- filterDEGs(levelsinc9)
DEGinc17 <- filterDEGs(levelsinc17)
DEGhatch <- filterDEGs(levelshatch)
DEGearly <- filterDEGs(levelsearly) %>%
  mutate(tissue = factor(tissue, levels = tissuelevels),
         sex = factor(sex , levels = sexlevels))





makefig3 <- function(tissue, label1){

  tissuetitle <- str_to_sentence(tissue, locale = "en")
  
  ylab <- paste(tissuetitle, "\nNo. of DEGs", sep = "")
  
  #ylab <- "Genes with increased expression"
  
  
  a <-  plot.volcano(tissue, "bldg_control") + 
    labs(subtitle = tissuetitle) 
  
  b <- makebargraphv5(DEGbldg, tissue, ylab, 
                      levelsbldgchar, labelsbldgchar) +
    labs(x = "Relative to nest-building controls,", subtitle = " ") 
  
  c <- makebargraphv5(DEGsequential, tissue, NULL,  
                      levelssequential, labelsssequential) +
    labs(x = "Relative to the previous stage,", subtitle = " ") 
  
  d <- makebargraphv5(DEGinc9, tissue, NULL,  
                      levelsinc9, labelsinc9) +
    labs(x = "... to inc.d9,", subtitle = " ") 
  
  
  e <- makebargraphv5(DEGinc17, tissue, NULL,  
                      levelsinc17, labelsinc17) +
    labs(x = "... to inc.d17,", subtitle = " ") 
  
  f <- makebargraphv5(DEGhatch, tissue, NULL,  
                      levelshatch, labelshatch) +
    labs(x = "... to hatch,", subtitle = " ") 
 
   g <- makebargraphv5(DEGearly, tissue, NULL,  
                      levelsearly, labelsearly) +
    labs(x = "Relative to early.", subtitle = " ") 
  
  fig <- plot_grid(b,c,d, e,f, nrow = 1, rel_widths = c(7,5.5,2.5,2.5,3.5),
                      labels = label1, label_size = 8, align = "h")
  
  return(fig)
}

ab <- makefig3("hypothalamus", c("A"))
cd <- makefig3("pituitary", c("B"))
ef <- makefig3("gonads", c("C"))

fig3 <- plot_grid(ab,cd,ef, nrow = 3)
fig3


png(file = "figures/fig3-1.png", width = 7, height = 7, 
    units = 'in', res = 300)
plot(fig3) 
dev.off()

pdf(file = "figures/fig3-1.pdf", width=7, height=7)
plot(fig3)
dev.off()




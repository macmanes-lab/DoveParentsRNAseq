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


###### DEGs

allDEG <- read_csv("results/03_allDEG.csv") 


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


makefigs234 <- function(tissue, dff, dfm, 
                        candidategene, comp1, comp2, figno){
  
  fsubtitle = paste("Female", tissue, sep = " ")
  msubtitle = paste("Male", tissue, sep = " ")

  a <- plotcandidatechar(dff, candidategene) + labs(subtitle = fsubtitle)
  b <- plotremoval(dff, candidategene) + labs(y = NULL) + labs(subtitle = " ")
  c <- plotreplacement(dff, candidategene) + labs(y = NULL) + labs(subtitle = " ")
  d <- plot.volcano(tissue, "female",  comp1) + labs(subtitle = " ")
  d2 <- plot.volcano(tissue, "female",  comp2) + labs(subtitle = " ") 
  
  abc <- plot_grid(a,b,c, nrow = 1)
  abcd <- plot_grid(abc,d,d2, nrow = 1, rel_widths = c(3,1,1), 
                    labels = c("A", "B"), label_size = 8)
  
  e <- makebargraphv4(DEGcharmanip, tissue, 
                      "No. of DEGs\n(-) decreased  increased (+)", 
                      levelscontrolcharmanip, labelscontrolcharmanip, 
                      "female") +
    labs(x = "Control versus all other reproductive and parental stages",
         subtitle = " ")
  
  f <- makebargraphv4(DEGbldg, tissue, 
                      "No. of DEGs\n (-)decreased  increased (+)", 
                      levelsbldgcharmanip, labelsbldgcharmanip, 
                      "female") +
    labs(x = "Nest-building versus all other parental stages",
         subtitle = " ")
  
  g <- makebargraphv4(DEGsequential, tissue, 
                      NULL,  
                      levelssequential, labelsssequential,
                      "female") +
    labs(x = "Comparison of sequential parental stages",
         subtitle = " ")
  
  g2 <- makebargraphv4(DEGmanip, tissue, 
                      NULL,  
                      levelsmanip, labelssmanip,
                      "female") +
    labs(x = "Comparison of manipulations",
         subtitle = " ")
  
  efg <- plot_grid(e,f,g,g2, nrow = 1, rel_widths = c(2,2,1,1),
                   labels = c("C"), label_size = 8)
  
  
  h <- plotcandidatechar(dfm, candidategene) + labs(subtitle = msubtitle)
  i <- plotremoval(dfm, candidategene) + labs(y = NULL) + labs(subtitle = " ")
  j <- plotreplacement(dfm, candidategene) + labs(y = NULL) + labs(subtitle = " ")
  l <- plot.volcano(tissue, "male",  comp1)  + labs(subtitle = " ")
  l2 <- plot.volcano(tissue, "male",  comp2) + labs(subtitle = " ") 
  
  hij <- plot_grid(h,i,j, nrow = 1)
  hijl <- plot_grid(hij,l,l2, nrow = 1, rel_widths = c(3,1,1),
                    labels = c("D", "E"), label_size = 8)
  
  l <- makebargraphv4(DEGcharmanip, tissue, 
                      "No. of DEGs\n(-) decreased  increased (+)", 
                      levelscontrolcharmanip, labelscontrolcharmanip, 
                      "male") +
    labs(x = "Control versus all other reproductive and parental stages",
         subtitle = " ") 

  
  m <- makebargraphv4(DEGbldg, tissue, 
                      "No. of DEGs\n (-)decreased  increased (+)", 
                      levelsbldgcharmanip, labelsbldgcharmanip, 
                      "male") +
    labs(x = "Nest-building versus all other parental stages",
         subtitle = " ")
  
  n <- makebargraphv4(DEGsequential, tissue, 
                      NULL,  
                      levelssequential, labelsssequential,
                      "male") +
    labs(x = "Comparison of sequential parental stages", 
         subtitle = " ")
  
  n2 <- makebargraphv4(DEGmanip, tissue, 
                       NULL,  
                       levelsmanip, labelssmanip,
                       "male") +
    labs(x = "Comparison of manipulations",
         subtitle = " ")
  
  lmn <- plot_grid(l,m,n,n2, nrow = 1, rel_widths = c(2,2,1,1),
                   labels = c("F"), label_size = 8)
  
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



makefigs234("hypothalamus", hypf,hypm,
            "HTR2C",
            "control_hatch", "hatch_early",
            "fig2")

makefigs234("pituitary", pitf, pitm,
            "PRL",
            "inc.d9_inc.d17", "hatch_early",
            "fig3")

makefigs234("gonads", gonf, gonm,
            "ESR1",
            "control_lay", "hatch_early",
            "fig4")
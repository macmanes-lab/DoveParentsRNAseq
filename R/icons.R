# individual icons
library(tidyverse)
library(png)
library(grid)

birds <- data.frame(iconpath=list.files("../figures/images/icons", full.names = T), stringsAsFactors = F)
birds$icons <- sapply(strsplit(as.character(birds$iconpath),'../figures/images/icons/DoveParentsRNAseq_'), "[", 2)

iconstreatment <- data.frame(treatment = c("control", "bldg", "lay", 
                                          "inc.d3", "inc.d9", "inc.d17", 
                                          "hatch", "n5", "n9",
                                          "m.inc.d3" ,  "m.inc.d9" , "m.inc.d17" , "m.n2",
                                          "m.inc.d8" , "prolong" , "extend") ,
                             icons = c("control.png", "bldg.png", "lay.png",
                                       "incubation.png", "incubation.png", "incubation.png",
                                       "hatch.png", "chicklings.png", "chicklings.png", 
                                       "bldg.png", "bldg.png", "bldg.png", "bldg.png",
                                       "hatch.png", "incubation.png","hatch.png"))
iconstreatment$music <- "https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcS9zgaht3nmVL4moYqz0uKBhYefpTgYrX69_qyL8kmznImDRk-KRw&s"

birds <- left_join(iconstreatment, birds, by = "icons")

control <- png::readPNG(birds$iconpath[1])
control <-  grid::rasterGrob(control, interpolate=TRUE)

bldg <- png::readPNG(birds$iconpath[2])
bldg <-  grid::rasterGrob(bldg, interpolate=TRUE)

lay <- png::readPNG(birds$iconpath[3])
lay <-  grid::rasterGrob(lay, interpolate=TRUE)

inc <- png::readPNG(birds$iconpath[4])
inc <-  grid::rasterGrob(inc, interpolate=TRUE)

hatch <- png::readPNG(birds$iconpath[7])
hatch <-  grid::rasterGrob(hatch, interpolate=TRUE)

nestling <- png::readPNG(birds$iconpath[8])
nestling <-  grid::rasterGrob(nestling, interpolate=TRUE)

# all characterization, for volcano plots
characterization <- png::readPNG("../figures/images/DoveParentsRNAseq_charicons.png")
characterization <- ggdraw() +  draw_image(characterization, scale = 1)
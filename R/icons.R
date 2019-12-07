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
                                       "removal-egg.png", "removal-egg.png", "removal-egg.png", "removal-chick.png",
                                       "manip-hatch.png", "manip-inc.png","manip-hatch.png"))
iconstreatment$music <- "https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcS9zgaht3nmVL4moYqz0uKBhYefpTgYrX69_qyL8kmznImDRk-KRw&s"

birds <- left_join(iconstreatment, birds, by = "icons")

control <-  grid::rasterGrob(png::readPNG(birds$iconpath[1]), interpolate=TRUE)
bldg <-  grid::rasterGrob(png::readPNG(birds$iconpath[2]), interpolate=TRUE)
lay <-  grid::rasterGrob(png::readPNG(birds$iconpath[3]), interpolate=TRUE)
inc <-  grid::rasterGrob(png::readPNG(birds$iconpath[4]), interpolate=TRUE)
hatch <-  grid::rasterGrob(png::readPNG(birds$iconpath[7]), interpolate=TRUE)
nestling <-  grid::rasterGrob(png::readPNG(birds$iconpath[8]), interpolate=TRUE)
removeegg <-  grid::rasterGrob(png::readPNG(birds$iconpath[10]), interpolate=TRUE)
removechick <-  grid::rasterGrob(png::readPNG(birds$iconpath[13]), interpolate=TRUE)
maniphatch <- grid::rasterGrob(png::readPNG(birds$iconpath[14]), interpolate=TRUE)
manipinc <- grid::rasterGrob(png::readPNG(birds$iconpath[15]), interpolate=TRUE)

# experimental design
characterization <- ggdraw() + draw_image(png::readPNG("../figures/images/DoveParentsRNAseq_timeline-char.png"), scale = 1)
removal <- ggdraw() + draw_image(png::readPNG("../figures/images/DoveParentsRNAseq_timeline-remove.png"), scale = 1)
timing <- ggdraw() + draw_image(png::readPNG("../figures/images/DoveParentsRNAseq_timeline-timing.png"), scale = 1)

removalonly <-  ggdraw() + draw_image(png::readPNG("../figures/images/DoveParentsRNAseq_volcanos-removal.png"), scale = 1)
timingonly <- ggdraw() + draw_image(png::readPNG("../figures/images/DoveParentsRNAseq_volcanos-timing.png"), scale = 1)

# individulal comparisons
early <-  ggdraw() + draw_image(png::readPNG("../figures/images/DoveParentsRNAseq_timeline-early.png"), scale = 1)
prolong <-  ggdraw() + draw_image(png::readPNG("../figures/images/DoveParentsRNAseq_timeline-prolong.png"), scale = 1)
extend <-  ggdraw() + draw_image(png::readPNG("../figures/images/DoveParentsRNAseq_timeline-delay.png"), scale = 1)



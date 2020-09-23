library(tidyverse)
library(cowplot)
library(Rtsne)

source("R/themes.R")
source("R/functions.R")

#  Experimental design, tSNE analysis figure

countData <- read_csv("results/00_counts.csv") %>%
  column_to_rownames(var = "X1")
countData <- as.data.frame(t(countData))
head(countData[1:3])

colData <- read_csv("metadata/00_colData.csv") %>%
  mutate(treatment = factor(treatment, levels = alllevels),
         tissue = factor(tissue, levels = tissuelevels)) %>% 
  column_to_rownames(var = "X1") 
head(colData)

# check ready for analysis
# row.names(countData) == row.names(colData)
head(row.names(countData) == row.names(colData))

## tsne

# prep for tsne useing count data and the custom `subsetmaketsne` function
hyptsne <- subsetmaketsne("hypothalamus", alllevelswithmanip, sexlevels)
pittsne <- subsetmaketsne("pituitary", alllevelswithmanip, sexlevels)
gontsne <- subsetmaketsne("gonads", alllevelswithmanip, sexlevels)

hyptsnechar <- subsetmaketsne("hypothalamus", charlevels, sexlevels)
pittsnechar <- subsetmaketsne("pituitary", charlevels, sexlevels)
gontsnechar <- subsetmaketsne("gonads", charlevels, sexlevels)


## Figure 1


ab <- png::readPNG("figures/images/fig_fig1a.png")
ab <- ggdraw() +  draw_image(ab, scale = 1)


f <- plottsne(hyptsnechar, hyptsnechar$treatment, allcolors) + 
  labs(subtitle = "hypothalamus", y = "Characterization\ntSNE2") 
g <- plottsne(pittsnechar, pittsnechar$treatment, allcolors ) + 
  labs(subtitle = "pituitary") 
theme(strip.text = element_blank())
h <- plottsne(gontsnechar, gontsnechar$treatment, allcolors ) + 
  labs(subtitle = "gonads") 

fgh <- plot_grid(f,g,h, ncol = 3, 
                 label_size = 8,
                 labels = c("C", "D", "E"),
                 rel_widths = c(1.2,1,1))



c <- plottsne(hyptsne, hyptsne$treatment, allcolors) + 
  labs(subtitle = "hypothalamus", y = "Manipulation\ntSNE2") 
d <- plottsne(pittsne, pittsne$treatment, allcolors ) + 
  labs(subtitle = "pituitary") 
theme(strip.text = element_blank())
e <- plottsne(gontsne, gontsne$treatment, allcolors ) + 
  labs(subtitle = "gonads") 

cde <- plot_grid(c,d,e, ncol = 3, 
                 labels = c("F", "G", "H"), 
                label_size = 8,
                rel_widths = c(1.2,1,1))



fig1 <- plot_grid(ab, fgh, cde, nrow = 3, rel_heights = c(2,1,1))
fig1

pdf(file="figures/fig1-1.pdf", width=7, height=6)
plot(fig1)
dev.off()

png("figures/fig1-1.png", width = 7, height = 6, 
    units = 'in', res = 300)
plot(fig1) 
dev.off()


sessionInfo()


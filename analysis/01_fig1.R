library(tidyverse)
library(cowplot)
library(Rtsne)


source("R/themes.R")
source("R/functions.R")


#  Experimental design, tSNE analysis, and bar chars


## import limma counts and sample info


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

# prep for tsne useing count data from limma and the custom `subsetmaketsne` function
hyptsne <- subsetmaketsne("hypothalamus", alllevels, sexlevels)
pittsne <- subsetmaketsne("pituitary", alllevels, sexlevels)
gontsne <- subsetmaketsne("gonads", alllevels, sexlevels)

## Figure 1

e1 <- plottsne(hyptsne, hyptsne$treatment, allcolors) + 
  labs(subtitle = "hypothalamus") 
e2 <- plottsne(pittsne, pittsne$treatment, allcolors ) + 
  labs(subtitle = "pituitary") 
  theme(strip.text = element_blank())
e3 <- plottsne(gontsne, gontsne$treatment, allcolors ) + 
  labs(subtitle = "gonads") 

e <- plot_grid(e1,e2,e3, ncol = 3, labels = c("D"), label_size = 8)

a <- png::readPNG("figures/images/fig_fig1a.png")
a <- ggdraw() +  draw_image(a, scale = 1)


fig1 <- plot_grid(a, e, nrow = 2, rel_heights = c(2.4,1.6))
#fig1

pdf(file="figures/fig1-1.pdf", width=7, height=4.5)
plot(fig1)
dev.off()

png("figures/fig1-1.png", width = 7, height = 4.5, 
    units = 'in', res = 300)
plot(fig1) 
dev.off()


sessionInfo()


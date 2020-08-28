# Fig 2

library(tidyverse)
library(cowplot)
library(ggsignif)

source("../R/themes.R")
source("../R/functions.R")
source("../R/genelists.R")




g <- plotcandidatechar(hypvsdf, "HTR2C") + labs(subtitle = "Hypothalamus")
h <- plotcandidatechar(hypvsdm, "HTR2C") + labs(y = NULL, x = "") + labs(subtitle = " ")
i <- plotremoval(hypvsdf, "HTR2C")+ labs(y = NULL) + labs(subtitle = " ")
j <- plotremoval(hypvsdm, "HTR2C") + labs(y = NULL, x = "")+ labs(subtitle = " ")
k <- plotreplacement(hypvsdf, "HTR2C") + labs(y = NULL)+ labs(subtitle = " ")
l <- plotreplacement(hypvsdm, "HTR2C") + labs(y = NULL, x = "")+ labs(subtitle = " ")

a <- plot.volcano("hypothalamus", sexlevels,  "control_bldg") + 
  facet_wrap(~sex) 

b <- makebargraphv3(DEGcontrolreplace, "hypothalamus","No. of DEGs\n(-) decreased  increased (+)", comparisonlabelscontrolreaplce) +
  labs(x = "Control versus all other reproductive and parental stages") +
  geom_rect(mapping=aes(xmin=0.5, xmax=1.5, ymin=-1000, ymax = 350, fill = F), color="black", alpha=0.5) +
  annotate("text", x = 1, y = -1075, label = "D", size = 2.5)   

c <- makebargraphv3(DEGbldg, "hypothalamus", "No. of DEGs\n (-)decreased  increased (+)", comparisonlabelsbldg) +
  labs(x = "Nest-building versus all other parental stages")
d <- makebargraphv3(DEGchar, "hypothalamus", NULL,  comparisonlabelscharnobldg) +
  labs(x = "Comparison of sequential parental stages")

e <- makebargraphv3(DEGremove, "hypothalamus", "No. of DEGs\n (-)decreased  increased (+)", comparisonlevelsremoval) +
  labs(x = "Offspring removal versus temporal control")
f <- makebargraphv3(DEGreplace, "hypothalamus", NULL,  comparisonlevelsreplace) +
  labs(x = "Offspring removal versus temporal or external control")


m <- png::readPNG("../figures/venn-eggs-hyp.png")
m <- ggdraw() +  draw_image(m, scale = 1)

n <- png::readPNG("../figures/venn-chicks-hyp.png")
n <- ggdraw() +  draw_image(n, scale = 1)

ghi <- plot_grid(g,h,k,l,nrow = 1,
                labels = c("A", "", "B"), label_size = 8, hjust = 0,
                rel_widths = c(9,9,6,6))

ab <- plot_grid(a,b,rel_widths = c(1,3), 
                labels = c("C", "D"), label_size = 8, hjust = 0)

cd <- plot_grid(c,d,rel_widths = c(1.1,1), align = "h",
                labels = c("E", "F"), label_size = 8, hjust = 0)

ef <- plot_grid(f,m,n, rel_widths = c(2,1,1), nrow =1,
                labels = c("G", "H", "I"), label_size = 8, hjust = 0)


fig2 <- plot_grid(ghi,ab,cd,ef,  ncol = 1)
fig2


a <- plot.volcano("pituitary", sexlevels,  "control_bldg") + 
  facet_wrap(~sex) 

b <- makebargraphv3(DEGcontrolreplace, "pituitary","No. of DEGs\n(-) decreased  increased (+)", comparisonlabelscontrolreaplce) +
  labs(x = "Control versus all other reproductive and parental stages") 

c <- makebargraphv3(DEGbldg, "pituitary", "No. of DEGs\n (-)decreased  increased (+)", comparisonlabelsbldg) +
  labs(x = "Nest-building versus all other parental stages")
d <- makebargraphv3(DEGchar, "pituitary", NULL,  comparisonlabelscharnobldg) +
  labs(x = "Comparison of sequential parental stages")

e <- makebargraphv3(DEGremove, "pituitary", "No. of DEGs\n (-)decreased  increased (+)", comparisonlevelsremoval) +
  labs(x = "Offspring removal versus temporal control")
f <- makebargraphv3(DEGreplace, "pituitary", NULL,  comparisonlablelssreplace) +
  labs(x = "Offspring removal versus temporal or external control")

g <- plotcandidatechar(pitvsdf, "PRL") + labs(subtitle = "Pituitary")
h <- plotcandidatechar(pitvsdm, "PRL") + labs(y = NULL, x = "") + labs(subtitle = " ")
i <- plotremoval(pitvsdf, "PRL")+ labs(y = NULL) + labs(subtitle = " ")
j <- plotremoval(pitvsdm, "PRL") + labs(y = NULL, x = "")+ labs(subtitle = " ")
k <- plotreplacement(pitvsdf, "PRL") + labs(y = NULL)+ labs(subtitle = " ")
l <- plotreplacement(pitvsdm, "PRL") + labs(y = NULL, x = "")+ labs(subtitle = " ")

m <- png::readPNG("../figures/venn-eggs-pit.png")
m <- ggdraw() +  draw_image(m, scale = 1)

n <- png::readPNG("../figures/venn-chicks-pit.png")
n <- ggdraw() +  draw_image(n, scale = 1)

ghi <- plot_grid(g,h,k,l,nrow = 1,
                labels = c("A", "", "B"), label_size = 8, hjust = 0,
                rel_widths = c(9,9,6,6))

ab <- plot_grid(a,b,rel_widths = c(1,2.5), 
                labels = c("C", "D"), label_size = 8, hjust = 0)

cd <- plot_grid(c,d,rel_widths = c(1.1,1), align = "h",
                labels = c("E", "F"), label_size = 8, hjust = 0)

ef <- plot_grid(f,m,n, rel_widths = c(2,1,1), nrow =1,
                labels = c("G", "H", "I"), label_size = 8, hjust = 0)


fig3 <- plot_grid(ghi,ab,cd,ef,  ncol = 1)
fig3


a <- plot.volcano("gonad", sexlevels,  "control_bldg") + 
  facet_wrap(~sex) 

b <- makebargraphv3(DEGcontrolreplace, "gonad","No. of DEGs\n(-) decreased  increased (+)", comparisonlevelscontrolreplace, comparisonlevelscontrolreplace) +
  labs(x = "Control versus all other reproductive and parental stages") 

c <- makebargraphv3(DEGbldg, "gonad", "No. of DEGs\n (-)decreased  increased (+)", comparisonlevelsbldg, comparisonlabelsbldg) +
  labs(x = "Nest-building versus all other parental stages")
d <- makebargraphv3(DEGchar, "gonad", NULL,  comparisonlevelscharnobldg, comparisonlabelscharnobldg) +
  labs(x = "Comparison of sequential parental stages")

e <- makebargraphv3(DEGremove, "gonad", "No. of DEGs\n (-)decreased  increased (+)", comparisonlevelsremoval, comparisonlevelsremoval) +
  labs(x = "Offspring removal versus temporal control")
f <- makebargraphv3(DEGreplace, "gonad", NULL,  comparisonlevelsreplace, comparisonlablelssreplace) +
  labs(x = "Offspring removal versus temporal or external control")

g <- plotcandidatechar(pitvsdf, "ESR1") + labs(subtitle = "Gonads")
h <- plotcandidatechar(pitvsdm, "ESR1") + labs(y = NULL, x = "") + labs(subtitle = " ")
i <- plotremoval(pitvsdf, "ESR1")+ labs(y = NULL) + labs(subtitle = " ")
j <- plotremoval(pitvsdm, "ESR1") + labs(y = NULL, x = "")+ labs(subtitle = " ")
k <- plotreplacement(pitvsdf, "ESR1") + labs(y = NULL)+ labs(subtitle = " ")
l <- plotreplacement(pitvsdm, "ESR1") + labs(y = NULL, x = "")+ labs(subtitle = " ")

m <- png::readPNG("../figures/venn-eggs-gon.png")
m <- ggdraw() +  draw_image(m, scale = 1)

n <- png::readPNG("../figures/venn-chicks-gon.png")
n <- ggdraw() +  draw_image(n, scale = 1)

ghi <- plot_grid(g,h,k,l,nrow = 1,
                labels = c("A", "", "B"), label_size = 8, hjust = 0,
                rel_widths = c(9,9,6,6))

ab <- plot_grid(a,b,rel_widths = c(1,2.5), 
                labels = c("C", "D"), label_size = 8, hjust = 0)

cd <- plot_grid(c,d,rel_widths = c(1.1,1), align = "h",
                labels = c("E", "F"), label_size = 8, hjust = 0)

ef <- plot_grid(f,m,n, rel_widths = c(2,1,1), nrow =1,
                labels = c("G", "H", "I"), label_size = 8, hjust = 0)


fig4 <- plot_grid(ghi,ab,cd,ef,  ncol = 1)
fig4


pdf(file="../figures/fig2-1.pdf", width=7, height=7)
plot(fig2)
dev.off()

png("../figures/fig2-1.png", width = 7, height = 7, 
    units = 'in', res = 300)
plot(fig2) 
dev.off()


pdf(file="../figures/fig3-1.pdf", width=7, height=7)
plot(fig3)
dev.off()

png("../figures/fig3-1.png", width = 7, height = 7, 
    units = 'in', res = 300)
plot(fig3) 
dev.off()

pdf(file="../figures/fig4-1.pdf", width=7, height=7)
plot(fig4)
dev.off()

png("../figures/fig4-1.png", width = 7, height = 7, 
    units = 'in', res = 300)
plot(fig4) 
dev.off()
```



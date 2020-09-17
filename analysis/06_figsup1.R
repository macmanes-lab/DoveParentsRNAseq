library(cowplot)

a <- png::readPNG("figures/fig_S1.png")
a <- ggdraw() +  draw_image(a, scale = 1)


b <- png::readPNG("figures/fig_S2.png")
b <- ggdraw() +  draw_image(b, scale = 1)


fig <- plot_grid(a,b, ncol = 1, rel_heights = c(1,1),
          labels = c("A", "B"), label_size = 8)
fig

pdf(file="figures/figsup1-1.pdf", width=5, height=5)
plot(fig)
dev.off()

png("figures/figsup1-1.png", width = 5, height = 5, 
    units = 'in', res = 300)
plot(fig) 
dev.off()

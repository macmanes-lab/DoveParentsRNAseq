library(colorblindr)
library(cowplot)

fig1a <- png::readPNG("../figures/images/fig_fig1b.png")
fig1a <- ggdraw() +  draw_image(fig1a, scale = 1)

fig1b <- png::readPNG("../figures/images/fig_fig1sup1.png")
fig1b <- ggdraw() +  draw_image(fig1b, scale = 1)

fig4a <- png::readPNG("../figures/images/fig_fig4a.png")
fig4a <- ggdraw() +  draw_image(fig4a, scale = 1)
 
fig1 <- plot_grid(fig1a, fig1b, 
                  nrow = 1, rel_widths = c(1,2.5))
fig <- plot_grid(fig1, fig4a, nrow = 2)

cvd_grid(fig)

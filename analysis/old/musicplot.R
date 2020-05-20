
musicplot <- function(whichgene, whichtissue){
  
  p <- candidatevsd %>%
    group_by(treatment, tissue, gene, sex)  %>% 
    summarize(median = median(counts, na.rm = T), 
              se = sd(counts,  na.rm = T)/sqrt(length(counts))) %>%
    dplyr::mutate(scaled = rescale(median, to = c(0, 7))) %>%
    dplyr::mutate(image = "../figures/images/musicnote.png")   %>%
    filter(
      gene %in% whichgene,
      tissue %in% whichtissue      ) %>% 
    collect() %>%
    drop_na() %>%
    ggplot( aes(x = treatment, y = median)) +
    geom_errorbar(aes(ymin = median - se, 
                      ymax = median + se, color = "white"),  width=0) +
    geom_image(aes(image=image), size = 0.1)+
    theme_B3() +
    theme(legend.position = "none",
          title = element_text(face = "italic"),
          strip.text = element_text(color = "white"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.line = element_blank(), axis.ticks = element_blank()) +
    scale_color_manual(values = allcolors) +
    labs(y = " ",x = NULL ) +
    facet_wrap(~sex, nrow = 2, scales = "free_y")
  return(p)
}

music <- png::readPNG("../figures/images/fig_music.png")
music <- ggdraw() +  draw_image(music, scale = 1)

musicalgenes <- png::readPNG("../figures/images/fig_musicalgenes.png")
musicalgenes <- ggdraw() +  draw_image(musicalgenes, scale = 1)


p1 <- musicplot(c("AVP", "OXT"), "hypothalamus") + labs(subtitle = "Hypothalamic AVP & OXT" ) + theme(axis.text.x = element_blank())
p2 <- musicplot(c("AVP", "OXT"), "pituitary") + labs(subtitle = "Pituitary AVP & OXT " ) + theme(axis.text.x = element_blank())
p3 <- musicplot(c("AVP", "OXT"), "gonad") + labs(subtitle = "Gonadal  AVP & OXT"  )

p123 <- plot_grid(music, p1, music, p3, ncol = 2, rel_widths = c(0.125,1))

plot_grid(musicalgenes, p123, nrow = 1, rel_widths = c(1.5,1))


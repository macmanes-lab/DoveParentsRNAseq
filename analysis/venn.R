

### removal overlap

```{r venn-rm}

filteredDEGs <- allDEG %>%
  filter(lfc > 0.14 | lfc < -0.14 ) %>%
  filter(tissue == "hypothalamus")

removalcomps <- c("inc.d3_m.inc.d3", "inc.d9_m.inc.d9", "inc.d17_m.inc.d17", "hatch_m.n2")

hatch_m.n2DEGs <- filteredDEGs %>%
  filter(comparison == "hatch_m.n2") %>%
  drop_na() %>%
  pull(gene)

inc.d17_m.inc.d17DEGs <- filteredDEGs %>%
  filter(comparison == "inc.d17_m.inc.d17") %>%
  drop_na() %>%
  pull(gene)

inc.d3_m.inc.d3DEGs <- filteredDEGs %>%
  filter(comparison == "inc.d3_m.inc.d3") %>%
  drop_na() %>%
  pull(gene)

inc.d9_m.inc.d9DEGs <- filteredDEGs %>%
  filter(comparison == "inc.d9_m.inc.d9") %>%
  drop_na() %>%
  pull(gene)


intersect(hatch_m.n2DEGs, inc.d17_m.inc.d17DEGs) %>%
  intersect(.,inc.d3_m.inc.d3DEGs) %>%
  intersect(.,inc.d9_m.inc.d9DEGs)

venn.diagram(
  x = list(inc.d3_m.inc.d3DEGs,hatch_m.n2DEGs,
           inc.d9_m.inc.d9DEGs, inc.d17_m.inc.d17DEGs),
  category.names = c("Inc 3" , "Hatch" , "Inc 9", " Inc 17"),
  filename = '../figures/venn-hyp.png',
  output=FALSE,
  print.mode = "raw",
  imagetype = "png",
  col=c("#CDCDCD", '#262625', '#959595','#959595'),
  fill = c("#CDCDCD", "#262625", "#959595", "#959595"),
  main = "Hypothalamus - Offspring Removal",
  sub = "Differentially expressed genes relative to internal controls"
)


```

## chicks / nesting care


```{r venn-chicks}

filteredDEGs <- allDEG %>%
  filter(lfc > 1.1 | lfc < -1.1 ) %>%
  filter(tissue == "hypothalamus")

chicks <- c("bldg_hatch", "bldg_n5", "bldg_n9",
            "bldg_extend", "bldg_early")

bldg_hatchDEGs <- filteredDEGs %>%
  filter(comparison == "bldg_hatch") %>%
  drop_na() %>%
  pull(gene)

bldg_n5DEGs <- filteredDEGs %>%
  filter(comparison == "bldg_n5") %>%
  drop_na() %>%
  pull(gene)

bldg_n9DEGs <- filteredDEGs %>%
  filter(comparison == "bldg_n9") %>%
  drop_na() %>%
  pull(gene)

bldg_extendDEGs <- filteredDEGs %>%
  filter(comparison == "bldg_extend") %>%
  drop_na() %>%
  pull(gene)

bldg_earlyDEGs <- filteredDEGs %>%
  filter(comparison == "bldg_early") %>%
  drop_na() %>%
  pull(gene)


# chicks
venn.diagram(
  x = list(bldg_n5DEGs, bldg_hatchDEGs, 
           bldg_earlyDEGs,bldg_extendDEGs,  bldg_n9DEGs),
  category.names = c( "N5" ,"Hatch" , "Early", " Extend", "N9"),
  filename = '../figures/venn-chicks-pit.png',
  output=FALSE,
  print.mode = "raw",
  imagetype = "png",
  height = 1500,
  width = 1750,
  resolution = 500,
  cex = 0.75,
  cat.cex = 0.75,
  col=c('#3182bd', "#6baed6",  '#cbc9e2', '#6a51a3', '#08519c'),
  fill = c("#3182bd", "#6baed6", "#cbc9e2", '#6a51a3', "#08519c" )
)

```


## eggs / incubation

```{r venn-eggs}

filteredDEGs <- allDEG %>%
  filter(lfc > 1.1 | lfc < -1.1 ) %>%
  filter(tissue == "hypothalamus") 

eggs <- c("bldg_lay", 
          "bldg_inc.d3", "bldg_inc.d9", "bldg_inc.d17",
          "bldg_prolong")

bldg_layDEGs <- filteredDEGs %>%
  filter(comparison == "bldg_lay") %>%
  drop_na() %>%
  pull(gene)

bldg_inc.d3DEGs <- filteredDEGs %>%
  filter(comparison == "bldg_inc.d3") %>%
  drop_na() %>%
  pull(gene)

bldg_inc.d9DEGs <- filteredDEGs %>%
  filter(comparison == "bldg_inc.d9") %>%
  drop_na() %>%
  pull(gene)

bldg_inc.d17DEGs <- filteredDEGs %>%
  filter(comparison == "bldg_inc.d17") %>%
  drop_na() %>%
  pull(gene)

bldg_prolongDEGs <- filteredDEGs %>%
  filter(comparison == "bldg_prolong") %>%
  drop_na() %>%
  pull(gene)

# chicks
venn.diagram(
  x = list(bldg_inc.d3DEGs, bldg_layDEGs,  bldg_prolongDEGs,
           bldg_inc.d17DEGs, bldg_inc.d9DEGs),
  category.names = c("Inc3" , "Lay" , "Prolong", " Inc17", "Inc9"),
  filename = '../figures/venn-eggs-pit.png',
  output=FALSE,
  print.mode = "raw",
  imagetype = "png",
  height = 1500,
  width = 1750,
  resolution = 500,
  cex = 0.75,
  cat.cex = 0.75,
  col=c( '#78c679', "#fed98e",  '#9e9ac8','#006837', '#31a354'),
  fill = c("#78c679", "#fed98e",  "#9e9ac8",'#006837', "#31a354")
)


m <- png::readPNG("figures/venn-eggs-hyp.png")
m <- ggdraw() +  draw_image(m, scale = 1)

n <- png::readPNG("figures/venn-chicks-hyp.png")
n <- ggdraw() +  draw_image(n, scale = 1)


```
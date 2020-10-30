library(tidyverse)
library(corrr)
library(cowplot)

source("R/genelists.R")
source("R/themes.R")

# sample information ----
samples <- read_csv("metadata/00_birds_sachecked.csv") %>%
  mutate(id = bird) %>%
  select(id, horm.id)
glimpse(samples)


# hormones ---
hormones <- read_csv("results/AllHorm_02042020_nomiss.csv") %>%
  mutate(id = str_replace_all(id, c("-" = ".", "/" = "."))) %>%
  mutate(prl = as.numeric(prl),
         cort = as.numeric(cort),
         p4 = as.numeric(p4),
         e2 = as.numeric(e2),
         t = as.numeric(t)) %>%
  full_join(., samples) %>%
  filter(treatment != "NA", treatment != ".") %>%
  mutate(e2t = ifelse(sex == "f", e2, t)) %>%
  mutate(sex = recode(sex, "f" = "female", "m" = "male" )) %>%
  select(id, treatment, sex, prl, cort, p4, e2t)  
glimpse(hormones)

# vsds (gene expression) ----
pitvsdAll <- read_csv("results/03_pitvsdAll.csv")
hypvsdAll <- read_csv("results/03_hypvsdAll.csv")
gonvsdAll <- read_csv("results/03_gonvsdAll.csv")

# join genes and hormones ----
genesnhomrmones <- function(vsds, whichsex){
  df <-  vsds %>% 
    inner_join(hormones, ., by = "id") %>%
    filter(sex %in% whichsex) 
  print(head(df)[1:10])
  return(df)
}

pitf <- genesnhomrmones(pitvsdAll, "female")
gonf <- genesnhomrmones(hypvsdAll, "female")
hypf <- genesnhomrmones(gonvsdAll, "female") 

pitm <- genesnhomrmones(pitvsdAll, "male")
gonm <- genesnhomrmones(hypvsdAll, "male")
hypm <- genesnhomrmones(gonvsdAll, "male") 

pit <- genesnhomrmones(pitvsdAll, c("female", "male"))
gon <- genesnhomrmones(hypvsdAll, c("female", "male"))
hyp <- genesnhomrmones(gonvsdAll, c("female", "male")) 

# correlations ---

plotcorrelation <- function(df, xvar, yvar, xlab, ylab){
  df %>%
    ggplot(aes(x = xvar, y = yvar)) +
    geom_point(aes(color = treatment)) +
    geom_smooth(method = "lm", aes(color = sex)) +
    scale_y_log10() +
    scale_x_log10() +
    labs(x = xlab, y = ylab) +
    scale_color_manual(values = allcolors) +
    theme(legend.position = "none") +
    facet_wrap(~sex)
}


a <- plotcorrelation(hyp,hyp$PRL, hyp$prl, 
                     "Hypothalamic prolactin gene expression", 
                     " ")

b <- plotcorrelation(pit,pit$PRL, pit$prl, 
                     "Pituitary prolactin gene expression", 
                     "Circulating prolactin (ng/mL)") +
  theme(strip.text = element_blank())

c <- plotcorrelation(gon,gon$PRL, gon$prl, 
                     "Gonadal prolactin gene expression", 
                     " ") +
  theme(strip.text = element_blank())

plot_grid(a,b,c, ncol = 1)



calculatecorrs <- function(df){
  df2 <- df %>%
    select(-id, -sex, -treatment) %>%
    correlate(.)
  print(head(df2)[1:6])
  return(df2)
}


hypfcorrs <- calculatecorrs(hypf)
hypmcorrs <- calculatecorrs(hypm)

pitfcorrs <- calculatecorrs(pitf)
pitmcorrs <- calculatecorrs(pitm)

gonfcorrs <- calculatecorrs(gonf)
gonmcorrs <- calculatecorrs(gonm)

gettop10PRL <- function(df){
  
  top5 <- df %>% focus(PRL) %>% arrange(desc(PRL)) %>% head(.,5)
  bottom5 <- df %>% focus(PRL) %>% arrange(PRL) %>% head(.,5)
  top10PRL <- rbind(top5, bottom5) 
  return(top10PRL)
}

a <- gettop10PRL(pitfcorrs) %>% rename("female" = rowname)
b <- gettop10PRL(pitmcorrs) %>% rename("male" = rowname)
cbind(a,b)


# candidate gene correlations
geneshormforcorr <- c(candidategenes, "prl", "cort", "p4", "e2t")
geneshormforcorr <- geneshormforcorr[! geneshormforcorr %in% 
                                       c("MC3R", "SERPINA4", "VIP")]


plotcorrrs <- function(df){
  p <- df %>% 
    filter(rowname %in% geneshormforcorr) %>%
    select(rowname, geneshormforcorr) %>% 
    rplot(colors = 
            c("skyblue1", "white","indianred2"),
          shape = 1) +
    theme_minimal(base_size = 9) +
    geom_point(size = 2) +
    theme(axis.text.x = element_text(angle = 90,
                                     face = "italic"),
          legend.position = "none") 
  return(p)
}


a <- plotcorrrs(hypfcorrs)
b <- plotcorrrs(hypmcorrs)
c <- plotcorrrs(pitfcorrs)
d <- plotcorrrs(pitmcorrs)
e <- plotcorrrs(gonfcorrs)
f <- plotcorrrs(gonmcorrs)

p <- plot_grid(a,b,c,d,e,f, ncol = 2)
p

png("~/Desktop/corrs.png", width = 10, height = 12, 
    units = 'in', res = 300)
plot(p) 
dev.off()


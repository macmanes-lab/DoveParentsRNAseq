
library(tidyverse)
library(randomForest)
library(caret)
require(caTools)
library(cowplot)

source("R/themes.R")


# https://towardsdatascience.com/random-forest-in-r-f66adf80ec9

genesnhomrmones <- full_join(hormones, candidatevsds, by = c("id", "sex", "treatment"))  %>%
  select(id, sex, tissue, treatment, 
         prl, cort, p4, e2t, everything()) %>%
  mutate(tissue = factor(tissue, levels = tissuelevels),
         treatment = factor(treatment, levels = alllevels)) %>% 
  mutate_if(is.character, as.factor) %>%
  drop_na() %>%
  mutate(external = fct_collapse(treatment,
                               eggs = c("lay", "inc.d3", "inc.d9", "inc.d17", "prolong"),
                               chicks = c("hatch", "n5", "n9", "extend" , "early" ),
                               loss = c("m.inc.d3",  "m.inc.d9",  "m.inc.d17",  "m.n2"),
                               controls = c("control", "bldg")),
       internal = fct_collapse(treatment,
                               earlier = c("lay", "inc.d3", "m.inc.d3",  "m.inc.d9",
                                         "inc.d9",  "early"),
                               later = c("hatch", "n5", "n9", "extend" ,  "prolong", "inc.d17",
                                        "m.inc.d17",  "m.n2"),
                               controls = c("control", "bldg")),
       hybrid = fct_collapse(treatment,
                               earlypresent = c("lay", "inc.d3", "inc.d9",  "early"),
                               earlyloss = c( "m.inc.d3",  "m.inc.d9"),
                               laterpresent = c("hatch", "n5", "n9", "extend" ,  "prolong", "inc.d17"),
                               laterloss = c("m.inc.d17",  "m.n2"),
                               controls = c("control", "bldg")),
       study = fct_collapse(treatment,
                             char = c("lay", "inc.d3", "inc.d9", "hatch", "n5", "n9", "inc.d17" ),
                             manip = c( "m.inc.d3",  "m.inc.d9", "m.inc.d17",  "m.n2",
                                        "early", "extend" ,  "prolong"),
                             controls = c("control", "bldg"))) %>%
  select(-id,  -MC3R, -CRH, -NPVF, -VIPR1) %>%
  select(treatment, external, internal, hybrid, study, sex, tissue, everything())
head(genesnhomrmones)




# need bird info to join with vsds 


birds <- read_csv("metadata/00_birds_fixed.csv") %>%
  mutate(external = fct_collapse(treatment,
                                 eggs = c("lay", "inc.d3", "inc.d9", "inc.d17", "prolong"),
                                 chicks = c("hatch", "n5", "n9", "extend" , "early" ),
                                 loss = c("m.inc.d3",  "m.inc.d9",  "m.inc.d17",  "m.n2"),
                                 controls = c("control", "bldg")),
         internal = fct_collapse(treatment,
                                 earlier = c("lay", "inc.d3", "m.inc.d3",  "m.inc.d9",
                                             "inc.d9",  "early"),
                                 later = c("hatch", "n5", "n9", "extend" ,  "prolong", "inc.d17",
                                           "m.inc.d17",  "m.n2"),
                                 controls = c("control", "bldg")),
         hybrid = fct_collapse(treatment,
                               earlypresent = c("lay", "inc.d3", "inc.d9",  "early"),
                               earlyloss = c( "m.inc.d3",  "m.inc.d9"),
                               laterpresent = c("hatch", "n5", "n9", "extend" ,  "prolong", "inc.d17"),
                               laterloss = c("m.inc.d17",  "m.n2"),
                               controls = c("control", "bldg")),
         study = fct_collapse(treatment,
                              char = c("lay", "inc.d3", "inc.d9", "hatch", "n5", "n9", "inc.d17" ),
                              manip = c( "m.inc.d3",  "m.inc.d9", "m.inc.d17",  "m.n2",
                                         "early", "extend" ,  "prolong"),
                              controls = c("control", "bldg")))

# now, get vsd for only DEGs

DEGlist <- read_csv("results/04_allDEG.csv") %>%
  filter(!grepl("control", comparison)) %>%
  select(gene) %>% arrange(gene) %>%
  filter(!grepl("LOC|\\.", gene)) %>%
  distinct(gene) %>% pull(gene)

slimDEGs <- DEGlist[5000:6000]

# select which genes to look at
savecols <- c("id", "group", slimDEGs)
savecols <- as.vector(savecols) 

vsd_path <- "results/"   # path to the data
vsd_files <- c("03_hypvsdAll.csv" , "03_pitvsdAll.csv", "03_gonvsdAll.csv")
vsd_pathfiles <- paste0(vsd_path, vsd_files)
vsd_pathfiles

DEGvsds <- vsd_pathfiles %>%
  map_dfr(~read_csv(.x), .id = "group") 
DEGvsds2 <- DEGvsds %>%
  dplyr::select(one_of(savecols))  %>%
  left_join(birds, ., by = "id") %>%
  mutate(group = factor(group),
         treatment = factor(treatment, levels = alllevels)) %>%
  mutate(tissue = fct_recode(group, "hypothathalamus" = "1",
                              "pituitary" = "2",
                              "gonads" = "3")) %>%
  drop_na(tissue) %>% 
  select(-id, -group) %>%
  select(treatment, external, internal, hybrid, study, sex, tissue, everything())  %>%
  select_if(~ !any(is.na(.))) %>%
  mutate_if(is.character, as.factor) 
head(DEGvsds2) 


sample = sample.split(DEGvsds2$treatment, SplitRatio = .75)
train = subset(DEGvsds2, sample == TRUE)
test  = subset(DEGvsds2, sample == FALSE)
dim(train)
dim(test)

rf1 <- randomForest(
  treatment ~ . * sex * tissue,
  data=train
)

rf1 # OOB estimate of  error rate: candidate genes 40.25%, Pgenes 35.22%
pred = predict(rf1, newdata=test[-1])
res1 <- test %>%
  mutate(pred = pred) %>%
  select(treatment, pred)
res1

rf2 <- randomForest(
  external ~ . * sex * tissue,
  data=train
)

rf2 # OOB estimate of  error rate: candidates: 1.24%, Pgens 2.03%
pred2 = predict(rf2, newdata=test[-2])
res2 <- test %>%
  mutate(pred = pred2) %>%
  select(external, pred)
res2


rf3 <- randomForest(
  internal ~ . * sex * tissue,
  data=train
)

rf3  # OOB estimate of  error rate: candidates: 1.24%, Pgenes 0%

pred3 = predict(rf3, newdata=test[-3])
res3 <- test %>%
  mutate(pred = pred3) %>%
  select(internal, pred)
res3



rf4 <- randomForest(
  hybrid ~ . * sex * tissue,
  data=train
)
rf4 

pred4 = predict(rf4, newdata=test[-4])
res4 <- test %>%
  mutate(pred = pred4) %>%
  select(hybrid, pred)
res4

rf5 <- randomForest(
  study ~ . * sex * tissue,
  data=train
)
rf5 

pred5 = predict(rf5, newdata=test[-5])
res5 <- test %>%
  mutate(pred = pred5) %>%
  select(study, pred)
res5

a <- ggplot(res1, aes(x = treatment, fill = pred)) +
  geom_bar(position = "fill") +
  theme_B3() +
  guides(fill=guide_legend(ncol=2)) +
  scale_fill_manual(values = allcolors) +
  labs(y = "percent predicted", 
       x = "parental stage or treatment",
       subtitle = "OOB estimate of  error rate: 40.25%",
       title = "Prediction based off 75/25 training/testing split") +
  expand_limits(y = 0)


b <- ggplot(res2, aes(x = external, fill = pred)) +
  geom_bar(position = "fill") +
  theme_B3() +
  scale_fill_manual(values = colorhypothesis) +
  labs(y = "percent predicted", x = "external stimuli",
       subtitle = "OOB estimate of  error rate: 2.03%") +
  expand_limits(y = 0)

c <- ggplot(res3, aes(x = internal, fill = pred)) +
  geom_bar(position = "fill") +
  theme_B3() +
  scale_fill_manual(values = colorhypothesis) +
  labs(y = "percent predicted", 
       x = "timing in parental care cycle",
       subtitle = "OOB estimate of  error rate: 1.24%") +
  expand_limits(y = 0)


e <- ggplot(res5, aes(x = study, fill = pred)) +
  geom_bar(position = "fill") +
  theme_B3() +
  #scale_fill_manual(values = colorhypothesis) +
  labs(y = "percent predicted", 
       x = "characterization versus manipulation",
       subtitle = "OOB estimate of  error rate: 5.22") +
  expand_limits(y = 0)

d <- ggplot(res4, aes(x = hybrid, fill = pred)) +
  geom_bar(position = "fill") +
  theme_B3() +
  #scale_fill_manual(values = colorhypothesis) +
  labs(y = "percent predicted", 
       x = "hybrid model internal and external",
       subtitle = "OOB estimate of  error rate: 0.72%") +
  expand_limits(y = 0)

bc <- plot_grid(b,c)
de <- plot_grid(e,d)
abc <- plot_grid(a,bc,de, nrow = 3,
                 rel_heights = c(1.2,1,1))

pdf(file="figures/fig-randomforest.pdf", width=7, height=6)
plot(abc)
dev.off()

png("figures/fig-randomforest.png", width = 7, height = 6, 
    units = 'in', res = 300)
plot(abc) 
dev.off()
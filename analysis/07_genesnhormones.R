library(tidyverse)

hormones <- read_csv("results/AllHorm_02042020_nomiss.csv") %>%
  select(id, prl) %>%
  mutate(id = str_replace_all(id, c("-" = ".", "/" = "."))) %>%
  mutate(prl = as.numeric(prl))

genes <- read_csv("results/03_pitvsdAll.csv") %>%
  mutate(id = sapply(strsplit(samples, '\\_'), "[", 1),
         id = str_replace_all(id, c(".ATLAS" = "")),
         id = tolower(id)) 

genesnhomrones <- full_join(genes, hormones, by = "id") %>%
  arrange(id) 
head(genesnhomrones)

genesnhomrones %>%
  select(prl, PRL) %>%
  drop_na() %>%
  ggplot(aes(x = prl, y = PRL)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "circulating prolactin (ng/mL)",
       y = "prolactin gene expression (counts)")

genesnhomrones %>%
  select(prl, STT3B) %>%
  drop_na() %>%
  ggplot(aes(x = prl, y = STT3B)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "circulating prolactin (ng/mL)",
       y = "CAMKK2 gene expression (counts)")

prlcor <- genesnhomrones %>%
  select(-samples, -id) %>%
  correlate(.) %>%
  focus(prl) 

prlcor %>% arrange(desc(prl))
prlcor %>% arrange(prl)

#write.csv(genesnhomrones, "~/Desktop/genesnhomrones.csv")
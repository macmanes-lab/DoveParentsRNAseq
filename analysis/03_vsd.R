library(tidyverse)
library(readr)

source("R/themes.R")
source("R/genelists.R")

### variance stabilized gene expression  (vsd) 


vsd_path <- "results/DEseq2/treatment/"   # path to the data
vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
vsd_pathfiles <- paste0(vsd_path, vsd_files)
vsd_files

allvsd <- vsd_pathfiles %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x), .id = "file_name")  %>% 
  dplyr::rename("gene" = "X1") 

## before pivoting, check names of df with
## head(names(allvsd)) and tail(names(allvsd))

allvsd <- allvsd %>% 
  pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
               names_to = "samples", values_to = "counts")  %>% 
  drop_na() 
allvsd <- as_tibble(allvsd)
head(allvsd[2:4])
tail(allvsd[2:4])

# for plotting gene expression over time
getcandidatevsdmanip <- function(whichgenes, whichtissue){
  
  candidates  <- allvsd %>%
    filter(gene %in% whichgenes) %>%
    dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1),
                  sextissue = sapply(strsplit(sextissue, 'results/DEseq2/treatment/'), "[", 2),
                  sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                  sex = factor(sex),
                  tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                  treatment = sapply(strsplit(samples, '\\_'), "[", 4),
                  treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1),
                  treatment = recode(treatment, "m.hatch" = "m.n2", "extend.hatch" = "m.n2",
                                     "inc.prolong" = "prolong", "m.inc.d8" = "early"),
                  treatment =  factor(treatment, levels = alllevels))  %>%
    dplyr::select(gene, counts, sex, tissue, treatment, samples)  %>%
    drop_na()
  return(candidates)
}

candidatevsd <- getcandidatevsdmanip(candidategenes)

hypvsdf <- candidatevsd %>% filter(tissue %in% "hypothalamus", sex == "female")  
pitvsdf <- candidatevsd %>% filter(tissue %in% "pituitary", sex == "female")
gonvsdf <- candidatevsd %>% filter(tissue %in% "gonads", sex == "female")  
hypvsdm <- candidatevsd %>% filter(tissue %in% "hypothalamus", sex == "male")  
pitvsdm <- candidatevsd %>% filter(tissue %in% "pituitary", sex == "male") 
gonvsdm <- candidatevsd %>% filter(tissue %in% "gonads", sex == "male") 

write.csv(hypvsdf, "results/03_hypvsdf.csv", row.names = F)
write.csv(pitvsdf, "results/03_pitvsdf.csv", row.names = F)
write.csv(gonvsdf, "results/03_gonvsdf.csv", row.names = F)
write.csv(hypvsdm, "results/03_hypvsdm.csv", row.names = F)
write.csv(pitvsdm, "results/03_pitvsdm.csv", row.names = F)
write.csv(gonvsdm, "results/03_gonvsdm.csv", row.names = F)
write.csv(candidatevsd, "results/03_candidatevsd.csv", row.names = F)

shinyvsd <- getcandidatevsdmanip(shinygenes)
write.csv(shinyvsd, "../musicalgenes/data/candidatecounts.csv", row.names = F)

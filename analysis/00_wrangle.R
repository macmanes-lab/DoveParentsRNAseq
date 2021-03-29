# Data wrangling
library(tidyverse)

source("R/themes.R")

## import Kallisto transcript data, make gene info file 

# import count data, set rows at entreziz

print("reading kallisto data")
kallistodata <- read.table("results/kallistocounts.txt", 
                           sep = ",", row.names = NULL) %>%
  dplyr::rename("NCBI" = "entrezid", "gene" = "Name") %>%
  select(-row.names) %>% 
  filter(gene != "") %>%  # rm sample without valid identifier
  as_tibble()
head(kallistodata)

## save gene info
geneinfo <- kallistodata %>%
  select(gene, geneid, NCBI)

# for example, these gene have 2 and 3 isoforms
geneinfo %>% filter(gene == "GRIN1")
geneinfo %>% filter(gene == "CACNA1C")

## check for genes that have multple transcripts expressed

print("creating isoforms df")
isoforms <- kallistodata %>%
  group_by(gene) %>%
  summarize(n = n()) %>%
  filter(n > 1) %>%
  group_by(n) %>% 
  summarize(genes = str_c(gene, collapse = ", ")) %>%
  arrange(desc(n)) %>%
  rename("n(isoforms)" = "n") %>% 
  mutate("n(counts)" = str_count(genes, ",") + 1)
head(isoforms)


# aggregate transcript counts to gene counts
print("aggregating count data")

countData <- kallistodata %>% 
  select(-geneid, -NCBI) %>% 
  pivot_longer(-gene, names_to = "samples", values_to = "counts") %>%  
  pivot_wider(
    names_from = samples, 
    values_from = counts,
    values_fn = list(counts = sum))  %>% 
  arrange(gene) %>% 
  select(-contains("R2XR"))  #remove duplicate samples
countData <- as.data.frame(countData)
row.names(countData) <- countData$gene ## make gene the row name
countData[1] <- NULL ## make gene the row name
countData <- round(countData) #round all value for deseq2
head(countData[1:2])

## remove duplicate bird
names(countData)

# print tolal num of genes and samples
dim(countData)

## wrangle colData 

print("creating colData file")

colData <- read.table(file.path( "metadata/kallistosamples.txt"),
                      header = F, stringsAsFactors = F) %>%
  # use strsplit to cut the filename into meaningful columns
  mutate(bird = sapply(strsplit(V1,'\\_'), "[", 1),
         sex = sapply(strsplit(V1,'\\_'), "[", 2),
         tissue = sapply(strsplit(V1,'\\_'), "[", 3),
         temp = sapply(strsplit(V1,'\\_'), "[", 4)) %>%
  mutate(treatmenttemp = sapply(strsplit(temp,'\\.'), "[", 1),
         NYNO = sapply(strsplit(temp,'\\.'), "[", 2)) %>%
  
  # rename variables
  mutate(treatment = ifelse(grepl("extend-hatch", treatmenttemp), "extend",
                            ifelse(grepl("inc-prolong", treatmenttemp), "prolong",
                                   ifelse(grepl("m.hatch", treatmenttemp), "m.n2",
                                          ifelse(grepl("m.inc.d8", treatmenttemp), "early",
                                                 treatmenttemp))))) %>%
  select(-temp, -NYNO, -treatmenttemp ) %>%
  # replace dashes with periods (to make deseq happy)
  mutate(bird = gsub("-", ".", bird),
         treatment = gsub("-", ".", treatment),
         V1 = gsub("-", ".", V1)) %>%
  
  ## fix sample assignments
  
  mutate(sex = ifelse(grepl("blu119.w84.x", bird), "male", sex),
         treatment = ifelse(grepl("y128.g23.x|blu108.w40.o158|x.blu43.g132|r37.w100.x", 
                                  bird), "m.inc.d9", treatment)) %>%
  filter(bird != "r.r.x.ATLAS.R2XR") %>%
  
  mutate(bird = tolower(bird),
         sex = factor(sex),
         tissue = factor(tissue, levels = c("hypothalamus", "pituitary", "gonad"))) %>% 
  mutate(tissue = fct_recode(tissue, "gonads" = "gonad"),
         treatment = factor(treatment, levels = alllevels)) %>%

  ## add external and internalhypothesis
  mutate(external = fct_collapse(treatment,
                                 eggs = c("lay", "inc.d3", "inc.d9", "inc.d17", "prolong"),
                                 chicks = c("hatch", "n5", "n9", "extend" , "early" ),
                                 nest = c("bldg", "m.inc.d3",  "m.inc.d9",  "m.inc.d17",  "m.n2"),
                                 control = c("control")),
         internal = fct_collapse(treatment,
                                 early = c("lay", "inc.d3", "bldg", "m.inc.d3",  "m.inc.d9",
                                           "inc.d9",  "early"),
                                 late = c("hatch", "n5", "n9", "extend" ,  "prolong", "inc.d17",
                                          "m.inc.d17",  "m.n2"),
                                 control = c("control"))) %>%
  mutate(group = paste(sex, tissue, treatment, sep = "."),
         study = ifelse(grepl("m.|extend|prolong", treatment), 
                        "manipulation", "charcterization")) %>%
  mutate(study = factor(study, levels = c("charcterization", "manipulation")))
str(colData)
head(colData)

## check that rownames and colnames match for DESeq
ncol(countData) == nrow(colData)

## summarize data

colData %>% select(sex, tissue, treatment, study)  %>%  summary()
table(colData$sex, colData$treatment, colData$tissue)

colData %>%
  count(sex, tissue, treatment, sort = TRUE)

length(unique(colData$bird))
table(colData$sex, colData$tissue)


## save bird data to join with 

print("creating sample file for future hamonoization")

samples <- read_csv("metadata/00_birds_sachecked.csv") %>%
  select(-treatment, -sex) %>%
  mutate(bird = tolower(bird)) %>%
  full_join(., colData, by = c("bird")) %>%
  dplyr::rename("rnaseq.id" = "V1") %>%
  select(bird, horm.id, rnaseq.id, sex, tissue, treatment, external, internal, group) %>%
  distinct()  %>%
  filter(rnaseq.id != "NA")
head(samples)

## save files for downstream use

print("saving files")

write.csv(countData, "results/00_countData.csv")
write.csv(colData, "metadata/00_colData.csv")
write.csv(samples, "metadata/00_samples.csv")
write.csv(geneinfo, "metadata/00_geneinfo.csv", row.names = TRUE)
write.csv(isoforms, "results/00_isoforms.csv", row.names = TRUE)




# Data wrangling
library(tidyverse)

## import Kallisto transcript data, make gene info file 

# import count data, set rows at entreziz
print("reading kallisto data")
kallistodata <- read.table("results/kallistocounts.txt", 
                           sep = ",", row.names = NULL) 
kallistodata <- kallistodata %>%
  dplyr::rename("NCBI" = "entrezid",
                "gene" = "Name")
head(kallistodata[1:3])

## check for genes that have multple transcripts expressed

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

# for example, these gene have 2 and 3 isoforms
geneinfo %>% filter(gene == "GRIN1")
geneinfo %>% filter(gene == "CACNA1C")

## group data by gene

# aggregate transcript counts to gene counts
countData <- kallistodata %>% 
  select(-row.names, -geneid, -NCBI) %>% 
  pivot_longer(-gene, names_to = "samples", values_to = "counts") %>%  
  pivot_wider(
    names_from = samples, 
    values_from = counts,
    values_fn = list(counts = sum))  %>% 
  arrange(gene) %>% 
  filter(gene != "")
countData <- as.data.frame(countData)
row.names(countData) <- countData$gene ## make gene the row name
countData[1] <- NULL ## make gene the row name
countData <- round(countData) #round all value to nearest 1s place
head(countData[13:15])

# print tolal num of genes and samples
dim(countData)

## wrangle colData 

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
         V1 = gsub("-", ".", V1))
# specify levels that will be factors
cols <- c("sex", "treatment", "tissue")
colData[cols] <- lapply(colData[cols], factor) 
colData <- colData %>%
  mutate(tissue = factor(tissue, levels = c("hypothalamus", "pituitary", "gonad"))) %>% 
  mutate(tissue = fct_recode(tissue, "gonads" = "gonad")) %>%
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
colData <- as.data.frame(colData)
row.names(colData) <- colData$V1
str(colData)
head(colData)

## check that rownames and colnames match for DESeq
ncol(countData) == nrow(colData)

## summarize data

colData %>% select(sex, tissue, treatment, study)  %>%  summary()
table(colData$sex, colData$treatment, colData$tissue)

length(unique(colData$bird))
table(colData$sex, colData$tissue)

## save bird data 

birds <- colData %>%
  select(bird, sex, treatment) %>%
  distinct()
head(birds)

## save files for downstream use

write.csv(hugopresent, "metadata/00_hugopresent.csv")
write.csv(hugoabsent, "metadata/00_hugoasent.csv")

write.csv(countData, "results/00_counts.csv")
write.csv(isoforms, "results/00_geneswithisoforms.csv", row.names = T)
write.csv(colData, "metadata/00_colData.csv")
write.csv(birds, "metadata/00_birds.csv")
write.csv(geneinfo, "metadata/00_geneinfo.csv", row.names = TRUE)
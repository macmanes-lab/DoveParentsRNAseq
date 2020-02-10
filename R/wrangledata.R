# data wrangle for figure 1

# Note, the input file is too big for storage on github :(

pseudocounts <- read_csv("../results/01_pseudo.counts.csv")
head(pseudocounts[1:3])
pseudocounts <- as.data.frame(pseudocounts)
row.names(pseudocounts) <- pseudocounts$X1
pseudocounts$X1 <- NULL
# prep count data for all samples
countData <- as.data.frame(t(pseudocounts))
head(countData[1:3])

# prep col data for all samples
colData <- read.csv("../metadata/00_samples.csv", header = T, row.names = 1)
colData$treatment <- factor(colData$treatment, levels = alllevels)
colData <- colData %>% mutate(tissue = fct_recode(tissue, "gonads" = "gonad"))
colData$tissue <- factor(colData$tissue, levels = tissuelevels)
row.names(colData) <- colData$V1

# check ready for analysis
#row.names(countData) == row.names(colData)


## PRL vsd

vsd_path <- "../results/DEseq2/"   # path to the data
vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
vsd_pathfiles <- paste0(vsd_path, vsd_files)
vsd_files

PRLvsd <- vsd_pathfiles %>%
  setNames(nm = .) %>% 
  map_df(~read_csv(.x), .id = "file_name")  %>% 
  dplyr::rename("gene" = "X1") %>% 
  pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
               names_to = "samples", values_to = "counts") %>%
  filter(gene == "PRL")  %>%
  dplyr::mutate(sextissue = sapply(strsplit(file_name, '_vsd.csv'), "[", 1)) %>%
  dplyr::mutate(sextissue = sapply(strsplit(sextissue, '../results/DEseq2/'), "[", 2)) %>%
  dplyr::mutate(sex = sapply(strsplit(sextissue, '\\_'), "[", 1),
                tissue = sapply(strsplit(sextissue, '\\_'), "[", 2),
                treatment = sapply(strsplit(samples, '\\_'), "[", 4)) %>%
  dplyr::mutate(treatment = sapply(strsplit(treatment, '.NYNO'), "[", 1)) %>%
  dplyr::select(sex, tissue, treatment, gene, samples, counts)
PRLvsd$treatment <- factor(PRLvsd$treatment, levels = alllevels)

PRLhyp <- PRLvsd %>% filter(tissue == "hypothalamus")
PRLpit <- PRLvsd %>% filter(tissue == "pituitary")
PRLgon <- PRLvsd %>% filter(tissue == "gonad")


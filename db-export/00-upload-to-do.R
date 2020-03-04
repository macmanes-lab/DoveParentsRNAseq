library(purrr)
library(dplyr)
library(tidyr)
library(data.table)
library(RPostgreSQL)

# allvsd ----

geneids <- fread("metadata/00_geneinfo.csv")

vsd_path <- "results/DEseq2/"   # path to the data
vsd_files <- dir(vsd_path, pattern = "*vsd.csv") # get file names
vsd_pathfiles <- paste0(vsd_path, vsd_files)

# vsd_files

allvsd <- vsd_pathfiles %>%
  setNames(nm = .) %>% 
  map_df(~fread(.x), .id = "file_name")  %>% 
  rename("gene" = "V1") %>% 
  pivot_longer(cols = L.G118_female_gonad_control:y98.o50.x_male_pituitary_inc.d3, 
               names_to = "samples", values_to = "counts")

# alldeg ----

DEG_path <- "results/DEseq2/"   # path to the data
DEG_files <- dir(DEG_path, pattern = "*DEGs") # get file names
DEG_pathfiles <- paste0(DEG_path, DEG_files)
#DEG_files

source("R/themes.R")

allDEG <- DEG_pathfiles %>%
  setNames(nm = .) %>% 
  map_df(~fread(.x), .id = "file_name") %>% 
  as_tibble() %>% 
  
  mutate(
    DEG = sapply(strsplit(as.character(file_name),'results/DEseq2/'), "[", 2),
    DEG = sapply(strsplit(as.character(DEG),'_diffexp.csv'), "[", 1),
    tissue = sapply(strsplit(as.character(DEG),'\\.'), "[", 1),
    down = sapply(strsplit(as.character(DEG),'\\_'), "[", 3),
    up = sapply(strsplit(as.character(DEG),'\\_'), "[", 4),
    comparison = paste(down,up, sep = "_"),
    sex = sapply(strsplit(as.character(sextissue),'\\_'), "[", 1),
    tissue = sapply(strsplit(as.character(sextissue),'\\_'), "[", 2)
  ) %>% 
  
  select(sex,tissue,comparison, direction, gene, lfc, padj, logpadj)
  
allDEG <- allDEG %>% 
  mutate(
    tissue = factor(tissue, levels = tissuelevel),
    comparison = factor(comparison , levels = comparisonlevels),
    direction = factor(direction, levels = charlevels)
  )

# upload to db ----

con <- DBI::dbConnect(
  DBI::dbDriver("PostgreSQL"),
  host = Sys.getenv("psql_rayna_ip"),
  port = Sys.getenv("psql_rayna_port"),
  user = Sys.getenv("psql_rayna_usr"),
  password = Sys.getenv("psql_rayna_pwd"),
  dbname = Sys.getenv("psql_rayna_db")
)

dbSendQuery(con, "DROP TABLE IF EXISTS public.allvsd")
dbSendQuery(con, "DROP TABLE IF EXISTS public.alldeg")

dbSendQuery(
  con,
  "CREATE TABLE public.allvsd
      (
      file_name varchar(255) DEFAULT NULL,
      gene varchar(255) NOT NULL,
      samples varchar(255) DEFAULT NULL,
      counts decimal(16,2) DEFAULT NULL,
      CONSTRAINT allvsd_pk PRIMARY KEY (file_name, gene, samples)
      )"
)

dbSendQuery(
  con,
  "CREATE TABLE public.alldeg
      (
      sex varchar(255) DEFAULT NULL,
      tissue varchar(255) NOT NULL,
      comparison varchar(255) DEFAULT NULL,
      direction varchar(255) DEFAULT NULL,
      gene varchar(255) DEFAULT NULL,
      lfc decimal(16,2) DEFAULT NULL,
      padj decimal(16,2) DEFAULT NULL,
      logpadj decimal(16,2) DEFAULT NULL,
      CONSTRAINT alldeg_pk PRIMARY KEY (sex, tissue, comparison, direction, gene)
      )"
)

obs_allvsd <- as.numeric(dbGetQuery(con, "SELECT COUNT(*) FROM public.allvsd"))

if (obs_allvsd == 0) {
  dbWriteTable(con, "allvsd", allvsd, append = TRUE, overwrite = FALSE, row.names = FALSE)
}

obs_alldeg <- as.numeric(dbGetQuery(con, "SELECT COUNT(*) FROM public.alldeg"))

if (obs_alldeg == 0) {
  dbWriteTable(con, "alldeg", mutate_if(allDEG, is.factor, as.character),
               append = TRUE, overwrite = FALSE, row.names = FALSE)
}

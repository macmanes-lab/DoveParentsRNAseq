library(purrr)
library(dplyr)
library(tidyr)
library(data.table)
library(RPostgreSQL)

source("R/themes.R")

###### variance statbilized data

# allvsd ----

geneids <- fread("metadata/00_geneinfo.csv")

vsd_path <- "results/"   # path to the data
vsd_files <- dir(vsd_path, pattern = "*vsd[fm].csv") # get file names
vsd_pathfiles <- paste0(vsd_path, vsd_files)

# vsd_files

allvsd <- vsd_pathfiles %>%
  setNames(nm = .) %>% 
  map_df(~fread(.x), .id = "file_name") 
# alldeg ----

allDEG <- read_csv("results/03_allDEG.csv")

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

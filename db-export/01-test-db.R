library(RPostgreSQL)
library(dplyr)
library(dbplyr)

con <- dbConnect(
  dbDriver("PostgreSQL"),
  host = Sys.getenv("psql_rayna_ip"),
  port = Sys.getenv("psql_rayna_port"),
  user = Sys.getenv("psql_rayna_usr"),
  password = Sys.getenv("psql_rayna_pwd"),
  dbname = Sys.getenv("psql_rayna_db")
)

tbl(con, in_schema("public", "alldeg")) %>% 
  filter(gene == "PRL")

df <- as_tibble(tbl(con, in_schema("public", "alldeg")))
head(df)
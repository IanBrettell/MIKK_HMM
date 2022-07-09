# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
IN = list(
    "/hps/nobackup/birney/users/ian/MIKK_all/0.08/F0.csv",
    "/hps/nobackup/birney/users/ian/MIKK_all/0.08/F2.csv",
    "/hps/nobackup/birney/users/ian/MIKK_all/0.08/Kiyosu_CC.csv"
)

## True
IN = snakemake@input
OUT = snakemake@output[[1]]

names(IN) = IN %>% 
  unlist() %>% 
  basename() %>% 
  stringr::str_remove(".csv")

# Read in files and combine

df = purrr::map_dfr(IN, readr::read_csv, .id = "dataset")

# Write to file

readr::write_csv(df, OUT)

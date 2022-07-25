# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug

#IN = list("/hps/nobackup/birney/users/ian/somites/gcta/mlma_loco/permuted/hdrr/None/5000/1/intercept/None/1.loco.mlma",
#          "/hps/nobackup/birney/users/ian/somites/gcta/mlma_loco/permuted/hdrr/None/5000/1/intercept/None/2.loco.mlma")

## True

IN = snakemake@input
OUT = snakemake@output[[1]]

# Read in files and process

names(IN) = IN %>% 
  unlist() %>% 
  basename(.) %>% 
  stringr::str_remove(".loco.mlma")

df = purrr::map_dfr(IN, readr::read_tsv, .id = "SEED") %>% 
  dplyr::group_by(SEED) %>% 
  dplyr::summarise(MIN_P = min(p))

# Write to file

readr::write_csv(df, OUT)

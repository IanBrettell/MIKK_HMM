# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug

IN = list("/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco/permuted/hdrr/5000/0.8/dge/notrans/novel_object/1/All/1.loco.mlma",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco/permuted/hdrr/5000/0.8/dge/notrans/novel_object/1/All/2.loco.mlma")

## True

IN = snakemake@input
OUT = snakemake@output[[1]]

# Read in files and process

names(IN) = IN %>% 
  unlist() %>% 
  basename(.) %>% 
  stringr::str_remove(".loco.mlma")

df = purrr::map_dfr(IN, function(FILE) {
  readr::read_tsv(FILE,
                  col_types = c("iciccdddd"))
  }, .id = "SEED") %>% 
  dplyr::group_by(SEED) %>% 
  dplyr::summarise(MIN_P = min(p, na.rm = T))

# Write to file

readr::write_csv(df, OUT)

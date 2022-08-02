# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug

IN = list("/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco_mean_speed_consol/permuted/hdrr/0.05/5000/0.8/sge/invnorm/novel_object/time-quadrant/1.mlma",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco_mean_speed_consol/permuted/hdrr/0.05/5000/0.8/sge/invnorm/novel_object/time-quadrant/2.mlma")

## True

IN = snakemake@input
OUT = snakemake@output[[1]]

# Read in files and process

names(IN) = IN %>% 
  unlist() %>% 
  basename(.) %>% 
  stringr::str_remove(".mlma")

df = purrr::map_dfr(IN, function(FILE) {
  readr::read_tsv(FILE,
                  col_types = c("iciccdddd"))
  }, .id = "SEED") %>% 
  dplyr::group_by(SEED) %>% 
  dplyr::summarise(MIN_P = min(p, na.rm = T))

# Write to file

readr::write_csv(df, OUT)

# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco_mean_speed/true/hdrr/0.05/5000/0.8/dge/invnorm/novel_object/time-quadrant",
                pattern = ".mlma",
                full.names = T)
## True
IN = snakemake@input
OUT = snakemake@output[[1]]

# Read in files

out = purrr::map_dfr(IN,
                     ~readr::read_tsv(.x,
                                      col_types = c("iciccdddd"))) %>% 
  # filter out NaNs
  dplyr::filter(!is.na(p))


# Write to file

readr::write_tsv(out, OUT)

# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

IN = list.files("/hps/nobackup/birney/users/ian/MIKK_HMM/coverage/hdrr/bwamem2",
                full.names = T)


dat_list = readr::read_tsv(IN[1])

# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
IN = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/1/dist_angle/5.csv"

## True
IN = snakemake@input[[1]]
OUT_F0 = snakemake@output[["F0"]]
OUT_F2 = snakemake@output[["F2"]]
OUT_Kiyosu_CC = snakemake@output[["Kiyosu_CC"]]

# Read in data

df = readr::read_csv(IN)

# Split by generation

out = split(df, f = df$dataset)

# Write to files

readr::write_csv(out[["F0"]], OUT_F0)
readr::write_csv(out[["F2"]], OUT_F2)
readr::write_csv(out[["Kiyosu_CC"]], OUT_Kiyosu_CC)

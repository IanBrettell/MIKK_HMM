# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = here::here("config/F2_samples_converted.csv")

## True
IN = snakemake@input[[1]]
OUT = snakemake@output[[1]]

# Get counts

df = readr::read_csv(IN) %>% 
  dplyr::count(pat_line, mat_line) %>% 
  readr::write_csv(OUT)

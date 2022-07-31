# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = list.files(here::here("config/samples_no_calls/5000/0.8"),
                full.names = T)

## True
IN = snakemake@input
OUT = snakemake@output[["csv"]]


names(IN) = IN %>% 
  unlist() %>% 
  basename() %>% 
  stringr::str_remove(".csv")

# Read in files

df = purrr::map_dfr(IN, function(CSV){
  readr::read_csv(CSV,
                  col_names = "SAMPLE",
                  col_types = "i")
}, .id = "CHROM") %>% 
  # get unique samples
  dplyr::distinct(SAMPLE) %>% 
  dplyr::arrange(SAMPLE) %>% 
  dplyr::select(SAMPLE)

# Write to file
readr::write_csv(df, OUT)

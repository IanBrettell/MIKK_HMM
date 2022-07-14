# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = "/nfs/research/birney/projects/indigene/raw_data/genecore/ian_behaviour/AAAMKYVHV/aspera_AAAMKYVHV"

## True
IN = snakemake@input[["in_dir"]]
OUT = snakemake@output[["out"]]

# List files

df = tibble::tibble(PATH = list.files(IN, full.names = T)) %>% 
  # get basename
  dplyr::mutate(BASENAME = basename(PATH)) %>% 
  # extract sample ID
  tidyr::separate(col = "BASENAME",
                  into =  c(rep(NA, 5), "ID", "PAIR", NA),
                  sep = "_") %>% 
  # remove 'lane1' string in `ID`
  dplyr::mutate(NEW_ID = stringr::str_remove(ID, "lane1")) %>% 
  # pivot wider to put PAIRS into separate rows
  tidyr::pivot_wider(names_from = "PAIR",
                     values_from = "PATH") %>% 
  # rename columns
  dplyr::select(sample = NEW_ID,
                fq1 = `1`,
                fq2 = `2`)

# Write to file

readr::write_tsv(df, OUT)

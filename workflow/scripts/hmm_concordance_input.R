# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
#IN_FILE = "/hps/nobackup/birney/users/ian/pilot/merged/1.csv"

## True
IN_FILE = snakemake@input[[1]]
OUT_A = snakemake@output[["A"]]
OUT_B = snakemake@output[["B"]]

# Read in file

df_list = readr::read_csv(IN_FILE) %>% 
  # create column with assay/date/time/quadrant/fish to isolate
  # each individual fish track
  tidyr::unite(col = adtqf,
               assay, date, time, quadrant, fish,
               sep = "_",
               remove = F) %>% 
  # split by `dtq`
  split(., f = .$adtqf)

# Randomise order of tracks

index = 1:length(df_list)
irand = sample(index)

df_rand = df_list[irand]

# Split into 2 (test/train and vice versa)

half = round(length(df_rand) / 2)

randA = df_rand[1:half] %>% 
  dplyr::bind_rows() %>% 
  # remove `adtqf` column
  dplyr::select(-adtqf)
randB = df_rand[(half + 1):length(df_rand)] %>% 
  dplyr::bind_rows() %>% 
  # remove `adtqf` column
  dplyr::select(-adtqf)

# Write to files

readr::write_csv(randA, OUT_A)
readr::write_csv(randB, OUT_B)


# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
PED = "/hps/nobackup/birney/users/ian/MIKK_HMM/peds/F2/hdrr/5000/0.8.ped"
SAMPLES_FILE = "config/F2_samples_converted.csv"
ASSAY = "open_field" %>% 
  stringr::str_replace('_', ' ')

## True
PED = snakemake@input[["ped"]]
SAMPLES_FILE = snakemake@input[["samples_file"]]
OUT_CAT = snakemake@output[["cat"]]
OUT_QUANT = snakemake@output[["quant"]]
ASSAY = snakemake@params[["assay"]] %>% 
  stringr::str_replace('_', ' ')

# Read in files

pedp = readr::read_tsv(PED,
                       col_names = "SAMPLE",
                       col_types = "i",
                       col_select = 1)

samples = readr::read_csv(SAMPLES_FILE) %>% 
#  # create `indiv` column
#  tidyr::unite(date, time, quadrant,
#               col = "indiv",
#               sep = "_",
#               remove = F) %>% 
  # rename `finclip_id` as SAMPLE to match `pedp`
  dplyr::rename(SAMPLE = finclip_id) %>%
  # select key columns
  dplyr::select(SAMPLE, date, time, quadrant, tank_side)

# Join

covar = pedp %>% 
  # join samples with ID and individual
  dplyr::left_join(.,
                   samples,
                   by = "SAMPLE") %>% 
  dplyr::mutate(IID = SAMPLE) %>% 
  dplyr::select(SAMPLE, IID, everything()) 

# Split into categorical and quantitative variables
covar_cat = covar %>% 
  dplyr::select(SAMPLE, IID, date, quadrant, tank_side)

covar_quant = covar %>% 
  dplyr::select(SAMPLE, IID, time)

readr::write_tsv(covar_cat, OUT_CAT, col_names = F)
readr::write_tsv(covar_quant, OUT_QUANT, col_names = F)

# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
PED = "/hps/nobackup/birney/users/ian/MIKK_HMM/peds/F2/hdrr/5000/0.8.ped"
PHENO_FILE = "/hps/nobackup/birney/users/ian/MIKK_HMM/state_freq_F2/dist_angle/0.05/15/dge/invnorm/1.csv"
SAMPLES_FILE = "config/F2_samples_converted.csv"
ASSAY = "open_field" %>% 
  stringr::str_replace('_', ' ')

## True
PED = snakemake@input[["ped"]]
PHENO_FILE = snakemake@input[["phenos"]]
SAMPLES_FILE = snakemake@input[["samples_file"]]
OUT = snakemake@output[[1]]
ASSAY = snakemake@params[["assay"]] %>% 
  stringr::str_replace('_', ' ')

# Read in files

pedp = readr::read_tsv(PED,
                       col_names = "SAMPLE",
                       col_types = "i",
                       col_select = 1)

phenos = readr::read_csv(PHENO_FILE) %>% 
  # filter for assay
  dplyr::filter(assay == ASSAY) %>% 
  # remove _ref or _test from `indiv`
  dplyr::mutate(indiv = stringr::str_replace_all(indiv,
                                                 c("_test" = "", "_ref" = "")))


samples = readr::read_csv(SAMPLES_FILE) %>% 
  # create `indiv` column
  tidyr::unite(date, time, quadrant,
               col = "indiv",
               sep = "_",
               remove = T) %>% 
  # rename `finclip_id` as SAMPLE to match `pedp`
  dplyr::rename(SAMPLE = finclip_id)

# Join

phen = pedp %>% 
  # join samples with ID and individual
  dplyr::left_join(.,
                   samples %>% 
                     dplyr::select(SAMPLE, indiv),
                   by = "SAMPLE") %>% 
  # join phenotype by individual
  dplyr::left_join(.,
                   phenos %>% 
                     dplyr::select(indiv, state_freq),
                   by = "indiv") %>% 
  dplyr::mutate(IID = SAMPLE) %>% 
  dplyr::select(SAMPLE, IID, state_freq) 

readr::write_tsv(phen, OUT, col_names = F)

# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
FAM = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8.fam"
SAMPLES_FILE = "config/F2_samples_converted.csv"
COVARS_QC = "config/covars_quant_cat.csv"
ASSAY = "open_field" %>% 
  stringr::str_replace('_', ' ')
COVARS = "time-quadrant" %>% 
  stringr::str_split("-", simplify = T) %>% 
  as.vector()

## True
FAM = snakemake@input[["fam"]]
SAMPLES_FILE = snakemake@input[["samples_file"]]
COVARS_QC = snakemake@input[["covars_qc"]]
OUT_CAT = snakemake@output[["cat"]]
OUT_QUANT = snakemake@output[["quant"]]
ASSAY = snakemake@params[["assay"]] %>% 
  stringr::str_replace('_', ' ')
COVARS = snakemake@params[["covars"]] %>% 
  stringr::str_split("-", simplify = T) %>% 
  as.vector()

# Split covars by quant vs cat
covars_qc = readr::read_csv(COVARS_QC)
covars_cat = covars_qc %>% 
  dplyr::filter(Q_C == "cat") %>% 
  dplyr::pull(COVAR)
covars_quant = covars_qc %>% 
  dplyr::filter(Q_C == "quant") %>% 
  dplyr::pull(COVAR)

COVARS_C = COVARS[COVARS %in% covars_cat]
COVARS_Q = COVARS[COVARS %in% covars_quant]

# Read in files

fam = genio::read_fam(FAM) %>% 
  dplyr::select(SAMPLE = id) %>% 
  dplyr::mutate(SAMPLE = as.numeric(SAMPLE))

samples = readr::read_csv(SAMPLES_FILE) %>% 
  # rename `finclip_id` as SAMPLE to match `fam`
  dplyr::rename(SAMPLE = finclip_id) %>%
  # select key columns
  dplyr::select(SAMPLE, dplyr::all_of(COVARS))

# Join

covar = fam %>% 
  # join samples with ID and individual
  dplyr::left_join(.,
                   samples,
                   by = "SAMPLE") %>% 
  dplyr::mutate(IID = SAMPLE) %>% 
  dplyr::select(SAMPLE, IID, everything()) 

# Split into categorical and quantitative variables
covar_cat = covar %>% 
  dplyr::select(SAMPLE, IID, dplyr::all_of(COVARS_C))

covar_quant = covar %>% 
  dplyr::select(SAMPLE, IID, dplyr::all_of(COVARS_Q))

readr::write_tsv(covar_cat, OUT_CAT, col_names = F)
readr::write_tsv(covar_quant, OUT_QUANT, col_names = F)

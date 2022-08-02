# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug

FAM = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8.fam"
IN_CAT = "/hps/nobackup/birney/users/ian/MIKK_HMM/covars_ms/true/time-quadrant.covar"
IN_QUANT = "/hps/nobackup/birney/users/ian/MIKK_HMM/covars_ms/true/time-quadrant.qcovar"
SEED = 1

## True
FAM = snakemake@input[["fam"]]
IN_CAT = snakemake@input[["covar_cat"]]
IN_QUANT = snakemake@input[["covar_quant"]]
SEED = snakemake@params[["seed"]] %>% 
  as.numeric()
OUT_CAT = snakemake@output[["cat"]]
OUT_QUANT = snakemake@output[["quant"]]

# Read in files

fam = genio::read_fam(FAM) %>% 
  dplyr::select(FID = id) %>% 
  dplyr::mutate(FID = as.numeric(FID))

df_cat = readr::read_tsv(IN_CAT, 
                         col_names = F) %>% 
  # rename first two columns
  dplyr::rename(FID = X1,
                IID = X2)

df_quant = readr::read_tsv(IN_QUANT, 
                          col_names = F) %>% 
  # rename first two columns
  dplyr::rename(FID = X1,
                IID = X2)

# Reorder sample

## Categorical covariates

set.seed(SEED)
df_cat$FID = sample(df_cat$FID)
set.seed(SEED)
df_cat$IID = sample(df_cat$IID)

## Quantitative covariates

set.seed(SEED)
df_quant$FID = sample(df_quant$FID)
set.seed(SEED)
df_quant$IID = sample(df_quant$IID)

# Arrange back to numerical order

out_cat = dplyr::left_join(fam,
                           df_cat,
                           by = "FID")

out_quant = dplyr::left_join(fam,
                             df_quant,
                             by = "FID")

# Save to file

readr::write_tsv(out_cat, OUT_CAT, col_names = F)
readr::write_tsv(out_quant, OUT_QUANT, col_names = F)

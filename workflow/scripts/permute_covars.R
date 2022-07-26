# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug

PED = "/hps/nobackup/birney/users/ian/MIKK_HMM/peds/F2/hdrr/5000/0.8.ped"
IN_CAT = "/hps/nobackup/birney/users/ian/MIKK_HMM/covars/true/All.covar"
IN_QUANT = "/hps/nobackup/birney/users/ian/MIKK_HMM/covars/true/All.qcovar"
SEED = 1

## True
PED = snakemake@input[["ped"]]
IN_CAT = snakemake@input[["covar_cat"]]
IN_QUANT = snakemake@input[["covar_quant"]]
SEED = snakemake@params[["seed"]] %>% 
  as.numeric()
OUT_CAT = snakemake@output[["cat"]]
OUT_QUANT = snakemake@output[["quant"]]

# Read in files

pedp = readr::read_tsv(PED,
                       col_names = "FID",
                       col_types = "i",
                       col_select = 1)

df_cat = readr::read_tsv(IN_CAT, 
                         col_names = c("FID", "IID", "DATE", "QUADRANT", "TANK_SIDE"))

df_quant = readr::read_tsv(IN_QUANT, 
                          col_names = c("FID", "IID", "TIME"))

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

out_cat = dplyr::left_join(pedp,
                           df_cat,
                           by = "FID")

out_quant = dplyr::left_join(pedp,
                             df_quant,
                             by = "FID")

# Save to file

readr::write_tsv(out_cat, OUT_CAT, col_names = F)
readr::write_tsv(out_quant, OUT_QUANT, col_names = F)

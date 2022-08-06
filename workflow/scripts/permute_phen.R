# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
FAM = "/hps/nobackup/birney/users/ian/MIKK_HMM/beds/F2/hdrr/5000/0.8.fam"
PHEN = "/hps/nobackup/birney/users/ian/MIKK_HMM/phens_ms/true/0.05/dge/invnorm/novel_object.phen"
SEED = 1

## True
FAM = snakemake@input[["fam"]]
PHEN = snakemake@input[["phen"]]
SEED = snakemake@params[["seed"]] %>% 
  as.numeric()
OUT = snakemake@output[[1]]

# Read in files

fam = genio::read_fam(FAM) %>% 
  dplyr::select(SAMPLE = id) %>% 
  dplyr::mutate(SAMPLE = as.numeric(SAMPLE))

df = readr::read_tsv(PHEN, 
                     col_names = c("FID", "IID", "PHENO"))

# Reorder sample

set.seed(SEED)
df$FID = sample(df$FID)
set.seed(SEED)
df$IID = sample(df$IID)

# Arrange back to original order

out = fam %>% 
  dplyr::left_join(df,
                   by = c("SAMPLE" = "FID")) 

# Save to file

readr::write_tsv(out, OUT, col_names = F)

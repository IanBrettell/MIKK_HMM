# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug

#IN = "/hps/nobackup/birney/users/ian/somites/covars/true/Microscope-reporter_pheno.covar"
#COVARS = "Microscope-reporter_pheno"
#SEED = 1

## True
IN = snakemake@input[[1]]
COVARS = snakemake@params[["covars"]]
SEED = snakemake@params[["seed"]] %>% 
  as.numeric()
OUT = snakemake@output[[1]]

# Parse covars

covs = stringr::str_split(COVARS, "-") %>% 
  unlist()

# Read in file

df = readr::read_tsv(IN, 
                     col_names = c("FID", "IID", covs))

# Reorder sample

set.seed(SEED)
df$FID = sample(df$FID)
set.seed(SEED)
df$IID = sample(df$IID)

# Arrange back to numerical order

out = df %>% 
  dplyr::arrange(FID)

# Save to file

readr::write_tsv(out, OUT, col_names = F)

# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Get variables

## Debug
#IN = "/hps/software/users/birney/ian/repos/somites/results/annotations_invnorm/hdrr/None/5000/0.8/intercept/snps.csv"

## True
IN = snakemake@input[[1]]
OUT = snakemake@output[["out"]]

# Read in file and process

df = readr::read_csv(IN) %>% 
  dplyr::select(CHROM, START = POS, END = POS, REF, ALT) %>% 
  # combine alleles
  tidyr::unite("REF_ALT", REF, ALT, sep = "/", remove = T) %>% 
  # add strand
  dplyr::mutate(STRAND = "+")

# Write to file

readr::write_tsv(df, OUT, col_names = F)

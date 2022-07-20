# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
IN = list("/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos/hdrr/5000/0.8/1_38-2_21-2.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos/hdrr/5000/0.8/2_38-2_21-2.csv",
          "/hps/nobackup/birney/users/ian/MIKK_HMM/F2_with_genos/hdrr/5000/0.8/99_38-2_40-1.csv")

## True
IN = snakemake@input
OUT_PED = snakemake@output[["ped"]]
OUT_MAP = snakemake@output[["map"]]

# Read in files
dat_list = purrr::map(IN, function(DF){
  df = readr::read_csv(DF)
  # Get sample
  ID = as.character(df$SAMPLE[1])
  # Remove sample column and rename GENO_NT column
  df = df %>% 
    dplyr::select(-SAMPLE, {{ID}} := GENO_NT)
})

# Join

df = dat_list %>%
  purrr::reduce(full_join, by=c("CHROM", "POS")) %>% 
  dplyr::arrange(CHROM, POS)

# Create .ped file

ped = df %>% 
  dplyr::mutate(SNP = paste(CHROM, POS, sep = ":")) %>% 
  dplyr::select(-c(CHROM, POS)) %>% 
  # pivot into 3 column
  tidyr::pivot_longer(cols = -SNP, names_to = "SAMPLE",values_to = "GT") %>% 
  tidyr::pivot_wider(id_cols = SAMPLE, names_from = SNP, values_from = GT)

# Create .map file

map = tibble::tibble(SNP = colnames(ped)) %>% 
  dplyr::filter(SNP != "SAMPLE") %>% 
  tidyr::separate(col = "SNP",into = c("CHROM", "POS"),sep = ":",remove = F) %>% 
  dplyr::select(CHROM, SNP, POS)

# Write to file

readr::write_tsv(ped, OUT_PED, col_names = F)
readr::write_tsv(map, OUT_MAP, col_names = F)

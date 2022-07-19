# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
F2 = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/F2/hdrr/hmmlearn_split/5000/0.8/1_38-2_21-2.csv"
F0 = "/hps/nobackup/birney/users/ian/MIKK_HMM/sites_files/F0/homo_divergent/38-2_21-2.txt"
BIN_LENGTH = 5000

## True
F2 = snakemake@input[["F2"]]
F0 = snakemake@input[["F0"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
OUT = snakemake@output[[1]]

# Read in data

f2 = readr::read_csv(F2)
f0 = readr::read_tsv(F0,
                     col_names = c("CHROM", "POS_1", "POS_2", "REF", "ALT", "PAT_GT", "MAT_GT")) %>% 
  dplyr::rename(POS = POS_1) %>% 
  dplyr::select(-POS_2) %>% 
  # remove CHROM == "MT"
  dplyr::filter(!CHROM == "MT") %>% 
  # convert to numberic
  dplyr::mutate(CHROM = as.numeric(CHROM)) %>% 
  # Add `BIN`
  dplyr::mutate(BIN = floor(POS / BIN_LENGTH)) %>% 
  # Get PAT and MAT alleles
  dplyr::mutate(PAT_AL = dplyr::if_else(PAT_GT == "0/0",
                                        REF,
                                        ALT),
                MAT_AL = dplyr::if_else(MAT_GT == "0/0",
                                        REF,
                                        ALT))

# Join `f0` and `f2`

final = dplyr::left_join(f0,
                         f2,
                         by = c("CHROM", "BIN")) %>% 
  # add SNP ID 
  dplyr::mutate(SNP = paste(CHROM, POS, sep = ":")) %>% 
  # get genotypes
  dplyr::mutate(GENO_NT = dplyr::case_when(STATE == 0 ~ paste(PAT_AL, PAT_AL, sep = ""),
                                           STATE == 1 ~ paste(REF, ALT, sep = ""),
                                           STATE == 2 ~ paste(MAT_AL, MAT_AL, sep = "")),
                GENO_AB = dplyr::case_when(GENO_NT == paste(REF, REF, sep = "") ~ "AA",
                                           GENO_NT == paste(REF, ALT, sep = "") ~ "AB",
                                           GENO_NT == paste(ALT, ALT, sep = "") ~ "BB")) %>% 
  # fill missing genotypes
  tidyr::fill(c(GENO_NT, GENO_AB, SAMPLE, PAT_LINE, MAT_LINE, STATE),
              .direction = "updown") %>% 
  # select key columns
  dplyr::select(SNP, CHROM, POS, REF, ALT, SAMPLE, PAT_LINE, PAT_AL, MAT_LINE, MAT_AL, STATE, GENO_NT, GENO_AB_REF_ALT = GENO_AB) %>% 
  # order by CHROM, POS
  dplyr::arrange(CHROM, POS)

# Select columns for creating .ped

out = final %>% 
  dplyr::select(SAMPLE, CHROM, POS, GENO_NT) 

# Save .ped and .map

readr::write_csv(out, OUT)

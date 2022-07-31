# Send log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug
F2 = "/hps/nobackup/birney/users/ian/MIKK_HMM/hmm_out/F2/hdrr/hmmlearn_split/5000/0.8/249_50-2_18-2.csv"
F0 = "/hps/nobackup/birney/users/ian/MIKK_HMM/genos/F0/biallelic_snps/50-2_18-2.txt"
BIN_LENGTH = 5000

## True
F2 = snakemake@input[["F2"]]
F0 = snakemake@input[["F0"]]
BIN_LENGTH = snakemake@params[["bin_length"]] %>% 
  as.numeric()
OUT_NT = snakemake@output[["nt"]]
OUT_AB = snakemake@output[["ab"]]

# Get OUT_DIRs
OUT_DIR_NT = OUT_NT[[1]] %>% 
  dirname()
OUT_DIR_AB = OUT_AB[[1]] %>% 
  dirname()

# Read in data

## F2
f2 = readr::read_csv(F2)

## F0
genos = readr::read_tsv(F0)

# Tidy up column names
colnames(genos) = colnames(genos) %>% 
  stringr::str_remove(pattern = "\\[.*\\]") %>% 
  stringr::str_remove(pattern = "# ") %>% 
  stringr::str_replace(pattern = ":", "_")
## Rename line names as `PAT` and `MAT`
col_locs = which(stringr::str_detect(colnames(genos), "GT"))
colnames(genos)[col_locs[1]] = "PAT_GT"
colnames(genos)[col_locs[2]] = "MAT_GT"

# Recode genotypes

f0 = genos %>% 
  # remove allele depth columns for Cab and Kaga
  dplyr::select(-dplyr::ends_with("AD")) %>% 
  # remove CHROM == "MT"
  dplyr::filter(!CHROM == "MT") %>% 
  # convert to numeric
  dplyr::mutate(CHROM = as.numeric(CHROM)) %>% 
  # Add `BIN`
  dplyr::mutate(BIN = floor(POS / BIN_LENGTH)) %>% 
  # recode alleles so they're in a consistent format, NA for missing or het
  dplyr::mutate(dplyr::across(dplyr::ends_with("GT"),
                              ~dplyr::case_when(. == "1|1" | . == "1/1" ~ "1/1",
                                                . == "0|0" | . == "0/0" ~ "0/0",
                                                . == "0|1" | . == "0/1" ~ NA_character_,
                                                . == ".|." | . == "./." ~ NA_character_,
                                                TRUE ~ NA_character_))) %>% 
  # filter out NAs
  tidyr::drop_na() %>% 
  # Get PAT and MAT alleles, NA if heterozygous
  dplyr::mutate(PAT_AL = dplyr::case_when(PAT_GT == "0/0" ~ REF,
                                          PAT_GT == "1/1" ~ ALT,
                                          TRUE ~ NA_character_),
                MAT_AL = dplyr::case_when(MAT_GT == "0/0" ~ REF,
                                          MAT_GT == "1/1" ~ ALT,
                                          TRUE ~ NA_character_)) %>% 
  # remove alleles with NA for either PAT_AL or MAT_AL
  dplyr::filter(!(is.na(PAT_AL) | is.na(MAT_AL)))

# Join `f0` and `f2` and fill missing genotypes

final = dplyr::inner_join(f0,
                          f2,
                          by = c("CHROM", "BIN")) %>% 
  # group by CHROM
  dplyr::group_by(CHROM) %>% 
  # fill missing states
  tidyr::fill(c(SAMPLE, PAT_LINE, MAT_LINE, STATE),
              .direction = "downup") %>%
  dplyr::ungroup() %>% 
  # add SNP ID 
  dplyr::mutate(SNP = paste(CHROM, POS, sep = ":")) %>% 
  # get genotypes
  dplyr::mutate(GENO_NT = dplyr::case_when(STATE == 0 ~ paste(PAT_AL, PAT_AL, sep = ""),
                                           # if the F2 is HET in that block but both parents share the same allele for a SNP, keep it HOM
                                           STATE == 1 & PAT_AL == MAT_AL ~ paste(PAT_AL, MAT_AL, sep = ""),
                                           # otherwise it's HET
                                           STATE == 1 & PAT_AL != MAT_AL ~ paste(REF, ALT, sep = ""),
                                           STATE == 2 ~ paste(MAT_AL, MAT_AL, sep = "")),
                GENO_AB = dplyr::case_when(GENO_NT == paste(REF, REF, sep = "") ~ "AA",
                                           GENO_NT == paste(REF, ALT, sep = "") ~ "AB",
                                           GENO_NT == paste(ALT, ALT, sep = "") ~ "BB")) %>% 
  # select key columns
  dplyr::select(SNP, CHROM, POS, REF, ALT, SAMPLE, PAT_LINE, PAT_AL, MAT_LINE, MAT_AL, STATE, GENO_NT, GENO_AB_REF_ALT = GENO_AB) %>% 
  # order by CHROM, POS
  dplyr::arrange(CHROM, POS)

# Select columns for creating .ped

out = final %>% 
  dplyr::select(SAMPLE, CHROM, POS, GENO_NT) %>% 
  # split by CHROM
  split(., f = .$CHROM)

# Write to file

## Create empty files so the job still completes even if there are no SNPs on a chrom
invisible(lapply(1:24, function(CHROM){
  file.create(file.path(OUT_DIR_NT, paste(CHROM, ".csv", sep = "")))
}))

counter = 0
invisible(lapply(out, function(DF){
  counter <<- counter + 1
  TARGET_CHROM = names(out)[counter]
  OUT_FILE = file.path(OUT_DIR_NT, paste(TARGET_CHROM, ".csv", sep = ""))
  readr::write_csv(DF, OUT_FILE)
}))


# Write AB as well

final_ab = final %>% 
  dplyr::select(SAMPLE, CHROM, POS, GENO_AB_REF_ALT) %>% 
  # split by CHROM
  split(., f = .$CHROM)

## Create empty files so the job still completes even if there are no SNPs on a chrom
invisible(lapply(1:24, function(CHROM){
  file.create(file.path(OUT_DIR_AB, paste(CHROM, ".csv", sep = "")))
}))

counter = 0
invisible(lapply(final_ab, function(DF){
  counter <<- counter + 1
  TARGET_CHROM = names(out)[counter]
  readr::write_csv(DF, file.path(OUT_DIR_AB, paste(TARGET_CHROM, ".csv", sep = "")))
}))


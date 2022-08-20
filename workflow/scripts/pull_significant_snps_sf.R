# Send output to log

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type = "message")

# Load libraries

library(tidyverse)

# Set variables

## Debug

DAT = list("/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco_state_freq_consol/true/hdrr/dist_angle/0.08/15/5000/0.8/dge/invnorm/novel_object/1/time-quadrant.mlma",
           "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco_state_freq_consol/true/hdrr/dist_angle/0.08/15/5000/0.8/dge/invnorm/novel_object/2/time-quadrant.mlma",
           "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco_state_freq_consol/true/hdrr/dist_angle/0.08/15/5000/0.8/sge/invnorm/novel_object/3/time-quadrant.mlma",
           "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco_state_freq_consol/true/hdrr/dist_angle/0.08/15/5000/0.8/sge/invnorm/novel_object/4/time-quadrant.mlma")

MIN_P = list("/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco/min_p/hdrr/dist_angle/0.08/15/5000/0.8/dge/invnorm/novel_object/1/time-quadrant.csv",
             "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco/min_p/hdrr/dist_angle/0.08/15/5000/0.8/dge/invnorm/novel_object/2/time-quadrant.csv",
             "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco/min_p/hdrr/dist_angle/0.08/15/5000/0.8/sge/invnorm/novel_object/3/time-quadrant.csv",
             "/hps/nobackup/birney/users/ian/MIKK_HMM/gcta/mlma_loco/min_p/hdrr/dist_angle/0.08/15/5000/0.8/sge/invnorm/novel_object/4/time-quadrant.csv")

ALLELES = "/hps/nobackup/birney/users/ian/MIKK_HMM/genos/F0/biallelic_snps_all/all.txt"


## True
DAT = snakemake@input[["res"]]
MIN_P = snakemake@input[["min_p"]]
ALLELES = snakemake@input[["alleles"]]
OUT_SIGS = snakemake@output[["sigs"]]

# Set names
SPLIT_DAT = purrr::map(DAT, function(FILE){
  STRING = FILE %>% 
    stringr::str_split("/", simplify = T)
  OUT = paste(STRING[, 17], STRING[, 19], STRING[,20], sep = "__")
}) %>% 
  unlist()

names(DAT) = SPLIT_DAT

SPLIT_MIN = purrr::map(MIN_P, function(FILE){
  STRING = FILE %>% 
    stringr::str_split("/", simplify = T)
  OUT = paste(STRING[, 17], STRING[, 19], STRING[,20], sep = "__")
}) %>% 
  unlist()

names(MIN_P) = SPLIT_MIN

# Data and min_p files need to correspond
# If not, produce an error
if (all(SPLIT_DAT == SPLIT_MIN) == F){
  stop()
}

# Read in alleles file
alleles = readr::read_tsv(ALLELES,
                          col_types = c("cicc"))
## Tidy up column names
colnames(alleles) = colnames(alleles) %>% 
  stringr::str_remove(pattern = "\\[.*\\]") %>% 
  stringr::str_remove(pattern = "# ") %>% 
  stringr::str_replace(pattern = ":", "_")
alleles = alleles %>% 
  # remove 'MT'
  dplyr::filter(CHROM != "MT") %>% 
  # convert CHROM to numeric
  dplyr::mutate(CHROM = as.numeric(CHROM))

# Read in full files

dat_list = purrr::map(DAT, function(MLMA){
  # Read in MLMA results
  readr::read_tsv(MLMA,
                  col_types = c("iciccdddd"))
})

# Pull out significant SNPs

sig_df = purrr::imap_dfr(dat_list, function(MLMA, i){
  # Pull minimum p-value from corresponding file
  MIN = readr::read_csv(MIN_P[[i]]) %>% 
    dplyr::pull(MIN_P) %>% 
    min(.)
  # Pull out significant SNPs from MLMA results
  OUT = MLMA %>% 
    dplyr::filter(p < MIN)
  return(OUT)
}, .id = "FILE") %>% 
  # Separate metadata from FILE column
  tidyr::separate(FILE, into = c("DGE_SGE", "ASSAY", "STATE"), sep = "__") %>% 
  # Bind with `alleles` to get REF and ALT
  dplyr::left_join(alleles %>% 
                     dplyr::select(CHROM, POS, REF, ALT),
                   by = c("Chr" = "CHROM",
                          "bp" = "POS")) %>% 
  # re-arrange column order
  dplyr::select(DGE_SGE, ASSAY, STATE, CHROM = Chr, POS = bp, REF, ALT, everything())

# Write to files

readr::write_csv(sig_df, OUT_SIGS)


